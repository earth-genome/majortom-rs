use geo::{BoundingRect, Coord, Intersects, Polygon, Rect};
use rayon::prelude::*;
use wide::{f64x4, CmpGe, CmpLe};

use crate::cell::{GridCell, GEOHASH_PRECISION};
use crate::error::GridError;
use crate::geohash::decode;

/// WGS84 equatorial radius in metres.
const EARTH_RADIUS: f64 = 6378137.0;

/// Epsilon used when expanding the row/column search bounds, matching the Go
/// and Python reference implementations.
const EPSILON: f64 = 1e-10;

/// Row count at which [`generate_grid_cells`](MajorTomGrid::generate_grid_cells)
/// switches from a sequential scan to rayon. Below this, thread-pool overhead
/// dominates for typical small AOIs.
const PARALLEL_ROW_THRESHOLD: i64 = 32;

/// AOI bounding box used to bound the row/column search and for the cheap AABB
/// prefilter. `max_lon` may be expanded by +360° when the AOI crosses the
/// antimeridian (matching the Go/Python ports).
#[derive(Clone, Copy)]
struct AoiBounds {
    min_lon: f64,
    min_lat: f64,
    max_lon: f64,
    max_lat: f64,
}

/// Returns true when the axis-aligned cell bbox overlaps the AOI bbox.
///
/// Matches the inline pre-check in the Go `mtgrid` implementation.
#[inline]
fn aabb_overlaps(min_lon: f64, min_lat: f64, max_lon: f64, max_lat: f64, aoi: AoiBounds) -> bool {
    !(max_lon < aoi.min_lon
        || min_lon > aoi.max_lon
        || max_lat < aoi.min_lat
        || min_lat > aoi.max_lat)
}

/// SIMD longitude AABB mask for four candidate cells that already share a
/// latitude range known to overlap the AOI.
///
/// Returns per-lane booleans (via the cmp mask's `to_array`), avoiding
/// `move_mask` bit-order ambiguity across AVX / non-AVX backends.
#[inline]
fn lon_aabb_lanes(min_lons: f64x4, max_lons: f64x4, aoi: AoiBounds) -> [bool; 4] {
    let ge = max_lons.cmp_ge(f64x4::splat(aoi.min_lon));
    let le = min_lons.cmp_le(f64x4::splat(aoi.max_lon));
    let bits = (ge & le).to_array();
    [
        bits[0].to_bits() != 0,
        bits[1].to_bits() != 0,
        bits[2].to_bits() != 0,
        bits[3].to_bits() != 0,
    ]
}

/// Returns true when the axis-aligned cell bbox intersects `polygon`.
///
/// Uses [`Rect`] so we can reject candidates without allocating a
/// [`GridCell`] (geohash encode + polygon ring).
#[inline]
fn cell_intersects(
    polygon: &Polygon<f64>,
    min_lon: f64,
    min_lat: f64,
    max_lon: f64,
    max_lat: f64,
) -> bool {
    Rect::new(
        Coord {
            x: min_lon,
            y: min_lat,
        },
        Coord {
            x: max_lon,
            y: max_lat,
        },
    )
    .intersects(polygon)
}

/// Push a cell when it passes the AOI AABB prefilter and a precise `Rect`
/// intersection against `polygon`.
#[inline]
fn maybe_push_cell(
    cells: &mut Vec<GridCell>,
    polygon: &Polygon<f64>,
    aoi: AoiBounds,
    min_lon: f64,
    min_lat: f64,
    max_lon: f64,
    max_lat: f64,
    is_primary: bool,
) {
    if !aabb_overlaps(min_lon, min_lat, max_lon, max_lat, aoi) {
        return;
    }
    if cell_intersects(polygon, min_lon, min_lat, max_lon, max_lat) {
        cells.push(GridCell::from_bbox(
            min_lon, min_lat, max_lon, max_lat, is_primary,
        ));
    }
}

/// An equal-area Major TOM grid.
///
/// Construct one with [`MajorTomGrid::new`] and then query it with
/// [`generate_grid_cells`](MajorTomGrid::generate_grid_cells),
/// [`cell_from_id`](MajorTomGrid::cell_from_id) or
/// [`migrate_cell_id`](MajorTomGrid::migrate_cell_id).
#[derive(Debug, Clone)]
pub struct MajorTomGrid {
    d: f64,
    overlap: bool,
    row_count: i64,
    lat_spacing: f64,
    lat_offset: f64,
}

impl MajorTomGrid {
    /// Creates a new grid with cell edge length `d` (in metres).
    ///
    /// Returns [`GridError::InvalidSpacing`] when `d == 0`.
    pub fn new(d: u64, overlap: bool) -> Result<Self, GridError> {
        if d == 0 {
            return Err(GridError::InvalidSpacing);
        }
        let d = d as f64;
        let row_count = (std::f64::consts::PI * EARTH_RADIUS / d).ceil().max(2.0) as i64;
        let lat_spacing = (180.0 / row_count as f64).min(89.0);
        let lat_offset = if row_count % 2 != 0 {
            lat_spacing / 2.0
        } else {
            0.0
        };
        Ok(MajorTomGrid {
            d,
            overlap,
            row_count,
            lat_spacing,
            lat_offset,
        })
    }

    /// Whether the grid emits half-spacing overlap cells.
    pub fn overlap(&self) -> bool {
        self.overlap
    }

    /// Number of latitude rows spanning the globe.
    pub fn row_count(&self) -> i64 {
        self.row_count
    }

    /// The latitude spacing between rows, in degrees.
    pub fn lat_spacing(&self) -> f64 {
        self.lat_spacing
    }

    /// The latitude centring offset, in degrees.
    pub fn lat_offset(&self) -> f64 {
        self.lat_offset
    }

    /// Returns the latitude (degrees) of the south edge of the given row.
    pub fn row_lat(&self, row_idx: i64) -> f64 {
        -90.0 + self.lat_offset + row_idx as f64 * self.lat_spacing
    }

    /// Returns the longitude spacing (degrees) for cells at the given latitude.
    pub fn lon_spacing(&self, lat: f64) -> f64 {
        let lat_rad = lat.clamp(-89.0, 89.0).to_radians();
        let circumference = 2.0 * std::f64::consts::PI * EARTH_RADIUS * lat_rad.cos();
        let n_cols = (circumference / self.d).ceil();
        360.0 / n_cols.max(1.0)
    }

    /// Returns the longitude centring offset for a given longitude spacing.
    pub fn lon_offset(&self, lon_spacing: f64) -> f64 {
        let n_cols = if lon_spacing > 0.0 {
            (360.0 / lon_spacing).round() as i64
        } else {
            0
        };
        if n_cols % 2 != 0 {
            lon_spacing / 2.0
        } else {
            0.0
        }
    }

    /// Returns the longitude (degrees) of the west edge of the given column.
    pub fn col_lon(&self, col_idx: i64, lon_spacing: f64, lon_offset: f64) -> f64 {
        -180.0 + lon_offset + col_idx as f64 * lon_spacing
    }

    /// Generates every primary (and, if enabled, overlap) cell that
    /// geometrically intersects `polygon`.
    ///
    /// Rows are scanned in parallel (via rayon) once the AOI spans at least
    /// [`PARALLEL_ROW_THRESHOLD`] latitude rows. Within each row, candidate
    /// longitudes are AABB-filtered in SIMD batches of four before a precise
    /// [`Rect`] intersection and [`GridCell`] construction.
    pub fn generate_grid_cells(&self, polygon: &Polygon<f64>) -> Vec<GridCell> {
        let bounds = match polygon.bounding_rect() {
            Some(b) => b,
            None => return Vec::new(),
        };
        let min_lon = bounds.min().x;
        let min_lat = bounds.min().y;
        let mut max_lon = bounds.max().x;
        let max_lat = bounds.max().y;
        if min_lon > max_lon {
            max_lon += 360.0;
        }
        let aoi = AoiBounds {
            min_lon,
            min_lat,
            max_lon,
            max_lat,
        };

        let mut start_row = ((min_lat + 90.0 - self.lat_offset) / self.lat_spacing).floor() as i64;
        let mut end_row = ((max_lat + 90.0 - self.lat_offset) / self.lat_spacing).ceil() as i64;
        while self.row_lat(start_row) > min_lat + EPSILON {
            start_row -= 1;
        }
        while self.row_lat(end_row) < max_lat - EPSILON {
            end_row += 1;
        }

        let half_lat_spacing = self.lat_spacing / 2.0;
        let row_count = end_row - start_row + 1;

        if row_count >= PARALLEL_ROW_THRESHOLD {
            (start_row..=end_row)
                .into_par_iter()
                .map(|row_idx| self.cells_for_row(row_idx, polygon, aoi, half_lat_spacing))
                .reduce(Vec::new, |mut acc, mut row| {
                    acc.append(&mut row);
                    acc
                })
        } else {
            let mut cells = Vec::with_capacity((row_count as usize).saturating_mul(2));
            for row_idx in start_row..=end_row {
                cells.extend(self.cells_for_row(row_idx, polygon, aoi, half_lat_spacing));
            }
            cells
        }
    }

    /// Collects intersecting cells for a single latitude row.
    fn cells_for_row(
        &self,
        row_idx: i64,
        polygon: &Polygon<f64>,
        aoi: AoiBounds,
        half_lat_spacing: f64,
    ) -> Vec<GridCell> {
        let lat = self.row_lat(row_idx);
        let lon_spacing = self.lon_spacing(lat);
        let lon_offset = self.lon_offset(lon_spacing);
        let half_lon_spacing = lon_spacing / 2.0;
        let cell_max_lat = lat + self.lat_spacing;

        let mut start_col = ((aoi.min_lon + 180.0 - lon_offset) / lon_spacing).floor() as i64;
        let mut end_col = ((aoi.max_lon + 180.0 - lon_offset) / lon_spacing).ceil() as i64;
        while self.col_lon(start_col, lon_spacing, lon_offset) > aoi.min_lon + EPSILON {
            start_col -= 1;
        }
        while self.col_lon(end_col, lon_spacing, lon_offset) < aoi.max_lon - EPSILON {
            end_col += 1;
        }

        let primary_lat_ok = cell_max_lat >= aoi.min_lat && lat <= aoi.max_lat;
        let overlap_lat = lat + half_lat_spacing;
        let overlap_max_lat = overlap_lat + self.lat_spacing;
        let overlap_lat_ok =
            self.overlap && overlap_max_lat >= aoi.min_lat && overlap_lat <= aoi.max_lat;

        let col_count = (end_col - start_col + 1).max(0) as usize;
        let mut cells =
            Vec::with_capacity(col_count.saturating_mul(if self.overlap { 2 } else { 1 }));

        let mut col_idx = start_col;
        while col_idx <= end_col {
            let remaining = end_col - col_idx + 1;
            if remaining >= 4 && (primary_lat_ok || overlap_lat_ok) {
                // West edges from col_lon (exact match to the scalar path).
                // SIMD is used only for the AABB reject mask; cell construction
                // stays scalar so centroids/geohashes stay bit-identical.
                let min_arr = [
                    self.col_lon(col_idx, lon_spacing, lon_offset),
                    self.col_lon(col_idx + 1, lon_spacing, lon_offset),
                    self.col_lon(col_idx + 2, lon_spacing, lon_offset),
                    self.col_lon(col_idx + 3, lon_spacing, lon_offset),
                ];
                let min_lons = f64x4::new(min_arr);
                let max_lons = min_lons + f64x4::splat(lon_spacing);

                if primary_lat_ok {
                    let lanes = lon_aabb_lanes(min_lons, max_lons, aoi);
                    for lane in 0..4 {
                        if !lanes[lane] {
                            continue;
                        }
                        let lon = min_arr[lane];
                        let cell_max_lon = lon + lon_spacing;
                        if cell_intersects(polygon, lon, lat, cell_max_lon, cell_max_lat) {
                            cells.push(GridCell::from_bbox(
                                lon,
                                lat,
                                cell_max_lon,
                                cell_max_lat,
                                true,
                            ));
                        }
                    }
                }

                if overlap_lat_ok {
                    let overlap_mins = min_lons + f64x4::splat(half_lon_spacing);
                    let overlap_maxs = overlap_mins + f64x4::splat(lon_spacing);
                    let lanes = lon_aabb_lanes(overlap_mins, overlap_maxs, aoi);
                    for lane in 0..4 {
                        if !lanes[lane] {
                            continue;
                        }
                        let lon = min_arr[lane] + half_lon_spacing;
                        let cell_max_lon = lon + lon_spacing;
                        if cell_intersects(polygon, lon, overlap_lat, cell_max_lon, overlap_max_lat)
                        {
                            cells.push(GridCell::from_bbox(
                                lon,
                                overlap_lat,
                                cell_max_lon,
                                overlap_max_lat,
                                false,
                            ));
                        }
                    }
                }

                col_idx += 4;
            } else {
                let lon = self.col_lon(col_idx, lon_spacing, lon_offset);
                let cell_max_lon = lon + lon_spacing;

                if primary_lat_ok {
                    maybe_push_cell(
                        &mut cells,
                        polygon,
                        aoi,
                        lon,
                        lat,
                        cell_max_lon,
                        cell_max_lat,
                        true,
                    );
                }

                if overlap_lat_ok {
                    let overlap_lon = lon + half_lon_spacing;
                    maybe_push_cell(
                        &mut cells,
                        polygon,
                        aoi,
                        overlap_lon,
                        overlap_lat,
                        overlap_lon + lon_spacing,
                        overlap_max_lat,
                        false,
                    );
                }

                col_idx += 1;
            }
        }

        cells
    }

    /// Reconstructs the [`GridCell`] identified by `id`.
    ///
    /// IDs longer than 11 characters are truncated. The cell is located by
    /// decoding the geohash centre and searching the 3×3 row/column
    /// neighbourhood to absorb floating-point edge cases (and to find overlap
    /// cells).
    pub fn cell_from_id(&self, id: &str) -> Result<GridCell, GridError> {
        let search_id: &str = if id.len() > GEOHASH_PRECISION {
            &id[..GEOHASH_PRECISION]
        } else {
            id
        };
        if search_id.len() != GEOHASH_PRECISION {
            return Err(GridError::InvalidCellId);
        }

        let (center_lat, center_lon) = decode(search_id)?;

        let half_lat = self.lat_spacing / 2.0;
        for row_offset in [0_i64, -1, 1] {
            let row_idx = ((center_lat + 90.0 - self.lat_offset) / self.lat_spacing).floor() as i64
                + row_offset;
            let row_lat = self.row_lat(row_idx);
            let lon_spacing = self.lon_spacing(row_lat);
            let lon_offset = self.lon_offset(lon_spacing);
            let half_lon = lon_spacing / 2.0;

            for col_offset in [0_i64, -1, 1] {
                let col_idx =
                    ((center_lon + 180.0 - lon_offset) / lon_spacing).floor() as i64 + col_offset;
                let cell_lon = self.col_lon(col_idx, lon_spacing, lon_offset);

                let primary = GridCell::from_bbox(
                    cell_lon,
                    row_lat,
                    cell_lon + lon_spacing,
                    row_lat + self.lat_spacing,
                    true,
                );
                if primary.id() == search_id {
                    return Ok(primary);
                }

                if self.overlap {
                    let overlap_lon = cell_lon + half_lon;
                    let overlap_lat = row_lat + half_lat;
                    let overlap = GridCell::from_bbox(
                        overlap_lon,
                        overlap_lat,
                        overlap_lon + lon_spacing,
                        overlap_lat + self.lat_spacing,
                        false,
                    );
                    if overlap.id() == search_id {
                        return Ok(overlap);
                    }
                }
            }
        }

        Err(GridError::CellNotFound(id.to_string()))
    }

    /// Maps a cell ID from a prior grid version onto the current grid.
    ///
    /// Decodes the old geohash to recover an approximate centroid and returns
    /// the current grid's **primary** cell that contains that point. The ID
    /// must be at least 11 characters.
    pub fn migrate_cell_id(&self, old_id: &str) -> Result<GridCell, GridError> {
        let search_id: &str = if old_id.len() > GEOHASH_PRECISION {
            &old_id[..GEOHASH_PRECISION]
        } else {
            old_id
        };
        if search_id.len() != GEOHASH_PRECISION {
            return Err(GridError::InvalidCellId);
        }

        let (lat, lon) = decode(search_id)?;

        let row_idx = ((lat + 90.0 - self.lat_offset) / self.lat_spacing).floor() as i64;
        let row_lat = self.row_lat(row_idx);
        let lon_spacing = self.lon_spacing(row_lat);
        let lon_offset = self.lon_offset(lon_spacing);
        let col_idx = ((lon + 180.0 - lon_offset) / lon_spacing).floor() as i64;
        let cell_lon = self.col_lon(col_idx, lon_spacing, lon_offset);

        Ok(GridCell::from_bbox(
            cell_lon,
            row_lat,
            cell_lon + lon_spacing,
            row_lat + self.lat_spacing,
            true,
        ))
    }
}
