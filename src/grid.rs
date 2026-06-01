use geo::{BoundingRect, Intersects, Polygon};

use crate::cell::{GridCell, GEOHASH_PRECISION};
use crate::error::GridError;
use crate::geohash::decode;

/// WGS84 equatorial radius in metres.
const EARTH_RADIUS: f64 = 6378137.0;

/// Epsilon used when expanding the row/column search bounds, matching the Go
/// and Python reference implementations.
const EPSILON: f64 = 1e-10;

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

        let mut start_row = ((min_lat + 90.0 - self.lat_offset) / self.lat_spacing).floor() as i64;
        let mut end_row = ((max_lat + 90.0 - self.lat_offset) / self.lat_spacing).ceil() as i64;
        while self.row_lat(start_row) > min_lat + EPSILON {
            start_row -= 1;
        }
        while self.row_lat(end_row) < max_lat - EPSILON {
            end_row += 1;
        }

        let half_lat_spacing = self.lat_spacing / 2.0;
        let mut cells = Vec::new();

        for row_idx in start_row..=end_row {
            let lat = self.row_lat(row_idx);
            let lon_spacing = self.lon_spacing(lat);
            let lon_offset = self.lon_offset(lon_spacing);
            let half_lon_spacing = lon_spacing / 2.0;

            let mut start_col = ((min_lon + 180.0 - lon_offset) / lon_spacing).floor() as i64;
            let mut end_col = ((max_lon + 180.0 - lon_offset) / lon_spacing).ceil() as i64;
            while self.col_lon(start_col, lon_spacing, lon_offset) > min_lon + EPSILON {
                start_col -= 1;
            }
            while self.col_lon(end_col, lon_spacing, lon_offset) < max_lon - EPSILON {
                end_col += 1;
            }

            for col_idx in start_col..=end_col {
                let lon = self.col_lon(col_idx, lon_spacing, lon_offset);

                let primary =
                    GridCell::from_bbox(lon, lat, lon + lon_spacing, lat + self.lat_spacing, true);
                if primary.geom.intersects(polygon) {
                    cells.push(primary);
                }

                if self.overlap {
                    let overlap_lon = lon + half_lon_spacing;
                    let overlap_lat = lat + half_lat_spacing;
                    let overlap = GridCell::from_bbox(
                        overlap_lon,
                        overlap_lat,
                        overlap_lon + lon_spacing,
                        overlap_lat + self.lat_spacing,
                        false,
                    );
                    if overlap.geom.intersects(polygon) {
                        cells.push(overlap);
                    }
                }
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
