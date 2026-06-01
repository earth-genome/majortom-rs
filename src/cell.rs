use geo::{Coord, LineString, Polygon};

use crate::geohash::encode;

/// The geohash precision used for every cell ID. Matches the Python
/// (`geolib`) and Go (`pierrre/geohash`) reference implementations.
pub(crate) const GEOHASH_PRECISION: usize = 11;

/// A single cell of the Major TOM grid.
///
/// A cell is an axis-aligned rectangle in geographic (lon/lat) coordinates and
/// carries a stable identifier: the precision-11 geohash of its centroid.
#[derive(Debug, Clone, PartialEq)]
pub struct GridCell {
    /// The cell geometry as a closed rectangular polygon.
    pub geom: Polygon<f64>,
    /// Whether this is a primary cell (`true`) or a half-spacing overlap cell
    /// (`false`).
    pub is_primary: bool,
    id: String,
}

impl GridCell {
    /// Builds a rectangular cell from its bounding box and pre-computes the
    /// geohash-11 ID from the bounding-box centre.
    ///
    /// The vertex order matches the Go and Python reference implementations
    /// (counter-clockwise starting at the south-west corner).
    pub(crate) fn from_bbox(
        min_lon: f64,
        min_lat: f64,
        max_lon: f64,
        max_lat: f64,
        is_primary: bool,
    ) -> Self {
        let center_lon = (min_lon + max_lon) / 2.0;
        let center_lat = (min_lat + max_lat) / 2.0;
        let id = encode(center_lat, center_lon, GEOHASH_PRECISION);

        let ring = LineString::new(vec![
            Coord {
                x: min_lon,
                y: min_lat,
            },
            Coord {
                x: max_lon,
                y: min_lat,
            },
            Coord {
                x: max_lon,
                y: max_lat,
            },
            Coord {
                x: min_lon,
                y: max_lat,
            },
            Coord {
                x: min_lon,
                y: min_lat,
            },
        ]);

        GridCell {
            geom: Polygon::new(ring, vec![]),
            is_primary,
            id,
        }
    }

    /// Returns the cell's stable identifier: the precision-11 geohash of its
    /// centroid.
    pub fn id(&self) -> &str {
        &self.id
    }
}
