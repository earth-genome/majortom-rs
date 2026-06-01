//! A Rust implementation of the ESA [Major TOM](https://github.com/ESA-PhiLab/Major-TOM)
//! equal-area grid.
//!
//! This crate ports the logic from the
//! [`earth-genome/majortom`](https://github.com/earth-genome/majortom) (Python) and
//! [`earth-genome/mtgrid`](https://github.com/earth-genome/mtgrid) (Go) libraries, with
//! the explicit goal of **behavioural parity**: for the same input polygon and
//! configuration it produces the exact same set of grid cells (geometry +
//! geohash ID), verified against the shared `cross_language_reference.json`
//! fixture.
//!
//! # Example
//!
//! ```
//! use geo::{Coord, LineString, Polygon};
//! use majortom::MajorTomGrid;
//!
//! // A 1/10th-degree square near (0, 0).
//! let aoi = Polygon::new(
//!     LineString::new(vec![
//!         Coord { x: 0.0, y: 0.0 },
//!         Coord { x: 0.0, y: 0.1 },
//!         Coord { x: 0.1, y: 0.1 },
//!         Coord { x: 0.1, y: 0.0 },
//!         Coord { x: 0.0, y: 0.0 },
//!     ]),
//!     vec![],
//! );
//!
//! let grid = MajorTomGrid::new(320, true).unwrap();
//! let cells = grid.generate_grid_cells(&aoi);
//! assert!(!cells.is_empty());
//! for cell in &cells {
//!     // Each cell carries a stable precision-11 geohash ID.
//!     assert_eq!(cell.id().len(), 11);
//! }
//! ```

mod cell;
mod error;
mod geohash;
mod grid;

pub use cell::GridCell;
pub use error::GridError;
pub use grid::MajorTomGrid;
