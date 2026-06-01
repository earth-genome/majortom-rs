# majortom-rs

A Rust implementation of the ESA [Major TOM](https://github.com/ESA-PhiLab/Major-TOM)
equal-area grid.

This crate ports the logic from the
[`earth-genome/majortom`](https://github.com/earth-genome/majortom) (Python) and
[`earth-genome/mtgrid`](https://github.com/earth-genome/mtgrid) (Go) libraries,
which in turn are based on the ESA-PhiLab
[grid implementation](https://github.com/ESA-PhiLab/Major-TOM/blob/main/src/grid.py).

The goal is **behavioural parity**: for the same input polygon and configuration
this crate produces the exact same set of grid cells (geometry + geohash ID) as
the Python and Go libraries. This is verified against the shared
`testdata/cross_language_reference.json` fixture (copied verbatim from `mtgrid`).

## Usage

```rust
use geo::{Coord, LineString, Polygon};
use majortom::MajorTomGrid;

fn main() {
    // A 1/10th-degree square near (0, 0).
    let aoi = Polygon::new(
        LineString::new(vec![
            Coord { x: 0.0, y: 0.0 },
            Coord { x: 0.0, y: 0.1 },
            Coord { x: 0.1, y: 0.1 },
            Coord { x: 0.1, y: 0.0 },
            Coord { x: 0.0, y: 0.0 },
        ]),
        vec![],
    );

    // A grid with 320 m cells and overlap cells enabled.
    let grid = MajorTomGrid::new(320, true).unwrap();

    for cell in grid.generate_grid_cells(&aoi) {
        println!("cell id: {}", cell.id());
    }
}
```

## API

- `MajorTomGrid::new(d: u64, overlap: bool)` — construct a grid with cell edge
  length `d` (metres). Returns an error when `d == 0`.
- `generate_grid_cells(&self, polygon: &geo::Polygon<f64>) -> Vec<GridCell>` —
  every primary (and, if enabled, overlap) cell that geometrically intersects
  the polygon.
- `cell_from_id(&self, id: &str) -> Result<GridCell, GridError>` — reconstruct a
  cell from its geohash ID (IDs longer than 11 characters are truncated).
- `migrate_cell_id(&self, old_id: &str) -> Result<GridCell, GridError>` — map a
  cell ID from a prior grid version onto the current grid's primary cell.

Each `GridCell` exposes its geometry (`geom`), an `is_primary` flag, and a
cached precision-11 geohash `id()`.

## Geohash parity

Cell IDs are the precision-11 geohash of the cell centroid. To match the Python
(`geolib`) and Go (`pierrre/geohash`) references byte-for-byte — including cells
straddling the equator and prime meridian, where the centroid sits
infinitesimally close to a cell boundary — this crate vendors a small
classic interval-bisection geohash encoder/decoder (`src/geohash.rs`) rather
than depending on the `geohash` crate, whose floating-point bit trick rounds
such sub-ULP coordinates onto the wrong side of the boundary.

## Development

```bash
cargo test          # unit + integration + cross-language conformance tests
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo bench         # optional criterion benchmarks
```

## License

MIT — see [LICENSE](LICENSE).
