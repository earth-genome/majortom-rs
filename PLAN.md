# Implementation Plan: `majortom-rs`

A Rust implementation of the ESA [Major TOM](https://github.com/ESA-PhiLab/Major-TOM)
equal-area grid, porting the logic from the existing
[`earth-genome/majortom`](https://github.com/earth-genome/majortom) (Python) and
[`earth-genome/mtgrid`](https://github.com/earth-genome/mtgrid) (Go) libraries.

The goal is **behavioural parity**: for the same input polygon and configuration,
this crate must produce the exact same set of grid cells (geometry + geohash ID)
as the Python and Go libraries, verified against the shared
`cross_language_reference.json` fixture.

---

## 1. What the reference libraries do

Both libraries implement the same algorithm. The Go version is the most recent and
includes a cross-language conformance fixture, so it is the primary reference; the
Python version is the canonical source of truth for numeric values.

### Core constants
- `radius = 6378137` (WGS84 equatorial radius, metres)
- `geohash precision = 11`

### Grid construction (`new(d, overlap)`)
Given a cell edge length `d` (metres) and an `overlap` flag:
- `row_count = max(2, ceil(pi * radius / d))`
- `lat_spacing = min(180 / row_count, 89)` (degrees)
- `lat_offset = lat_spacing / 2` if `row_count` is odd, else `0`
  (this centring makes the equator land on a grid line)

### Per-row / per-column geometry
- `row_lat(row_idx) = -90 + lat_offset + row_idx * lat_spacing`
- `lon_spacing(lat)`:
  - `lat_rad = radians(clamp(lat, -89, 89))`
  - `circumference = 2 * pi * radius * cos(lat_rad)`
  - `n_cols = ceil(circumference / d)`
  - `lon_spacing = 360 / max(n_cols, 1)`
- `lon_offset(lon_spacing)`:
  - `n_cols = round(360 / lon_spacing)` (when `lon_spacing > 0`)
  - `lon_spacing / 2` if `n_cols` is odd, else `0`
- `col_lon(col_idx, lon_spacing, lon_offset) = -180 + lon_offset + col_idx * lon_spacing`

### `generate_grid_cells(polygon)`
1. Take the polygon's bounding box `(min_lon, min_lat, max_lon, max_lat)`.
   If `min_lon > max_lon`, add `360` to `max_lon` (antimeridian wrap).
2. Compute `start_row`/`end_row` from `(lat + 90 - lat_offset) / lat_spacing`
   (floor / ceil), then expand outward with `while` loops until the row latitudes
   fully bracket `[min_lat, max_lat]` (epsilon `1e-10`).
3. For each row, compute `lon_spacing`, `lon_offset`, and `start_col`/`end_col`
   analogously, expanding with `while` loops.
4. For each `(row, col)` build the **primary** cell rectangle
   `[lon, lat] -> [lon+lon_spacing, lat+lat_spacing]`.
5. If `overlap`, also build an **overlap** cell shifted by
   `(+lon_spacing/2, +lat_spacing/2)`.
6. Emit a cell only if it geometrically **intersects** the input polygon
   (true polygon intersection, not just bbox overlap).

### `GridCell` ID
The cell ID is the **geohash (precision 11)** of the cell centroid. For an
axis-aligned rectangle the centroid equals the bounding-box centre, which is what
Go uses (`Bound().Center()`); Python uses Shapely's polygon centroid — identical
for these rectangles.

### `cell_from_id(id)`
- Truncate IDs longer than 11 chars to 11.
- Decode the geohash to its centre `(lat, lon)`.
- Compute the nominal `(row, col)` directly, then search the 3×3 neighbourhood
  (`row_offset, col_offset ∈ {0, -1, +1}`), reconstructing primary (and, if
  `overlap`, overlap) cells and returning the one whose geohash ID matches.
  The neighbour search handles floating-point edge cases.

### `migrate_cell_id(old_id)`
Decode an old-grid geohash to a point, then return the **primary** cell of the
current grid that contains that point (single direct row/col computation, no search).

### Key parity nuances to preserve
- **Iteration bounds differ between Python and Go** (Python uses an inclusive
  `range(start, end+1)`; Go uses exclusive `rowIdx < endRow`). Both yield the same
  *emitted* set because the intersection filter drops the boundary cells. The Rust
  port must reproduce the **emitted set**, validated against the fixture — not a
  particular loop convention.
- Latitude clamping to `±89` inside `lon_spacing` only.
- `max(n_cols, 1)` and `max(2, ...)` guards.
- Epsilon `1e-10` in the row/col expansion loops.
- Antimeridian `max_lon += 360` handling.

---

## 2. Target dependencies (Rust crates)

| Concern | Crate | Notes |
|---|---|---|
| Geometry types & predicates | [`geo`](https://crates.io/crates/geo) (`geo-types` + `Intersects`) | `Polygon<f64>`, `Rect`, `Coord`; `Intersects` trait for precise polygon intersection (Shapely/orb-predicates equivalent). |
| Geohash encode/decode | [`geohash`](https://crates.io/crates/geohash) | `encode(Coord, len)` and `decode(&str) -> (Coord, lon_err, lat_err)`; `decode` returns the cell centre directly. Standard base32 geohash — must verify byte-for-byte parity vs `geolib`/`pierrre`. |
| Errors | [`thiserror`](https://crates.io/crates/thiserror) | Typed error enum for the public API. |
| (dev) GeoJSON / fixtures | [`serde`](https://crates.io/crates/serde) + [`serde_json`](https://crates.io/crates/serde_json) | Parse `cross_language_reference.json` for conformance tests. |

**Verification gate:** before building anything else, write a tiny spike that
encodes a few known `(lat, lon)` points at precision 11 and confirms the output
matches the IDs in `cross_language_reference.json` (e.g. `dr18yzvf3fv`). If the
`geohash` crate's base32 ordering or rounding differs, we vendor a small geohash
implementation matching `geolib`/`pierrre` exactly. This is the single biggest
cross-language risk and must be retired first.

---

## 3. Proposed crate layout

```
majortom-rs/
├── Cargo.toml
├── README.md                         # usage, ported from Go/Python READMEs
├── PLAN.md                           # this file
├── src/
│   ├── lib.rs                        # public re-exports + crate docs
│   ├── grid.rs                       # MajorTomGrid: new, spacing/offset/row/col helpers,
│   │                                 #   generate_grid_cells, cell_from_id, migrate_cell_id
│   ├── cell.rs                       # GridCell: geometry + cached geohash id, is_primary
│   └── error.rs                      # GridError (invalid spacing, cell not found, bad id)
├── tests/
│   ├── grid_tests.rs                 # ported unit tests (counts, overlap on/off, lookup, edge cases)
│   ├── esa_alignment_tests.rs        # latitude/longitude alignment vs ESA linspace+mod reference
│   ├── migrate_tests.rs              # migrate_cell_id behaviour
│   └── cross_language_tests.rs       # parity vs testdata/cross_language_reference.json
├── testdata/
│   └── cross_language_reference.json # copied verbatim from earth-genome/mtgrid
└── benches/
    └── grid_bench.rs                 # optional criterion benchmarks (generate / cell_from_id / id)
```

## 4. Public API (proposed)

```rust
pub struct GridCell {
    pub geom: geo::Polygon<f64>,
    pub is_primary: bool,
    // cached geohash-11 id
}
impl GridCell {
    pub fn id(&self) -> &str;
}

pub struct MajorTomGrid { /* d, overlap, row_count, lat_spacing, lat_offset */ }

impl MajorTomGrid {
    pub fn new(d: u64, overlap: bool) -> Result<Self, GridError>; // d must be > 0

    /// Returns all primary (+ overlap) cells intersecting `polygon`.
    pub fn generate_grid_cells(&self, polygon: &geo::Polygon<f64>)
        -> Vec<GridCell>;

    pub fn cell_from_id(&self, id: &str) -> Result<GridCell, GridError>;
    pub fn migrate_cell_id(&self, old_id: &str) -> Result<GridCell, GridError>;
}
```

API/idiom decisions (chosen, not blocking):
- `generate_grid_cells` returns a `Vec<GridCell>` to mirror Go. A streaming
  `Iterator` variant (mirroring Python's generator) can be added later if needed.
- `d: u64` matches Go's `uint64`; `new` returns `Err` on `d == 0` (Python raises on
  `d <= 0`).
- Accept `&geo::Polygon<f64>`; a convenience that also accepts `Rect`/bbox can be
  layered on. No need for the Go `derefGeometry` pointer-unwrapping hack.

## 5. Implementation steps

1. **Scaffold**: `cargo init --lib`, set crate name `majortom`, edition 2021, add
   dependencies, fill `Cargo.toml` metadata (license MIT, repo URL).
2. **Geohash parity spike** (Section 2 gate): confirm the `geohash` crate matches
   the fixture IDs at precision 11; decide crate-vs-vendor.
3. **`error.rs`**: `GridError` enum (`InvalidSpacing`, `InvalidCellId`, `CellNotFound`,
   `GeohashError`).
4. **`cell.rs`**: `GridCell` with rectangle constructor that computes centroid →
   geohash-11 id once and caches it; `is_primary` flag; `id()` accessor.
5. **`grid.rs` — construction & helpers**: `new`, `row_lat`, `lon_spacing`,
   `lon_offset`, `col_lon` exactly per Section 1. Keep helpers `pub(crate)` (or
   `pub`) so alignment tests can call them like the Go/Python tests do.
6. **`grid.rs` — `generate_grid_cells`**: bbox extraction (antimeridian wrap), row
   then column expansion loops, primary + overlap cell construction, precise
   `Intersects` filtering. (Single-threaded first for correctness; optional
   parallelism via `rayon` later — Go parallelises rows but ordering doesn't matter
   since results are compared as sets.)
7. **`grid.rs` — `cell_from_id`**: truncate-to-11, decode, 3×3 neighbour search.
8. **`grid.rs` — `migrate_cell_id`**: decode + single direct row/col → primary cell;
   validate ID length (`>= 11`, error otherwise).
9. **`lib.rs`**: re-export `MajorTomGrid`, `GridCell`, `GridError`; crate-level docs
   with the usage example.
10. **README**: port usage from the Go/Python READMEs to Rust.

## 6. Testing & verification strategy

Port the existing test suites (Python `tests/test_major_tom_grid.py`, Go
`mtgrid_test.go`) and add the cross-language fixture as the authoritative gate:

- **Unit / behaviour** (`grid_tests.rs`):
  - `bigSouthampton` polygon with `d=320, overlap=true` → **exactly 225 cells**
    (matches Python `test_generate_grid_cells_basic`); all cells intersect the AOI.
  - `overlap=false` yields strictly fewer cells than `overlap=true`.
  - Tiny polygon near `(0,0)` → ≥ 1 cell.
  - High-latitude polygon `(170..171, 80..81)` → ≥ 1 cell.
  - Round-trip: every generated cell's id resolves via `cell_from_id` back to an
    equal geometry.
  - Edge-case overlap IDs found only via neighbour search: `6r32gxpn0w4`
    (Python) and `gcp0yqzxpk4` (Go) → `is_primary == false` where applicable.
- **ESA alignment** (`esa_alignment_tests.rs`): reproduce the ESA `linspace + mod`
  reference for latitudes and longitudes at `d ∈ {5,7,10,13,50,100} km` and assert
  the grid lines match within `1e-10`; assert equator and prime meridian fall on
  grid lines.
- **Migrate** (`migrate_tests.rs`): migrated cell contains the decoded old centroid
  for `["dr18zj1ntew","dr19n8zgg4e","s000003037z","6r32gxpn0w4"]`; result is
  primary; id length 11; `"short"` errors.
- **Cross-language conformance** (`cross_language_tests.rs`): for every case in
  `testdata/cross_language_reference.json`
  (`southampton_320_overlap`, `southampton_320_no_overlap`, `equator_5000_overlap`,
  `equator_10000_no_overlap`, `high_lat_320_overlap`): generated **cell count**
  matches, every reference cell **id** is present, and the corresponding cell's
  bbox matches the reference `coords` within `1e-8`. **This is the definition of
  done for parity.**

Tooling: `cargo test`, `cargo fmt --check`, `cargo clippy -- -D warnings`. Add a
GitHub Actions workflow mirroring `mtgrid`'s `smoke.yaml` (fmt + clippy + test on
stable).

## 7. Risks & mitigations
- **Geohash parity** (highest risk): retired up front by the Section 2 spike;
  fallback is a vendored geohash matching `geolib`/`pierrre` base32 + bit order.
- **Floating-point drift** vs Python/Go: use `f64` throughout and identical
  operation order; the `1e-8`/`1e-10` tolerances in the fixtures absorb the rest.
- **Intersection-predicate semantics** (boundary touches): `geo`'s `Intersects`
  must treat edge contact the way Shapely/orb-predicates do; the fixed cell counts
  (e.g. 225, 135) in the conformance tests will catch any discrepancy.
- **Antimeridian / pole handling**: covered by the `high_lat` and equator fixtures.

## 8. Out of scope (initial release)
- Row-level parallelism (add later behind a feature flag if benchmarks justify it).
- Python/Go-style packaging/publishing (crates.io release) — can follow once parity
  is green.
- `check_overlaps.py`-style auxiliary tooling.
