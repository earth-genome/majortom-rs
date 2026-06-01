use geo::{BoundingRect, Coord, LineString, Polygon};
use majortom::MajorTomGrid;

fn polygon(coords: &[(f64, f64)]) -> Polygon<f64> {
    let ring: Vec<Coord<f64>> = coords.iter().map(|&(x, y)| Coord { x, y }).collect();
    Polygon::new(LineString::new(ring), vec![])
}

/// The `bigSouthampton` AOI used by the Python and Go test suites.
fn big_southampton() -> Polygon<f64> {
    polygon(&[
        (-76.35673421721803, 39.55614384974018),
        (-76.35673421721803, 39.53123810591927),
        (-76.3131967920373, 39.53123810591927),
        (-76.3131967920373, 39.55614384974018),
        (-76.35673421721803, 39.55614384974018),
    ])
}

#[test]
fn generate_grid_cells_basic_count() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let aoi = big_southampton();
    let cells = grid.generate_grid_cells(&aoi);
    assert_eq!(
        cells.len(),
        225,
        "bigSouthampton d=320 overlap should yield 225 cells"
    );
}

#[test]
fn all_generated_cells_intersect_aoi() {
    use geo::Intersects;
    let grid = MajorTomGrid::new(320, true).unwrap();
    let aoi = big_southampton();
    for cell in grid.generate_grid_cells(&aoi) {
        assert!(
            cell.geom.intersects(&aoi),
            "every returned cell must intersect the AOI"
        );
    }
}

#[test]
fn overlap_off_yields_fewer_cells() {
    let aoi = big_southampton();
    let with_overlap = MajorTomGrid::new(320, true)
        .unwrap()
        .generate_grid_cells(&aoi);
    let no_overlap = MajorTomGrid::new(320, false)
        .unwrap()
        .generate_grid_cells(&aoi);
    assert!(
        no_overlap.len() < with_overlap.len(),
        "overlap=false ({}) should yield fewer cells than overlap=true ({})",
        no_overlap.len(),
        with_overlap.len()
    );
}

#[test]
fn tiny_polygon_near_origin_yields_cells() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let aoi = polygon(&[
        (-0.0001, -0.0001),
        (-0.0001, 0.0001),
        (0.0001, 0.0001),
        (0.0001, -0.0001),
        (-0.0001, -0.0001),
    ]);
    let cells = grid.generate_grid_cells(&aoi);
    assert!(
        !cells.is_empty(),
        "a tiny polygon at the equator should yield at least one cell"
    );
}

#[test]
fn high_latitude_polygon_yields_cells() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let aoi = polygon(&[
        (170.0, 80.0),
        (171.0, 80.0),
        (171.0, 81.0),
        (170.0, 81.0),
        (170.0, 80.0),
    ]);
    let cells = grid.generate_grid_cells(&aoi);
    assert!(
        !cells.is_empty(),
        "high-latitude polygon should yield cells"
    );
}

#[test]
fn round_trip_cell_from_id() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let aoi = big_southampton();
    for cell in grid.generate_grid_cells(&aoi) {
        let found = grid
            .cell_from_id(cell.id())
            .expect("cell id should resolve");
        assert_eq!(found.id(), cell.id());
        let a = found.geom.bounding_rect().unwrap();
        let b = cell.geom.bounding_rect().unwrap();
        let tol = 1e-9;
        assert!(
            (a.min().x - b.min().x).abs() < tol
                && (a.min().y - b.min().y).abs() < tol
                && (a.max().x - b.max().x).abs() < tol
                && (a.max().y - b.max().y).abs() < tol,
            "resolved cell geometry should match the original for id {}",
            cell.id()
        );
    }
}

#[test]
fn edge_case_overlap_cell_found_via_neighbor_search() {
    // From the Python suite: this overlap cell is only found via the ±1
    // neighbour search and must be reported as non-primary.
    let grid = MajorTomGrid::new(320, true).unwrap();
    let cell = grid
        .cell_from_id("6r32gxpn0w4")
        .expect("edge-case id should resolve");
    assert_eq!(cell.id(), "6r32gxpn0w4");
    assert!(
        !cell.is_primary,
        "6r32gxpn0w4 should be an overlap (non-primary) cell"
    );
}

#[test]
fn odd_tile_id_and_truncation() {
    // From the Go suite: a tricky cell, plus the truncate-to-11 behaviour.
    let grid = MajorTomGrid::new(320, true).unwrap();
    let cell = grid
        .cell_from_id("gcp0yqzxpk4")
        .expect("gcp0yqzxpk4 should resolve");
    assert_eq!(cell.id(), "gcp0yqzxpk4");

    let truncated = grid
        .cell_from_id("gcp0yqzxpk4t24vzxu52")
        .expect("over-long id should be truncated and resolve");
    assert_eq!(truncated.id(), "gcp0yqzxpk4");
}

#[test]
fn new_rejects_zero_spacing() {
    assert!(MajorTomGrid::new(0, true).is_err());
}

#[test]
fn ids_are_unique_within_an_aoi() {
    use std::collections::HashSet;
    let grid = MajorTomGrid::new(320, true).unwrap();
    let cells = grid.generate_grid_cells(&big_southampton());
    let mut seen = HashSet::new();
    for cell in &cells {
        assert!(
            seen.insert(cell.id().to_string()),
            "duplicate id {}",
            cell.id()
        );
    }
}
