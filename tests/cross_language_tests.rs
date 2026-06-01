use std::collections::HashMap;

use geo::{BoundingRect, Coord, LineString, Polygon};
use majortom::MajorTomGrid;
use serde::Deserialize;

#[derive(Deserialize)]
struct CrossLangCell {
    id: String,
    coords: Vec<Vec<f64>>,
}

#[derive(Deserialize)]
struct Config {
    d: u64,
    overlap: bool,
}

#[derive(Deserialize)]
struct CrossLangCase {
    count: usize,
    cells: Vec<CrossLangCell>,
    config: Config,
    polygon: Vec<Vec<f64>>,
}

fn polygon_from_points(points: &[Vec<f64>]) -> Polygon<f64> {
    let ring: Vec<Coord<f64>> = points.iter().map(|p| Coord { x: p[0], y: p[1] }).collect();
    Polygon::new(LineString::new(ring), vec![])
}

#[test]
fn cross_language_conformance() {
    let data = include_str!("../testdata/cross_language_reference.json");
    let cases: HashMap<String, CrossLangCase> =
        serde_json::from_str(data).expect("parse reference JSON");

    assert!(!cases.is_empty(), "fixture should contain cases");

    for (name, case) in &cases {
        let polygon = polygon_from_points(&case.polygon);
        let grid = MajorTomGrid::new(case.config.d, case.config.overlap).unwrap();
        let cells = grid.generate_grid_cells(&polygon);

        assert_eq!(
            cells.len(),
            case.count,
            "[{name}] cell count mismatch: rust={} reference={}",
            cells.len(),
            case.count
        );

        let by_id: HashMap<&str, &_> = cells.iter().map(|c| (c.id(), c)).collect();

        for ref_cell in &case.cells {
            let cell = by_id.get(ref_cell.id.as_str()).unwrap_or_else(|| {
                panic!(
                    "[{name}] reference cell id {} not found in output",
                    ref_cell.id
                )
            });

            let bounds = cell.geom.bounding_rect().unwrap();
            let ref_min_lon = ref_cell.coords[0][0];
            let ref_min_lat = ref_cell.coords[0][1];
            let ref_max_lon = ref_cell.coords[2][0];
            let ref_max_lat = ref_cell.coords[2][1];

            let tol = 1e-8;
            assert!(
                (bounds.min().x - ref_min_lon).abs() < tol
                    && (bounds.min().y - ref_min_lat).abs() < tol
                    && (bounds.max().x - ref_max_lon).abs() < tol
                    && (bounds.max().y - ref_max_lat).abs() < tol,
                "[{name}] cell {} bbox differs: rust=[{},{}]-[{},{}] reference=[{},{}]-[{},{}]",
                ref_cell.id,
                bounds.min().x,
                bounds.min().y,
                bounds.max().x,
                bounds.max().y,
                ref_min_lon,
                ref_min_lat,
                ref_max_lon,
                ref_max_lat
            );
        }
    }
}
