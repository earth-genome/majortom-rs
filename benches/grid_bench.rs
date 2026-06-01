use criterion::{criterion_group, criterion_main, Criterion};
use geo::{Coord, LineString, Polygon};
use majortom::MajorTomGrid;
use std::hint::black_box;

fn big_southampton() -> Polygon<f64> {
    Polygon::new(
        LineString::new(vec![
            Coord {
                x: -76.35673421721803,
                y: 39.55614384974018,
            },
            Coord {
                x: -76.35673421721803,
                y: 39.53123810591927,
            },
            Coord {
                x: -76.3131967920373,
                y: 39.53123810591927,
            },
            Coord {
                x: -76.3131967920373,
                y: 39.55614384974018,
            },
            Coord {
                x: -76.35673421721803,
                y: 39.55614384974018,
            },
        ]),
        vec![],
    )
}

fn bench_generate(c: &mut Criterion) {
    let poly = big_southampton();
    let grid = MajorTomGrid::new(320, true).unwrap();
    c.bench_function("generate_grid_cells_overlap", |b| {
        b.iter(|| grid.generate_grid_cells(black_box(&poly)))
    });

    let grid_no = MajorTomGrid::new(320, false).unwrap();
    c.bench_function("generate_grid_cells_no_overlap", |b| {
        b.iter(|| grid_no.generate_grid_cells(black_box(&poly)))
    });
}

fn bench_cell_from_id(c: &mut Criterion) {
    let grid = MajorTomGrid::new(320, true).unwrap();
    c.bench_function("cell_from_id", |b| {
        b.iter(|| grid.cell_from_id(black_box("dr19n8f7v6e")).unwrap())
    });
}

criterion_group!(benches, bench_generate, bench_cell_from_id);
criterion_main!(benches);
