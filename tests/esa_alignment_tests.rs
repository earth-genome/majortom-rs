use majortom::MajorTomGrid;
use std::f64::consts::PI;

const EARTH_RADIUS_KM: f64 = 6378.137;
const EARTH_RADIUS_M: f64 = 6378137.0;

/// ESA reference latitude grid lines via the `linspace + mod` approach.
fn esa_latitudes(dist_km: f64) -> Vec<f64> {
    let num_divisions = (PI * EARTH_RADIUS_KM / dist_km).ceil() as usize;
    let step = 180.0 / num_divisions as f64;
    let mut lats: Vec<f64> = (0..num_divisions)
        .map(|i| {
            let v = -90.0 + i as f64 * step;
            let mut v = v.rem_euclid(180.0);
            if v < 0.0 {
                v += 180.0;
            }
            v - 90.0
        })
        .collect();
    lats.sort_by(|a, b| a.partial_cmp(b).unwrap());
    lats
}

/// ESA reference longitude grid lines for a given latitude.
fn esa_longitudes(lat: f64, dist_km: f64) -> Vec<f64> {
    let circumference = 2.0 * PI * EARTH_RADIUS_KM * (lat * PI / 180.0).cos();
    let num_divisions = (circumference / dist_km).ceil() as usize;
    let step = 360.0 / num_divisions as f64;
    let mut lons: Vec<f64> = (0..num_divisions)
        .map(|i| {
            let v = -180.0 + i as f64 * step;
            let mut v = v.rem_euclid(360.0);
            if v < 0.0 {
                v += 360.0;
            }
            v - 180.0
        })
        .collect();
    lons.sort_by(|a, b| a.partial_cmp(b).unwrap());
    lons
}

fn eg_latitudes(grid: &MajorTomGrid) -> Vec<f64> {
    (0..grid.row_count()).map(|i| grid.row_lat(i)).collect()
}

fn eg_longitudes(grid: &MajorTomGrid, lat: f64, d: f64) -> Vec<f64> {
    let lon_spacing = grid.lon_spacing(lat);
    let lon_offset = grid.lon_offset(lon_spacing);
    let lat_rad = lat.clamp(-89.0, 89.0).to_radians();
    let n_cols = (2.0 * PI * EARTH_RADIUS_M * lat_rad.cos() / d).ceil() as i64;
    (0..n_cols)
        .map(|i| grid.col_lon(i, lon_spacing, lon_offset))
        .collect()
}

#[test]
fn latitude_alignment() {
    for dist_km in [5.0, 10.0, 50.0, 100.0] {
        let dist_m = (dist_km * 1000.0) as u64;
        let grid = MajorTomGrid::new(dist_m, false).unwrap();
        let esa = esa_latitudes(dist_km);
        let eg = eg_latitudes(&grid);
        assert_eq!(
            esa.len(),
            eg.len(),
            "dist={dist_km}km: latitude count mismatch"
        );
        for (i, (a, b)) in eg.iter().zip(esa.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-10,
                "dist={dist_km}km: latitude[{i}] differs: eg={a} esa={b}"
            );
        }
    }
}

#[test]
fn longitude_alignment() {
    for dist_km in [5.0, 10.0, 50.0, 100.0] {
        let dist_m = (dist_km * 1000.0) as u64;
        let grid = MajorTomGrid::new(dist_m, false).unwrap();
        for test_lat in [0.0, 30.0, 45.0, 60.0] {
            let esa = esa_longitudes(test_lat, dist_km);
            let eg = eg_longitudes(&grid, test_lat, dist_m as f64);
            assert_eq!(
                esa.len(),
                eg.len(),
                "dist={dist_km}km, lat={test_lat}: longitude count mismatch"
            );
            for (i, (a, b)) in eg.iter().zip(esa.iter()).enumerate() {
                assert!(
                    (a - b).abs() < 1e-10,
                    "dist={dist_km}km, lat={test_lat}: longitude[{i}] differs: eg={a} esa={b}"
                );
            }
        }
    }
}

#[test]
fn equator_on_grid_line() {
    for dist_km in [5.0, 7.0, 10.0, 13.0, 50.0, 100.0] {
        let dist_m = (dist_km * 1000.0) as u64;
        let grid = MajorTomGrid::new(dist_m, false).unwrap();
        let lats = eg_latitudes(&grid);
        assert!(
            lats.contains(&0.0),
            "dist={dist_km}km: equator (0.0) should be a grid line"
        );
    }
}

#[test]
fn prime_meridian_on_grid_line() {
    for dist_km in [5.0, 7.0, 10.0, 13.0, 50.0, 100.0] {
        let dist_m = (dist_km * 1000.0) as u64;
        let grid = MajorTomGrid::new(dist_m, false).unwrap();
        for test_lat in [0.0, 30.0, 60.0] {
            let lons = eg_longitudes(&grid, test_lat, dist_m as f64);
            assert!(
                lons.contains(&0.0),
                "dist={dist_km}km, lat={test_lat}: prime meridian (0.0) should be a grid line"
            );
        }
    }
}
