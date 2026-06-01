use geo::BoundingRect;
use majortom::MajorTomGrid;

/// Independent classic-bisection geohash decode, returning the cell centre
/// `(lat, lon)`. Mirrors `geolib`/`pierrre` and lets us verify migration
/// containment without relying on the crate's internals.
fn decode_center(hash: &str) -> (f64, f64) {
    let (mut lat_min, mut lat_max) = (-90.0f64, 90.0f64);
    let (mut lon_min, mut lon_max) = (-180.0f64, 180.0f64);
    let mut even = true;
    let base32 = b"0123456789bcdefghjkmnpqrstuvwxyz";
    for c in hash.bytes() {
        let idx = base32.iter().position(|&b| b == c).expect("valid base32") as u8;
        for n in (0..5).rev() {
            let bit = (idx >> n) & 1;
            if even {
                let mid = (lon_min + lon_max) / 2.0;
                if bit == 1 {
                    lon_min = mid;
                } else {
                    lon_max = mid;
                }
            } else {
                let mid = (lat_min + lat_max) / 2.0;
                if bit == 1 {
                    lat_min = mid;
                } else {
                    lat_max = mid;
                }
            }
            even = !even;
        }
    }
    ((lat_min + lat_max) / 2.0, (lon_min + lon_max) / 2.0)
}

#[test]
fn migrated_cell_contains_old_centroid() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    for old_id in ["dr18zj1ntew", "dr19n8zgg4e", "s000003037z", "6r32gxpn0w4"] {
        let cell = grid
            .migrate_cell_id(old_id)
            .expect("migrate should succeed");
        let (lat, lon) = decode_center(old_id);
        let b = cell.geom.bounding_rect().unwrap();
        assert!(
            lon >= b.min().x && lon <= b.max().x && lat >= b.min().y && lat <= b.max().y,
            "migrated cell for {old_id} must contain decoded centroid ({lat:.6}, {lon:.6})"
        );
    }
}

#[test]
fn migrated_cell_is_primary_with_valid_id() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let cell = grid.migrate_cell_id("dr18zj1ntew").unwrap();
    assert!(cell.is_primary, "migrated cell should be primary");
    assert_eq!(cell.id().len(), 11, "migrated cell id should be 11 chars");
}

#[test]
fn migrate_rejects_short_id() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    assert!(grid.migrate_cell_id("short").is_err());
}

#[test]
fn migrate_truncates_long_id() {
    let grid = MajorTomGrid::new(320, true).unwrap();
    let a = grid.migrate_cell_id("dr18zj1ntew").unwrap();
    let b = grid.migrate_cell_id("dr18zj1ntewt24vzxu52").unwrap();
    assert_eq!(
        a.id(),
        b.id(),
        "over-long ids should be truncated to 11 chars"
    );
}
