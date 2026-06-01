//! A minimal geohash encoder/decoder using the classic interval-bisection
//! algorithm.
//!
//! This is vendored (rather than using the `geohash` crate) to guarantee
//! byte-for-byte parity with the Python (`geolib`) and Go (`pierrre/geohash`)
//! reference implementations. In particular, those libraries bisect the
//! latitude/longitude intervals and compare the coordinate against each
//! midpoint with `>=`. Because every midpoint is an exact dyadic multiple of
//! 90/180, the comparison is exact in `f64`, so a coordinate that is
//! infinitesimally below a cell boundary (e.g. a centroid of `-5e-15` at the
//! equator) is correctly placed in the lower cell. The `geohash` crate instead
//! uses a floating-point bit trick (`value / 180 + 1.5`) that rounds such
//! sub-ULP values to the wrong side of the boundary, which breaks parity for
//! cells straddling the equator or the prime meridian.

use crate::error::GridError;

/// Standard geohash base32 alphabet (note: excludes `a`, `i`, `l`, `o`).
const BASE32: &[u8; 32] = b"0123456789bcdefghjkmnpqrstuvwxyz";

/// Reverse lookup table for [`BASE32`]; `0xff` marks invalid characters.
fn base32_index(c: u8) -> Option<u8> {
    BASE32.iter().position(|&b| b == c).map(|p| p as u8)
}

/// Encodes a latitude/longitude pair into a geohash of the given precision.
///
/// Uses the same interval-bisection scheme (with a `>=` midpoint comparison) as
/// `geolib` and `pierrre/geohash`.
pub(crate) fn encode(lat: f64, lon: f64, precision: usize) -> String {
    let mut lat_min = -90.0f64;
    let mut lat_max = 90.0f64;
    let mut lon_min = -180.0f64;
    let mut lon_max = 180.0f64;

    let mut even = true;
    let mut bit = 0u8;
    let mut index = 0usize;
    let mut out = String::with_capacity(precision);

    while out.len() < precision {
        if even {
            let mid = (lon_min + lon_max) / 2.0;
            if lon >= mid {
                index = index * 2 + 1;
                lon_min = mid;
            } else {
                index *= 2;
                lon_max = mid;
            }
        } else {
            let mid = (lat_min + lat_max) / 2.0;
            if lat >= mid {
                index = index * 2 + 1;
                lat_min = mid;
            } else {
                index *= 2;
                lat_max = mid;
            }
        }
        even = !even;
        bit += 1;
        if bit == 5 {
            out.push(BASE32[index] as char);
            bit = 0;
            index = 0;
        }
    }

    out
}

/// Decodes a geohash to the centre `(lat, lon)` of its cell.
///
/// Mirrors the bisection used by the reference libraries and returns the
/// midpoint of the resulting latitude/longitude interval.
pub(crate) fn decode(hash: &str) -> Result<(f64, f64), GridError> {
    if hash.is_empty() {
        return Err(GridError::InvalidCellId);
    }

    let mut lat_min = -90.0f64;
    let mut lat_max = 90.0f64;
    let mut lon_min = -180.0f64;
    let mut lon_max = 180.0f64;
    let mut even = true;

    for c in hash.bytes() {
        let index = base32_index(c).ok_or(GridError::InvalidCellId)?;
        for n in (0..5).rev() {
            let bit = (index >> n) & 1;
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

    Ok(((lat_min + lat_max) / 2.0, (lon_min + lon_max) / 2.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_known_points() {
        assert_eq!(
            encode(39.54980995879779, -76.33601540899066, 11),
            "dr18yzvf3fv"
        );
        assert_eq!(
            encode(39.55555893480778, -76.34979642927647, 11),
            "dr19n8ggbf9"
        );
    }

    #[test]
    fn encodes_boundary_straddlers() {
        // Centroid just below the equator must land in the lower cell.
        assert_eq!(
            encode(-4.850286838831153e-15, -0.17964071856288427, 11),
            "7zzzgzvpypf"
        );
        // Centroid just west of the prime meridian.
        assert_eq!(
            encode(-0.49401197604790553, -5.053249485520439e-15, 11),
            "7zzvpypfpvr"
        );
    }

    #[test]
    fn round_trips() {
        for id in ["dr18yzvf3fv", "6r32gxpn0w4", "gcp0yqzxpk4", "s000003037z"] {
            let (lat, lon) = decode(id).unwrap();
            assert_eq!(encode(lat, lon, 11), id);
        }
    }

    #[test]
    fn rejects_invalid_characters() {
        assert!(decode("dr18yzvf3fa").is_err()); // 'a' is not in the alphabet
        assert!(decode("").is_err());
    }
}
