//! n-vector geometry on the WGS84 ellipsoid (contracts §13.3).
//!
//! A point is `(n, h)`: `n` the outward geodetic normal (unit up vector,
//! `(cos φ cos λ, cos φ sin λ, sin φ)`) and `h` the ellipsoid height. Unlike
//! latitude and longitude, `n` has no singularity at the poles.

use crate::GeoError;
use crate::geodetic::{Ecef, WGS84_A_M, WGS84_F};
use crate::vec::Vec3;

/// First eccentricity squared, e² = f(2 − f).
pub(crate) fn e2() -> f64 {
    WGS84_F * (2.0 - WGS84_F)
}

/// W = √(1 − e² sin² φ), with sin φ = n_z.
fn w_factor(n: Vec3) -> f64 {
    libm::sqrt(1.0 - e2() * (n.z * n.z))
}

/// Prime-vertical (N) and meridian (M) radii of curvature at normal `n`, metres.
pub(crate) fn radii(n: Vec3) -> (f64, f64) {
    let w = w_factor(n);
    let big_n = WGS84_A_M / w;
    let m = WGS84_A_M * (1.0 - e2()) / (w * w * w);
    (big_n, m)
}

/// The ECEF point at normal `n` and ellipsoid height `h`:
/// `(N n_x, N n_y, N(1 − e²) n_z) + h n`, the geodetic forward map.
pub(crate) fn point(n: Vec3, h: f64) -> Vec3 {
    let (big_n, _) = radii(n);
    Vec3::new(big_n * n.x, big_n * n.y, big_n * (1.0 - e2()) * n.z) + n * h
}

/// The unit normal and ellipsoid height of an ECEF point (geodetic inverse).
pub(crate) fn normal_of(p: Ecef) -> Result<(Vec3, f64), GeoError> {
    let g = p.to_geodetic()?;
    Ok((normal(g.lat_rad, g.lon_rad), g.h_m))
}

/// The unit normal at latitude and longitude (radians).
pub(crate) fn normal(lat_rad: f64, lon_rad: f64) -> Vec3 {
    let (sp, cp) = libm::sincos(lat_rad);
    let (sl, cl) = libm::sincos(lon_rad);
    Vec3::new(cp * cl, cp * sl, sp)
}

/// The rate of change of the normal for horizontal velocity `x` (⊥ n) at height h:
/// `K x = x/(N + h) + c (w·x) w`, with `w = ẑ − n_z n` (|w| = cos φ, along north) and
/// `c = a e² / (W³ (M + h)(N + h))`, which equals `(1/(M+h) − 1/(N+h))/cos² φ` and is
/// smooth at the poles. Northward motion turns n at `v/(M + h)`, eastward at
/// `v/(N + h)`.
pub(crate) fn curvature(n: Vec3, h: f64, x: Vec3) -> Vec3 {
    let w = w_factor(n);
    let (big_n, m) = radii(n);
    let c = WGS84_A_M * e2() / ((w * w * w) * ((m + h) * (big_n + h)));
    let wv = Vec3::new(-(n.z * n.x), -(n.z * n.y), 1.0 - n.z * n.z);
    x * (1.0 / (big_n + h)) + wv * (c * wv.dot(x))
}

/// Local north and east unit vectors at normal `n`, or `None` within `min_cos` of a
/// pole (cos φ below it).
pub(crate) fn north_east(n: Vec3, min_cos: f64) -> Option<(Vec3, Vec3)> {
    let east = Vec3::Z.cross(n);
    let cos_lat = east.norm();
    if cos_lat < min_cos {
        return None;
    }
    let east = east * (1.0 / cos_lat);
    Some((n.cross(east), east))
}

/// Smallest cos φ at which a track angle is defined (about 2·10⁻⁸ ° from a pole).
pub(crate) const TRACK_MIN_COS: f64 = 1.0e-9;

/// Track angle of horizontal direction `d` at normal `n`, radians clockwise from north
/// in (−π, π].
pub(crate) fn track_angle(n: Vec3, d: Vec3) -> Result<f64, GeoError> {
    let (north, east) = north_east(n, TRACK_MIN_COS)
        .ok_or_else(|| GeoError::Trajectory("a track angle is undefined at a pole".to_owned()))?;
    Ok(libm::atan2(d.dot(east), d.dot(north)))
}

/// The horizontal unit direction with track angle `track` at normal `n`.
pub(crate) fn direction(n: Vec3, track: f64) -> Result<Vec3, GeoError> {
    let (north, east) = north_east(n, TRACK_MIN_COS)
        .ok_or_else(|| GeoError::Trajectory("a heading is undefined at a pole".to_owned()))?;
    let (s, c) = libm::sincos(track);
    Ok(north * c + east * s)
}

/// `a` wrapped to (−π, π].
pub(crate) fn wrap_pi(a: f64) -> f64 {
    let two_pi = 2.0 * core::f64::consts::PI;
    let mut r = libm::remainder(a, two_pi);
    if r <= -core::f64::consts::PI {
        r += two_pi;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodetic::Geodetic;

    #[test]
    fn point_matches_the_geodetic_forward_map() {
        for (lat, lon, h) in [(0.0, 0.0, 0.0), (51.5, -0.1, 10_000.0), (-89.9, 170.0, 5.0)] {
            let g = Geodetic::from_degrees(lat, lon, h);
            let p = point(normal(g.lat_rad, g.lon_rad), h);
            assert!(p.distance(g.to_ecef().vec()) < 1e-8, "{lat} {lon}");
        }
    }

    #[test]
    fn curvature_rates_match_finite_differences() {
        // Move 10 m north, east and diagonally; the change in n over that step,
        // divided by the step, is K applied to the unit direction.
        for lat in [0.0_f64, 45.0, 89.0, 89.999_999] {
            let h = 9_000.0;
            let g = Geodetic::from_degrees(lat, 10.0, h);
            let n = normal(g.lat_rad, g.lon_rad);
            let (north, east) = north_east(n, 0.0).unwrap();
            for dir in [north, east, (north + east).normalized()] {
                let step = 10.0;
                let p = point(n, h) + dir * step;
                let (n2, _) = normal_of(Ecef::from_vec(p)).unwrap();
                let fd = (n2 - n) * (1.0 / step);
                let k = curvature(n, h, dir);
                assert!(
                    (fd - k).norm() < 1e-6 * k.norm(),
                    "lat {lat}: {fd:?} vs {k:?}"
                );
            }
        }
    }

    #[test]
    fn track_and_direction_invert() {
        let n = normal(0.7, -2.0);
        for t in [-3.0, -1.0, 0.0, 0.5, 3.1] {
            let d = direction(n, t).unwrap();
            assert!((track_angle(n, d).unwrap() - t).abs() < 1e-15);
        }
        assert!(track_angle(Vec3::Z, Vec3::X).is_err());
        assert_eq!(wrap_pi(3.0 * core::f64::consts::PI), core::f64::consts::PI);
        assert!((wrap_pi(-0.5) + 0.5).abs() < 1e-16);
    }
}
