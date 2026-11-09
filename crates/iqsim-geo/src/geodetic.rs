//! WGS84 geodetic coordinates and ECEF positions (contracts §12.2).
//!
//! **Forward** (geodetic → ECEF) is the closed form
//! `N = a/√(1 − e²·sin²φ)`,
//! `(x, y, z) = ((N + h)·cos φ·cos λ, (N + h)·cos φ·sin λ, (N·(1 − e²) + h)·sin φ)`.
//!
//! **Inverse** (ECEF → geodetic) is Bowring's (1976) iteration on the parametric
//! latitude β:
//!
//! 1. `p = √(x² + y²)`, `β₀ = atan2(z, (1 − f)·p)`.
//! 2. `φ = atan2(z + e′²·b·sin³β, p − e²·a·cos³β)`, then
//!    `β = atan2((1 − f)·sin φ, cos φ)`.
//! 3. Repeat 2 until `|φₖ − φₖ₋₁| ≤ 10⁻¹⁵ rad` (about 4 ulp of φ near ±π/2), at most
//!    [`MAX_INVERSE_ITERATIONS`] times; not converging is an error.
//! 4. `h = p·cos φ + z·sin φ − a·√(1 − e²·sin²φ)`, which is well conditioned at every
//!    latitude, including the poles.
//!
//! The error after k steps falls roughly as the cube of the previous one, so a point
//! near the surface converges after 2 steps and a lunar-distance point after 2–3.
//! The domain is [`MIN_INVERSE_RADIUS_M`] (1000 km) `≤ |r| ≤`
//! [`MAX_INVERSE_RADIUS_M`] (10¹² m): far outside the ellipsoid's evolute (radius
//! ≈ 43 km), where geodetic latitude is not unique, and far below f64 overflow.
//! On the polar axis (`p = 0`) the latitude is ±90°, the longitude is 0 by convention
//! and `h = |z| − b`. Longitude is `atan2(y, x)` in (−180°, 180°].
//!
//! Every transcendental goes through `libm`, and every product and sum is a separate
//! operation, so the bits are the same on every machine.

use crate::GeoError;
use crate::vec::Vec3;

/// Speed of light in vacuum, m/s (exact, SI).
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// WGS84 semi-major axis a, metres.
pub const WGS84_A_M: f64 = 6_378_137.0;

/// WGS84 flattening f.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// Smallest distance from the Earth's centre the inverse accepts, metres.
pub const MIN_INVERSE_RADIUS_M: f64 = 1.0e6;

/// Largest distance from the Earth's centre the inverse accepts, metres (10¹² m,
/// about 7 AU). Beyond it `hypot` and the height could overflow; nothing in the
/// simulator's scope is that far.
pub const MAX_INVERSE_RADIUS_M: f64 = 1.0e12;

/// Iteration cap of the inverse (it never needs more than 4 inside the domain).
pub const MAX_INVERSE_ITERATIONS: u32 = 8;

/// Inverse stopping threshold on the latitude update, radians.
const INVERSE_TOL_RAD: f64 = 1.0e-15;

const DEG_TO_RAD: f64 = core::f64::consts::PI / 180.0;

/// WGS84 first eccentricity squared, `e² = f·(2 − f)`.
fn e2() -> f64 {
    WGS84_F * (2.0 - WGS84_F)
}

/// WGS84 semi-minor axis `b = a·(1 − f)`, metres.
pub fn wgs84_b_m() -> f64 {
    WGS84_A_M * (1.0 - WGS84_F)
}

/// A WGS84 Earth-centred, Earth-fixed position, metres.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ecef {
    /// x, metres.
    pub x: f64,
    /// y, metres.
    pub y: f64,
    /// z, metres.
    pub z: f64,
}

impl Ecef {
    /// The position (x, y, z) metres.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The position from a vector.
    pub const fn from_vec(v: Vec3) -> Self {
        Self::new(v.x, v.y, v.z)
    }

    /// The position as a vector.
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    /// The ECEF position of WGS84 geodetic latitude and longitude (degrees) and height
    /// above the ellipsoid (metres). Same bits as
    /// `Geodetic::from_degrees(lat, lon, h).to_ecef()`.
    pub fn from_geodetic(lat_deg: f64, lon_deg: f64, h_m: f64) -> Self {
        Geodetic::from_degrees(lat_deg, lon_deg, h_m).to_ecef()
    }

    /// The geodetic coordinates of this position (module docs: Bowring's iteration).
    ///
    /// Errors: a non-finite coordinate, or a point closer than
    /// [`MIN_INVERSE_RADIUS_M`] to the Earth's centre.
    pub fn to_geodetic(self) -> Result<Geodetic, GeoError> {
        to_geodetic_counted(self).map(|(g, _)| g)
    }

    /// [`Self::to_geodetic`] and the number of Bowring steps it took (0 on the polar
    /// axis), for convergence diagnostics and tests.
    pub fn to_geodetic_with_steps(self) -> Result<(Geodetic, u32), GeoError> {
        to_geodetic_counted(self)
    }

    /// Straight-line distance to `other`, metres (IEEE sqrt, fixed operation order).
    pub fn distance(self, other: Ecef) -> f64 {
        let dx = other.x - self.x;
        let dy = other.y - self.y;
        let dz = other.z - self.z;
        let s = dx * dx + dy * dy;
        let s = s + dz * dz;
        s.sqrt()
    }
}

/// WGS84 geodetic coordinates: latitude and longitude in radians, height above the
/// ellipsoid in metres.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Geodetic {
    /// Geodetic latitude φ, radians, in [−π/2, π/2].
    pub lat_rad: f64,
    /// Longitude λ, radians (east positive).
    pub lon_rad: f64,
    /// Height above the WGS84 ellipsoid, metres.
    pub h_m: f64,
}

impl Geodetic {
    /// From radians and metres.
    pub const fn new(lat_rad: f64, lon_rad: f64, h_m: f64) -> Self {
        Self {
            lat_rad,
            lon_rad,
            h_m,
        }
    }

    /// From degrees and metres (`deg · (π/180)`, one rounding each).
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, h_m: f64) -> Self {
        Self::new(lat_deg * DEG_TO_RAD, lon_deg * DEG_TO_RAD, h_m)
    }

    /// Latitude in degrees.
    pub fn lat_deg(self) -> f64 {
        self.lat_rad / DEG_TO_RAD
    }

    /// Longitude in degrees.
    pub fn lon_deg(self) -> f64 {
        self.lon_rad / DEG_TO_RAD
    }

    /// The ECEF position (module docs: closed form).
    pub fn to_ecef(self) -> Ecef {
        let (sp, cp) = libm::sincos(self.lat_rad);
        let (sl, cl) = libm::sincos(self.lon_rad);
        let e2 = e2();
        let n = WGS84_A_M / (1.0 - e2 * (sp * sp)).sqrt();
        let r = (n + self.h_m) * cp;
        Ecef {
            x: r * cl,
            y: r * sl,
            z: (n * (1.0 - e2) + self.h_m) * sp,
        }
    }
}

/// The inverse and the number of Bowring steps it took (0 on the polar axis).
fn to_geodetic_counted(e: Ecef) -> Result<(Geodetic, u32), GeoError> {
    if !(e.x.is_finite() && e.y.is_finite() && e.z.is_finite()) {
        return Err(GeoError::NonFinite("ECEF position"));
    }
    let p = libm::hypot(e.x, e.y);
    let r = libm::hypot(p, e.z);
    if r.is_nan() || r > MAX_INVERSE_RADIUS_M {
        return Err(GeoError::OutsideDomain { radius_m: r });
    }
    if r < MIN_INVERSE_RADIUS_M {
        return Err(GeoError::InsideDomain { radius_m: r });
    }
    let a = WGS84_A_M;
    let b = wgs84_b_m();
    let e2 = e2();
    let ep2 = e2 / (1.0 - e2);
    let one_minus_f = 1.0 - WGS84_F;

    if p == 0.0 {
        let lat = if e.z > 0.0 {
            core::f64::consts::FRAC_PI_2
        } else {
            -core::f64::consts::FRAC_PI_2
        };
        return Ok((Geodetic::new(lat, 0.0, e.z.abs() - b), 0));
    }

    let mut beta = libm::atan2(e.z, one_minus_f * p);
    let mut lat = f64::NAN;
    let mut steps = 0;
    loop {
        if steps == MAX_INVERSE_ITERATIONS {
            return Err(GeoError::NoConvergence {
                what: "geodetic inverse",
                iterations: steps,
            });
        }
        steps += 1;
        let (sb, cb) = libm::sincos(beta);
        let sb3 = sb * sb * sb;
        let cb3 = cb * cb * cb;
        let num = e.z + ep2 * b * sb3;
        let den = p - e2 * a * cb3;
        let next = libm::atan2(num, den);
        let done = (next - lat).abs() <= INVERSE_TOL_RAD;
        lat = next;
        if done {
            break;
        }
        let (sp, cp) = libm::sincos(lat);
        beta = libm::atan2(one_minus_f * sp, cp);
    }

    let (sp, cp) = libm::sincos(lat);
    let w = (1.0 - e2 * (sp * sp)).sqrt();
    let h = p * cp + e.z * sp;
    let h = h - a * w;
    let mut lon = libm::atan2(e.y, e.x);
    if lon == -core::f64::consts::PI {
        lon = core::f64::consts::PI;
    }
    Ok((Geodetic::new(lat, lon, h), steps))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geodetic_matches_reference_points() {
        // Exact points: equator at 0° and 90° east, and the pole.
        let e = Ecef::from_geodetic(0.0, 0.0, 0.0);
        assert!((e.x - WGS84_A_M).abs() < 1e-9 && e.y.abs() < 1e-9 && e.z.abs() < 1e-9);
        let e = Ecef::from_geodetic(0.0, 90.0, 5.0);
        assert!(e.x.abs() < 1e-9 && (e.y - WGS84_A_M - 5.0).abs() < 1e-9);
        let b = wgs84_b_m();
        let p = Ecef::from_geodetic(90.0, 0.0, 100.0);
        assert!(p.x.abs() < 1e-6 && (p.z - (b + 100.0)).abs() < 1e-6);
        // Elsewhere: h = 0 lies on the ellipsoid, and height moves along the geodetic
        // normal (cos φ cos λ, cos φ sin λ, sin φ).
        for (lat, lon) in [(38.8898, -77.0365), (-33.8688, 151.2093), (71.0, 25.0)] {
            let s0 = Ecef::from_geodetic(lat, lon, 0.0);
            let on = (s0.x * s0.x + s0.y * s0.y) / (WGS84_A_M * WGS84_A_M) + s0.z * s0.z / (b * b);
            assert!((on - 1.0).abs() < 1e-14, "{on}");
            let s1 = Ecef::from_geodetic(lat, lon, 1000.0);
            let (sp, cp) = libm::sincos(lat.to_radians());
            let (sl, cl) = libm::sincos(lon.to_radians());
            let normal = Ecef::new(cp * cl, cp * sl, sp);
            let d = Ecef::new(s1.x - s0.x, s1.y - s0.y, s1.z - s0.z);
            assert!(
                d.distance(Ecef::new(
                    1000.0 * normal.x,
                    1000.0 * normal.y,
                    1000.0 * normal.z
                )) < 1e-6
            );
        }
    }

    #[test]
    fn distance_is_euclidean() {
        let a = Ecef::new(1.0, 2.0, 3.0);
        assert_eq!(a.distance(Ecef::new(4.0, 6.0, 3.0)), 5.0);
        assert_eq!(a.distance(a), 0.0);
    }

    #[test]
    fn inverse_axis_and_domain() {
        let b = wgs84_b_m();
        let n = Ecef::new(0.0, 0.0, b + 10.0).to_geodetic().unwrap();
        assert_eq!(n.lat_rad, core::f64::consts::FRAC_PI_2);
        assert_eq!(n.lon_rad, 0.0);
        assert!((n.h_m - 10.0).abs() < 1e-9);
        let s = Ecef::new(0.0, 0.0, -(b - 5.0)).to_geodetic().unwrap();
        assert_eq!(s.lat_rad, -core::f64::consts::FRAC_PI_2);
        assert!((s.h_m + 5.0).abs() < 1e-9);
        assert!(matches!(
            Ecef::new(1.0, 2.0, 3.0).to_geodetic(),
            Err(GeoError::InsideDomain { .. })
        ));
        // Review regression: a finite but huge point overflowed hypot to an infinite
        // height that was returned as Ok.
        assert!(matches!(
            Ecef::new(1.7e308, 1.7e308, 0.0).to_geodetic(),
            Err(GeoError::OutsideDomain { .. })
        ));
        assert!(Ecef::new(0.0, 0.0, 9.9e11).to_geodetic().is_ok());
        assert!(matches!(
            Ecef::new(f64::NAN, 0.0, 7e6).to_geodetic(),
            Err(GeoError::NonFinite(_))
        ));
        // −180° maps to +180°.
        let w = Ecef::new(-7e6, -0.0, 0.0).to_geodetic().unwrap();
        assert_eq!(w.lon_rad, core::f64::consts::PI);
    }
}
