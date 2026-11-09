//! Reference frames, trajectories, and propagation geometry.
//! Per the plan sections "World model" and "Antennas and propagation", this crate
//! owns WGS84/ECEF/ECI/NED transforms, C2 motion, attitude, geoid altitude,
//! and light-time geometry that retains fractional acquisition time.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]

/// Speed of light in vacuum, m/s (exact, SI).
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// WGS84 semi-major axis a, metres.
pub const WGS84_A_M: f64 = 6_378_137.0;

/// WGS84 flattening f.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;

/// A WGS84 Earth-centred, Earth-fixed position, metres.
///
/// Prompt 5 needs only fixed positions and straight-line distance. Frames, light time
/// and motion arrive in phase 3.
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

    /// The ECEF position of WGS84 geodetic latitude and longitude (degrees) and height
    /// above the ellipsoid (metres): `N = a/√(1 − e²·sin²φ)`,
    /// `(x, y, z) = ((N + h)·cos φ·cos λ, (N + h)·cos φ·sin λ, (N·(1 − e²) + h)·sin φ)`.
    /// Trigonometry goes through `libm` so the bits match on every platform.
    pub fn from_geodetic(lat_deg: f64, lon_deg: f64, h_m: f64) -> Self {
        let d2r = core::f64::consts::PI / 180.0;
        let (sp, cp) = libm::sincos(lat_deg * d2r);
        let (sl, cl) = libm::sincos(lon_deg * d2r);
        let e2 = WGS84_F * (2.0 - WGS84_F);
        let n = WGS84_A_M / (1.0 - e2 * (sp * sp)).sqrt();
        let r = (n + h_m) * cp;
        Self {
            x: r * cl,
            y: r * sl,
            z: (n * (1.0 - e2) + h_m) * sp,
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_is_euclidean() {
        let a = Ecef::new(1.0, 2.0, 3.0);
        assert_eq!(a.distance(Ecef::new(4.0, 6.0, 3.0)), 5.0);
        assert_eq!(a.distance(a), 0.0);
    }

    #[test]
    fn geodetic_matches_reference_points() {
        // Exact points: equator at 0° and 90° east, and the pole.
        let e = Ecef::from_geodetic(0.0, 0.0, 0.0);
        assert!((e.x - WGS84_A_M).abs() < 1e-9 && e.y.abs() < 1e-9 && e.z.abs() < 1e-9);
        let e = Ecef::from_geodetic(0.0, 90.0, 5.0);
        assert!(e.x.abs() < 1e-9 && (e.y - WGS84_A_M - 5.0).abs() < 1e-9);
        let b = WGS84_A_M * (1.0 - WGS84_F);
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
}
