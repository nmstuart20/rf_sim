//! Local north-east-down frame at a geodetic reference point (contracts §12.3).
//!
//! The NED axes at geodetic latitude φ and longitude λ, in ECEF coordinates, are the
//! columns of `C_ecef_ned`:
//!
//! ```text
//! north = (−sin φ·cos λ, −sin φ·sin λ,  cos φ)
//! east  = (−sin λ,        cos λ,         0    )
//! down  = (−cos φ·cos λ, −cos φ·sin λ, −sin φ)
//! ```
//!
//! "Down" is along the inward ellipsoid normal (geodetic, not geocentric), so it is
//! orthogonal to the ellipsoid at the reference point. The NED origin is the
//! reference point's ECEF position, including its height.

use crate::geodetic::{Ecef, Geodetic};
use crate::quat::Quat;
use crate::vec::{Mat3, Vec3};

/// A local NED frame (module docs).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LocalNed {
    reference: Geodetic,
    origin: Ecef,
    /// `C_ecef_ned`: NED coordinates → ECEF coordinates.
    c_ecef_ned: Mat3,
}

impl LocalNed {
    /// The NED frame at `reference`.
    pub fn new(reference: Geodetic) -> Self {
        let (sp, cp) = libm::sincos(reference.lat_rad);
        let (sl, cl) = libm::sincos(reference.lon_rad);
        let c_ecef_ned = Mat3([
            [-(sp * cl), -sl, -(cp * cl)],
            [-(sp * sl), cl, -(cp * sl)],
            [cp, 0.0, -sp],
        ]);
        Self {
            reference,
            origin: reference.to_ecef(),
            c_ecef_ned,
        }
    }

    /// The reference point.
    pub const fn reference(&self) -> Geodetic {
        self.reference
    }

    /// The origin in ECEF.
    pub const fn origin(&self) -> Ecef {
        self.origin
    }

    /// `C_ecef_ned` (columns are the north, east and down axes in ECEF).
    pub const fn c_ecef_ned(&self) -> Mat3 {
        self.c_ecef_ned
    }

    /// `q_ecef_ned`, the quaternion of [`Self::c_ecef_ned`].
    pub fn q_ecef_ned(&self) -> Quat {
        Quat::from_mat3(&self.c_ecef_ned)
    }

    /// An ECEF direction (or displacement) in NED coordinates.
    pub fn vector_to_ned(&self, v_ecef: Vec3) -> Vec3 {
        self.c_ecef_ned.transpose().mul_vec(v_ecef)
    }

    /// A NED direction (or displacement) in ECEF coordinates.
    pub fn vector_to_ecef(&self, v_ned: Vec3) -> Vec3 {
        self.c_ecef_ned.mul_vec(v_ned)
    }

    /// An ECEF point relative to the origin, in NED.
    pub fn point_to_ned(&self, p: Ecef) -> Vec3 {
        self.vector_to_ned(p.vec() - self.origin.vec())
    }

    /// A NED point (relative to the origin) in ECEF.
    pub fn point_to_ecef(&self, p_ned: Vec3) -> Ecef {
        Ecef::from_vec(self.origin.vec() + self.vector_to_ecef(p_ned))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodetic::WGS84_A_M;

    fn close(a: Vec3, b: Vec3, tol: f64) -> bool {
        a.distance(b) <= tol
    }

    #[test]
    fn hand_computed_axes() {
        // At (0°, 0°): north = +z, east = +y, down = −x.
        let n = LocalNed::new(Geodetic::new(0.0, 0.0, 0.0));
        assert!(close(n.vector_to_ecef(Vec3::X), Vec3::Z, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Y), Vec3::Y, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Z), -Vec3::X, 1e-16));
        // At (0°, 90° E): north = +z, east = −x, down = −y.
        let n = LocalNed::new(Geodetic::from_degrees(0.0, 90.0, 0.0));
        assert!(close(n.vector_to_ecef(Vec3::X), Vec3::Z, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Y), -Vec3::X, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Z), -Vec3::Y, 1e-16));
        // At the north pole on the 0° meridian: north = −x, east = +y, down = −z.
        let n = LocalNed::new(Geodetic::from_degrees(90.0, 0.0, 0.0));
        assert!(close(n.vector_to_ecef(Vec3::X), -Vec3::X, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Y), Vec3::Y, 1e-16));
        assert!(close(n.vector_to_ecef(Vec3::Z), -Vec3::Z, 1e-16));
        // At 45° N, 0°: north = (−√½, 0, √½), down = (−√½, 0, −√½).
        let n = LocalNed::new(Geodetic::from_degrees(45.0, 0.0, 0.0));
        let h = core::f64::consts::FRAC_1_SQRT_2;
        assert!(close(
            n.vector_to_ecef(Vec3::X),
            Vec3::new(-h, 0.0, h),
            1e-15
        ));
        assert!(close(
            n.vector_to_ecef(Vec3::Z),
            Vec3::new(-h, 0.0, -h),
            1e-15
        ));
    }

    #[test]
    fn down_is_the_ellipsoid_normal_and_points_round_trip() {
        let g = Geodetic::from_degrees(38.9, -77.0, 120.0);
        let n = LocalNed::new(g);
        // A point 1 km "up" is the same geodetic position 1 km higher.
        let up = n.point_to_ecef(Vec3::new(0.0, 0.0, -1000.0));
        let want = Geodetic::new(g.lat_rad, g.lon_rad, g.h_m + 1000.0).to_ecef();
        assert!(up.distance(want) < 1e-8);
        // Orthonormal, right-handed.
        let c = n.c_ecef_ned();
        let p = c.mul(&c.transpose());
        for i in 0..3 {
            for j in 0..3 {
                let e = if i == j { 1.0 } else { 0.0 };
                assert!((p.0[i][j] - e).abs() < 1e-15);
            }
        }
        let (nn, ee, dd) = (
            n.vector_to_ecef(Vec3::X),
            n.vector_to_ecef(Vec3::Y),
            n.vector_to_ecef(Vec3::Z),
        );
        assert!(close(nn.cross(ee), dd, 1e-15));
        // Points round-trip, and the quaternion matches the matrix.
        let p = Ecef::new(WGS84_A_M * 0.6, -WGS84_A_M * 0.5, 4.0e6);
        assert!(n.point_to_ecef(n.point_to_ned(p)).distance(p) < 1e-8);
        let q = n.q_ecef_ned();
        for v in [Vec3::X, Vec3::Y, Vec3::new(3.0, -1.0, 2.0)] {
            assert!(close(q.rotate(v), n.vector_to_ecef(v), 1e-14));
        }
    }
}
