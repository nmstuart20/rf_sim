//! Rotation quaternions (contracts §12.4).
//!
//! **Convention.** Hamilton quaternions (`i² = j² = k² = ijk = −1`), stored scalar
//! first as `(w, x, y, z)`. A quaternion named `q_a_b` takes the coordinates of a
//! vector in frame *b* to its coordinates in frame *a*:
//!
//! ```text
//! v_a = q_a_b ⊗ v_b ⊗ q_a_b*        (rotate)
//! q_a_c = q_a_b ⊗ q_b_c             (compose: apply q_b_c first)
//! ```
//!
//! The same quaternion has two readings, and both are correct:
//!
//! - **active**: it rotates a vector, in a fixed frame, by the rotation that turns
//!   frame *a*'s axes onto frame *b*'s axes;
//! - **passive**: it re-expresses a fixed vector from frame *b* coordinates in frame
//!   *a* coordinates.
//!
//! For example the attitude `q_ned_body` maps body-axis coordinates (x forward,
//! y right, z down) to local NED coordinates. [`Quat::to_mat3`] gives the matching
//! direction-cosine matrix `C_a_b` with `v_a = C_a_b · v_b`.
//!
//! [`Quat::from_euler_zyx`] builds the aerospace 3-2-1 sequence
//! `q = q_z(yaw) ⊗ q_y(pitch) ⊗ q_x(roll)`: positive yaw turns the nose right
//! (east of north), positive pitch raises the nose, positive roll lowers the right
//! wing. Every operation is a separate multiply or add in a fixed order.

use core::ops::Mul;

use crate::vec::{Mat3, Vec3};

/// `cos(pitch)` below which [`Quat::to_euler_zyx`] treats the attitude as gimbal
/// locked. Off the lock, yaw and roll carry errors of about 10⁻¹⁶/cos(pitch); in the
/// locked branch the dropped split costs about cos(pitch). 10⁻⁸ balances the two, so
/// the rebuilt rotation is within ~2·10⁻⁸ everywhere.
const GIMBAL_LOCK_COS: f64 = 1e-8;

/// A Hamilton quaternion `w + x·i + y·j + z·k` (module docs).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quat {
    /// Scalar part.
    pub w: f64,
    /// i component.
    pub x: f64,
    /// j component.
    pub y: f64,
    /// k component.
    pub z: f64,
}

impl Default for Quat {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Quat {
    /// The identity rotation.
    pub const IDENTITY: Self = Self::new(1.0, 0.0, 0.0, 0.0);

    /// `(w, x, y, z)`, taken as given (not normalized).
    pub const fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }

    /// The vector part.
    pub const fn vector(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }

    /// Rotation by `angle_rad` about `axis` (right-handed; the axis is normalized
    /// here): `(cos(θ/2), sin(θ/2)·axis)`.
    pub fn from_axis_angle(axis: Vec3, angle_rad: f64) -> Self {
        let a = axis.normalized();
        let (s, c) = libm::sincos(angle_rad * 0.5);
        Self::new(c, a.x * s, a.y * s, a.z * s)
    }

    /// Aerospace 3-2-1 Euler angles (radians): `q_z(yaw) ⊗ q_y(pitch) ⊗ q_x(roll)`.
    pub fn from_euler_zyx(yaw_rad: f64, pitch_rad: f64, roll_rad: f64) -> Self {
        let qz = Self::from_axis_angle(Vec3::Z, yaw_rad);
        let qy = Self::from_axis_angle(Vec3::Y, pitch_rad);
        let qx = Self::from_axis_angle(Vec3::X, roll_rad);
        qz * qy * qx
    }

    /// The 3-2-1 Euler angles `(yaw, pitch, roll)` of a unit quaternion, radians:
    /// yaw and roll in (−π, π], pitch in [−π/2, π/2]. Rebuilding with
    /// [`Self::from_euler_zyx`] gives the same rotation (to ~2·10⁻⁸ within 10⁻⁸ rad
    /// of the lock, to ~10⁻¹⁵ away from it).
    ///
    /// Near pitch ±90° (gimbal lock, `cos(pitch) < 10⁻⁸`) only yaw ∓ roll is defined;
    /// roll is then 0 and yaw carries the whole rotation about the vertical:
    /// `yaw = atan2(−C₀₁, C₁₁)` (both locks give the same formula).
    pub fn to_euler_zyx(self) -> (f64, f64, f64) {
        let m = self.to_mat3().0;
        // atan2 instead of asin(−C₂₀): asin is ill-conditioned next to ±1, where a
        // 10⁻¹⁶ error in C₂₀ moves the angle by 10⁻⁸ rad.
        let cp = libm::hypot(m[0][0], m[1][0]);
        let pitch = libm::atan2(-m[2][0], cp);
        if cp < GIMBAL_LOCK_COS {
            let yaw = libm::atan2(-m[0][1], m[1][1]);
            return (yaw, pitch, 0.0);
        }
        let yaw = libm::atan2(m[1][0], m[0][0]);
        let roll = libm::atan2(m[2][1], m[2][2]);
        (yaw, pitch, roll)
    }

    /// The conjugate `(w, −x, −y, −z)`: the inverse rotation of a unit quaternion
    /// (`q_b_a = q_a_b*`).
    pub fn conj(self) -> Self {
        Self::new(self.w, -self.x, -self.y, -self.z)
    }

    /// `|q|`.
    pub fn norm(self) -> f64 {
        let s = self.w * self.w + self.x * self.x;
        let s = s + self.y * self.y;
        (s + self.z * self.z).sqrt()
    }

    /// `q/|q|`, with the sign chosen so that `w ≥ 0` (q and −q are the same rotation).
    pub fn normalized(self) -> Self {
        let n = self.norm();
        let n = if self.w < 0.0 { -n } else { n };
        let k = 1.0 / n;
        Self::new(self.w * k, self.x * k, self.y * k, self.z * k)
    }

    /// `q ⊗ (0, v) ⊗ q*` for a unit quaternion, as
    /// `t = 2·(q_v × v)`, `v + w·t + q_v × t`.
    pub fn rotate(self, v: Vec3) -> Vec3 {
        let qv = self.vector();
        let t = qv.cross(v) * 2.0;
        v + t * self.w + qv.cross(t)
    }

    /// The direction-cosine matrix `C` with `C·v = self.rotate(v)` (unit quaternion).
    pub fn to_mat3(self) -> Mat3 {
        let Self { w, x, y, z } = self;
        let (xx, yy, zz) = (x * x, y * y, z * z);
        let (xy, xz, yz) = (x * y, x * z, y * z);
        let (wx, wy, wz) = (w * x, w * y, w * z);
        Mat3([
            [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy)],
            [2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
            [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)],
        ])
    }

    /// The unit quaternion (w ≥ 0) of a proper rotation matrix (Shepperd's method:
    /// the largest of 4w², 4x², 4y², 4z² picks the stable formula).
    pub fn from_mat3(m: &Mat3) -> Self {
        let m = &m.0;
        let tr = m[0][0] + m[1][1];
        let tr = tr + m[2][2];
        let q = if tr >= m[0][0] && tr >= m[1][1] && tr >= m[2][2] {
            let s = (1.0 + tr).sqrt() * 2.0;
            Self::new(
                0.25 * s,
                (m[2][1] - m[1][2]) / s,
                (m[0][2] - m[2][0]) / s,
                (m[1][0] - m[0][1]) / s,
            )
        } else if m[0][0] >= m[1][1] && m[0][0] >= m[2][2] {
            let s = ((1.0 + m[0][0]) - m[1][1] - m[2][2]).sqrt() * 2.0;
            Self::new(
                (m[2][1] - m[1][2]) / s,
                0.25 * s,
                (m[0][1] + m[1][0]) / s,
                (m[0][2] + m[2][0]) / s,
            )
        } else if m[1][1] >= m[2][2] {
            let s = ((1.0 + m[1][1]) - m[0][0] - m[2][2]).sqrt() * 2.0;
            Self::new(
                (m[0][2] - m[2][0]) / s,
                (m[0][1] + m[1][0]) / s,
                0.25 * s,
                (m[1][2] + m[2][1]) / s,
            )
        } else {
            let s = ((1.0 + m[2][2]) - m[0][0] - m[1][1]).sqrt() * 2.0;
            Self::new(
                (m[1][0] - m[0][1]) / s,
                (m[0][2] + m[2][0]) / s,
                (m[1][2] + m[2][1]) / s,
                0.25 * s,
            )
        };
        q.normalized()
    }
}

impl Mul for Quat {
    type Output = Self;

    /// The Hamilton product `self ⊗ o`.
    fn mul(self, o: Self) -> Self {
        let (a, b) = (self, o);
        let w = a.w * b.w - a.x * b.x;
        let w = w - a.y * b.y;
        let w = w - a.z * b.z;
        let x = a.w * b.x + a.x * b.w;
        let x = x + a.y * b.z;
        let x = x - a.z * b.y;
        let y = a.w * b.y - a.x * b.z;
        let y = y + a.y * b.w;
        let y = y + a.z * b.x;
        let z = a.w * b.z + a.x * b.y;
        let z = z - a.y * b.x;
        let z = z + a.z * b.w;
        Self::new(w, x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const D: f64 = core::f64::consts::PI / 180.0;

    fn close(a: Vec3, b: Vec3, tol: f64) -> bool {
        a.distance(b) <= tol
    }

    #[test]
    fn hamilton_units() {
        let i = Quat::new(0.0, 1.0, 0.0, 0.0);
        let j = Quat::new(0.0, 0.0, 1.0, 0.0);
        let k = Quat::new(0.0, 0.0, 0.0, 1.0);
        let m1 = Quat::new(-1.0, 0.0, 0.0, 0.0);
        assert_eq!(i * i, m1);
        assert_eq!(j * j, m1);
        assert_eq!(k * k, m1);
        assert_eq!(i * j, k);
        assert_eq!(j * k, i);
        assert_eq!(k * i, j);
        assert_eq!(i * j * k, m1);
        assert_eq!(j * i, Quat::new(0.0, 0.0, 0.0, -1.0));
    }

    #[test]
    fn active_right_handed_rotation() {
        // +90° about z takes x to y (counter-clockwise seen from +z).
        let q = Quat::from_axis_angle(Vec3::Z, 90.0 * D);
        assert!(close(q.rotate(Vec3::X), Vec3::Y, 1e-15));
        assert!(close(q.rotate(Vec3::Y), -Vec3::X, 1e-15));
        // +90° about x takes y to z; about y takes z to x.
        let q = Quat::from_axis_angle(Vec3::X, 90.0 * D);
        assert!(close(q.rotate(Vec3::Y), Vec3::Z, 1e-15));
        let q = Quat::from_axis_angle(Vec3::Y, 90.0 * D);
        assert!(close(q.rotate(Vec3::Z), Vec3::X, 1e-15));
    }

    #[test]
    fn aerospace_euler_cases() {
        // Heading 090: the nose points east.
        let q = Quat::from_euler_zyx(90.0 * D, 0.0, 0.0);
        assert!(close(q.rotate(Vec3::X), Vec3::Y, 1e-15));
        // Pitch +10°: the nose is above the horizon, so its down component is negative.
        let q = Quat::from_euler_zyx(0.0, 10.0 * D, 0.0);
        let (s, c) = libm::sincos(10.0 * D);
        assert!(close(q.rotate(Vec3::X), Vec3::new(c, 0.0, -s), 1e-15));
        // Roll +30°: the right wing goes down.
        let q = Quat::from_euler_zyx(0.0, 0.0, 30.0 * D);
        let (s, c) = libm::sincos(30.0 * D);
        assert!(close(q.rotate(Vec3::Y), Vec3::new(0.0, c, s), 1e-15));
        // Combined, hand-computed: heading 090, pitch 10°: nose = (0, cos10°, −sin10°).
        let q = Quat::from_euler_zyx(90.0 * D, 10.0 * D, 0.0);
        let (s, c) = libm::sincos(10.0 * D);
        assert!(close(q.rotate(Vec3::X), Vec3::new(0.0, c, -s), 1e-15));
        // Heading 090, roll 30°: right wing = (−cos30°, 0, sin30°) (points south, down).
        let q = Quat::from_euler_zyx(90.0 * D, 0.0, 30.0 * D);
        let (s, c) = libm::sincos(30.0 * D);
        assert!(close(q.rotate(Vec3::Y), Vec3::new(-c, 0.0, s), 1e-15));
    }

    #[test]
    fn composition_order() {
        // q_a_c = q_a_b ⊗ q_b_c applies q_b_c first.
        let q_ab = Quat::from_axis_angle(Vec3::Z, 90.0 * D);
        let q_bc = Quat::from_axis_angle(Vec3::X, 90.0 * D);
        let q_ac = q_ab * q_bc;
        let v = Vec3::new(0.3, -1.2, 2.5);
        assert!(close(q_ac.rotate(v), q_ab.rotate(q_bc.rotate(v)), 1e-14));
        // y →(x 90°) z →(z 90°) z.
        assert!(close(q_ac.rotate(Vec3::Y), Vec3::Z, 1e-15));
        // The conjugate inverts.
        assert!(close(q_ac.conj().rotate(q_ac.rotate(v)), v, 1e-14));
    }

    /// Review regression: at exactly pitch +90° the split used to give
    /// (0, π/2, 0) for a quaternion with a 90° yaw-minus-roll, a different rotation.
    #[test]
    fn euler_at_gimbal_lock_keeps_the_rotation() {
        let mut cases = vec![
            Quat::new(0.5, 0.5, 0.5, -0.5),
            Quat::new(0.5, -0.5, 0.5, 0.5),
        ];
        for &(y, p, r) in &[
            (30.0, 90.0, 10.0),
            (-120.0, -90.0, 45.0),
            (170.0, 90.0, -170.0),
            (10.0, 89.999_999_9, 20.0),
            (10.0, -89.999_999_99, -20.0),
            (10.0, 89.999_999_5, 20.0),
            (0.0, 90.0 - 1e-12, 0.0),
        ] {
            cases.push(Quat::from_euler_zyx(y * D, p * D, r * D));
        }
        for q in cases {
            let q = q.normalized();
            let (y, p, r) = q.to_euler_zyx();
            let back = Quat::from_euler_zyx(y, p, r);
            for v in [Vec3::X, Vec3::Y, Vec3::Z] {
                assert!(
                    close(back.rotate(v), q.rotate(v), 5e-8),
                    "{q:?}: ({y}, {p}, {r})"
                );
            }
        }
    }

    #[test]
    fn matrix_and_euler_round_trips() {
        for &(y, p, r) in &[
            (0.0, 0.0, 0.0),
            (30.0, 10.0, -20.0),
            (-170.0, -80.0, 179.0),
            (123.0, 45.0, 60.0),
            (-45.0, 89.0, 5.0),
        ] {
            let q = Quat::from_euler_zyx(y * D, p * D, r * D);
            assert!((q.norm() - 1.0).abs() < 1e-15);
            let m = q.to_mat3();
            for v in [Vec3::X, Vec3::Y, Vec3::Z, Vec3::new(1.0, -2.0, 0.5)] {
                assert!(close(m.mul_vec(v), q.rotate(v), 1e-14));
            }
            let back = Quat::from_mat3(&m);
            let qn = q.normalized();
            for (a, b) in [
                (back.w, qn.w),
                (back.x, qn.x),
                (back.y, qn.y),
                (back.z, qn.z),
            ] {
                assert!((a - b).abs() < 1e-14, "{back:?} {qn:?}");
            }
            let (y2, p2, r2) = q.to_euler_zyx();
            assert!((y2 - y * D).abs() < 1e-9 && (p2 - p * D).abs() < 1e-9);
            assert!((r2 - r * D).abs() < 1e-9);
        }
    }
}
