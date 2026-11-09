//! Body and antenna frames (contracts §12.3, §12.4).
//!
//! - **Body**: x forward (nose), y right (starboard wing), z down. A platform's
//!   attitude is `q_ned_body` relative to the local NED frame at the platform's own
//!   position; [`crate::Quat::from_euler_zyx`] builds it from heading, pitch and roll.
//! - **Antenna**: the mount orientation `q_body_ant` maps antenna-frame coordinates to
//!   body coordinates. Pattern conventions in the antenna frame (boresight axis,
//!   θ/φ) belong to `iqsim-antenna` (prompt 10).
//! - A mount's **lever arm** is the antenna phase centre's offset from the platform
//!   origin in body coordinates, metres. [`PlatformPose::antenna_position`] applies
//!   it; the light-time solver starts using phase centres in prompt 10.

use crate::geodetic::{Ecef, Geodetic};
use crate::ned::LocalNed;
use crate::quat::Quat;
use crate::vec::Vec3;

/// A platform's position and attitude at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlatformPose {
    ned: LocalNed,
    q_ned_body: Quat,
    q_ecef_body: Quat,
}

impl PlatformPose {
    /// The pose at geodetic `position` with attitude `q_ned_body` (normalized here).
    pub fn new(position: Geodetic, q_ned_body: Quat) -> Self {
        let ned = LocalNed::new(position);
        let q_ned_body = q_ned_body.normalized();
        let q_ecef_body = (ned.q_ecef_ned() * q_ned_body).normalized();
        Self {
            ned,
            q_ned_body,
            q_ecef_body,
        }
    }

    /// The local NED frame at the platform.
    pub const fn ned(&self) -> &LocalNed {
        &self.ned
    }

    /// The attitude relative to local NED.
    pub const fn q_ned_body(&self) -> Quat {
        self.q_ned_body
    }

    /// `q_ecef_body = q_ecef_ned ⊗ q_ned_body`.
    pub const fn q_ecef_body(&self) -> Quat {
        self.q_ecef_body
    }

    /// An ECEF direction in body coordinates.
    pub fn ecef_to_body(&self, v_ecef: Vec3) -> Vec3 {
        self.q_ecef_body.conj().rotate(v_ecef)
    }

    /// A body direction in ECEF coordinates.
    pub fn body_to_ecef(&self, v_body: Vec3) -> Vec3 {
        self.q_ecef_body.rotate(v_body)
    }

    /// An ECEF direction in the antenna frame of `mount`.
    pub fn ecef_to_antenna(&self, mount: &Mount, v_ecef: Vec3) -> Vec3 {
        mount.q_body_ant.conj().rotate(self.ecef_to_body(v_ecef))
    }

    /// An antenna-frame direction of `mount` in ECEF coordinates.
    pub fn antenna_to_ecef(&self, mount: &Mount, v_ant: Vec3) -> Vec3 {
        self.body_to_ecef(mount.q_body_ant.rotate(v_ant))
    }

    /// The ECEF position of `mount`'s phase centre (origin + lever arm).
    pub fn antenna_position(&self, mount: &Mount) -> Ecef {
        Ecef::from_vec(self.ned.origin().vec() + self.body_to_ecef(mount.lever_arm_body_m))
    }
}

/// An antenna mount: phase-centre lever arm and orientation in the body frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mount {
    /// Phase-centre offset from the platform origin, body coordinates, metres.
    pub lever_arm_body_m: Vec3,
    /// Antenna-frame → body-frame rotation.
    pub q_body_ant: Quat,
}

impl Mount {
    /// A mount with the given lever arm and orientation (normalized here).
    pub fn new(lever_arm_body_m: Vec3, q_body_ant: Quat) -> Self {
        Self {
            lever_arm_body_m,
            q_body_ant: q_body_ant.normalized(),
        }
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
    fn body_axes_at_the_equator_hand_computed() {
        // At (0°, 0°) NED = (+z, +y, −x) in ECEF. Heading 090, level: nose east (+y),
        // right wing south (−z), belly down (−x).
        let pose = PlatformPose::new(
            Geodetic::new(0.0, 0.0, 10_000.0),
            Quat::from_euler_zyx(90.0 * D, 0.0, 0.0),
        );
        assert!(close(pose.body_to_ecef(Vec3::X), Vec3::Y, 1e-15));
        assert!(close(pose.body_to_ecef(Vec3::Y), -Vec3::Z, 1e-15));
        assert!(close(pose.body_to_ecef(Vec3::Z), -Vec3::X, 1e-15));
        // Banked 30° right while heading 090: the right wing (south) goes down, so the
        // belly tilts north: NED (sin30°, 0, cos30°) = ECEF (−cos30°, 0, sin30°).
        let pose = PlatformPose::new(
            Geodetic::new(0.0, 0.0, 10_000.0),
            Quat::from_euler_zyx(90.0 * D, 0.0, 30.0 * D),
        );
        let (s, c) = libm::sincos(30.0 * D);
        assert!(close(
            pose.body_to_ecef(Vec3::Z),
            Vec3::new(-c, 0.0, s),
            1e-15
        ));
        // A direction round-trips.
        let v = Vec3::new(0.2, -0.7, 0.4);
        assert!(close(pose.ecef_to_body(pose.body_to_ecef(v)), v, 1e-15));
    }

    #[test]
    fn antenna_frame_and_lever_arm() {
        // A belly antenna whose antenna z axis is the body z axis (down), rotated 90°
        // about body z so antenna x points along the right wing.
        let mount = Mount::new(
            Vec3::new(5.0, 0.0, 1.5),
            Quat::from_axis_angle(Vec3::Z, 90.0 * D),
        );
        let pose = PlatformPose::new(Geodetic::new(0.0, 0.0, 0.0), Quat::IDENTITY);
        // Heading north at (0°, 0°): body x = +z ECEF, y = +y, z = −x.
        assert!(close(pose.antenna_to_ecef(&mount, Vec3::X), Vec3::Y, 1e-15));
        assert!(close(
            pose.antenna_to_ecef(&mount, Vec3::Z),
            -Vec3::X,
            1e-15
        ));
        assert!(close(pose.ecef_to_antenna(&mount, Vec3::Y), Vec3::X, 1e-15));
        // Lever arm: 5 m forward (+z) and 1.5 m down (−x) from the origin.
        let o = pose.ned().origin();
        let p = pose.antenna_position(&mount);
        assert!(close(p.vec() - o.vec(), Vec3::new(-1.5, 0.0, 5.0), 1e-9));
    }
}
