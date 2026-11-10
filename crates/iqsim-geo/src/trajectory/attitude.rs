//! Attitude derived from the kinematics (plan "World model": Derived attitude;
//! contracts §13.5).
//!
//! From position, velocity, acceleration and jerk (ECEF), with `n` the geodetic up
//! normal at the position:
//!
//! - vertical rate `ḣ = v·n`, horizontal velocity `v_h = v − ḣn`, `ṅ = K v_h`
//!   (curvature map), `ḧ = a·n + v·ṅ`;
//! - air-relative horizontal velocity `u = v_h − w` for a constant wind `w` (NED);
//!   the **heading frame** H has x along `u` (heading; the wind crab is the angle from
//!   the track), y to the right, z down;
//! - **pitch** `θ = atan2(ḣ, |u|) + α` (flight-path angle plus angle of attack);
//! - **bank** `φ = atan2(a·y_H, g₀)`: the lateral acceleration a·y_H is `V ω` in a
//!   coordinated turn, so `tan φ = V ω / g₀` with g₀ = 9.806 65 m/s²;
//! - `q_ecef_body = q_ecef_H ⊗ q_y(θ) ⊗ q_x(φ)` (heading, then pitch, then bank: the
//!   aerospace 3-2-1 sequence relative to NED);
//! - angular rate (body relative to ECEF) `Ω_H + θ̇ y_H + φ̇ x_body`, with
//!   `Ω_H = ½ Σ hᵢ × ḣᵢ` over the H axes; every derivative is in closed form, using
//!   the jerk for φ̇.
//!
//! [`AttitudeMode::Ground`] has no bank and no angle of attack (pitch follows the
//! along-track slope). [`AttitudeMode::Fixed`] holds a constant attitude relative to
//! NED and is only for platforms that do not move.

use crate::GeoError;
use crate::geodetic::Ecef;
use crate::ned::LocalNed;
use crate::quat::Quat;
use crate::trajectory::earth::{curvature, north_east, radii};
use crate::trajectory::piece::Kinematics;
use crate::vec::{Mat3, Vec3};

/// Standard gravity g₀, m/s², used in the bank formula.
pub const STANDARD_GRAVITY_M_S2: f64 = 9.806_65;

/// Smallest air-relative horizontal speed with a defined heading, m/s.
pub const MIN_HEADING_SPEED_M_S: f64 = 0.1;

/// Wind needs north and east: rejected closer to a pole than this cos φ (≈ 0.01°).
const WIND_MIN_COS: f64 = 1.745e-4;

/// How a platform's attitude follows its motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AttitudeMode {
    /// Coordinated flight: heading from the air-relative track, pitch from the
    /// flight-path angle plus `angle_of_attack_rad`, bank from the lateral
    /// acceleration. `wind_ned_m_s` is the constant air velocity (towards), NED, m/s.
    Derived {
        /// Angle of attack added to the pitch, radians.
        angle_of_attack_rad: f64,
        /// Wind velocity (where the air moves to), NED, m/s. Down is ignored.
        wind_ned_m_s: Vec3,
    },
    /// Ground vehicle: heading along the track, pitch along the slope, no bank.
    Ground,
    /// A constant attitude relative to the local NED frame (static platforms).
    Fixed(Quat),
}

impl AttitudeMode {
    /// Coordinated flight with no wind and no angle of attack.
    pub const COORDINATED: Self = Self::Derived {
        angle_of_attack_rad: 0.0,
        wind_ned_m_s: Vec3::ZERO,
    };
}

/// Attitude and angular rate (body axes, relative to ECEF) for kinematics `k`.
pub fn derive(mode: &AttitudeMode, k: &Kinematics) -> Result<(Quat, Vec3), GeoError> {
    let (alpha, wind_ned, banks) = match *mode {
        AttitudeMode::Fixed(q) => {
            let ned = LocalNed::new(Ecef::from_vec(k.p).to_geodetic()?);
            return Ok(((ned.q_ecef_ned() * q).normalized(), Vec3::ZERO));
        }
        AttitudeMode::Ground => (0.0, Vec3::ZERO, false),
        AttitudeMode::Derived {
            angle_of_attack_rad,
            wind_ned_m_s,
        } => (angle_of_attack_rad, wind_ned_m_s, true),
    };
    let g = Ecef::from_vec(k.p).to_geodetic()?;
    let (sp, cp) = libm::sincos(g.lat_rad);
    let (sl, cl) = libm::sincos(g.lon_rad);
    let n = Vec3::new(cp * cl, cp * sl, sp);
    let h = g.h_m;

    let hdot = k.v.dot(n);
    let v_h = k.v - n * hdot;
    let ndot = curvature(n, h, v_h);
    let hddot = k.a.dot(n) + k.v.dot(ndot);
    let v_h_dot = k.a - n * hddot - ndot * hdot;

    // Wind in ECEF and its rate (the NED frame turns as the platform moves).
    let (w, w_dot) = if wind_ned.x == 0.0 && wind_ned.y == 0.0 {
        (Vec3::ZERO, Vec3::ZERO)
    } else {
        let (north, east) = north_east(n, WIND_MIN_COS).ok_or_else(|| {
            GeoError::Trajectory("wind is undefined within 0.01° of a pole".to_owned())
        })?;
        let (big_n, m) = radii(n);
        let lat_dot = v_h.dot(north) / (m + h);
        let lon_dot = v_h.dot(east) / ((big_n + h) * cp);
        let rho = Vec3::Z * lon_dot - east * lat_dot;
        let w = north * wind_ned.x + east * wind_ned.y;
        (w, rho.cross(w))
    };

    let u = v_h - w;
    let u_dot = v_h_dot - w_dot;
    let speed = u.norm();
    if speed.is_nan() || speed < MIN_HEADING_SPEED_M_S {
        return Err(GeoError::Trajectory(format!(
            "attitude is undefined at an air-relative horizontal speed of {speed} m/s (minimum {MIN_HEADING_SPEED_M_S} m/s)"
        )));
    }
    let h1 = u * (1.0 / speed);
    let speed_dot = h1.dot(u_dot);
    let h1_dot = (u_dot - h1 * speed_dot) * (1.0 / speed);
    let h3 = -n;
    let h3_dot = -ndot;
    let h2 = h3.cross(h1);
    let h2_dot = h3_dot.cross(h1) + h3.cross(h1_dot);
    let omega_h = (h1.cross(h1_dot) + h2.cross(h2_dot) + h3.cross(h3_dot)) * 0.5;

    let gamma = libm::atan2(hdot, speed);
    let gamma_dot = (hddot * speed - hdot * speed_dot) / (speed * speed + hdot * hdot);
    let theta = gamma + alpha;

    let (bank, bank_dot) = if banks {
        let g0 = STANDARD_GRAVITY_M_S2;
        let lat = k.a.dot(h2);
        let lat_dot = k.j.dot(h2) + k.a.dot(h2_dot);
        (libm::atan2(lat, g0), g0 * lat_dot / (g0 * g0 + lat * lat))
    } else {
        (0.0, 0.0)
    };

    let c_ecef_h = Mat3::from_rows(
        Vec3::new(h1.x, h2.x, h3.x),
        Vec3::new(h1.y, h2.y, h3.y),
        Vec3::new(h1.z, h2.z, h3.z),
    );
    let q_h = Quat::from_mat3(&c_ecef_h);
    let q_hp = q_h * Quat::from_axis_angle(Vec3::Y, theta);
    let q = (q_hp * Quat::from_axis_angle(Vec3::X, bank)).normalized();
    let x_body = q_hp.rotate(Vec3::X);
    let omega_ecef = omega_h + h2 * gamma_dot + x_body * bank_dot;
    Ok((q, q.conj().rotate(omega_ecef)))
}
