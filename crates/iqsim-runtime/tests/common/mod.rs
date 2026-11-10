//! Shared test references.

#![allow(dead_code)]

/// Speed of light, m/s.
pub const C: f64 = 299_792_458.0;
/// WGS84 Earth rotation rate, rad/s.
pub const OMEGA_E: f64 = 7.292_115e-5;

/// Light time between two Earth-fixed points (ECEF metres), receiver first, solved in
/// the rotating Earth's inertial frame independently of the simulator's solver.
///
/// In ECI the receiver at t is `R_z(θ(t))·p_rx` and the transmitter at t − τ is
/// `R_z(θ(t) − ω_E·τ)·p_tx`, so `c·τ = |p_rx − R_z(−ω_E·τ)·p_tx|` for every t (the
/// Sagnac effect). Solved by fixed-point iteration to convergence in f64.
pub fn eci_delay(rx: [f64; 3], tx: [f64; 3]) -> f64 {
    let mut tau = 0.0_f64;
    for _ in 0..50 {
        let (s, c) = libm::sincos(-OMEGA_E * tau);
        let t = [c * tx[0] - s * tx[1], s * tx[0] + c * tx[1], tx[2]];
        let d = [rx[0] - t[0], rx[1] - t[1], rx[2] - t[2]];
        tau = libm::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / C;
    }
    tau
}

use iqsim_core::{AcqInstant, TruthTime};
use iqsim_geo::trajectory::{Span, State};
use iqsim_geo::{EarthOrientation, Eci, GeoError, Quat, Trajectory, Vec3};

/// GM of the Earth (WGS84), m³/s².
pub const MU: f64 = 3.986_004_418e14;

/// Motion defined in ECI in closed form, for analytic light-time references.
#[derive(Clone, Copy, Debug)]
pub enum EciMotion {
    /// `p0 + u·t`.
    Uniform {
        /// Position at t = 0, m.
        p0: [f64; 3],
        /// Velocity, m/s.
        u: [f64; 3],
    },
    /// A circular equatorial orbit of radius `a`, angle `phase0 + n·t` (prograde).
    Circular {
        /// Radius, m.
        a: f64,
        /// Angle at t = 0, rad.
        phase0: f64,
    },
}

/// Seconds of an instant (f64; the tests use small times).
pub fn seconds(t: AcqInstant) -> f64 {
    (t.ps().as_ps() as f64 + t.frac() as f64 / 18_446_744_073_709_551_616.0) / 1e12
}

impl EciMotion {
    /// ECI position, velocity and acceleration at `t` seconds.
    pub fn eci(&self, t: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        match *self {
            Self::Uniform { p0, u } => (
                [p0[0] + u[0] * t, p0[1] + u[1] * t, p0[2] + u[2] * t],
                u,
                [0.0; 3],
            ),
            Self::Circular { a, phase0 } => {
                let n = libm::sqrt(MU / (a * a * a));
                let (s, c) = libm::sincos(phase0 + n * t);
                (
                    [a * c, a * s, 0.0],
                    [-a * n * s, a * n * c, 0.0],
                    [-a * n * n * c, -a * n * n * s, 0.0],
                )
            }
        }
    }
}

impl Trajectory for EciMotion {
    fn state(&self, t: AcqInstant) -> Result<State, GeoError> {
        let (r, v, a) = self.eci(seconds(t));
        let o = EarthOrientation::at(t);
        let w = Vec3::new(0.0, 0.0, OMEGA_E);
        let (r, v, a) = (
            Vec3::new(r[0], r[1], r[2]),
            Vec3::new(v[0], v[1], v[2]),
            Vec3::new(a[0], a[1], a[2]),
        );
        let p = o.eci_to_ecef(Eci::from_vec(r));
        let rp = p.vec();
        let vp = o.vector_to_ecef(v) - w.cross(rp);
        let ap = o.vector_to_ecef(a) - w.cross(vp) * 2.0 - w.cross(w.cross(rp));
        Ok(State {
            position: p,
            velocity: vp,
            acceleration: ap,
            q_ecef_body: Quat::IDENTITY,
            angular_rate_body: Vec3::ZERO,
        })
    }

    fn span(&self) -> Span {
        Span {
            start: None,
            end: None,
        }
    }

    fn knots(&self, _from: TruthTime, _to: TruthTime) -> Vec<TruthTime> {
        Vec::new()
    }
}

/// Reference light time from an Earth-fixed receiver `rx` (ECEF) at receive time `t`
/// seconds to an emitter moving in ECI: `c·τ = |R_z(ω_E·t)·rx − r(t − τ)|`, iterated to
/// convergence in f64. Independent of the simulator's solver and Earth-angle routine.
pub fn reference_delay(rx: [f64; 3], tx: &EciMotion, t: f64) -> f64 {
    let (s, c) = libm::sincos(OMEGA_E * t);
    let r = [c * rx[0] - s * rx[1], s * rx[0] + c * rx[1], rx[2]];
    let mut tau = 0.0_f64;
    for _ in 0..30 {
        let (p, _, _) = tx.eci(t - tau);
        let d = [r[0] - p[0], r[1] - p[1], r[2] - p[2]];
        tau = libm::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / C;
    }
    tau
}
