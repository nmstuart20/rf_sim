//! Platform motion: the [`Trajectory`] trait and its sources (contracts §13).
//!
//! A trajectory is a pure function of time. [`Trajectory::state`] returns position,
//! velocity, acceleration, attitude and angular rate at any instant in its
//! [`Span`], and the same call always returns the same bits.
//!
//! Sources:
//!
//! - [`Static`]: a fixed position and attitude.
//! - [`Flight`]: everything that moves. It covers maneuver primitives
//!   ([`Flight::maneuvers`]), waypoint routes ([`Flight::route`]), recorded tracks
//!   ([`Flight::track`]) and ground-clamped motion ([`AttitudeMode::Ground`] with a
//!   [`Ground`] surface). [`Flight::append`] adds primitives from the current state
//!   at a given time (runtime maneuvers).
//!
//! Every [`Flight`] stores its kinematics as quintic Hermite pieces in ECEF, with knots
//! at whole picoseconds ([`piece`]). Velocity, acceleration and jerk are the exact
//! derivatives of position, and position, velocity and acceleration agree at every
//! knot, so motion is C2. Attitude is derived from the kinematics ([`attitude`]).
//!
//! All arithmetic is f64 through `libm`, with each product and sum a separate
//! operation in a fixed order.

pub mod attitude;
mod earth;
mod flight;
pub mod piece;
pub mod profile;
mod route;
pub mod spline;

use std::fmt;

use iqsim_core::{AcqInstant, TruthTime};
use iqsim_units::{Altitude, AltitudeRef};

use crate::GeoError;
use crate::altitude::AltitudeModel;
use crate::frames::PlatformPose;
use crate::geodetic::{Ecef, Geodetic};
use crate::ned::LocalNed;
use crate::quat::Quat;
use crate::vec::Vec3;

pub use attitude::{AttitudeMode, STANDARD_GRAVITY_M_S2};
pub use flight::{
    ClimbEnd, ClimbRate, Flight, Limits, Primitive, Start, TrackPoint, TurnDirection, TurnEnd,
    TurnRate, Waypoint,
};

/// The kinematic state of a platform at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    /// Platform origin, ECEF metres.
    pub position: Ecef,
    /// Earth-relative velocity in ECEF axes, m/s.
    pub velocity: Vec3,
    /// Earth-relative acceleration (second derivative of the ECEF position), m/s².
    pub acceleration: Vec3,
    /// Attitude: body-frame coordinates → ECEF coordinates (contracts §12.4).
    pub q_ecef_body: Quat,
    /// Angular velocity of the body frame relative to ECEF, in body axes, rad/s.
    pub angular_rate_body: Vec3,
}

impl State {
    /// The attitude relative to the local NED frame at the platform.
    pub fn q_ned_body(&self) -> Result<Quat, GeoError> {
        let ned = LocalNed::new(self.position.to_geodetic()?);
        Ok((ned.q_ecef_ned().conj() * self.q_ecef_body).normalized())
    }

    /// The pose (position and attitude) for antenna-frame transforms.
    pub fn pose(&self) -> Result<PlatformPose, GeoError> {
        let g = self.position.to_geodetic()?;
        Ok(PlatformPose::new(g, self.q_ned_body()?))
    }
}

/// The interval on which a trajectory is defined (both ends inclusive; `None` is
/// unbounded).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span {
    /// First instant, if bounded.
    pub start: Option<TruthTime>,
    /// Last instant, if bounded.
    pub end: Option<TruthTime>,
}

impl Span {
    /// Whether `t` lies in the span.
    pub fn contains(&self, t: AcqInstant) -> bool {
        self.start
            .is_none_or(|s| AcqInstant::from_truth_time(s) <= t)
            && self.end.is_none_or(|e| t <= AcqInstant::from_truth_time(e))
    }
}

/// A platform's motion and attitude as a pure function of time.
pub trait Trajectory: Send + Sync + fmt::Debug {
    /// The state at acquisition instant `t` (fractional picoseconds allowed). Outside
    /// [`span`](Self::span) it is [`GeoError::OutsideSpan`].
    fn state(&self, t: AcqInstant) -> Result<State, GeoError>;

    /// The state at whole-picosecond instant `t`.
    fn state_at(&self, t: TruthTime) -> Result<State, GeoError> {
        self.state(AcqInstant::from_truth_time(t))
    }

    /// Where the trajectory is defined.
    fn span(&self) -> Span;

    /// Knots in `[from, to]`, ascending: the instants where jerk (and so angular rate)
    /// may be discontinuous. Polynomial fits should split there (prompt 9).
    fn knots(&self, from: TruthTime, to: TruthTime) -> Vec<TruthTime>;
}

/// `t − t0` in seconds: an exact Q64.64 picosecond difference, rounded once to f64
/// picoseconds, then divided by 10¹². Callers keep `t − t0` small (one piece), so no
/// large absolute time ever meets a floating-point product.
pub(crate) fn seconds_since(t: AcqInstant, t0: TruthTime) -> f64 {
    let ps = i128::from(t.ps().as_ps()) - i128::from(t0.as_ps());
    let frac = t.frac() as f64 * (1.0 / 18_446_744_073_709_551_616.0);
    (ps as f64 + frac) / 1.0e12
}

/// The ground under a ground-clamped platform: its ellipsoid height at a point.
///
/// Placeholder until terrain (prompt 22): [`FlatGround`] is a constant height. A
/// terrain model must be C2 along the path, or the clamped motion is not.
pub trait Ground: Send + Sync + fmt::Debug {
    /// Ellipsoid height of the ground at (φ, λ) radians, metres.
    fn height_m(&self, lat_rad: f64, lon_rad: f64) -> f64;
}

/// Ground at a constant ellipsoid height everywhere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatGround {
    /// Ground height above the ellipsoid, metres.
    pub height_m: f64,
}

impl FlatGround {
    /// The ellipsoid itself, which is how AGL heights resolve until terrain exists
    /// (contracts §11.4).
    pub const ELLIPSOID: Self = Self { height_m: 0.0 };
}

impl Ground for FlatGround {
    fn height_m(&self, _lat_rad: f64, _lon_rad: f64) -> f64 {
        self.height_m
    }
}

/// How scenario altitudes become ellipsoid heights: the geoid and atmosphere for MSL
/// and flight levels, the ground for AGL.
#[derive(Clone, Copy, Debug)]
pub struct Heights<'a> {
    /// Geoid and atmosphere.
    pub model: &'a AltitudeModel,
    /// Ground for AGL heights.
    pub ground: &'a dyn Ground,
}

impl Heights<'_> {
    /// The ellipsoid height of `alt` at (φ, λ) radians.
    pub fn ellipsoid_height_m(
        &self,
        lat_rad: f64,
        lon_rad: f64,
        alt: Altitude,
    ) -> Result<f64, GeoError> {
        match alt.reference() {
            AltitudeRef::Agl => Ok(self.ground.height_m(lat_rad, lon_rad) + alt.meters()),
            _ => self.model.ellipsoid_height_m(lat_rad, lon_rad, alt),
        }
    }
}

/// A platform that does not move: fixed position and attitude relative to the local
/// NED frame, zero velocity, acceleration and angular rate (relative to ECEF).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Static {
    state: State,
}

impl Static {
    /// A static platform at `position` with attitude `q_ned_body`.
    pub fn new(position: Geodetic, q_ned_body: Quat) -> Self {
        let pose = PlatformPose::new(position, q_ned_body);
        Self {
            state: State {
                position: position.to_ecef(),
                velocity: Vec3::ZERO,
                acceleration: Vec3::ZERO,
                q_ecef_body: pose.q_ecef_body(),
                angular_rate_body: Vec3::ZERO,
            },
        }
    }
}

impl Trajectory for Static {
    fn state(&self, _t: AcqInstant) -> Result<State, GeoError> {
        Ok(self.state)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seconds_since_is_exact_for_small_spans() {
        let t0 = TruthTime::from_ps(5_000_000_000_000_000_000);
        let t = AcqInstant::new(TruthTime::from_ps(5_000_000_000_250_000_000), 1 << 63);
        assert_eq!(seconds_since(t, t0), 250_000_000.5 / 1.0e12);
        assert_eq!(seconds_since(AcqInstant::from_truth_time(t0), t0), 0.0);
        let early = AcqInstant::from_truth_time(TruthTime::from_ps(-1_000_000_000_000));
        assert_eq!(seconds_since(early, TruthTime::EPOCH), -1.0);
    }

    #[test]
    fn static_platform_is_constant() {
        let g = Geodetic::from_degrees(51.0, -1.0, 100.0);
        let s = Static::new(g, Quat::from_euler_zyx(0.5, 0.0, 0.0));
        let a = s.state_at(TruthTime::EPOCH).unwrap();
        let b = s.state_at(TruthTime::from_ps(i64::MAX)).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.velocity, Vec3::ZERO);
        let (yaw, _, _) = a.q_ned_body().unwrap().to_euler_zyx();
        assert!((yaw - 0.5).abs() < 1e-12);
        assert!(
            s.span()
                .contains(AcqInstant::from_truth_time(TruthTime::MIN))
        );
    }
}
