//! Moving platforms: maneuver primitives, routes, tracks and appends (contracts §13).
//!
//! A [`Flight`] is a [`Path`] of quintic pieces plus an [`AttitudeMode`]. Maneuvers
//! are flown segment by segment. Each segment drives three channels ([`profile`]):
//! turn rate ω, speed V and vertical rate ḣ (h is its exact integral). The horizontal
//! motion integrates on the ellipsoid in n-vector form (the `earth` module) with classic RK4:
//!
//! ```text
//! ṅ = K(n, h)·V_h d                 (K: curvature map, V_h = √(V² − ḣ²))
//! ḋ = −ω (n × d) − (d·ṅ) n          (ω > 0 turns right; ω = 0 is parallel transport)
//! ```
//!
//! so "straight" is a geodesic of the height surface with zero lateral acceleration,
//! and the right-hand lateral acceleration is V_h ω. Quintic Hermite nodes sit at the
//! segment ends, at every channel breakpoint and at most [`MAX_NODE_SPACING_PS`] apart;
//! node position, velocity and acceleration come from the ODE state and the channels
//! in closed form:
//!
//! ```text
//! p = point(n, h)    v = V_h d + ḣ n    a = V̇_h d + V_h ḋ + ḧ n + ḣ ṅ
//! ```
//!
//! [`profile`]: crate::trajectory::profile

use iqsim_core::{AcqInstant, TruthDuration, TruthTime};
use iqsim_units::Altitude;

use crate::GeoError;
use crate::geodetic::{Ecef, Geodetic};
use crate::trajectory::attitude::{AttitudeMode, STANDARD_GRAVITY_M_S2, derive};
use crate::trajectory::earth::{
    curvature, direction, normal, normal_of, point, track_angle, wrap_pi,
};
use crate::trajectory::piece::{Kinematics, Path, Piece};
use crate::trajectory::profile::{
    Channel, RampShape, ScurveShape, hold_for, move_by, ramp, scurve,
};
use crate::trajectory::spline::CubicSpline3;
use crate::trajectory::{Heights, Span, State, Trajectory};
use crate::vec::Vec3;

/// Largest spacing of Hermite nodes, ps (1 s).
pub const MAX_NODE_SPACING_PS: i64 = 1_000_000_000_000;
/// Nodes closer than this (1 ms) are merged; a channel breakpoint then falls just
/// inside a piece.
pub const MIN_NODE_SPACING_PS: i64 = 1_000_000_000;
/// Largest RK4 step, seconds.
pub const MAX_RK_STEP_S: f64 = 0.125;
/// Largest turn per RK4 step, radians (the step shrinks for fast turns).
pub const MAX_RK_TURN_RAD: f64 = 0.02;
/// Largest turn between Hermite nodes, radians (the spacing shrinks for fast turns).
pub const MAX_NODE_TURN_RAD: f64 = 0.1;
/// Shortest segment, ps (1 ms). Shorter transitions are held to this length; an
/// explicit shorter duration is an error.
pub const MIN_SEGMENT_PS: i64 = MIN_NODE_SPACING_PS;
/// Smallest horizontal ground speed of a maneuvering platform, m/s.
pub const MIN_SPEED_M_S: f64 = 1.0;
/// Largest bank angle a turn may demand, radians (80°).
pub const MAX_BANK_RAD: f64 = 80.0 * core::f64::consts::PI / 180.0;

/// Performance limits and route defaults, SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Limits {
    /// Bank angle of route turns, radians (default 25°).
    pub bank_rad: f64,
    /// Roll rate, rad/s (default 10°/s); the turn rate ramps at g₀·roll_rate/V.
    pub roll_rate_rad_s: f64,
    /// Largest longitudinal acceleration in speed changes, m/s² (default 1).
    pub accel_m_s2: f64,
    /// Longitudinal jerk, m/s³ (default 0.5).
    pub jerk_m_s3: f64,
    /// Largest vertical acceleration, m/s² (default 0.5).
    pub vertical_accel_m_s2: f64,
    /// Vertical jerk, m/s³ (default 0.25).
    pub vertical_jerk_m_s3: f64,
    /// Route climb and descent rate, m/s (default 1500 ft/min = 7.62 m/s).
    pub climb_rate_m_s: f64,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            bank_rad: 25.0_f64.to_radians(),
            roll_rate_rad_s: 10.0_f64.to_radians(),
            accel_m_s2: 1.0,
            jerk_m_s3: 0.5,
            vertical_accel_m_s2: 0.5,
            vertical_jerk_m_s3: 0.25,
            climb_rate_m_s: 7.62,
        }
    }
}

impl Limits {
    fn check(&self) -> Result<(), GeoError> {
        let all = [
            self.bank_rad,
            self.roll_rate_rad_s,
            self.accel_m_s2,
            self.jerk_m_s3,
            self.vertical_accel_m_s2,
            self.vertical_jerk_m_s3,
            self.climb_rate_m_s,
        ];
        if !all.iter().all(|v| v.is_finite() && *v > 0.0) {
            return plan_err("every limit must be positive and finite");
        }
        if self.bank_rad >= MAX_BANK_RAD {
            return plan_err("the route bank angle must be below 80°");
        }
        Ok(())
    }
}

/// Turn direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TurnDirection {
    /// Counter-clockwise seen from above.
    Left,
    /// Clockwise seen from above.
    Right,
}

impl TurnDirection {
    fn sign(self) -> f64 {
        match self {
            Self::Left => -1.0,
            Self::Right => 1.0,
        }
    }
}

/// How hard a turn is (magnitudes).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TurnRate {
    /// Turn rate, rad/s.
    Rate(f64),
    /// Bank angle, radians; the rate is g₀ tan φ / V at the turn's start speed.
    Bank(f64),
}

/// When a turn ends (the turn rate is back to zero).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TurnEnd {
    /// After this long in total (roll-in, hold, roll-out). Needs a direction.
    Duration(TruthDuration),
    /// On this track angle (radians from north); the shortest way unless a direction
    /// is given.
    Heading(f64),
    /// After turning this much relative to a straight path (radians, magnitude). Needs
    /// a direction.
    By(f64),
}

/// How steep a climb is. Signed for [`ClimbEnd::Duration`] (negative descends); for
/// [`ClimbEnd::Altitude`] only the magnitude counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ClimbRate {
    /// Vertical rate, m/s.
    Rate(f64),
    /// Flight-path angle, radians; the rate is V sin γ at the climb's start speed.
    Angle(f64),
}

/// When a climb ends (level again).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ClimbEnd {
    /// Level at this altitude, converted at the climb's start point.
    Altitude(Altitude),
    /// After this long in total.
    Duration(TruthDuration),
}

/// A maneuver primitive. Each one commands all three channels: those it does not
/// name return to rest (wings level, level flight, constant speed), so every
/// primitive ends in steady straight-and-level flight unless it says otherwise. A
/// primitive lasts until its transitions complete, or its duration if that is later.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Primitive {
    /// Straight and level for `duration`.
    Straight {
        /// How long.
        duration: TruthDuration,
    },
    /// A coordinated, level turn.
    Turn {
        /// Turn rate or bank.
        rate: TurnRate,
        /// Direction (required for duration and angle ends).
        direction: Option<TurnDirection>,
        /// End condition.
        end: TurnEnd,
    },
    /// A straight climb or descent.
    Climb {
        /// Vertical rate or flight-path angle.
        rate: ClimbRate,
        /// End condition.
        end: ClimbEnd,
    },
    /// A straight, level change to a new speed.
    Speed {
        /// Target speed, m/s.
        speed_m_s: f64,
    },
}

/// The initial state of a maneuvering platform: straight and level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Start {
    /// Truth time of the initial state; the trajectory starts here.
    pub time: TruthTime,
    /// Initial position (ellipsoid height).
    pub position: Geodetic,
    /// Initial speed, m/s.
    pub speed_m_s: f64,
    /// Initial track angle, radians from north. Routes default to the bearing of the
    /// first waypoint.
    pub track_rad: Option<f64>,
}

/// A route waypoint, flown over.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Waypoint {
    /// Latitude, radians.
    pub lat_rad: f64,
    /// Longitude, radians.
    pub lon_rad: f64,
    /// Altitude at the waypoint (reached by it).
    pub altitude: Altitude,
    /// Speed reached by the waypoint, m/s (default: keep the current speed).
    pub speed_m_s: Option<f64>,
}

/// One recorded track position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TrackPoint {
    /// When the platform was there.
    pub time: TruthTime,
    /// Where it was (ellipsoid height).
    pub position: Geodetic,
}

/// A moving platform (module docs).
#[derive(Clone, Debug, PartialEq)]
pub struct Flight {
    path: Path,
    attitude: AttitudeMode,
    limits: Limits,
}

pub(crate) fn plan_err<T>(m: impl Into<String>) -> Result<T, GeoError> {
    Err(GeoError::Trajectory(m.into()))
}

/// The continuous state at a segment boundary.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Air {
    /// Up normal.
    pub n: Vec3,
    /// Horizontal unit track direction.
    pub d: Vec3,
    /// Ellipsoid height, m.
    pub h: f64,
    /// Speed |v|, m/s.
    pub speed: f64,
    /// dV/dt.
    pub speed_dot: f64,
    /// Vertical rate.
    pub hdot: f64,
    /// Vertical acceleration.
    pub hddot: f64,
    /// Turn rate (right positive).
    pub omega: f64,
}

impl Air {
    /// The state implied by position, velocity and acceleration (exact inverse of the
    /// node formulas, with n and h from the geodetic inverse).
    pub(crate) fn from_kinematics(p: Vec3, v: Vec3, a: Vec3) -> Result<Self, GeoError> {
        let (n, h) = normal_of(Ecef::from_vec(p))?;
        let hdot = v.dot(n);
        let v_h = v - n * hdot;
        let vh = v_h.norm();
        if vh.is_nan() || vh < MIN_SPEED_M_S {
            return plan_err(format!(
                "cannot maneuver from a horizontal speed of {vh} m/s (minimum {MIN_SPEED_M_S} m/s)"
            ));
        }
        let d = v_h * (1.0 / vh);
        let speed = v.norm();
        let ndot = curvature(n, h, v_h);
        let right = d.cross(n);
        Ok(Self {
            n,
            d,
            h,
            speed,
            speed_dot: v.dot(a) / speed,
            hdot,
            hddot: a.dot(n) + v.dot(ndot),
            omega: (a.dot(right) - hdot * ndot.dot(right)) / vh,
        })
    }

    /// Latitude and longitude of the normal, radians.
    pub(crate) fn lat_lon(&self) -> (f64, f64) {
        let n = self.n;
        (
            libm::atan2(n.z, libm::hypot(n.x, n.y)),
            libm::atan2(n.y, n.x),
        )
    }

    /// Track angle, radians.
    pub(crate) fn track(&self) -> Result<f64, GeoError> {
        track_angle(self.n, self.d)
    }

    /// Position.
    pub(crate) fn position(&self) -> Vec3 {
        point(self.n, self.h)
    }
}

/// One planned segment: the three channels from its start.
#[derive(Clone, Debug)]
pub(crate) struct Segment {
    pub omega: Channel,
    pub speed: Channel,
    pub climb: Channel,
    /// Height at the start.
    pub h0: f64,
    /// Duration, ps (rounded from seconds).
    pub duration_ps: i64,
}

impl Segment {
    /// A segment lasting until its channels finish, or `at_least_s` if later.
    pub(crate) fn new(
        omega: Channel,
        speed: Channel,
        climb: Channel,
        h0: f64,
        at_least_s: f64,
    ) -> Result<Self, GeoError> {
        let s = omega
            .duration()
            .max(speed.duration())
            .max(climb.duration())
            .max(at_least_s);
        if !(s.is_finite() && s < 9.0e6) {
            return plan_err(format!("a maneuver lasting {s} s is too long"));
        }
        let mut duration_ps = libm::round(s * 1.0e12) as i64;
        if duration_ps > 0 && duration_ps < MIN_SEGMENT_PS {
            if at_least_s * 1.0e12 > 0.5 {
                return plan_err(format!(
                    "a maneuver duration of {at_least_s} s is below the 1 ms minimum"
                ));
            }
            duration_ps = MIN_SEGMENT_PS;
        }
        Ok(Self {
            omega,
            speed,
            climb,
            h0,
            duration_ps,
        })
    }

    /// The largest |ω| anywhere in the segment, rad/s.
    fn max_turn_rate(&self) -> f64 {
        self.omega.max_abs()
    }

    fn breakpoints(&self) -> impl Iterator<Item = f64> + '_ {
        self.omega
            .breakpoints()
            .chain(self.speed.breakpoints())
            .chain(self.climb.breakpoints())
    }

    /// (V_h, V̇_h, ḣ, ḧ, ω, h) at local seconds `s`.
    fn channels(&self, s: f64) -> Result<(f64, f64, f64, f64, f64, f64), GeoError> {
        let (v, vd, _, _) = self.speed.at(s);
        let (hd, hdd, _, dh) = self.climb.at(s);
        let (w, _, _, _) = self.omega.at(s);
        let vh2 = v * v - hd * hd;
        let vh = libm::sqrt(vh2.max(0.0));
        if vh.is_nan() || vh < MIN_SPEED_M_S {
            return plan_err(format!(
                "the horizontal speed would fall to {vh} m/s (speed {v} m/s, vertical rate {hd} m/s; minimum {MIN_SPEED_M_S} m/s)"
            ));
        }
        let vhd = (v * vd - hd * hdd) / vh;
        Ok((vh, vhd, hd, hdd, w, self.h0 + dh))
    }

    /// ODE right-hand side: (ṅ, ḋ).
    fn rhs(&self, s: f64, n: Vec3, d: Vec3) -> Result<(Vec3, Vec3), GeoError> {
        let (vh, _, _, _, w, h) = self.channels(s)?;
        let ndot = curvature(n, h, d * vh);
        let ddot = -(n.cross(d) * w) - n * d.dot(ndot);
        Ok((ndot, ddot))
    }

    /// Position, velocity and acceleration at `s` for ODE state (n, d).
    fn kinematics(&self, s: f64, n: Vec3, d: Vec3) -> Result<[Vec3; 3], GeoError> {
        let (vh, vhd, hd, hdd, w, h) = self.channels(s)?;
        let ndot = curvature(n, h, d * vh);
        let ddot = -(n.cross(d) * w) - n * d.dot(ndot);
        Ok([
            point(n, h),
            d * vh + n * hd,
            d * vhd + ddot * vh + n * hdd + ndot * hd,
        ])
    }
}

/// Node offsets (ps from the segment start) for a segment.
fn nodes(seg: &Segment) -> Vec<i64> {
    let end = seg.duration_ps;
    let mut cand: Vec<i64> = seg
        .breakpoints()
        .map(|s| libm::round(s * 1.0e12) as i64)
        .filter(|&p| 0 < p && p < end)
        .collect();
    cand.push(end);
    cand.sort_unstable();
    cand.dedup();
    let mut kept = vec![0_i64];
    for c in cand {
        let last = *kept.last().unwrap_or(&0);
        if c - last >= MIN_NODE_SPACING_PS {
            kept.push(c);
        } else if c == end {
            if kept.len() > 1 {
                kept.pop();
            }
            kept.push(end);
        }
    }
    // Fast turns get closer nodes (at most MAX_NODE_TURN_RAD apart), but never
    // closer than the merge distance.
    let w = seg.max_turn_rate();
    let spacing = if w > 0.0 {
        let s = libm::floor(MAX_NODE_TURN_RAD / w * 1.0e12);
        if s < MAX_NODE_SPACING_PS as f64 {
            (s as i64).max(MIN_NODE_SPACING_PS)
        } else {
            MAX_NODE_SPACING_PS
        }
    } else {
        MAX_NODE_SPACING_PS
    };
    let mut out = vec![0_i64];
    for w in kept.windows(2) {
        let (a, b) = (w[0], w[1]);
        let m = (b - a + spacing - 1) / spacing;
        for k in 1..m {
            out.push(a + ((i128::from(b - a) * i128::from(k)) / i128::from(m)) as i64);
        }
        out.push(b);
    }
    out
}

/// Normalizes n and keeps d a horizontal unit vector.
fn renormalize(n: Vec3, d: Vec3) -> (Vec3, Vec3) {
    let n = n.normalized();
    let d = (d - n * d.dot(n)).normalized();
    (n, d)
}

/// The result of flying one segment.
pub(crate) struct Flown {
    pub pieces: Vec<Piece>,
    pub end: Air,
    pub end_kin: [Vec3; 3],
    pub end_time: TruthTime,
}

/// Flies `seg` from `t0`, ODE state `air`, with exact start kinematics `k0` (the
/// first piece starts from them, so the join is C2 whatever `air` is).
pub(crate) fn fly(
    seg: &Segment,
    t0: TruthTime,
    air: &Air,
    k0: [Vec3; 3],
) -> Result<Flown, GeoError> {
    let end_time =
        TruthTime::from_ps(t0.as_ps().checked_add(seg.duration_ps).ok_or_else(|| {
            GeoError::Trajectory("the maneuver runs past the end of truth time".to_owned())
        })?);
    if seg.duration_ps == 0 {
        // Nothing to do (e.g. a speed change to the current speed).
        return Ok(Flown {
            pieces: Vec::new(),
            end: *air,
            end_kin: k0,
            end_time,
        });
    }
    let offsets = nodes(seg);
    let w_max = seg.max_turn_rate();
    let max_step = if w_max > 0.0 {
        MAX_RK_STEP_S.min(MAX_RK_TURN_RAD / w_max)
    } else {
        MAX_RK_STEP_S
    };
    let (mut n, mut d) = (air.n, air.d);
    let mut prev = k0;
    let mut pieces = Vec::with_capacity(offsets.len());
    for w in offsets.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (sa, sb) = (a as f64 / 1.0e12, b as f64 / 1.0e12);
        let steps = libm::ceil((sb - sa) / max_step).max(1.0);
        let dt = (sb - sa) / steps;
        for k in 0..steps as u32 {
            let s = sa + dt * f64::from(k);
            let (n1, d1) = seg.rhs(s, n, d)?;
            let (n2, d2) = seg.rhs(s + 0.5 * dt, n + n1 * (0.5 * dt), d + d1 * (0.5 * dt))?;
            let (n3, d3) = seg.rhs(s + 0.5 * dt, n + n2 * (0.5 * dt), d + d2 * (0.5 * dt))?;
            let (n4, d4) = seg.rhs(s + dt, n + n3 * dt, d + d3 * dt)?;
            let sixth = dt / 6.0;
            n = n + (n1 + n2 * 2.0 + n3 * 2.0 + n4) * sixth;
            d = d + (d1 + d2 * 2.0 + d3 * 2.0 + d4) * sixth;
            (n, d) = renormalize(n, d);
        }
        let k1 = seg.kinematics(sb, n, d)?;
        pieces.push(Piece::hermite(
            TruthTime::from_ps(t0.as_ps() + a),
            TruthTime::from_ps(t0.as_ps() + b),
            prev,
            k1,
        )?);
        prev = k1;
    }
    let s_end = seg.duration_ps as f64 / 1.0e12;
    let (vh, _, hd, hdd, w, h) = seg.channels(s_end)?;
    let (v, vd, _, _) = seg.speed.at(s_end);
    let _ = vh;
    Ok(Flown {
        pieces,
        end: Air {
            n,
            d,
            h,
            speed: v,
            speed_dot: vd,
            hdot: hd,
            hddot: hdd,
            omega: w,
        },
        end_kin: prev,
        end_time,
    })
}

/// Plans and flies primitives (shared by construction and append).
pub(crate) struct Planner<'a> {
    pub limits: Limits,
    pub ground: bool,
    pub heights: Heights<'a>,
}

impl Planner<'_> {
    fn ramp_rate(&self, speed: f64) -> f64 {
        STANDARD_GRAVITY_M_S2 * self.limits.roll_rate_rad_s / speed
    }

    pub(crate) fn omega_shape(&self, air: &Air) -> RampShape {
        RampShape {
            x0: air.omega,
            rate: self.ramp_rate(air.speed),
        }
    }

    pub(crate) fn climb_shape(&self, air: &Air) -> ScurveShape {
        ScurveShape {
            x0: air.hdot,
            d0: air.hddot,
            max_d: self.limits.vertical_accel_m_s2,
            jerk: self.limits.vertical_jerk_m_s3,
        }
    }

    pub(crate) fn wings_level(&self, air: &Air) -> Result<Channel, GeoError> {
        ramp(air.omega, 0.0, self.ramp_rate(air.speed))
    }

    pub(crate) fn speed_to(&self, air: &Air, target: Option<f64>) -> Result<Channel, GeoError> {
        let (j, a) = (self.limits.jerk_m_s3, self.limits.accel_m_s2);
        let target = match target {
            Some(v) => {
                if !(v >= MIN_SPEED_M_S && v.is_finite()) {
                    return plan_err(format!(
                        "a speed of {v} m/s is below the minimum {MIN_SPEED_M_S} m/s"
                    ));
                }
                v
            }
            // Hold: bring the acceleration to zero.
            None => air.speed + air.speed_dot * air.speed_dot.abs() / (2.0 * j),
        };
        let steps = scurve(air.speed, air.speed_dot, target, a, j)?;
        Ok(Channel::from_steps(air.speed, air.speed_dot, &steps))
    }

    pub(crate) fn level(&self, air: &Air) -> Result<Channel, GeoError> {
        let steps = scurve(
            air.hdot,
            air.hddot,
            0.0,
            self.limits.vertical_accel_m_s2,
            self.limits.vertical_jerk_m_s3,
        )?;
        Ok(Channel::from_steps(air.hdot, air.hddot, &steps))
    }

    fn turn_peak(&self, air: &Air, rate: TurnRate) -> Result<f64, GeoError> {
        let w = match rate {
            TurnRate::Rate(w) => w,
            TurnRate::Bank(b) => {
                if b.is_nan() || b.abs() >= MAX_BANK_RAD {
                    return plan_err("a turn bank angle must be below 80°");
                }
                STANDARD_GRAVITY_M_S2 * libm::tan(b.abs()) / air.speed
            }
        };
        let bank = libm::atan(air.speed * w / STANDARD_GRAVITY_M_S2);
        if !(w.is_finite() && w > 0.0) || bank >= MAX_BANK_RAD {
            return plan_err(format!(
                "a turn rate of {w} rad/s at {} m/s needs a bank of {}° (limit 80°)",
                air.speed,
                bank.to_degrees()
            ));
        }
        Ok(w.abs())
    }

    /// Plans and flies one primitive from (`t`, `air`, `k`).
    pub(crate) fn primitive(
        &self,
        prim: &Primitive,
        t: TruthTime,
        air: &Air,
        k: [Vec3; 3],
    ) -> Result<Flown, GeoError> {
        let secs = |d: TruthDuration| -> Result<f64, GeoError> {
            if d.as_ps() <= 0 {
                return plan_err("a maneuver duration must be positive");
            }
            Ok(d.as_ps() as f64 / 1.0e12)
        };
        let seg = match *prim {
            Primitive::Straight { duration } => Segment::new(
                self.wings_level(air)?,
                self.speed_to(air, None)?,
                self.level(air)?,
                air.h,
                secs(duration)?,
            )?,
            Primitive::Speed { speed_m_s } => Segment::new(
                self.wings_level(air)?,
                self.speed_to(air, Some(speed_m_s))?,
                self.level(air)?,
                air.h,
                0.0,
            )?,
            Primitive::Climb { rate, end } => {
                if self.ground {
                    return plan_err(
                        "a ground platform cannot climb (its height follows the ground)",
                    );
                }
                let c = match rate {
                    ClimbRate::Rate(c) => c,
                    ClimbRate::Angle(g) => air.speed * libm::sin(g),
                };
                if !(c.is_finite() && c != 0.0) {
                    return plan_err("a climb rate must be non-zero");
                }
                let climb = match end {
                    ClimbEnd::Duration(d) => {
                        hold_for(&self.climb_shape(air), c, secs(d)?, "climb")?
                    }
                    ClimbEnd::Altitude(alt) => {
                        let (lat, lon) = air.lat_lon();
                        let target = self.heights.ellipsoid_height_m(lat, lon, alt)?;
                        move_by(&self.climb_shape(air), target - air.h, c.abs())?
                    }
                };
                Segment::new(
                    self.wings_level(air)?,
                    self.speed_to(air, None)?,
                    climb,
                    air.h,
                    0.0,
                )?
            }
            Primitive::Turn {
                rate,
                direction,
                end,
            } => {
                let peak = self.turn_peak(air, rate)?;
                let need_dir = || {
                    direction.ok_or_else(|| {
                        GeoError::Trajectory(
                            "a turn by duration or by angle needs a direction".to_owned(),
                        )
                    })
                };
                let omega = match end {
                    TurnEnd::Duration(d) => hold_for(
                        &self.omega_shape(air),
                        need_dir()?.sign() * peak,
                        secs(d)?,
                        "turn",
                    )?,
                    TurnEnd::By(angle) => move_by(
                        &self.omega_shape(air),
                        need_dir()?.sign() * angle.abs(),
                        peak,
                    )?,
                    TurnEnd::Heading(target) => {
                        return self.turn_to_heading(air, t, k, peak, target, direction);
                    }
                };
                Segment::new(
                    omega,
                    self.speed_to(air, None)?,
                    self.level(air)?,
                    air.h,
                    0.0,
                )?
            }
        };
        fly(&seg, t, air, k)
    }

    /// Turn to a track angle: shooting on the turn angle, because the track also
    /// changes as the platform moves over the curved Earth.
    fn turn_to_heading(
        &self,
        air: &Air,
        t: TruthTime,
        k: [Vec3; 3],
        peak: f64,
        target: f64,
        direction: Option<TurnDirection>,
    ) -> Result<Flown, GeoError> {
        let tau = 2.0 * core::f64::consts::PI;
        let diff = wrap_pi(target - air.track()?);
        let mut delta = match direction {
            None => diff,
            Some(TurnDirection::Right) => {
                if diff < 0.0 {
                    diff + tau
                } else {
                    diff
                }
            }
            Some(TurnDirection::Left) => {
                if diff > 0.0 {
                    diff - tau
                } else {
                    diff
                }
            }
        };
        let speed = self.speed_to(air, None)?;
        let level = self.level(air)?;
        for _ in 0..12 {
            let omega = move_by(&self.omega_shape(air), delta, peak)?;
            let seg = Segment::new(omega, speed.clone(), level.clone(), air.h, 0.0)?;
            let flown = fly(&seg, t, air, k)?;
            let err = wrap_pi(target - flown.end.track()?);
            if err.abs() <= 1.0e-12 {
                return Ok(flown);
            }
            delta += err;
        }
        plan_err(format!(
            "a turn to heading {}° did not converge",
            target.to_degrees()
        ))
    }
}

impl Flight {
    fn empty(attitude: AttitudeMode, limits: Limits) -> Result<Self, GeoError> {
        limits.check()?;
        if matches!(attitude, AttitudeMode::Fixed(_)) {
            return plan_err("a fixed attitude is only for static platforms");
        }
        Ok(Self {
            path: Path::default(),
            attitude,
            limits,
        })
    }

    fn planner<'a>(&self, heights: Heights<'a>) -> Planner<'a> {
        Planner {
            limits: self.limits,
            ground: self.attitude == AttitudeMode::Ground,
            heights,
        }
    }

    /// The starting state for maneuvers or a route.
    fn start_state(start: &Start, track: f64) -> Result<(Air, [Vec3; 3]), GeoError> {
        let g = start.position;
        if !(g.lat_rad.is_finite()
            && g.lon_rad.is_finite()
            && g.h_m.is_finite()
            && track.is_finite())
        {
            return Err(GeoError::NonFinite("trajectory start"));
        }
        if !(start.speed_m_s >= MIN_SPEED_M_S && start.speed_m_s.is_finite()) {
            return plan_err(format!(
                "the start speed must be at least {MIN_SPEED_M_S} m/s"
            ));
        }
        let n = normal(g.lat_rad, g.lon_rad);
        let air = Air {
            n,
            d: direction(n, track)?,
            h: g.h_m,
            speed: start.speed_m_s,
            speed_dot: 0.0,
            hdot: 0.0,
            hddot: 0.0,
            omega: 0.0,
        };
        let seg = Segment::new(
            Channel::constant(0.0),
            Channel::constant(air.speed),
            Channel::constant(0.0),
            air.h,
            0.0,
        )?;
        let k = seg.kinematics(0.0, air.n, air.d)?;
        Ok((air, k))
    }

    fn push_all(&mut self, pieces: Vec<Piece>) {
        for p in pieces {
            self.path.push(p);
        }
    }

    fn run(
        &mut self,
        t: TruthTime,
        air: Air,
        k: [Vec3; 3],
        primitives: &[Primitive],
        heights: Heights<'_>,
    ) -> Result<(), GeoError> {
        let planner = self.planner(heights);
        let (mut t, mut air, mut k) = (t, air, k);
        for (i, prim) in primitives.iter().enumerate() {
            let flown = planner.primitive(prim, t, &air, k).map_err(|e| match e {
                GeoError::Trajectory(m) => GeoError::Trajectory(format!("maneuver {}: {m}", i + 1)),
                other => other,
            })?;
            self.push_all(flown.pieces);
            (t, air, k) = (flown.end_time, flown.end, flown.end_kin);
        }
        Ok(())
    }

    /// Maneuver primitives flown in order from a straight-and-level start.
    pub fn maneuvers(
        start: Start,
        limits: Limits,
        attitude: AttitudeMode,
        primitives: &[Primitive],
        heights: Heights<'_>,
    ) -> Result<Self, GeoError> {
        let track = start
            .track_rad
            .ok_or_else(|| GeoError::Trajectory("maneuvers need an initial heading".to_owned()))?;
        if primitives.is_empty() {
            return plan_err("maneuvers need at least one primitive");
        }
        let mut f = Self::empty(attitude, limits)?;
        let (air, k) = Self::start_state(&start, track)?;
        f.run(start.time, air, k, primitives, heights)?;
        Ok(f)
    }

    /// A waypoint route flown over each waypoint in order (contracts §13.6).
    pub fn route(
        start: Start,
        limits: Limits,
        attitude: AttitudeMode,
        waypoints: &[Waypoint],
        heights: Heights<'_>,
    ) -> Result<Self, GeoError> {
        let Some(first) = waypoints.first() else {
            return plan_err("a route needs at least one waypoint");
        };
        let mut f = Self::empty(attitude, limits)?;
        let track = match start.track_rad {
            Some(t) => t,
            None => {
                let n = normal(start.position.lat_rad, start.position.lon_rad);
                let p0 = point(n, start.position.h_m);
                let pw = point(normal(first.lat_rad, first.lon_rad), start.position.h_m);
                let rel = pw - p0;
                let rel = rel - n * rel.dot(n);
                if rel.norm() < 1.0 {
                    return plan_err("the first waypoint is at the start position; give a heading");
                }
                track_angle(n, rel.normalized())?
            }
        };
        let (air, k) = Self::start_state(&start, track)?;
        let planner = f.planner(heights);
        let pieces = crate::trajectory::route::fly_route(&planner, start.time, air, k, waypoints)?;
        f.push_all(pieces);
        Ok(f)
    }

    /// A recorded track smoothed to C2 (contracts §13.7): a cubic smoothing spline per
    /// ECEF axis with time constant `smoothing_s` (0 interpolates the points).
    pub fn track(
        points: &[TrackPoint],
        smoothing_s: f64,
        limits: Limits,
        attitude: AttitudeMode,
    ) -> Result<Self, GeoError> {
        let mut f = Self::empty(attitude, limits)?;
        let Some(first) = points.first() else {
            return plan_err("a track needs at least two points");
        };
        let origin = first.position.to_ecef().vec();
        let times: Vec<TruthTime> = points.iter().map(|p| p.time).collect();
        let rel: Vec<Vec3> = points
            .iter()
            .map(|p| p.position.to_ecef().vec() - origin)
            .collect();
        let spline = CubicSpline3::fit(&times, &rel, smoothing_s)
            .map_err(|e| GeoError::Trajectory(format!("track: {e}")))?;
        for (i, c) in spline.coeffs.iter().enumerate() {
            f.path.push(Piece::cubic(
                spline.knots[i],
                spline.knots[i + 1],
                origin,
                *c,
            ));
        }
        Ok(f)
    }

    /// Adds primitives from the state at `t` (runtime maneuvers). Everything planned
    /// after `t` is replaced; every instant before `t` keeps its exact bits. The new
    /// motion starts from the exact position, velocity and acceleration at `t`, so the
    /// join is C2. On error the flight is unchanged.
    pub fn append(
        &mut self,
        t: TruthTime,
        primitives: &[Primitive],
        heights: Heights<'_>,
    ) -> Result<(), GeoError> {
        let (air, k) = self.state_for_append(t)?;
        let mut next = self.clone();
        next.path.truncate_at(t);
        next.run(t, air, k, primitives, heights)?;
        *self = next;
        Ok(())
    }

    /// Replaces the plan after `t` with a waypoint route from the state at `t` ("go
    /// to"). Same guarantees as [`append`](Self::append).
    pub fn append_route(
        &mut self,
        t: TruthTime,
        waypoints: &[Waypoint],
        heights: Heights<'_>,
    ) -> Result<(), GeoError> {
        let (air, k) = self.state_for_append(t)?;
        let mut next = self.clone();
        next.path.truncate_at(t);
        let planner = next.planner(heights);
        let pieces = crate::trajectory::route::fly_route(&planner, t, air, k, waypoints)?;
        next.push_all(pieces);
        *self = next;
        Ok(())
    }

    /// Extends the flight with straight-and-level flight until at least `t` (no-op if
    /// it already reaches `t`).
    pub fn extend_to(&mut self, t: TruthTime, heights: Heights<'_>) -> Result<(), GeoError> {
        let end = self.path.end().unwrap_or(TruthTime::MIN);
        if t <= end {
            return Ok(());
        }
        let gap = t.as_ps().checked_sub(end.as_ps()).ok_or_else(|| {
            GeoError::Trajectory("the extension is longer than the truth-time range".to_owned())
        })?;
        let duration = TruthDuration::from_ps(gap.max(MIN_SEGMENT_PS));
        self.append(end, &[Primitive::Straight { duration }], heights)
    }

    fn state_for_append(&self, t: TruthTime) -> Result<(Air, [Vec3; 3]), GeoError> {
        let kin = self.path.eval(AcqInstant::from_truth_time(t))?;
        let air = Air::from_kinematics(kin.p, kin.v, kin.a)?;
        Ok((air, [kin.p, kin.v, kin.a]))
    }

    /// The kinematic pieces.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The attitude mode.
    pub fn attitude(&self) -> AttitudeMode {
        self.attitude
    }

    /// Kinematics (with jerk) at `t`.
    pub fn kinematics(&self, t: AcqInstant) -> Result<Kinematics, GeoError> {
        self.path.eval(t)
    }
}

impl Trajectory for Flight {
    fn state(&self, t: AcqInstant) -> Result<State, GeoError> {
        let k = self.path.eval(t)?;
        let (q, w) = derive(&self.attitude, &k)?;
        Ok(State {
            position: Ecef::from_vec(k.p),
            velocity: k.v,
            acceleration: k.a,
            q_ecef_body: q,
            angular_rate_body: w,
        })
    }

    fn span(&self) -> Span {
        Span {
            start: self.path.start(),
            end: self.path.end(),
        }
    }

    fn knots(&self, from: TruthTime, to: TruthTime) -> Vec<TruthTime> {
        self.path.knots(from, to)
    }
}
