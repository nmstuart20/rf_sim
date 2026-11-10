//! Waypoint routes (contracts §13.6).
//!
//! Each waypoint is flown over. A leg from the current state to the next waypoint is
//! one segment whose channels run together from the leg start:
//!
//! - turn by Δ at the route bank (rate g₀ tan φ / V at the leg's start speed);
//! - climb or descend at the route climb rate to the waypoint altitude (converted at
//!   the waypoint), level by the waypoint;
//! - change to the waypoint speed, if given;
//! - straight flight for the rest of the leg.
//!
//! Two-variable shooting finds Δ and the leg duration T. From the end state, the
//! along-track miss `e_a` and the cross-track miss `e_c` (right positive) update
//! `T += e_a / V_h` and `Δ += e_c / L`, with L the distance flown after the turn.
//! It stops when both are below [`ARRIVAL_TOLERANCE_M`], and fails after
//! [`MAX_ITERATIONS`] or when a transition cannot finish before the waypoint.

use iqsim_core::TruthTime;

use crate::GeoError;
use crate::trajectory::STANDARD_GRAVITY_M_S2;
use crate::trajectory::earth::{normal, point, radii};
use crate::trajectory::flight::{Air, Planner, Segment, Waypoint, fly, plan_err};
use crate::trajectory::piece::Piece;
use crate::trajectory::profile::move_by;
use crate::vec::Vec3;

/// Largest miss at a waypoint when shooting stops, metres.
pub const ARRIVAL_TOLERANCE_M: f64 = 1.0e-4;
/// Shooting iterations per leg.
pub const MAX_ITERATIONS: u32 = 40;

pub(crate) fn fly_route(
    planner: &Planner<'_>,
    t: TruthTime,
    air: Air,
    k: [Vec3; 3],
    waypoints: &[Waypoint],
) -> Result<Vec<Piece>, GeoError> {
    let (mut t, mut air, mut k) = (t, air, k);
    let mut pieces = Vec::new();
    for (i, wp) in waypoints.iter().enumerate() {
        let flown = leg(planner, t, &air, k, wp).map_err(|e| match e {
            GeoError::Trajectory(m) => GeoError::Trajectory(format!("waypoint {}: {m}", i + 1)),
            other => other,
        })?;
        pieces.extend(flown.pieces);
        (t, air, k) = (flown.end_time, flown.end, flown.end_kin);
    }
    Ok(pieces)
}

/// Angle between normals beyond which the miss is measured along the surface.
const NEAR_RAD: f64 = 1.0e-3;

/// Along-track and cross-track miss (right positive) from `state` to the waypoint.
/// Near the waypoint it is the chord projected on the local horizontal. Farther away
/// it is the surface arc `angle(n, n_w)·(N + h)` along the initial direction to the
/// waypoint, so a far (even antipodal-side) waypoint is never mistaken for an arrival.
fn miss(state: &Air, n_w: Vec3, target: Vec3) -> Result<(f64, f64), GeoError> {
    let right = state.d.cross(state.n);
    let angle = libm::atan2(state.n.cross(n_w).norm(), state.n.dot(n_w));
    let dir = if angle < NEAR_RAD {
        let rel = target - state.position();
        rel - state.n * rel.dot(state.n)
    } else {
        let u = n_w - state.n * n_w.dot(state.n);
        if u.norm() < 1.0e-9 {
            return plan_err(
                "the waypoint is antipodal to the current position; add a waypoint in between",
            );
        }
        let (big_n, _) = radii(state.n);
        u.normalized() * (angle * (big_n + state.h))
    };
    Ok((dir.dot(state.d), dir.dot(right)))
}

fn leg(
    planner: &Planner<'_>,
    t: TruthTime,
    air: &Air,
    k: [Vec3; 3],
    wp: &Waypoint,
) -> Result<crate::trajectory::flight::Flown, GeoError> {
    let h_w = planner
        .heights
        .ellipsoid_height_m(wp.lat_rad, wp.lon_rad, wp.altitude)?;
    let target = point(normal(wp.lat_rad, wp.lon_rad), h_w);
    let dh = h_w - air.h;
    // Ground platforms change height the same way, between the waypoints' AGL
    // heights over the placeholder ground.
    let climb = if dh.abs() > 1.0e-9 {
        move_by(&planner.climb_shape(air), dh, planner.limits.climb_rate_m_s)?
    } else {
        planner.level(air)?
    };
    let speed = planner.speed_to(air, wp.speed_m_s)?;
    let omega_max = STANDARD_GRAVITY_M_S2 * libm::tan(planner.limits.bank_rad) / air.speed;

    let n_w = normal(wp.lat_rad, wp.lon_rad);
    let (along0, cross0) = miss(air, n_w, target)?;
    let mut delta = libm::atan2(cross0, along0);
    let mut duration = libm::hypot(along0, cross0) / air.speed;
    for _ in 0..MAX_ITERATIONS {
        let omega = move_by(&planner.omega_shape(air), delta, omega_max)?;
        let turn_end = omega.duration();
        let min = turn_end.max(climb.duration()).max(speed.duration());
        let used = duration.max(min);
        let seg = Segment::new(omega, speed.clone(), climb.clone(), air.h, used)?;
        let flown = fly(&seg, t, air, k)?;
        let e = &flown.end;
        let (along, cross) = miss(e, n_w, target)?;
        let distance = target.distance(e.position());
        if along.abs() <= ARRIVAL_TOLERANCE_M
            && cross.abs() <= ARRIVAL_TOLERANCE_M
            && distance <= 10.0 * ARRIVAL_TOLERANCE_M
        {
            return Ok(flown);
        }
        let vh = libm::sqrt((e.speed * e.speed - e.hdot * e.hdot).max(1.0));
        if used <= min && along < -ARRIVAL_TOLERANCE_M && duration < min {
            return plan_err(format!(
                "the waypoint is too close: the turn, climb and speed change need {min} s"
            ));
        }
        duration = used + along / vh;
        let lever = (vh * (used - turn_end)).max(vh);
        delta += cross / lever;
    }
    plan_err("the route leg did not converge (is the waypoint inside the turn circle?)")
}
