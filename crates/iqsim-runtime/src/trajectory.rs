//! Scenario trajectories → `iqsim-geo` trajectories (contracts §11.4, §13).
//!
//! [`build`] turns a platform's `[platforms.trajectory]` table into a
//! [`Trajectory`]: a [`Static`] platform or a [`Flight`]. Scenario altitudes resolve
//! through the run's [`AltitudeModel`] (geoid, atmosphere), and AGL heights over the
//! ellipsoid, as everywhere until terrain exists. [`build_until`] also extends a
//! flight whose plan ends before the run does, straight and level
//! ([`Flight::extend_to`]), so every platform is defined for the whole run.

use std::sync::Arc;

use iqsim_core::TruthTime;
use iqsim_geo::altitude::AltitudeModel;
use iqsim_geo::trajectory::{
    AttitudeMode, ClimbEnd, ClimbRate, FlatGround, Heights, Limits, Primitive, Start, TrackPoint,
    TurnDirection, TurnEnd, TurnRate, Waypoint,
};
use iqsim_geo::{Flight, Geodetic, Quat, Static, Trajectory, Vec3};
use iqsim_scenario::{DirectionSpec, ManeuverKind, ManeuverSpec, Platform, TrajectoryKind};
use iqsim_units::{Altitude, Position};

use crate::plan::ResolveError;

fn geo(owner: &str) -> impl Fn(iqsim_geo::GeoError) -> ResolveError + '_ {
    move |e| ResolveError(format!("{owner}: {e}"))
}

fn bad(owner: &str, m: &str) -> ResolveError {
    ResolveError(format!("{owner}: {m}"))
}

fn geodetic(p: &Position, model: &AltitudeModel, owner: &str) -> Result<Geodetic, ResolveError> {
    let (lat, lon) = (
        p.latitude().degrees().to_radians(),
        p.longitude().degrees().to_radians(),
    );
    let h = model
        .ellipsoid_height_m(lat, lon, p.altitude())
        .map_err(geo(owner))?;
    Ok(Geodetic::new(lat, lon, h))
}

fn limits(p: &Platform) -> Limits {
    let l = &p.trajectory.limits;
    let d = Limits::default();
    Limits {
        bank_rad: l
            .bank
            .as_ref()
            .map_or(d.bank_rad, |x| x.get_ref().radians()),
        roll_rate_rad_s: l
            .roll_rate
            .as_ref()
            .map_or(d.roll_rate_rad_s, |x| x.get_ref().radians_per_second()),
        accel_m_s2: l
            .accel
            .as_ref()
            .map_or(d.accel_m_s2, |x| x.get_ref().meters_per_second_squared()),
        jerk_m_s3: l
            .jerk
            .as_ref()
            .map_or(d.jerk_m_s3, |x| x.get_ref().meters_per_second_cubed()),
        vertical_accel_m_s2: l
            .vertical_accel
            .as_ref()
            .map_or(d.vertical_accel_m_s2, |x| {
                x.get_ref().meters_per_second_squared()
            }),
        vertical_jerk_m_s3: l.vertical_jerk.as_ref().map_or(d.vertical_jerk_m_s3, |x| {
            x.get_ref().meters_per_second_cubed()
        }),
        climb_rate_m_s: l
            .climb_rate
            .as_ref()
            .map_or(d.climb_rate_m_s, |x| x.get_ref().meters_per_second()),
    }
}

fn attitude(p: &Platform) -> AttitudeMode {
    let t = &p.trajectory;
    if t.kind == TrajectoryKind::Ground {
        return AttitudeMode::Ground;
    }
    let wind = t.wind.as_ref().map_or(Vec3::ZERO, |w| {
        let w = w.get_ref();
        let (s, c) = libm::sincos(w.from.radians());
        let v = w.speed.meters_per_second();
        // The air moves away from `from`.
        Vec3::new(-(v * c), -(v * s), 0.0)
    });
    AttitudeMode::Derived {
        angle_of_attack_rad: t
            .angle_of_attack
            .as_ref()
            .map_or(0.0, |a| a.get_ref().radians()),
        wind_ned_m_s: wind,
    }
}

fn primitive(m: &ManeuverSpec, owner: &str) -> Result<Primitive, ResolveError> {
    let dur = |m: &ManeuverSpec| m.duration.as_ref().map(|d| d.get_ref().get());
    let missing = || bad(owner, "the maneuver is incomplete (run `sim lint`)");
    Ok(match m.kind {
        ManeuverKind::Straight => Primitive::Straight {
            duration: dur(m).ok_or_else(missing)?,
        },
        ManeuverKind::Speed => Primitive::Speed {
            speed_m_s: m
                .speed
                .as_ref()
                .ok_or_else(missing)?
                .get_ref()
                .meters_per_second(),
        },
        ManeuverKind::Turn => {
            let rate = match (&m.rate, &m.bank) {
                (Some(r), None) => TurnRate::Rate(r.get_ref().radians_per_second()),
                (None, Some(b)) => TurnRate::Bank(b.get_ref().radians()),
                _ => return Err(missing()),
            };
            let end = match (dur(m), &m.heading, &m.by) {
                (Some(d), None, None) => TurnEnd::Duration(d),
                (None, Some(h), None) => TurnEnd::Heading(h.get_ref().radians()),
                (None, None, Some(b)) => TurnEnd::By(b.get_ref().radians()),
                _ => return Err(missing()),
            };
            Primitive::Turn {
                rate,
                direction: m.direction.map(|d| match d {
                    DirectionSpec::Left => TurnDirection::Left,
                    DirectionSpec::Right => TurnDirection::Right,
                }),
                end,
            }
        }
        ManeuverKind::Climb => {
            let rate = match (&m.vertical_rate, &m.angle) {
                (Some(r), None) => ClimbRate::Rate(r.get_ref().meters_per_second()),
                (None, Some(a)) => ClimbRate::Angle(a.get_ref().radians()),
                _ => return Err(missing()),
            };
            let end = match (&m.altitude, dur(m)) {
                (Some(a), None) => ClimbEnd::Altitude(*a.get_ref()),
                (None, Some(d)) => ClimbEnd::Duration(d),
                _ => return Err(missing()),
            };
            Primitive::Climb { rate, end }
        }
    })
}

fn waypoint(lat: f64, lon: f64, alt: Altitude, speed: Option<f64>) -> Waypoint {
    Waypoint {
        lat_rad: lat.to_radians(),
        lon_rad: lon.to_radians(),
        altitude: alt,
        speed_m_s: speed,
    }
}

/// The platform's trajectory as planned (module docs). Expects a lint-clean scenario.
pub fn build(p: &Platform, model: &AltitudeModel) -> Result<Arc<dyn Trajectory>, ResolveError> {
    build_inner(p, model, None)
}

/// The platform's trajectory, extended straight and level to at least `until` if its
/// plan ends earlier.
pub fn build_until(
    p: &Platform,
    model: &AltitudeModel,
    until: TruthTime,
) -> Result<Arc<dyn Trajectory>, ResolveError> {
    build_inner(p, model, Some(until))
}

fn build_inner(
    p: &Platform,
    model: &AltitudeModel,
    until: Option<TruthTime>,
) -> Result<Arc<dyn Trajectory>, ResolveError> {
    let owner = format!("platform \"{}\"", p.name.as_str());
    let owner = owner.as_str();
    let t = &p.trajectory;
    let heights = Heights {
        model,
        ground: &FlatGround::ELLIPSOID,
    };
    let position = || {
        p.position
            .as_ref()
            .ok_or_else(|| bad(owner, "needs a position"))
            .and_then(|pos| geodetic(pos, model, owner))
    };
    let start = || -> Result<Start, ResolveError> {
        Ok(Start {
            time: TruthTime::from_ps(t.start.as_ref().map_or(0, |s| s.get_ref().get().as_ps())),
            position: position()?,
            speed_m_s: t
                .speed
                .as_ref()
                .ok_or_else(|| bad(owner, "needs a speed"))?
                .get_ref()
                .meters_per_second(),
            track_rad: t.heading.as_ref().map(|h| h.get_ref().radians()),
        })
    };
    let maneuvers = || -> Result<Vec<Primitive>, ResolveError> {
        t.maneuvers
            .iter()
            .enumerate()
            .map(|(i, m)| primitive(m.get_ref(), &format!("{owner}, maneuver {}", i + 1)))
            .collect()
    };
    let waypoints = || -> Vec<Waypoint> {
        t.waypoints
            .iter()
            .map(|w| {
                let w = w.get_ref();
                waypoint(
                    w.lat.degrees(),
                    w.lon.degrees(),
                    *w.alt.get_ref(),
                    w.speed.as_ref().map(|s| s.get_ref().meters_per_second()),
                )
            })
            .collect()
    };
    let mode = attitude(p);
    let flight = match t.kind {
        TrajectoryKind::Static => {
            let a = t.attitude.as_ref().map(|a| a.get_ref());
            let angle = |x: Option<iqsim_units::Angle>| x.map_or(0.0, |a| a.radians());
            let q = a.map_or(Quat::IDENTITY, |a| {
                Quat::from_euler_zyx(angle(a.heading), angle(a.pitch), angle(a.roll))
            });
            return Ok(Arc::new(Static::new(position()?, q)));
        }
        TrajectoryKind::Maneuvers => {
            Flight::maneuvers(start()?, limits(p), mode, &maneuvers()?, heights)
        }
        TrajectoryKind::Route => Flight::route(start()?, limits(p), mode, &waypoints(), heights),
        TrajectoryKind::Ground if !t.waypoints.is_empty() => {
            Flight::route(start()?, limits(p), mode, &waypoints(), heights)
        }
        TrajectoryKind::Ground => {
            Flight::maneuvers(start()?, limits(p), mode, &maneuvers()?, heights)
        }
        TrajectoryKind::Track => {
            let points = t
                .points
                .iter()
                .map(|q| {
                    let q = q.get_ref();
                    let (lat, lon) = (q.lat.degrees().to_radians(), q.lon.degrees().to_radians());
                    let h = model
                        .ellipsoid_height_m(lat, lon, q.alt)
                        .map_err(geo(owner))?;
                    Ok(TrackPoint {
                        time: TruthTime::from_ps(q.time.get_ref().get().as_ps()),
                        position: Geodetic::new(lat, lon, h),
                    })
                })
                .collect::<Result<Vec<_>, ResolveError>>()?;
            let tau = t
                .smoothing
                .as_ref()
                .map_or(0.0, |s| s.get_ref().get().as_ps() as f64 / 1.0e12);
            Flight::track(&points, tau, limits(p), mode)
        }
    };
    let mut flight = flight.map_err(geo(owner))?;
    if let Some(t) = until {
        flight.extend_to(t, heights).map_err(geo(owner))?;
    }
    Ok(Arc::new(flight))
}
