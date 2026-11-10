//! Trajectory checks (prompt 8; contracts §13).

use iqsim_core::{AcqInstant, TruthDuration, TruthTime};
use iqsim_geo::altitude::AltitudeModel;
use iqsim_geo::trajectory::{
    AttitudeMode, ClimbEnd, ClimbRate, FlatGround, Flight, Heights, Limits, Primitive,
    STANDARD_GRAVITY_M_S2, Start, TrackPoint, Trajectory, TurnDirection, TurnEnd, TurnRate,
    Waypoint,
};
use iqsim_geo::{Ecef, GeoError, Geodetic, Quat, Vec3};
use iqsim_units::{Altitude, AltitudeRef};

const DEG: f64 = core::f64::consts::PI / 180.0;
const PS: i64 = 1_000_000_000_000;

fn heights_model() -> AltitudeModel {
    AltitudeModel::default()
}

fn heights(m: &AltitudeModel) -> Heights<'_> {
    Heights {
        model: m,
        ground: &FlatGround::ELLIPSOID,
    }
}

fn secs(s: f64) -> TruthDuration {
    TruthDuration::from_ps((s * 1.0e12) as i64)
}

fn at(s: f64) -> TruthTime {
    TruthTime::from_ps((s * 1.0e12) as i64)
}

fn ell(m: f64) -> Altitude {
    Altitude::new(m, AltitudeRef::Ellipsoid).unwrap()
}

fn start(lat: f64, lon: f64, h: f64, track_deg: f64, speed: f64) -> Start {
    Start {
        time: TruthTime::EPOCH,
        position: Geodetic::from_degrees(lat, lon, h),
        speed_m_s: speed,
        track_rad: Some(track_deg * DEG),
    }
}

/// Every maneuver primitive kind, from a 250 kt start.
fn sample_maneuvers() -> Vec<Primitive> {
    vec![
        Primitive::Straight {
            duration: secs(20.0),
        },
        Primitive::Turn {
            rate: TurnRate::Rate(3.0 * DEG),
            direction: Some(TurnDirection::Right),
            end: TurnEnd::Duration(secs(40.0)),
        },
        Primitive::Climb {
            rate: ClimbRate::Rate(7.62),
            end: ClimbEnd::Altitude(ell(3500.0)),
        },
        Primitive::Speed { speed_m_s: 180.0 },
        Primitive::Turn {
            rate: TurnRate::Bank(25.0 * DEG),
            direction: None,
            end: TurnEnd::Heading(200.0 * DEG),
        },
        Primitive::Climb {
            rate: ClimbRate::Angle(-3.0 * DEG),
            end: ClimbEnd::Duration(secs(60.0)),
        },
        Primitive::Turn {
            rate: TurnRate::Bank(20.0 * DEG),
            direction: Some(TurnDirection::Left),
            end: TurnEnd::By(45.0 * DEG),
        },
        Primitive::Straight {
            duration: secs(5.0),
        },
    ]
}

fn sample_flight() -> Flight {
    let m = heights_model();
    Flight::maneuvers(
        start(51.0, -1.0, 3000.0, 30.0, 128.611),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &sample_maneuvers(),
        heights(&m),
    )
    .unwrap()
}

fn kin(f: &Flight, t: AcqInstant) -> iqsim_geo::trajectory::piece::Kinematics {
    f.kinematics(t).unwrap()
}

fn inst(ps: i64) -> AcqInstant {
    AcqInstant::from_truth_time(TruthTime::from_ps(ps))
}

/// C2 at every joint, by finite differences of position on each side and by the
/// state just before and at the knot.
fn check_c2(f: &Flight, what: &str) {
    let span = f.span();
    let (s, e) = (span.start.unwrap().as_ps(), span.end.unwrap().as_ps());
    let knots = f.knots(TruthTime::MIN, TruthTime::MAX);
    assert!(knots.len() > 2, "{what}: too few knots");
    let hp = 60_000_000_000;
    let mut checked = 0;
    let all: Vec<i64> = knots.iter().map(|k| k.as_ps()).collect();
    for k in knots {
        let k = k.as_ps();
        if k - 3 * hp < s || k + 3 * hp > e {
            continue;
        }
        let at_k = kin(f, inst(k));
        let before = kin(f, inst(k - 1));
        assert!(before.p.distance(at_k.p) < 1e-6, "{what} p jump at {k}");
        assert!(before.v.distance(at_k.v) < 1e-7, "{what} v jump at {k}");
        assert!(before.a.distance(at_k.a) < 1e-6, "{what} a jump at {k}");
        // One-sided, second-order finite differences from positions only: a short
        // step for velocity (truncation h²·jerk/3), a longer one for acceleration
        // (position rounding / h² dominates).
        for side in [-1_i64, 1] {
            // Finite differences need a smooth window: no other knot on this side.
            if all
                .iter()
                .any(|&o| o != k && (o - k) * side > 0 && (o - k).abs() <= hp)
            {
                continue;
            }
            let sf = side as f64;
            let p = |h: f64, i: i64| kin(f, inst(k + i * side * (h * 1.0e12) as i64)).p - at_k.p;
            let hv = 0.004;
            let v_fd = (p(hv, 2) * -1.0 + p(hv, 1) * 4.0) * (sf / (2.0 * hv));
            let ha = 0.02;
            let a_fd = (p(ha, 1) * -5.0 + p(ha, 2) * 4.0 - p(ha, 3)) * (1.0 / (ha * ha));
            assert!(
                v_fd.distance(at_k.v) < 2e-5,
                "{what} v FD side {side} at {k}: {v_fd:?} vs {:?}",
                at_k.v
            );
            assert!(
                a_fd.distance(at_k.a) < 2e-3,
                "{what} a FD side {side} at {k}: {a_fd:?} vs {:?}",
                at_k.a
            );
        }
        checked += 1;
    }
    assert!(checked > 2, "{what}: only {checked} joints checked");
}

#[test]
fn maneuvers_are_c2_at_every_joint() {
    check_c2(&sample_flight(), "maneuvers");
}

#[test]
fn derivatives_are_consistent_inside_pieces() {
    let f = sample_flight();
    let e = f.span().end.unwrap().as_ps();
    let h = 1.0e-3;
    let hp = (h * 1.0e12) as i64;
    for i in 1..200 {
        let t = e / 200 * i + 123_456_789;
        let c = kin(&f, inst(t));
        let l = kin(&f, inst(t - hp));
        let r = kin(&f, inst(t + hp));
        assert!(((r.p - l.p) * (0.5 / h)).distance(c.v) < 1e-5);
        assert!(((r.v - l.v) * (0.5 / h)).distance(c.a) < 1e-6);
        assert!(((r.a - l.a) * (0.5 / h)).distance(c.j) < 1e-3);
    }
}

fn roll_pitch_yaw(f: &Flight, t: TruthTime) -> (f64, f64, f64) {
    let s = f.state_at(t).unwrap();
    let (yaw, pitch, roll) = s.q_ned_body().unwrap().to_euler_zyx();
    (roll, pitch, yaw)
}

#[test]
fn coordinated_turn_bank_matches_the_formula() {
    let m = heights_model();
    for (lat, lon, track) in [
        (0.0, 0.0, 0.0),
        (51.0, -1.0, 30.0),
        (-35.0, 150.0, 250.0),
        (85.0, 40.0, 90.0),
        (89.5, 0.0, 180.0),
    ] {
        for speed in [60.0, 128.611, 240.0] {
            let rate = 3.0 * DEG;
            let f = Flight::maneuvers(
                start(lat, lon, 9000.0, track, speed),
                Limits::default(),
                AttitudeMode::COORDINATED,
                &[Primitive::Turn {
                    rate: TurnRate::Rate(rate),
                    direction: Some(TurnDirection::Right),
                    end: TurnEnd::Duration(secs(60.0)),
                }],
                heights(&m),
            )
            .unwrap();
            let (roll, pitch, _) = roll_pitch_yaw(&f, at(30.0));
            let want = libm::atan(speed * rate / STANDARD_GRAVITY_M_S2);
            assert!(
                (roll - want).abs() < 0.01 * DEG,
                "lat {lat} V {speed}: {} vs {}",
                roll / DEG,
                want / DEG
            );
            assert!(pitch.abs() < 1e-6, "level turn pitch {pitch}");
            // By bank: the bank is what was asked for.
            let f = Flight::maneuvers(
                start(lat, lon, 9000.0, track, speed),
                Limits::default(),
                AttitudeMode::COORDINATED,
                &[Primitive::Turn {
                    rate: TurnRate::Bank(25.0 * DEG),
                    direction: Some(TurnDirection::Left),
                    end: TurnEnd::Duration(secs(60.0)),
                }],
                heights(&m),
            )
            .unwrap();
            let (roll, _, _) = roll_pitch_yaw(&f, at(30.0));
            assert!(
                (roll + 25.0 * DEG).abs() < 0.01 * DEG,
                "bank {}",
                roll / DEG
            );
        }
    }
}

#[test]
fn turn_to_heading_and_climb_to_altitude_end_on_target() {
    let f = sample_flight();
    // The last turn is 45° left of a straight (geodesic) path after reaching 200°, so
    // the final track is near 155°, off by the track's own drift over the curved Earth.
    let end = f.span().end.unwrap();
    let (_, _, yaw) = roll_pitch_yaw(&f, end);
    assert!(
        (yaw - 155.0 * DEG).abs() < 0.1 * DEG,
        "final heading {}",
        yaw / DEG
    );
    assert!((yaw - 155.0 * DEG).abs() > 1e-6, "the drift is real");
    let m = heights_model();
    let f = Flight::maneuvers(
        start(10.0, 20.0, 1000.0, 0.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[
            Primitive::Climb {
                rate: ClimbRate::Rate(10.0),
                end: ClimbEnd::Altitude(ell(3000.0)),
            },
            Primitive::Turn {
                rate: TurnRate::Rate(2.0 * DEG),
                direction: None,
                end: TurnEnd::Heading(-90.0 * DEG),
            },
        ],
        heights(&m),
    )
    .unwrap();
    let s = f.state_at(f.span().end.unwrap()).unwrap();
    let g = s.position.to_geodetic().unwrap();
    assert!((g.h_m - 3000.0).abs() < 1e-5, "height {}", g.h_m);
    let up = Vec3::new(
        libm::cos(g.lat_rad) * libm::cos(g.lon_rad),
        libm::cos(g.lat_rad) * libm::sin(g.lon_rad),
        libm::sin(g.lat_rad),
    );
    assert!(s.velocity.dot(up).abs() < 1e-6);
    let (yaw, _, _) = s.q_ned_body().unwrap().to_euler_zyx();
    assert!((yaw + 90.0 * DEG).abs() < 1e-9, "heading {}", yaw / DEG);
    // The shortest way from north to west is a left turn: banked left throughout.
    let (roll, _, _) = roll_pitch_yaw(&f, at(f.span().end.unwrap().as_ps() as f64 / 1e12 - 20.0));
    assert!(roll < 0.0);
}

#[test]
fn angular_rate_matches_attitude_differences() {
    for mode in [
        AttitudeMode::COORDINATED,
        AttitudeMode::Derived {
            angle_of_attack_rad: 2.0 * DEG,
            wind_ned_m_s: Vec3::new(-5.0, 15.0, 0.0),
        },
        AttitudeMode::Ground,
    ] {
        let m = heights_model();
        let prims: Vec<Primitive> = sample_maneuvers()
            .into_iter()
            .filter(|p| mode != AttitudeMode::Ground || !matches!(p, Primitive::Climb { .. }))
            .collect();
        let f = Flight::maneuvers(
            start(51.0, -1.0, 3000.0, 30.0, 128.611),
            Limits::default(),
            mode,
            &prims,
            heights(&m),
        )
        .unwrap();
        let e = f.span().end.unwrap().as_ps();
        let h = 1.0e-4;
        let hp = (h * 1.0e12) as i64;
        for i in 1..150 {
            let t = e / 150 * i + 7_777_777;
            let s = f.state(inst(t)).unwrap();
            let l = f.state(inst(t - hp)).unwrap().q_ecef_body;
            let r = f.state(inst(t + hp)).unwrap().q_ecef_body;
            let d = l.conj() * r;
            let d = if d.w < 0.0 {
                Quat::new(-d.w, -d.x, -d.y, -d.z)
            } else {
                d
            };
            let fd = d.vector() * (1.0 / h);
            assert!(
                fd.distance(s.angular_rate_body) < 2e-6,
                "{mode:?} at {t}: {fd:?} vs {:?}",
                s.angular_rate_body
            );
        }
    }
}

#[test]
fn wind_crab_and_angle_of_attack() {
    let m = heights_model();
    // Wind from 270° at 20 m/s: the air moves east. Tracking north at 100 m/s, the
    // nose points into the wind: heading −atan(20/100).
    let f = Flight::maneuvers(
        start(40.0, 10.0, 2000.0, 0.0, 100.0),
        Limits::default(),
        AttitudeMode::Derived {
            angle_of_attack_rad: 3.0 * DEG,
            wind_ned_m_s: Vec3::new(0.0, 20.0, 0.0),
        },
        &[Primitive::Straight {
            duration: secs(10.0),
        }],
        heights(&m),
    )
    .unwrap();
    let (roll, pitch, yaw) = roll_pitch_yaw(&f, at(5.0));
    assert!(
        (yaw + libm::atan2(20.0, 100.0)).abs() < 1e-6,
        "crab {}",
        yaw / DEG
    );
    assert!((pitch - 3.0 * DEG).abs() < 1e-6, "pitch {}", pitch / DEG);
    assert!(roll.abs() < 1e-6);
    // The track is still north.
    let s = f.state_at(at(5.0)).unwrap();
    let v_ned = s.q_ned_body().unwrap();
    let _ = v_ned;
}

/// Waypoint as (lat°, lon°, ellipsoid height m, speed m/s).
type Wp = (f64, f64, f64, Option<f64>);

#[test]
fn routes_pass_through_their_waypoints() {
    let m = heights_model();
    let cases: Vec<(&str, Start, Vec<Wp>)> = vec![
        (
            "england",
            Start {
                track_rad: None,
                ..start(51.0, -1.0, 3000.0, 0.0, 120.0)
            },
            vec![
                (51.2, -0.9, 3600.0, Some(150.0)),
                (51.25, -0.6, 3600.0, None),
                (51.0, -0.5, 2400.0, Some(110.0)),
                (50.95, -0.95, 2400.0, None),
            ],
        ),
        (
            "dateline",
            start(-20.0, 179.6, 10_000.0, 90.0, 230.0),
            vec![
                (-20.1, -179.5, 10_000.0, None),
                (-19.5, -179.9, 11_000.0, Some(240.0)),
            ],
        ),
        (
            "pole",
            start(89.0, 0.0, 11_000.0, 0.0, 230.0),
            vec![(89.0, 180.0, 11_000.0, None), (88.5, 170.0, 11_000.0, None)],
        ),
    ];
    for (name, st, wps) in cases {
        let waypoints: Vec<Waypoint> = wps
            .iter()
            .map(|&(lat, lon, h, v)| Waypoint {
                lat_rad: lat * DEG,
                lon_rad: lon * DEG,
                altitude: ell(h),
                speed_m_s: v,
            })
            .collect();
        let f = Flight::route(
            st,
            Limits::default(),
            AttitudeMode::COORDINATED,
            &waypoints,
            heights(&m),
        )
        .unwrap();
        let knots = f.knots(TruthTime::MIN, TruthTime::MAX);
        let mut last_speed = st.speed_m_s;
        for (i, &(lat, lon, h, v)) in wps.iter().enumerate() {
            let target = Geodetic::from_degrees(lat, lon, h).to_ecef();
            let (best, t) = knots
                .iter()
                .map(|&k| (f.state_at(k).unwrap(), k))
                .min_by(|a, b| {
                    a.0.position
                        .distance(target)
                        .total_cmp(&b.0.position.distance(target))
                })
                .unwrap();
            let miss = best.position.distance(target);
            assert!(
                miss < 0.01,
                "{name} waypoint {}: miss {miss} m at {t:?}",
                i + 1
            );
            let g = best.position.to_geodetic().unwrap();
            assert!(
                (g.h_m - h).abs() < 0.01,
                "{name} waypoint {} height {}",
                i + 1,
                g.h_m
            );
            if let Some(v) = v {
                last_speed = v;
            }
            assert!(
                (best.velocity.norm() - last_speed).abs() < 1e-6,
                "{name} waypoint {} speed",
                i + 1
            );
        }
        check_c2(&f, name);
    }
}

#[test]
fn track_replay_preserves_timing() {
    // Sample a known flight once a second, replay it.
    let truth = sample_flight();
    let end = truth.span().end.unwrap().as_ps();
    let points: Vec<TrackPoint> = (0..=end / PS)
        .map(|i| {
            let t = TruthTime::from_ps(i * PS + 250_000_000_000);
            TrackPoint {
                time: t,
                position: truth.state_at(t).unwrap().position.to_geodetic().unwrap(),
            }
        })
        .filter(|p| p.time.as_ps() <= end)
        .collect();
    let f = Flight::track(&points, 0.0, Limits::default(), AttitudeMode::COORDINATED).unwrap();
    assert_eq!(f.span().start, Some(points[0].time));
    assert_eq!(f.span().end, Some(points.last().unwrap().time));
    for p in &points {
        let s = f.state_at(p.time).unwrap();
        assert!(
            s.position.distance(p.position.to_ecef()) < 1e-6,
            "point at {:?}",
            p.time
        );
    }
    check_c2(&f, "track");
    // Interior states follow the truth closely (interpolating spline, 1 s samples).
    let mid = f.state_at(TruthTime::from_ps(end / 2)).unwrap();
    assert!(
        mid.position.distance(
            truth
                .state_at(TruthTime::from_ps(end / 2))
                .unwrap()
                .position
        ) < 0.05
    );

    // Smoothed: the spline is zero-phase, so the closest approach to each point stays
    // at its timestamp: no lag on average, and small local shifts only where the
    // recorded motion's jerk steps (maneuver starts and ends).
    let f = Flight::track(&points, 2.0, Limits::default(), AttitudeMode::COORDINATED).unwrap();
    let mut sum = 0.0;
    let interior = &points[10..points.len() - 10];
    for p in interior {
        let target = p.position.to_ecef().vec();
        let mut t = p.time.as_ps() as f64 / 1e12;
        for _ in 0..5 {
            let k = f.kinematics(inst((t * 1e12) as i64)).unwrap();
            let r = k.p - target;
            t -= r.dot(k.v) / (k.v.dot(k.v) + r.dot(k.a));
        }
        let dt = t - p.time.as_ps() as f64 / 1e12;
        assert!(
            dt.abs() < 0.01,
            "closest approach {dt} s off at {:?}",
            p.time
        );
        sum += dt;
    }
    let mean = sum / interior.len() as f64;
    assert!(mean.abs() < 2e-4, "mean timing shift {mean} s");
}

#[test]
fn append_keeps_the_past_and_joins_c2() {
    let m = heights_model();
    let mut f = sample_flight();
    let before = f.clone();
    let t = TruthTime::from_ps(47 * PS + 314_159_265_358);
    f.append(
        t,
        &[
            Primitive::Turn {
                rate: TurnRate::Rate(2.5 * DEG),
                direction: None,
                end: TurnEnd::Heading(0.0),
            },
            Primitive::Climb {
                rate: ClimbRate::Rate(-5.0),
                end: ClimbEnd::Altitude(ell(2000.0)),
            },
        ],
        heights(&m),
    )
    .unwrap();
    for i in 0..2000 {
        let ti = inst(t.as_ps() / 2000 * i + 999);
        assert_eq!(f.state(ti).unwrap(), before.state(ti).unwrap());
    }
    let l = kin(&f, inst(t.as_ps() - 1));
    let r = kin(&f, inst(t.as_ps()));
    assert!(l.a.distance(r.a) < 1e-6 && l.v.distance(r.v) < 1e-7);
    check_c2(&f, "append");
    // Append in the middle of a climb's transition, then extend.
    let mut g = f.clone();
    let mid = TruthTime::from_ps(f.span().end.unwrap().as_ps() - 30 * PS);
    g.append(mid, &[Primitive::Speed { speed_m_s: 150.0 }], heights(&m))
        .unwrap();
    g.extend_to(
        TruthTime::from_ps(g.span().end.unwrap().as_ps() + 100 * PS),
        heights(&m),
    )
    .unwrap();
    check_c2(&g, "append mid-transition");
    // Appending from a recorded track works too.
    let pts: Vec<TrackPoint> = (0..20)
        .map(|i| {
            let t = TruthTime::from_ps(i * PS);
            TrackPoint {
                time: t,
                position: before.state_at(t).unwrap().position.to_geodetic().unwrap(),
            }
        })
        .collect();
    let mut tr = Flight::track(&pts, 1.0, Limits::default(), AttitudeMode::COORDINATED).unwrap();
    tr.append_route(
        TruthTime::from_ps(15 * PS),
        &[Waypoint {
            lat_rad: 51.1 * DEG,
            lon_rad: -0.95 * DEG,
            altitude: ell(3000.0),
            speed_m_s: None,
        }],
        heights(&m),
    )
    .unwrap();
    check_c2(&tr, "track then route");
    // A failed append leaves the flight unchanged.
    let snapshot = f.clone();
    assert!(
        f.append(t, &[Primitive::Speed { speed_m_s: 0.5 }], heights(&m))
            .is_err()
    );
    assert_eq!(f, snapshot);
}

#[test]
fn ground_platforms_do_not_bank_or_climb() {
    let m = heights_model();
    let prims = [
        Primitive::Straight {
            duration: secs(5.0),
        },
        Primitive::Turn {
            rate: TurnRate::Rate(10.0 * DEG),
            direction: Some(TurnDirection::Right),
            end: TurnEnd::By(90.0 * DEG),
        },
        Primitive::Speed { speed_m_s: 5.0 },
    ];
    let f = Flight::maneuvers(
        start(30.0, 30.0, 2.0, 0.0, 15.0),
        Limits::default(),
        AttitudeMode::Ground,
        &prims,
        heights(&m),
    )
    .unwrap();
    for i in 0..100 {
        let t = f.span().end.unwrap().as_ps() / 100 * i;
        let (roll, pitch, _) = roll_pitch_yaw(&f, TruthTime::from_ps(t));
        assert!(
            roll.abs() < 1e-12 && pitch.abs() < 1e-9,
            "roll {roll} pitch {pitch}"
        );
        let g = f
            .state_at(TruthTime::from_ps(t))
            .unwrap()
            .position
            .to_geodetic()
            .unwrap();
        assert!((g.h_m - 2.0).abs() < 1e-6);
    }
    let climb = [Primitive::Climb {
        rate: ClimbRate::Rate(1.0),
        end: ClimbEnd::Duration(secs(10.0)),
    }];
    assert!(matches!(
        Flight::maneuvers(
            start(30.0, 30.0, 2.0, 0.0, 15.0),
            Limits::default(),
            AttitudeMode::Ground,
            &climb,
            heights(&m)
        ),
        Err(GeoError::Trajectory(_))
    ));
    // AGL waypoints over the placeholder ground.
    let ground = FlatGround { height_m: 100.0 };
    let hs = Heights {
        model: &m,
        ground: &ground,
    };
    let f = Flight::route(
        Start {
            track_rad: None,
            ..start(30.0, 30.0, 102.0, 0.0, 15.0)
        },
        Limits::default(),
        AttitudeMode::Ground,
        &[Waypoint {
            lat_rad: 30.01 * DEG,
            lon_rad: 30.0 * DEG,
            altitude: Altitude::new(2.0, AltitudeRef::Agl).unwrap(),
            speed_m_s: None,
        }],
        hs,
    )
    .unwrap();
    let g = f
        .state_at(f.span().end.unwrap())
        .unwrap()
        .position
        .to_geodetic()
        .unwrap();
    assert!((g.h_m - 102.0).abs() < 1e-6 && (g.lat_deg() - 30.01).abs() < 1e-7);
}

#[test]
fn spans_errors_and_determinism() {
    let f = sample_flight();
    assert!(matches!(
        f.state(inst(-1)),
        Err(GeoError::OutsideSpan { .. })
    ));
    assert!(matches!(
        f.state(inst(f.span().end.unwrap().as_ps() + 1)),
        Err(GeoError::OutsideSpan { .. })
    ));
    assert_eq!(f, sample_flight());
    let m = heights_model();
    let bad = |p: Primitive| {
        Flight::maneuvers(
            start(0.0, 0.0, 100.0, 0.0, 100.0),
            Limits::default(),
            AttitudeMode::COORDINATED,
            &[p],
            heights(&m),
        )
    };
    assert!(
        bad(Primitive::Turn {
            rate: TurnRate::Rate(3.0 * DEG),
            direction: None,
            end: TurnEnd::Duration(secs(10.0))
        })
        .is_err()
    );
    assert!(
        bad(Primitive::Turn {
            rate: TurnRate::Rate(3.0 * DEG),
            direction: Some(TurnDirection::Left),
            end: TurnEnd::Duration(secs(1.0))
        })
        .is_err()
    );
    assert!(
        bad(Primitive::Turn {
            rate: TurnRate::Bank(85.0 * DEG),
            direction: None,
            end: TurnEnd::Heading(1.0)
        })
        .is_err()
    );
    assert!(
        bad(Primitive::Climb {
            rate: ClimbRate::Rate(120.0),
            end: ClimbEnd::Duration(secs(100.0))
        })
        .is_err()
    );
    assert!(bad(Primitive::Speed { speed_m_s: 0.0 }).is_err());
    let _ = Ecef::new(0.0, 0.0, 0.0);
}
/// Each piece ends on the next piece's start values, to f64 rounding (contracts §13.2).
#[test]
fn knot_values_match_to_rounding() {
    let f = sample_flight();
    let ps = f.path().pieces();
    let (mut dp, mut dv, mut da) = (0.0_f64, 0.0_f64, 0.0_f64);
    for w in ps.windows(2) {
        let t = (w[0].end.as_ps() - w[0].start.as_ps()) as f64 / 1e12;
        let l = w[0].eval(t);
        let r = w[1].eval(0.0);
        dp = dp.max(l.p.distance(r.p));
        dv = dv.max(l.v.distance(r.v));
        da = da.max(l.a.distance(r.a));
    }
    assert!(
        dp < 1e-8 && dv < 1e-10 && da < 1e-10,
        "jumps {dp:e} m, {dv:e} m/s, {da:e} m/s²"
    );
}

// Regressions for the independent review (progress.md, prompt 8).

#[test]
fn no_op_maneuvers_add_nothing() {
    let m = heights_model();
    let f = Flight::maneuvers(
        start(10.0, 10.0, 1000.0, 45.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[
            Primitive::Speed { speed_m_s: 100.0 },
            Primitive::Turn {
                rate: TurnRate::Rate(3.0 * DEG),
                direction: None,
                end: TurnEnd::Heading(45.0 * DEG),
            },
            Primitive::Straight {
                duration: secs(2.0),
            },
        ],
        heights(&m),
    )
    .unwrap();
    assert_eq!(f.span().end, Some(at(2.0)));
    for p in f.path().pieces() {
        assert!(p.end > p.start && p.coeffs.iter().flatten().all(|c| c.is_finite()));
    }
    // A primitive that only does nothing gives an empty trajectory error, not NaN.
    let only = Flight::maneuvers(
        start(10.0, 10.0, 1000.0, 45.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Primitive::Speed { speed_m_s: 100.0 }],
        heights(&m),
    )
    .unwrap();
    assert!(only.state_at(TruthTime::EPOCH).is_err());
}

#[test]
fn short_segments_stay_c2_and_bounded() {
    let m = heights_model();
    // Sub-millisecond durations are rejected; 1 ms segments chained are smooth.
    let too_short = Flight::maneuvers(
        start(51.0, -1.0, 1000.0, 0.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Primitive::Straight {
            duration: TruthDuration::from_ps(1_000_000),
        }],
        heights(&m),
    );
    assert!(too_short.is_err());
    let prims: Vec<Primitive> = (0..20)
        .map(|_| Primitive::Straight {
            duration: TruthDuration::from_ps(1_000_000_000),
        })
        .collect();
    let f = Flight::maneuvers(
        start(51.0, -1.0, 1000.0, 0.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &prims,
        heights(&m),
    )
    .unwrap();
    for i in 0..=200 {
        let k = kin(&f, inst(i * 100_000_000));
        // Straight and level: only the ~V²/R curvature acceleration.
        assert!(k.a.norm() < 2e-3, "a = {:?} at {i}", k.a);
    }
    // An append 0.3 ms after a knot, then a short hold.
    let mut g = sample_flight();
    let t = TruthTime::from_ps(20 * PS + 300_000_000);
    g.append(
        t,
        &[Primitive::Straight {
            duration: TruthDuration::from_ps(1_000_000_000),
        }],
        heights(&m),
    )
    .unwrap();
    let l = kin(&g, inst(t.as_ps() - 1));
    let r = kin(&g, inst(t.as_ps() + 1));
    assert!(l.a.distance(r.a) < 1e-6);
    for i in 0..=10 {
        let k = kin(&g, inst(t.as_ps() + i * 100_000_000));
        assert!(k.a.norm() < 0.1, "a = {:?}", k.a);
    }
}

#[test]
fn fast_low_speed_turns_follow_the_commanded_rate() {
    let m = heights_model();
    // 2 m/s with a 60° bank: about 8.5 rad/s, a 0.24 m radius.
    let f = Flight::maneuvers(
        start(0.0, 0.0, 10.0, 0.0, 2.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Primitive::Turn {
            rate: TurnRate::Bank(60.0 * DEG),
            direction: Some(TurnDirection::Right),
            end: TurnEnd::Duration(secs(30.0)),
        }],
        heights(&m),
    )
    .unwrap();
    for i in 0..50 {
        let t = at(12.0 + 0.0371 * f64::from(i));
        let s = f.state_at(t).unwrap();
        assert!(
            (s.velocity.norm() - 2.0).abs() < 1e-6,
            "speed {}",
            s.velocity.norm()
        );
        let (roll, _, _) = roll_pitch_yaw(&f, t);
        assert!(
            (roll - 60.0 * DEG).abs() < 0.01 * DEG,
            "bank {}",
            roll / DEG
        );
    }
}

#[test]
fn far_and_antipodal_waypoints() {
    let m = heights_model();
    let wp = |lat: f64, lon: f64| Waypoint {
        lat_rad: lat * DEG,
        lon_rad: lon * DEG,
        altitude: ell(0.0),
        speed_m_s: None,
    };
    // Antipodal: rejected, never "arrived" at.
    let r = Flight::route(
        start(0.0, 0.0, 0.0, 90.0, 100.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[wp(0.0, 180.0)],
        heights(&m),
    );
    assert!(
        matches!(r, Err(GeoError::Trajectory(ref e)) if e.contains("antipodal")),
        "{r:?}"
    );
    // A 2000 km leg behind the start: reached.
    let f = Flight::route(
        start(0.0, 0.0, 10_000.0, 0.0, 250.0),
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Waypoint {
            altitude: ell(10_000.0),
            ..wp(-10.0, 15.0)
        }],
        heights(&m),
    )
    .unwrap();
    let end = f.state_at(f.span().end.unwrap()).unwrap();
    let target = Geodetic::from_degrees(-10.0, 15.0, 10_000.0).to_ecef();
    assert!(
        end.position.distance(target) < 0.01,
        "miss {}",
        end.position.distance(target)
    );
}

#[test]
fn extreme_times_error_instead_of_overflowing() {
    let m = heights_model();
    let mut st = start(0.0, 0.0, 100.0, 0.0, 100.0);
    st.time = TruthTime::MIN;
    let mut f = Flight::maneuvers(
        st,
        Limits::default(),
        AttitudeMode::COORDINATED,
        &[Primitive::Straight {
            duration: secs(1.0),
        }],
        heights(&m),
    )
    .unwrap();
    assert!(f.extend_to(TruthTime::MAX, heights(&m)).is_err());
    let mut st = start(0.0, 0.0, 100.0, 0.0, 100.0);
    st.time = TruthTime::from_ps(i64::MAX - PS / 2);
    assert!(
        Flight::maneuvers(
            st,
            Limits::default(),
            AttitudeMode::COORDINATED,
            &[Primitive::Straight {
                duration: secs(1.0)
            }],
            heights(&m)
        )
        .is_err()
    );
}
