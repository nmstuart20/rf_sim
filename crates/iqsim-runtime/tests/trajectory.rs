//! Scenario trajectories (prompt 8): every kind in `examples/trajectories.toml` builds
//! into an `iqsim-geo` trajectory with the expected behavior, and `sim run` resolves
//! and fingerprints every one of them (prompt 9).

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};

use iqsim_core::TruthTime;
use iqsim_geo::{AltitudeModel, Ecef, Geodetic};
use iqsim_runtime::{RunOptions, run_offline, trajectory};

fn example() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/trajectories.toml")
}

fn scenario() -> iqsim_scenario::Scenario {
    let source = std::fs::read_to_string(example()).unwrap();
    iqsim_scenario::parse(&source).unwrap()
}

#[test]
fn every_example_trajectory_builds() {
    let s = scenario();
    let model = AltitudeModel::default();
    let build = |name: &str| {
        let p = s
            .platforms
            .iter()
            .find(|p| p.name.as_str() == name)
            .unwrap();
        trajectory::build(p, &model).unwrap()
    };
    let deg = core::f64::consts::PI / 180.0;

    // Static, heading 045.
    let site = build("site");
    let st = site.state_at(TruthTime::from_ps(i64::MAX)).unwrap();
    let (yaw, _, _) = st.q_ned_body().unwrap().to_euler_zyx();
    assert!((yaw - 45.0 * deg).abs() < 1e-12);
    assert!(site.span().end.is_none());

    // Maneuvers: the second primitive turns to 180°; the flight ends descending at
    // 280 kt after a 90° left turn.
    let ac1 = build("ac1");
    let end = ac1.span().end.unwrap();
    let st = ac1.state_at(end).unwrap();
    assert!((st.velocity.norm() - 280.0 * 1852.0 / 3600.0).abs() < 1e-6);
    let g = st.position.to_geodetic().unwrap();
    let up = Ecef::from_geodetic(g.lat_deg(), g.lon_deg(), 1.0).vec() - st.position.vec();
    assert!(st.velocity.dot(up).abs() < 1e-6, "levelled off at the end");
    assert!(ac1.knots(TruthTime::MIN, TruthTime::MAX).len() > 100);

    // Route: passes over each waypoint (MSL = ellipsoid with no geoid grid).
    let ac2 = build("ac2");
    let knots = ac2.knots(TruthTime::MIN, TruthTime::MAX);
    for (lat, lon, ft) in [(50.9, -1.3, 4000.0), (51.05, -1.25, 4000.0)] {
        let target = Geodetic::from_degrees(lat, lon, ft * 0.3048).to_ecef();
        let miss = knots
            .iter()
            .map(|&k| ac2.state_at(k).unwrap().position.distance(target))
            .fold(f64::INFINITY, f64::min);
        assert!(miss < 0.01, "route miss {miss} m");
    }

    // Track: defined exactly over its points.
    let ac3 = build("ac3");
    assert_eq!(ac3.span().start, Some(TruthTime::EPOCH));
    assert_eq!(ac3.span().end, Some(TruthTime::from_ps(20_000_000_000_000)));

    // Ground: 2 m above the (ellipsoid) ground throughout, never banked.
    let truck = build("truck");
    let e = truck.span().end.unwrap().as_ps();
    for i in 0..=50 {
        let st = truck.state_at(TruthTime::from_ps(e / 50 * i)).unwrap();
        let g = st.position.to_geodetic().unwrap();
        assert!((g.h_m - 2.0).abs() < 1e-6);
        let (_, _, roll) = st.q_ned_body().unwrap().to_euler_zyx();
        assert!(roll.abs() < 1e-12);
    }
}

#[test]
fn sim_run_resolves_every_moving_platform() {
    // Prompt 8 refused moving platforms here; since prompt 9 they are built until the
    // run ends and fingerprinted in the resolved scenario.
    let out = std::env::temp_dir().join(format!("iqsim-traj-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);
    let m = run_offline(&example(), &out, RunOptions::default()).unwrap();
    let platforms = m.scenario.resolved["platforms"].as_array().unwrap();
    let kinds: Vec<&str> = platforms
        .iter()
        .map(|p| p["kind"].as_str().unwrap())
        .collect();
    assert_eq!(kinds, ["static", "maneuvers", "route", "track", "ground"]);
    for p in platforms {
        assert_eq!(p["trajectory_sha256"].as_str().unwrap().len(), 64);
    }
    assert!(platforms[0]["span_start_ps"].is_null());
    assert_eq!(platforms[1]["span_start_ps"], 0);
    let _ = std::fs::remove_dir_all(out);
}
