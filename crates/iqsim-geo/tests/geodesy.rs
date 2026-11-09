//! Geodetic ↔ ECEF and Earth rotation against independent references
//! (`tests/data/gen_oracle.py`: PROJ 9.8.1 and a 60-digit mpmath reference).
//!
//! Tolerances (contracts §12.2):
//! - ECEF position and height: `max(10⁻⁶ m, 4 ulp(|r|))`. 4 ulp(|r|) exceeds 1 µm
//!   only beyond 2³⁰ m ≈ 1.07·10⁹ m, so the whole tested range (to 4·10⁸ m) is held
//!   to 1 µm.
//! - Latitude and longitude as arc lengths: `|Δφ|·(N + h)` and `|Δλ|·p`, each within
//!   the same metre tolerance, so the angular bound scales with distance and the
//!   longitude is not checked on the polar axis.

#![forbid(unsafe_code)]

use core::f64::consts::{PI, TAU};

use iqsim_core::{AcqInstant, TruthTime};
use iqsim_geo::geodetic::{MAX_INVERSE_ITERATIONS, wgs84_b_m};
use iqsim_geo::{Ecef, Geodetic, WGS84_A_M, WGS84_F, earth_rotation_angle};
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;

const POS_TOL_M: f64 = 1e-6;

fn ulp(x: f64) -> f64 {
    let x = x.abs();
    f64::from_bits(x.to_bits() + 1) - x
}

fn tol(r: f64) -> f64 {
    POS_TOL_M.max(4.0 * ulp(r))
}

fn rows(text: &str) -> Vec<Vec<f64>> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect()
}

/// Arc-length differences (lat, lon, h) between two geodetic points, metres.
fn geodetic_diff(a: Geodetic, lat_deg: f64, lon_deg: f64, h: f64) -> (f64, f64, f64) {
    let e2 = WGS84_F * (2.0 - WGS84_F);
    let (sp, cp) = libm::sincos(a.lat_rad);
    let n = WGS84_A_M / (1.0 - e2 * sp * sp).sqrt();
    let dlat = (a.lat_rad - lat_deg.to_radians()).abs() * (n + a.h_m).abs();
    let mut dl = (a.lon_rad - lon_deg.to_radians()).abs();
    if dl > PI {
        dl = TAU - dl;
    }
    let p = ((n + a.h_m) * cp).abs();
    (dlat, dl * p, (a.h_m - h).abs())
}

#[test]
fn forward_matches_proj_and_mpmath() {
    let data = rows(include_str!("data/geodetic_forward.txt"));
    assert!(data.len() > 500);
    let mut worst = 0.0_f64;
    for v in &data {
        let e = Ecef::from_geodetic(v[0], v[1], v[2]);
        let r = e.vec().norm();
        for (k, name) in [(3, "PROJ"), (6, "mpmath")] {
            let d = e.distance(Ecef::new(v[k], v[k + 1], v[k + 2]));
            worst = worst.max(d);
            assert!(d <= tol(r), "{name}: {v:?} -> {e:?}: {d} m");
        }
    }
    eprintln!("forward: worst {worst:e} m over {} points", data.len());
}

#[test]
fn inverse_matches_mpmath() {
    let data = rows(include_str!("data/geodetic_inverse.txt"));
    assert!(data.len() > 700);
    let (mut wl, mut wn, mut wh, mut steps) = (0.0_f64, 0.0_f64, 0.0_f64, 0);
    for v in &data {
        let e = Ecef::new(v[0], v[1], v[2]);
        let r = e.vec().norm();
        let (g, k) = e.to_geodetic_with_steps().unwrap();
        steps = steps.max(k);
        let (dl, dn, dh) = geodetic_diff(g, v[6], v[7], v[8]);
        wl = wl.max(dl);
        wn = wn.max(dn);
        wh = wh.max(dh);
        assert!(
            dl <= tol(r) && dn <= tol(r) && dh <= tol(r),
            "{v:?} -> {g:?}: {dl} {dn} {dh}"
        );
    }
    assert!(steps <= 4, "{steps} Bowring steps");
    eprintln!("inverse vs mpmath: lat {wl:e} m, lon {wn:e} m, h {wh:e} m, max {steps} steps");
}

/// PROJ's geocentric inverse is a fixed-step approximation: about 10⁻⁵ m at the
/// surface, degrading to centimetres at 2000 km (measured when the data was made).
/// Near the surface it still confirms the conventions (signs, longitude range).
#[test]
fn inverse_agrees_with_proj_near_the_surface() {
    let data = rows(include_str!("data/geodetic_inverse.txt"));
    let mut n = 0;
    for v in &data {
        let e = Ecef::new(v[0], v[1], v[2]);
        if e.vec().norm() > 6.4e6 {
            continue;
        }
        n += 1;
        let g = e.to_geodetic().unwrap();
        let (dl, dn, dh) = geodetic_diff(g, v[3], v[4], v[5]);
        assert!(
            dl <= 2e-5 && dn <= 2e-5 && dh <= 2e-5,
            "{v:?}: {dl} {dn} {dh}"
        );
    }
    assert!(n > 100);
}

#[test]
fn poles_dateline_and_extremes() {
    let b = wgs84_b_m();
    for h in [-11_000.0, 0.0, 400_000.0, 35_786_000.0, 3.844e8] {
        for lat in [90.0, -90.0] {
            let g = Ecef::from_geodetic(lat, 0.0, h).to_geodetic().unwrap();
            assert!((g.lat_deg() - lat).abs() * (b + h) * PI / 180.0 <= tol(b + h));
            assert!((g.h_m - h).abs() <= tol(b + h), "{lat} {h}: {}", g.h_m);
        }
        // Both sides of the antimeridian come back on their own side; ±180° itself
        // comes back as +180°.
        for (lon, want) in [
            (179.999_999_9, 179.999_999_9),
            (-179.999_999_9, -179.999_999_9),
            (180.0, 180.0),
            (-180.0, 180.0),
        ] {
            let g = Ecef::from_geodetic(12.0, lon, h).to_geodetic().unwrap();
            let r = WGS84_A_M + h;
            let mut d = (g.lon_deg() - want).abs();
            if d > 180.0 {
                d = 360.0 - d;
            }
            assert!(d * PI / 180.0 * r <= tol(r), "{lon} {h}: {}", g.lon_deg());
        }
    }
}

/// 2048 cases unless `PROPTEST_CASES` is set (AGENTS.md); failures persist to an
/// explicit path because integration tests have no lib.rs beside them.
fn config() -> ProptestConfig {
    let mut config = ProptestConfig {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "proptest-regressions/geodesy.txt",
        ))),
        ..ProptestConfig::default()
    };
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = 2048;
    }
    config
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn geodetic_round_trip(
        lat in -90.0_f64..=90.0,
        lon in -180.0_f64..180.0,
        h in prop_oneof![-11_000.0_f64..20_000.0, 20_000.0_f64..2.0e6, 2.0e6_f64..4.0e8],
    ) {
        let e = Ecef::from_geodetic(lat, lon, h);
        let r = e.vec().norm();
        let (g, steps) = e.to_geodetic_with_steps().unwrap();
        prop_assert!(steps <= 4, "{} steps", steps);
        let (dl, dn, dh) = geodetic_diff(g, lat, lon, h);
        // At the pole, longitude is arbitrary and |Δλ|·p is ~0 anyway.
        prop_assert!(dl <= tol(r) && dn <= tol(r) && dh <= tol(r), "{:?}: {} {} {}", g, dl, dn, dh);
        prop_assert!(g.to_ecef().distance(e) <= tol(r));
    }

    #[test]
    fn ecef_round_trip(
        dir in (-1.0_f64..1.0, -1.0_f64..1.0, -1.0_f64..1.0),
        r in prop_oneof![1.0e6_f64..6.3e6, 6.3e6_f64..6.4e6, 6.4e6_f64..4.0e8],
    ) {
        let (x, y, z) = dir;
        let n = (x * x + y * y + z * z).sqrt();
        prop_assume!(n > 1e-3);
        let e = Ecef::new(x / n * r, y / n * r, z / n * r);
        let (g, steps) = e.to_geodetic_with_steps().unwrap();
        prop_assert!(steps < MAX_INVERSE_ITERATIONS, "{} steps", steps);
        let back = g.to_ecef();
        prop_assert!(back.distance(e) <= tol(r), "{:?} -> {:?}: {}", e, back, back.distance(e));
    }
}

#[test]
fn earth_rotation_angle_matches_mpmath() {
    let text = include_str!("data/eci_angles.txt");
    let mut n = 0;
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split_whitespace().collect();
        let t = AcqInstant::new(
            TruthTime::from_ps(f[0].parse().unwrap()),
            f[1].parse().unwrap(),
        );
        let want: f64 = f[2].parse().unwrap();
        let got = earth_rotation_angle(t);
        assert!((0.0..=TAU).contains(&got), "{got}");
        let mut d = (got - want).abs();
        if d > PI {
            d = TAU - d;
        }
        // ≤ 0.5 ulp + 3.5·10⁻¹⁹ rad from the exact value (contracts §12.1), and the
        // reference is itself rounded to 0.5 ulp. (Review regression: the earlier
        // f64(r)·unit conversion reached 1.52 ulp, e.g. at ps = 5858921581365495166.)
        assert!(d <= ulp(want).max(3.5e-19), "{l}: {got} vs {want} ({d:e})");
        n += 1;
    }
    assert!(n > 200);
}
