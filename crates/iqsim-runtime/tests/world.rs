//! The Earth model in scenario resolution, manifests and replay (prompt 7): the
//! explicit zero-geoid fallback and its warning, geoid grids moving MSL positions and
//! registered as assets, replay rejecting a changed grid, and flight-level positions
//! resolved through the pressure-altitude chain.

#![forbid(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use iqsim_contracts::AssetKind;
use iqsim_geo::{AltitudeModel, Atmosphere, Ecef, Geoid, GeoidGrid};
use iqsim_runtime::offline::load_scenario;
use iqsim_runtime::plan::ZERO_GEOID_WARNING;
use iqsim_runtime::{RunError, RunOptions, replay, run_offline};

fn scratch(tag: &str) -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "iqsim-world-{}-{tag}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// `examples/minimal.toml` shortened to 1 ms, with `extra` appended.
fn minimal(extra: &str) -> String {
    let base = fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/minimal.toml"),
    )
    .unwrap();
    let base = base.replace("duration = \"10 ms\"", "duration = \"1 ms\"");
    format!("{base}\n{extra}")
}

/// A 2.5°-spaced GeographicLib PGM grid whose heights vary with position
/// (−40 m … +40 m).
fn pgm() -> Vec<u8> {
    let width = 144;
    let height = width / 2 + 1;
    let mut out = format!(
        "P5\n# Description synthetic test geoid\n# Offset -40\n# Scale 0.002\n{width} {height}\n65535\n"
    )
    .into_bytes();
    for i in 0..height {
        for j in 0..width {
            let v = (20_000 + 7 * i * j + 113 * j + 251 * i) % 40_000;
            out.extend_from_slice(&(v as u16).to_be_bytes());
        }
    }
    out
}

fn stream_position(file: &Path, id: &str) -> Ecef {
    let l = load_scenario(file).unwrap();
    let s = l.resolved.streams.iter().find(|s| s.id == id).unwrap();
    Ecef::new(
        s.position_ecef_m[0],
        s.position_ecef_m[1],
        s.position_ecef_m[2],
    )
}

fn expect_replay_error(manifest: &Path, what: &str, needle: &str) {
    match replay(manifest, RunOptions { threads: 1 }) {
        Err(RunError::Replay(msg)) => assert!(msg.contains(needle), "{what}: {msg}"),
        other => panic!("{what}: expected a replay rejection, got {other:?}"),
    }
}

#[test]
fn zero_geoid_is_explicit_and_warned() {
    let dir = scratch("zero");
    let file = dir.join("s.toml");
    fs::write(&file, minimal("")).unwrap();
    let out = dir.join("out");
    let m = run_offline(&file, &out, RunOptions::default()).unwrap();
    assert_eq!(m.warnings, vec![ZERO_GEOID_WARNING.to_owned()]);
    assert_eq!(m.scenario.resolved["world"]["geoid"]["model"], "zero");
    assert_eq!(m.scenario.resolved["world"]["pressure_offset_pa"], 0.0);
    assert!(m.assets.iter().all(|a| a.kind != AssetKind::GeoidGrid));
    let manifest = out.join("run/manifest.json");
    replay(&manifest, RunOptions::default()).unwrap();

    // A manifest whose warning was edited away is rejected.
    let text = fs::read_to_string(&manifest).unwrap();
    let edited = text.replace(ZERO_GEOID_WARNING, "nothing to see");
    assert_ne!(text, edited);
    fs::write(&manifest, edited).unwrap();
    expect_replay_error(&manifest, "edited warning", "warnings differ");
}

#[test]
fn geoid_grid_moves_msl_positions_and_replay_checks_it() {
    let dir = scratch("grid");
    let bytes = pgm();
    fs::write(dir.join("test-geoid.pgm"), &bytes).unwrap();
    let zero = dir.join("zero.toml");
    fs::write(&zero, minimal("")).unwrap();
    let file = dir.join("s.toml");
    fs::write(
        &file,
        minimal("[world]\ngeoid = { path = \"test-geoid.pgm\" }\n"),
    )
    .unwrap();

    // The site is "25 m MSL": with the grid it sits N higher along the normal.
    let p0 = stream_position(&zero, "rx.ch0");
    let p1 = stream_position(&file, "rx.ch0");
    let g = p0.to_geodetic().unwrap();
    assert!((g.h_m - 25.0).abs() < 1e-6);
    let grid = GeoidGrid::from_pgm_bytes(&bytes).unwrap();
    let n = grid.height_m(g.lat_rad, g.lon_rad).unwrap();
    assert!(n.abs() > 1.0, "the test grid should move the site: N = {n}");
    let want = Ecef::from_geodetic(g.lat_deg(), g.lon_deg(), 25.0 + n);
    assert!(p1.distance(want) < 1e-6, "{}", p1.distance(want));
    assert!((p1.distance(p0) - n.abs()).abs() < 1e-6);
    // The beacon is "10 m AGL", still taken as ellipsoidal (no terrain yet), so it
    // does not move.
    let l0 = load_scenario(&zero).unwrap();
    let l1 = load_scenario(&file).unwrap();
    assert_eq!(
        l0.resolved.emitters[0].position_ecef_m,
        l1.resolved.emitters[0].position_ecef_m
    );

    let out = dir.join("out");
    let m = run_offline(&file, &out, RunOptions::default()).unwrap();
    assert!(m.warnings.is_empty(), "{:?}", m.warnings);
    let a = m
        .assets
        .iter()
        .find(|a| a.kind == AssetKind::GeoidGrid)
        .unwrap();
    assert_eq!(a.logical_name, "geoid");
    assert_eq!(a.format, "geographiclib-pgm");
    assert_eq!(a.byte_len, bytes.len() as u64);
    assert!(a.metadata_sha256.is_some());
    assert_eq!(m.scenario.resolved["world"]["geoid"]["model"], "grid");
    assert_eq!(m.scenario.resolved["world"]["geoid"]["width"], 144);
    let manifest = out.join("run/manifest.json");
    replay(&manifest, RunOptions::default()).unwrap();

    // Review regression: a manifest with the geoid asset removed, or with its
    // metadata hash edited, is rejected even though every listed file verifies.
    let text = fs::read_to_string(&manifest).unwrap();
    let mut v: serde_json::Value = serde_json::from_str(&text).unwrap();
    let assets = v["assets"].as_array().unwrap().clone();
    v["assets"] = serde_json::Value::Array(Vec::new());
    fs::write(&manifest, serde_json::to_string_pretty(&v).unwrap()).unwrap();
    expect_replay_error(&manifest, "omitted asset", "asset list differs");
    let mut edited = assets.clone();
    edited[0]["metadata_sha256"] = serde_json::Value::String("0".repeat(64));
    v["assets"] = serde_json::Value::Array(edited);
    fs::write(&manifest, serde_json::to_string_pretty(&v).unwrap()).unwrap();
    expect_replay_error(&manifest, "edited metadata hash", "asset list differs");
    fs::write(&manifest, &text).unwrap();
    replay(&manifest, RunOptions::default()).unwrap();

    // One changed sample byte, then a missing file: rejected before rendering.
    let path = dir.join("test-geoid.pgm");
    let mut changed = bytes.clone();
    let last = changed.len() - 1;
    changed[last] ^= 1;
    fs::write(&path, &changed).unwrap();
    expect_replay_error(
        &manifest,
        "changed byte",
        "changed since the run was recorded",
    );
    fs::remove_file(&path).unwrap();
    expect_replay_error(&manifest, "missing grid", "geoid");
    fs::write(&path, &bytes).unwrap();
    replay(&manifest, RunOptions::default()).unwrap();

    // A malformed grid is a resolution error.
    fs::write(&path, &bytes[..bytes.len() - 2]).unwrap();
    match run_offline(&file, &dir.join("out2"), RunOptions::default()) {
        Err(e) => assert!(e.to_string().contains("world.geoid"), "{e}"),
        Ok(_) => panic!("a truncated geoid grid was accepted"),
    }
}

#[test]
fn flight_level_positions_resolve_through_the_pressure_chain() {
    let dir = scratch("fl");
    let bytes = pgm();
    fs::write(dir.join("g.pgm"), &bytes).unwrap();
    let text = minimal("[world]\ngeoid = { path = \"g.pgm\" }\npressure_offset = \"-20 hPa\"\n")
        .replace("alt = \"10 m AGL\"", "alt = \"FL350\"");
    let file = dir.join("s.toml");
    fs::write(&file, text).unwrap();
    let l = load_scenario(&file).unwrap();
    let e = l.resolved.emitters[0].position_ecef_m;
    let g = Ecef::new(e[0], e[1], e[2]).to_geodetic().unwrap();
    // The forward Mode C chain on the resolved position gives back FL350.
    let model = AltitudeModel::new(
        Geoid::Grid(GeoidGrid::from_pgm_bytes(&bytes).unwrap()),
        Atmosphere::new(-2_000.0).unwrap(),
    );
    let ft = model.pressure_altitude_ft(g).unwrap();
    assert!((ft - 35_000.0).abs() < 1e-6, "{ft}");
    assert_eq!(model.mode_c_code(g).unwrap().altitude_ft(), 35_000);
    // Low pressure puts FL350 lower than on a standard day.
    let msl = model.msl_height_m(g).unwrap();
    let std = Atmosphere::STANDARD
        .msl_height_for_pressure_altitude_m(10_668.0)
        .unwrap();
    assert!(msl < std - 100.0, "{msl} vs {std}");
    assert_eq!(l.resolved.world.pressure_offset_pa, -2_000.0);
    // And the run renders.
    run_offline(&file, &dir.join("out"), RunOptions::default()).unwrap();
}
