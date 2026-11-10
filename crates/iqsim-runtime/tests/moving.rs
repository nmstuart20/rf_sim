//! `examples/moving.toml` end to end (prompt 9): an aircraft CW emitter rendered by the
//! light-time geometry engine through `sim run`, with the Doppler shift in the
//! antenna-port samples matching truth, the transmitter knot split logged in truth,
//! look-back before the flight starts, golden hashes, and replay.

#![forbid(unsafe_code)]

use std::f64::consts::PI;
use std::fs;
use std::path::{Path, PathBuf};

use iqsim_runtime::{RunOptions, replay, run_offline};

fn workspace_file(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(relative)
}

#[test]
fn moving_example_renders_doppler_logs_its_splits_and_replays() {
    let out = std::env::temp_dir().join(format!("iqsim-moving-{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    let m = run_offline(
        &workspace_file("examples/moving.toml"),
        &out,
        RunOptions { threads: 3 },
    )
    .unwrap();
    assert_eq!(m.outputs.len(), 2);
    let truth: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("run/truth.json")).unwrap()).unwrap();
    assert_eq!(truth["truth_schema_version"], 2);
    let link = truth["links"]
        .as_array()
        .unwrap()
        .iter()
        .find(|l| l["stream"] == "port")
        .unwrap()
        .clone();
    // Closing at about 250 kt·cos(10°): Doppler ≈ +460 Hz at 1090.2 MHz.
    let doppler = link["doppler_hz"].as_f64().unwrap();
    let range_rate = link["range_rate_m_s"].as_f64().unwrap();
    assert!((-130.0..-120.0).contains(&range_rate), "{range_rate}");
    assert!((doppler + 1090.2e6 * range_rate / 299_792_458.0).abs() < 0.01);
    // The first sample hears an emission from before the flight's start (t = 0).
    assert!(link["emission_time"][0].as_i64().unwrap() < 0);
    // The turn starts 5 ms in: one transmitter-knot split, in block 5.
    let fit = &link["geometry_fit"];
    assert_eq!(fit["blocks"], 20);
    let splits = fit["splits"].as_array().unwrap();
    assert!(
        splits
            .iter()
            .any(|s| s["reason"] == "tx_knot" && s["block"] == 5),
        "{splits:?}"
    );
    for (k, bound) in [("delay_samples", 2.5e-4), ("phase_rad", 2.5e-3)] {
        assert!(fit["max_residuals"][k].as_f64().unwrap() <= bound);
    }

    // The antenna-port tone sits at +200 kHz + Doppler: estimate it from the sc16
    // samples over the first 4 ms.
    let port = m.outputs.iter().find(|o| o.stream_id == "port").unwrap();
    let bytes = fs::read(out.join(&port.path)).unwrap();
    let iq: Vec<(f64, f64)> = bytes
        .chunks(4)
        .map(|c| {
            (
                f64::from(i16::from_le_bytes([c[0], c[1]])),
                f64::from(i16::from_le_bytes([c[2], c[3]])),
            )
        })
        .collect();
    let (mut re, mut im) = (0.0, 0.0);
    for k in 1..16_000 {
        let (a, b) = (iq[k], iq[k - 1]);
        re += a.0 * b.0 + a.1 * b.1;
        im += a.1 * b.0 - a.0 * b.1;
    }
    let f = libm::atan2(im, re) * 4e6 / (2.0 * PI);
    eprintln!("moving example: tone at {f:.3} Hz, Doppler {doppler:.3} Hz");
    assert!((f - (200e3 + doppler)).abs() < 1.0, "tone {f} Hz");

    for o in &m.outputs {
        let bytes = fs::read(out.join(&o.path)).unwrap();
        iqsim_testkit::assert_golden(
            &format!("runtime/moving_{}", o.stream_id.replace('.', "_")),
            &bytes,
        );
    }
    iqsim_testkit::assert_golden(
        "runtime/moving_truth",
        &fs::read(out.join("run/truth.json")).unwrap(),
    );
    let report = replay(&out.join("run/manifest.json"), RunOptions { threads: 1 }).unwrap();
    assert_eq!(report.outputs_matched, 2);
    let _ = fs::remove_dir_all(out);
}
