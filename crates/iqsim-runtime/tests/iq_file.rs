//! Prompt 6a acceptance tests: IQ-file emitters through the two-stage resampler.
//!
//! Fixtures are generated here with the test's own f64/libm code. References are
//! independent of the stage-1 kernel and the Farrow kernel: tone levels use the
//! responses' defining |H(f)|, and the chirp reference is a direct sum over native
//! samples of a kernel computed by time-domain Simpson quadrature of the test's own
//! sinc·Gaussian impulse responses.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use iqsim_contracts::canonical_json;
use iqsim_core::amplitude::LOSSLESS_CULL_LSB;
use iqsim_core::{
    BlockSpec, ChannelTimebase, IdealClockMap, ResponseSpec, SampleRate, TapKind, TapResponse,
    TruthTime,
};
use iqsim_geo::Ecef;
use iqsim_render::{Backend, CpuBackend};
use iqsim_runtime::static_scene::{
    ContinuousSpec, SceneError, SharedSource, StaticEmitter, StaticScene, StaticSceneBuilder,
    StaticTarget, StaticWaveform,
};
use iqsim_runtime::{RunError, RunOptions, replay, run_offline};
use iqsim_waveform::{EndPolicy, MemorySource, SourceTimeline};
use serde_json::json;
use sha2::{Digest, Sha256};

const FS: f64 = 16e6;
/// One block: 1 ms, so every 1 kHz-multiple tone has whole cycles in it.
const LEN: u64 = 16_000;
const FN_HZ: u64 = 2_000_000;
const FN: f64 = 2e6;
const C: f64 = 299_792_458.0;
const CENTER_UHZ: i64 = 1_090_000_000_000_000;
/// Not a whole number of native or target samples of delay.
const RANGE_M: f64 = 30_000.123;
const POWER_W: f64 = 1000.0;
const PI: f64 = core::f64::consts::PI;

fn mhz(x: f64) -> i64 {
    libm::round(x * 1e12) as i64
}

fn rate() -> SampleRate {
    SampleRate::from_hz(16_000_000).unwrap()
}

fn native() -> SampleRate {
    SampleRate::from_hz(FN_HZ).unwrap()
}

fn clock() -> IdealClockMap {
    IdealClockMap::new(ChannelTimebase::new(TruthTime::from_ps(0), rate()))
}

/// φ: 0.8 / 1.0 MHz full widths at 2 MS/s (the scenario defaults), 80 dB.
fn phi_spec() -> ResponseSpec {
    ResponseSpec {
        passband_edge_uhz: mhz(0.8),
        stopband_edge_uhz: mhz(1.0),
        stopband_atten_mdb: 80_000,
        group_delay_ps: 0,
    }
}

fn phi() -> TapResponse {
    TapResponse::new(TapKind::SourceReconstruction, phi_spec(), native(), 0).unwrap()
}

/// H_rx: 6 / 7.5 MHz, 80 dB, 1 µs group delay.
fn virtual_rx(cull: f64) -> StaticTarget {
    StaticTarget {
        kind: TapKind::VirtualReceiver,
        center_uhz: CENTER_UHZ,
        rate: rate(),
        response: ResponseSpec {
            passband_edge_uhz: mhz(6.0),
            stopband_edge_uhz: mhz(7.5),
            stopband_atten_mdb: 80_000,
            group_delay_ps: 1_000_000,
        },
        full_scale_volts: 1.0,
        cull_threshold_lsb: cull,
    }
}

/// B_out: 7 / 8 MHz, 60 dB, zero group delay, its own full scale.
fn antenna_port(cull: f64) -> StaticTarget {
    StaticTarget {
        kind: TapKind::AntennaPort,
        center_uhz: CENTER_UHZ,
        rate: rate(),
        response: ResponseSpec {
            passband_edge_uhz: mhz(7.0),
            stopband_edge_uhz: mhz(8.0),
            stopband_atten_mdb: 60_000,
            group_delay_ps: 0,
        },
        full_scale_volts: 0.5,
        cull_threshold_lsb: cull,
    }
}

fn targets() -> [StaticTarget; 2] {
    [
        virtual_rx(LOSSLESS_CULL_LSB),
        antenna_port(LOSSLESS_CULL_LSB),
    ]
}

fn response_of(t: &StaticTarget) -> TapResponse {
    TapResponse::new(t.kind, t.response, t.rate, t.center_uhz).unwrap()
}

fn tone(n: usize, parts: &[(f64, f64)]) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let mut v = [0.0, 0.0];
            for &(a, f) in parts {
                let (s, c) = libm::sincos(2.0 * PI * f * k as f64 / FN);
                v[0] += a * c;
                v[1] += a * s;
            }
            v
        })
        .collect()
}

fn continuous(x: Vec<[f64; 2]>, policy: EndPolicy, start_ps: i64) -> StaticWaveform {
    let peak = x
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
        .fold(0.0, f64::max);
    StaticWaveform::Continuous(ContinuousSpec {
        timeline: SourceTimeline {
            len: x.len() as u64,
            policy,
        },
        source: SharedSource::new(MemorySource(x)),
        rate: native(),
        start: TruthTime::from_ps(start_ps),
        response: phi_spec(),
        scale: 1.0,
        peak,
    })
}

fn scene(carrier_uhz: i64, waveform: StaticWaveform) -> StaticScene {
    StaticScene {
        receiver: Ecef::new(6_378_137.0, 0.0, 0.0),
        emitters: vec![StaticEmitter {
            position: Ecef::new(6_378_137.0, RANGE_M, 0.0),
            carrier_uhz,
            power_w: POWER_W,
            waveform,
        }],
    }
}

fn friis_volts(carrier_hz: f64) -> f64 {
    (POWER_W * 50.0).sqrt() * (C / carrier_hz) / (4.0 * PI * RANGE_M)
}

fn builder<'c>(
    scene: &StaticScene,
    target: &StaticTarget,
    clock: &'c IdealClockMap,
    len: u64,
) -> StaticSceneBuilder<'c> {
    StaticSceneBuilder::new(scene, target, clock, BlockSpec::new(len).unwrap()).unwrap()
}

fn render_acc(b: &StaticSceneBuilder<'_>, block: u64, backend: CpuBackend) -> Vec<i64> {
    let work = b.build_block(block).unwrap();
    let mut acc = vec![0_i64; 2 * work.len as usize];
    backend.render(&work, &mut acc).unwrap();
    acc
}

fn volts(b: &StaticSceneBuilder<'_>, block: u64) -> Vec<(f64, f64)> {
    let q = b.scale().lsb_volts();
    render_acc(b, block, CpuBackend::new(4))
        .chunks(2)
        .map(|c| (c[0] as f64 * q, c[1] as f64 * q))
        .collect()
}

/// Mean of y·e^{−j2πf(k0 + k)/F_s}: the complex amplitude of a tone at f.
fn project(y: &[(f64, f64)], k0: u64, f_hz: f64) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    for (k, &(a, b)) in y.iter().enumerate() {
        let (s, c) = libm::sincos(-2.0 * PI * f_hz * (k0 + k as u64) as f64 / FS);
        re += a * c - b * s;
        im += a * s + b * c;
    }
    (re / y.len() as f64, im / y.len() as f64)
}

fn mag(v: (f64, f64)) -> f64 {
    (v.0 * v.0 + v.1 * v.1).sqrt()
}

#[test]
fn tones_land_at_the_right_frequency_and_level_without_native_rate_images() {
    // A 2 MS/s source translated +5 MHz into a 16 MS/s channel: the translation is
    // the carrier phase at 16 MS/s, so no 2 MHz images may appear.
    let clock = clock();
    let offset = 5e6;
    for f_t in [-650_000.0, 0.0, 300_000.0, 700_000.0] {
        let x = tone(2000, &[(0.9, f_t)]);
        let sc = scene(CENTER_UHZ + mhz(5.0), continuous(x, EndPolicy::Loop, 0));
        for target in targets() {
            let resp = response_of(&target);
            let b = builder(&sc, &target, &clock, LEN);
            let y = volts(&b, 1);
            let f = offset + f_t;
            let a = project(&y, LEN, f);
            let expect = friis_volts(1095e6) * 0.9 * phi().magnitude(f_t) * resp.magnitude(f);
            let db = 20.0 * libm::log10(mag(a) / expect);
            assert!(db.abs() < 0.1, "{:?} f_t {f_t}: {db} dB", target.kind);
            assert!(
                db.abs() < 1e-3,
                "{:?} f_t {f_t}: {db} dB (precision)",
                target.kind
            );
            let mut worst_img = 0.0_f64;
            // Frequency from the phase advance between the two half blocks.
            let h = y.len() / 2;
            let p0 = project(&y[..h], LEN, f);
            let p1 = project(&y[h..], LEN + h as u64, f);
            let dphi = libm::atan2(p0.0 * p1.1 - p0.1 * p1.0, p0.0 * p1.0 + p0.1 * p1.1);
            let f_err = dphi / (2.0 * PI * (h as f64 / FS));
            assert!(f_err.abs() < 1.0, "{:?} f_t {f_t}: {f_err} Hz", target.kind);
            // Native-rate images at f ± 2k MHz inside the target band: ≤ −80 dB.
            for k in [-3.0, -2.0, -1.0, 1.0, 2.0, 3.0] {
                let g = f + k * FN;
                if g.abs() < 7.9e6 {
                    let img = mag(project(&y, LEN, g));
                    worst_img = worst_img.max(img / expect);
                    assert!(
                        img <= 1e-4 * expect,
                        "{:?} f_t {f_t}: image at {g} is {}",
                        target.kind,
                        img / expect
                    );
                }
            }
            eprintln!(
                "tone {:?} f_t {f_t}: level {db:.2e} dB, freq {f_err:.2e} Hz, worst image {worst_img:.1e}",
                target.kind
            );
        }
    }
}

#[test]
fn straddling_source_matches_the_response_and_leaves_no_alias() {
    // Emitter +8 MHz: a tone at +0.7 MHz lands at 8.7 MHz, past Nyquist (its naive alias
    // is −7.3 MHz); one at −0.7 MHz lands at 7.3 MHz, inside the transition bands.
    let clock = clock();
    let x = tone(2000, &[(0.5, 700_000.0), (0.4, -700_000.0)]);
    let sc = scene(CENTER_UHZ + mhz(8.0), continuous(x, EndPolicy::Loop, 0));
    for target in [virtual_rx(0.0), antenna_port(0.0)] {
        let resp = response_of(&target);
        let b = builder(&sc, &target, &clock, LEN);
        let y = volts(&b, 1);
        let v = friis_volts(1098e6);
        let inband = mag(project(&y, LEN, 7.3e6));
        let expect = v * 0.4 * phi().magnitude(-700_000.0) * resp.magnitude(7.3e6);
        let db = 20.0 * libm::log10(inband / expect);
        assert!(db.abs() < 0.1, "{:?}: in-band {db} dB", target.kind);
        let alias = mag(project(&y, LEN, -7.3e6));
        eprintln!(
            "straddle {:?}: in-band {db:.2e} dB, alias {:.1e} of unfiltered",
            target.kind,
            alias / (v * 0.5)
        );
        let bound = libm::pow(
            10.0,
            -f64::from(target.response.stopband_atten_mdb) / 20_000.0,
        );
        assert!(
            alias <= bound * v * 0.5,
            "{:?}: alias {} of the unfiltered level",
            target.kind,
            alias / (v * 0.5)
        );
    }
}

#[test]
fn out_of_band_iq_emitter_is_culled_or_exactly_zero() {
    // 1030 MHz into the 1090 MHz channel: the bands do not overlap, so K is zero.
    let clock = clock();
    let x = tone(2000, &[(1.0, 100_000.0)]);
    let sc = scene(1_030_000_000_000_000, continuous(x, EndPolicy::Loop, 0));
    for target in targets() {
        let b = builder(&sc, &target, &clock, LEN);
        assert!(
            b.build_block(1).unwrap().paths.is_empty(),
            "{:?}",
            target.kind
        );
    }
    for target in [virtual_rx(0.0), antenna_port(0.0)] {
        let b = builder(&sc, &target, &clock, LEN);
        let acc = render_acc(&b, 1, CpuBackend::new(2));
        assert!(acc.iter().all(|&v| v == 0), "{:?}", target.kind);
    }
}

/// The test's own truncated sinc·Gaussian impulse response (1/s).
fn imp(r: &TapResponse, t: f64) -> f64 {
    let x = t - r.group_delay_s();
    let th = (13.0 / (2.0 * PI * PI * r.sigma_hz() * r.sigma_hz())).sqrt();
    if x.abs() > th {
        return 0.0;
    }
    let two_fc = 2.0 * r.cutoff_hz();
    let arg = PI * two_fc * x;
    let sinc = if arg == 0.0 {
        1.0
    } else {
        libm::sin(arg) / arg
    };
    let s = r.sigma_hz();
    two_fc * sinc * libm::exp(-2.0 * PI * PI * s * s * x * x)
}

/// K(t) in native samples: T_n·∫ φ(u)·h(tT_n − u)·e^{−j2πΔf(tT_n − u)} du, Simpson.
fn kernel_ref(phi: &TapResponse, h: &TapResponse, offset_hz: f64, t: f64) -> (f64, f64) {
    let tn = 1.0 / FN;
    let ts = t * tn;
    let tphi = (13.0 / (2.0 * PI * PI * phi.sigma_hz() * phi.sigma_hz())).sqrt();
    let th = (13.0 / (2.0 * PI * PI * h.sigma_hz() * h.sigma_hz())).sqrt();
    let lo = (-tphi).max(ts - h.group_delay_s() - th);
    let hi = tphi.min(ts - h.group_delay_s() + th);
    if hi <= lo {
        return (0.0, 0.0);
    }
    let n = 800;
    let dx = (hi - lo) / n as f64;
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..=n {
        let u = lo + i as f64 * dx;
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let v = ts - u;
        let g = imp(phi, u) * imp(h, v);
        let (s, c) = libm::sincos(-2.0 * PI * offset_hz * v);
        re += w * g * c;
        im += w * g * s;
    }
    (re * dx / 3.0 * tn, im * dx / 3.0 * tn)
}

fn chirp(n: usize) -> Vec<[f64; 2]> {
    (0..n)
        .map(|k| {
            let k = k as f64;
            let ph = 2.0 * PI * (-0.35 * k + 0.35 * k * k / n as f64);
            let (s, c) = libm::sincos(ph);
            [0.8 * c, 0.8 * s]
        })
        .collect()
}

#[test]
fn chirp_matches_an_independent_direct_sum_at_fractional_delays() {
    let clock = clock();
    let n = 4096;
    let x = chirp(n);
    let offset_hz = 5e6;
    let carrier_hz = 1095e6;
    // Two start times move the native coordinate's fraction (by 0.37 native and 0.3
    // target samples).
    for start_ps in [0_i64, 185_000] {
        let sc = scene(
            CENTER_UHZ + mhz(5.0),
            continuous(x.clone(), EndPolicy::ZeroPad, start_ps),
        );
        let tau = RANGE_M / C;
        for target in targets() {
            let h = response_of(&target);
            let b = builder(&sc, &target, &clock, LEN);
            let y = volts(&b, 1);
            let v = friis_volts(carrier_hz);
            let k0 = 3000usize;
            let m = 34usize;
            let mut env = Vec::new();
            let mut reference = Vec::new();
            for (k, &(yr, yi)) in y.iter().enumerate().skip(k0).take(m) {
                let t = (LEN as usize + k) as f64 / FS;
                let w = FN * (t - tau - start_ps as f64 * 1e-12);
                let (mut rr, mut ri) = (0.0, 0.0);
                let lo = libm::floor(w - 75.0).max(0.0) as usize;
                let hi = (libm::ceil(w + 70.0) as usize).min(n - 1);
                for (j, xj) in x.iter().enumerate().take(hi + 1).skip(lo) {
                    let kv = kernel_ref(&phi(), &h, offset_hz, w - j as f64);
                    rr += xj[0] * kv.0 - xj[1] * kv.1;
                    ri += xj[0] * kv.1 + xj[1] * kv.0;
                }
                // Remove the carrier rotation e^{j2π(Δf·t − f_e·τ)} from the output.
                let th = 2.0 * PI * (offset_hz * t - carrier_hz * tau);
                let (s, c) = libm::sincos(-th);
                env.push(((yr * c - yi * s) / v, (yr * s + yi * c) / v));
                reference.push((rr, ri));
            }
            let err = env
                .iter()
                .zip(&reference)
                .map(|(a, r)| mag((a.0 - r.0, a.1 - r.1)))
                .fold(0.0, f64::max);
            assert!(
                err < 2e-5 * 0.8,
                "{:?} start {start_ps}: {err}",
                target.kind
            );
            // Least-squares envelope delay against the reference: ≤ 0.001 samples.
            let (mut num, mut den) = (0.0, 0.0);
            for i in 1..m - 1 {
                let d = (
                    (reference[i + 1].0 - reference[i - 1].0) / 2.0,
                    (reference[i + 1].1 - reference[i - 1].1) / 2.0,
                );
                let e = (env[i].0 - reference[i].0, env[i].1 - reference[i].1);
                num += d.0 * e.0 + d.1 * e.1;
                den += d.0 * d.0 + d.1 * d.1;
            }
            let delay = num / den;
            eprintln!(
                "chirp {:?} start {start_ps}: max error {:.2e} of peak, delay {delay:.1e} samples",
                target.kind,
                err / 0.8
            );
            assert!(
                delay.abs() < 1e-3,
                "{:?}: delay error {delay} samples",
                target.kind
            );
        }
    }
}

#[test]
fn output_is_independent_of_threads_chunks_and_block_split() {
    let clock = clock();
    let x = chirp(4096);
    let sc = scene(
        CENTER_UHZ + mhz(2.5),
        continuous(x, EndPolicy::Loop, 333_000),
    );
    for target in targets() {
        let b = builder(&sc, &target, &clock, LEN);
        let reference = render_acc(&b, 2, CpuBackend::new(1));
        for backend in [
            CpuBackend::new(2),
            CpuBackend::new(8),
            CpuBackend::new(3).with_chunk_len(1),
            CpuBackend::new(4).with_chunk_len(777),
        ] {
            assert_eq!(render_acc(&b, 2, backend), reference, "{:?}", target.kind);
        }
        // The same samples from blocks of 4000 and 1777: continuous across block
        // boundaries, within f32 rounding of the coordinate.
        let full: Vec<(f64, f64)> = [1, 2].iter().flat_map(|&k| volts(&b, k)).collect();
        let peak = full.iter().map(|&v| mag(v)).fold(0.0, f64::max);
        assert!(peak > 0.0);
        for split in [4000_u64, 1777] {
            let bs = builder(&sc, &target, &clock, split);
            let first = LEN;
            let last = 3 * LEN;
            let mut got = Vec::new();
            let mut blk = first / split;
            while blk * split < last {
                let y = volts(&bs, blk);
                for (i, v) in y.into_iter().enumerate() {
                    let n = blk * split + i as u64;
                    if (first..last).contains(&n) {
                        got.push(v);
                    }
                }
                blk += 1;
            }
            let diff = got
                .iter()
                .zip(&full)
                .map(|(a, b)| mag((a.0 - b.0, a.1 - b.1)))
                .fold(0.0, f64::max);
            eprintln!(
                "split {:?} {split}: {:.1e} of peak",
                target.kind,
                diff / peak
            );
            assert!(
                diff <= 1e-6 * peak,
                "{:?} split {split}: {}",
                target.kind,
                diff / peak
            );
        }
    }
}

#[test]
fn end_policies_are_deterministic() {
    let clock = clock();
    // 2000 samples = 1 ms, a whole number of 5 MHz carrier cycles: a looped source
    // repeats bit for bit, block after block.
    let x = tone(2000, &[(0.6, 250_000.0), (0.3, -410_000.0)]);
    let looped = scene(
        CENTER_UHZ + mhz(5.0),
        continuous(x.clone(), EndPolicy::Loop, 0),
    );
    for target in targets() {
        let b = builder(&looped, &target, &clock, LEN);
        let a1 = render_acc(&b, 1, CpuBackend::new(2));
        assert!(a1.iter().any(|&v| v != 0));
        assert_eq!(
            render_acc(&b, 2, CpuBackend::new(3)),
            a1,
            "{:?}",
            target.kind
        );
        assert_eq!(
            render_acc(&b, 5, CpuBackend::new(1)),
            a1,
            "{:?}",
            target.kind
        );
    }
    // Zero pad: after the last sample and the filter tail, exactly nothing.
    let padded = scene(
        CENTER_UHZ + mhz(5.0),
        continuous(x.clone(), EndPolicy::ZeroPad, 0),
    );
    for target in targets() {
        let b = builder(&padded, &target, &clock, LEN);
        assert!(
            render_acc(&b, 0, CpuBackend::new(2))
                .iter()
                .any(|&v| v != 0)
        );
        let late = b.build_block(3).unwrap();
        assert!(late.paths.is_empty(), "{:?}", target.kind);
    }
    // Stop: the run may not outlast the file. Block 0 ends at 1 ms − 100 µs of
    // retarded time and is fine; block 1 needs samples past the end.
    let stopped = scene(CENTER_UHZ + mhz(5.0), continuous(x, EndPolicy::Stop, 0));
    for target in targets() {
        let b = builder(&stopped, &target, &clock, LEN);
        b.build_block(0).unwrap();
        assert!(matches!(
            b.build_block(1),
            Err(SceneError::SourceExhausted { index: 0, block: 1 })
        ));
    }
}

// ------------------------------------------------------------------ files and runs

fn scratch(tag: &str) -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "iqsim-iqfile-{}-{tag}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::Relaxed)
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn workspace_file(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(relative)
}

/// Writes a SigMF pair `<dir>/<stem>.sigmf-{meta,data}`.
fn write_sigmf(dir: &Path, stem: &str, cf32: bool, rate_hz: f64, x: &[[f64; 2]]) {
    let mut data = Vec::new();
    for v in x {
        for c in v {
            if cf32 {
                data.extend_from_slice(&(*c as f32).to_le_bytes());
            } else {
                let q = (c * 32768.0).round_ties_even().clamp(-32768.0, 32767.0) as i16;
                data.extend_from_slice(&q.to_le_bytes());
            }
        }
    }
    let meta = json!({
        "global": {
            "core:datatype": if cf32 { "cf32_le" } else { "ci16_le" },
            "core:sample_rate": rate_hz,
            "core:version": "1.2.0",
        },
        "captures": [{ "core:sample_start": 0 }],
        "annotations": [],
    });
    fs::write(dir.join(format!("{stem}.sigmf-data")), data).unwrap();
    fs::write(
        dir.join(format!("{stem}.sigmf-meta")),
        canonical_json(&meta),
    )
    .unwrap();
}

/// A scenario with one IQ-file emitter 3 km from a site with a channel and a cf32
/// antenna port.
fn scenario_toml(iq_file: &str, power: &str, duration: &str) -> String {
    format!(
        r#"schema_version = 1
epoch = 2026-01-01T00:00:00Z
seed = 7
run_mode = "offline"

[metadata]
name = "iq-test"

[render]
duration = "{duration}"
block = "1 ms"

[[platforms]]
name = "src"
position = {{ lat = 38.9, lon = -77.0, alt = "10 m HAE" }}

[[platforms]]
name = "site"
position = {{ lat = 38.927, lon = -77.0, alt = "10 m HAE" }}

[[emitters]]
name = "file"
platform = "src"
frequency = "1092 MHz"
power = "{power}"
behavior = "always_on"
waveform = "iq_file"
iq_file = {iq_file}

[[receivers]]
name = "rx"
platform = "site"

[[receivers.channels]]
name = "ch0"
center_frequency = "1090 MHz"
sample_rate = "10 MS/s"
bandwidth = "8 MHz"
front_end = {{ gain = "30 dB", noise_figure = "5 dB" }}
adc = {{ full_scale = "-30 dBm", bits = 16 }}
sinks = ["rec"]

[[antenna_ports]]
name = "port"
platform = "site"
center_frequency = "1090 MHz"
sample_rate = "10 MS/s"
bandwidth = "8 MHz"
full_scale = "-40 dBm"
sinks = ["port_rec"]

[[sinks]]
name = "rec"
kind = "sigmf"

[[sinks]]
name = "port_rec"
kind = "sigmf"
format = "cf32_le"
"#
    )
}

/// Range of the scenario's link, from the same geodesy the run uses.
fn link_range() -> f64 {
    let a = Ecef::from_geodetic(38.9, -77.0, 10.0);
    let b = Ecef::from_geodetic(38.927, -77.0, 10.0);
    a.distance(b)
}

#[test]
fn port_power_matches_the_link_budget_for_both_formats_and_normalizations() {
    let f_t = 400_000.0;
    let a = 0.5;
    let x = tone(20_000, &[(a, f_t)]);
    let dir = scratch("power");
    for cf32 in [false, true] {
        let stem = if cf32 { "tone_cf32" } else { "tone_ci16" };
        write_sigmf(&dir, stem, cf32, FN, &x);
        for (norm, level) in [("full_scale", a * a), ("rms", 1.0)] {
            let toml = scenario_toml(
                &format!(r#"{{ path = "{stem}", end = "loop", normalization = "{norm}" }}"#),
                "30 dBm",
                "3 ms",
            );
            let file = dir.join(format!("{stem}-{norm}.toml"));
            fs::write(&file, toml).unwrap();
            let out = dir.join(format!("out-{stem}-{norm}"));
            let m = run_offline(&file, &out, RunOptions { threads: 2 }).unwrap();
            let port = m.outputs.iter().find(|o| o.stream_id == "port").unwrap();
            let bytes = fs::read(out.join(&port.path)).unwrap();
            let fs_eff = iqsim_core::TargetScale::new(
                iqsim_core::amplitude::rms_volts(1e-7), // −40 dBm
            )
            .unwrap()
            .full_scale_effective_volts();
            let y: Vec<(f64, f64)> = bytes
                .chunks(8)
                .map(|c| {
                    let re = f32::from_le_bytes(c[0..4].try_into().unwrap());
                    let im = f32::from_le_bytes(c[4..8].try_into().unwrap());
                    (f64::from(re) * fs_eff, f64::from(im) * fs_eff)
                })
                .collect();
            // Steady state: the last 2 ms (20 000 samples at 10 MS/s, whole cycles).
            let k0 = 10_000;
            let rate = 10e6;
            let f = 2e6 + f_t;
            let (mut re, mut im) = (0.0, 0.0);
            for (k, v) in y[k0..].iter().enumerate() {
                let (s, c) = libm::sincos(-2.0 * PI * f * (k0 + k) as f64 / rate);
                re += v.0 * c - v.1 * s;
                im += v.0 * s + v.1 * c;
            }
            let nk = (y.len() - k0) as f64;
            let p = (re * re + im * im) / (nk * nk) / 50.0;
            let path = (C / 1092e6) / (4.0 * PI * link_range());
            let bout = TapResponse::new(
                TapKind::AntennaPort,
                ResponseSpec {
                    passband_edge_uhz: mhz(4.0),
                    stopband_edge_uhz: mhz(5.0),
                    stopband_atten_mdb: 80_000,
                    group_delay_ps: 0,
                },
                SampleRate::from_hz(10_000_000).unwrap(),
                CENTER_UHZ,
            )
            .unwrap();
            let g = phi().magnitude(f_t) * bout.magnitude(f);
            let expect = 1.0 * level * path * path * g * g;
            let db = 10.0 * libm::log10(p / expect);
            eprintln!("port power {stem} {norm}: {db:.4} dB");
            assert!(db.abs() < 0.1, "{stem} {norm}: {db} dB");
        }
    }
}

#[test]
fn replay_rejects_changed_iq_bytes_and_metadata() {
    let dir = scratch("replay");
    write_sigmf(&dir, "rec", false, FN, &chirp(8192));
    let file = dir.join("s.toml");
    fs::write(
        &file,
        scenario_toml(
            r#"{ path = "rec.sigmf-meta", start = "0.25 ms", end = "loop", normalization = "rms" }"#,
            "20 dBm",
            "2 ms",
        ),
    )
    .unwrap();
    let out = dir.join("out");
    let m = run_offline(&file, &out, RunOptions { threads: 3 }).unwrap();
    assert_eq!(m.assets.len(), 2);
    assert!(
        m.assets
            .iter()
            .all(|a| a.metadata_sha256.is_some() || a.format == "sigmf-meta")
    );
    let manifest = out.join("run/manifest.json");
    assert_eq!(
        replay(&manifest, RunOptions { threads: 1 })
            .unwrap()
            .outputs_matched,
        2
    );

    let expect = |what: &str, needle: &str| match replay(&manifest, RunOptions { threads: 2 }) {
        Err(RunError::Replay(msg)) => assert!(msg.contains(needle), "{what}: {msg}"),
        other => panic!("{what}: expected a replay rejection, got {other:?}"),
    };
    let data = dir.join("rec.sigmf-data");
    let original = fs::read(&data).unwrap();
    let mut changed = original.clone();
    changed[1001] ^= 1;
    fs::write(&data, &changed).unwrap();
    expect("one IQ byte", "changed since the run was recorded");
    fs::write(&data, &original).unwrap();
    let meta = dir.join("rec.sigmf-meta");
    let meta_text = fs::read(&meta).unwrap();
    let mut edited = meta_text.clone();
    edited.push(b'\n');
    fs::write(&meta, &edited).unwrap();
    expect("metadata", "changed since the run was recorded");
    fs::write(&meta, &meta_text).unwrap();
    fs::remove_file(&data).unwrap();
    expect("missing data", "missing");
    fs::write(&data, &original).unwrap();
    replay(&manifest, RunOptions { threads: 2 }).unwrap();
}

#[test]
fn stop_policy_fails_a_run_that_outlasts_the_recording() {
    let dir = scratch("stop");
    write_sigmf(&dir, "short", true, FN, &chirp(2000)); // 1 ms
    let file = dir.join("s.toml");
    fs::write(
        &file,
        scenario_toml(
            r#"{ path = "short", end = "stop", normalization = "full_scale" }"#,
            "20 dBm",
            "3 ms",
        ),
    )
    .unwrap();
    let e = run_offline(&file, &dir.join("out"), RunOptions { threads: 1 }).unwrap_err();
    assert!(e.to_string().contains("ends before block"), "{e}");
}

#[test]
fn malformed_recordings_are_rejected_at_resolution() {
    let dir = scratch("bad");
    write_sigmf(&dir, "ok", false, FN, &chirp(100));
    let run = |iq: &str, tag: &str| {
        let file = dir.join(format!("{tag}.toml"));
        fs::write(&file, scenario_toml(iq, "0 dBm", "1 ms")).unwrap();
        run_offline(
            &file,
            &dir.join(format!("out-{tag}")),
            RunOptions { threads: 1 },
        )
        .unwrap_err()
        .to_string()
    };
    // Truncated data.
    let mut d = fs::read(dir.join("ok.sigmf-data")).unwrap();
    d.pop();
    fs::write(dir.join("trunc.sigmf-data"), &d).unwrap();
    fs::copy(dir.join("ok.sigmf-meta"), dir.join("trunc.sigmf-meta")).unwrap();
    let e = run(
        r#"{ path = "trunc", end = "loop", normalization = "rms" }"#,
        "trunc",
    );
    assert!(e.contains("truncated"), "{e}");
    // Unsupported datatype.
    let meta = fs::read_to_string(dir.join("ok.sigmf-meta")).unwrap();
    fs::write(
        dir.join("be.sigmf-meta"),
        meta.replace("ci16_le", "ci16_be"),
    )
    .unwrap();
    fs::copy(dir.join("ok.sigmf-data"), dir.join("be.sigmf-data")).unwrap();
    let e = run(
        r#"{ path = "be", end = "loop", normalization = "rms" }"#,
        "be",
    );
    assert!(e.contains("unsupported core:datatype"), "{e}");
    // A stated rate that differs from the recording.
    let e = run(
        r#"{ path = "ok", end = "loop", normalization = "rms", sample_rate = "2.5 MS/s" }"#,
        "rate",
    );
    assert!(e.contains("differs from the recording"), "{e}");
    // Missing file.
    let e = run(
        r#"{ path = "nope", end = "loop", normalization = "rms" }"#,
        "missing",
    );
    assert!(e.contains("nope.sigmf-meta"), "{e}");
    // Non-finite cf32.
    let mut nan = Vec::new();
    for _ in 0..10 {
        nan.extend_from_slice(&0.5_f32.to_le_bytes());
    }
    nan.extend_from_slice(&f32::NAN.to_le_bytes());
    nan.extend_from_slice(&0.0_f32.to_le_bytes());
    fs::write(dir.join("nan.sigmf-data"), &nan).unwrap();
    fs::write(
        dir.join("nan.sigmf-meta"),
        meta.replace("ci16_le", "cf32_le"),
    )
    .unwrap();
    let e = run(
        r#"{ path = "nan", end = "loop", normalization = "rms" }"#,
        "nan",
    );
    assert!(e.contains("non-finite"), "{e}");
}

// ------------------------------------------------------------------ example

/// The example recording: 10 ms at 2 MS/s of a ±0.6 MHz chirp plus a 0.25 MHz tone,
/// ci16. Regenerated here so its provenance is checked.
fn example_fixture() -> (Vec<u8>, Vec<u8>) {
    let n = 20_000usize;
    let mut data = Vec::with_capacity(4 * n);
    for k in 0..n {
        let t = k as f64;
        let ph = 2.0 * PI * (-0.3 * t + 0.3 * t * t / n as f64);
        let (s1, c1) = libm::sincos(ph);
        let (s2, c2) = libm::sincos(2.0 * PI * 0.125 * t);
        for v in [0.5 * c1 + 0.25 * c2, 0.5 * s1 + 0.25 * s2] {
            let q = (v * 32768.0).round_ties_even() as i16;
            data.extend_from_slice(&q.to_le_bytes());
        }
    }
    let meta = json!({
        "global": {
            "core:datatype": "ci16_le",
            "core:sample_rate": 2_000_000.0,
            "core:version": "1.2.0",
            "core:description": "iqsim example: 10 ms, chirp ±0.6 MHz at 0.5 FS plus a 0.25 MHz tone at 0.25 FS",
        },
        "captures": [{ "core:sample_start": 0, "core:frequency": 1_092_000_000.0 }],
        "annotations": [],
    });
    let mut meta = canonical_json(&meta);
    meta.push(b'\n');
    (meta, data)
}

#[test]
fn example_recording_matches_its_generator() {
    let (meta, data) = example_fixture();
    let dir = workspace_file("examples/iq");
    if std::env::var_os("IQSIM_WRITE_FIXTURES").is_some() {
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("chirp_2msps.sigmf-meta"), &meta).unwrap();
        fs::write(dir.join("chirp_2msps.sigmf-data"), &data).unwrap();
    }
    assert_eq!(fs::read(dir.join("chirp_2msps.sigmf-meta")).unwrap(), meta);
    assert_eq!(fs::read(dir.join("chirp_2msps.sigmf-data")).unwrap(), data);
}

/// Golden hashes for examples/iq_file.toml, and replay of it.
#[test]
fn iq_file_example_golden_hashes() {
    let dir = scratch("example");
    let out = dir.join("out");
    let m = run_offline(
        &workspace_file("examples/iq_file.toml"),
        &out,
        RunOptions { threads: 3 },
    )
    .unwrap();
    assert_eq!(m.outputs.len(), 2);
    for o in &m.outputs {
        let bytes = fs::read(out.join(&o.path)).unwrap();
        assert_eq!(hex::encode(Sha256::digest(&bytes)), o.sha256);
        assert_eq!(o.adc_clipped_components, 0);
        iqsim_testkit::assert_golden(
            &format!("runtime/iq_file_{}", o.stream_id.replace('.', "_")),
            &bytes,
        );
    }
    let truth = fs::read(out.join("run/truth.json")).unwrap();
    iqsim_testkit::assert_golden("runtime/iq_file_truth", &truth);
    let report = replay(&out.join("run/manifest.json"), RunOptions { threads: 1 }).unwrap();
    assert_eq!(report.outputs_matched, 2);
}

// ------------------------------------------------------------------ review regressions

fn spec_of(x: Vec<[f64; 2]>, policy: EndPolicy, scale: f64) -> ContinuousSpec {
    let StaticWaveform::Continuous(mut c) = continuous(x, policy, 0) else {
        unreachable!()
    };
    c.scale = scale;
    c
}

#[test]
fn rms_scaling_of_a_tiny_recording_matches_a_unit_one() {
    // Review P1: a 1e-30 source with a 1e30 factor was culled as silent.
    let clock = clock();
    let unit = tone(2000, &[(0.9, 300_000.0)]);
    let tiny: Vec<[f64; 2]> = unit.iter().map(|v| [v[0] * 1e-30, v[1] * 1e-30]).collect();
    for target in targets() {
        let a = builder(
            &scene(
                CENTER_UHZ + mhz(5.0),
                StaticWaveform::Continuous(spec_of(unit.clone(), EndPolicy::Loop, 1.0)),
            ),
            &target,
            &clock,
            LEN,
        );
        let b = builder(
            &scene(
                CENTER_UHZ + mhz(5.0),
                StaticWaveform::Continuous(spec_of(tiny.clone(), EndPolicy::Loop, 1e30)),
            ),
            &target,
            &clock,
            LEN,
        );
        let ya = render_acc(&a, 1, CpuBackend::new(2));
        let yb = render_acc(&b, 1, CpuBackend::new(2));
        let peak = ya.iter().map(|v| v.unsigned_abs()).max().unwrap();
        assert!(peak > 1 << 20, "{:?}", target.kind);
        let diff = ya
            .iter()
            .zip(&yb)
            .map(|(p, q)| (p - q).unsigned_abs())
            .max()
            .unwrap();
        assert!(
            diff as f64 <= 1e-6 * peak as f64,
            "{:?}: {diff} of {peak}",
            target.kind
        );
    }
}

#[test]
fn one_source_with_two_end_policies_keeps_them_apart() {
    // Review P2: links sharing a source were grouped without their end policy.
    let clock = clock();
    let x = tone(2000, &[(0.5, 200_000.0)]);
    let zp = spec_of(x, EndPolicy::ZeroPad, 1.0);
    let mut lp = zp.clone();
    lp.timeline.policy = EndPolicy::Loop;
    assert_eq!(zp.source, lp.source);
    for order in [[zp.clone(), lp.clone()], [lp.clone(), zp.clone()]] {
        let mut sc = scene(CENTER_UHZ + mhz(5.0), StaticWaveform::Cw);
        let e = sc.emitters[0].clone();
        sc.emitters = order
            .iter()
            .map(|c| StaticEmitter {
                waveform: StaticWaveform::Continuous(c.clone()),
                ..e.clone()
            })
            .collect();
        for target in targets() {
            let b = builder(&sc, &target, &clock, LEN);
            // Block 3 is past the end and the tail: only the looping link remains.
            let work = b.build_block(3).unwrap();
            assert_eq!(work.paths.len(), 1, "{:?}", target.kind);
            assert_eq!(work.windows.headers.len(), 1);
        }
    }
}

#[test]
fn oversized_windows_are_rejected_before_they_are_built() {
    // Review P2: a 16 MS/s source needs U = 16; a 4.8 M-sample block would need a
    // 77 M-point window (2²⁶ limit). It must fail before allocating anything.
    let clock = clock();
    let n16 = SampleRate::from_hz(16_000_000).unwrap();
    let mut c = spec_of(vec![[0.5, 0.0]; 1000], EndPolicy::Loop, 1.0);
    c.rate = n16;
    c.response = ResponseSpec {
        passband_edge_uhz: mhz(6.4),
        stopband_edge_uhz: mhz(8.0),
        stopband_atten_mdb: 80_000,
        group_delay_ps: 0,
    };
    let sc = scene(CENTER_UHZ, StaticWaveform::Continuous(c));
    let b = builder(&sc, &antenna_port(LOSSLESS_CULL_LSB), &clock, 4_800_000);
    match b.build_block(0) {
        Err(SceneError::WindowTooLong { index: 0, points }) => assert!(points > 1 << 26),
        other => panic!("expected WindowTooLong, got {other:?}"),
    }
}

#[test]
fn stop_ignores_the_discarded_tail_of_a_final_partial_block() {
    // Review P2: a 1 ms run in 2 ms blocks of a 1 ms recording is covered.
    let clock = clock();
    let x = tone(2000, &[(0.5, 200_000.0)]);
    let sc = scene(CENTER_UHZ + mhz(5.0), continuous(x, EndPolicy::Stop, 0));
    for target in targets() {
        let full = builder(&sc, &target, &clock, 2 * LEN);
        assert!(matches!(
            full.build_block(0),
            Err(SceneError::SourceExhausted { .. })
        ));
        let limited = builder(&sc, &target, &clock, 2 * LEN).with_sample_limit(LEN);
        limited.build_block(0).unwrap();
    }
}
