//! Prompt 6 validation checks (plan "Validation strategy"), through the full offline
//! pipeline: scenario → static geometry → both render targets → CPU backend → receiver
//! chain or antenna-port conversion → SigMF → run manifest → replay.
//!
//! References are computed here, independently of the code under test: WGS84 → ECEF and
//! the Friis link budget from their textbook formulas, kTBF from SI constants, a
//! radix-2 FFT for tone levels and the Welch PSD. The selected tap response's |H(f)|
//! comes from `TapResponse::magnitude`, which is the definition of the response
//! (checked against the impulse response's Fourier transform in prompt 5).

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use iqsim_contracts::RunManifest;
use iqsim_core::{
    BlockSpec, ChannelTimebase, IdealClockMap, ResponseSpec, SampleIndex, SampleRate, TapKind,
    TapResponse, TargetScale, TruthTime,
};
use iqsim_io::sigmf::{decode_cf32, decode_ci16, read_recording};
use iqsim_receiver::{
    AdcConfig, ChainConfig, FrontEnd, OutputFormat, ReceiverChannelState, Samples,
};
use iqsim_render::{Backend, CpuBackend};
use iqsim_runtime::static_scene::{
    PulseSpec, StaticEmitter, StaticScene, StaticSceneBuilder, StaticTarget, StaticWaveform,
};
use iqsim_runtime::{RunError, RunOptions, replay, run_offline};
use sha2::{Digest, Sha256};

const PI: f64 = core::f64::consts::PI;
const C: f64 = 299_792_458.0;
const K_B: f64 = 1.380_649e-23;
const R0: f64 = 50.0;

// ------------------------------------------------------------------ helpers

/// x² (std `powi` is banned workspace-wide).
trait Sq {
    fn sq(self) -> f64;
}

impl Sq for f64 {
    fn sq(self) -> f64 {
        self * self
    }
}

fn scratch(tag: &str) -> PathBuf {
    static N: AtomicU32 = AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "iqsim-validation-{}-{tag}-{}",
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

/// Writes `toml` and runs it offline; returns the run directory and manifest.
fn run_toml(tag: &str, toml: &str) -> (PathBuf, RunManifest) {
    let dir = scratch(tag);
    let file = dir.join("scenario.toml");
    fs::write(&file, toml).unwrap();
    let out = dir.join("out");
    let m = run_offline(&file, &out, RunOptions { threads: 4 }).unwrap_or_else(|e| panic!("{e}"));
    (out, m)
}

/// Complex samples of a recording, scaled to the stream's full scale (1.0 = V_fs).
fn read_stream(out: &Path, sink: &str, stem: &str) -> Vec<(f64, f64)> {
    let r = read_recording(&out.join("consumer").join(sink), stem).unwrap();
    let dt = r.meta["global"]["core:datatype"]
        .as_str()
        .unwrap()
        .to_owned();
    let v: Vec<f64> = if dt == "ci16_le" {
        decode_ci16(&r.data)
            .unwrap()
            .into_iter()
            .map(|c| f64::from(c) / 32768.0)
            .collect()
    } else {
        decode_cf32(&r.data)
            .unwrap()
            .into_iter()
            .map(f64::from)
            .collect()
    };
    v.chunks(2).map(|c| (c[0], c[1])).collect()
}

/// WGS84 geodetic → ECEF, written out here from the textbook formula.
fn ecef(lat_deg: f64, lon_deg: f64, h: f64) -> [f64; 3] {
    let a = 6_378_137.0;
    let f = 1.0 / 298.257_223_563;
    let e2 = f * (2.0 - f);
    let (sp, cp) = libm::sincos(lat_deg.to_radians());
    let (sl, cl) = libm::sincos(lon_deg.to_radians());
    let n = a / (1.0 - e2 * sp * sp).sqrt();
    [
        (n + h) * cp * cl,
        (n + h) * cp * sl,
        (n * (1.0 - e2) + h) * sp,
    ]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).sq() + (a[1] - b[1]).sq() + (a[2] - b[2]).sq()).sqrt()
}

/// Receiver at 0°N 0°E, emitters at 0°N 0.05°E (≈ 5.6 km), both on the ellipsoid.
fn range_m() -> f64 {
    dist(ecef(0.0, 0.0, 0.0), ecef(0.0, 0.05, 0.0))
}

/// Friis: received power for `p_w` at `f_hz` over `range_m()`, isotropic antennas.
fn friis_w(p_w: f64, f_hz: f64) -> f64 {
    let l = C / f_hz / (4.0 * PI * range_m());
    p_w * l * l
}

fn db(x: f64) -> f64 {
    10.0 * libm::log10(x)
}

/// In-place radix-2 FFT (forward, e^{−j2πkn/N}).
fn fft(x: &mut [(f64, f64)]) {
    let n = x.len();
    assert!(n.is_power_of_two());
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            x.swap(i, j);
        }
    }
    let mut len = 2;
    while len <= n {
        let ang = -2.0 * PI / len as f64;
        for start in (0..n).step_by(len) {
            for k in 0..len / 2 {
                let (s, c) = libm::sincos(ang * k as f64);
                let (ar, ai) = x[start + k];
                let (br, bi) = x[start + k + len / 2];
                let (tr, ti) = (br * c - bi * s, br * s + bi * c);
                x[start + k] = (ar + tr, ai + ti);
                x[start + k + len / 2] = (ar - tr, ai - ti);
            }
        }
        len <<= 1;
    }
}

/// Power (in full-scale² units) of the tone in DFT bin `k` of `x` (rectangular window,
/// coherent sampling): |X[k]|²/N².
fn tone_power(x: &[(f64, f64)], k: i64) -> f64 {
    let n = x.len();
    let (mut re, mut im) = (0.0, 0.0);
    for (i, (a, b)) in x.iter().enumerate() {
        let ang = -2.0 * PI * ((k as i128 * i as i128).rem_euclid(n as i128)) as f64 / n as f64;
        let (s, c) = libm::sincos(ang);
        re += a * c - b * s;
        im += a * s + b * c;
    }
    (re * re + im * im) / (n as f64 * n as f64)
}

fn mean_power(x: &[(f64, f64)]) -> f64 {
    x.iter().map(|(a, b)| a * a + b * b).sum::<f64>() / x.len() as f64
}

fn dbm_to_w(dbm: f64) -> f64 {
    libm::pow(10.0, dbm / 10.0) * 1e-3
}

/// Full-scale voltage of a "-x dBm" full scale.
fn v_fs(dbm: f64) -> f64 {
    (dbm_to_w(dbm) * R0).sqrt()
}

/// A scenario around a 1090 MHz, 16 MS/s receiver channel and antenna port.
fn scenario(
    emitters: &[(f64, f64)],
    channel_extra: &str,
    port_extra: &str,
    duration: &str,
) -> String {
    let mut s = format!(
        r#"schema_version = 1
epoch = 2026-03-01T12:00:00Z
seed = 42
run_mode = "offline"

[metadata]
name = "validation"

[render]
duration = "{duration}"
block = "1 ms"

[[platforms]]
name = "site"
position = {{ lat = 0, lon = 0, alt = "0 m HAE" }}

[[platforms]]
name = "far"
position = {{ lat = 0, lon = 0.05, alt = "0 m HAE" }}

[[receivers]]
name = "rx"
platform = "site"

[[receivers.channels]]
name = "ch0"
center_frequency = "1090 MHz"
sample_rate = "16 MS/s"
bandwidth = "12 MHz"
sinks = ["rec"]
{channel_extra}

[[antenna_ports]]
name = "port"
platform = "site"
center_frequency = "1090 MHz"
sample_rate = "16 MS/s"
bandwidth = "14 MHz"
response = {{ attenuation = "60 dB" }}
full_scale = "-30 dBm"
sinks = ["port_rec"]
{port_extra}

[[sinks]]
name = "rec"
kind = "sigmf"

[[sinks]]
name = "port_rec"
kind = "sigmf"
format = "cf32_le"
"#
    );
    for (i, (f_mhz, p_dbm)) in emitters.iter().enumerate() {
        s.push_str(&format!(
            "\n[[emitters]]\nname = \"e{i}\"\nplatform = \"far\"\nfrequency = \"{f_mhz} MHz\"\npower = \"{p_dbm} dBm\"\nbehavior = \"always_on\"\nwaveform = \"cw\"\n"
        ));
    }
    s
}

fn rate16() -> SampleRate {
    SampleRate::from_hz(16_000_000).unwrap()
}

/// The channel's H_rx and the port's B_out as resolved from `scenario` (bandwidth/2 and
/// rate/2 edges, 80 dB and 60 dB).
fn responses() -> (TapResponse, TapResponse) {
    let mk = |kind, fp, a| {
        TapResponse::new(
            kind,
            ResponseSpec {
                passband_edge_uhz: fp,
                stopband_edge_uhz: 8_000_000_000_000,
                stopband_atten_mdb: a,
                group_delay_ps: 0,
            },
            rate16(),
            1_090_000_000_000_000,
        )
        .unwrap()
    };
    (
        mk(TapKind::VirtualReceiver, 6_000_000_000_000, 80_000),
        mk(TapKind::AntennaPort, 7_000_000_000_000, 60_000),
    )
}

// ------------------------------------------------------------------ checks

/// Power and in-band level: every tone lands at its frequency with the closed-form
/// link-budget power times |H(f)|², within 0.1 dB, on both taps.
#[test]
fn tone_levels_match_the_link_budget_and_selected_response_on_both_taps() {
    // Offsets on exact 1 kHz bins of a 16 000-sample (1 ms) record.
    let offsets_mhz = [0.0, 1.5, -3.25, 5.5, -6.6, 7.2];
    let emitters: Vec<(f64, f64)> = offsets_mhz.iter().map(|o| (1090.0 + o, 30.0)).collect();
    let (out, _) = run_toml(
        "levels",
        &scenario(
            &emitters,
            "front_end = { gain = \"20 dB\" }\nthermal_noise = false\nadc = { full_scale = \"-40 dBm\" }",
            "",
            "1 ms",
        ),
        // Power at the receiver input ≈ −58 dBm per tone; +20 dB gain against a −40 dBm
        // ADC full scale leaves ≈ 3 dB of crest headroom over six tones' sum.
    );
    let (h_rx, b_out) = responses();
    let rx = read_stream(&out, "rec", "rx.ch0");
    let port = read_stream(&out, "port_rec", "port");
    assert_eq!(rx.len(), 16_000);
    let gain = 100.0; // 20 dB
    for (o, (f_mhz, p_dbm)) in offsets_mhz.iter().zip(&emitters) {
        let k = (o * 1000.0).round() as i64;
        let budget = friis_w(dbm_to_w(*p_dbm), f_mhz * 1e6);
        let f_off = o * 1e6;
        // Virtual receiver: code/32768 = v/V_adc at the ADC input.
        let v = v_fs(-40.0);
        let measured = tone_power(&rx, k) * v * v / R0 / gain;
        let expected = budget * h_rx.magnitude(f_off).sq();
        let err = db(measured) - db(expected);
        eprintln!(
            "rx   {o:+5} MHz: {:.3} dBm, error {err:+.4} dB",
            db(measured * 1e3)
        );
        assert!(err.abs() < 0.1, "rx {o} MHz: {err} dB");
        // Antenna port: cf32 = v/V_fs.
        let v = v_fs(-30.0);
        let measured = tone_power(&port, k) * v * v / R0;
        let expected = budget * b_out.magnitude(f_off).sq();
        let err = db(measured) - db(expected);
        eprintln!(
            "port {o:+5} MHz: {:.3} dBm, error {err:+.4} dB",
            db(measured * 1e3)
        );
        assert!(err.abs() < 0.1, "port {o} MHz: {err} dB");
    }
}

/// Out-of-band rejection: a 1030 MHz interrogation into the 1090 MHz, 16 MS/s channel
/// would fold to +4 MHz. Nothing above the selected stopband appears there, while an
/// in-band reference tone keeps its level.
#[test]
fn a_1030_mhz_emitter_leaves_no_aliased_image() {
    let (out, _) = run_toml(
        "oob",
        &scenario(
            &[(1030.0, 60.0), (1091.0, 0.0)],
            "front_end = { gain = \"20 dB\" }\nthermal_noise = false\nadc = { full_scale = \"-40 dBm\" }",
            "",
            "1 ms",
        ),
    );
    let (h_rx, b_out) = responses();
    let alias_bin = 4000; // −60 MHz mod 16 MHz = +4 MHz
    let strong = friis_w(dbm_to_w(60.0), 1030e6);
    for (stream, sink, scale, gain, resp, atten_db) in [
        ("rx.ch0", "rec", v_fs(-40.0), 100.0, h_rx, 80.0),
        ("port", "port_rec", v_fs(-30.0), 1.0, b_out, 60.0),
    ] {
        let x = read_stream(&out, sink, stream);
        let image = tone_power(&x, alias_bin) * scale * scale / R0 / gain;
        eprintln!(
            "{stream}: image {:.1} dB below the unfiltered 1030 MHz level",
            db(strong) - db(image.max(1e-300))
        );
        assert!(
            image <= strong * libm::pow(10.0, -atten_db / 10.0),
            "{stream}: aliased image"
        );
        let reference = tone_power(&x, 1000) * scale * scale / R0 / gain;
        let expected = friis_w(dbm_to_w(0.0), 1091e6) * resp.magnitude(1e6).sq();
        assert!(
            (db(reference) - db(expected)).abs() < 0.1,
            "{stream}: in-band level"
        );
    }
}

/// Noise floor: with no emitters, the channel's noise is kT₀·F·G over the noise
/// bandwidth ∫|H_rx|² df (total power) and kT₀·F·G·|H(f)|² per hertz (Welch PSD in the
/// passband), each within 0.2 dB.
#[test]
fn noise_floor_matches_ktbf_over_the_noise_bandwidth() {
    // 2¹⁷ samples = 8.192 ms at 16 MS/s.
    let (out, _) = run_toml(
        "noise",
        &scenario(
            &[],
            "front_end = { gain = \"40 dB\", noise_figure = \"6 dB\" }\nadc = { full_scale = \"-40 dBm\" }",
            "",
            "8.192 ms",
        ),
    );
    let (h_rx, _) = responses();
    let x = read_stream(&out, "rec", "rx.ch0");
    assert_eq!(x.len(), 1 << 17);
    let v = v_fs(-40.0);
    let g = 1e4;
    let f = libm::pow(10.0, 0.6);
    let n0 = K_B * 290.0 * f * g; // W/Hz at the ADC input
    // Noise bandwidth by the trapezoid rule over ±F_s/2 (|H| ≤ 10⁻⁴ beyond).
    let steps = 160_000;
    let df = 16e6 / steps as f64;
    let mut bn = 0.0;
    for i in 0..=steps {
        let ff = -8e6 + i as f64 * df;
        let w = if i == 0 || i == steps { 0.5 } else { 1.0 };
        bn += w * h_rx.magnitude(ff).sq() * df;
    }
    let measured = mean_power(&x) * v * v / R0;
    let err = db(measured) - db(n0 * bn);
    eprintln!("noise power error {err:+.4} dB (B_n = {:.4} MHz)", bn / 1e6);
    assert!(err.abs() < 0.2, "total noise power {err} dB");
    // Welch: 1024-point Hann segments, 50 % overlap.
    let seg = 1024;
    let win: Vec<f64> = (0..seg)
        .map(|i| 0.5 - 0.5 * libm::cos(2.0 * PI * i as f64 / seg as f64))
        .collect();
    let u: f64 = win.iter().map(|w| w * w).sum();
    let mut psd = vec![0.0; seg];
    let mut count = 0;
    let mut start = 0;
    while start + seg <= x.len() {
        let mut buf: Vec<(f64, f64)> = (0..seg)
            .map(|i| (x[start + i].0 * win[i], x[start + i].1 * win[i]))
            .collect();
        fft(&mut buf);
        for (p, b) in psd.iter_mut().zip(&buf) {
            *p += b.0 * b.0 + b.1 * b.1;
        }
        count += 1;
        start += seg / 2;
    }
    let mut worst: f64 = 0.0;
    let mut band = (0.0, 0.0);
    for (k, p) in psd.iter().enumerate() {
        let ff = if k < seg / 2 { k } else { k as i64 as usize } as f64;
        let ff = if k < seg / 2 { ff } else { ff - seg as f64 } * 16e6 / seg as f64;
        if ff.abs() > 4.5e6 {
            continue;
        }
        let measured = p / count as f64 / (16e6 * u) * v * v / R0;
        let expected = n0 * h_rx.magnitude(ff).sq();
        band.0 += measured;
        band.1 += expected;
        worst = worst.max((db(measured) - db(expected)).abs());
    }
    let err = db(band.0) - db(band.1);
    eprintln!("Welch passband PSD error {err:+.4} dB (worst single bin {worst:.2} dB)");
    assert!(err.abs() < 0.2, "passband PSD {err} dB");
}

/// Delay: a pulse's time of arrival through the full virtual-receiver chain (gain and
/// 16-bit ADC) is within 0.001 samples of truth after removing H_rx's group delay.
#[test]
fn pulse_time_of_arrival_survives_the_receiver_chain() {
    let fs = 16e6;
    let rate = rate16();
    let resp = ResponseSpec {
        passband_edge_uhz: 6_000_000_000_000,
        stopband_edge_uhz: 7_500_000_000_000,
        stopband_atten_mdb: 80_000,
        group_delay_ps: 1_000_000,
    };
    let rx = [6_378_137.0, 0.0, 0.0];
    let em = [6_378_137.0, 30_000.0, 0.0];
    let tau = dist(rx, em) / C;
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::from_ps(0), rate));
    let v_adc = 1e-3;
    let target = StaticTarget {
        kind: TapKind::VirtualReceiver,
        center_uhz: 1_090_000_000_000_000,
        rate,
        response: resp,
        full_scale_volts: v_adc / 10.0,
        cull_threshold_lsb: iqsim_core::amplitude::LOSSLESS_CULL_LSB,
    };
    for frac in [0.0, 0.3, 0.77] {
        let emission_ps = 300_000_000 + (frac * 62_500.0) as i64;
        let scene = StaticScene {
            receiver: iqsim_geo::Ecef::new(rx[0], rx[1], rx[2]),
            emitters: vec![StaticEmitter {
                position: iqsim_geo::Ecef::new(em[0], em[1], em[2]),
                carrier_uhz: 1_090_000_000_000_000,
                power_w: 100.0,
                waveform: StaticWaveform::Pulses {
                    emission: TruthTime::from_ps(emission_ps),
                    pulses: vec![PulseSpec {
                        start_ps: 0,
                        rise_ps: 100_000,
                        fall_start_ps: 450_000,
                        fall_ps: 100_000,
                        amplitude: 1.0,
                    }],
                },
            }],
        };
        let builder =
            StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(16_384).unwrap())
                .unwrap();
        let work = builder.build_block(0).unwrap();
        let mut acc = vec![0_i64; 2 * work.len as usize];
        CpuBackend::new(4).render(&work, &mut acc).unwrap();
        let chain = ChainConfig {
            rate,
            response: *builder.response(),
            input_scale: TargetScale::new(target.full_scale_volts).unwrap(),
            adc_full_scale_volts: v_adc,
            output: OutputFormat::Cf32,
            front_end: Some(FrontEnd {
                gain_db: 20.0,
                noise_figure_db: 3.0,
            }),
            noise: None,
            lo: None,
            iq: None,
            decimation: None,
            agc: None,
            adc: Some(AdcConfig { bits: 16 }),
        };
        let mut state = ReceiverChannelState::new(chain).unwrap();
        let o = state.process(0, SampleIndex::new(0), &acc, &clock).unwrap();
        let Samples::Cf32(y) = o.samples else {
            panic!()
        };
        assert_eq!(o.clipped, 0);
        let (mut num, mut den) = (0.0, 0.0);
        for (k, c) in y.chunks(2).enumerate() {
            let e = f64::from(c[0]).sq() + f64::from(c[1]).sq();
            num += k as f64 * e;
            den += e;
        }
        let measured = num / den;
        let expected = (emission_ps as f64 * 1e-12 + tau + 1e-6 + 275e-9) * fs;
        let err = measured - expected;
        eprintln!("TOA frac {frac}: error {err:.2e} samples");
        assert!(err.abs() < 1e-3, "TOA error {err} samples");
    }
}

/// Tap isolation: changing H_rx, gain, noise, clock, LO, ADC, IQ, decimation, AGC or
/// the seed changes the virtual-receiver output but never the antenna-port output.
#[test]
fn antenna_port_hash_ignores_every_virtual_receiver_setting() {
    let base_channel = "front_end = { gain = \"20 dB\", noise_figure = \"5 dB\" }\nadc = { full_scale = \"-40 dBm\" }";
    let emitters = [(1090.0, 30.0), (1093.0, 20.0)];
    let hashes = |tag: &str, channel: &str, seed: Option<u64>| {
        let mut s = scenario(&emitters, channel, "", "1 ms");
        if let Some(seed) = seed {
            s = s.replace("seed = 42", &format!("seed = {seed}"));
        }
        let (_, m) = run_toml(tag, &s);
        let get = |id: &str| {
            m.outputs
                .iter()
                .find(|o| o.stream_id == id)
                .unwrap()
                .sha256
                .clone()
        };
        (get("rx.ch0"), get("port"))
    };
    let (rx0, port0) = hashes("iso-base", base_channel, None);
    let variants = [
        (
            "bandwidth",
            base_channel.to_owned() + "\nresponse = { stopband = \"15 MHz\" }",
        ),
        ("gain", base_channel.replace("20 dB", "25 dB")),
        ("nf", base_channel.replace("5 dB", "8 dB")),
        (
            "noise-off",
            base_channel.to_owned() + "\nthermal_noise = false",
        ),
        (
            "clock",
            base_channel.to_owned() + "\nclock = { model = \"fixed_offset\", offset = \"13 ns\" }",
        ),
        (
            "lo",
            base_channel.to_owned() + "\nlo = { model = \"fixed_offset\", offset = \"1.5 kHz\" }",
        ),
        (
            "adc",
            base_channel.replace("\"-40 dBm\" }", "\"-40 dBm\", bits = 10 }"),
        ),
        (
            "iq",
            base_channel.to_owned()
                + "\niq = { gain_imbalance = \"0.5 dB\", phase_skew = \"2 deg\", dc_i = 0.01 }",
        ),
        (
            "decimation",
            base_channel.to_owned() + "\ndecimation = { factor = 2 }",
        ),
        (
            "agc",
            base_channel.to_owned() + "\nagc = { target = \"-20 dB\" }",
        ),
    ];
    for (tag, channel) in &variants {
        let (rx, port) = hashes(tag, channel, None);
        assert_eq!(port, port0, "{tag}: antenna-port output changed");
        assert_ne!(rx, rx0, "{tag}: virtual-receiver output did not change");
    }
    let (rx, port) = hashes("seed", base_channel, Some(43));
    assert_eq!(port, port0, "seed: antenna-port output changed");
    assert_ne!(rx, rx0, "seed: noise did not change");
    // Changing the antenna port's own configuration does change it.
    let base = scenario(&emitters, base_channel, "", "1 ms");
    for (tag, s) in [
        (
            "port-fs",
            base.replace("full_scale = \"-30 dBm\"", "full_scale = \"-31 dBm\""),
        ),
        (
            "port-gd",
            base.replace(
                "attenuation = \"60 dB\" }",
                "attenuation = \"60 dB\", group_delay = \"100 ns\" }",
            ),
        ),
    ] {
        let (_, m) = run_toml(tag, &s);
        assert_ne!(m.outputs[1].sha256, port0, "{tag}");
    }
}

/// Golden hashes for examples/minimal.toml: both recordings and the truth record. The
/// manifest's hashes match the files, and consumer files carry no truth or manifest.
#[test]
fn minimal_example_golden_hashes() {
    let dir = scratch("minimal");
    let out = dir.join("out");
    let m = run_offline(
        &workspace_file("examples/minimal.toml"),
        &out,
        RunOptions { threads: 3 },
    )
    .unwrap();
    assert_eq!(m.outputs.len(), 2);
    for o in &m.outputs {
        let bytes = fs::read(out.join(&o.path)).unwrap();
        assert_eq!(hex::encode(Sha256::digest(&bytes)), o.sha256);
        assert_eq!(bytes.len() as u64, 4 * o.sample_count);
        assert_eq!(o.adc_clipped_components, 0);
        iqsim_testkit::assert_golden(
            &format!("runtime/minimal_{}", o.stream_id.replace('.', "_")),
            &bytes,
        );
    }
    let truth = fs::read(out.join("run/truth.json")).unwrap();
    assert_eq!(hex::encode(Sha256::digest(&truth)), m.truth.sha256);
    iqsim_testkit::assert_golden("runtime/minimal_truth", &truth);
    // Nothing privileged under consumer/.
    for sink in fs::read_dir(out.join("consumer")).unwrap() {
        for f in fs::read_dir(sink.unwrap().path()).unwrap() {
            let p = f.unwrap().path();
            let name = p.file_name().unwrap().to_string_lossy().into_owned();
            assert!(
                name.ends_with(".sigmf-data") || name.ends_with(".sigmf-meta"),
                "{name}"
            );
            if name.ends_with(".sigmf-meta") {
                let text = fs::read_to_string(&p).unwrap();
                for word in ["truth", "manifest", "seed", "range", "delay", "emitter"] {
                    assert!(!text.contains(word), "{name} mentions {word}");
                }
            }
        }
    }
}

/// Replay re-renders identically (with a different thread count) and rejects every
/// changed input before rendering.
#[test]
fn replay_matches_and_rejects_changed_inputs() {
    let dir = scratch("replay");
    let file = dir.join("minimal.toml");
    let original = fs::read_to_string(workspace_file("examples/minimal.toml")).unwrap();
    fs::write(&file, &original).unwrap();
    let out = dir.join("out");
    let m = run_offline(&file, &out, RunOptions { threads: 4 }).unwrap();
    let manifest = out.join("run/manifest.json");
    let report = replay(&manifest, RunOptions { threads: 1 }).unwrap();
    assert_eq!(report.outputs_matched, 2);
    assert_eq!(report.deterministic_sha256, m.deterministic_sha256());

    let expect = |what: &str, needle: &str| match replay(&manifest, RunOptions { threads: 2 }) {
        Err(RunError::Replay(msg)) => assert!(msg.contains(needle), "{what}: {msg}"),
        other => panic!("{what}: expected a replay rejection, got {other:?}"),
    };
    // The scenario file changed (even a comment).
    fs::write(&file, original.clone() + "\n# edited\n").unwrap();
    expect("scenario edit", "changed since the run");
    fs::remove_file(&file).unwrap();
    expect("scenario removed", "cannot be read");
    fs::write(&file, &original).unwrap();
    replay(&manifest, RunOptions { threads: 2 }).unwrap();

    let text = fs::read_to_string(&manifest).unwrap();
    let edit = |f: &dyn Fn(&mut RunManifest)| {
        let mut m = RunManifest::from_json(&text).unwrap();
        f(&mut m);
        fs::write(&manifest, m.to_json()).unwrap();
    };
    edit(&|m| m.numeric_contract_version += 1);
    expect("numeric contract", "numeric contract version");
    edit(&|m| m.build.cargo_lock_sha256 = "00".repeat(32));
    expect("build", "different build");
    // A source edit with the same version and lockfile is a different build too.
    edit(&|m| m.build.source_sha256 = "00".repeat(32));
    expect("source fingerprint", "different build");
    assert_eq!(m.build.source_sha256.len(), 64);
    edit(&|m| m.scenario.resolved["seed"] = 7.into());
    expect("resolved scenario", "does not match its hash");
    edit(&|m| {
        m.scenario.resolved["seed"] = 7.into();
        m.scenario.resolved_sha256 = iqsim_contracts::canonical_sha256(&m.scenario.resolved);
    });
    expect("resolved scenario rehashed", "resolves differently");
    let asset = dir.join("pattern.bin");
    fs::write(&asset, b"pattern v1").unwrap();
    let mut reg = iqsim_contracts::AssetRegistry::new();
    reg.register(
        iqsim_contracts::AssetKind::AntennaPattern,
        "mast",
        asset.to_str().unwrap(),
        b"pattern v1",
        "raw",
        None,
    )
    .unwrap();
    let records = reg.records();
    edit(&|m| m.assets = records.clone());
    // Since prompt 7 the asset list must be exactly what the scenario resolves to, so
    // an extra record is rejected; file verification still comes first.
    expect("extra asset", "asset list differs");
    fs::write(&asset, b"pattern v2").unwrap();
    expect("asset changed", "changed since the run was recorded");
    fs::remove_file(&asset).unwrap();
    expect("asset missing", "is missing");
    edit(&|m| m.outputs[0].sha256 = "11".repeat(32));
    expect("output hash", "re-rendered differently");
    edit(&|m| m.truth.sha256 = "22".repeat(32));
    expect("truth hash", "truth record differs");
    edit(&|m| {
        m.event_log = Some(iqsim_contracts::EventLogRecord {
            path: "log".into(),
            sha256: String::new(),
            record_count: 0,
            first_sequence: 0,
            last_sequence: 0,
            final_commit_ps: 0,
        })
    });
    expect("event log", "event-log replay");
    // Unknown manifest fields are rejected.
    fs::write(
        &manifest,
        text.replace("\"seed\":", "\"extra\": 1,\n  \"seed\":"),
    )
    .unwrap();
    assert!(replay(&manifest, RunOptions { threads: 1 }).is_err());
}

/// A run refuses to write into a non-empty directory and reports unsupported scenario
/// features instead of guessing.
#[test]
fn runs_reject_unsafe_outputs_and_unsupported_features() {
    let dir = scratch("reject");
    let file = dir.join("s.toml");
    let base = scenario(&[(1090.0, 0.0)], "", "", "1 ms");
    fs::write(&file, &base).unwrap();
    let out = dir.join("out");
    fs::create_dir_all(&out).unwrap();
    fs::write(out.join("keep"), b"x").unwrap();
    assert!(matches!(
        run_offline(&file, &out, RunOptions::default()),
        Err(RunError::Io(_))
    ));
    for (edit, needle) in [
        (
            base.replace("[render]\nduration = \"1 ms\"\nblock = \"1 ms\"\n", ""),
            "[render]",
        ),
        (
            base.replace("waveform = \"cw\"", "waveform = \"mode_s\""),
            "built-in waveform",
        ),
        (
            base.replace(
                "kind = \"sigmf\"\n\n[[sinks]]",
                "kind = \"vita49\"\n\n[[sinks]]",
            ),
            "only SigMF",
        ),
        (
            base.replace("block = \"1 ms\"", "block = \"1.00001 ms\""),
            "whole number of samples",
        ),
        // Flight levels resolve since prompt 7; a missing geoid grid does not.
        (
            format!("{base}\n[world]\ngeoid = {{ path = \"missing-egm96-5.pgm\" }}\n"),
            "world.geoid",
        ),
    ] {
        fs::write(&file, &edit).unwrap();
        let out = dir.join(format!("o-{}", needle.len()));
        match run_offline(&file, &out, RunOptions::default()) {
            Err(e) => assert!(e.to_string().contains(needle), "{needle}: {e}"),
            Ok(_) => panic!("{needle}: accepted"),
        }
    }
}

/// Truth maps each stream's indices to time (review fix): decimation changes the
/// output count and rate and reports its delay; a clock offset moves the first
/// sample's acquisition instant and reported label; the antenna port keeps its own.
#[test]
fn truth_records_each_streams_index_and_time_mapping() {
    let channel = "decimation = { factor = 2 }\nclock = { model = \"fixed_offset\", offset = \"13 ns\" }\nthermal_noise = false";
    let (out, m) = run_toml(
        "truth-map",
        &scenario(&[(1090.0, 0.0)], channel, "", "1 ms"),
    );
    let truth: serde_json::Value =
        serde_json::from_slice(&fs::read(out.join("run/truth.json")).unwrap()).unwrap();
    let rx = &truth["streams"][0];
    assert_eq!(rx["stream"], "rx.ch0");
    assert_eq!(rx["render_sample_count"], 16_000);
    assert_eq!(rx["output_sample_count"], 8_000);
    assert_eq!(rx["output_rate"], serde_json::json!([8_000_000, 1]));
    assert_eq!(rx["output_index_to_render_index_factor"], 2);
    assert!(rx["decimation_delay_render_samples"].as_u64().unwrap() > 0);
    assert_eq!(
        rx["first_sample_acquisition"],
        serde_json::json!([13_000, 0])
    );
    assert_eq!(rx["first_sample_reported"], serde_json::json!([13_000, 0]));
    let port = &truth["streams"][1];
    assert_eq!(port["output_sample_count"], 16_000);
    assert_eq!(port["decimation_delay_render_samples"], 0);
    assert_eq!(port["first_sample_acquisition"], serde_json::json!([0, 0]));
    let fs_eff = port["effective_full_scale_volts"].as_f64().unwrap();
    assert!((fs_eff / v_fs(-30.0) - 1.0).abs() < 1e-7);
    // The recording agrees with truth.
    assert_eq!(m.outputs[0].sample_count, 8_000);
    assert_eq!(
        (m.outputs[0].rate_num, m.outputs[0].rate_den),
        (8_000_000, 1)
    );
}
