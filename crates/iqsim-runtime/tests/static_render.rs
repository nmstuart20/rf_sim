//! Prompt 5 acceptance tests: static geometry → work list → CPU render.
//!
//! Reference values here are computed independently of the builder and renderer:
//! Friis and light time from the test's own geometry, the response's defining
//! formula, and a direct quadrature of the band-limited pulse (no edge tables).

mod common;

use iqsim_core::amplitude::LOSSLESS_CULL_LSB;
use iqsim_core::{
    BlockSpec, ChannelTimebase, IdealClockMap, ResponseSpec, SampleRate, TapKind, TapResponse,
    TruthTime,
};
use iqsim_geo::Ecef;
use iqsim_render::{Backend, CpuBackend, RenderError};
use iqsim_runtime::scene::{
    PulseSpec, SceneError, StaticEmitter, StaticScene, StaticSceneBuilder, StaticTarget,
    StaticWaveform,
};

const FS: f64 = 16e6;
const LEN: u64 = 16_384;
const C: f64 = 299_792_458.0;
const CENTER_UHZ: i64 = 1_090_000_000_000_000;
const RANGE_M: f64 = 30_000.0;
const POWER_W: f64 = 1000.0;
const PI: f64 = core::f64::consts::PI;

fn rate() -> SampleRate {
    SampleRate::from_hz(16_000_000).unwrap()
}

fn clock() -> IdealClockMap {
    IdealClockMap::new(ChannelTimebase::new(TruthTime::from_ps(0), rate()))
}

fn mhz(x: f64) -> i64 {
    (x * 1e12) as i64
}

/// H_rx: 6 / 7.5 MHz, 80 dB, 1 µs group delay.
fn virtual_rx() -> StaticTarget {
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
        cull_threshold_lsb: LOSSLESS_CULL_LSB,
    }
}

/// B_out: 7 / 8 MHz, 60 dB, zero group delay, its own full scale.
fn antenna_port() -> StaticTarget {
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
        cull_threshold_lsb: LOSSLESS_CULL_LSB,
    }
}

/// The light time of the test link in ECI, Sagnac included (independent reference).
fn link_tau() -> f64 {
    common::eci_delay([6_378_137.0, 0.0, 0.0], [6_378_137.0, RANGE_M, 0.0])
}

fn receiver() -> Ecef {
    Ecef::new(6_378_137.0, 0.0, 0.0)
}

fn emitter(carrier_uhz: i64, waveform: StaticWaveform) -> StaticEmitter {
    StaticEmitter {
        position: Ecef::new(6_378_137.0, RANGE_M, 0.0),
        carrier_uhz,
        power_w: POWER_W,
        waveform,
    }
}

/// Independent Friis voltage at the receiver.
fn friis_volts(carrier_hz: f64) -> f64 {
    (POWER_W * 50.0).sqrt() * (C / carrier_hz) / (4.0 * PI * (C * link_tau()))
}

fn trapezoid(emission_ps: i64) -> StaticWaveform {
    StaticWaveform::Pulses {
        emission: TruthTime::from_ps(emission_ps),
        pulses: vec![PulseSpec {
            start_ps: 0,
            rise_ps: 100_000,
            fall_start_ps: 450_000,
            fall_ps: 100_000,
            amplitude: 1.0,
        }],
    }
}

/// Renders block `b` and returns complex volts.
fn render_volts(builder: &StaticSceneBuilder<'_>, b: u64) -> Vec<(f64, f64)> {
    let acc = render_acc(builder, b, CpuBackend::new(4));
    let q = builder.scale().lsb_volts();
    acc.chunks(2)
        .map(|c| (c[0] as f64 * q, c[1] as f64 * q))
        .collect()
}

fn render_acc(builder: &StaticSceneBuilder<'_>, b: u64, backend: CpuBackend) -> Vec<i64> {
    let work = builder.build_block(b).unwrap();
    let mut acc = vec![0_i64; 2 * work.len as usize];
    backend.render(&work, &mut acc).unwrap();
    acc
}

#[test]
fn pulse_time_of_arrival_is_within_a_thousandth_of_a_sample() {
    let tau = link_tau();
    let clock = clock();
    for target in [virtual_rx(), antenna_port()] {
        let gd = target.response.group_delay_ps as f64 * 1e-12;
        for frac in [0.0, 0.1, 0.25, 0.5, 0.731, 0.999] {
            let emission_ps = 300_000_000 + (frac * 62_500.0) as i64;
            let scene = StaticScene {
                receiver: receiver(),
                emitters: vec![emitter(CENTER_UHZ, trapezoid(emission_ps))],
            };
            let builder =
                StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(LEN).unwrap())
                    .unwrap();
            let y = render_volts(&builder, 0);
            // Energy centroid: exact for a band-limited signal sampled above twice its
            // bandwidth. The pulse is symmetric about the trapezoid's centre (275 ns).
            let (mut num, mut den) = (0.0, 0.0);
            for (k, (i, q)) in y.iter().enumerate() {
                let e = i * i + q * q;
                num += k as f64 * e;
                den += e;
            }
            let measured = num / den;
            let expected = (emission_ps as f64 * 1e-12 + tau + gd + 275e-9) * FS;
            let err = measured - expected;
            eprintln!("TOA {:?} frac {frac}: error {err:.3e} samples", target.kind);
            assert!(
                err.abs() < 1e-3,
                "{:?} frac {frac}: TOA error {err:.2e} samples",
                target.kind
            );
        }
    }
}

/// Direct quadrature of ∫ e(τ − x)·h(x)·e^{−j2πf x} dx (composite Simpson over ±6 µs
/// around the group delay, 0.25 ns panels): an untabulated, untruncated reference.
fn reference_edge(resp: &TapResponse, f_off: f64, rise_s: f64, tau: f64) -> (f64, f64) {
    let gd = resp.group_delay_s();
    let fc = resp.cutoff_hz();
    let a = 2.0 * PI * PI * resp.sigma_hz() * resp.sigma_hz();
    let half = 6e-6;
    let n = 48_000;
    let dx = 2.0 * half / f64::from(n);
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..=n {
        let x = gd - half + f64::from(i) * dx;
        let u = tau - x;
        let e = if u <= 0.0 {
            0.0
        } else if rise_s == 0.0 || u >= rise_s {
            1.0
        } else {
            u / rise_s
        };
        if e == 0.0 {
            continue;
        }
        let t = x - gd;
        let arg = PI * 2.0 * fc * t;
        let sinc = if arg == 0.0 {
            1.0
        } else {
            libm::sin(arg) / arg
        };
        let h = 2.0 * fc * sinc * libm::exp(-(a * t * t));
        let w = if i == 0 || i == n {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let (s, c) = libm::sincos(-2.0 * PI * f_off * x);
        re += w * e * h * c;
        im += w * e * h * s;
    }
    (re * dx / 3.0, im * dx / 3.0)
}

#[test]
fn pulses_match_an_independent_band_limited_reference() {
    let tau = link_tau();
    let clock = clock();
    let target = virtual_rx();
    // In band (+3 MHz), and 1030 MHz interrogations 60 MHz below the channel: what
    // remains is the pulse's own spectral splatter inside H_rx, never an alias.
    for carrier_uhz in [CENTER_UHZ + mhz(3.0), mhz(1030.0)] {
        let emission_ps = 200_000_123;
        let scene = StaticScene {
            receiver: receiver(),
            emitters: vec![emitter(carrier_uhz, trapezoid(emission_ps))],
        };
        let builder =
            StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(LEN).unwrap()).unwrap();
        for h in &builder.build_block(0).unwrap().tables.headers {
            eprintln!("table: 2^{} points/sample, {} points", h.log2_res, h.len);
        }
        let y = render_volts(&builder, 0);
        let resp = builder.response();
        let carrier_hz = carrier_uhz as f64 * 1e-6;
        let f_off = (carrier_uhz - CENTER_UHZ) as f64 * 1e-6;
        let v = friis_volts(carrier_hz);
        let t_e = emission_ps as f64 * 1e-12;
        let (mut worst, mut peak) = (0.0_f64, 0.0_f64);
        let first = ((t_e + tau) * FS) as usize - 100;
        for (k, yk) in y.iter().enumerate().skip(first).take(260) {
            let t = k as f64 / FS;
            let r = reference_edge(resp, f_off, 100e-9, t - tau - t_e);
            let f = reference_edge(resp, f_off, 100e-9, t - tau - t_e - 450e-9);
            let (gr, gi) = (r.0 - f.0, r.1 - f.1);
            let phi = 2.0 * PI * (f_off * t - carrier_hz * tau);
            let (s, c) = libm::sincos(phi);
            let er = v * (gr * c - gi * s);
            let ei = v * (gr * s + gi * c);
            let d = ((yk.0 - er) * (yk.0 - er) + (yk.1 - ei) * (yk.1 - ei)).sqrt();
            worst = worst.max(d);
            peak = peak.max((er * er + ei * ei).sqrt());
        }
        eprintln!(
            "pulse offset {f_off}: worst {:.3e}·v, peak {:.3e}·v",
            worst / v,
            peak / v
        );
        // Outside the window the output must be (numerically) zero.
        let quiet = y[..first].iter().chain(&y[first + 260..]);
        let stray = quiet
            .map(|(i, q)| (i * i + q * q).sqrt())
            .fold(0.0, f64::max);
        assert!(
            worst <= 1e-5 * v,
            "offset {f_off}: max error {:.2e} V^-1 relative (peak {:.2e})",
            worst / v,
            peak / v
        );
        assert!(
            stray <= 1e-5 * v,
            "offset {f_off}: stray output {stray:.2e} V"
        );
        assert!(peak > 0.0);
    }
}

#[test]
fn cw_tone_has_the_expected_frequency_level_and_phase_for_each_tap() {
    let tau = link_tau();
    let clock = clock();
    for target in [virtual_rx(), antenna_port()] {
        for off_mhz in [0.0, 2.5, -4.0, 6.9, -7.2] {
            let carrier_uhz = CENTER_UHZ + mhz(off_mhz);
            let scene = StaticScene {
                receiver: receiver(),
                emitters: vec![emitter(carrier_uhz, StaticWaveform::Cw)],
            };
            let builder =
                StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(LEN).unwrap())
                    .unwrap();
            let y = render_volts(&builder, 1);
            let f = off_mhz * 1e6;
            let resp = builder.response();
            let expected = friis_volts(carrier_uhz as f64 * 1e-6) * resp.magnitude(f);
            let n0 = LEN as f64;
            let (mut sr, mut si, mut pr, mut pi) = (0.0, 0.0, 0.0, 0.0);
            for (k, yk) in y.iter().enumerate() {
                let t = (n0 + k as f64) / FS;
                let (s, c) = libm::sincos(-2.0 * PI * f * t);
                sr += yk.0 * c - yk.1 * s;
                si += yk.0 * s + yk.1 * c;
                if k > 0 {
                    let p = y[k - 1];
                    pr += yk.0 * p.0 + yk.1 * p.1;
                    pi += yk.1 * p.0 - yk.0 * p.1;
                }
            }
            let level = (sr * sr + si * si).sqrt() / LEN as f64;
            let rel = level / expected - 1.0;
            eprintln!(
                "CW {:?} {off_mhz} MHz: |H| {:.6e} level err {rel:.2e}",
                target.kind,
                resp.magnitude(f)
            );
            // 1e-5 relative is 8.7e-5 dB, far inside the 0.01 dB contract tolerance.
            assert!(
                rel.abs() < 1e-5,
                "{:?} {off_mhz} MHz: level error {rel:.2e}",
                target.kind
            );
            let f_est = libm::atan2(pi, pr) * FS / (2.0 * PI);
            assert!((f_est - f).abs() < 1.0, "{off_mhz} MHz: frequency {f_est}");
            // Phase: −f_c·τ − f·τ_g, mod one cycle.
            let gd = target.response.group_delay_ps as f64 * 1e-12;
            let want = -(carrier_uhz as f64 * 1e-6 * tau) - f * gd;
            let got = libm::atan2(si, sr) / (2.0 * PI);
            let mut d = got - want;
            d -= libm::round(d);
            assert!(
                d.abs() * 2.0 * PI < 1e-4,
                "{off_mhz} MHz: phase error {d} cycles"
            );
        }
    }
}

#[test]
fn out_of_band_emitters_leave_no_alias_above_the_stopband() {
    let clock = clock();
    for culled in [true, false] {
        let mut target = virtual_rx();
        if !culled {
            target.cull_threshold_lsb = 0.0;
        }
        // 1030 MHz is 60 MHz below the 1090 MHz channel: at 16 MS/s it would fold to
        // −4 MHz. 1099 MHz is beyond Nyquist and would fold to −7 MHz.
        for (carrier_uhz, alias_hz) in [(mhz(1030.0), -4e6), (mhz(1099.0), -7e6)] {
            let scene = StaticScene {
                receiver: receiver(),
                emitters: vec![emitter(carrier_uhz, StaticWaveform::Cw)],
            };
            let builder =
                StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(LEN).unwrap())
                    .unwrap();
            let work = builder.build_block(0).unwrap();
            assert_eq!(work.paths.is_empty(), culled);
            let y = render_volts(&builder, 0);
            let (mut sr, mut si) = (0.0, 0.0);
            for (k, yk) in y.iter().enumerate() {
                let (s, c) = libm::sincos(-2.0 * PI * alias_hz * k as f64 / FS);
                sr += yk.0 * c - yk.1 * s;
                si += yk.0 * s + yk.1 * c;
            }
            let image = (sr * sr + si * si).sqrt() / LEN as f64;
            let unfiltered = friis_volts(carrier_uhz as f64 * 1e-6);
            assert!(image <= 1e-4 * unfiltered, "image {image:.2e} V");
        }
    }
    // A case whose image is not quantized away: 8.3 MHz above the antenna-port centre
    // (B_out 7/8 MHz, 60 dB) folds to −7.7 MHz at about −130 dB, which is above one
    // LSB, so it renders. It must stay below the declared 60 dB stopband.
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![emitter(CENTER_UHZ + mhz(8.3), StaticWaveform::Cw)],
    };
    let builder = StaticSceneBuilder::new(
        &scene,
        &antenna_port(),
        &clock,
        BlockSpec::new(LEN).unwrap(),
    )
    .unwrap();
    assert_eq!(builder.build_block(0).unwrap().paths.len(), 1);
    let y = render_volts(&builder, 0);
    let (mut sr, mut si) = (0.0, 0.0);
    for (k, yk) in y.iter().enumerate() {
        let (s, c) = libm::sincos(-2.0 * PI * -7.7e6 * k as f64 / FS);
        sr += yk.0 * c - yk.1 * s;
        si += yk.0 * s + yk.1 * c;
    }
    let image = (sr * sr + si * si).sqrt() / LEN as f64;
    let unfiltered = friis_volts((CENTER_UHZ + mhz(8.3)) as f64 * 1e-6);
    eprintln!("8.3 MHz image {:.3e} of unfiltered", image / unfiltered);
    assert!(
        image > 0.0 && image <= 1e-3 * unfiltered,
        "image {image:.2e} V"
    );
}

#[test]
fn lossless_culling_never_changes_output_bits() {
    let clock = clock();
    let mut weak = emitter(CENTER_UHZ + mhz(1.0), StaticWaveform::Cw);
    weak.power_w = 1e-16; // ≈ 0.03 LSB at the receiver.
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![emitter(CENTER_UHZ + mhz(2.0), StaticWaveform::Cw), weak],
    };
    let spec = BlockSpec::new(LEN).unwrap();
    let mut keep_all = virtual_rx();
    keep_all.cull_threshold_lsb = 0.0;
    let culled = StaticSceneBuilder::new(&scene, &virtual_rx(), &clock, spec).unwrap();
    let full = StaticSceneBuilder::new(&scene, &keep_all, &clock, spec).unwrap();
    assert_eq!(culled.build_block(0).unwrap().paths.len(), 1);
    assert_eq!(full.build_block(0).unwrap().paths.len(), 2);
    assert_eq!(
        render_acc(&culled, 0, CpuBackend::new(1)),
        render_acc(&full, 0, CpuBackend::new(1))
    );
}

fn busy_scene() -> StaticScene {
    let block_end_ps = 1_024_000_000;
    let tau_ps = (link_tau() * 1e12) as i64;
    let train = |emission_ps: i64| StaticWaveform::Pulses {
        emission: TruthTime::from_ps(emission_ps),
        pulses: (0..4)
            .map(|i| PulseSpec {
                start_ps: i * 1_450_000,
                rise_ps: 50_000 + i * 25_000,
                fall_start_ps: i * 1_450_000 + 450_000,
                fall_ps: 100_000,
                amplitude: 1.0 - 0.1 * i as f32,
            })
            .collect(),
    };
    StaticScene {
        receiver: receiver(),
        emitters: vec![
            emitter(CENTER_UHZ + mhz(1.25), StaticWaveform::Cw),
            emitter(CENTER_UHZ - mhz(3.333_333), StaticWaveform::Cw),
            // Straddles the block 0 / block 1 boundary.
            emitter(
                CENTER_UHZ + mhz(1.3),
                train(block_end_ps - tau_ps - 2_000_000),
            ),
            emitter(mhz(1030.0), train(500_000_000)),
        ],
    }
}

#[test]
fn output_is_identical_across_threads_and_chunking() {
    let clock = clock();
    let scene = busy_scene();
    let spec = BlockSpec::new(LEN).unwrap();
    let mut bytes = Vec::new();
    for target in [virtual_rx(), antenna_port()] {
        let builder = StaticSceneBuilder::new(&scene, &target, &clock, spec).unwrap();
        for b in [0, 1] {
            // Both blocks must contain pulse trains, not only the CW paths.
            assert!(!builder.build_block(b).unwrap().trains.is_empty());
            let reference = render_acc(&builder, b, CpuBackend::new(1));
            assert!(reference.iter().any(|&v| v != 0));
            for backend in [
                CpuBackend::new(2),
                CpuBackend::new(8),
                CpuBackend::new(3).with_chunk_len(777),
                CpuBackend::new(1).with_chunk_len(1),
            ] {
                assert_eq!(render_acc(&builder, b, backend), reference, "{backend:?}");
            }
            for v in reference {
                bytes.extend_from_slice(&v.to_le_bytes());
            }
        }
    }
    iqsim_testkit::assert_golden("render/static_scene", &bytes);
}

#[test]
fn unsafe_headroom_is_rejected_and_overflow_is_an_error() {
    let clock = clock();
    let spec = BlockSpec::new(LEN).unwrap();
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![emitter(CENTER_UHZ, StaticWaveform::Cw)],
    };
    // A 1e-20 V full scale puts the 1.6e-4 V tone at ~2^94 LSB.
    let mut tiny = virtual_rx();
    tiny.full_scale_volts = 1e-20;
    let builder = StaticSceneBuilder::new(&scene, &tiny, &clock, spec).unwrap();
    assert!(matches!(
        builder.build_block(0),
        Err(SceneError::Amplitude(
            iqsim_core::AmplitudeError::Headroom { .. }
        ))
    ));

    // A hand-made work list that bypasses the proof must fail loudly, never wrap.
    let ok = StaticSceneBuilder::new(&scene, &virtual_rx(), &clock, spec).unwrap();
    let mut work = ok.build_block(0).unwrap();
    let lsb = f64::from(work.target.lsb_per_volt);
    work.paths[0].amplitude.c[0] = (2f64.sqrt() * 3.3e18 / lsb) as f32; // |v| ≈ 2^62.3 LSB
    work.paths.push(work.paths[0]);
    work.paths.push(work.paths[0]);
    let mut acc = vec![0_i64; 2 * LEN as usize];
    let e = CpuBackend::new(2).render(&work, &mut acc).unwrap_err();
    assert!(matches!(e, RenderError::Overflow { .. }), "{e:?}");
    work.paths.truncate(1);
    work.paths[0].amplitude.c[0] = (1e19 * 2.0 / lsb) as f32;
    let mut acc = vec![0_i64; 2 * LEN as usize];
    let e = CpuBackend::new(1).render(&work, &mut acc).unwrap_err();
    assert!(matches!(e, RenderError::Contribution { .. }), "{e:?}");
}

#[test]
fn antenna_port_is_isolated_from_virtual_receiver_settings() {
    let clock = clock();
    let scene = busy_scene();
    let spec = BlockSpec::new(LEN).unwrap();
    let port = StaticSceneBuilder::new(&scene, &antenna_port(), &clock, spec).unwrap();
    let port_hash = port.build_block(1).unwrap().content_hash();
    let port_out = render_acc(&port, 1, CpuBackend::new(2));
    let mut outputs = Vec::new();
    for (pass, fs_mhz, gd, vfs) in [(6.0, 7.5, 1_000_000, 1.0), (5.0, 8.0, 0, 0.1)] {
        let mut rx = virtual_rx();
        rx.response.passband_edge_uhz = mhz(pass);
        rx.response.stopband_edge_uhz = mhz(fs_mhz);
        rx.response.group_delay_ps = gd;
        rx.full_scale_volts = vfs;
        let b = StaticSceneBuilder::new(&scene, &rx, &clock, spec).unwrap();
        outputs.push(render_acc(&b, 1, CpuBackend::new(2)));
        // Rebuilding the antenna-port target after any virtual-receiver change gives
        // the same work list and the same samples.
        let again = StaticSceneBuilder::new(&scene, &antenna_port(), &clock, spec).unwrap();
        assert_eq!(again.build_block(1).unwrap().content_hash(), port_hash);
        assert_eq!(render_acc(&again, 1, CpuBackend::new(2)), port_out);
    }
    // The virtual-receiver output itself does change with its settings.
    assert_ne!(outputs[0], outputs[1]);
    // The same response definition on the two taps has different identities.
    let same = antenna_port().response;
    let a = TapResponse::new(TapKind::AntennaPort, same, rate(), CENTER_UHZ).unwrap();
    let v = TapResponse::new(TapKind::VirtualReceiver, same, rate(), CENTER_UHZ).unwrap();
    assert_ne!(a.id(), v.id());
}

#[test]
fn a_rise_that_outlasts_the_fall_renders_the_same_in_any_chunking() {
    // Regression (prompt 5 review): rise 10 µs, fall start 1 µs, fall 0. The pulse is
    // non-zero until the rise table ends, after the fall table; small chunks used to
    // drop that tail.
    let clock = clock();
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![emitter(
            CENTER_UHZ,
            StaticWaveform::Pulses {
                emission: TruthTime::from_ps(300_000_000),
                pulses: vec![PulseSpec {
                    start_ps: 0,
                    rise_ps: 10_000_000,
                    fall_start_ps: 1_000_000,
                    fall_ps: 0,
                    amplitude: 1.0,
                }],
            },
        )],
    };
    let builder =
        StaticSceneBuilder::new(&scene, &virtual_rx(), &clock, BlockSpec::new(LEN).unwrap())
            .unwrap();
    let whole = render_acc(&builder, 0, CpuBackend::new(1));
    for chunk in [1, 7, 64, 777] {
        assert_eq!(
            render_acc(&builder, 0, CpuBackend::new(2).with_chunk_len(chunk)),
            whole,
            "chunk {chunk}"
        );
    }
    // The tail between the fall table's end and the rise table's end is non-zero.
    let tau = link_tau();
    let k = ((300e-6 + tau + 1e-6 + 6e-6) * FS) as usize;
    assert!(whole[2 * k] != 0 || whole[2 * k + 1] != 0);
}

#[test]
fn headroom_covers_the_scaled_amplitude_intermediate() {
    // Regression (prompt 5 review): a strong path with a tiny pulse amplitude passed
    // the old proof (bound ≈ 10⁹ LSB) while a·lsb_per_volt overflowed f32 in the kernel.
    let clock = clock();
    let mut e = emitter(
        CENTER_UHZ,
        StaticWaveform::Pulses {
            emission: TruthTime::from_ps(300_000_000),
            pulses: vec![PulseSpec {
                start_ps: 0,
                rise_ps: 100_000,
                fall_start_ps: 450_000,
                fall_ps: 100_000,
                amplitude: 1e-30,
            }],
        },
    );
    e.power_w = 1e12; // ≈ 5 V at the receiver.
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![e],
    };
    let mut target = virtual_rx();
    target.full_scale_volts = 1.1e-26; // lsb_per_volt ≈ 10³⁸
    target.cull_threshold_lsb = 0.0;
    let builder =
        StaticSceneBuilder::new(&scene, &target, &clock, BlockSpec::new(LEN).unwrap()).unwrap();
    assert!(matches!(
        builder.build_block(0),
        Err(SceneError::Amplitude(
            iqsim_core::AmplitudeError::Headroom { .. }
        ))
    ));
}

#[test]
fn links_beyond_the_light_time_range_fail_the_block() {
    // Prompt 5 rejected f_c·R/c above 2⁴⁰ cycles because f64 could not resolve its
    // fraction; the geometry engine computes that fraction exactly (`carrier_phase`).
    // A 10¹⁵ m link is still rejected, by the light-time range (|τ| ≤ 10⁶ s), when the
    // block is built.
    let clock = clock();
    let mut far = emitter(CENTER_UHZ, StaticWaveform::Cw);
    far.position = Ecef::new(6_378_137.0, 1e15, 0.0);
    let scene = StaticScene {
        receiver: receiver(),
        emitters: vec![far],
    };
    let builder =
        StaticSceneBuilder::new(&scene, &virtual_rx(), &clock, BlockSpec::new(LEN).unwrap())
            .unwrap();
    assert!(matches!(
        builder.build_block(0),
        Err(SceneError::Fit {
            index: 0,
            block: 0,
            ..
        })
    ));
}
