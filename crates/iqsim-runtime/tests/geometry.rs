//! The light-time geometry engine (prompt 9): a moving CW emitter's frequency and
//! phase trajectory against an analytic solution, a LEO-like mover that forces fit
//! refinement and still meets every bound at dense points, splits at trajectory knots
//! mapped from emission to receive time, and output that does not depend on threads or
//! chunking.

#![forbid(unsafe_code)]

mod common;

use std::f64::consts::PI;
use std::sync::Arc;

use common::{C, EciMotion, reference_delay};
use iqsim_core::amplitude::LOSSLESS_CULL_LSB;
use iqsim_core::{
    AcqInstant, BlockSpec, ChannelTimebase, ClockMap, Delay, IdealClockMap, Phase, ResponseSpec,
    SampleIndex, SampleRate, TapKind, TruthDuration, TruthTime,
};
use iqsim_geo::light_time;
use iqsim_geo::trajectory::{AttitudeMode, FlatGround, Heights, Limits, Primitive, Start};
use iqsim_geo::{AltitudeModel, Ecef, Flight, Geodetic, Quat, Static, Trajectory};
use iqsim_render::{Backend, CpuBackend};
use iqsim_runtime::geometry::{
    AMPLITUDE_BUDGET_REL, DELAY_BUDGET_SAMPLES, PHASE_BUDGET_RAD, SplitReason,
};
use iqsim_runtime::scene::{
    PulseSpec, Scene, SceneBuilder, SceneEmitter, StaticTarget, StaticWaveform,
};
use iqsim_worklist::{BlockWorkList, PathSegment};

const R_E: f64 = 6_378_137.0;
const RX: [f64; 3] = [R_E, 0.0, 0.0];

fn mhz(x: f64) -> i64 {
    (x * 1e12) as i64
}

fn target(center_uhz: i64, rate: SampleRate, half_band_hz: f64) -> StaticTarget {
    StaticTarget {
        kind: TapKind::AntennaPort,
        center_uhz,
        rate,
        response: ResponseSpec {
            passband_edge_uhz: (half_band_hz * 1e6) as i64,
            stopband_edge_uhz: (half_band_hz * 1.25e6) as i64,
            stopband_atten_mdb: 80_000,
            group_delay_ps: 0,
        },
        full_scale_volts: 1e-3,
        cull_threshold_lsb: LOSSLESS_CULL_LSB,
    }
}

fn receiver() -> Arc<dyn Trajectory> {
    Arc::new(Static::from_ecef(
        Ecef::new(RX[0], RX[1], RX[2]),
        Quat::IDENTITY,
    ))
}

fn emitter(t: Arc<dyn Trajectory>, carrier_uhz: i64, waveform: StaticWaveform) -> SceneEmitter {
    SceneEmitter {
        trajectory: t,
        carrier_uhz,
        power_w: 10.0,
        waveform,
    }
}

/// The segments of emitter `link` in a work list, in sample order.
fn segments(w: &BlockWorkList, link: u32) -> Vec<PathSegment> {
    let mut v: Vec<PathSegment> = w
        .paths
        .iter()
        .filter(|p| p.coord_source == link)
        .copied()
        .collect();
    v.sort_by_key(|p| p.first);
    v
}

fn wrap_rad(cycles: f64) -> f64 {
    (cycles - libm::round(cycles)).abs() * 2.0 * PI
}

/// Moving CW emitter: the rendered phase follows f_off·t − f_e·τ(t) for an emitter in
/// uniform ECI motion (closed-form light time), and the frequency follows the analytic
/// Doppler, through closest approach.
#[test]
fn moving_cw_emitter_matches_the_analytic_phase_and_frequency() {
    let fs = 2_000_000.0;
    let rate = SampleRate::from_hz(2_000_000).unwrap();
    let carrier_uhz = mhz(2200.0);
    let f_e = 2.2e9;
    let center = carrier_uhz - mhz(0.2);
    let f_off = 0.2e6;
    // 7.5 km/s, passing 50 km over the receiver at t ≈ 4 s, with 30 m/s up.
    let motion = EciMotion::Uniform {
        p0: [R_E + 50e3, -30e3, 2e3],
        u: [30.0, 7_500.0, 0.0],
    };
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 20_000; // 10 ms blocks
    let scene = Scene {
        receiver: receiver(),
        emitters: vec![emitter(Arc::new(motion), carrier_uhz, StaticWaveform::Cw)],
    };
    let tgt = target(center, rate, 800e3);
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let mut worst_phase = 0.0_f64;
    let mut worst_freq = 0.0_f64;
    let mut doppler_span = (f64::MAX, f64::MIN);
    // Blocks before, through and after closest approach (t ≈ 4 s).
    for block in [0_u64, 200, 399, 400, 401, 600] {
        let work = b.build_block(block).unwrap();
        let mut acc = vec![0_i64; 2 * len as usize];
        CpuBackend::new(4).render(&work, &mut acc).unwrap();
        let mut phases = Vec::with_capacity(len as usize);
        for k in 0..len as usize {
            let n = block * len + k as u64;
            let t = n as f64 / fs;
            let (i, q) = (acc[2 * k] as f64, acc[2 * k + 1] as f64);
            let measured = libm::atan2(q, i) / (2.0 * PI);
            let tau = reference_delay(RX, &motion, t);
            let want = f_off * t - f_e * tau;
            worst_phase = worst_phase.max(wrap_rad(measured - want));
            phases.push(measured);
        }
        // Frequency over 2000-sample windows (1 ms) from the unwrapped phase.
        for w0 in (0..len as usize - 2000).step_by(4000) {
            let mut acc_phase = 0.0;
            for k in w0..w0 + 2000 {
                let mut d = phases[k + 1] - phases[k];
                d -= libm::round(d);
                acc_phase += d;
            }
            let f_meas = acc_phase * fs / 2000.0;
            let (t0, t1) = (
                (block * len + w0 as u64) as f64 / fs,
                (block * len + w0 as u64 + 2000) as f64 / fs,
            );
            let f_ref = f_off
                - f_e * (reference_delay(RX, &motion, t1) - reference_delay(RX, &motion, t0))
                    / (t1 - t0);
            worst_freq = worst_freq.max((f_meas - f_ref).abs());
            let dop = f_ref - f_off;
            doppler_span = (doppler_span.0.min(dop), doppler_span.1.max(dop));
        }
    }
    eprintln!(
        "moving CW: worst phase {worst_phase:.3e} rad, worst frequency {worst_freq:.3e} Hz, \
         Doppler {:.1} .. {:.1} Hz",
        doppler_span.0, doppler_span.1
    );
    assert!(doppler_span.0 < -10e3 && doppler_span.1 > 20e3);
    assert!(worst_phase <= PHASE_BUDGET_RAD, "phase {worst_phase}");
    // 0.01 rad over a 1 ms window is 1.6 Hz; the estimate is far tighter.
    assert!(worst_freq < 0.05, "frequency {worst_freq}");
}

/// A LEO-like circular orbit passing overhead, at 10 GHz, in 4 s blocks: a single
/// cubic per block misses the budget, the fit refines, and every segment then meets
/// every bound at dense points against an independent light-time reference.
#[test]
fn a_leo_pass_forces_refinement_and_meets_every_bound() {
    let fs = 1_000_000.0;
    let rate = SampleRate::from_hz(1_000_000).unwrap();
    let carrier_uhz = mhz(10_000.0);
    let f_e = 1e10;
    let center = carrier_uhz - mhz(0.1);
    let f_off = 0.1e6;
    let a = R_E + 550e3;
    // Overhead (angle 0 in ECI, where the receiver also is at t ≈ 0) at t ≈ 6 s.
    let n = libm::sqrt(common::MU / (a * a * a));
    let motion = EciMotion::Circular {
        a,
        phase0: -n * 6.0,
    };
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 4_000_000; // 4 s blocks
    let emission = TruthTime::from_ps(5_000_000_000_000);
    let pulses = StaticWaveform::Pulses {
        emission,
        pulses: vec![PulseSpec {
            start_ps: 0,
            rise_ps: 1_000_000,
            fall_start_ps: 20_000_000,
            fall_ps: 1_000_000,
            amplitude: 1.0,
        }],
    };
    let scene = Scene {
        receiver: receiver(),
        emitters: vec![
            emitter(Arc::new(motion), carrier_uhz, StaticWaveform::Cw),
            emitter(Arc::new(motion), carrier_uhz, pulses),
        ],
    };
    let tgt = target(center, rate, 400e3);
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let resp = *b.response();
    let mut refined = 0usize;
    let mut worst = [0.0_f64; 4];
    for block in [0_u64, 1, 2] {
        let (work, log) = b.build_block_logged(block).unwrap();
        for p in &log.paths {
            refined += p
                .fit
                .splits
                .iter()
                .filter(|s| s.reason == SplitReason::Refine)
                .count();
        }
        let cw = segments(&work, 0);
        let pl = segments(&work, 1);
        assert!(cw.len() > 1, "block {block}: {} segment(s)", cw.len());
        // Segments tile the block.
        assert_eq!(cw[0].first, 0);
        for w in cw.windows(2) {
            assert_eq!(w[0].first + w[0].len, w[1].first);
        }
        assert_eq!(
            u64::from(cw.last().unwrap().first + cw.last().unwrap().len),
            len
        );
        for (link, segs) in [(0, &cw), (1, &pl)] {
            for s in segs.iter() {
                let sl = u64::from(s.len);
                let delay = s.delay.to_poly::<Delay>(sl).unwrap();
                let phase = s.phase.to_poly::<Phase>(sl).unwrap();
                let coord = s.coord.to_poly::<Delay>(sl).unwrap();
                let step = (sl / 97).max(1);
                for j in (0..sl).step_by(step as usize).chain([sl - 1]) {
                    let k = block * len + u64::from(s.first) + j;
                    let t = k as f64 / fs;
                    let tau = reference_delay(RX, &motion, t);
                    let d = (delay.eval(j).unwrap().to_samples_f64() - tau * fs).abs();
                    let ph = phase.eval(j).unwrap().to_cycles_f64();
                    let p = wrap_rad(ph - (f_off * t - f_e * tau));
                    worst[0] = worst[0].max(d);
                    worst[1] = worst[1].max(p);
                    if link == 0 {
                        // Friis · |H(f_inst)|.
                        let h = 1e-4;
                        let rate_tau = (reference_delay(RX, &motion, t + h)
                            - reference_delay(RX, &motion, t - h))
                            / (2.0 * h);
                        let want = libm::sqrt(10.0 * 50.0) * (C / f_e) / (4.0 * PI * C * tau)
                            * resp.magnitude(f_off - f_e * rate_tau);
                        let got = f64::from(s.amplitude.eval(j));
                        worst[3] = worst[3].max((got / want - 1.0).abs());
                    } else {
                        // Retarded coordinate F_s·(t − τ − t_ref).
                        let w = coord.eval(j).unwrap().to_samples_f64();
                        let want = fs * (t - tau - 5.0);
                        worst[2] = worst[2].max((w - want).abs());
                    }
                }
            }
        }
    }
    eprintln!(
        "LEO: {refined} refinement splits; worst delay {:.3e} samples, phase {:.3e} rad, \
         coordinate {:.3e} samples, amplitude {:.3e}",
        worst[0], worst[1], worst[2], worst[3]
    );
    assert!(refined > 0);
    assert!(worst[0] <= DELAY_BUDGET_SAMPLES);
    assert!(worst[1] <= PHASE_BUDGET_RAD);
    assert!(worst[2] <= DELAY_BUDGET_SAMPLES);
    assert!(worst[3] <= AMPLITUDE_BUDGET_REL);
}

/// An aircraft whose maneuvers start at t = 0 (light time looks before its start) and
/// whose turn-in and roll-out knots fall inside blocks: the fit splits there, each
/// split is the first sample whose emission time is at or after the knot, and the
/// rendered output does not depend on threads or chunking.
#[test]
fn transmitter_knots_split_at_their_receive_time_and_output_is_thread_independent() {
    let rate = SampleRate::from_hz(4_000_000).unwrap();
    let carrier_uhz = mhz(1090.0);
    let model = AltitudeModel::default();
    let heights = Heights {
        model: &model,
        ground: &FlatGround::ELLIPSOID,
    };
    let start = Start {
        time: TruthTime::EPOCH,
        position: Geodetic::from_degrees(0.3, 0.2, 9000.0),
        speed_m_s: 250.0,
        track_rad: Some(0.0),
    };
    let flight = Flight::maneuvers(
        start,
        Limits::default(),
        AttitudeMode::Derived {
            angle_of_attack_rad: 0.0,
            wind_ned_m_s: iqsim_geo::Vec3::ZERO,
        },
        &[
            Primitive::Straight {
                duration: TruthDuration::from_ms(3).unwrap(),
            },
            Primitive::Turn {
                rate: iqsim_geo::trajectory::TurnRate::Rate(1f64.to_radians()),
                direction: Some(iqsim_geo::trajectory::TurnDirection::Right),
                end: iqsim_geo::trajectory::TurnEnd::By(0.5f64.to_radians()),
            },
        ],
        heights,
    )
    .unwrap();
    let mut flight = flight;
    flight
        .extend_to(TruthTime::from_ps(5_000_000_000_000), heights)
        .unwrap();
    let tx: Arc<dyn Trajectory> = Arc::new(flight);
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 8_000; // 2 ms blocks
    let scene = Scene {
        receiver: Arc::new(Static::new(
            Geodetic::from_degrees(0.0, 0.0, 10.0),
            Quat::IDENTITY,
        )),
        emitters: vec![emitter(Arc::clone(&tx), carrier_uhz, StaticWaveform::Cw)],
    };
    let tgt = target(carrier_uhz - mhz(0.5), rate, 1.6e6);
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let geometry = b.path_geometry(0).unwrap().clone();
    let knots = tx.knots(TruthTime::EPOCH, TruthTime::from_ps(4_000_000_000_000));
    assert!(knots.len() >= 3, "{knots:?}");
    let emit = |n: u64| {
        let t = clock.instant(SampleIndex::new(n)).unwrap();
        geometry.at(t).unwrap().tx.instant
    };
    // Block 0 hears emissions from before the flight starts (the look-back).
    assert!(emit(0) < AcqInstant::from_truth_time(TruthTime::EPOCH));
    let mut blocks = vec![0_u64];
    for &kn in &knots {
        let lt = light_time::solve_emit(
            &*geometry.rx,
            &*geometry.tx,
            AcqInstant::from_truth_time(kn),
        )
        .unwrap();
        let n = clock.index_at_or_after(lt.rx.instant).unwrap().get();
        blocks.push(n / len);
    }
    blocks.dedup();
    let mut found = 0;
    let mut hashes = Vec::new();
    for &block in &blocks {
        let (work, log) = b.build_block_logged(block).unwrap();
        for s in log.paths[0]
            .fit
            .splits
            .iter()
            .filter(|s| s.reason == SplitReason::TransmitterKnot)
        {
            found += 1;
            let k = block * len + s.at;
            let hit = knots.iter().any(|&kn| {
                let kn = AcqInstant::from_truth_time(kn);
                emit(k - 1) < kn && kn <= emit(k)
            });
            assert!(hit, "block {block}: split at {} is not at a knot", s.at);
        }
        let mut outs = Vec::new();
        for backend in [
            CpuBackend::new(1),
            CpuBackend::new(3),
            CpuBackend::new(8).with_chunk_len(777),
            CpuBackend::new(2).with_chunk_len(1),
        ] {
            let mut acc = vec![0_i64; 2 * len as usize];
            backend.render(&work, &mut acc).unwrap();
            outs.push(acc);
        }
        assert!(outs.windows(2).all(|w| w[0] == w[1]), "block {block}");
        assert!(outs[0].iter().any(|&v| v != 0));
        hashes.push(work.content_hash());
    }
    // Every knot whose receive time is inside a block (not on a block boundary).
    assert!(
        found + 1 >= knots.len(),
        "{found} transmitter-knot splits for {knots:?}"
    );
    // Rebuilding gives the same work lists.
    let b2 = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    assert_eq!(b2.build_block(0).unwrap().content_hash(), hashes[0]);
}

// ------------------------------------------------------------------ review regressions

/// Review P1: a Doppler sweep through a sharp response edge could vanish between
/// check points. f_inst may now move only σ/8 between checks, so the fit refines and
/// every sample's amplitude matches Friis·|H(f_inst)|.
#[test]
fn a_doppler_sweep_through_a_sharp_edge_is_resolved() {
    let fs = 1000.0;
    let rate = SampleRate::from_hz(1000).unwrap();
    let carrier_uhz = mhz(31_000.0);
    let f_e = 31e9;
    // 100 m/s relative to the (Earth-rotating) receiver, passing 1 km from it at
    // t ≈ 0.6 s: Doppler sweeps ≈ 1 kHz/s through 0.
    let motion = EciMotion::Uniform {
        p0: [R_E + 1000.0, -60.0, 0.0],
        u: [0.0, 100.0 + common::OMEGA_E * R_E, 0.0],
    };
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 1000;
    let scene = Scene {
        receiver: receiver(),
        emitters: vec![emitter(Arc::new(motion), carrier_uhz, StaticWaveform::Cw)],
    };
    let tgt = target(carrier_uhz, rate, 100.0);
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let resp = *b.response();
    let (work, log) = b.build_block_logged(0).unwrap();
    let segs = log.paths[0].fit.segments.len();
    let v = |tau: f64| libm::sqrt(500.0) * (C / f_e) / (4.0 * PI * C * tau);
    let mut worst = 0.0_f64;
    let mut seen_in_band = false;
    let culled_or = segments(&work, 0);
    for k in 0..len {
        let t = k as f64 / fs;
        let tau = reference_delay(RX, &motion, t);
        let h = 1e-4;
        let rate_tau =
            (reference_delay(RX, &motion, t + h) - reference_delay(RX, &motion, t - h)) / (2.0 * h);
        let want = v(tau) * resp.magnitude(-f_e * rate_tau);
        seen_in_band |= resp.magnitude(-f_e * rate_tau) > 0.5;
        let got = culled_or
            .iter()
            .find(|s| u64::from(s.first) <= k && k < u64::from(s.first + s.len))
            .map_or(0.0, |s| f64::from(s.amplitude.eval(k - u64::from(s.first))));
        // 0.01 dB above the −40 dB floor, 10⁻⁵·V_ref below it (contracts §2).
        let allowed = (AMPLITUDE_BUDGET_REL * want).max(1e-5 * v(tau)) + 1e-12;
        worst = worst.max((got - want).abs() / allowed);
    }
    eprintln!("sharp edge: {segs} segments, worst amplitude error {worst:.3} of budget");
    assert!(seen_in_band);
    assert!(segs > 16, "{segs}");
    assert!(worst <= 1.0, "{worst}");
}

/// Review P1: the CW group delay is applied exactly (the geometry at t − τ_g). A
/// 31 GHz pass at 7.5 km/s and 550 km with τ_g = 1 ms is 0.03 rad away from the
/// first-order fold at closest approach; the rendered phase follows the exact one.
#[test]
fn cw_group_delay_follows_the_delayed_trajectory() {
    let fs = 2_000_000.0;
    let rate = SampleRate::from_hz(2_000_000).unwrap();
    let carrier_uhz = mhz(31_000.0);
    let f_e = 31e9;
    let motion = EciMotion::Uniform {
        p0: [R_E + 550e3, -7_500.0, 0.0],
        u: [0.0, 7_500.0, 0.0],
    };
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 20_000;
    let mut tgt = target(carrier_uhz, rate, 800e3);
    tgt.response.group_delay_ps = 1_000_000_000;
    let scene = Scene {
        receiver: receiver(),
        emitters: vec![emitter(Arc::new(motion), carrier_uhz, StaticWaveform::Cw)],
    };
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let tg = 1e-3;
    let (mut worst, mut first_order) = (0.0_f64, 0.0_f64);
    for block in [49_u64, 50] {
        let work = b.build_block(block).unwrap();
        let mut acc = vec![0_i64; 2 * len as usize];
        CpuBackend::new(4).render(&work, &mut acc).unwrap();
        for k in (0..len as usize).step_by(37) {
            let t = (block * len + k as u64) as f64 / fs;
            let measured = libm::atan2(acc[2 * k + 1] as f64, acc[2 * k] as f64) / (2.0 * PI);
            let exact = -f_e * reference_delay(RX, &motion, t - tg);
            let h = 1e-4;
            let f_inst = -f_e
                * (reference_delay(RX, &motion, t + h) - reference_delay(RX, &motion, t - h))
                / (2.0 * h);
            let approx = -f_e * reference_delay(RX, &motion, t) - f_inst * tg;
            worst = worst.max(wrap_rad(measured - exact));
            first_order = first_order.max(wrap_rad(exact - approx));
        }
    }
    eprintln!(
        "group delay: worst {worst:.3e} rad; first-order fold would miss {first_order:.3e} rad"
    );
    assert!(first_order > PHASE_BUDGET_RAD);
    assert!(worst <= PHASE_BUDGET_RAD, "{worst}");
}

/// Review P1: pulses are filtered at the Doppler-shifted offset. A 31 GHz emitter
/// closing at 7.5 km/s is shifted by ≈ +775.5 kHz; with its nominal offset at
/// −775.5 kHz the received pulse sits at the channel centre, and its plateau has the
/// full Friis level (at the nominal offset it would have been out of band).
#[test]
fn pulses_are_filtered_at_their_doppler_shifted_offset() {
    let fs = 1_000_000.0;
    let rate = SampleRate::from_hz(1_000_000).unwrap();
    let f_e = 31e9;
    let doppler = f_e * 7_500.0 / C;
    let center = mhz(31_000.0) + (doppler * 1e6) as i64;
    // Closing head-on along x.
    let motion = EciMotion::Uniform {
        p0: [R_E + 300e3, 0.0, 0.0],
        u: [-7_500.0, 0.0, 0.0],
    };
    let clock = IdealClockMap::new(ChannelTimebase::new(TruthTime::EPOCH, rate));
    let len = 4_000;
    let pulses = StaticWaveform::Pulses {
        emission: TruthTime::EPOCH,
        pulses: vec![PulseSpec {
            start_ps: 0,
            rise_ps: 10_000_000,
            fall_start_ps: 3_000_000_000,
            fall_ps: 10_000_000,
            amplitude: 1.0,
        }],
    };
    let scene = Scene {
        receiver: receiver(),
        emitters: vec![emitter(Arc::new(motion), mhz(31_000.0), pulses)],
    };
    let tgt = target(center, rate, 400e3);
    let b = SceneBuilder::new(&scene, &tgt, &clock, BlockSpec::new(len).unwrap()).unwrap();
    let work = b.build_block(0).unwrap();
    let mut acc = vec![0_i64; 2 * len as usize];
    CpuBackend::new(2).render(&work, &mut acc).unwrap();
    let q = b.scale().lsb_volts();
    // Mid-pulse: emitted 1.5 ms after 0, received ≈ 1 ms later.
    let k = 2_500;
    let t = k as f64 / fs;
    let tau = reference_delay(RX, &motion, t);
    let want = libm::sqrt(500.0) * (C / f_e) / (4.0 * PI * C * tau);
    let got = libm::hypot(acc[2 * k] as f64 * q, acc[2 * k + 1] as f64 * q);
    let db = 20.0 * libm::log10(got / want);
    eprintln!("Doppler-shifted pulse: level {db:.2e} dB");
    assert!(db.abs() < 0.01, "{db} dB");
}

/// Review P2: a moving receiver whose plan ends inside the final block, with a clock
/// offset, renders: trajectories are extended past the last fitted instant.
#[test]
fn trajectories_cover_the_whole_final_block() {
    let toml = r#"
schema_version = 1
epoch = 2026-01-01T00:00:00Z
seed = 1
run_mode = "offline"
[metadata]
name = "tail"
[render]
duration = "1.5 ms"
block = "1 ms"
[[platforms]]
name = "beacon"
position = { lat = 51.0, lon = -1.0, alt = "10 m AGL" }
[[platforms]]
name = "car"
position = { lat = 51.01, lon = -1.0, alt = "2 m AGL" }
[platforms.trajectory]
kind = "ground"
heading = "0 deg"
speed = "20 m/s"
[[platforms.trajectory.maneuvers]]
kind = "straight"
duration = "1 ms"
[[emitters]]
name = "cw"
platform = "beacon"
frequency = "1090 MHz"
power = "30 dBm"
behavior = "always_on"
waveform = "cw"
[[receivers]]
name = "rx"
platform = "car"
[[receivers.channels]]
name = "ch0"
center_frequency = "1090 MHz"
sample_rate = "4 MS/s"
bandwidth = "3 MHz"
clock = { model = "fixed_offset", offset = "100 us" }
sinks = ["rec"]
[[sinks]]
name = "rec"
kind = "sigmf"
"#;
    let dir = std::env::temp_dir().join(format!("iqsim-tail-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("tail.toml");
    std::fs::write(&path, toml).unwrap();
    let m = iqsim_runtime::run_offline(
        &path,
        &dir.join("out"),
        iqsim_runtime::RunOptions::default(),
    )
    .unwrap();
    assert_eq!(m.outputs[0].sample_count, 6000);
    let _ = std::fs::remove_dir_all(dir);
}
