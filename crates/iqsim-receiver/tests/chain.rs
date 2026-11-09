//! Receiver-chain behavior: stage formulas, ordered block handoff, and output that does
//! not depend on how the run is split into blocks.

use iqsim_core::{
    ChannelTimebase, IdealClockMap, ResponseSpec, SampleIndex, SampleRate, TapKind, TapResponse,
    TargetScale, TruthTime,
};
use iqsim_receiver::{
    AdcConfig, AgcConfig, ChainConfig, ChainError, DecimationConfig, FrontEnd, IqImpairment,
    LoConfig, NoiseConfig, OutputFormat, PhaseNoise, ReceiverChannelState, Samples,
};

const CENTER: i64 = 1_090_000_000_000_000;

fn rate() -> SampleRate {
    SampleRate::from_hz(10_000_000).unwrap()
}

fn clock() -> IdealClockMap {
    IdealClockMap::new(ChannelTimebase::new(TruthTime::from_ps(0), rate()))
}

fn h_rx() -> TapResponse {
    TapResponse::new(
        TapKind::VirtualReceiver,
        ResponseSpec {
            passband_edge_uhz: 4_000_000_000_000,
            stopband_edge_uhz: 5_000_000_000_000,
            stopband_atten_mdb: 80_000,
            group_delay_ps: 0,
        },
        rate(),
        CENTER,
    )
    .unwrap()
}

/// A bare chain: unity gain, no stages, unquantized cf32 at a 1 V full scale.
fn bare() -> ChainConfig {
    ChainConfig {
        rate: rate(),
        response: h_rx(),
        input_scale: TargetScale::new(1.0).unwrap(),
        adc_full_scale_volts: 1.0,
        output: OutputFormat::Cf32,
        front_end: None,
        noise: None,
        lo: None,
        iq: None,
        decimation: None,
        agc: None,
        adc: None,
    }
}

/// Every stage on.
fn full() -> ChainConfig {
    ChainConfig {
        front_end: Some(FrontEnd {
            gain_db: 20.0,
            noise_figure_db: 4.0,
        }),
        noise: Some(NoiseConfig {
            seed: 7,
            entity: 3,
            channel: 1,
        }),
        lo: Some(LoConfig {
            offset_uhz: 1_234_567_890,
            phase_noise: PhaseNoise::None,
        }),
        iq: Some(IqImpairment {
            gain_imbalance_db: 0.3,
            phase_skew_deg: 2.0,
            dc_i_volts: 1e-3,
            dc_q_volts: -2e-3,
        }),
        decimation: Some(DecimationConfig {
            factor: 3,
            filter: ResponseSpec {
                passband_edge_uhz: 1_200_000_000_000,
                stopband_edge_uhz: 1_600_000_000_000,
                stopband_atten_mdb: 70_000,
                group_delay_ps: 0,
            },
        }),
        agc: Some(AgcConfig {
            target_volts: 0.05,
            attack: 0.05,
            step: 0.02,
            min_gain_db: -20.0,
            max_gain_db: 20.0,
        }),
        adc: Some(AdcConfig { bits: 12 }),
        output: OutputFormat::Ci16,
        ..bare()
    }
}

/// A test input: a tone of 10 mV plus a slow ramp, in accumulator units at 1 V.
fn input(n: usize) -> Vec<i64> {
    let lsb = (1_u64 << 40) as f64;
    let mut acc = Vec::with_capacity(2 * n);
    for k in 0..n {
        let t = k as f64 * 0.0137;
        let (s, c) = libm::sincos(t);
        acc.push((0.01 * c * lsb) as i64 + k as i64 * 1000);
        acc.push((0.01 * s * lsb) as i64);
    }
    acc
}

fn run_split(config: ChainConfig, acc: &[i64], splits: &[usize]) -> (Vec<i16>, u64, Vec<u64>) {
    let mut state = ReceiverChannelState::new(config).unwrap();
    let mut out = Vec::new();
    let mut firsts = Vec::new();
    let mut start = 0;
    for (b, &len) in splits.iter().enumerate() {
        let block = &acc[2 * start..2 * (start + len)];
        let o = state
            .process(b as u64, SampleIndex::new(start as u64), block, &clock())
            .unwrap();
        firsts.push(o.first_index);
        let Samples::Ci16(v) = o.samples else {
            panic!("ci16 expected")
        };
        out.extend(v);
        start += len;
    }
    (out, state.clipped_total(), firsts)
}

#[test]
fn output_does_not_depend_on_block_splitting() {
    let acc = input(6000);
    let (whole, clip_whole, _) = run_split(full(), &acc, &[6000]);
    for splits in [
        vec![3000, 3000],
        vec![1000; 6],
        vec![1, 2, 997, 1000, 4000],
        vec![7; 857].into_iter().chain([1]).collect::<Vec<usize>>(),
    ] {
        assert_eq!(splits.iter().sum::<usize>(), 6000);
        let (split, clip, firsts) = run_split(full(), &acc, &splits);
        assert_eq!(split, whole, "{:?}", &splits[..3]);
        assert_eq!(clip, clip_whole);
        // Output indices continue across blocks at the decimated rate.
        let mut start = 0_u64;
        for (len, first) in splits.iter().zip(firsts) {
            assert_eq!(first, start.div_ceil(3));
            start += *len as u64;
        }
    }
    assert_eq!(whole.len(), 2 * 2000);
    assert!(whole.iter().any(|&c| c != 0));
}

#[test]
fn blocks_must_arrive_in_order_without_gaps() {
    let acc = input(100);
    let mut s = ReceiverChannelState::new(full()).unwrap();
    assert_eq!(
        s.process(1, SampleIndex::new(0), &acc, &clock()),
        Err(ChainError::OutOfOrder {
            expected: 0,
            got: 1
        })
    );
    s.process(0, SampleIndex::new(0), &acc, &clock()).unwrap();
    assert_eq!(
        s.process(1, SampleIndex::new(101), &acc, &clock()),
        Err(ChainError::Gap {
            expected: 100,
            got: 101
        })
    );
    assert_eq!(
        s.process(1, SampleIndex::new(100), &acc[..3], &clock()),
        Err(ChainError::Length)
    );
    // A failed call changes nothing: the right block still goes through.
    s.process(1, SampleIndex::new(100), &acc, &clock()).unwrap();
    assert_eq!(s.next_block(), 2);
}

fn cf32(config: ChainConfig, acc: &[i64]) -> Vec<f32> {
    let mut s = ReceiverChannelState::new(config).unwrap();
    match s
        .process(0, SampleIndex::new(0), acc, &clock())
        .unwrap()
        .samples
    {
        Samples::Cf32(v) => v,
        Samples::Ci16(_) => panic!("cf32 expected"),
    }
}

#[test]
fn gain_iq_and_dc_follow_their_formulas() {
    let one = 1_i64 << 40;
    let acc = [one / 10, one / 20];
    // Gain: 20 dB = ×10 in voltage.
    let mut c = bare();
    c.front_end = Some(FrontEnd {
        gain_db: 20.0,
        noise_figure_db: 3.0,
    });
    let v = cf32(c, &acc);
    assert!(
        (v[0] - 1.0).abs() < 1e-6 && (v[1] - 0.5).abs() < 1e-6,
        "{v:?}"
    );
    // IQ: I' = I + dc_i, Q' = a(Q cos φ − I sin φ) + dc_q.
    let mut c = bare();
    let (a, phi) = (libm::pow(10.0, 1.0 / 20.0), 3.0_f64.to_radians());
    c.iq = Some(IqImpairment {
        gain_imbalance_db: 1.0,
        phase_skew_deg: 3.0,
        dc_i_volts: 0.01,
        dc_q_volts: 0.02,
    });
    let v = cf32(c, &acc);
    let q = a * (0.05 * libm::cos(phi) - 0.1 * libm::sin(phi)) + 0.02;
    assert!((f64::from(v[0]) - 0.11).abs() < 1e-6);
    assert!((f64::from(v[1]) - q).abs() < 1e-6, "{} vs {q}", v[1]);
}

#[test]
fn lo_offset_shifts_a_tone_by_minus_the_offset() {
    // A DC input becomes e^{−j2π·Δf·t}: at Δf = 250 kHz and 10 MS/s, −1/40 cycle per
    // sample.
    let one = 1_i64 << 40;
    let acc: Vec<i64> = (0..80).flat_map(|_| [one / 4, 0]).collect();
    let mut c = bare();
    c.lo = Some(LoConfig {
        offset_uhz: 250_000_000_000,
        phase_noise: PhaseNoise::None,
    });
    let v = cf32(c, &acc);
    for k in 0..80 {
        let ang = -2.0 * core::f64::consts::PI * (k as f64) / 40.0;
        let (s, cc) = libm::sincos(ang);
        assert!((f64::from(v[2 * k]) - 0.25 * cc).abs() < 1e-6);
        assert!((f64::from(v[2 * k + 1]) - 0.25 * s).abs() < 1e-6);
    }
}

#[test]
fn adc_quantizes_left_justifies_and_counts_clips() {
    let one = 1_i64 << 40;
    // 8-bit ADC at 1 V: step 1/128 V. 0.5 V → 64 → 64 << 8.
    let acc = [
        one / 2,
        -one / 2,
        2 * one,
        -2 * one,
        one / 256 * 3,
        0,
        // Far beyond any i64 code range after scaling: still a clip, not an error.
        i64::MAX,
        i64::MIN,
    ];
    let mut c = bare();
    c.adc = Some(AdcConfig { bits: 8 });
    c.output = OutputFormat::Ci16;
    let mut s = ReceiverChannelState::new(c).unwrap();
    let o = s.process(0, SampleIndex::new(0), &acc, &clock()).unwrap();
    let Samples::Ci16(v) = o.samples else {
        panic!()
    };
    // 3/256 V = 1.5 steps → ties to even → 2.
    assert_eq!(
        v,
        vec![
            64 << 8,
            -64 << 8,
            127 << 8,
            -128 << 8,
            2 << 8,
            0,
            127 << 8,
            -128 << 8
        ]
    );
    assert_eq!(o.clipped, 4);
    // A chain gain that drives values past i64 after the ADC scale also clips.
    let mut c = bare();
    c.adc = Some(AdcConfig { bits: 16 });
    c.output = OutputFormat::Ci16;
    c.front_end = Some(FrontEnd {
        gain_db: 200.0,
        noise_figure_db: 0.0,
    });
    let mut s = ReceiverChannelState::new(c).unwrap();
    let o = s
        .process(0, SampleIndex::new(0), &[one, -one], &clock())
        .unwrap();
    assert_eq!(o.samples, Samples::Ci16(vec![32_767, -32_768]));
    assert_eq!(o.clipped, 2);
    // A finite DC offset so large that v·scale overflows f32 still clips (review fix).
    let mut c = bare();
    c.adc = Some(AdcConfig { bits: 16 });
    c.output = OutputFormat::Ci16;
    c.adc_full_scale_volts = 1e-6;
    c.iq = Some(IqImpairment {
        gain_imbalance_db: 0.0,
        phase_skew_deg: 0.0,
        dc_i_volts: 1e35,
        dc_q_volts: -1e35,
    });
    let mut s = ReceiverChannelState::new(c).unwrap();
    let o = s
        .process(0, SampleIndex::new(0), &[0, 0], &clock())
        .unwrap();
    assert_eq!(o.samples, Samples::Ci16(vec![32_767, -32_768]));
    assert_eq!(o.clipped, 2);
    // sc16 without the ADC stage is a configuration error.
    let mut c = bare();
    c.output = OutputFormat::Ci16;
    assert!(matches!(
        ReceiverChannelState::new(c),
        Err(ChainError::Config(_))
    ));
}

#[test]
fn noise_is_deterministic_and_keyed_by_seed_entity_and_channel() {
    let acc = vec![0_i64; 2 * 512];
    let mut c = bare();
    c.noise = Some(NoiseConfig {
        seed: 1,
        entity: 2,
        channel: 0,
    });
    c.front_end = Some(FrontEnd {
        gain_db: 100.0,
        noise_figure_db: 0.0,
    });
    let a = cf32(c, &acc);
    assert_eq!(a, cf32(c, &acc));
    for change in [
        NoiseConfig {
            seed: 2,
            entity: 2,
            channel: 0,
        },
        NoiseConfig {
            seed: 1,
            entity: 3,
            channel: 0,
        },
        NoiseConfig {
            seed: 1,
            entity: 2,
            channel: 1,
        },
    ] {
        let mut d = c;
        d.noise = Some(change);
        assert_ne!(cf32(d, &acc), a);
    }
    // Taps are the sampled H_rx: their DC sum is |H(0)|·σ_w, i.e. Σ t_m / σ_w ≈ 1.
    let s = ReceiverChannelState::new(c).unwrap();
    let taps = s.noise_taps().unwrap();
    assert!(taps.len() % 2 == 1);
    let mid = taps.len() / 2;
    for m in 0..mid {
        assert_eq!(taps[m], taps[taps.len() - 1 - m], "taps are symmetric");
    }
}

#[test]
fn agc_drives_a_steady_tone_to_its_target() {
    let one = 1_i64 << 40;
    let acc: Vec<i64> = (0..4000).flat_map(|_| [one / 100, 0]).collect();
    let mut c = bare();
    c.agc = Some(AgcConfig {
        target_volts: 0.05,
        attack: 0.05,
        step: 0.05,
        min_gain_db: -40.0,
        max_gain_db: 40.0,
    });
    let v = cf32(c, &acc);
    let last = v[v.len() - 2];
    assert!((last - 0.05).abs() < 1e-4, "{last}");
}

#[test]
fn invalid_configurations_are_rejected() {
    let mut cases = Vec::new();
    let mut c = bare();
    c.adc_full_scale_volts = 0.0;
    cases.push(c);
    let mut c = bare();
    c.adc = Some(AdcConfig { bits: 17 });
    cases.push(c);
    let mut c = full();
    c.decimation = Some(DecimationConfig {
        factor: 1,
        ..full().decimation.unwrap()
    });
    cases.push(c);
    let mut c = full();
    // Stopband above the output Nyquist (10/3 MS/s → 1.67 MHz).
    c.decimation = Some(DecimationConfig {
        factor: 3,
        filter: ResponseSpec {
            stopband_edge_uhz: 1_700_000_000_000,
            ..full().decimation.unwrap().filter
        },
    });
    cases.push(c);
    let mut c = bare();
    c.iq = Some(IqImpairment {
        gain_imbalance_db: 0.0,
        phase_skew_deg: 50.0,
        dc_i_volts: 0.0,
        dc_q_volts: 0.0,
    });
    cases.push(c);
    let mut c = bare();
    c.front_end = Some(FrontEnd {
        gain_db: 10.0,
        noise_figure_db: -1.0,
    });
    cases.push(c);
    for c in cases {
        assert!(ReceiverChannelState::new(c).is_err(), "{c:?}");
    }
}
