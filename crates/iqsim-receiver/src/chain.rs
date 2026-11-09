//! The virtual-receiver chain (plan "Receiver model", contracts §5.4).
//!
//! It runs once per channel, after the complete reduced path sum of the channel's
//! virtual-receiver target. H_rx is already applied to every path before sampling
//! (contracts §5.1), so the chain never filters the signal by H_rx again. Stages, in
//! order, each switched on by its `Option` in [`ChainConfig`]:
//!
//! 0. Accumulator → volts at the receiver input: `x = f32(acc) / lsb_per_volt`.
//! 1. Front-end gain ([`FrontEnd`]): `x ·= g`, g the voltage gain. Its noise figure
//!    sets the level of stage 2.
//! 2. Thermal noise ([`NoiseConfig`]): `x += Σ_m t_m·w[n + m]`, where `w` is the
//!    unit-power Philox noise of the channel's `ReceiverThermalNoise` stream at the
//!    channel sample index and `t_m` is the digital equivalent of H_rx scaled to the
//!    output-referred density `N₀ = k·T₀·F·g²`. The total is `N₀·B_n`, with B_n the
//!    noise bandwidth `∫|H_rx|² df`.
//! 3. LO frequency offset ([`LoConfig`]): `x ·= e^{−j2π·Δf·t_acq(n)}`, with the exact
//!    modular phase of `OffsetFrequency::phase_at` at the clock map's acquisition
//!    instant. Phase noise is a stub ([`PhaseNoise::None`]).
//! 4. IQ imbalance and DC offset ([`IqImpairment`]).
//! 5. Decimation ([`DecimationConfig`]): FIR lowpass and keep every M-th sample.
//! 6. AGC ([`AgcConfig`]), off by default.
//! 7. ADC ([`AdcConfig`]): quantization to `bits` and clipping, left-justified in sc16.
//!
//! Every stage is f32 with each operation written separately in a fixed order (no FMA,
//! no reassociation), so a later GPU port can reproduce it. Constants are computed once
//! on the control plane in f64 with `libm` and rounded to f32.
//!
//! # State ownership and block handoff
//!
//! One [`ReceiverChannelState`] owns a channel's stateful processes: the decimation
//! filter history and phase, and the AGC gain and power estimate. Blocks must arrive in
//! order (block b+1 after block b, contiguous sample indices); anything else is an error,
//! never a reset. Noise needs no carried state: it is a correlation over counter-indexed
//! white noise, so each block draws its own `L − 1` sample lookahead from the counter.
//! LO phase is computed from the acquisition instant of every sample. The output is
//! therefore identical however the run is split into blocks.

use std::fmt;

use iqsim_core::{
    ClockMap, NoiseStream, NumericError, OffsetFrequency, Phase, ResponseError, ResponseSpec,
    SampleIndex, SampleRate, StreamId, TapKind, TapResponse, TargetScale, TimeError,
    f32_to_i64_ties_even, sincos,
};

use crate::{BOLTZMANN_J_PER_K, REFERENCE_TEMPERATURE_K};

/// Front-end gain and noise figure (stage 1).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrontEnd {
    /// Power gain from the receiver input to the ADC input, dB.
    pub gain_db: f64,
    /// Noise figure, dB (≥ 0). Used by the noise stage.
    pub noise_figure_db: f64,
}

/// Thermal-noise stream identity (stage 2). The level comes from [`FrontEnd`]: with the
/// front end off, the gain is 0 dB and the noise figure 0 dB (F = 1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NoiseConfig {
    /// Scenario seed.
    pub seed: u64,
    /// Receiver entity id.
    pub entity: u32,
    /// Channel index within the receiver (the `ReceiverThermalNoise` stream instance).
    pub channel: u16,
}

/// LO phase-noise model. Only the stub exists; prompt 15 adds a Philox-driven process
/// with its own registered stream kind.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PhaseNoise {
    /// No phase noise.
    #[default]
    None,
}

/// LO frequency error (stage 3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LoConfig {
    /// LO frequency minus nominal, µHz. The output is shifted by −offset.
    pub offset_uhz: i64,
    /// Phase noise (stub).
    pub phase_noise: PhaseNoise,
}

/// IQ imbalance and DC offset (stage 4):
/// `I' = I + dc_i`, `Q' = a·(Q·cos φ − I·sin φ) + dc_q`, with `a = 10^(gain/20)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IqImpairment {
    /// Q-branch gain relative to I, dB.
    pub gain_imbalance_db: f64,
    /// Quadrature skew φ, degrees (|φ| < 45).
    pub phase_skew_deg: f64,
    /// DC offset on I at the ADC input, volts.
    pub dc_i_volts: f64,
    /// DC offset on Q at the ADC input, volts.
    pub dc_q_volts: f64,
}

/// Digital decimation (stage 5): FIR lowpass at the input rate, then keep the samples
/// whose input index is a multiple of the factor. Output sample j is input sample j·M
/// filtered causally, so the filter adds `(L − 1)/2` input samples of group delay.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DecimationConfig {
    /// Decimation factor M ≥ 2.
    pub factor: u32,
    /// Lowpass definition at the input rate. Its group delay must be 0 (the FIR is
    /// centred), and its stopband edge at most the output Nyquist frequency.
    pub filter: ResponseSpec,
}

/// Automatic gain control (stage 6), per output sample:
/// `y = g·x`, `p += α·(|y|² − p)`, `g += β·g·(T² − p)/T²`, then g is clamped.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AgcConfig {
    /// Target RMS magnitude T at the ADC input, volts.
    pub target_volts: f64,
    /// Power-tracker coefficient α, 0 < α ≤ 1.
    pub attack: f64,
    /// Gain step β, 0 < β ≤ 1.
    pub step: f64,
    /// Smallest AGC gain, dB.
    pub min_gain_db: f64,
    /// Largest AGC gain, dB.
    pub max_gain_db: f64,
}

/// ADC quantization and clipping (stage 7).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AdcConfig {
    /// Resolution, 2 to 16 bits. Codes are left-justified in sc16.
    pub bits: u32,
}

/// Output sample format of a channel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    /// sc16 ADC codes (needs the ADC stage).
    Ci16,
    /// f32 normalized to the ADC full scale: the code / 32768 with the ADC on, or the
    /// unquantized `v / V_fs` with it off.
    Cf32,
}

/// The configuration of one channel's chain. Each `Option` is a stage switch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChainConfig {
    /// Channel sample rate (the virtual target's render rate).
    pub rate: SampleRate,
    /// The channel's H_rx, used to shape the noise.
    pub response: TapResponse,
    /// The virtual target's accumulator scale.
    pub input_scale: TargetScale,
    /// ADC full scale: the per-component clip voltage. A tone of power
    /// `V_fs²/R₀` just reaches it.
    pub adc_full_scale_volts: f64,
    /// Output format.
    pub output: OutputFormat,
    /// Stage 1.
    pub front_end: Option<FrontEnd>,
    /// Stage 2.
    pub noise: Option<NoiseConfig>,
    /// Stage 3.
    pub lo: Option<LoConfig>,
    /// Stage 4.
    pub iq: Option<IqImpairment>,
    /// Stage 5.
    pub decimation: Option<DecimationConfig>,
    /// Stage 6.
    pub agc: Option<AgcConfig>,
    /// Stage 7.
    pub adc: Option<AdcConfig>,
}

/// A chain failure. Every variant is a hard error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChainError {
    /// A configuration value is invalid.
    Config(&'static str),
    /// A filter response is invalid.
    Response(ResponseError),
    /// A numeric conversion failed.
    Numeric(NumericError),
    /// A clock or index operation failed.
    Time(TimeError),
    /// A block arrived out of order.
    OutOfOrder {
        /// The block the state expected.
        expected: u64,
        /// The block it got.
        got: u64,
    },
    /// A block's first sample does not continue the previous block.
    Gap {
        /// The sample index the state expected.
        expected: u64,
        /// The sample index it got.
        got: u64,
    },
    /// The accumulator slice has odd length.
    Length,
    /// A stage produced a non-finite value.
    NonFinite,
}

impl fmt::Display for ChainError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Config(what) => write!(f, "invalid receiver chain: {what}"),
            Self::Response(e) => write!(f, "invalid receiver filter: {e}"),
            Self::Numeric(e) => write!(f, "receiver numeric error: {e}"),
            Self::Time(e) => write!(f, "receiver clock error: {e}"),
            Self::OutOfOrder { expected, got } => {
                write!(f, "receiver block {got} arrived, expected block {expected}")
            }
            Self::Gap { expected, got } => {
                write!(
                    f,
                    "receiver block starts at sample {got}, expected {expected}"
                )
            }
            Self::Length => f.write_str("receiver input must be interleaved I/Q pairs"),
            Self::NonFinite => f.write_str("receiver stage produced a non-finite value"),
        }
    }
}

impl std::error::Error for ChainError {}

impl From<NumericError> for ChainError {
    fn from(e: NumericError) -> Self {
        Self::Numeric(e)
    }
}

impl From<TimeError> for ChainError {
    fn from(e: TimeError) -> Self {
        Self::Time(e)
    }
}

impl From<ResponseError> for ChainError {
    fn from(e: ResponseError) -> Self {
        Self::Response(e)
    }
}

/// Output samples of one block, interleaved I/Q.
#[derive(Clone, Debug, PartialEq)]
pub enum Samples {
    /// sc16.
    Ci16(Vec<i16>),
    /// cf32.
    Cf32(Vec<f32>),
}

impl Samples {
    /// Number of complex samples.
    pub fn len(&self) -> usize {
        match self {
            Self::Ci16(v) => v.len() / 2,
            Self::Cf32(v) => v.len() / 2,
        }
    }

    /// True when there are no samples.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// One block of channel output.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelOutput {
    /// Output-rate index of the first sample.
    pub first_index: u64,
    /// The samples.
    pub samples: Samples,
    /// Components clipped by the ADC in this block.
    pub clipped: u64,
}

fn db_to_voltage(db: f64) -> f64 {
    libm::exp10(db / 20.0)
}

fn finite_f32(x: f64, what: &'static str) -> Result<f32, ChainError> {
    let y = x as f32;
    if y.is_finite() {
        Ok(y)
    } else {
        Err(ChainError::Config(what))
    }
}

/// The digital equivalent of a response: `t_m = scale·h(τ_g + (m − c)·T)·T` for
/// m = 0..L, with L = 2·⌈T_h/T⌉ + 1 and c = (L − 1)/2, so its DTFT is H(f) up to the
/// response's own truncation and the aliasing of |H| beyond F_s/2 (below the stopband).
fn sampled_response(resp: &TapResponse, rate_hz: f64, scale: f64) -> Result<Vec<f32>, ChainError> {
    let t = 1.0 / rate_hz;
    let half = libm::ceil(resp.half_width_s() * rate_hz);
    if !(half.is_finite() && half <= 1_000_000.0) {
        return Err(ChainError::Config("filter is longer than 2·10⁶ taps"));
    }
    let half = half as i64;
    let mut taps = Vec::with_capacity((2 * half + 1) as usize);
    for m in -half..=half {
        let tm = resp.group_delay_s() + (m as f64) * t;
        let v = resp.impulse(tm) * t;
        let v = v * scale;
        taps.push(finite_f32(v, "filter tap is not finite")?);
    }
    Ok(taps)
}

/// Prepared constants of a chain.
#[derive(Clone, Debug)]
struct Prepared {
    lsb_per_volt: f32,
    gain: Option<f32>,
    noise: Option<(NoiseStream, Vec<f32>)>,
    lo: Option<OffsetFrequency>,
    iq: Option<[f32; 5]>,
    decim: Option<(u32, Vec<f32>)>,
    agc: Option<[f32; 5]>,
    adc: Option<(u32, f32)>,
    inv_fs: f32,
    output: OutputFormat,
    output_rate: SampleRate,
}

/// The owner of one channel's receiver state. Blocks go through [`Self::process`] in
/// order.
#[derive(Clone, Debug)]
pub struct ReceiverChannelState {
    config: ChainConfig,
    p: Prepared,
    next_block: u64,
    next_index: u64,
    /// Decimation history: the last L − 1 inputs, oldest first.
    history: Vec<[f32; 2]>,
    agc_gain: f32,
    agc_power: f32,
    clipped_total: u64,
}

impl ReceiverChannelState {
    /// Validates the configuration and prepares every constant. The state starts at
    /// block 0, sample 0, with zero filter history and the AGC at unity gain.
    pub fn new(config: ChainConfig) -> Result<Self, ChainError> {
        let rate_hz = config.rate.num() as f64 / config.rate.den() as f64;
        let v_fs = config.adc_full_scale_volts;
        if !(v_fs.is_finite() && v_fs > 0.0) {
            return Err(ChainError::Config("ADC full scale must be positive"));
        }
        if config.output == OutputFormat::Ci16 && config.adc.is_none() {
            return Err(ChainError::Config("sc16 output needs the ADC stage"));
        }
        let (g_lin, f_lin) = match config.front_end {
            Some(fe) => {
                if !(fe.gain_db.is_finite() && fe.gain_db.abs() <= 200.0) {
                    return Err(ChainError::Config("front-end gain must be within ±200 dB"));
                }
                if !(fe.noise_figure_db.is_finite() && (0.0..=100.0).contains(&fe.noise_figure_db))
                {
                    return Err(ChainError::Config("noise figure must be 0 to 100 dB"));
                }
                (
                    db_to_voltage(fe.gain_db),
                    libm::exp10(fe.noise_figure_db / 10.0),
                )
            }
            None => (1.0, 1.0),
        };
        let gain = match config.front_end {
            Some(_) => Some(finite_f32(g_lin, "front-end gain")?),
            None => None,
        };
        let noise = match config.noise {
            Some(n) => {
                // Output-referred density N₀ = k·T₀·F·g², white over F_s before shaping.
                let n0 = BOLTZMANN_J_PER_K * REFERENCE_TEMPERATURE_K * f_lin * (g_lin * g_lin);
                let sigma_w =
                    (n0 * rate_hz * iqsim_core::amplitude::REFERENCE_IMPEDANCE_OHMS).sqrt();
                let taps = sampled_response(&config.response, rate_hz, sigma_w)?;
                let stream = NoiseStream::new(
                    n.seed,
                    n.entity,
                    StreamId::receiver_thermal_noise(n.channel),
                );
                Some((stream, taps))
            }
            None => None,
        };
        let lo = match config.lo {
            Some(lo) => {
                let PhaseNoise::None = lo.phase_noise;
                Some(OffsetFrequency::from_uhz(lo.offset_uhz)?)
            }
            None => None,
        };
        let iq = match config.iq {
            Some(iq) => {
                if !(iq.phase_skew_deg.is_finite() && iq.phase_skew_deg.abs() < 45.0) {
                    return Err(ChainError::Config("IQ phase skew must be below 45°"));
                }
                if !(iq.gain_imbalance_db.is_finite() && iq.gain_imbalance_db.abs() <= 20.0) {
                    return Err(ChainError::Config(
                        "IQ gain imbalance must be within ±20 dB",
                    ));
                }
                let a = db_to_voltage(iq.gain_imbalance_db);
                let phi = iq.phase_skew_deg * (core::f64::consts::PI / 180.0);
                let (s, c) = libm::sincos(phi);
                Some([
                    finite_f32(a * c, "IQ coefficient")?,
                    finite_f32(a * s, "IQ coefficient")?,
                    finite_f32(iq.dc_i_volts, "DC offset")?,
                    finite_f32(iq.dc_q_volts, "DC offset")?,
                    0.0,
                ])
            }
            None => None,
        };
        let (decim, output_rate) = match config.decimation {
            Some(d) => {
                if d.factor < 2 {
                    return Err(ChainError::Config("decimation factor must be at least 2"));
                }
                if d.filter.group_delay_ps != 0 {
                    return Err(ChainError::Config(
                        "decimation filter group delay must be 0",
                    ));
                }
                let out = SampleRate::from_ratio(
                    config.rate.num(),
                    config
                        .rate
                        .den()
                        .checked_mul(u64::from(d.factor))
                        .ok_or(ChainError::Config("decimated rate out of range"))?,
                )
                .map_err(|_| ChainError::Config("decimated rate out of range"))?;
                // Stopband edge ≤ output Nyquist: 2·f_s·den_out ≤ 10⁶·num_out, exactly.
                if 2 * i128::from(d.filter.stopband_edge_uhz) * i128::from(out.den())
                    > 1_000_000 * i128::from(out.num())
                {
                    return Err(ChainError::Config(
                        "decimation stopband edge is above the output Nyquist frequency",
                    ));
                }
                let resp = TapResponse::new(
                    TapKind::VirtualReceiver,
                    d.filter,
                    config.rate,
                    config.response.center_uhz(),
                )?;
                (
                    Some((d.factor, sampled_response(&resp, rate_hz, 1.0)?)),
                    out,
                )
            }
            None => (None, config.rate),
        };
        let agc = match config.agc {
            Some(a) => {
                let ok = a.target_volts.is_finite()
                    && a.target_volts > 0.0
                    && a.attack > 0.0
                    && a.attack <= 1.0
                    && a.step > 0.0
                    && a.step <= 1.0
                    && a.min_gain_db.is_finite()
                    && a.max_gain_db.is_finite()
                    && a.min_gain_db <= a.max_gain_db
                    && a.max_gain_db.abs() <= 200.0
                    && a.min_gain_db.abs() <= 200.0;
                if !ok {
                    return Err(ChainError::Config("invalid AGC settings"));
                }
                let t2 = a.target_volts * a.target_volts;
                Some([
                    finite_f32(t2, "AGC target")?,
                    finite_f32(a.attack, "AGC attack")?,
                    finite_f32(a.step / t2, "AGC step")?,
                    finite_f32(db_to_voltage(a.min_gain_db), "AGC gain")?,
                    finite_f32(db_to_voltage(a.max_gain_db), "AGC gain")?,
                ])
            }
            None => None,
        };
        let adc = match config.adc {
            Some(adc) => {
                if !(2..=16).contains(&adc.bits) {
                    return Err(ChainError::Config("ADC resolution must be 2 to 16 bits"));
                }
                let codes = f64::from(1_u32 << (adc.bits - 1));
                Some((adc.bits, finite_f32(codes / v_fs, "ADC scale")?))
            }
            None => None,
        };
        let p = Prepared {
            lsb_per_volt: config.input_scale.lsb_per_volt(),
            gain,
            noise,
            lo,
            iq,
            decim,
            agc,
            adc,
            inv_fs: finite_f32(1.0 / v_fs, "ADC full scale")?,
            output: config.output,
            output_rate,
        };
        let hist = p.decim.as_ref().map_or(0, |(_, t)| t.len() - 1);
        let agc_power = p.agc.map_or(0.0, |a| a[0]);
        Ok(Self {
            config,
            p,
            next_block: 0,
            next_index: 0,
            history: vec![[0.0, 0.0]; hist],
            agc_gain: 1.0,
            agc_power,
            clipped_total: 0,
        })
    }

    /// The configuration.
    pub fn config(&self) -> &ChainConfig {
        &self.config
    }

    /// The output sample rate (the channel rate divided by any decimation factor).
    pub fn output_rate(&self) -> SampleRate {
        self.p.output_rate
    }

    /// The noise-shaping taps (unit-power white noise in, volts out), if noise is on.
    pub fn noise_taps(&self) -> Option<&[f32]> {
        self.p.noise.as_ref().map(|(_, t)| t.as_slice())
    }

    /// The decimation filter taps, if decimation is on.
    pub fn decimation_taps(&self) -> Option<&[f32]> {
        self.p.decim.as_ref().map(|(_, t)| t.as_slice())
    }

    /// Components clipped so far.
    pub fn clipped_total(&self) -> u64 {
        self.clipped_total
    }

    /// The next block this state accepts.
    pub fn next_block(&self) -> u64 {
        self.next_block
    }

    /// Processes block `block`, whose first channel sample is `first`, from the
    /// complete reduced accumulator sum `acc` (I/Q interleaved). `clock` is the
    /// channel's ClockMap, which gives the LO stage its acquisition instants.
    pub fn process(
        &mut self,
        block: u64,
        first: SampleIndex,
        acc: &[i64],
        clock: &dyn ClockMap,
    ) -> Result<ChannelOutput, ChainError> {
        if block != self.next_block {
            return Err(ChainError::OutOfOrder {
                expected: self.next_block,
                got: block,
            });
        }
        if first.get() != self.next_index {
            return Err(ChainError::Gap {
                expected: self.next_index,
                got: first.get(),
            });
        }
        if !acc.len().is_multiple_of(2) {
            return Err(ChainError::Length);
        }
        let n = acc.len() / 2;
        let p = &self.p;

        // Stage 0: volts at the receiver input.
        let mut x: Vec<[f32; 2]> = acc
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| [c[0] as f32 / p.lsb_per_volt, c[1] as f32 / p.lsb_per_volt])
            .collect();

        // Stage 1: gain.
        if let Some(g) = p.gain {
            for v in &mut x {
                v[0] *= g;
                v[1] *= g;
            }
        }

        // Stage 2: shaped thermal noise, a correlation over counter-indexed white noise.
        if let Some((stream, taps)) = &p.noise {
            let l = taps.len();
            let mut w = Vec::with_capacity(n + l - 1);
            for j in 0..(n + l - 1) as u64 {
                let (a, b) = stream.sample(first.checked_add(j)?);
                w.push([a, b]);
            }
            for (k, v) in x.iter_mut().enumerate() {
                let (mut si, mut sq) = (0.0_f32, 0.0_f32);
                for (m, t) in taps.iter().enumerate() {
                    let s = w[k + m];
                    si += *t * s[0];
                    sq += *t * s[1];
                }
                v[0] += si;
                v[1] += sq;
            }
        }

        // Stage 3: LO offset at the actual acquisition instants.
        if let Some(off) = p.lo {
            for (k, v) in x.iter_mut().enumerate() {
                let t = clock.instant(first.checked_add(k as u64)?)?;
                let phi = Phase::from_bits(0).wrapping_sub(off.phase_at(t));
                let (s, c) = sincos(phi);
                let re = v[0] * c - v[1] * s;
                let im = v[0] * s + v[1] * c;
                *v = [re, im];
            }
        }

        // Stage 4: IQ imbalance, then DC.
        if let Some([ac, as_, dci, dcq, _]) = p.iq {
            for v in &mut x {
                let q = v[1] * ac - v[0] * as_;
                v[0] += dci;
                v[1] = q + dcq;
            }
        }

        // Stage 5: decimation with carried history.
        let mut first_out = first.get();
        if let Some((m, taps)) = &p.decim {
            let m = u64::from(*m);
            let h = self.history.len();
            let mut buf = Vec::with_capacity(h + n);
            buf.extend_from_slice(&self.history);
            buf.extend_from_slice(&x);
            first_out = first.get().div_ceil(m);
            let mut out = Vec::new();
            let mut idx = first_out
                .checked_mul(m)
                .ok_or(ChainError::Time(TimeError::IndexOverflow))?;
            let end = first.checked_add(n as u64)?.get();
            while idx < end {
                // Input idx sits at buf[h + idx − first]; tap m' multiplies x[idx − m'].
                let pos = h + (idx - first.get()) as usize;
                let (mut si, mut sq) = (0.0_f32, 0.0_f32);
                for (k, t) in taps.iter().enumerate() {
                    let s = buf[pos - k];
                    si += *t * s[0];
                    sq += *t * s[1];
                }
                out.push([si, sq]);
                idx += m;
            }
            self.history.copy_from_slice(&buf[buf.len() - h..]);
            x = out;
        }

        // Stage 6: AGC.
        if let Some([t2, alpha, beta, gmin, gmax]) = p.agc {
            let (mut g, mut pw) = (self.agc_gain, self.agc_power);
            for v in &mut x {
                let yi = v[0] * g;
                let yq = v[1] * g;
                let e = yi * yi + yq * yq;
                pw += alpha * (e - pw);
                let corr = beta * (t2 - pw);
                g += g * corr;
                g = g.clamp(gmin, gmax);
                *v = [yi, yq];
            }
            self.agc_gain = g;
            self.agc_power = pw;
        }

        // Stage 7: ADC, then the output format.
        let mut clipped = 0_u64;
        let samples = match p.adc {
            Some((bits, scale)) => {
                let hi = (1_i64 << (bits - 1)) - 1;
                let lo = -(1_i64 << (bits - 1));
                let shift = 16 - bits;
                let mut codes = Vec::with_capacity(2 * x.len());
                for v in &x {
                    for c in v {
                        if !c.is_finite() {
                            return Err(ChainError::NonFinite);
                        }
                        // A finite input times a finite positive scale is never NaN; it
                        // can overflow to ±∞, which the clamp below turns into a clip.
                        let y = *c * scale;
                        // Clamp to one code past each end first (exact in f32), so a
                        // value far out of range clips instead of failing conversion.
                        let y = y.clamp((lo - 1) as f32, (hi + 1) as f32);
                        let q = f32_to_i64_ties_even(y)?;
                        let q = if q > hi {
                            clipped += 1;
                            hi
                        } else if q < lo {
                            clipped += 1;
                            lo
                        } else {
                            q
                        };
                        codes.push((q << shift) as i16);
                    }
                }
                match p.output {
                    OutputFormat::Ci16 => Samples::Ci16(codes),
                    OutputFormat::Cf32 => Samples::Cf32(
                        codes
                            .iter()
                            .map(|&c| f32::from(c) * (1.0 / 32_768.0))
                            .collect(),
                    ),
                }
            }
            None => {
                let mut out = Vec::with_capacity(2 * x.len());
                for v in &x {
                    for c in v {
                        let y = *c * p.inv_fs;
                        if !y.is_finite() {
                            return Err(ChainError::NonFinite);
                        }
                        out.push(y);
                    }
                }
                Samples::Cf32(out)
            }
        };

        self.next_block += 1;
        self.next_index = first.checked_add(n as u64)?.get();
        self.clipped_total += clipped;
        Ok(ChannelOutput {
            first_index: first_out,
            samples,
            clipped,
        })
    }
}
