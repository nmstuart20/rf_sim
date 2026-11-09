//! The scenario file's top-level types.
//!
//! Every table rejects unknown keys. Fields marked as placeholders hold only the
//! choices that exist today (for example `clock = { model = "ideal" }`); later
//! prompts add variants and fields without changing the surrounding layout.

use iqsim_units::{DataRate, Duration, Frequency, Gain, Position, Power, SampleRate};
use schemars::JsonSchema;
use serde::Deserialize;
use toml::Spanned;

use crate::epoch::Epoch;
use crate::ident::{Ident, LibraryRef};

/// The scenario format version this build reads.
pub const SCENARIO_SCHEMA_VERSION: u32 = 1;

/// A scenario file.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Scenario {
    /// Scenario format version. This build reads version 1.
    #[schemars(with = "u32", extend("const" = SCENARIO_SCHEMA_VERSION))]
    pub schema_version: Spanned<u32>,
    /// Descriptive information.
    pub metadata: Metadata,
    /// UTC instant of truth time zero.
    pub epoch: Epoch,
    /// Seed for every random stream (0 to 2⁶³−1, the TOML integer range).
    pub seed: u64,
    /// How the run is paced.
    pub run_mode: RunMode,
    /// Render settings for `sim run` (duration and block length).
    pub render: Option<Render>,
    /// Where library items are looked up.
    #[serde(default)]
    pub library: Library,
    /// Platforms: anything with a position that carries emitters or receivers.
    #[serde(default)]
    pub platforms: Vec<Platform>,
    /// RF emitters.
    #[serde(default)]
    pub emitters: Vec<Emitter>,
    /// Virtual receivers and their channels.
    #[serde(default)]
    pub receivers: Vec<Receiver>,
    /// Antenna-port outputs (the RF field at an antenna, before any receiver).
    #[serde(default)]
    pub antenna_ports: Vec<AntennaPort>,
    /// Output sinks that channels and antenna ports stream to.
    #[serde(default)]
    pub sinks: Vec<Sink>,
    /// Shared hardware resources (NIC, disk, USB, PCIe) that sinks use.
    #[serde(default)]
    pub resources: Vec<Resource>,
}

/// Descriptive information about the scenario.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Metadata {
    /// Short scenario name.
    #[schemars(with = "String")]
    pub name: Spanned<String>,
    /// Free-text description.
    pub description: Option<String>,
}

/// How a run is paced.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum RunMode {
    /// As fast as possible; no deadlines.
    Offline,
    /// Paced to the wall clock, with sink deadlines and resource budgets.
    Realtime,
}

/// Render settings (`sim run`).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Render {
    /// How much truth time to render, from the epoch.
    #[schemars(with = "Duration")]
    pub duration: Spanned<Duration>,
    /// Block length (default 1 ms). It must be a whole number of samples at every
    /// channel and antenna-port rate.
    #[schemars(with = "Option<Duration>")]
    pub block: Option<Spanned<Duration>>,
}

/// Details of a pre-sampling tap response (H_rx or B_out), beyond its bandwidth.
///
/// The response is a linear-phase lowpass with passband edge `bandwidth / 2` and
/// stopband edge `stopband / 2` around the centre (contracts §5.1).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ResponseDetail {
    /// Full width at which the stopband begins (default: the complex sample rate).
    /// It must exceed the bandwidth and be at most the sample rate.
    #[schemars(with = "Option<Frequency>")]
    pub stopband: Option<Spanned<Frequency>>,
    /// Minimum stopband attenuation (default 80 dB, 20 to 200 dB).
    #[schemars(with = "Option<Gain>")]
    pub attenuation: Option<Spanned<Gain>>,
    /// Group delay (default 0, at most 1 ms).
    #[schemars(with = "Option<Duration>")]
    pub group_delay: Option<Spanned<Duration>>,
}

/// Front-end gain and noise figure (receiver chain stage 1).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FrontEndSpec {
    /// Power gain from the receiver input to the ADC input.
    pub gain: Gain,
    /// Noise figure (default 0 dB).
    pub noise_figure: Option<Gain>,
}

/// IQ imbalance and DC offset (receiver chain stage 4).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IqSpec {
    /// Q-branch gain relative to I (default 0 dB).
    pub gain_imbalance: Option<Gain>,
    /// Quadrature skew (default 0°, below 45°).
    pub phase_skew: Option<iqsim_units::Angle>,
    /// DC offset on I as a fraction of the ADC full scale (default 0).
    #[serde(default)]
    pub dc_i: f64,
    /// DC offset on Q as a fraction of the ADC full scale (default 0).
    #[serde(default)]
    pub dc_q: f64,
}

/// Digital decimation (receiver chain stage 5).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct DecimationSpec {
    /// Decimation factor, at least 2.
    pub factor: u32,
    /// Passband full width (default 0.8 × the output rate).
    pub bandwidth: Option<Frequency>,
    /// Stopband full width (default the output rate).
    pub stopband: Option<Frequency>,
    /// Stopband attenuation (default 80 dB).
    pub attenuation: Option<Gain>,
}

/// Automatic gain control (receiver chain stage 6, off unless given).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AgcSpec {
    /// Target RMS level relative to the ADC full scale (default −12 dB).
    pub target: Option<Gain>,
    /// Power-tracker coefficient, 0 < attack ≤ 1 (default 0.01).
    pub attack: Option<f64>,
    /// Gain step, 0 < step ≤ 1 (default 0.01).
    pub step: Option<f64>,
    /// Smallest AGC gain (default −40 dB).
    pub min_gain: Option<Gain>,
    /// Largest AGC gain (default 40 dB).
    pub max_gain: Option<Gain>,
}

/// The virtual ADC (receiver chain stage 7).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AdcSpec {
    /// Full scale: the power of a CW tone whose components just reach the clip level
    /// (default 0 dBm).
    pub full_scale: Option<Power>,
    /// Resolution, 2 to 16 bits (default 16).
    pub bits: Option<u32>,
    /// Quantize and clip (default true). Off, the channel writes unquantized cf32.
    pub enabled: Option<bool>,
}

/// The library search path.
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Library {
    /// Directories searched in order for library items, relative to the scenario
    /// file unless absolute.
    #[serde(default)]
    pub paths: Vec<String>,
}

/// A platform: a body with a position, trajectory and antenna mounts.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Platform {
    /// Unique entity name.
    pub name: Ident,
    /// Initial position.
    pub position: Position,
    /// Motion model (placeholder: only `static`).
    #[serde(default)]
    pub trajectory: Trajectory,
    /// Antenna mounts on the platform body.
    #[serde(default)]
    pub mounts: Vec<Mount>,
}

/// A platform motion model (placeholder).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Trajectory {
    /// The trajectory kind.
    pub kind: TrajectoryKind,
}

/// Available motion models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum TrajectoryKind {
    /// The platform stays at its position.
    #[default]
    Static,
}

/// An antenna mount on a platform (placeholder: lever arm and orientation come later).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Mount {
    /// Name, unique within the platform.
    pub name: Ident,
    /// Antenna pattern from the library.
    pub antenna: LibraryRef,
}

/// An RF emitter on a platform.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Emitter {
    /// Unique entity name.
    pub name: Ident,
    /// The platform carrying the emitter.
    pub platform: Ident,
    /// The mount the emitter radiates from (default: the platform origin).
    pub mount: Option<Ident>,
    /// Carrier frequency.
    #[schemars(with = "Frequency")]
    pub frequency: Spanned<Frequency>,
    /// Transmit power at the antenna input.
    pub power: Power,
    /// Behavior from the library (placeholder reference).
    pub behavior: LibraryRef,
    /// Waveform from the library (placeholder reference).
    pub waveform: LibraryRef,
}

/// A virtual receiver: an antenna on a platform feeding one or more channels.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Receiver {
    /// Unique entity name.
    pub name: Ident,
    /// The platform carrying the receiver antenna.
    pub platform: Ident,
    /// The mount the receiver antenna sits on (default: the platform origin).
    pub mount: Option<Ident>,
    /// Receiver channels (at least one).
    pub channels: Vec<Channel>,
}

/// One tuned channel of a virtual receiver.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Channel {
    /// Name, unique within the receiver.
    pub name: Ident,
    /// Tuned centre frequency.
    #[schemars(with = "Frequency")]
    pub center_frequency: Spanned<Frequency>,
    /// Complex sample rate.
    #[schemars(with = "SampleRate")]
    pub sample_rate: Spanned<SampleRate>,
    /// Passband width of the receiver response H_rx (below its stopband).
    #[schemars(with = "Frequency")]
    pub bandwidth: Spanned<Frequency>,
    /// H_rx details beyond the bandwidth.
    #[serde(default)]
    pub response: ResponseDetail,
    /// Sampling clock model.
    #[serde(default)]
    pub clock: ClockModelSpec,
    /// Local-oscillator model.
    #[serde(default)]
    pub lo: LoModelSpec,
    /// Front-end gain and noise figure (off when absent: 0 dB gain, 0 dB noise figure).
    pub front_end: Option<FrontEndSpec>,
    /// Thermal noise kTBF (default true).
    pub thermal_noise: Option<bool>,
    /// IQ imbalance and DC offset (off when absent).
    pub iq: Option<IqSpec>,
    /// Digital decimation (off when absent).
    pub decimation: Option<DecimationSpec>,
    /// Automatic gain control (off when absent).
    pub agc: Option<AgcSpec>,
    /// The virtual ADC (defaults when absent: 0 dBm full scale, 16 bits, on).
    pub adc: Option<AdcSpec>,
    /// Sinks this channel streams to.
    #[serde(default)]
    pub sinks: Vec<Ident>,
}

/// A sampling-clock model.
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ClockModelSpec {
    /// The clock model.
    pub model: ClockModel,
    /// Sampling-instant offset for `fixed_offset`: every sample is acquired this much
    /// later than nominal.
    #[schemars(with = "Option<Duration>")]
    pub offset: Option<Spanned<Duration>>,
}

/// Available sampling-clock models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ClockModel {
    /// Samples exactly at the nominal instants (`IdealClockMap`).
    #[default]
    Ideal,
    /// Nominal instants shifted by a constant `offset` (an ideal map with a later epoch).
    FixedOffset,
}

/// A local-oscillator model.
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LoModelSpec {
    /// The LO model.
    pub model: LoModel,
    /// LO frequency error for `fixed_offset` (LO minus nominal).
    #[schemars(with = "Option<Frequency>")]
    pub offset: Option<Spanned<Frequency>>,
}

/// Available local-oscillator models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum LoModel {
    /// No frequency or phase error.
    #[default]
    Ideal,
    /// A constant frequency error `offset`; phase noise is not modelled yet.
    FixedOffset,
}

/// An antenna-port output: the field at an antenna, band-limited by B_out, with
/// no virtual-receiver effects (contracts §5).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AntennaPort {
    /// Unique entity name.
    pub name: Ident,
    /// The platform carrying the antenna.
    pub platform: Ident,
    /// The mount the antenna sits on (default: the platform origin).
    pub mount: Option<Ident>,
    /// Output centre frequency.
    #[schemars(with = "Frequency")]
    pub center_frequency: Spanned<Frequency>,
    /// Complex output sample rate.
    #[schemars(with = "SampleRate")]
    pub sample_rate: Spanned<SampleRate>,
    /// Bandwidth of the output response B_out (below its stopband).
    #[schemars(with = "Frequency")]
    pub bandwidth: Spanned<Frequency>,
    /// B_out details beyond the bandwidth.
    #[serde(default)]
    pub response: ResponseDetail,
    /// Output full scale: the power of a CW tone at sc16 code 32768 (default 0 dBm).
    pub full_scale: Option<Power>,
    /// Sinks this output streams to.
    #[serde(default)]
    pub sinks: Vec<Ident>,
}

/// An output sink (placeholder: endpoint and file settings come with the sinks).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Sink {
    /// Unique entity name.
    pub name: Ident,
    /// Sink type.
    pub kind: SinkKind,
    /// Sample format (default `ci16_le`).
    #[serde(default)]
    pub format: SampleFormatSpec,
    /// Shared resources this sink's traffic crosses.
    #[serde(default)]
    pub resources: Vec<Ident>,
}

/// Sample formats a sink can write.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum SampleFormatSpec {
    /// Complex int16, little-endian (sc16).
    #[default]
    Ci16Le,
    /// Complex float32, little-endian.
    Cf32Le,
}

/// Output sink types (plan "Interfaces", IQ outputs).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum SinkKind {
    /// SigMF recording.
    Sigmf,
    /// Strict, acknowledged framed ZMQ stream.
    ZmqFramed,
    /// VITA-49 over UDP.
    Vita49,
    /// Best-effort raw ZMQ PUB preview.
    ZmqPreview,
}

/// A shared hardware resource with an optional usable capacity (contracts §8).
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Resource {
    /// Unique entity name.
    pub name: Ident,
    /// Resource type.
    pub kind: ResourceKind,
    /// Usable capacity, e.g. "1 Gb/s". Omitted means unknown, reported as unverified.
    pub capacity: Option<DataRate>,
}

/// Shared resource types.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ResourceKind {
    /// A network interface.
    Nic,
    /// A disk or storage volume.
    Disk,
    /// A USB controller.
    Usb,
    /// A PCIe link or root port.
    Pcie,
}
