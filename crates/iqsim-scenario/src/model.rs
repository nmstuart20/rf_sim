//! The scenario file's top-level types.
//!
//! Every table rejects unknown keys. Fields marked as placeholders hold only the
//! choices that exist today (for example `clock = { model = "ideal" }`); later
//! prompts add variants and fields without changing the surrounding layout.

use iqsim_units::{
    Acceleration, Altitude, Angle, AngularRate, DataRate, Duration, Frequency, Gain, Jerk,
    Latitude, Longitude, Position, Power, Pressure, SampleRate, Speed,
};
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
    /// Earth model: geoid grid and local atmosphere (MSL heights, flight levels and
    /// Mode C altitude).
    #[serde(default)]
    pub world: World,
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

/// Earth model settings (contracts §12.5, §12.6).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct World {
    /// Geoid grid. Without one, the geoid height is 0 everywhere (MSL is the
    /// ellipsoid) and the run manifest carries a warning.
    pub geoid: Option<GeoidSpec>,
    /// Local sea-level pressure minus 1013.25 hPa (QNH − standard), default 0. Shifts
    /// pressure altitude (Mode C) and the geometric height of flight levels.
    #[schemars(with = "Option<Pressure>")]
    pub pressure_offset: Option<Spanned<Pressure>>,
}

/// A geoid grid file.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GeoidSpec {
    /// A GeographicLib PGM geoid grid (e.g. `egm96-5.pgm`), relative to the scenario
    /// file unless absolute.
    #[schemars(with = "String")]
    pub path: Spanned<String>,
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
    /// Initial position (required, except for `track` trajectories, which start at
    /// their first point).
    pub position: Option<Position>,
    /// Motion model (default `static`).
    #[serde(default)]
    pub trajectory: Trajectory,
    /// Antenna mounts on the platform body.
    #[serde(default)]
    pub mounts: Vec<Mount>,
}

/// A platform motion model (contracts §13). Which fields apply depends on `kind`;
/// `sim lint` rejects the others.
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Trajectory {
    /// The trajectory kind (default `static`).
    #[serde(default)]
    pub kind: TrajectoryKind,
    /// Truth time of the initial state (maneuvers, route, ground; default 0 s).
    #[schemars(with = "Option<Duration>")]
    pub start: Option<Spanned<Duration>>,
    /// Initial track angle from true north (maneuvers and ground: required; route:
    /// default towards the first waypoint).
    #[schemars(with = "Option<Angle>")]
    pub heading: Option<Spanned<Angle>>,
    /// Initial speed (maneuvers, route, ground: required).
    #[schemars(with = "Option<Speed>")]
    pub speed: Option<Spanned<Speed>>,
    /// Fixed attitude of a static platform (default: nose north, level).
    #[schemars(with = "Option<AttitudeSpec>")]
    pub attitude: Option<Spanned<AttitudeSpec>>,
    /// Angle of attack added to the derived pitch (airborne kinds; default 0°).
    #[schemars(with = "Option<Angle>")]
    pub angle_of_attack: Option<Spanned<Angle>>,
    /// Constant wind, which crabs the derived heading (airborne kinds; default none).
    #[schemars(with = "Option<WindSpec>")]
    pub wind: Option<Spanned<WindSpec>>,
    /// Performance limits and route defaults.
    #[serde(default)]
    pub limits: LimitsSpec,
    /// Maneuver primitives, flown in order (maneuvers; ground).
    #[serde(default)]
    #[schemars(with = "Vec<ManeuverSpec>")]
    pub maneuvers: Vec<Spanned<ManeuverSpec>>,
    /// Waypoints, flown over in order (route; ground).
    #[serde(default)]
    #[schemars(with = "Vec<WaypointSpec>")]
    pub waypoints: Vec<Spanned<WaypointSpec>>,
    /// Recorded positions, in increasing time order (track; at least two).
    #[serde(default)]
    #[schemars(with = "Vec<TrackPointSpec>")]
    pub points: Vec<Spanned<TrackPointSpec>>,
    /// Track smoothing time constant τ (default 0 s, which interpolates the points).
    /// Motion faster than 1/τ rad/s is attenuated by 1/(1 + (ωτ)⁴).
    #[schemars(with = "Option<Duration>")]
    pub smoothing: Option<Spanned<Duration>>,
}

/// Available motion models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum TrajectoryKind {
    /// The platform stays at its position.
    #[default]
    Static,
    /// Maneuver primitives from an initial heading and speed.
    Maneuvers,
    /// A waypoint route.
    Route,
    /// A recorded track, smoothed to C2.
    Track,
    /// A ground vehicle, clamped to the ground (the ellipsoid until terrain exists):
    /// maneuvers without climbs, or waypoints at AGL heights.
    Ground,
}

/// A fixed attitude relative to local north-east-down.
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AttitudeSpec {
    /// Heading (yaw) from true north (default 0°).
    pub heading: Option<Angle>,
    /// Pitch, nose up positive (default 0°).
    pub pitch: Option<Angle>,
    /// Roll, right wing down positive (default 0°).
    pub roll: Option<Angle>,
}

/// A constant wind.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WindSpec {
    /// Wind speed.
    pub speed: Speed,
    /// Direction the wind blows from, from true north.
    pub from: Angle,
}

/// Performance limits (contracts §13.4 has the defaults).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LimitsSpec {
    /// Bank angle of route turns (default 25°, below 80°).
    #[schemars(with = "Option<Angle>")]
    pub bank: Option<Spanned<Angle>>,
    /// Roll rate (default 10 deg/s).
    #[schemars(with = "Option<AngularRate>")]
    pub roll_rate: Option<Spanned<AngularRate>>,
    /// Largest acceleration in speed changes (default 1 m/s²).
    #[schemars(with = "Option<Acceleration>")]
    pub accel: Option<Spanned<Acceleration>>,
    /// Jerk in speed changes (default 0.5 m/s³).
    #[schemars(with = "Option<Jerk>")]
    pub jerk: Option<Spanned<Jerk>>,
    /// Largest vertical acceleration (default 0.5 m/s²).
    #[schemars(with = "Option<Acceleration>")]
    pub vertical_accel: Option<Spanned<Acceleration>>,
    /// Vertical jerk (default 0.25 m/s³).
    #[schemars(with = "Option<Jerk>")]
    pub vertical_jerk: Option<Spanned<Jerk>>,
    /// Route climb and descent rate (default 1500 ft/min).
    #[schemars(with = "Option<Speed>")]
    pub climb_rate: Option<Spanned<Speed>>,
}

/// One maneuver primitive.
///
/// - `straight`: `duration`.
/// - `turn`: `rate` or `bank`; then `duration`, `heading` or `by`; `direction` is
///   required with `duration` and `by` (with `heading` the default is the shortest
///   way).
/// - `climb`: `vertical_rate` or `angle` (negative descends); then `altitude` or
///   `duration`.
/// - `speed`: `speed`.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ManeuverSpec {
    /// The primitive.
    pub kind: ManeuverKind,
    /// Total duration (straight; turn and climb by duration).
    #[schemars(with = "Option<Duration>")]
    pub duration: Option<Spanned<Duration>>,
    /// Turn rate.
    #[schemars(with = "Option<AngularRate>")]
    pub rate: Option<Spanned<AngularRate>>,
    /// Turn bank angle.
    #[schemars(with = "Option<Angle>")]
    pub bank: Option<Spanned<Angle>>,
    /// Turn direction.
    pub direction: Option<DirectionSpec>,
    /// Turn to this track angle.
    #[schemars(with = "Option<Angle>")]
    pub heading: Option<Spanned<Angle>>,
    /// Turn by this angle relative to a straight path.
    #[schemars(with = "Option<Angle>")]
    pub by: Option<Spanned<Angle>>,
    /// Climb or descent rate.
    #[schemars(with = "Option<Speed>")]
    pub vertical_rate: Option<Spanned<Speed>>,
    /// Climb or descent flight-path angle.
    #[schemars(with = "Option<Angle>")]
    pub angle: Option<Spanned<Angle>>,
    /// Level off at this altitude.
    #[schemars(with = "Option<Altitude>")]
    pub altitude: Option<Spanned<Altitude>>,
    /// Target speed.
    #[schemars(with = "Option<Speed>")]
    pub speed: Option<Spanned<Speed>>,
}

/// Maneuver primitives.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ManeuverKind {
    /// Straight and level.
    Straight,
    /// A coordinated, level turn.
    Turn,
    /// A straight climb or descent.
    Climb,
    /// A straight, level speed change.
    Speed,
}

/// Turn directions.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum DirectionSpec {
    /// Counter-clockwise seen from above.
    Left,
    /// Clockwise seen from above.
    Right,
}

/// A route waypoint.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct WaypointSpec {
    /// Latitude.
    pub lat: Latitude,
    /// Longitude.
    pub lon: Longitude,
    /// Altitude reached by the waypoint (AGL for ground platforms).
    #[schemars(with = "Altitude")]
    pub alt: Spanned<Altitude>,
    /// Speed reached by the waypoint (default: unchanged).
    #[schemars(with = "Option<Speed>")]
    pub speed: Option<Spanned<Speed>>,
}

/// A recorded track position.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TrackPointSpec {
    /// Truth time of the fix.
    #[schemars(with = "Duration")]
    pub time: Spanned<Duration>,
    /// Latitude.
    pub lat: Latitude,
    /// Longitude.
    pub lon: Longitude,
    /// Altitude.
    pub alt: Altitude,
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
    /// Waveform: the built-in `cw`, the built-in `iq_file` (with [`Self::iq_file`]), or
    /// a library reference (not resolved yet).
    pub waveform: LibraryRef,
    /// The recording an `iq_file` waveform plays (required exactly when
    /// `waveform = "iq_file"`).
    pub iq_file: Option<IqFileSpec>,
}

/// An IQ-file emitter: a SigMF ci16_le or cf32_le recording played on the emitter's
/// carrier (contracts §5.6). The emitter `frequency` is the RF frequency of the file's
/// 0 Hz, and its native rate is the recording's `core:sample_rate`.
#[derive(Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct IqFileSpec {
    /// The recording: its `.sigmf-meta` or `.sigmf-data` file, or the common stem,
    /// relative to the scenario file unless absolute.
    #[schemars(with = "String")]
    pub path: Spanned<String>,
    /// Truth time at which sample 0 is emitted (default 0 s, not negative).
    #[schemars(with = "Option<Duration>")]
    pub start: Option<Spanned<Duration>>,
    /// What happens after the last sample.
    pub end: EndPolicySpec,
    /// How sample values map to the emitter `power`.
    pub normalization: NormalizationSpec,
    /// The native rate, if stated: it must equal the recording's `core:sample_rate`.
    #[schemars(with = "Option<SampleRate>")]
    pub sample_rate: Option<Spanned<SampleRate>>,
    /// The source reconstruction response (defaults: bandwidth 0.8 × and stopband 1 ×
    /// the native rate, 80 dB).
    #[serde(default)]
    pub response: SourceResponseSpec,
}

/// End-of-file behavior of an IQ-file emitter.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum EndPolicySpec {
    /// The recording must cover the run: rendering past its last sample is an error.
    Stop,
    /// Silence after the last sample.
    ZeroPad,
    /// The recording repeats; its last sample is followed directly by its first.
    Loop,
}

/// Level normalization of an IQ-file emitter.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum NormalizationSpec {
    /// `power` is the power of a full-scale sample (|x| = 1: cf32 magnitude 1, ci16
    /// magnitude 32768).
    FullScale,
    /// `power` is the recording's mean power over all its samples.
    Rms,
}

/// The source reconstruction response of an IQ-file emitter (full widths).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct SourceResponseSpec {
    /// Passband full width (default 0.8 × the native rate).
    #[schemars(with = "Option<Frequency>")]
    pub bandwidth: Option<Spanned<Frequency>>,
    /// Stopband full width (default the native rate, at most the native rate).
    #[schemars(with = "Option<Frequency>")]
    pub stopband: Option<Spanned<Frequency>>,
    /// Stopband attenuation (default 80 dB, 20 to 200 dB).
    #[schemars(with = "Option<Gain>")]
    pub attenuation: Option<Spanned<Gain>>,
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
