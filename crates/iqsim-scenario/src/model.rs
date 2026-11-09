//! The scenario file's top-level types.
//!
//! Every table rejects unknown keys. Fields marked as placeholders hold only the
//! choices that exist today (for example `clock = { model = "ideal" }`); later
//! prompts add variants and fields without changing the surrounding layout.

use iqsim_units::{DataRate, Frequency, Position, Power, SampleRate};
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
    /// Usable bandwidth of the receiver response H_rx (at most the sample rate).
    #[schemars(with = "Frequency")]
    pub bandwidth: Spanned<Frequency>,
    /// Sampling clock model (placeholder: only `ideal`).
    #[serde(default)]
    pub clock: ClockModelSpec,
    /// Local-oscillator model (placeholder: only `ideal`).
    #[serde(default)]
    pub lo: LoModelSpec,
    /// Sinks this channel streams to.
    #[serde(default)]
    pub sinks: Vec<Ident>,
}

/// A sampling-clock model (placeholder).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ClockModelSpec {
    /// The clock model.
    pub model: ClockModel,
}

/// Available sampling-clock models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ClockModel {
    /// Samples exactly at the nominal instants (`IdealClockMap`).
    #[default]
    Ideal,
}

/// A local-oscillator model (placeholder).
#[derive(Debug, Default, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LoModelSpec {
    /// The LO model.
    pub model: LoModel,
}

/// Available local-oscillator models.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum LoModel {
    /// No frequency or phase error.
    #[default]
    Ideal,
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
    /// Bandwidth of the output response B_out (at most the sample rate).
    #[schemars(with = "Frequency")]
    pub bandwidth: Spanned<Frequency>,
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
    /// Shared resources this sink's traffic crosses.
    #[serde(default)]
    pub resources: Vec<Ident>,
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
