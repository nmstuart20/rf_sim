//! Scenario resolution: a parsed scenario becomes a [`ResolvedScenario`] (every default
//! filled in, serializable, hashed into the run manifest) and then a [`RunPlan`] of
//! runtime objects.
//!
//! What a v1 offline run supports, and how each scenario item maps:
//!
//! - Platforms are static. A position becomes WGS84 ECEF. An MSL height adds the
//!   geoid height N(φ, λ) from `[world] geoid` (N = 0 without a grid, with a manifest
//!   warning). A flight level is a pressure altitude: it goes through the inverse ISA
//!   chain with `[world] pressure_offset`, then the geoid (contracts §12.6). AGL
//!   heights are still taken as heights above the ellipsoid, because there is no
//!   terrain model yet (phase 8). Mounts are ignored: antennas are isotropic at the
//!   platform origin.
//! - A geoid grid is read once here; its file SHA-256, length and header values go
//!   into the resolved scenario (so a changed grid changes the resolved hash) and the
//!   manifest registers it as a `geoid_grid` asset.
//! - Emitters must use `behavior = "always_on"` and the built-in `waveform = "cw"` or
//!   `waveform = "iq_file"`. Library resolution comes later.
//! - An `iq_file` path is relative to the scenario file. Resolution opens the recording
//!   (one streaming pass: validation, hashes, peak and mean power) and records its
//!   native rate (the shortest decimal of `core:sample_rate`, parsed exactly), length,
//!   start, end policy, normalization factor, reconstruction response and the SHA-256
//!   of both files, so a changed recording changes the resolved hash.
//! - Each receiver channel is one virtual-receiver stream (`<receiver>.<channel>`), and
//!   each antenna port one antenna-port stream (`<port>`). Only SigMF sinks are written.
//! - A response's passband edge is `bandwidth / 2` and its stopband edge `stopband / 2`
//!   (default: the sample rate), both floored to whole µHz.
//! - The virtual target's accumulator full scale is the ADC full scale referred to the
//!   receiver input (`V_adc / g`), so the accumulator and the ADC use the same headroom.
//! - Receiver entity ids are the first four bytes (little-endian) of
//!   SHA-256(`"iqsim.entity.v1\0"` ‖ name), so they do not depend on file order.

use std::fmt;
use std::path::{Path, PathBuf};

use iqsim_core::amplitude::{LOSSLESS_CULL_LSB, rms_volts};
use iqsim_core::{
    ChannelTimebase, IdealClockMap, MAX_BLOCK_LEN, ResponseSpec, SampleRate, TapKind, TapResponse,
    TargetScale, TruthTime,
};
use iqsim_geo::{AltitudeModel, Atmosphere, Ecef, Geoid, GeoidGrid};
use iqsim_io::iq_input::IqInput;
use iqsim_receiver::{
    AdcConfig, AgcConfig, ChainConfig, DecimationConfig, FrontEnd, IqImpairment, LoConfig,
    NoiseConfig, OutputFormat, PhaseNoise,
};
use iqsim_scenario::{
    ClockModel, EndPolicySpec, IqFileSpec, LoModel, NormalizationSpec, ResponseDetail,
    SampleFormatSpec, Scenario, SinkKind,
};
use iqsim_units::{AltitudeRef, Frequency, Position};
use iqsim_waveform::continuous::source_response;
use iqsim_waveform::{EndPolicy, SourceTimeline};
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::static_scene::{
    ContinuousSpec, SharedSource, StaticEmitter, StaticScene, StaticTarget, StaticWaveform,
};

/// Pages an IQ-file reader may cache (64 MiB at 1 MiB per page).
pub const IQ_CACHE_PAGES: usize = 64;

/// Default block length, ps (1 ms).
pub const DEFAULT_BLOCK_PS: i64 = 1_000_000_000;

/// Default stopband attenuation, milli-dB (80 dB).
pub const DEFAULT_ATTENUATION_MDB: u32 = 80_000;

/// A scenario that cannot be run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolveError(pub String);

impl fmt::Display for ResolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ResolveError {}

fn err<T>(msg: impl Into<String>) -> Result<T, ResolveError> {
    Err(ResolveError(msg.into()))
}

/// A pre-sampling response definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct ResolvedResponse {
    /// Passband edge, µHz.
    pub passband_edge_uhz: i64,
    /// Stopband edge, µHz.
    pub stopband_edge_uhz: i64,
    /// Stopband attenuation, milli-dB.
    pub stopband_atten_mdb: u32,
    /// Group delay, ps.
    pub group_delay_ps: i64,
}

impl ResolvedResponse {
    fn spec(self) -> ResponseSpec {
        ResponseSpec {
            passband_edge_uhz: self.passband_edge_uhz,
            stopband_edge_uhz: self.stopband_edge_uhz,
            stopband_atten_mdb: self.stopband_atten_mdb,
            group_delay_ps: self.group_delay_ps,
        }
    }
}

/// Resolved front end.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct ResolvedFrontEnd {
    /// Gain, dB.
    pub gain_db: f64,
    /// Noise figure, dB.
    pub noise_figure_db: f64,
}

/// Resolved IQ impairment.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct ResolvedIq {
    /// Q gain relative to I, dB.
    pub gain_imbalance_db: f64,
    /// Quadrature skew, degrees.
    pub phase_skew_deg: f64,
    /// DC on I, volts at the ADC input.
    pub dc_i_volts: f64,
    /// DC on Q, volts at the ADC input.
    pub dc_q_volts: f64,
}

/// Resolved decimation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct ResolvedDecimation {
    /// Factor.
    pub factor: u32,
    /// Filter at the input rate (group delay 0).
    pub filter: ResolvedResponse,
}

/// Resolved AGC.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct ResolvedAgc {
    /// Target RMS, volts at the ADC input.
    pub target_volts: f64,
    /// Attack coefficient.
    pub attack: f64,
    /// Gain step.
    pub step: f64,
    /// Smallest gain, dB.
    pub min_gain_db: f64,
    /// Largest gain, dB.
    pub max_gain_db: f64,
}

/// Resolved receiver chain of one channel.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
pub struct ResolvedChain {
    /// Stage 1.
    pub front_end: Option<ResolvedFrontEnd>,
    /// Stage 2 on.
    pub thermal_noise: bool,
    /// Stage 3: LO offset, µHz.
    pub lo_offset_uhz: Option<i64>,
    /// Stage 4.
    pub iq: Option<ResolvedIq>,
    /// Stage 5.
    pub decimation: Option<ResolvedDecimation>,
    /// Stage 6.
    pub agc: Option<ResolvedAgc>,
    /// ADC full scale, volts per component.
    pub adc_full_scale_volts: f64,
    /// Stage 7: ADC bits when quantizing.
    pub adc_bits: Option<u32>,
}

/// One resolved emitter.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResolvedEmitter {
    /// Name.
    pub name: String,
    /// ECEF position, m.
    pub position_ecef_m: [f64; 3],
    /// Carrier, µHz.
    pub carrier_uhz: i64,
    /// Power, W.
    pub power_w: f64,
    /// Built-in waveform: `cw` or `iq_file`.
    pub waveform: String,
    /// The recording of an `iq_file` emitter.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<ResolvedIqSource>,
}

/// A resolved IQ-file source.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResolvedIqSource {
    /// Absolute path of the `.sigmf-meta` file.
    pub meta_path: String,
    /// Absolute path of the `.sigmf-data` file.
    pub data_path: String,
    /// `ci16_le` or `cf32_le`.
    pub datatype: String,
    /// Native rate numerator.
    pub rate_num: u64,
    /// Native rate denominator.
    pub rate_den: u64,
    /// Complex samples N.
    pub sample_count: u64,
    /// Emission time of sample 0, ps.
    pub start_ps: i64,
    /// `stop`, `zero_pad` or `loop`.
    pub end: String,
    /// `full_scale` or `rms`.
    pub normalization: String,
    /// Factor applied to the normalized samples (1, or 1/RMS).
    pub scale: f64,
    /// max |x| of the normalized samples.
    pub peak: f64,
    /// Source reconstruction response φ.
    pub response: ResolvedResponse,
    /// Data file length, bytes.
    pub data_byte_len: u64,
    /// SHA-256 of the data file.
    pub data_sha256: String,
    /// Metadata file length, bytes.
    pub meta_byte_len: u64,
    /// SHA-256 of the metadata file.
    pub meta_sha256: String,
}

/// The `.sigmf-meta` and `.sigmf-data` paths for a scenario `iq_file.path`.
fn recording_paths(base: &Path, text: &str) -> (PathBuf, PathBuf) {
    let stem = text
        .strip_suffix(".sigmf-meta")
        .or_else(|| text.strip_suffix(".sigmf-data"))
        .unwrap_or(text);
    let stem = base.join(stem);
    let with = |ext: &str| {
        let mut p = stem.clone().into_os_string();
        p.push(ext);
        PathBuf::from(p)
    };
    (with(".sigmf-meta"), with(".sigmf-data"))
}

/// The exact rate a SigMF `core:sample_rate` names: its shortest round-trip decimal,
/// parsed as an exact quantity.
fn exact_rate(hz: f64, owner: &str) -> Result<SampleRate, ResolveError> {
    let text = format!("{hz} S/s");
    text.parse::<iqsim_units::SampleRate>()
        .map(|r| r.get())
        .map_err(|e| ResolveError(format!("{owner}: core:sample_rate {hz}: {e}")))
}

fn resolve_iq(q: &IqFileSpec, base: &Path, owner: &str) -> Result<ResolvedIqSource, ResolveError> {
    let (meta, data) = recording_paths(base, q.path.get_ref());
    let canon = |p: &Path| {
        std::fs::canonicalize(p).map_err(|e| ResolveError(format!("{owner}: {}: {e}", p.display())))
    };
    let (meta, data) = (canon(&meta)?, canon(&data)?);
    let input = IqInput::open(&meta, &data).map_err(|e| ResolveError(format!("{owner}: {e}")))?;
    let m = input.meta();
    let sum = input.summary();
    let rate = exact_rate(m.sample_rate_hz, owner)?;
    if let Some(stated) = &q.sample_rate
        && stated.get_ref().get() != rate
    {
        return err(format!(
            "{owner}: iq_file.sample_rate {} differs from the recording's {}",
            stated.get_ref(),
            iqsim_units::SampleRate::new(rate)
        ));
    }
    let rate_uhz = rate_uhz(rate);
    let r = &q.response;
    let bw = r
        .bandwidth
        .as_ref()
        .map_or(rate_uhz / 5 * 4, |b| b.get_ref().as_micro_hz());
    let stop = r
        .stopband
        .as_ref()
        .map_or(rate_uhz, |b| b.get_ref().as_micro_hz());
    let atten = match &r.attenuation {
        Some(a) => {
            let mdb = libm::round(a.get_ref().db() * 1000.0);
            if !(20_000.0..=200_000.0).contains(&mdb) {
                return err(format!("{owner}: attenuation must be 20 to 200 dB"));
            }
            mdb as u32
        }
        None => DEFAULT_ATTENUATION_MDB,
    };
    let response = ResolvedResponse {
        passband_edge_uhz: bw / 2,
        stopband_edge_uhz: stop / 2,
        stopband_atten_mdb: atten,
        group_delay_ps: 0,
    };
    source_response(rate, response.spec())
        .map_err(|e| ResolveError(format!("{owner}: iq_file response: {e}")))?;
    let scale = match q.normalization {
        NormalizationSpec::FullScale => 1.0,
        NormalizationSpec::Rms => {
            let mean = sum.sum_power / sum.sample_count as f64;
            if mean.is_nan() || mean <= 0.0 {
                return err(format!(
                    "{owner}: an all-zero recording cannot be RMS-normalized"
                ));
            }
            1.0 / mean.sqrt()
        }
    };
    Ok(ResolvedIqSource {
        meta_path: meta.display().to_string(),
        data_path: data.display().to_string(),
        datatype: m.datatype.name().to_owned(),
        rate_num: rate.num(),
        rate_den: rate.den(),
        sample_count: sum.sample_count,
        start_ps: q.start.as_ref().map_or(0, |s| s.get_ref().get().as_ps()),
        end: match q.end {
            EndPolicySpec::Stop => "stop",
            EndPolicySpec::ZeroPad => "zero_pad",
            EndPolicySpec::Loop => "loop",
        }
        .to_owned(),
        normalization: match q.normalization {
            NormalizationSpec::FullScale => "full_scale",
            NormalizationSpec::Rms => "rms",
        }
        .to_owned(),
        scale,
        peak: sum.peak_magnitude,
        response,
        data_byte_len: sum.byte_len,
        data_sha256: hex::encode(sum.sha256),
        meta_byte_len: m.meta_byte_len,
        meta_sha256: hex::encode(m.meta_sha256),
    })
}

/// One resolved output stream (a tap of a channel).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResolvedStream {
    /// Stream id.
    pub id: String,
    /// `virtual_receiver` or `antenna_port`.
    pub tap: String,
    /// Receiver or antenna-port name.
    pub entity: String,
    /// Receiver entity id (noise key), 0 for antenna ports.
    pub entity_id: u32,
    /// Channel index within the receiver.
    pub channel_index: u16,
    /// ECEF position of the antenna, m.
    pub position_ecef_m: [f64; 3],
    /// Centre, µHz.
    pub center_uhz: i64,
    /// Render rate numerator.
    pub rate_num: u64,
    /// Render rate denominator.
    pub rate_den: u64,
    /// Block length, samples at the render rate.
    pub block_len: u64,
    /// Samples at the render rate.
    pub sample_count: u64,
    /// H_rx or B_out.
    pub response: ResolvedResponse,
    /// Accumulator full scale, volts.
    pub full_scale_volts: f64,
    /// Sampling-instant offset, ps (virtual receiver only).
    pub clock_offset_ps: i64,
    /// Receiver chain (virtual receiver only).
    pub chain: Option<ResolvedChain>,
    /// `ci16_le` or `cf32_le`.
    pub format: String,
    /// SigMF sinks.
    pub sinks: Vec<String>,
}

/// The resolved geoid model.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "model", rename_all = "snake_case")]
pub enum ResolvedGeoid {
    /// No grid configured: N = 0 (MSL is the ellipsoid). The manifest carries
    /// [`ZERO_GEOID_WARNING`].
    Zero,
    /// A GeographicLib PGM grid.
    Grid {
        /// Absolute path it was read from.
        path: String,
        /// File length, bytes.
        byte_len: u64,
        /// File SHA-256, hex.
        sha256: String,
        /// Columns.
        width: usize,
        /// Rows.
        height: usize,
        /// Header `Offset`, m.
        offset_m: f64,
        /// Header `Scale`, m per count.
        scale_m: f64,
        /// Interpolation (`bilinear`).
        interpolation: String,
    },
}

/// The resolved Earth model.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResolvedWorld {
    /// Geoid.
    pub geoid: ResolvedGeoid,
    /// Sea-level pressure offset from 1013.25 hPa, Pa.
    pub pressure_offset_pa: f64,
}

/// The run-manifest warning when no geoid grid is configured.
pub const ZERO_GEOID_WARNING: &str =
    "no geoid grid configured: geoid height is 0 everywhere (MSL heights are ellipsoid heights)";

/// Warnings a run with this resolved scenario records in its manifest.
pub fn resolution_warnings(resolved: &ResolvedScenario) -> Vec<String> {
    let mut w = Vec::new();
    if resolved.world.geoid == ResolvedGeoid::Zero {
        w.push(ZERO_GEOID_WARNING.to_owned());
    }
    w
}

/// A scenario with every default filled in. Its canonical JSON hash identifies the run
/// configuration in the manifest.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ResolvedScenario {
    /// Resolution rules version.
    pub resolver_version: u32,
    /// Scenario name.
    pub name: String,
    /// Epoch, POSIX seconds.
    pub epoch_unix_seconds: i64,
    /// Epoch nanoseconds.
    pub epoch_nanos: u32,
    /// Seed.
    pub seed: u64,
    /// Rendered duration, ps.
    pub duration_ps: i64,
    /// Block length, ps.
    pub block_ps: i64,
    /// Earth model.
    pub world: ResolvedWorld,
    /// Emitters, in file order.
    pub emitters: Vec<ResolvedEmitter>,
    /// Streams: receiver channels in file order, then antenna ports.
    pub streams: Vec<ResolvedStream>,
}

/// Receiver entity id from its name (module docs).
pub fn entity_id(name: &str) -> u32 {
    let mut h = Sha256::new();
    h.update(b"iqsim.entity.v1\0");
    h.update(name.as_bytes());
    let d = h.finalize();
    u32::from_le_bytes([d[0], d[1], d[2], d[3]])
}

fn ecef(position: &Position, model: &AltitudeModel, what: &str) -> Result<[f64; 3], ResolveError> {
    let alt = position.altitude();
    let (lat, lon) = (
        position.latitude().degrees(),
        position.longitude().degrees(),
    );
    let (lat_rad, lon_rad) = (lat.to_radians(), lon.to_radians());
    let geo = |e: iqsim_geo::GeoError| ResolveError(format!("{what}: {e}"));
    let h = match alt.reference() {
        AltitudeRef::Ellipsoid | AltitudeRef::Agl => alt.meters(),
        AltitudeRef::Msl => model
            .ellipsoid_height_from_msl_m(lat_rad, lon_rad, alt.meters())
            .map_err(geo)?,
        AltitudeRef::PressureStd => model
            .ellipsoid_height_from_pressure_altitude_m(lat_rad, lon_rad, alt.meters())
            .map_err(geo)?,
    };
    let e = Ecef::from_geodetic(lat, lon, h);
    Ok([e.x, e.y, e.z])
}

/// `num·10⁶/den` floored: a rate in µHz.
fn rate_uhz(rate: SampleRate) -> i64 {
    (i128::from(rate.num()) * 1_000_000 / i128::from(rate.den())) as i64
}

/// Samples in `ps` at `rate`, exactly, or an error.
fn exact_samples(ps: i64, rate: SampleRate, what: &str) -> Result<u64, ResolveError> {
    let n = i128::from(ps) * i128::from(rate.num());
    let d = i128::from(rate.den()) * 1_000_000_000_000;
    if n % d != 0 {
        return err(format!(
            "{what}: the block is not a whole number of samples"
        ));
    }
    u64::try_from(n / d).map_err(|_| ResolveError(format!("{what}: block too long")))
}

fn response(
    bandwidth: Frequency,
    rate: SampleRate,
    detail: &ResponseDetail,
) -> Result<ResolvedResponse, ResolveError> {
    let stop = detail
        .stopband
        .as_ref()
        .map_or(rate_uhz(rate), |s| s.get_ref().as_micro_hz());
    let atten = match &detail.attenuation {
        Some(a) => {
            let mdb = libm::round(a.get_ref().db() * 1000.0);
            if !(20_000.0..=200_000.0).contains(&mdb) {
                return err("stopband attenuation must be 20 to 200 dB");
            }
            mdb as u32
        }
        None => DEFAULT_ATTENUATION_MDB,
    };
    Ok(ResolvedResponse {
        passband_edge_uhz: bandwidth.as_micro_hz() / 2,
        stopband_edge_uhz: stop / 2,
        stopband_atten_mdb: atten,
        group_delay_ps: detail
            .group_delay
            .as_ref()
            .map_or(0, |d| d.get_ref().get().as_ps()),
    })
}

fn format_of(
    scenario: &Scenario,
    sinks: &[iqsim_scenario::Ident],
    owner: &str,
) -> Result<(String, Vec<String>), ResolveError> {
    let mut format: Option<SampleFormatSpec> = None;
    let mut names = Vec::new();
    for s in sinks {
        let Some(sink) = scenario
            .sinks
            .iter()
            .find(|k| k.name.as_str() == s.as_str())
        else {
            return err(format!("{owner}: unknown sink \"{}\"", s.as_str()));
        };
        if sink.kind != SinkKind::Sigmf {
            return err(format!(
                "{owner}: sink \"{}\" is {:?}; only SigMF sinks are implemented so far",
                s.as_str(),
                sink.kind
            ));
        }
        if format.is_some_and(|f| f != sink.format) {
            return err(format!(
                "{owner}: its sinks ask for different sample formats"
            ));
        }
        format = Some(sink.format);
        names.push(s.as_str().to_owned());
    }
    let f = match format.unwrap_or_default() {
        SampleFormatSpec::Ci16Le => "ci16_le",
        SampleFormatSpec::Cf32Le => "cf32_le",
    };
    Ok((f.to_owned(), names))
}

fn platform_position(
    scenario: &Scenario,
    model: &AltitudeModel,
    name: &str,
    owner: &str,
) -> Result<[f64; 3], ResolveError> {
    let p = scenario
        .platforms
        .iter()
        .find(|p| p.name.as_str() == name)
        .ok_or_else(|| ResolveError(format!("{owner}: unknown platform \"{name}\"")))?;
    ecef(&p.position, model, owner)
}

/// Resolves a parsed, lint-clean scenario. Relative input paths are taken from
/// `base`, the scenario file's directory.
pub fn resolve(scenario: &Scenario, base: &Path) -> Result<ResolvedScenario, ResolveError> {
    let Some(render) = &scenario.render else {
        return err("`sim run` needs a [render] table with a duration");
    };
    let duration_ps = render.duration.get_ref().get().as_ps();
    let block_ps = render
        .block
        .as_ref()
        .map_or(DEFAULT_BLOCK_PS, |b| b.get_ref().get().as_ps());
    if duration_ps <= 0 || block_ps <= 0 {
        return err("render duration and block must be positive");
    }
    let (world, model) = resolve_world(scenario, base)?;

    let mut emitters = Vec::new();
    for e in &scenario.emitters {
        let owner = format!("emitter \"{}\"", e.name.as_str());
        let waveform = e.waveform.as_str();
        if !matches!(waveform, "cw" | "iq_file") || e.behavior.as_str() != "always_on" {
            return err(format!(
                "{owner}: only the built-in waveforms \"cw\" and \"iq_file\" with behavior \"always_on\" are available so far"
            ));
        }
        let source = match (&e.iq_file, waveform) {
            (Some(q), "iq_file") => Some(resolve_iq(q, base, &owner)?),
            (None, "cw") => None,
            _ => {
                return err(format!(
                    "{owner}: iq_file must be given exactly for waveform \"iq_file\""
                ));
            }
        };
        emitters.push(ResolvedEmitter {
            name: e.name.as_str().to_owned(),
            position_ecef_m: platform_position(scenario, &model, e.platform.as_str(), &owner)?,
            carrier_uhz: e.frequency.get_ref().as_micro_hz(),
            power_w: e.power.watts(),
            waveform: waveform.to_owned(),
            source,
        });
    }

    let mut streams = Vec::new();
    let mut ids: Vec<(u32, String)> = Vec::new();
    for r in &scenario.receivers {
        let id = entity_id(r.name.as_str());
        if let Some((_, other)) = ids.iter().find(|(i, _)| *i == id) {
            return err(format!(
                "receivers \"{other}\" and \"{}\" hash to the same entity id; rename one",
                r.name.as_str()
            ));
        }
        ids.push((id, r.name.as_str().to_owned()));
        let position = platform_position(
            scenario,
            &model,
            r.platform.as_str(),
            &format!("receiver \"{}\"", r.name.as_str()),
        )?;
        for (ci, c) in r.channels.iter().enumerate() {
            let sid = format!("{}.{}", r.name.as_str(), c.name.as_str());
            let owner = format!("channel \"{sid}\"");
            let rate = c.sample_rate.get_ref().get();
            let resp = response(*c.bandwidth.get_ref(), rate, &c.response)?;
            let front_end = c.front_end.as_ref().map(|f| ResolvedFrontEnd {
                gain_db: f.gain.db(),
                noise_figure_db: f.noise_figure.map_or(0.0, |n| n.db()),
            });
            let adc_fs_w = c
                .adc
                .as_ref()
                .and_then(|a| a.full_scale)
                .map_or(1e-3, |p| p.watts());
            let v_adc = rms_volts(adc_fs_w);
            let quantize = c.adc.as_ref().and_then(|a| a.enabled).unwrap_or(true);
            let adc_bits = quantize.then(|| c.adc.as_ref().and_then(|a| a.bits).unwrap_or(16));
            let decimation = match &c.decimation {
                Some(d) => {
                    if d.factor < 2 {
                        return err(format!("{owner}: decimation factor must be at least 2"));
                    }
                    let out_uhz = rate_uhz(rate) / i64::from(d.factor);
                    let bw = d.bandwidth.map_or(out_uhz * 4 / 5, |b| b.as_micro_hz());
                    let stop = d.stopband.map_or(out_uhz, |b| b.as_micro_hz());
                    let atten = d.attenuation.map_or(DEFAULT_ATTENUATION_MDB, |a| {
                        libm::round(a.db() * 1000.0) as u32
                    });
                    Some(ResolvedDecimation {
                        factor: d.factor,
                        filter: ResolvedResponse {
                            passband_edge_uhz: bw / 2,
                            stopband_edge_uhz: stop / 2,
                            stopband_atten_mdb: atten,
                            group_delay_ps: 0,
                        },
                    })
                }
                None => None,
            };
            let chain = ResolvedChain {
                front_end,
                thermal_noise: c.thermal_noise.unwrap_or(true),
                lo_offset_uhz: match c.lo.model {
                    LoModel::Ideal => None,
                    LoModel::FixedOffset => Some(
                        c.lo.offset
                            .as_ref()
                            .ok_or_else(|| ResolveError(format!("{owner}: LO offset missing")))?
                            .get_ref()
                            .as_micro_hz(),
                    ),
                },
                iq: c.iq.as_ref().map(|q| ResolvedIq {
                    gain_imbalance_db: q.gain_imbalance.map_or(0.0, |g| g.db()),
                    phase_skew_deg: q.phase_skew.map_or(0.0, |a| a.degrees()),
                    dc_i_volts: q.dc_i * v_adc,
                    dc_q_volts: q.dc_q * v_adc,
                }),
                decimation,
                agc: c.agc.as_ref().map(|a| ResolvedAgc {
                    target_volts: v_adc * libm::exp10(a.target.map_or(-12.0, |g| g.db()) / 20.0),
                    attack: a.attack.unwrap_or(0.01),
                    step: a.step.unwrap_or(0.01),
                    min_gain_db: a.min_gain.map_or(-40.0, |g| g.db()),
                    max_gain_db: a.max_gain.map_or(40.0, |g| g.db()),
                }),
                adc_full_scale_volts: v_adc,
                adc_bits,
            };
            let g = front_end.map_or(1.0, |f| libm::exp10(f.gain_db / 20.0));
            let (format, sinks) = format_of(scenario, &c.sinks, &owner)?;
            if format == "ci16_le" && adc_bits.is_none() {
                return err(format!("{owner}: ci16_le output needs the ADC enabled"));
            }
            let clock_offset_ps = match c.clock.model {
                ClockModel::Ideal => 0,
                ClockModel::FixedOffset => c
                    .clock
                    .offset
                    .as_ref()
                    .ok_or_else(|| ResolveError(format!("{owner}: clock offset missing")))?
                    .get_ref()
                    .get()
                    .as_ps(),
            };
            streams.push(stream(
                sid,
                "virtual_receiver",
                r.name.as_str(),
                id,
                ci,
                position,
                c.center_frequency.get_ref().as_micro_hz(),
                rate,
                block_ps,
                duration_ps,
                resp,
                v_adc / g,
                clock_offset_ps,
                Some(chain),
                format,
                sinks,
            )?);
        }
    }
    for a in &scenario.antenna_ports {
        let owner = format!("antenna port \"{}\"", a.name.as_str());
        let rate = a.sample_rate.get_ref().get();
        let resp = response(*a.bandwidth.get_ref(), rate, &a.response)?;
        let (format, sinks) = format_of(scenario, &a.sinks, &owner)?;
        let v_fs = rms_volts(a.full_scale.map_or(1e-3, |p| p.watts()));
        streams.push(stream(
            a.name.as_str().to_owned(),
            "antenna_port",
            a.name.as_str(),
            0,
            0,
            platform_position(scenario, &model, a.platform.as_str(), &owner)?,
            a.center_frequency.get_ref().as_micro_hz(),
            rate,
            block_ps,
            duration_ps,
            resp,
            v_fs,
            0,
            None,
            format,
            sinks,
        )?);
    }
    Ok(ResolvedScenario {
        resolver_version: RESOLVER_VERSION,
        name: scenario.metadata.name.get_ref().clone(),
        epoch_unix_seconds: scenario.epoch.unix_seconds(),
        epoch_nanos: scenario.epoch.nanosecond(),
        seed: scenario.seed,
        duration_ps,
        block_ps,
        world,
        emitters,
        streams,
    })
}

/// Resolution rules version (2 since prompt 7: the Earth model).
pub const RESOLVER_VERSION: u32 = 2;

/// Reads the geoid grid (if any) and builds the altitude model.
fn resolve_world(
    scenario: &Scenario,
    base: &Path,
) -> Result<(ResolvedWorld, AltitudeModel), ResolveError> {
    let pressure_offset_pa = scenario
        .world
        .pressure_offset
        .as_ref()
        .map_or(0.0, |p| p.get_ref().pascals());
    let atmosphere = Atmosphere::new(pressure_offset_pa)
        .map_err(|e| ResolveError(format!("world.pressure_offset: {e}")))?;
    let (geoid, resolved) = match &scenario.world.geoid {
        None => (Geoid::Zero, ResolvedGeoid::Zero),
        Some(g) => {
            let rel = Path::new(g.path.get_ref());
            let path = if rel.is_absolute() {
                rel.to_path_buf()
            } else {
                base.join(rel)
            };
            let bytes = std::fs::read(&path)
                .map_err(|e| ResolveError(format!("world.geoid: {}: {e}", path.display())))?;
            let grid = GeoidGrid::from_pgm_bytes(&bytes)
                .map_err(|e| ResolveError(format!("world.geoid: {}: {e}", path.display())))?;
            let m = grid.metadata();
            let resolved = ResolvedGeoid::Grid {
                path: std::path::absolute(&path)
                    .unwrap_or(path)
                    .display()
                    .to_string(),
                byte_len: bytes.len() as u64,
                sha256: hex::encode(Sha256::digest(&bytes)),
                width: m.width,
                height: m.height,
                offset_m: m.offset_m,
                scale_m: m.scale_m,
                interpolation: m.interpolation.to_owned(),
            };
            (Geoid::Grid(grid), resolved)
        }
    };
    Ok((
        ResolvedWorld {
            geoid: resolved,
            pressure_offset_pa,
        },
        AltitudeModel::new(geoid, atmosphere),
    ))
}

#[allow(clippy::too_many_arguments)]
fn stream(
    id: String,
    tap: &str,
    entity: &str,
    entity_id: u32,
    channel_index: usize,
    position_ecef_m: [f64; 3],
    center_uhz: i64,
    rate: SampleRate,
    block_ps: i64,
    duration_ps: i64,
    response: ResolvedResponse,
    full_scale_volts: f64,
    clock_offset_ps: i64,
    chain: Option<ResolvedChain>,
    format: String,
    sinks: Vec<String>,
) -> Result<ResolvedStream, ResolveError> {
    let owner = format!("stream \"{id}\"");
    let block_len = exact_samples(block_ps, rate, &owner)?;
    if block_len == 0 || block_len > MAX_BLOCK_LEN {
        return err(format!(
            "{owner}: block length {block_len} samples is out of range"
        ));
    }
    let tb = ChannelTimebase::new(TruthTime::from_ps(0), rate);
    let sample_count = tb
        .index_at_or_after(TruthTime::from_ps(duration_ps))
        .map_err(|e| ResolveError(format!("{owner}: {e}")))?
        .get();
    Ok(ResolvedStream {
        id,
        tap: tap.to_owned(),
        entity: entity.to_owned(),
        entity_id,
        channel_index: u16::try_from(channel_index)
            .map_err(|_| ResolveError(format!("{owner}: too many channels")))?,
        position_ecef_m,
        center_uhz,
        rate_num: rate.num(),
        rate_den: rate.den(),
        block_len,
        sample_count,
        response,
        full_scale_volts,
        clock_offset_ps,
        chain,
        format,
        sinks,
    })
}

/// Runtime objects for one stream.
#[derive(Clone, Debug)]
pub struct StreamPlan {
    /// The resolved stream.
    pub resolved: ResolvedStream,
    /// The scene seen from this stream's antenna.
    pub scene: StaticScene,
    /// The render target.
    pub target: StaticTarget,
    /// The stream's ClockMap.
    pub clock: IdealClockMap,
    /// The receiver chain (virtual receiver only).
    pub chain: Option<ChainConfig>,
    /// Output format.
    pub format: OutputFormat,
}

impl StreamPlan {
    /// The render rate.
    pub fn rate(&self) -> SampleRate {
        self.target.rate
    }
}

/// Everything an offline run needs.
#[derive(Clone, Debug)]
pub struct RunPlan {
    /// The resolved scenario.
    pub resolved: ResolvedScenario,
    /// One plan per stream.
    pub streams: Vec<StreamPlan>,
}

impl RunPlan {
    /// Builds runtime objects from a resolved scenario.
    pub fn new(resolved: ResolvedScenario) -> Result<Self, ResolveError> {
        let mut emitters: Vec<StaticEmitter> = Vec::new();
        for e in &resolved.emitters {
            let waveform = match &e.source {
                None => StaticWaveform::Cw,
                Some(src) => StaticWaveform::Continuous(open_source(src, &e.name)?),
            };
            emitters.push(StaticEmitter {
                position: Ecef::new(
                    e.position_ecef_m[0],
                    e.position_ecef_m[1],
                    e.position_ecef_m[2],
                ),
                carrier_uhz: e.carrier_uhz,
                power_w: e.power_w,
                waveform,
            });
        }
        let mut streams = Vec::new();
        for s in &resolved.streams {
            let owner = format!("stream \"{}\"", s.id);
            let rate = SampleRate::from_ratio(s.rate_num, s.rate_den)
                .map_err(|e| ResolveError(format!("{owner}: {e}")))?;
            let kind = if s.chain.is_some() {
                TapKind::VirtualReceiver
            } else {
                TapKind::AntennaPort
            };
            let target = StaticTarget {
                kind,
                center_uhz: s.center_uhz,
                rate,
                response: s.response.spec(),
                full_scale_volts: s.full_scale_volts,
                cull_threshold_lsb: LOSSLESS_CULL_LSB,
            };
            let clock = IdealClockMap::new(ChannelTimebase::new(
                TruthTime::from_ps(s.clock_offset_ps),
                rate,
            ));
            let format = if s.format == "ci16_le" {
                OutputFormat::Ci16
            } else {
                OutputFormat::Cf32
            };
            let chain = match &s.chain {
                Some(c) => Some(
                    chain_config(&resolved, s, c, &target, format)
                        .map_err(|e| ResolveError(format!("{owner}: {e}")))?,
                ),
                None => None,
            };
            let p = s.position_ecef_m;
            streams.push(StreamPlan {
                resolved: s.clone(),
                scene: StaticScene {
                    receiver: Ecef::new(p[0], p[1], p[2]),
                    emitters: emitters.clone(),
                },
                target,
                clock,
                chain,
                format,
            });
        }
        Ok(Self { resolved, streams })
    }
}

/// The [`iqsim_waveform::ContinuousSource`] view of an IQ-file reader.
struct IqFileSource {
    len: u64,
    reader: iqsim_io::iq_input::IqFileReader,
}

impl iqsim_waveform::ContinuousSource for IqFileSource {
    fn len(&self) -> u64 {
        self.len
    }
    fn fill(&mut self, first: u64, out: &mut [[f64; 2]]) -> Result<(), String> {
        self.reader.read(first, out).map_err(|e| e.to_string())
    }
}

/// Opens a resolved IQ-file source for rendering. The recording must still be the one
/// resolution hashed.
fn open_source(src: &ResolvedIqSource, name: &str) -> Result<ContinuousSpec, ResolveError> {
    let owner = format!("emitter \"{name}\"");
    let input = IqInput::open(Path::new(&src.meta_path), Path::new(&src.data_path))
        .map_err(|e| ResolveError(format!("{owner}: {e}")))?;
    if hex::encode(input.summary().sha256) != src.data_sha256
        || hex::encode(input.meta().meta_sha256) != src.meta_sha256
    {
        return err(format!(
            "{owner}: the recording changed since it was resolved"
        ));
    }
    let rate = SampleRate::from_ratio(src.rate_num, src.rate_den)
        .map_err(|e| ResolveError(format!("{owner}: {e}")))?;
    let policy = match src.end.as_str() {
        "stop" => EndPolicy::Stop,
        "zero_pad" => EndPolicy::ZeroPad,
        "loop" => EndPolicy::Loop,
        other => return err(format!("{owner}: unknown end policy {other}")),
    };
    Ok(ContinuousSpec {
        source: SharedSource::new(IqFileSource {
            len: src.sample_count,
            reader: input.reader(IQ_CACHE_PAGES),
        }),
        timeline: SourceTimeline {
            len: src.sample_count,
            policy,
        },
        rate,
        start: TruthTime::from_ps(src.start_ps),
        response: src.response.spec(),
        scale: src.scale,
        peak: src.peak,
    })
}

fn chain_config(
    resolved: &ResolvedScenario,
    s: &ResolvedStream,
    c: &ResolvedChain,
    target: &StaticTarget,
    output: OutputFormat,
) -> Result<ChainConfig, String> {
    let response = TapResponse::new(target.kind, target.response, target.rate, target.center_uhz)
        .map_err(|e| e.to_string())?;
    let input_scale = TargetScale::new(target.full_scale_volts).map_err(|e| e.to_string())?;
    Ok(ChainConfig {
        rate: target.rate,
        response,
        input_scale,
        adc_full_scale_volts: c.adc_full_scale_volts,
        output,
        front_end: c.front_end.map(|f| FrontEnd {
            gain_db: f.gain_db,
            noise_figure_db: f.noise_figure_db,
        }),
        noise: c.thermal_noise.then_some(NoiseConfig {
            seed: resolved.seed,
            entity: s.entity_id,
            channel: s.channel_index,
        }),
        lo: c.lo_offset_uhz.map(|o| LoConfig {
            offset_uhz: o,
            phase_noise: PhaseNoise::None,
        }),
        iq: c.iq.map(|q| IqImpairment {
            gain_imbalance_db: q.gain_imbalance_db,
            phase_skew_deg: q.phase_skew_deg,
            dc_i_volts: q.dc_i_volts,
            dc_q_volts: q.dc_q_volts,
        }),
        decimation: c.decimation.map(|d| DecimationConfig {
            factor: d.factor,
            filter: d.filter.spec(),
        }),
        agc: c.agc.map(|a| AgcConfig {
            target_volts: a.target_volts,
            attack: a.attack,
            step: a.step,
            min_gain_db: a.min_gain_db,
            max_gain_db: a.max_gain_db,
        }),
        adc: c.adc_bits.map(|bits| AdcConfig { bits }),
    })
}
