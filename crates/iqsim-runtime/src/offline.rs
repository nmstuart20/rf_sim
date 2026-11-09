//! Offline runs and replay (`sim run --offline`, `sim replay`).
//!
//! # Run directory
//!
//! ```text
//! <out>/consumer/<sink>/<stream>.sigmf-data   consumer-visible recordings
//! <out>/consumer/<sink>/<stream>.sigmf-meta
//! <out>/run/manifest.json                      privileged: run manifest
//! <out>/run/truth.json                         privileged: canonical truth record
//! ```
//!
//! Consumer files hold only samples and SigMF core metadata. The manifest and the truth
//! record never go into `consumer/`.
//!
//! # Pipeline, per stream
//!
//! Static work list for block b → CPU render into zeroed i64 accumulators → the
//! complete sum goes once through the virtual-receiver chain ([`ReceiverChannelState`])
//! or, for an antenna port, the [`PortConverter`] → every SigMF sink of the stream.
//! Blocks run in order. The last block is rendered whole and truncated to the stream's
//! sample count before the chain sees it.
//!
//! # Replay
//!
//! [`replay`] checks, before rendering anything: the manifest schema, the numeric
//! contract version, the build identity, the scenario file hash, the resolved-scenario
//! hash, every asset and the absence of an event log (live commands arrive in prompt
//! 16). It then re-renders and requires every output's count and SHA-256 and the truth
//! hash to match.

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use iqsim_contracts::{
    AssetKind, AssetRegistry, BackendRecord, BuildIdentity, MANIFEST_SCHEMA_VERSION,
    NondeterministicRecord, OutputRecord, RunManifest, ScenarioRecord, TruthRecord, canonical_json,
    canonical_sha256, sha256_hex,
};
use iqsim_core::{
    BlockSpec, ClockMap, IdealTimestampClock, NUMERIC_CONTRACT_VERSION, SampleIndex, TapKind,
    TapResponse, TimestampClock,
};
use iqsim_io::IoError;
use iqsim_io::sigmf::{CheckedSampleSink, SampleFormat, SigmfCapture, SigmfGlobal, SigmfWriter};
use iqsim_receiver::{
    ChainError, OutputFormat, PortConverter, PortError, ReceiverChannelState, Samples,
};
use iqsim_render::{Backend, CpuBackend, RenderError};
use serde_json::{Value, json};

use crate::plan::{
    ResolveError, ResolvedGeoid, ResolvedScenario, RunPlan, StreamPlan, resolution_warnings,
    resolve,
};
use crate::static_scene::{SceneError, StaticSceneBuilder};

/// A failed run or replay.
#[derive(Debug)]
pub enum RunError {
    /// The scenario file could not be read, parsed or validated (rendered diagnostics).
    Scenario(String),
    /// The scenario cannot be run.
    Resolve(ResolveError),
    /// Work-list construction failed.
    Scene(SceneError),
    /// Rendering failed.
    Render(RenderError),
    /// The receiver chain failed.
    Chain(ChainError),
    /// The antenna-port conversion failed.
    Port(PortError),
    /// A file operation failed.
    Io(String),
    /// Replay rejected the manifest or found different output.
    Replay(String),
}

impl fmt::Display for RunError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Scenario(s) => f.write_str(s.trim_end()),
            Self::Resolve(e) => write!(f, "cannot run the scenario: {e}"),
            Self::Scene(e) => write!(f, "work-list construction failed: {e}"),
            Self::Render(e) => write!(f, "render failed: {e}"),
            Self::Chain(e) => write!(f, "receiver chain failed: {e}"),
            Self::Port(e) => write!(f, "antenna-port output failed: {e}"),
            Self::Io(e) => f.write_str(e),
            Self::Replay(e) => write!(f, "replay rejected: {e}"),
        }
    }
}

impl std::error::Error for RunError {}

macro_rules! from_err {
    ($($t:ty => $v:ident),*) => {$(
        impl From<$t> for RunError {
            fn from(e: $t) -> Self { Self::$v(e) }
        }
    )*};
}
from_err!(ResolveError => Resolve, SceneError => Scene, RenderError => Render,
    ChainError => Chain, PortError => Port);

impl From<IoError> for RunError {
    fn from(e: IoError) -> Self {
        Self::Io(e.to_string())
    }
}

fn io_err(path: &Path, e: std::io::Error) -> RunError {
    RunError::Io(format!("{}: {e}", path.display()))
}

/// Run options that cannot change output bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RunOptions {
    /// CPU render threads (≥ 1).
    pub threads: usize,
}

impl Default for RunOptions {
    fn default() -> Self {
        let threads = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
        Self { threads }
    }
}

/// The build identity of this binary (contracts §6).
pub fn build_identity() -> BuildIdentity {
    BuildIdentity {
        package_version: env!("CARGO_PKG_VERSION").to_owned(),
        rustc: env!("IQSIM_RUSTC_RELEASE").to_owned(),
        target: env!("IQSIM_TARGET").to_owned(),
        profile: env!("IQSIM_PROFILE").to_owned(),
        features: Vec::new(),
        cargo_lock_sha256: sha256_hex(include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../Cargo.lock"
        ))),
        source_sha256: env!("IQSIM_SOURCE_SHA256").to_owned(),
    }
}

/// A scenario file read, linted and resolved.
#[derive(Clone, Debug)]
pub struct LoadedScenario {
    /// Absolute path.
    pub path: PathBuf,
    /// File bytes.
    pub bytes: Vec<u8>,
    /// Resolved form.
    pub resolved: ResolvedScenario,
}

/// Reads, lints and resolves a scenario file.
pub fn load_scenario(path: &Path) -> Result<LoadedScenario, RunError> {
    let bytes = fs::read(path).map_err(|e| io_err(path, e))?;
    let abs = fs::canonicalize(path).map_err(|e| io_err(path, e))?;
    let shown = path.display().to_string();
    let source = String::from_utf8(bytes.clone())
        .map_err(|_| RunError::Scenario(format!("{shown}: error: the file is not UTF-8")))?;
    let scenario = iqsim_scenario::parse(&source)
        .map_err(|d| RunError::Scenario(d.render(&shown, &source)))?;
    let diagnostics = iqsim_scenario::validate(&scenario);
    if !diagnostics.is_empty() {
        let text: String = diagnostics
            .iter()
            .map(|d| d.render(&shown, &source))
            .collect();
        return Err(RunError::Scenario(text));
    }
    let resolved = resolve(&scenario, abs.parent().unwrap_or(Path::new("/")))?;
    Ok(LoadedScenario {
        path: abs,
        bytes,
        resolved,
    })
}

/// Output samples of one stream after the chain or port conversion, written to sinks
/// as they come.
fn render_stream(
    plan: &StreamPlan,
    backend: &CpuBackend,
    emit: &mut dyn FnMut(&Samples) -> Result<(), RunError>,
) -> Result<u64, RunError> {
    let r = &plan.resolved;
    let blocks = BlockSpec::new(r.block_len).map_err(SceneError::Time)?;
    let builder = StaticSceneBuilder::new(&plan.scene, &plan.target, &plan.clock, blocks)?
        .with_sample_limit(r.sample_count);
    let mut chain = match plan.chain {
        Some(c) => Some(ReceiverChannelState::new(c)?),
        None => None,
    };
    let port = PortConverter::new(plan.format);
    let total = r.sample_count;
    let nblocks = total.div_ceil(r.block_len);
    for b in 0..nblocks {
        let work = builder.build_block(b)?;
        let mut acc = vec![0_i64; 2 * work.len as usize];
        backend.render(&work, &mut acc)?;
        let first = b * r.block_len;
        let take = (total - first).min(r.block_len);
        acc.truncate(2 * take as usize);
        let samples = match &mut chain {
            Some(state) => {
                state
                    .process(b, SampleIndex::new(first), &acc, &plan.clock)?
                    .samples
            }
            None => port.convert(&acc)?,
        };
        emit(&samples)?;
    }
    Ok(chain.map_or(0, |c| c.clipped_total()))
}

fn sigmf_format(f: OutputFormat) -> SampleFormat {
    match f {
        OutputFormat::Ci16 => SampleFormat::Ci16Le,
        OutputFormat::Cf32 => SampleFormat::Cf32Le,
    }
}

/// The output rate of a stream (after decimation).
fn output_rate(plan: &StreamPlan) -> Result<iqsim_core::SampleRate, RunError> {
    match plan.chain {
        Some(c) => Ok(ReceiverChannelState::new(c)?.output_rate()),
        None => Ok(plan.rate()),
    }
}

fn tap_name(plan: &StreamPlan) -> &'static str {
    match plan.target.kind {
        TapKind::AntennaPort => "antenna_port",
        TapKind::VirtualReceiver => "virtual_receiver",
        TapKind::SourceReconstruction => unreachable!("streams are taps"),
    }
}

/// The canonical truth record: per stream and emitter, the link geometry and level,
/// mapped to the stream's sample indices. It is a privileged export (never in consumer
/// files).
pub fn truth_record(plan: &RunPlan) -> Result<Value, RunError> {
    let mut links = Vec::new();
    let mut streams = Vec::new();
    for s in &plan.streams {
        let t = &s.target;
        let resp = TapResponse::new(t.kind, t.response, t.rate, t.center_uhz)
            .map_err(SceneError::Response)?;
        let scale =
            iqsim_core::TargetScale::new(t.full_scale_volts).map_err(SceneError::Amplitude)?;
        let rate_hz = t.rate.num() as f64 / t.rate.den() as f64;
        // Index and time mapping: render index n ↔ acquisition instant (ClockMap) ↔
        // reported label (timestamp clock); output index j ↔ render index j·M, plus the
        // decimation filter's (L − 1)/2 render samples of delay.
        let (factor, decim_delay, out_rate) = match s.chain {
            Some(c) => {
                let state = ReceiverChannelState::new(c)?;
                let factor = c.decimation.map_or(1, |d| u64::from(d.factor));
                let delay = state
                    .decimation_taps()
                    .map_or(0, |t| (t.len() as u64 - 1) / 2);
                (factor, delay, state.output_rate())
            }
            None => (1, 0, t.rate),
        };
        let t0 = s
            .clock
            .instant(SampleIndex::ZERO)
            .map_err(SceneError::Time)?;
        let reported = IdealTimestampClock.label(t0).map_err(SceneError::Time)?;
        streams.push(json!({
            "stream": s.resolved.id,
            "tap": tap_name(s),
            "render_rate": [t.rate.num(), t.rate.den()],
            "render_sample_count": s.resolved.sample_count,
            "output_rate": [out_rate.num(), out_rate.den()],
            "output_sample_count": s.resolved.sample_count.div_ceil(factor),
            "output_index_to_render_index_factor": factor,
            "decimation_delay_render_samples": decim_delay,
            "first_sample_acquisition": [t0.ps().as_ps(), t0.frac()],
            "first_sample_reported": [reported.ps, reported.frac],
            "effective_full_scale_volts": scale.full_scale_effective_volts(),
        }));
        for (e, re) in s.scene.emitters.iter().zip(&plan.resolved.emitters) {
            let range = s.scene.receiver.distance(e.position);
            let carrier_hz = e.carrier_uhz as f64 * 1e-6;
            let lambda = iqsim_geo::SPEED_OF_LIGHT_M_S / carrier_hz;
            let path = lambda / (4.0 * core::f64::consts::PI * range);
            let offset_hz = (e.carrier_uhz - t.center_uhz) as f64 * 1e-6;
            let gain = resp.magnitude(offset_hz);
            let delay_s = range / iqsim_geo::SPEED_OF_LIGHT_M_S;
            let mut link = json!({
                "stream": s.resolved.id,
                "emitter": re.name,
                "range_m": range,
                "delay_s": delay_s,
                "delay_render_samples": delay_s * rate_hz,
                "carrier_offset_hz": offset_hz,
                "response_gain": gain,
                "power_at_target_w": e.power_w * (path * path) * (gain * gain),
            });
            if let Some(src) = &re.source {
                // A wideband source has no single response gain: report the level of
                // its reference (a full-scale or RMS tone) before band limiting, and
                // where its samples land.
                let native_hz = src.rate_num as f64 / src.rate_den as f64;
                let map = link.as_object_mut().expect("link is an object");
                map.remove("response_gain");
                map.remove("power_at_target_w");
                map.insert(
                    "reference_power_at_target_w".to_owned(),
                    json!(e.power_w * (path * path)),
                );
                map.insert(
                    "iq_file".to_owned(),
                    json!({
                        "start_ps": src.start_ps,
                        "native_rate": [src.rate_num, src.rate_den],
                        "delay_native_samples": delay_s * native_hz,
                        "sample_count": src.sample_count,
                        "end": src.end,
                        "group_delay_s": resp.group_delay_s(),
                    }),
                );
            }
            links.push(link);
        }
    }
    Ok(json!({ "truth_schema_version": 1, "streams": streams, "links": links }))
}

/// Renders every stream that has a sink. `make_sink` returns one sample sink per
/// (stream, sink name); `finish` closes it and returns (relative path, count, hash).
fn render_all<S>(
    plan: &RunPlan,
    options: RunOptions,
    mut make_sink: impl FnMut(&StreamPlan, &str) -> Result<S, RunError>,
    mut write: impl FnMut(&mut S, &Samples) -> Result<(), RunError>,
    mut finish: impl FnMut(S) -> Result<(String, u64, String), RunError>,
) -> Result<Vec<OutputRecord>, RunError> {
    let backend = CpuBackend::new(options.threads.max(1));
    let mut outputs = Vec::new();
    for s in &plan.streams {
        if s.resolved.sinks.is_empty() {
            continue;
        }
        let mut sinks = Vec::new();
        for name in &s.resolved.sinks {
            sinks.push((name.clone(), make_sink(s, name)?));
        }
        let clipped = render_stream(s, &backend, &mut |samples| {
            for (_, sink) in &mut sinks {
                write(sink, samples)?;
            }
            Ok(())
        })?;
        let rate = output_rate(s)?;
        for (name, sink) in sinks {
            let (path, count, sha) = finish(sink)?;
            outputs.push(OutputRecord {
                stream_id: s.resolved.id.clone(),
                sink: name,
                tap: tap_name(s).to_owned(),
                sample_format: s.resolved.format.clone(),
                rate_num: rate.num(),
                rate_den: rate.den(),
                center_uhz: s.resolved.center_uhz,
                first_sample_index: 0,
                sample_count: count,
                sha256: sha,
                path,
                finalized: true,
                adc_clipped_components: clipped,
            });
        }
    }
    Ok(outputs)
}

fn write_samples_sigmf(w: &mut SigmfWriter, s: &Samples) -> Result<(), RunError> {
    match s {
        Samples::Ci16(v) => w.write_ci16(v)?,
        Samples::Cf32(v) => w.write_cf32(v)?,
    }
    Ok(())
}

/// Every external input the resolved scenario names (contracts §6): per IQ-file
/// emitter its data file (with the parameters that interpret it as metadata) and its
/// `.sigmf-meta` file.
pub fn asset_registry(resolved: &ResolvedScenario) -> Result<AssetRegistry, RunError> {
    let mut reg = AssetRegistry::new();
    let sha = |hex_text: &str| -> Result<[u8; 32], RunError> {
        let mut out = [0u8; 32];
        hex::decode_to_slice(hex_text, &mut out)
            .map_err(|e| RunError::Io(format!("bad asset hash: {e}")))?;
        Ok(out)
    };
    for e in &resolved.emitters {
        let Some(src) = &e.source else { continue };
        let meta = json!({
            "datatype": src.datatype,
            "rate": [src.rate_num, src.rate_den],
            "sample_count": src.sample_count,
            "start_ps": src.start_ps,
            "end": src.end,
            "normalization": src.normalization,
            "scale": src.scale,
            "response": serde_json::to_value(src.response).expect("response serializes"),
            "meta_sha256": src.meta_sha256,
        });
        reg.register_digest(
            AssetKind::IqInput,
            &format!("{}.data", e.name),
            &src.data_path,
            src.data_byte_len,
            sha(&src.data_sha256)?,
            &format!("sigmf/{}", src.datatype),
            Some(&meta),
        )
        .map_err(|e| RunError::Io(e.to_string()))?;
        reg.register_digest(
            AssetKind::IqInput,
            &format!("{}.meta", e.name),
            &src.meta_path,
            src.meta_byte_len,
            sha(&src.meta_sha256)?,
            "sigmf-meta",
            None,
        )
        .map_err(|e| RunError::Io(e.to_string()))?;
    }
    if let ResolvedGeoid::Grid {
        path,
        byte_len,
        sha256,
        width,
        height,
        offset_m,
        scale_m,
        interpolation,
    } = &resolved.world.geoid
    {
        let meta = json!({
            "width": width,
            "height": height,
            "offset_m": offset_m,
            "scale_m": scale_m,
            "interpolation": interpolation,
        });
        reg.register_digest(
            AssetKind::GeoidGrid,
            "geoid",
            path,
            *byte_len,
            sha(sha256)?,
            iqsim_geo::geoid::GEOID_FORMAT,
            Some(&meta),
        )
        .map_err(|e| RunError::Io(e.to_string()))?;
    }
    Ok(reg)
}

/// Length and SHA-256 of a file, read in bounded chunks.
fn file_digest(path: &str) -> Result<(u64, [u8; 32]), std::io::Error> {
    use sha2::{Digest, Sha256};
    use std::io::Read;
    let mut f = fs::File::open(path)?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut len = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        len += n as u64;
        h.update(&buf[..n]);
    }
    Ok((len, h.finalize().into()))
}

fn unix_ms() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_millis() as u64)
}

/// Runs a scenario offline into `out` (which must not exist or be empty) and returns
/// the manifest it wrote to `out/run/manifest.json`.
pub fn run_offline(
    scenario: &Path,
    out: &Path,
    options: RunOptions,
) -> Result<RunManifest, RunError> {
    let start_ms = unix_ms();
    let started = Instant::now();
    let loaded = load_scenario(scenario)?;
    let plan = RunPlan::new(loaded.resolved.clone())?;
    if out.exists() {
        let mut entries = fs::read_dir(out).map_err(|e| io_err(out, e))?;
        if entries.next().is_some() {
            return Err(RunError::Io(format!(
                "{}: the output directory is not empty",
                out.display()
            )));
        }
    }
    let run_dir = out.join("run");
    fs::create_dir_all(&run_dir).map_err(|e| io_err(&run_dir, e))?;
    let epoch = (plan.resolved.epoch_unix_seconds, plan.resolved.epoch_nanos);
    let name = plan.resolved.name.clone();

    let outputs = render_all(
        &plan,
        options,
        |s, sink| {
            let dir = out.join("consumer").join(sink);
            fs::create_dir_all(&dir).map_err(|e| io_err(&dir, e))?;
            let t0 = s
                .clock
                .instant(SampleIndex::ZERO)
                .map_err(SceneError::Time)?;
            let reported = IdealTimestampClock.label(t0).map_err(SceneError::Time)?;
            let datetime = iqsim_io::utc::format_datetime(epoch.0, epoch.1, reported)?;
            let w = SigmfWriter::create(
                &dir,
                &s.resolved.id,
                sigmf_format(s.format),
                SigmfGlobal {
                    sample_rate: output_rate(s)?,
                    description: Some(format!("{name}: {}", s.resolved.id)),
                    recorder: format!("iqsim {}", env!("CARGO_PKG_VERSION")),
                },
                vec![SigmfCapture {
                    sample_start: 0,
                    frequency_hz: s.resolved.center_uhz as f64 / 1e6,
                    datetime: Some(datetime),
                }],
            )?;
            Ok((format!("consumer/{sink}/{}.sigmf-data", s.resolved.id), w))
        },
        |(_, w), samples| write_samples_sigmf(w, samples),
        |(path, w)| {
            let summary = w.finish()?;
            Ok((path, summary.sample_count, summary.sha256_hex()))
        },
    )?;

    let truth = canonical_json(&truth_record(&plan)?);
    let truth_path = run_dir.join("truth.json");
    fs::write(&truth_path, &truth).map_err(|e| io_err(&truth_path, e))?;

    let resolved = serde_json::to_value(&plan.resolved).expect("resolved scenario serializes");
    let manifest = RunManifest {
        manifest_schema_version: MANIFEST_SCHEMA_VERSION,
        numeric_contract_version: NUMERIC_CONTRACT_VERSION,
        build: build_identity(),
        scenario: ScenarioRecord {
            source_path: loaded.path.display().to_string(),
            source_sha256: sha256_hex(&loaded.bytes),
            resolved_sha256: canonical_sha256(&resolved),
            resolved,
        },
        seed: plan.resolved.seed,
        run_mode: "offline".to_owned(),
        block_ps: plan.resolved.block_ps,
        duration_ps: plan.resolved.duration_ps,
        assets: asset_registry(&plan.resolved)?.records(),
        warnings: resolution_warnings(&plan.resolved),
        event_log: None,
        outputs,
        truth: TruthRecord {
            path: "run/truth.json".to_owned(),
            sha256: sha256_hex(&truth),
        },
        backend: BackendRecord {
            name: "cpu".to_owned(),
            threads: options.threads.max(1) as u32,
        },
        nondeterministic: NondeterministicRecord {
            wall_clock_start_unix_ms: start_ms,
            wall_clock_end_unix_ms: unix_ms(),
            host: None,
            render_seconds: started.elapsed().as_secs_f64(),
        },
    };
    let manifest_path = run_dir.join("manifest.json");
    fs::write(&manifest_path, manifest.to_json()).map_err(|e| io_err(&manifest_path, e))?;
    Ok(manifest)
}

/// What replay verified.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReplayReport {
    /// Streams re-rendered and matched.
    pub outputs_matched: usize,
    /// The manifest's deterministic hash.
    pub deterministic_sha256: String,
}

/// Checks a run manifest's inputs, re-renders the run and compares every output and
/// the truth record. Nothing is rendered unless every input check passes.
pub fn replay(manifest_path: &Path, options: RunOptions) -> Result<ReplayReport, RunError> {
    let text = fs::read_to_string(manifest_path).map_err(|e| io_err(manifest_path, e))?;
    let m = RunManifest::from_json(&text).map_err(RunError::Replay)?;
    if m.numeric_contract_version != NUMERIC_CONTRACT_VERSION {
        return Err(RunError::Replay(format!(
            "numeric contract version {} differs from this build's {NUMERIC_CONTRACT_VERSION}",
            m.numeric_contract_version
        )));
    }
    let build = build_identity();
    if m.build != build {
        return Err(RunError::Replay(format!(
            "the run was made by a different build ({:?}); this is {:?}. Compare against stored goldens instead",
            m.build, build
        )));
    }
    if m.event_log.is_some() {
        return Err(RunError::Replay(
            "runs with live commands need event-log replay (prompt 16)".to_owned(),
        ));
    }
    let source = Path::new(&m.scenario.source_path);
    let bytes = fs::read(source).map_err(|e| {
        RunError::Replay(format!("scenario {} cannot be read: {e}", source.display()))
    })?;
    if sha256_hex(&bytes) != m.scenario.source_sha256 {
        return Err(RunError::Replay(format!(
            "scenario {} changed since the run",
            source.display()
        )));
    }
    if canonical_sha256(&m.scenario.resolved) != m.scenario.resolved_sha256 {
        return Err(RunError::Replay(
            "the manifest's resolved scenario does not match its hash".to_owned(),
        ));
    }
    AssetRegistry::verify_digests(&m.assets, file_digest)
        .map_err(|e| RunError::Replay(e.to_string()))?;
    let loaded = load_scenario(source)?;
    let resolved = serde_json::to_value(&loaded.resolved).expect("resolved scenario serializes");
    if canonical_sha256(&resolved) != m.scenario.resolved_sha256 {
        return Err(RunError::Replay(
            "the scenario resolves differently in this build".to_owned(),
        ));
    }
    if m.seed != loaded.resolved.seed
        || m.block_ps != loaded.resolved.block_ps
        || m.duration_ps != loaded.resolved.duration_ps
    {
        return Err(RunError::Replay(
            "seed, block or duration differ from the resolved scenario".to_owned(),
        ));
    }
    if m.warnings != resolution_warnings(&loaded.resolved) {
        return Err(RunError::Replay(
            "the manifest's warnings differ from the resolved scenario's".to_owned(),
        ));
    }
    // Every asset the resolved scenario needs, with the same identity and metadata
    // (review fix: an omitted or edited record used to pass, because only the
    // listed files were verified).
    if asset_registry(&loaded.resolved)?.records() != m.assets {
        return Err(RunError::Replay(
            "the manifest's asset list differs from what the scenario resolves to".to_owned(),
        ));
    }

    // Every input checks out: re-render into hashing sinks.
    let plan = RunPlan::new(loaded.resolved)?;
    let outputs = render_all(
        &plan,
        options,
        |s, sink| {
            Ok((
                format!("consumer/{sink}/{}.sigmf-data", s.resolved.id),
                CheckedSampleSink::new(std::io::sink(), sigmf_format(s.format)),
            ))
        },
        |(_, w), samples| {
            match samples {
                Samples::Ci16(v) => w.write_ci16(v)?,
                Samples::Cf32(v) => w.write_cf32(v)?,
            }
            Ok(())
        },
        |(path, w)| Ok((path, w.sample_count(), hex::encode(w.sha256()))),
    )?;
    if outputs != m.outputs {
        let which = outputs
            .iter()
            .zip(&m.outputs)
            .find(|(a, b)| a != b)
            .map_or("the stream list".to_owned(), |(a, _)| {
                format!("stream \"{}\"", a.stream_id)
            });
        return Err(RunError::Replay(format!("{which} re-rendered differently")));
    }
    let truth = canonical_json(&truth_record(&plan)?);
    if sha256_hex(&truth) != m.truth.sha256 {
        return Err(RunError::Replay("the truth record differs".to_owned()));
    }
    Ok(ReplayReport {
        outputs_matched: outputs.len(),
        deterministic_sha256: m.deterministic_sha256(),
    })
}
