//! The versioned run manifest, the asset registry and canonical JSON (contracts §6).
//!
//! A manifest is a privileged run artifact. It is written beside the run, never into a
//! consumer-visible file, and `sim replay` checks it before rendering anything.
//!
//! # Canonical JSON (iqsim canonical JSON v1)
//!
//! Hashed documents (the resolved scenario, the truth record, the deterministic part of
//! a manifest) are serialized by [`canonical_json`]: no whitespace, object keys sorted by
//! their UTF-8 bytes, strings escaped as `serde_json` escapes them, integers in decimal,
//! and floats in `serde_json`'s shortest round-trip form. It follows RFC 8785 except for
//! float formatting, which only matters for floats outside 10⁻⁶..10²¹; exact quantities
//! (µHz, ps, rates) are integers, so they never reach that difference.
//!
//! # Deterministic and nondeterministic content
//!
//! [`RunManifest::deterministic_sha256`] hashes the canonical JSON of the manifest without
//! its `backend` and `nondeterministic` sections. Wall-clock times, host, thread counts
//! and timing live only there, so they never change a deterministic hash.

use std::collections::BTreeMap;
use std::fmt;

use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};

/// Manifest format version this build writes and reads.
pub const MANIFEST_SCHEMA_VERSION: u32 = 1;

/// Writes `value` as iqsim canonical JSON v1 (see the module docs).
pub fn canonical_json(value: &Value) -> Vec<u8> {
    let mut out = Vec::new();
    write_canonical(value, &mut out);
    out
}

fn write_canonical(value: &Value, out: &mut Vec<u8>) {
    match value {
        Value::Object(map) => {
            let mut keys: Vec<&String> = map.keys().collect();
            keys.sort_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
            out.push(b'{');
            for (i, k) in keys.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(
                    serde_json::to_string(k)
                        .expect("a string serializes")
                        .as_bytes(),
                );
                out.push(b':');
                write_canonical(&map[*k], out);
            }
            out.push(b'}');
        }
        Value::Array(items) => {
            out.push(b'[');
            for (i, v) in items.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                write_canonical(v, out);
            }
            out.push(b']');
        }
        other => out.extend_from_slice(
            serde_json::to_string(other)
                .expect("a scalar serializes")
                .as_bytes(),
        ),
    }
}

/// SHA-256 of `bytes` as lower-case hex.
pub fn sha256_hex(bytes: &[u8]) -> String {
    hex::encode(Sha256::digest(bytes))
}

/// SHA-256 of the canonical JSON of `value`, as hex.
pub fn canonical_sha256(value: &Value) -> String {
    sha256_hex(&canonical_json(value))
}

/// The build that produced a run. Replay requires every field to match.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BuildIdentity {
    /// Workspace package version.
    pub package_version: String,
    /// `rustc -vV` release line of the compiler.
    pub rustc: String,
    /// Target triple.
    pub target: String,
    /// Cargo profile (`debug` or `release`).
    pub profile: String,
    /// Enabled build features, sorted.
    pub features: Vec<String>,
    /// SHA-256 of the workspace `Cargo.lock` the binary was built from.
    pub cargo_lock_sha256: String,
    /// Fingerprint of the workspace sources the binary was built from (see
    /// `iqsim-runtime/build.rs`), so a code change is a different build.
    pub source_sha256: String,
}

/// The scenario a run used.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ScenarioRecord {
    /// Absolute path of the scenario file.
    pub source_path: String,
    /// SHA-256 of the scenario file bytes.
    pub source_sha256: String,
    /// The resolved configuration (every default filled in).
    pub resolved: Value,
    /// SHA-256 of the canonical JSON of `resolved`.
    pub resolved_sha256: String,
}

/// Kinds of external input (contracts §6). Appended to, never renumbered.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AssetKind {
    /// A library item (antenna, waveform, behavior, receiver profile, …).
    LibraryItem,
    /// A terrain tile.
    TerrainTile,
    /// A geoid grid.
    GeoidGrid,
    /// An antenna pattern file.
    AntennaPattern,
    /// An IQ input file.
    IqInput,
    /// A calibration table.
    CalibrationTable,
    /// A finalized live recording.
    LiveRecording,
    /// A finalized live timing/gap record.
    LiveTimingRecord,
}

/// One external input and its content identity.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssetRecord {
    /// What kind of input it is.
    pub kind: AssetKind,
    /// The name the scenario uses for it.
    pub logical_name: String,
    /// Where it was read from.
    pub resolved_path: String,
    /// Length in bytes.
    pub byte_len: u64,
    /// SHA-256 of the bytes.
    pub sha256: String,
    /// Format and version, e.g. `sigmf-1.2.0/ci16_le`.
    pub format: String,
    /// SHA-256 of the canonical JSON of the metadata that changes how the bytes are
    /// interpreted (rate, centre, scaling, …), or `None` when there is none.
    pub metadata_sha256: Option<String>,
}

/// An asset failed verification.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AssetMismatch {
    /// The asset could not be read.
    Missing {
        /// Its logical name.
        logical_name: String,
        /// The read error.
        error: String,
    },
    /// Its length or hash differs from the record.
    Changed {
        /// Its logical name.
        logical_name: String,
    },
    /// Two registrations of one (kind, name) disagree.
    Conflict {
        /// Its logical name.
        logical_name: String,
    },
}

impl fmt::Display for AssetMismatch {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing {
                logical_name,
                error,
            } => write!(f, "asset \"{logical_name}\" is missing: {error}"),
            Self::Changed { logical_name } => {
                write!(
                    f,
                    "asset \"{logical_name}\" changed since the run was recorded"
                )
            }
            Self::Conflict { logical_name } => {
                write!(
                    f,
                    "asset \"{logical_name}\" was registered twice with different content"
                )
            }
        }
    }
}

impl std::error::Error for AssetMismatch {}

/// The registry of every external input a run read, ordered by (kind, name).
///
/// Loaders register each input as they read it. The registry holds no I/O: callers
/// pass the bytes in, and [`Self::verify`] takes a reader.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AssetRegistry {
    assets: BTreeMap<(AssetKind, String), AssetRecord>,
}

impl AssetRegistry {
    /// An empty registry.
    pub fn new() -> Self {
        Self::default()
    }

    /// Registers `bytes`, read from `path`, as asset (`kind`, `logical_name`).
    /// Registering the same content twice is fine; different content is a conflict.
    pub fn register(
        &mut self,
        kind: AssetKind,
        logical_name: &str,
        path: &str,
        bytes: &[u8],
        format: &str,
        metadata: Option<&Value>,
    ) -> Result<(), AssetMismatch> {
        let record = AssetRecord {
            kind,
            logical_name: logical_name.to_owned(),
            resolved_path: path.to_owned(),
            byte_len: bytes.len() as u64,
            sha256: sha256_hex(bytes),
            format: format.to_owned(),
            metadata_sha256: metadata.map(canonical_sha256),
        };
        match self.assets.get(&(kind, logical_name.to_owned())) {
            Some(existing) if *existing != record => Err(AssetMismatch::Conflict {
                logical_name: logical_name.to_owned(),
            }),
            Some(_) => Ok(()),
            None => {
                self.assets.insert((kind, logical_name.to_owned()), record);
                Ok(())
            }
        }
    }

    /// Registers an asset whose length and SHA-256 were computed by a streaming pass
    /// (large files), with the same conflict rule as [`Self::register`].
    #[allow(clippy::too_many_arguments)]
    pub fn register_digest(
        &mut self,
        kind: AssetKind,
        logical_name: &str,
        path: &str,
        byte_len: u64,
        sha256: [u8; 32],
        format: &str,
        metadata: Option<&Value>,
    ) -> Result<(), AssetMismatch> {
        let record = AssetRecord {
            kind,
            logical_name: logical_name.to_owned(),
            resolved_path: path.to_owned(),
            byte_len,
            sha256: hex::encode(sha256),
            format: format.to_owned(),
            metadata_sha256: metadata.map(canonical_sha256),
        };
        self.insert(record)
    }

    fn insert(&mut self, record: AssetRecord) -> Result<(), AssetMismatch> {
        let key = (record.kind, record.logical_name.clone());
        match self.assets.get(&key) {
            Some(existing) if *existing != record => Err(AssetMismatch::Conflict {
                logical_name: record.logical_name,
            }),
            Some(_) => Ok(()),
            None => {
                self.assets.insert(key, record);
                Ok(())
            }
        }
    }

    /// Checks every record against the (length, SHA-256) that `digest` computes for its
    /// path, so large files can be verified without holding them in memory.
    pub fn verify_digests<E: fmt::Display>(
        records: &[AssetRecord],
        mut digest: impl FnMut(&str) -> Result<(u64, [u8; 32]), E>,
    ) -> Result<(), AssetMismatch> {
        for r in records {
            let (len, sha) = digest(&r.resolved_path).map_err(|e| AssetMismatch::Missing {
                logical_name: r.logical_name.clone(),
                error: e.to_string(),
            })?;
            if len != r.byte_len || hex::encode(sha) != r.sha256 {
                return Err(AssetMismatch::Changed {
                    logical_name: r.logical_name.clone(),
                });
            }
        }
        Ok(())
    }

    /// The records, in (kind, name) order.
    pub fn records(&self) -> Vec<AssetRecord> {
        self.assets.values().cloned().collect()
    }

    /// Checks every record against the bytes `read` returns for its path.
    pub fn verify<E: fmt::Display>(
        records: &[AssetRecord],
        mut read: impl FnMut(&str) -> Result<Vec<u8>, E>,
    ) -> Result<(), AssetMismatch> {
        for r in records {
            let bytes = read(&r.resolved_path).map_err(|e| AssetMismatch::Missing {
                logical_name: r.logical_name.clone(),
                error: e.to_string(),
            })?;
            if bytes.len() as u64 != r.byte_len || sha256_hex(&bytes) != r.sha256 {
                return Err(AssetMismatch::Changed {
                    logical_name: r.logical_name.clone(),
                });
            }
        }
        Ok(())
    }
}

/// The durable event log of a run with live commands (prompt 16).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EventLogRecord {
    /// Path of the log.
    pub path: String,
    /// SHA-256 of the log bytes.
    pub sha256: String,
    /// Number of records.
    pub record_count: u64,
    /// First ingress sequence.
    pub first_sequence: u64,
    /// Last ingress sequence.
    pub last_sequence: u64,
    /// Final commit boundary C, ps.
    pub final_commit_ps: i64,
}

/// One output stream of a run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputRecord {
    /// Stream id, e.g. `rx.ch0` or `port`.
    pub stream_id: String,
    /// The sink it was written to.
    pub sink: String,
    /// `antenna_port` or `virtual_receiver`.
    pub tap: String,
    /// `ci16_le` or `cf32_le`.
    pub sample_format: String,
    /// Output sample rate numerator, S/s.
    pub rate_num: u64,
    /// Output sample rate denominator.
    pub rate_den: u64,
    /// Centre frequency, µHz.
    pub center_uhz: i64,
    /// First sample index.
    pub first_sample_index: u64,
    /// Number of complex samples.
    pub sample_count: u64,
    /// SHA-256 of the sample bytes.
    pub sha256: String,
    /// Data file path relative to the run directory.
    pub path: String,
    /// The recording was finalized (count and hash checked at close).
    pub finalized: bool,
    /// Components clipped by the virtual ADC (a reported receiver effect; always 0 on
    /// an antenna port, where clipping is an error).
    pub adc_clipped_components: u64,
}

/// The canonical truth record of a run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TruthRecord {
    /// Path relative to the run directory.
    pub path: String,
    /// SHA-256 of the canonical truth bytes.
    pub sha256: String,
}

/// Rendering setup, recorded for diagnostics only. By contract it cannot change a
/// sample or truth hash, so it is excluded from the deterministic hash.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BackendRecord {
    /// Backend name (`cpu`).
    pub name: String,
    /// Worker threads.
    pub threads: u32,
}

/// Measurements that differ between identical runs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NondeterministicRecord {
    /// Wall-clock start, ms since the Unix epoch.
    pub wall_clock_start_unix_ms: u64,
    /// Wall-clock end, ms since the Unix epoch.
    pub wall_clock_end_unix_ms: u64,
    /// Host name, when known.
    pub host: Option<String>,
    /// Render wall time, seconds.
    pub render_seconds: f64,
}

/// A run manifest (contracts §6).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunManifest {
    /// [`MANIFEST_SCHEMA_VERSION`].
    pub manifest_schema_version: u32,
    /// `iqsim_core::NUMERIC_CONTRACT_VERSION` of the build.
    pub numeric_contract_version: u32,
    /// The build.
    pub build: BuildIdentity,
    /// The scenario.
    pub scenario: ScenarioRecord,
    /// Scenario seed.
    pub seed: u64,
    /// `offline` or `realtime`.
    pub run_mode: String,
    /// Block length, ps.
    pub block_ps: i64,
    /// Rendered duration, ps.
    pub duration_ps: i64,
    /// External inputs, in (kind, name) order.
    pub assets: Vec<AssetRecord>,
    /// Event log, when the run had live commands.
    pub event_log: Option<EventLogRecord>,
    /// Output streams, in stream order.
    pub outputs: Vec<OutputRecord>,
    /// Canonical truth record.
    pub truth: TruthRecord,
    /// Diagnostic backend setup (not hashed).
    pub backend: BackendRecord,
    /// Wall-clock measurements (not hashed).
    pub nondeterministic: NondeterministicRecord,
}

impl RunManifest {
    /// The manifest without `backend` and `nondeterministic`, as JSON.
    pub fn deterministic_view(&self) -> Value {
        let mut v = serde_json::to_value(self).expect("a manifest serializes");
        if let Value::Object(map) = &mut v {
            map.remove("backend");
            map.remove("nondeterministic");
        }
        v
    }

    /// SHA-256 of the canonical deterministic view.
    pub fn deterministic_sha256(&self) -> String {
        canonical_sha256(&self.deterministic_view())
    }

    /// Pretty JSON with a trailing newline, as written to `manifest.json`.
    pub fn to_json(&self) -> String {
        let mut s = serde_json::to_string_pretty(self).expect("a manifest serializes");
        s.push('\n');
        s
    }

    /// Parses a manifest, rejecting unknown fields and other schema versions.
    pub fn from_json(text: &str) -> Result<Self, String> {
        let m: Self = serde_json::from_str(text).map_err(|e| e.to_string())?;
        if m.manifest_schema_version != MANIFEST_SCHEMA_VERSION {
            return Err(format!(
                "manifest schema version {} is not supported; this build reads {MANIFEST_SCHEMA_VERSION}",
                m.manifest_schema_version
            ));
        }
        Ok(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_json_sorts_keys_and_drops_whitespace() {
        let v: Value =
            serde_json::from_str(r#"{ "b": [1, 2.5, {"z": null, "a": true}], "a": "x\n" }"#)
                .unwrap();
        assert_eq!(
            canonical_json(&v),
            br#"{"a":"x\n","b":[1,2.5,{"a":true,"z":null}]}"#.to_vec()
        );
        assert_eq!(
            canonical_sha256(&v),
            sha256_hex(br#"{"a":"x\n","b":[1,2.5,{"a":true,"z":null}]}"#)
        );
    }

    #[test]
    fn registry_orders_detects_conflicts_and_verifies() {
        let mut r = AssetRegistry::new();
        r.register(AssetKind::IqInput, "b", "/b", b"bb", "raw", None)
            .unwrap();
        r.register(AssetKind::LibraryItem, "a", "/a", b"a", "toml", None)
            .unwrap();
        r.register(AssetKind::IqInput, "b", "/b", b"bb", "raw", None)
            .unwrap();
        assert!(matches!(
            r.register(AssetKind::IqInput, "b", "/b", b"bc", "raw", None),
            Err(AssetMismatch::Conflict { .. })
        ));
        let recs = r.records();
        assert_eq!(recs[0].logical_name, "a");
        assert_eq!(recs[1].byte_len, 2);
        let ok = |p: &str| -> Result<Vec<u8>, String> {
            Ok(if p == "/a" {
                b"a".to_vec()
            } else {
                b"bb".to_vec()
            })
        };
        AssetRegistry::verify(&recs, ok).unwrap();
        let changed = |p: &str| -> Result<Vec<u8>, String> {
            Ok(if p == "/a" {
                b"A".to_vec()
            } else {
                b"bb".to_vec()
            })
        };
        assert_eq!(
            AssetRegistry::verify(&recs, changed),
            Err(AssetMismatch::Changed {
                logical_name: "a".into()
            })
        );
        let missing = |_: &str| -> Result<Vec<u8>, String> { Err("gone".into()) };
        assert!(matches!(
            AssetRegistry::verify(&recs, missing),
            Err(AssetMismatch::Missing { .. })
        ));
    }
}
