//! Versioned shared interfaces without handlers or I/O.
//! Per the plan section "Interfaces", this crate owns message envelopes,
//! command/ack/status schemas, IQ frame headers and context, and run manifests.
//! It also owns the handler traits implemented
//! by runtime and consumed by adapters, keeping control independent of runtime.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]

pub mod manifest;

pub use manifest::{
    AssetKind, AssetMismatch, AssetRecord, AssetRegistry, BackendRecord, BuildIdentity,
    EventLogRecord, MANIFEST_SCHEMA_VERSION, NondeterministicRecord, OutputRecord, RunManifest,
    ScenarioRecord, TruthRecord, canonical_json, canonical_sha256, sha256_hex,
};
