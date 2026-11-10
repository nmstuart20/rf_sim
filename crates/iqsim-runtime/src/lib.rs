//! Stateful CPU orchestration and ordered block execution.
//! Per the plan section "Architecture", this crate connects resolved scenarios,
//! behavior, geometry, work-list construction, rendering, receiver state, and sinks.
//! It owns shared commit-boundary admission, event logging, and verified replay
//! without depending on transport adapters or hardware integration.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts, iqsim-scenario, iqsim-geo, iqsim-antenna, iqsim-terrain, iqsim-behavior, iqsim-waveform, iqsim-worklist, iqsim-render, iqsim-receiver, iqsim-io.

#![forbid(unsafe_code)]

pub mod geometry;
pub mod offline;
pub mod plan;
pub mod scene;
pub mod trajectory;

pub use offline::{
    GeometryLog, LoadedScenario, PathFitSummary, ReplayReport, RunError, RunOptions,
    TRUTH_SCHEMA_VERSION, build_identity, load_scenario, replay, run_offline, truth_record,
};
pub use plan::{ResolveError, ResolvedScenario, RunPlan, StreamPlan, resolve};
