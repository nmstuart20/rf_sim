//! Stateful CPU orchestration and ordered block execution.
//! Per the plan section "Architecture", this crate connects resolved scenarios,
//! behavior, geometry, work-list construction, rendering, receiver state, and sinks.
//! It owns shared commit-boundary admission, event logging, and verified replay
//! without depending on transport adapters or hardware integration.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts, iqsim-scenario, iqsim-geo, iqsim-antenna, iqsim-terrain, iqsim-behavior, iqsim-waveform, iqsim-worklist, iqsim-render, iqsim-receiver, iqsim-io.

#![forbid(unsafe_code)]

pub mod offline;
pub mod plan;
pub mod static_scene;
pub mod trajectory;

pub use offline::{
    LoadedScenario, ReplayReport, RunError, RunOptions, build_identity, load_scenario, replay,
    run_offline, truth_record,
};
pub use plan::{ResolveError, ResolvedScenario, RunPlan, StreamPlan, resolve};
