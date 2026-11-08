//! Stateful emitter behavior and deterministic event scheduling.
//! Per the plan section "Emitters, behaviors and waveforms", this crate owns
//! transmission decisions and protocol timing, including SSR and transponder events.
//! Events and accepted commands follow deterministic ordering and commit boundaries.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts, iqsim-scenario, iqsim-geo.

#![forbid(unsafe_code)]
