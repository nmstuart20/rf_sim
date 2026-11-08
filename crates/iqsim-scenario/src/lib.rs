//! Scenario schema, validation, and resolution.
//! Per the plan section "Interfaces" (Scenario files), this crate owns strict
//! unit-aware configuration and named references to platforms, emitters, receivers,
//! libraries, output sinks, and shared resource capacities.
//!
//! Dependencies: iqsim-core, iqsim-units, iqsim-contracts.

#![forbid(unsafe_code)]
