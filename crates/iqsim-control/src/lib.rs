//! Transport adapters for the shared control interface.
//! Per the plan section "Interfaces" (Control plane), this crate adapts commands,
//! acknowledgements, and status to transports and dispatches to a contracts handler trait.
//! It defines no simulation logic and has no dependency on runtime.
//!
//! Dependencies: iqsim-core, iqsim-contracts.

#![forbid(unsafe_code)]
