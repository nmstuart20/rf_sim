//! IQ sinks, recordings, and separate privileged exports.
//! Per the plan section "Interfaces" (IQ outputs and Exports), this crate owns
//! strict framed delivery, VITA-49, preview, SigMF, truth, and knowledge exports.
//! Consumer-visible IQ metadata excludes simulation truth and run manifests.
//!
//! Dependencies: iqsim-core, iqsim-contracts.

#![forbid(unsafe_code)]
