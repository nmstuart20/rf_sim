//! IQ sinks, recordings, and separate privileged exports.
//! Per the plan section "Interfaces" (IQ outputs and Exports), this crate owns
//! strict framed delivery, VITA-49, preview, SigMF, truth, and knowledge exports.
//! Consumer-visible IQ metadata excludes simulation truth and run manifests.
//!
//! [`sigmf`] provides checked recording writes and readback helpers; [`utc`]
//! formats receiver-reported timestamps using integer calendar arithmetic.
//! [`iq_input`] validates SigMF emitter inputs and reads them through a bounded,
//! hash-verified page cache.
//!
//! Dependencies: iqsim-core, iqsim-contracts.

#![forbid(unsafe_code)]

mod error;
pub mod iq_input;
pub mod sigmf;
pub mod utc;

pub use error::IoError;
