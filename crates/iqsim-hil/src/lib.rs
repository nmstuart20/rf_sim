//! Open-loop SDR integration and device calibration.
//! Per the plan section "Hardware in the loop", this crate owns antenna-port TX
//! and recorded RX injection, device pacing, synchronization, and calibration.
//! Virtual receiver effects remain outside the hardware antenna-port path.
//! Unsafe will be allowed only in narrowly scoped, documented FFI modules (libuhd).
//!
//! Dependencies: iqsim-core, iqsim-contracts, iqsim-io.

#![deny(unsafe_code)]
