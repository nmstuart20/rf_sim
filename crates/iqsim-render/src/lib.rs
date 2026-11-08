//! Render backend interface and CPU reference backend.
//! Per the plan section "Compute backends and scaling", this crate consumes flat
//! block work lists with stateless kernels and checked i64 accumulation.
//! CPU output supplies the bit-exact reference for every alternate backend.
//!
//! Dependencies: iqsim-core, iqsim-worklist, iqsim-waveform.

#![forbid(unsafe_code)]
