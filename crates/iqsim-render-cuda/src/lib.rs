//! CUDA host integration for deterministic render kernels.
//! Per the plan section "Compute backends and scaling", this crate owns CUDA
//! backend loading and partition execution matching the CPU reference bit for bit.
//! Emitter partitions return raw i64 sums before receiver effects or finalization.
//! Unsafe will be allowed only in narrowly scoped, documented FFI modules (cudarc).
//!
//! Dependencies: iqsim-core, iqsim-worklist, iqsim-render.

#![deny(unsafe_code)]
