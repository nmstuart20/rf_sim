//! Typed physical quantities and unit-aware parsing.
//! Per the plan section "Interfaces" (Scenario files), this crate owns quantities
//! such as frequency, power, distance, angle, and duration, including unit checks.
//! It provides the unit vocabulary used when resolving scenario inputs.
//!
//! Dependencies: iqsim-core.
//!
//! # Writing quantities
//!
//! Scenario values are strings with a case-sensitive unit: `"1090 MHz"`,
//! `"30 dBm"`, `"250 kt"`, `"35000 ft MSL"`. A bare number where a unit is required
//! is rejected with a message that suggests the spelling. The number is read as an
//! exact decimal (`[+-]digits[.digits][e[+-]digits]`, a space before the unit is
//! optional), and contracts §11 lists each representation:
//!
//! | Type | Units | Stored as |
//! | --- | --- | --- |
//! | [`Frequency`] | Hz, kHz, MHz, GHz | i64 µHz, exact |
//! | [`SampleRate`] | S/s, kS/s, MS/s, GS/s | core rational, exact |
//! | [`Duration`] | s, ms, us/µs, ns, ps | core i64 ps, exact |
//! | [`DataRate`] | b/s … Gb/s, B/s … GB/s | u64 bit/s, exact |
//! | [`Power`] | W, mW, dBm, dBW | f64 W |
//! | [`Gain`] | dB, dBi | f64 dB |
//! | [`Distance`] | m, km, ft, nmi | f64 m |
//! | [`Speed`] | m/s, km/h, kt | f64 m/s |
//! | [`Angle`] | deg, °, rad, DMS | f64 deg |
//! | [`Altitude`] | m, km, ft, nmi + MSL/AGL/ellipsoid, or FL | f64 m + reference |
//! | [`Position`] | lat/lon (decimal or DMS) or [`Mgrs`], plus altitude | f64 deg |
//!
//! Exact types reject values finer than their resolution instead of rounding them.
//! Float types round once to the nearest f64; transcendental conversions (dB, the
//! UTM inverse) use `libm`, so results match on every platform.

#![forbid(unsafe_code)]

mod angle;
mod decimal;
mod error;
mod exact;
mod float;
pub mod mgrs;
mod position;
mod text;
pub mod utm;

pub use angle::{Angle, Latitude, Longitude};
pub use error::{QuantityKind, UnitError};
pub use exact::{DataRate, Duration, Frequency, SampleRate};
pub use float::{Distance, Gain, Power, Speed};
pub use mgrs::{Mgrs, MgrsError, MgrsErrorKind};
pub use position::{Altitude, AltitudeRef, Position, PositionError};
