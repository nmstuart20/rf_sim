//! Reference frames, trajectories, and propagation geometry.
//! Per the plan sections "World model" and "Antennas and propagation", this crate
//! owns WGS84/ECEF/ECI/NED transforms, C2 motion, attitude, geoid altitude,
//! and light-time geometry that retains fractional acquisition time.
//!
//! Since prompt 7 (contracts §12):
//!
//! - [`vec`](mod@vec): `Vec3` and `Mat3` with a fixed operation order.
//! - [`geodetic`]: WGS84 geodetic ↔ ECEF (closed form forward, Bowring's iteration
//!   inverse).
//! - [`eci`]: the pseudo-inertial ECI frame, Earth rotation with an integer-reduced
//!   angle.
//! - [`quat`]: Hamilton quaternions (`q_a_b` maps frame-b coordinates to frame a).
//! - [`ned`] and [`frames`]: local NED, body and antenna frames.
//! - [`geoid`]: GeographicLib PGM geoid grids with bilinear interpolation.
//! - [`atmosphere`], [`gillham`] and [`altitude`]: ISA pressure altitude, Gillham
//!   (Mode C) encoding and the ellipsoid → MSL → pressure-altitude → code chain.
//!
//! Since prompt 8 (contracts §13): [`trajectory`], the [`Trajectory`] trait with static
//! platforms, maneuver primitives, waypoint routes, track replay, ground-clamped
//! motion, derived attitude and runtime appends.
//!
//! All arithmetic is f64 with transcendentals from `libm`, so results are the same
//! bits on every machine.
//!
//! Dependencies: iqsim-core, iqsim-units.

#![forbid(unsafe_code)]

pub mod altitude;
pub mod atmosphere;
pub mod eci;
pub mod frames;
pub mod geodetic;
pub mod geoid;
pub mod gillham;
pub mod ned;
pub mod quat;
pub mod trajectory;
pub mod vec;

use std::fmt;

pub use altitude::AltitudeModel;
pub use atmosphere::Atmosphere;
pub use eci::{EARTH_ROTATION_RATE_RAD_S, EarthOrientation, Eci, earth_rotation_angle};
pub use frames::{Mount, PlatformPose};
pub use geodetic::{Ecef, Geodetic, SPEED_OF_LIGHT_M_S, WGS84_A_M, WGS84_F};
pub use geoid::{Geoid, GeoidGrid, GeoidGridMetadata};
pub use gillham::GillhamCode;
pub use ned::LocalNed;
pub use quat::Quat;
pub use trajectory::{Flight, Span, State, Static, Trajectory};
pub use vec::{Mat3, Vec3};

/// A geometry, geoid or altitude-chain failure.
#[derive(Clone, Debug, PartialEq)]
pub enum GeoError {
    /// An input was NaN or infinite.
    NonFinite(&'static str),
    /// The geodetic inverse was asked for a point too close to the Earth's centre.
    InsideDomain {
        /// Distance from the centre, metres.
        radius_m: f64,
    },
    /// The geodetic inverse was asked for a point too far from the Earth's centre.
    OutsideDomain {
        /// Distance from the centre, metres.
        radius_m: f64,
    },
    /// An iteration did not converge within its cap.
    NoConvergence {
        /// Which iteration.
        what: &'static str,
        /// Steps taken.
        iterations: u32,
    },
    /// A geoid grid file is malformed.
    GeoidFormat(String),
    /// A file could not be read.
    Io(String),
    /// A value is outside the modelled atmosphere.
    AtmosphereRange {
        /// The value.
        value: f64,
        /// What it is, with units.
        what: &'static str,
    },
    /// A pressure altitude has no Gillham code.
    Gillham(gillham::GillhamError),
    /// A trajectory was evaluated outside its span.
    OutsideSpan {
        /// The instant asked for (floor), ps.
        t_ps: i64,
        /// Span start, ps.
        start_ps: i64,
        /// Span end, ps.
        end_ps: i64,
    },
    /// A trajectory could not be built or evaluated (motion plan, attitude).
    Trajectory(String),
}

impl fmt::Display for GeoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite(what) => write!(f, "{what} is not finite"),
            Self::InsideDomain { radius_m } => write!(
                f,
                "a point {radius_m} m from the Earth's centre is inside the geodetic inverse's domain limit ({} m)",
                geodetic::MIN_INVERSE_RADIUS_M
            ),
            Self::OutsideDomain { radius_m } => write!(
                f,
                "a point {radius_m} m from the Earth's centre is beyond the geodetic inverse's domain limit ({} m)",
                geodetic::MAX_INVERSE_RADIUS_M
            ),
            Self::NoConvergence { what, iterations } => {
                write!(f, "{what} did not converge in {iterations} iterations")
            }
            Self::GeoidFormat(m) => write!(f, "geoid grid: {m}"),
            Self::Io(m) => write!(f, "{m}"),
            Self::AtmosphereRange { value, what } => {
                write!(f, "{what} {value} is outside the modelled ISA range")
            }
            Self::Gillham(e) => write!(f, "Mode C: {e}"),
            Self::OutsideSpan {
                t_ps,
                start_ps,
                end_ps,
            } => write!(
                f,
                "t = {t_ps} ps is outside the trajectory span [{start_ps}, {end_ps}] ps"
            ),
            Self::Trajectory(m) => write!(f, "trajectory: {m}"),
        }
    }
}

impl std::error::Error for GeoError {}
