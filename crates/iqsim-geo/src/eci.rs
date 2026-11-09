//! Earth-centred inertial frame and Earth rotation (contracts §12.1).
//!
//! The simulator's ECI frame is **pseudo-inertial**: it coincides with ECEF at the
//! scenario epoch (truth time 0), and ECEF then spins about the shared z axis at the
//! constant WGS84 rate ω_E = 7 292 115·10⁻¹¹ rad/s. Precession, nutation, polar motion
//! and UT1 are below the design's fidelity (plan, "World model": Frames), so this is
//! not GCRF; it is the frame the light-time solver works in.
//!
//! ```text
//! r_eci = R_z(θ)·r_ecef,   v_eci = R_z(θ)·(v_ecef + ω × r_ecef),   θ = ω_E·t
//! R_z(θ) = [[cos θ, −sin θ, 0], [sin θ, cos θ, 0], [0, 0, 1]]
//! ```
//!
//! **Angle reduction in integers.** θ grows to ~670 rad over the 106-day truth-time
//! range, where one f64 ulp is 1.1·10⁻¹³ rad (5 µm at GEO radius). Instead, for an
//! [`AcqInstant`] `t = ps + frac·2⁻⁶⁴` ps:
//!
//! 1. `S = ps·7 292 115·2⁴¹ + ⌊frac·7 292 115 / 2²³⌋`, exact in i128, in units of
//!    10⁻²³·2⁻⁴¹ rad (`7 292 115·10⁻¹¹ rad/s × 10⁻¹² s/ps = 7 292 115·10⁻²³ rad/ps`).
//!    `|S| ≤ 1.48·10³⁸ < 2¹²⁷` over the whole i64 ps range.
//! 2. `r = S mod T` with `T = round(2π·10²³·2⁴¹)` (error 0.05 units), so `r ∈ [0, T)`.
//! 3. `F = ⌊r·2⁶⁴/T⌋`, the angle as a 64-bit fraction of a turn, by exact long
//!    division (eleven 6-bit steps; every intermediate stays below 2¹²⁷).
//! 4. `θ = F·2⁻⁶⁴·2π` with 2π as a double-double `(τ_hi, τ_lo)`, F split into its top
//!    53 bits `F_h` and low 11 bits `F_l` (both exact in f64), and Dekker's exact
//!    two-product (no FMA): `(p, e) = F_h ⊗ τ_hi`, then
//!    `θ = (p + (e + (F_l·τ_hi + F_h·τ_lo)))·2⁻⁶⁴`.
//!
//! θ is in [0, 2π] and within 0.5 ulp + 3.5·10⁻¹⁹ rad of the exact angle at any time
//! in the truth-time range: the floor in step 3 (≤ 2⁻⁶⁴ turn = 3.4·10⁻¹⁹ rad, always
//! downward) dominates the other reduction errors (T, the floor in step 1: < 10⁻³⁰ rad).
//! Near 2π that is ≤ 0.5004 ulp, i.e. 4.5·10⁻¹⁶ rad (1.9·10⁻⁸ m at GEO radius,
//! 1.7·10⁻⁷ m at lunar distance). A naive f64 `ω·t` has a 1.1·10⁻¹³ rad ulp at 100 days.

use iqsim_core::{AcqInstant, TruthTime};

use crate::geodetic::Ecef;
use crate::vec::Vec3;

/// ω_E, rad/s (WGS84: 7 292 115·10⁻¹¹ exactly; this f64 is its nearest value).
pub const EARTH_ROTATION_RATE_RAD_S: f64 = 7.292_115e-5;

/// ω_E numerator over 10⁻¹¹ rad/s.
const OMEGA_E11: i128 = 7_292_115;

/// Fraction bits of the fixed-point angle.
const ANGLE_FRAC_BITS: u32 = 41;

/// `round(2π·10²³·2⁴¹)`: one turn in units of 10⁻²³·2⁻⁴¹ rad (mpmath, 80 digits).
const TURN: i128 = 1_381_687_060_943_054_741_360_535_823_246_976_831;

/// 2π as a double-double: the nearest f64 and the f64 nearest the remainder
/// (`0x1.921fb54442d18p+2`, `0x1.1a62633145c07p-52`, mpmath).
const TAU_HI: f64 = core::f64::consts::TAU;
const TAU_LO: f64 = 2.449_293_598_294_706_4e-16;

/// `2⁻⁶⁴`.
const TWO_M64: f64 = 1.0 / 18_446_744_073_709_551_616.0;

/// Dekker's split constant `2²⁷ + 1`.
const SPLIT: f64 = 134_217_729.0;

/// `a = hi + lo` with each half 26 bits (Dekker), separate operations only.
fn split(a: f64) -> (f64, f64) {
    let c = SPLIT * a;
    let hi = c - (c - a);
    (hi, a - hi)
}

/// `a·b = p + e` exactly (Dekker's two-product, no FMA).
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    (p, e)
}

/// `⌊r·2⁶⁴/TURN⌋` for `0 ≤ r < TURN`, by long division in 6-bit steps.
fn turn_fraction(r: i128) -> u64 {
    let mut rem = r;
    let mut q: u128 = 0;
    for _ in 0..11 {
        rem <<= 6;
        let d = rem / TURN;
        rem -= d * TURN;
        q = (q << 6) | d as u128;
    }
    // 66 quotient bits; keep the top 64 (floor).
    (q >> 2) as u64
}

/// A position in the pseudo-inertial ECI frame (module docs), metres.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Eci {
    /// x, metres.
    pub x: f64,
    /// y, metres.
    pub y: f64,
    /// z, metres.
    pub z: f64,
}

impl Eci {
    /// The position (x, y, z) metres.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// The position from a vector.
    pub const fn from_vec(v: Vec3) -> Self {
        Self::new(v.x, v.y, v.z)
    }

    /// The position as a vector.
    pub const fn vec(self) -> Vec3 {
        Vec3::new(self.x, self.y, self.z)
    }
}

/// The Earth rotation angle θ = ω_E·t at `t`, reduced to [0, 2π) (module docs).
pub fn earth_rotation_angle(t: AcqInstant) -> f64 {
    let whole = i128::from(t.ps().as_ps()) * OMEGA_E11 * (1_i128 << ANGLE_FRAC_BITS);
    let part = (i128::from(t.frac()) * OMEGA_E11) >> 23;
    let s = whole + part;
    let r = s.rem_euclid(TURN);
    let f = turn_fraction(r);
    let fh = (f & !0x7ff) as f64;
    let fl = (f & 0x7ff) as f64;
    let (p, e) = two_prod(fh, TAU_HI);
    let t = fl * TAU_HI + fh * TAU_LO;
    (p + (e + t)) * TWO_M64
}

/// The ECEF ↔ ECI rotation at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EarthOrientation {
    theta: f64,
    sin: f64,
    cos: f64,
}

impl EarthOrientation {
    /// The orientation at acquisition instant `t`.
    pub fn at(t: AcqInstant) -> Self {
        Self::from_angle(earth_rotation_angle(t))
    }

    /// The orientation at truth time `t`.
    pub fn at_time(t: TruthTime) -> Self {
        Self::at(AcqInstant::from_truth_time(t))
    }

    /// The orientation for a given rotation angle θ, radians.
    pub fn from_angle(theta: f64) -> Self {
        let (sin, cos) = libm::sincos(theta);
        Self { theta, sin, cos }
    }

    /// θ, radians.
    pub const fn angle(&self) -> f64 {
        self.theta
    }

    fn rz(&self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.cos * v.x - self.sin * v.y,
            self.sin * v.x + self.cos * v.y,
            v.z,
        )
    }

    fn rz_t(&self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.cos * v.x + self.sin * v.y,
            self.cos * v.y - self.sin * v.x,
            v.z,
        )
    }

    /// `r_eci = R_z(θ)·r_ecef`.
    pub fn ecef_to_eci(&self, p: Ecef) -> Eci {
        Eci::from_vec(self.rz(p.vec()))
    }

    /// `r_ecef = R_z(θ)ᵀ·r_eci`.
    pub fn eci_to_ecef(&self, p: Eci) -> Ecef {
        Ecef::from_vec(self.rz_t(p.vec()))
    }

    /// A free vector (direction or displacement, not a velocity) from ECEF to ECI.
    pub fn vector_to_eci(&self, v: Vec3) -> Vec3 {
        self.rz(v)
    }

    /// A free vector from ECI to ECEF.
    pub fn vector_to_ecef(&self, v: Vec3) -> Vec3 {
        self.rz_t(v)
    }

    /// Position and velocity, ECEF → ECI: `v_eci = R_z(θ)·(v_ecef + ω × r_ecef)`.
    pub fn state_to_eci(&self, p: Ecef, v_ecef: Vec3) -> (Eci, Vec3) {
        let w = Vec3::new(0.0, 0.0, EARTH_ROTATION_RATE_RAD_S);
        (self.ecef_to_eci(p), self.rz(v_ecef + w.cross(p.vec())))
    }

    /// Position and velocity, ECI → ECEF: `v_ecef = R_z(θ)ᵀ·v_eci − ω × r_ecef`.
    pub fn state_to_ecef(&self, p: Eci, v_eci: Vec3) -> (Ecef, Vec3) {
        let w = Vec3::new(0.0, 0.0, EARTH_ROTATION_RATE_RAD_S);
        let r = self.eci_to_ecef(p);
        (r, self.rz_t(v_eci) - w.cross(r.vec()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants() {
        assert_eq!(TAU_HI.to_bits(), 0x4019_21fb_5444_2d18);
        // mpmath: 2π − TAU_HI rounds to 0x1.1a62633145c07p-52.
        assert_eq!(TAU_LO.to_bits(), 0x3cb1_a626_3314_5c07);
        // Two-product is exact on a case f64 multiplication rounds.
        let (p, e) = two_prod(1.0 + f64::EPSILON, 1.0 + f64::EPSILON);
        assert_eq!(
            (p, e),
            (1.0 + 2.0 * f64::EPSILON, f64::EPSILON * f64::EPSILON)
        );
        // Turn fractions: 0, ½ and just below a full turn.
        assert_eq!(turn_fraction(0), 0);
        // TURN is odd: TURN/2 is just under half a turn, (TURN + 1)/2 just over.
        assert_eq!(turn_fraction(TURN / 2), (1 << 63) - 1);
        assert_eq!(turn_fraction((TURN + 1) / 2), 1 << 63);
        assert_eq!(turn_fraction(TURN - 1), u64::MAX);
        assert_eq!(
            earth_rotation_angle(AcqInstant::from_truth_time(TruthTime::EPOCH)),
            0.0
        );
    }

    #[test]
    fn ecef_point_moves_eastward_in_eci() {
        // A quarter of a sidereal turn later, a point on the x axis is on the y axis.
        let quarter_s = core::f64::consts::FRAC_PI_2 / EARTH_ROTATION_RATE_RAD_S;
        let t = TruthTime::from_ps((quarter_s * 1e12) as i64);
        let o = EarthOrientation::at_time(t);
        let p = o.ecef_to_eci(Ecef::new(7.0e6, 0.0, 1.0e6));
        assert!((p.x).abs() < 1e-3 && (p.y - 7.0e6).abs() < 1e-3 && p.z == 1.0e6);
        // Round trip.
        let q = Ecef::new(-1.2e6, 6.1e6, 3.3e5);
        assert!(o.eci_to_ecef(o.ecef_to_eci(q)).distance(q) < 1e-8);
    }

    #[test]
    fn fixed_ecef_point_has_omega_cross_r_velocity() {
        let o = EarthOrientation::at_time(TruthTime::from_ps(123_456_789_000_000));
        let r = Ecef::new(6.0e6, -2.0e6, 1.0e6);
        let (p, v) = o.state_to_eci(r, Vec3::ZERO);
        let w = Vec3::new(0.0, 0.0, EARTH_ROTATION_RATE_RAD_S);
        let want = w.cross(p.vec());
        assert!(v.distance(want) < 1e-9);
        let (r2, v2) = o.state_to_ecef(p, v);
        assert!(r2.distance(r) < 1e-8 && v2.norm() < 1e-9);
    }
}
