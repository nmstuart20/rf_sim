//! Physical amplitude normalization, the per-target accumulator scale and the headroom
//! proof (contracts §2).
//!
//! # Normalization
//!
//! A sample is the complex-envelope voltage at the target's reference plane across
//! R₀ = [`REFERENCE_IMPEDANCE_OHMS`] (50 Ω), RMS convention: `P = |v|²/R₀`, so a tone of
//! power P has `|v| = √(P·R₀)` ([`rms_volts`]). Path amplitude polynomials are in volts.
//!
//! # Accumulator scale
//!
//! Each render target declares a full-scale voltage V_fs ([`TargetScale`]). Its
//! accumulator unit is `q = V_fs·2⁻⁴⁰` ([`ACC_SCALE_BITS`] = 40), and the render path
//! multiplies volts by the f32 factor [`TargetScale::lsb_per_volt`] = f32(2⁴⁰/V_fs).
//! **That f32 factor is the definition of the scale**: q is exactly its reciprocal, so
//! every backend uses the same bits. The antenna-port target and the virtual receiver
//! declare their own V_fs; neither depends on the other or on any virtual ADC setting.
//!
//! # Headroom
//!
//! Each accumulator is an i64 per I and Q. Before a block renders, the control plane
//! sums, over every path of the target,
//!
//! ```text
//! units_p = ceil(max|a_p| · max|s_p| · lsb_per_volt)
//! ```
//!
//! and requires `Σ units_p ≤ L` = [`HEADROOM_LIMIT`] = 2⁶² ([`Headroom`]). A component
//! of one contribution is at most `|a|·|s|·|rot|·factor` before rounding, and rounding
//! adds at most ½, so `ceil` bounds it. The guard bit between 2⁶² and 2⁶³ absorbs the f32
//! rounding of the per-sample products (relative ≤ 2⁻²¹) and sincos's |rot| ≤ 1 + 6·10⁻⁸.
//! The supported peak overlap is [`FULL_SCALE_PATHS`] = 2²² simultaneous full-scale
//! contributions per sample. A configuration that fails is rejected
//! ([`AmplitudeError::Headroom`]). The CPU backend also uses checked adds.
//!
//! # Culling
//!
//! A path whose peak contribution is below [`LOSSLESS_CULL_LSB`] = ¼ LSB in every
//! sample rounds to zero in every sample: the computed component is below
//! ¼·(1 + 2⁻²⁰) < ½ LSB, and ties-to-even maps it to 0. Culling at or below that
//! threshold is therefore lossless (output bits are identical with and without it).
//! A target may configure a larger threshold, which then is a documented approximation.

use std::fmt;

/// Reference impedance R₀, ohms.
pub const REFERENCE_IMPEDANCE_OHMS: f64 = 50.0;

/// Accumulator scale exponent S: one LSB is V_fs·2⁻ˢ (−240 dBFS).
pub const ACC_SCALE_BITS: u32 = 40;

/// Headroom limit L on the per-target sum of contribution bounds, in LSB.
pub const HEADROOM_LIMIT: u64 = 1 << 62;

/// Number of simultaneous full-scale contributions that fit under [`HEADROOM_LIMIT`].
pub const FULL_SCALE_PATHS: u64 = HEADROOM_LIMIT >> ACC_SCALE_BITS;

/// Default and largest lossless cull threshold, in LSB.
pub const LOSSLESS_CULL_LSB: f64 = 0.25;

/// RMS voltage of a power `p_w` watts across R₀: `√(P·R₀)`.
pub fn rms_volts(p_w: f64) -> f64 {
    (p_w * REFERENCE_IMPEDANCE_OHMS).sqrt()
}

/// A rejected scale or a failed headroom proof.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AmplitudeError {
    /// V_fs is not finite and positive, or its scale factor is not a normal f32.
    FullScale,
    /// An amplitude or envelope bound was negative or not finite.
    Bound,
    /// The sum of contribution bounds exceeds [`HEADROOM_LIMIT`].
    Headroom {
        /// Sum of contribution bounds that was needed, LSB (saturated at `u64::MAX`).
        needed: u64,
    },
}

impl fmt::Display for AmplitudeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FullScale => f.write_str(
                "full-scale voltage must be finite and positive with a normal f32 scale factor",
            ),
            Self::Bound => f.write_str("amplitude bound must be finite and non-negative"),
            Self::Headroom { needed } => write!(
                f,
                "accumulator headroom exceeded: paths need {needed} LSB, limit is 2^62 \
                 ({FULL_SCALE_PATHS} full-scale contributions)"
            ),
        }
    }
}

impl std::error::Error for AmplitudeError {}

/// A render target's integer scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TargetScale {
    v_fs: f64,
    lsb_per_volt: f32,
}

impl TargetScale {
    /// The scale for full-scale voltage `v_fs` (volts, complex-envelope RMS).
    pub fn new(v_fs: f64) -> Result<Self, AmplitudeError> {
        if !v_fs.is_finite() || v_fs <= 0.0 {
            return Err(AmplitudeError::FullScale);
        }
        // 2⁴⁰ / V_fs, one correctly rounded f64 division, then one rounding to f32.
        let factor = ((1u64 << ACC_SCALE_BITS) as f64 / v_fs) as f32;
        if !factor.is_normal() {
            return Err(AmplitudeError::FullScale);
        }
        Ok(Self {
            v_fs,
            lsb_per_volt: factor,
        })
    }

    /// Declared full-scale voltage.
    pub const fn full_scale_volts(&self) -> f64 {
        self.v_fs
    }

    /// The f32 volts → LSB factor used by every backend.
    pub const fn lsb_per_volt(&self) -> f32 {
        self.lsb_per_volt
    }

    /// The effective full scale `2⁴⁰ / lsb_per_volt`: the voltage that is exactly 2⁴⁰
    /// accumulator LSB. It differs from the declared V_fs only by the f32 rounding of
    /// the scale factor (≤ 2⁻²⁴ relative). Conversions that work on the integer
    /// accumulator (the antenna port's sc16 code 32768, cf32 1.0) are relative to it.
    pub fn full_scale_effective_volts(&self) -> f64 {
        (1_u64 << ACC_SCALE_BITS) as f64 / f64::from(self.lsb_per_volt)
    }

    /// One LSB in volts: exactly `1 / lsb_per_volt`.
    pub fn lsb_volts(&self) -> f64 {
        1.0 / f64::from(self.lsb_per_volt)
    }

    /// Peak contribution of a path, in LSB, for amplitude bound `max_amp_volts` and
    /// envelope bound `max_env`.
    pub fn peak_lsb(&self, max_amp_volts: f64, max_env: f64) -> Result<f64, AmplitudeError> {
        if !(max_amp_volts.is_finite() && max_env.is_finite())
            || max_amp_volts < 0.0
            || max_env < 0.0
        {
            return Err(AmplitudeError::Bound);
        }
        Ok(max_amp_volts * max_env * f64::from(self.lsb_per_volt))
    }
}

/// Running headroom sum for one target and block.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct Headroom {
    used: u64,
}

impl Headroom {
    /// An empty sum.
    pub const fn new() -> Self {
        Self { used: 0 }
    }

    /// Adds a path whose peak contribution is `peak_lsb` (from
    /// [`TargetScale::peak_lsb`]). Fails once the sum exceeds [`HEADROOM_LIMIT`].
    pub fn add(&mut self, peak_lsb: f64) -> Result<(), AmplitudeError> {
        if !peak_lsb.is_finite() || peak_lsb < 0.0 {
            return Err(AmplitudeError::Bound);
        }
        let units = peak_lsb.ceil();
        // 2⁶² is exact in f64; anything at or above 2⁶⁴ saturates the report.
        let needed_f = units + self.used as f64;
        let units_u = if units >= 18_446_744_073_709_551_616.0 {
            u64::MAX
        } else {
            units as u64
        };
        let needed = self.used.saturating_add(units_u);
        if needed > HEADROOM_LIMIT || needed_f > HEADROOM_LIMIT as f64 {
            return Err(AmplitudeError::Headroom { needed });
        }
        self.used = needed;
        Ok(())
    }

    /// Sum so far, LSB.
    pub const fn used(&self) -> u64 {
        self.used
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scale_is_defined_by_the_f32_factor() {
        let s = TargetScale::new(1.0).unwrap();
        assert_eq!(s.lsb_per_volt(), 1_099_511_627_776.0_f32);
        assert_eq!(s.lsb_volts(), 1.0 / 1_099_511_627_776.0);
        assert_eq!(FULL_SCALE_PATHS, 1 << 22);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-40] {
            assert_eq!(
                TargetScale::new(bad),
                Err(AmplitudeError::FullScale),
                "{bad}"
            );
        }
        assert!((rms_volts(1e-3) - 0.223_606_797_749_979).abs() < 1e-15);
    }

    #[test]
    fn headroom_accepts_up_to_the_limit_and_rejects_beyond() {
        let s = TargetScale::new(1.0).unwrap();
        let full = s.peak_lsb(1.0, 1.0).unwrap();
        let mut h = Headroom::new();
        for _ in 0..FULL_SCALE_PATHS {
            h.add(full).unwrap();
        }
        assert_eq!(h.used(), HEADROOM_LIMIT);
        assert!(matches!(h.add(0.5), Err(AmplitudeError::Headroom { .. })));
        let mut h = Headroom::new();
        assert!(matches!(
            h.add(1e30),
            Err(AmplitudeError::Headroom { needed: u64::MAX })
        ));
        assert_eq!(Headroom::new().add(f64::NAN), Err(AmplitudeError::Bound));
        assert_eq!(s.peak_lsb(-1.0, 1.0), Err(AmplitudeError::Bound));
    }
}
