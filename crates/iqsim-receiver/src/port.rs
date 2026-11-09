//! Antenna-port output conversion (contracts §2 "Output conversion", §5.5).
//!
//! The antenna-port target is rendered through B_out with its own full scale V_fs and
//! the ideal output-device timebase. It bypasses every virtual-receiver effect: no
//! gain, noise, LO, IQ, decimation, AGC or ADC. Its accumulator is converted once,
//! with integers only where possible:
//!
//! - sc16: `code = round_ties_even(acc / 2²⁵)`, so code 32768 is 2⁴⁰ LSB, the target's
//!   *effective* full scale `2⁴⁰ / lsb_per_volt`
//!   (`TargetScale::full_scale_effective_volts`). That differs from the declared V_fs
//!   only by the f32 rounding of the scale factor (≤ 2⁻²⁴ relative), and the truth
//!   record reports it. A code outside [−32768, 32767] is an error, because clipping on
//!   the antenna port would corrupt the hardware drive signal silently.
//! - cf32: `f32(acc)·2⁻⁴⁰`, i.e. volts / effective full scale. The power-of-two scaling
//!   is exact, so the only rounding is the i64 → f32 conversion.

use std::fmt;

use crate::chain::{OutputFormat, Samples};

/// Accumulator LSBs per sc16 code: 2⁴⁰ / 2¹⁵.
const CODE_SHIFT: u32 = iqsim_core::amplitude::ACC_SCALE_BITS - 15;

/// An antenna-port conversion failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PortError {
    /// A component exceeded the sc16 range at this complex sample of the block.
    Clip {
        /// Sample index within the block.
        k: usize,
    },
    /// The accumulator slice has odd length.
    Length,
}

impl fmt::Display for PortError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Clip { k } => write!(
                f,
                "antenna-port output clips at block sample {k}: raise its full scale"
            ),
            Self::Length => f.write_str("antenna-port input must be interleaved I/Q pairs"),
        }
    }
}

impl std::error::Error for PortError {}

/// Converts antenna-port accumulators to output samples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PortConverter {
    format: OutputFormat,
}

/// `round_ties_even(a / 2^CODE_SHIFT)` in integers.
fn round_shift(a: i64) -> i64 {
    let q = a >> CODE_SHIFT;
    let r = a & ((1_i64 << CODE_SHIFT) - 1);
    let half = 1_i64 << (CODE_SHIFT - 1);
    if r > half || (r == half && q & 1 == 1) {
        q + 1
    } else {
        q
    }
}

impl PortConverter {
    /// A converter to `format`.
    pub const fn new(format: OutputFormat) -> Self {
        Self { format }
    }

    /// The output format.
    pub const fn format(&self) -> OutputFormat {
        self.format
    }

    /// Converts one block of I/Q-interleaved accumulators.
    pub fn convert(&self, acc: &[i64]) -> Result<Samples, PortError> {
        if !acc.len().is_multiple_of(2) {
            return Err(PortError::Length);
        }
        match self.format {
            OutputFormat::Ci16 => {
                let mut out = Vec::with_capacity(acc.len());
                for (i, &a) in acc.iter().enumerate() {
                    let c = round_shift(a);
                    let c = i16::try_from(c).map_err(|_| PortError::Clip { k: i / 2 })?;
                    out.push(c);
                }
                Ok(Samples::Ci16(out))
            }
            OutputFormat::Cf32 => {
                let scale = 1.0 / (1_u64 << iqsim_core::amplitude::ACC_SCALE_BITS) as f32;
                Ok(Samples::Cf32(
                    acc.iter().map(|&a| a as f32 * scale).collect(),
                ))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sc16_rounds_ties_to_even_and_rejects_clipping() {
        let u = 1_i64 << CODE_SHIFT;
        let h = u / 2;
        let conv = PortConverter::new(OutputFormat::Ci16);
        let acc = [0, h, 3 * h, -h, -3 * h, h + 1, u * 32_767, -u * 32_768];
        let Samples::Ci16(c) = conv.convert(&acc).unwrap() else {
            panic!()
        };
        assert_eq!(c, vec![0, 0, 2, 0, -2, 1, 32_767, -32_768]);
        assert_eq!(
            conv.convert(&[0, 0, u * 32_767 + h, 0]),
            Err(PortError::Clip { k: 1 })
        );
        assert_eq!(conv.convert(&[0]), Err(PortError::Length));
        // A value just below the tie rounds down.
        let Samples::Ci16(c) = conv
            .convert(&[u * 32_767 + h - 1, -u * 32_768 - h])
            .unwrap()
        else {
            panic!()
        };
        assert_eq!(c, vec![32_767, -32_768]);
    }

    #[test]
    fn cf32_is_volts_over_full_scale() {
        let conv = PortConverter::new(OutputFormat::Cf32);
        let one = 1_i64 << 40;
        let Samples::Cf32(c) = conv.convert(&[one, -one / 2, 3, 0]).unwrap() else {
            panic!()
        };
        assert_eq!(c, vec![1.0, -0.5, 3.0 / (1_u64 << 40) as f32, 0.0]);
    }
}
