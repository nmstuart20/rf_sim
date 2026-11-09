//! Gillham (Mode C) pressure-altitude encoding in 100 ft steps.
//!
//! For an exact altitude `alt`, the 500 ft band index is
//! `n500 = (alt + 1200) / 500`. Its reflected binary Gray code,
//! `n500 ^ (n500 >> 1)`, occupies D2 D4 A1 A2 A4 B1 B2 B4, most
//! significant bit first. Within the band, `m = ((alt + 1200) / 100) % 5`
//! selects C1 C2 C4 from `001, 011, 010, 110, 100`; odd bands reverse
//! that sequence (`m` becomes `4 - m`). C patterns `000`, `101`, `111`
//! and any code with D1 set are illegal. The 256 bands and five legal
//! C patterns cover exactly -1200 through 126700 ft, inclusive.
//!
//! Decoding reverses the Gray transform and the odd-band C reversal.
//! The mapping is checked against all 1280 rows of the AirSport Avionics
//! table by Darryl Phillips (1998/2012), checked in as
//! `tests/data/airsport_modecascii.txt`.
//!
//! [`GillhamCode::encode_pressure_altitude_ft`] quantizes pressure altitude
//! to the nearest 100 ft, with ties toward positive infinity: 949.9 ft
//! reports 900 ft, 950.0 ft reports 1000 ft, and -1250.0 ft reports -1200 ft.
//! This encoder quantization is distinct from numeric ties-to-even conversion.

use std::fmt;

/// Lowest pressure altitude with a Gillham code, feet.
pub const GILLHAM_MIN_FT: i32 = -1200;

/// Highest pressure altitude with a Gillham code, feet.
pub const GILLHAM_MAX_FT: i32 = 126_700;

const C_PATTERNS: [u16; 5] = [0b001, 0b011, 0b010, 0b110, 0b100];

/// A validated 12-bit Gillham (Mode C) code. D1 is always zero.
///
/// Internally, bits 11 through 0 are D1 D2 D4 A1 A2 A4 B1 B2 B4 C1 C2 C4,
/// the same layout returned by [`Self::wire_order`]. Bits above 11 are zero.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GillhamCode(u16);

impl GillhamCode {
    /// Encodes an exact multiple of 100 ft in the inclusive Gillham range.
    ///
    /// Returns [`GillhamError::OutOfRange`] outside the range, or
    /// [`GillhamError::NotMultipleOf100`] for an in-range nonmultiple.
    pub fn from_altitude_ft(alt_100ft_multiple: i32) -> Result<Self, GillhamError> {
        if !(GILLHAM_MIN_FT..=GILLHAM_MAX_FT).contains(&alt_100ft_multiple) {
            return Err(GillhamError::OutOfRange {
                altitude_ft: f64::from(alt_100ft_multiple),
            });
        }
        if alt_100ft_multiple % 100 != 0 {
            return Err(GillhamError::NotMultipleOf100 {
                altitude_ft: alt_100ft_multiple,
            });
        }
        // Range validation bounds the offset to 0..=127900, so all arithmetic fits.
        let offset = alt_100ft_multiple - GILLHAM_MIN_FT;
        let n500 = (offset / 500) as u16;
        let m = (offset / 100 % 5) as usize;
        let c_index = if n500 & 1 == 0 { m } else { 4 - m };
        let gray = n500 ^ (n500 >> 1);
        Ok(Self((gray << 3) | C_PATTERNS[c_index]))
    }

    /// Decodes this validated code to its exact pressure altitude, feet.
    pub fn altitude_ft(self) -> i32 {
        let mut n500 = self.0 >> 3;
        // Invert the eight-bit reflected Gray code by prefix XOR.
        n500 ^= n500 >> 1;
        n500 ^= n500 >> 2;
        n500 ^= n500 >> 4;
        let c_index = match self.0 & 0b111 {
            0b001 => 0,
            0b011 => 1,
            0b010 => 2,
            0b110 => 3,
            0b100 => 4,
            _ => unreachable!("GillhamCode constructors validate the C pattern"),
        };
        let m = if n500 & 1 == 0 { c_index } else { 4 - c_index };
        i32::from(n500) * 500 + m * 100 + GILLHAM_MIN_FT
    }

    /// Encodes pressure altitude rounded to the nearest 100 ft, ties toward +∞.
    ///
    /// Thus 950.0 ft reports 1000 ft and -1250.0 ft reports -1200 ft.
    /// Returns [`GillhamError::NonFinite`] for NaN or infinity, and
    /// [`GillhamError::OutOfRange`] when the rounded altitude is outside
    /// [`GILLHAM_MIN_FT`] through [`GILLHAM_MAX_FT`].
    pub fn encode_pressure_altitude_ft(pressure_alt_ft: f64) -> Result<Self, GillhamError> {
        if !pressure_alt_ft.is_finite() {
            return Err(GillhamError::NonFinite);
        }
        let lower = (pressure_alt_ft / 100.0).floor() * 100.0;
        // Compare the remainder instead of adding 50 first: that addition can
        // round a value immediately below a transition up to the transition.
        let rounded = if pressure_alt_ft - lower >= 50.0 {
            lower + 100.0
        } else {
            lower
        };
        if !(f64::from(GILLHAM_MIN_FT)..=f64::from(GILLHAM_MAX_FT)).contains(&rounded) {
            return Err(GillhamError::OutOfRange {
                altitude_ft: pressure_alt_ft,
            });
        }
        // The rounded value is an integral multiple of 100 within the i32 range.
        Self::from_altitude_ft(rounded as i32)
    }

    /// Validates named bits and constructs a code.
    ///
    /// D1 must be zero and C1 C2 C4 must be one of `001, 011, 010, 110, 100`.
    /// Every combination of the other eight bits and a legal C pattern decodes
    /// within the Gillham altitude range.
    pub fn from_bits(bits: GillhamBits) -> Result<Self, GillhamError> {
        let wire = (u16::from(bits.d1) << 11)
            | (u16::from(bits.d2) << 10)
            | (u16::from(bits.d4) << 9)
            | (u16::from(bits.a1) << 8)
            | (u16::from(bits.a2) << 7)
            | (u16::from(bits.a4) << 6)
            | (u16::from(bits.b1) << 5)
            | (u16::from(bits.b2) << 4)
            | (u16::from(bits.b4) << 3)
            | (u16::from(bits.c1) << 2)
            | (u16::from(bits.c2) << 1)
            | u16::from(bits.c4);
        Self::from_wire_order(wire)
    }

    /// Returns the twelve named bits.
    pub fn bits(self) -> GillhamBits {
        GillhamBits {
            a1: self.0 & (1 << 8) != 0,
            a2: self.0 & (1 << 7) != 0,
            a4: self.0 & (1 << 6) != 0,
            b1: self.0 & (1 << 5) != 0,
            b2: self.0 & (1 << 4) != 0,
            b4: self.0 & (1 << 3) != 0,
            c1: self.0 & (1 << 2) != 0,
            c2: self.0 & (1 << 1) != 0,
            c4: self.0 & 1 != 0,
            d1: false,
            d2: self.0 & (1 << 10) != 0,
            d4: self.0 & (1 << 9) != 0,
        }
    }

    /// Returns the octal ABCD value as an integer, for example `0o0620` at 0 ft.
    ///
    /// Each digit is `4 * bit4 + 2 * bit2 + bit1` for its letter group.
    pub fn squawk_octal(self) -> u16 {
        let bits = self.bits();
        let a = 4 * u16::from(bits.a4) + 2 * u16::from(bits.a2) + u16::from(bits.a1);
        let b = 4 * u16::from(bits.b4) + 2 * u16::from(bits.b2) + u16::from(bits.b1);
        let c = 4 * u16::from(bits.c4) + 2 * u16::from(bits.c2) + u16::from(bits.c1);
        let d = 4 * u16::from(bits.d4) + 2 * u16::from(bits.d2);
        (a << 9) | (b << 6) | (c << 3) | d
    }

    /// Returns four octal ABCD digits, including leading zeros, e.g. `"0620"`.
    pub fn squawk_string(self) -> String {
        format!("{:04o}", self.squawk_octal())
    }

    /// Returns bits 11..0 as D1 D2 D4 A1 A2 A4 B1 B2 B4 C1 C2 C4.
    pub fn wire_order(self) -> u16 {
        self.0
    }

    /// Validates a wire-order value and constructs a code.
    ///
    /// Rejects bits above bit 11, D1 set, and illegal C patterns, in that order.
    /// The remaining 1280 combinations decode within the Gillham range.
    pub fn from_wire_order(bits: u16) -> Result<Self, GillhamError> {
        if bits & !0x0fff != 0 {
            return Err(GillhamError::WireBitsOutOfRange { bits });
        }
        if bits & (1 << 11) != 0 {
            return Err(GillhamError::D1Set);
        }
        let c_bits = bits & 0b111;
        if !C_PATTERNS.contains(&c_bits) {
            return Err(GillhamError::IllegalC { c_bits });
        }
        Ok(Self(bits))
    }
}

/// Named Gillham bits; [`GillhamCode::from_bits`] validates their combination.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GillhamBits {
    /// A group's weight-1 bit.
    pub a1: bool,
    /// A group's weight-2 bit.
    pub a2: bool,
    /// A group's weight-4 bit.
    pub a4: bool,
    /// B group's weight-1 bit.
    pub b1: bool,
    /// B group's weight-2 bit.
    pub b2: bool,
    /// B group's weight-4 bit.
    pub b4: bool,
    /// C group's weight-1 bit.
    pub c1: bool,
    /// C group's weight-2 bit.
    pub c2: bool,
    /// C group's weight-4 bit.
    pub c4: bool,
    /// D group's weight-1 bit; always zero in a valid Mode C code.
    pub d1: bool,
    /// D group's weight-2 bit.
    pub d2: bool,
    /// D group's weight-4 bit.
    pub d4: bool,
}

/// An invalid altitude or Gillham bit combination.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GillhamError {
    /// The exact or rounded altitude is outside the inclusive Gillham range.
    OutOfRange {
        /// Input altitude, feet, before any encoder rounding.
        altitude_ft: f64,
    },
    /// An exact-altitude input was not a multiple of 100 ft.
    NotMultipleOf100 {
        /// Input altitude, feet.
        altitude_ft: i32,
    },
    /// Pressure altitude was NaN or infinite.
    NonFinite,
    /// C1 C2 C4 was `000`, `101` or `111`.
    IllegalC {
        /// C1 C2 C4 in bits 2..0.
        c_bits: u16,
    },
    /// D1 was set, which is forbidden for Mode C altitude.
    D1Set,
    /// A wire-order value had bits set above bit 11.
    WireBitsOutOfRange {
        /// The rejected wire-order value.
        bits: u16,
    },
}

impl fmt::Display for GillhamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OutOfRange { altitude_ft } => write!(
                f,
                "altitude {altitude_ft} ft is outside the Gillham range {GILLHAM_MIN_FT}..={GILLHAM_MAX_FT} ft after rounding"
            ),
            Self::NotMultipleOf100 { altitude_ft } => {
                write!(f, "altitude {altitude_ft} ft is not a multiple of 100 ft")
            }
            Self::NonFinite => f.write_str("pressure altitude must be finite"),
            Self::IllegalC { c_bits } => write!(f, "illegal Gillham C1 C2 C4 pattern {c_bits:03b}"),
            Self::D1Set => f.write_str("Gillham D1 must be zero"),
            Self::WireBitsOutOfRange { bits } => {
                write!(
                    f,
                    "Gillham wire value {bits:#06x} has bits above bit 11 set"
                )
            }
        }
    }
}

impl std::error::Error for GillhamError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hand_checked_reference_rows() {
        for (altitude, squawk, wire) in [
            (-1200, 0o0040, 0x001),
            (-1000, 0o0020, 0x002),
            (0, 0o0620, 0x01a),
            (30_700, 0o1040, 0x101),
            (30_800, 0o1044, 0x301),
            (62_700, 0o0044, 0x201),
            (62_800, 0o0046, 0x601),
            (126_700, 0o0042, 0x401),
        ] {
            let code = GillhamCode::from_altitude_ft(altitude).unwrap();
            assert_eq!(code.squawk_octal(), squawk, "{altitude} ft");
            assert_eq!(code.wire_order(), wire, "{altitude} ft");
            assert_eq!(
                GillhamCode::from_wire_order(wire).unwrap().altitude_ft(),
                altitude
            );
        }
    }
}
