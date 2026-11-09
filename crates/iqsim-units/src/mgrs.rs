//! Military Grid Reference System (MGRS) references on WGS84.
//!
//! A reference such as `18SUJ2337106519` names a UTM zone (18), a latitude band (S),
//! a 100 km square (column U, row J) and an easting/northing inside it at 1 m to
//! 100 km resolution (0–5 digits each). WGS84 uses the "AA" lettering scheme: the
//! column letters cycle through three sets of eight by zone, and the row letters
//! cycle every 2000 km, starting at F in even zones.
//!
//! A reference denotes a grid **square**, not a point. [`Mgrs::sw_corner`] gives the
//! square's south-west corner, which is what scenario positions use (the NGA reading
//! convention), and [`Mgrs::center`] its centre (GeographicLib's representative
//! point).
//!
//! The polar UPS regions (bands A, B, Y, Z) are not supported.

use std::fmt;

use crate::utm::{FALSE_NORTHING_SOUTH, Hemisphere, utm_to_geodetic};

const BANDS: &[u8; 20] = b"CDEFGHJKLMNPQRSTUVWX";
const COLUMNS: &[u8; 24] = b"ABCDEFGHJKLMNPQRSTUVWXYZ";
const ROWS: &[u8; 20] = b"ABCDEFGHJKLMNPQRSTUV";
const SQUARE_M: u32 = 100_000;
const ROW_CYCLE_M: u32 = 2_000_000;
/// Slack in degrees when checking that a square meets its latitude band.
const BAND_SLACK_DEG: f64 = 0.01;

/// What is wrong with an MGRS reference.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MgrsErrorKind {
    /// The reference does not start with a zone number.
    MissingZone,
    /// The zone number is outside 1..=60.
    ZoneOutOfRange(u32),
    /// A UPS (polar) band letter.
    Polar(char),
    /// Not a latitude band letter.
    BadBand(char),
    /// One of the non-existent Svalbard zones 32X, 34X, 36X.
    NonexistentZone(u8),
    /// The reference ended before both 100 km square letters.
    MissingSquare,
    /// The column letter is not valid for this zone.
    BadColumn(char),
    /// The row letter is not a valid row letter.
    BadRow(char),
    /// Trailing characters that are not digits.
    BadDigits,
    /// An odd number of digits, or more than ten.
    DigitCount(usize),
    /// The 100 km square does not meet the latitude band.
    NotInBand,
}

/// A rejected MGRS reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MgrsError {
    /// The reference as written.
    pub text: String,
    /// What is wrong.
    pub kind: MgrsErrorKind,
}

impl fmt::Display for MgrsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid MGRS reference \"{}\": ", self.text)?;
        match &self.kind {
            MgrsErrorKind::MissingZone => {
                write!(
                    f,
                    "it must start with a UTM zone number 1–60, e.g. \"18SUJ2337106519\""
                )
            }
            MgrsErrorKind::ZoneOutOfRange(z) => write!(f, "zone {z} is outside 1–60"),
            MgrsErrorKind::Polar(c) => write!(
                f,
                "band {c} is a polar (UPS) band, which is not supported; give latitude and longitude instead"
            ),
            MgrsErrorKind::BadBand(c) => write!(
                f,
                "'{c}' is not a latitude band letter (C–X, without I and O)"
            ),
            MgrsErrorKind::NonexistentZone(z) => {
                write!(
                    f,
                    "zone {z}X does not exist (Svalbard uses 31X, 33X, 35X and 37X)"
                )
            }
            MgrsErrorKind::MissingSquare => {
                write!(
                    f,
                    "the two 100 km square letters are missing after the band"
                )
            }
            MgrsErrorKind::BadColumn(c) => {
                write!(f, "'{c}' is not a valid 100 km column letter for this zone")
            }
            MgrsErrorKind::BadRow(c) => write!(
                f,
                "'{c}' is not a valid 100 km row letter (A–V, without I and O)"
            ),
            MgrsErrorKind::BadDigits => {
                write!(f, "only digits may follow the 100 km square letters")
            }
            MgrsErrorKind::DigitCount(n) => write!(
                f,
                "{n} digits after the square letters; expected an even count up to 10 (easting then northing)"
            ),
            MgrsErrorKind::NotInBand => {
                write!(f, "the 100 km square does not lie in the latitude band")
            }
        }
    }
}

impl std::error::Error for MgrsError {}

/// A parsed, validated MGRS reference (non-polar).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Mgrs {
    zone: u8,
    band: char,
    hemisphere: Hemisphere,
    /// UTM easting of the square's west edge, metres.
    easting: u32,
    /// UTM northing (with false northing) of the square's south edge, metres.
    northing: u32,
    /// Digits per axis, 0..=5.
    precision: u8,
    text: String,
}

fn band_limits(band_index: usize) -> (f64, f64) {
    let lo = -80.0 + 8.0 * band_index as f64;
    let hi = if band_index == BANDS.len() - 1 {
        84.0
    } else {
        lo + 8.0
    };
    (lo, hi)
}

impl Mgrs {
    /// Parses a reference. Whitespace is ignored and letters may be lower case.
    pub fn parse(text: &str) -> Result<Self, MgrsError> {
        let fail = |kind| MgrsError {
            text: text.to_owned(),
            kind,
        };
        let compact: String = text
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .collect();
        let b = compact.as_bytes();

        let zone_len = b.iter().take(2).take_while(|c| c.is_ascii_digit()).count();
        if zone_len == 0 {
            // UPS references have no zone: they start with the band letter.
            return Err(fail(match b.first() {
                Some(&c @ (b'A' | b'B' | b'Y' | b'Z')) => MgrsErrorKind::Polar(char::from(c)),
                _ => MgrsErrorKind::MissingZone,
            }));
        }
        let zone: u32 = compact[..zone_len].parse().expect("ASCII digits");
        if !(1..=60).contains(&zone) {
            return Err(fail(MgrsErrorKind::ZoneOutOfRange(zone)));
        }
        let zone = u8::try_from(zone).expect("zone <= 60");

        let Some(&band) = b.get(zone_len) else {
            return Err(fail(MgrsErrorKind::MissingSquare));
        };
        if matches!(band, b'A' | b'B' | b'Y' | b'Z') {
            return Err(fail(MgrsErrorKind::Polar(char::from(band))));
        }
        let Some(band_index) = BANDS.iter().position(|&c| c == band) else {
            return Err(fail(MgrsErrorKind::BadBand(display_char(
                &compact, zone_len,
            ))));
        };
        if band == b'X' && matches!(zone, 32 | 34 | 36) {
            return Err(fail(MgrsErrorKind::NonexistentZone(zone)));
        }

        let (Some(&column), Some(&row)) = (b.get(zone_len + 1), b.get(zone_len + 2)) else {
            return Err(fail(MgrsErrorKind::MissingSquare));
        };
        let set = usize::from((zone - 1) % 3);
        let column_index = COLUMNS[8 * set..8 * set + 8]
            .iter()
            .position(|&c| c == column)
            .ok_or_else(|| {
                fail(MgrsErrorKind::BadColumn(display_char(
                    &compact,
                    zone_len + 1,
                )))
            })?;
        let row_letter = ROWS
            .iter()
            .position(|&c| c == row)
            .ok_or_else(|| fail(MgrsErrorKind::BadRow(display_char(&compact, zone_len + 2))))?;
        let row_offset = if zone % 2 == 0 { 5 } else { 0 };
        let row_index = (row_letter + ROWS.len() - row_offset) % ROWS.len();

        let digits = &compact[zone_len + 3..];
        if !digits.bytes().all(|c| c.is_ascii_digit()) {
            return Err(fail(MgrsErrorKind::BadDigits));
        }
        if !digits.len().is_multiple_of(2) || digits.len() > 10 {
            return Err(fail(MgrsErrorKind::DigitCount(digits.len())));
        }
        let precision = digits.len() / 2;
        let unit = 10u32.pow(5 - u32::try_from(precision).expect("<= 5"));
        let axis = |s: &str| -> u32 {
            if s.is_empty() {
                0
            } else {
                s.parse::<u32>().expect("ASCII digits") * unit
            }
        };
        let (e_digits, n_digits) = digits.split_at(precision);

        let easting = (u32::try_from(column_index).expect("< 8") + 1) * SQUARE_M + axis(e_digits);
        let hemisphere = if band < b'N' {
            Hemisphere::South
        } else {
            Hemisphere::North
        };

        // The row letters repeat every 2000 km; the band picks the cycle. The test uses
        // the whole 100 km block, so it does not depend on the digit precision (as in
        // GeographicLib). A southern block must lie south of the equator.
        let (lo, hi) = band_limits(band_index);
        let block_easting = easting - axis(e_digits);
        let block_northing = (0..=5u32)
            .map(|k| u32::try_from(row_index).expect("< 20") * SQUARE_M + k * ROW_CYCLE_M)
            .find(|&n| {
                if hemisphere == Hemisphere::South && f64::from(n + SQUARE_M) > FALSE_NORTHING_SOUTH
                {
                    return false;
                }
                let (min, max) =
                    square_latitude_range(zone, hemisphere, block_easting, n, SQUARE_M);
                min < hi + BAND_SLACK_DEG && max > lo - BAND_SLACK_DEG
            })
            .ok_or_else(|| fail(MgrsErrorKind::NotInBand))?;
        let northing = block_northing + axis(n_digits);

        Ok(Self {
            zone,
            band: char::from(band),
            hemisphere,
            easting,
            northing,
            precision: u8::try_from(precision).expect("<= 5"),
            text: text.to_owned(),
        })
    }

    /// UTM zone, 1..=60.
    pub fn zone(&self) -> u8 {
        self.zone
    }

    /// Latitude band letter.
    pub fn band(&self) -> char {
        self.band
    }

    /// UTM hemisphere implied by the band.
    pub fn hemisphere(&self) -> Hemisphere {
        self.hemisphere
    }

    /// Digits per axis, 0 (100 km square) to 5 (1 m square).
    pub fn precision(&self) -> u8 {
        self.precision
    }

    /// Side of the grid square in metres.
    pub fn square_size_m(&self) -> u32 {
        10u32.pow(5 - u32::from(self.precision))
    }

    /// UTM easting and northing (with false northing) of the south-west corner, metres.
    pub fn utm_sw_corner(&self) -> (u32, u32) {
        (self.easting, self.northing)
    }

    /// The reference as written.
    pub fn as_str(&self) -> &str {
        &self.text
    }

    /// Latitude and longitude of the square's south-west corner, degrees.
    pub fn sw_corner(&self) -> (f64, f64) {
        utm_to_geodetic(
            self.zone,
            self.hemisphere,
            f64::from(self.easting),
            f64::from(self.northing),
        )
    }

    /// Latitude and longitude of the square's centre, degrees.
    pub fn center(&self) -> (f64, f64) {
        let half = f64::from(self.square_size_m()) / 2.0;
        utm_to_geodetic(
            self.zone,
            self.hemisphere,
            f64::from(self.easting) + half,
            f64::from(self.northing) + half,
        )
    }
}

/// The original character at byte `at` of the compacted reference, for messages.
fn display_char(compact: &str, at: usize) -> char {
    compact[at..].chars().next().unwrap_or('?')
}

/// Approximate latitude extent of a grid square: its corners, plus its points on the
/// central meridian where lines of constant northing reach their extreme latitude.
fn square_latitude_range(
    zone: u8,
    hemisphere: Hemisphere,
    easting: u32,
    northing: u32,
    size: u32,
) -> (f64, f64) {
    let (e0, e1) = (f64::from(easting), f64::from(easting + size));
    let (n0, n1) = (f64::from(northing), f64::from(northing + size));
    let mut eastings = vec![e0, e1];
    if (e0..=e1).contains(&crate::utm::FALSE_EASTING) {
        eastings.push(crate::utm::FALSE_EASTING);
    }
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &e in &eastings {
        for n in [n0, n1] {
            let (lat, _) = utm_to_geodetic(zone, hemisphere, e, n);
            min = min.min(lat);
            max = max.max(lat);
        }
    }
    (min, max)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn kind(s: &str) -> MgrsErrorKind {
        Mgrs::parse(s).unwrap_err().kind
    }

    #[test]
    fn decodes_letters_and_digits() {
        let m = Mgrs::parse("18SUJ2337106519").unwrap();
        assert_eq!((m.zone(), m.band(), m.precision()), (18, 'S', 5));
        assert_eq!(m.utm_sw_corner(), (323_371, 4_306_519));
        assert_eq!(m.hemisphere(), Hemisphere::North);
        let m = Mgrs::parse("18s uj 2337 0651").unwrap();
        assert_eq!(m.utm_sw_corner(), (323_370, 4_306_510));
        assert_eq!(m.square_size_m(), 10);
        let m = Mgrs::parse("30UYC").unwrap();
        assert_eq!(
            (m.utm_sw_corner(), m.square_size_m()),
            ((700_000, 5_700_000), 100_000)
        );
        assert_eq!(Mgrs::parse("4QFJ1234567890").unwrap().zone(), 4);
    }

    #[test]
    fn rejects_malformed_references() {
        assert_eq!(kind(""), MgrsErrorKind::MissingZone);
        assert_eq!(kind("SUJ23"), MgrsErrorKind::MissingZone);
        assert_eq!(kind("61SUJ"), MgrsErrorKind::ZoneOutOfRange(61));
        assert_eq!(kind("0SUJ"), MgrsErrorKind::ZoneOutOfRange(0));
        assert_eq!(kind("18ZUJ"), MgrsErrorKind::Polar('Z'));
        assert_eq!(kind("ZAH1234512345"), MgrsErrorKind::Polar('Z'));
        assert_eq!(kind("bAN"), MgrsErrorKind::Polar('B'));
        assert_eq!(kind("18IUJ"), MgrsErrorKind::BadBand('I'));
        assert_eq!(kind("32XNA"), MgrsErrorKind::NonexistentZone(32));
        assert_eq!(kind("18S"), MgrsErrorKind::MissingSquare);
        assert_eq!(kind("18SU"), MgrsErrorKind::MissingSquare);
        // Zone 18 uses the S–Z column set.
        assert_eq!(kind("18SAJ"), MgrsErrorKind::BadColumn('A'));
        assert_eq!(kind("18SUW"), MgrsErrorKind::BadRow('W'));
        assert_eq!(kind("18SUJ233"), MgrsErrorKind::DigitCount(3));
        assert_eq!(kind("18SUJ233710651900"), MgrsErrorKind::DigitCount(12));
        assert_eq!(kind("18SUJ23A1"), MgrsErrorKind::BadDigits);
        // Row J in zone 18 starts 300 km (mod 2000 km) north of the false origin;
        // band C spans roughly 1120–2010 km of southern northing.
        assert_eq!(kind("18CUJ"), MgrsErrorKind::NotInBand);
        assert!(Mgrs::parse("18CUA").is_ok());
        // Found by review: band validity must not depend on the digits. Block 31NEJ
        // (800–900 km north) meets band N even where a 1 m square inside it does not.
        assert_eq!(
            Mgrs::parse("31NEJ0000090000").unwrap().utm_sw_corner(),
            (500_000, 890_000)
        );
        assert!(Mgrs::parse("31NEJ").is_ok());
        // Band M's last row ends at the equator; a southern block at 10 000 km would
        // lie north of it.
        assert_eq!(kind("31MEA0000000000"), MgrsErrorKind::NotInBand);
        assert_eq!(kind("31MEA"), MgrsErrorKind::NotInBand);
        assert_eq!(
            Mgrs::parse("31MEV").unwrap().utm_sw_corner(),
            (500_000, 9_900_000)
        );
    }

    #[test]
    fn messages_name_the_reference() {
        let message = Mgrs::parse("18ZUJ").unwrap_err().to_string();
        assert!(message.contains("\"18ZUJ\""));
        assert!(message.contains("UPS"));
    }
}
