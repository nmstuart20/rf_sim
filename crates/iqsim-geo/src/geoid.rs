//! Geoid height from a gridded model (contracts §12.5).
//!
//! **File format: GeographicLib PGM** (the `egm96-5.pgm`, `egm2008-1.pgm`, … files
//! that GeographicLib distributes). A binary PGM:
//!
//! ```text
//! P5
//! # Description WGS84 EGM96, 5-minute grid
//! # Offset -108
//! # Scale 0.003
//! # …other comment lines…
//! 4320 2161
//! 65535
//! <width·height big-endian u16 samples>
//! ```
//!
//! Sample `(row i, column j)` is the geoid height `Offset + Scale·v` metres at
//! latitude `90° − i·Δ` and longitude `j·Δ` east, with `Δ = 360°/width`, so rows run
//! from the north pole to the south pole inclusive (`height = width/2 + 1`) and
//! columns from 0° to `360° − Δ`. `Offset` and `Scale` are required; other comment
//! lines are kept only as text for the description. Anything else (another magic,
//! maxval ≠ 65535, an inconsistent shape, missing or extra data bytes) is an error.
//!
//! **Interpolation: bilinear**, in the same order as GeographicLib's bilinear mode:
//! longitude is reduced to [0°, 360°), `fx = λ·width/360`, `fy = (90° − φ)·(height−1)/180`,
//! the cell is `(⌊fy⌋, ⌊fx⌋)` with the column after the last wrapping to 0 and the
//! last row folded into the cell above it, and
//! `a = (1−fx)·v₀₀ + fx·v₀₁`, `b = (1−fx)·v₁₀ + fx·v₁₁`, `c = (1−fy)·a + fy·b`,
//! `N = Offset + Scale·c` (with fx, fy the fractional parts). The result is a
//! continuous surface; its error against the full model is the file's own
//! documented bilinear error (EGM96 5′: 0.140 m max).
//!
//! When no grid is configured, [`Geoid::Zero`] gives N = 0 everywhere: MSL equals the
//! ellipsoid. The runtime records that choice and a warning in the run manifest.

use std::path::Path;

use crate::GeoError;

/// The format name recorded in the manifest.
pub const GEOID_FORMAT: &str = "geographiclib-pgm";

/// The interpolation name recorded in the manifest.
pub const GEOID_INTERPOLATION: &str = "bilinear";

/// A gridded geoid model (module docs).
#[derive(Clone, Debug, PartialEq)]
pub struct GeoidGrid {
    width: usize,
    height: usize,
    offset: f64,
    scale: f64,
    description: String,
    data: Vec<u16>,
}

/// Header values that change how a grid's bytes are interpreted.
#[derive(Clone, Debug, PartialEq)]
pub struct GeoidGridMetadata {
    /// [`GEOID_FORMAT`].
    pub format: &'static str,
    /// Columns (longitudes).
    pub width: usize,
    /// Rows (latitudes), pole to pole inclusive.
    pub height: usize,
    /// `Offset`, metres.
    pub offset_m: f64,
    /// `Scale`, metres per count.
    pub scale_m: f64,
    /// [`GEOID_INTERPOLATION`].
    pub interpolation: &'static str,
    /// The `# Description` text, or empty.
    pub description: String,
}

fn bad(msg: impl Into<String>) -> GeoError {
    GeoError::GeoidFormat(msg.into())
}

/// Skips a `#` comment at `pos` up to (not including) its CR or LF terminator, and
/// passes its text to `comment` (netpbm comment rule).
fn skip_comment(
    bytes: &[u8],
    pos: &mut usize,
    comment: &mut impl FnMut(&str),
) -> Result<(), GeoError> {
    let start = *pos + 1;
    while *pos < bytes.len() && bytes[*pos] != b'\n' && bytes[*pos] != b'\r' {
        *pos += 1;
    }
    let line = std::str::from_utf8(&bytes[start..*pos])
        .map_err(|_| bad("a header comment is not UTF-8"))?;
    comment(line.trim());
    Ok(())
}

/// Reads one whitespace-delimited header token, skipping `#` comment lines, which are
/// passed to `comment`.
fn token<'a>(
    bytes: &'a [u8],
    pos: &mut usize,
    comment: &mut impl FnMut(&str),
) -> Result<&'a str, GeoError> {
    loop {
        while *pos < bytes.len() && bytes[*pos].is_ascii_whitespace() {
            *pos += 1;
        }
        if *pos >= bytes.len() {
            return Err(bad("the header ends early"));
        }
        if bytes[*pos] == b'#' {
            skip_comment(bytes, pos, comment)?;
            continue;
        }
        let start = *pos;
        while *pos < bytes.len() && !bytes[*pos].is_ascii_whitespace() && bytes[*pos] != b'#' {
            *pos += 1;
        }
        return std::str::from_utf8(&bytes[start..*pos])
            .map_err(|_| bad("a header token is not UTF-8"));
    }
}

impl GeoidGrid {
    /// Parses a GeographicLib PGM grid (module docs).
    pub fn from_pgm_bytes(bytes: &[u8]) -> Result<Self, GeoError> {
        let mut pos = 0;
        let mut offset = None;
        let mut scale = None;
        let mut description = String::new();
        let mut comment_err = None;
        let mut on_comment = |line: &str| {
            let (key, value) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
            let value = value.trim();
            let parse = |v: &str| v.parse::<f64>().ok().filter(|x| x.is_finite());
            match key {
                "Offset" => match parse(value) {
                    Some(v) if offset.is_none() => offset = Some(v),
                    Some(_) => comment_err = Some(bad("Offset appears twice")),
                    None => comment_err = Some(bad(format!("bad Offset {value:?}"))),
                },
                "Scale" => match parse(value) {
                    Some(v) if scale.is_none() && v > 0.0 => scale = Some(v),
                    Some(_) => comment_err = Some(bad("Scale appears twice or is not positive")),
                    None => comment_err = Some(bad(format!("bad Scale {value:?}"))),
                },
                "Description" => description = value.to_string(),
                _ => {}
            }
        };
        if token(bytes, &mut pos, &mut on_comment)? != "P5" {
            return Err(bad("not a binary PGM (magic P5)"));
        }
        let mut num = |what: &str| -> Result<usize, GeoError> {
            let t = token(bytes, &mut pos, &mut on_comment)?;
            t.parse::<usize>()
                .map_err(|_| bad(format!("bad {what} {t:?}")))
        };
        let width = num("width")?;
        let height = num("height")?;
        let maxval = num("maxval")?;
        // A comment may follow maxval directly; its CR or LF terminator is then the
        // single whitespace byte before the raster.
        if pos < bytes.len() && bytes[pos] == b'#' {
            skip_comment(bytes, &mut pos, &mut on_comment)?;
        }
        drop(on_comment);
        if let Some(e) = comment_err {
            return Err(e);
        }
        let (Some(offset), Some(scale)) = (offset, scale) else {
            return Err(bad("the header has no Offset or no Scale comment"));
        };
        if maxval != 65_535 {
            return Err(bad(format!("maxval is {maxval}, not 65535")));
        }
        if width < 4 || width % 2 != 0 || height != width / 2 + 1 {
            return Err(bad(format!(
                "a {width}×{height} grid is not a global pole-to-pole grid (height = width/2 + 1)"
            )));
        }
        // Exactly one whitespace byte separates maxval from the data.
        if pos >= bytes.len() || !bytes[pos].is_ascii_whitespace() {
            return Err(bad("no data after the header"));
        }
        pos += 1;
        let n = width
            .checked_mul(height)
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| bad("the grid is too large"))?;
        let data = &bytes[pos..];
        if data.len() != n {
            return Err(bad(format!(
                "expected {n} data bytes for {width}×{height}, found {}",
                data.len()
            )));
        }
        let data = data
            .as_chunks::<2>()
            .0
            .iter()
            .map(|c| u16::from_be_bytes(*c))
            .collect();
        Ok(Self {
            width,
            height,
            offset,
            scale,
            description,
            data,
        })
    }

    /// Reads and parses a grid file.
    pub fn from_pgm_file(path: &Path) -> Result<Self, GeoError> {
        let bytes =
            std::fs::read(path).map_err(|e| GeoError::Io(format!("{}: {e}", path.display())))?;
        Self::from_pgm_bytes(&bytes)
    }

    /// The header values that define how the samples are read.
    pub fn metadata(&self) -> GeoidGridMetadata {
        GeoidGridMetadata {
            format: GEOID_FORMAT,
            width: self.width,
            height: self.height,
            offset_m: self.offset,
            scale_m: self.scale,
            interpolation: GEOID_INTERPOLATION,
            description: self.description.clone(),
        }
    }

    fn at(&self, row: usize, col: usize) -> f64 {
        f64::from(self.data[row * self.width + col])
    }

    /// Geoid height N at geodetic latitude/longitude (radians), metres (module docs).
    pub fn height_m(&self, lat_rad: f64, lon_rad: f64) -> Result<f64, GeoError> {
        if !(lat_rad.is_finite() && lon_rad.is_finite()) {
            return Err(GeoError::NonFinite("geoid latitude or longitude"));
        }
        let rad_to_deg = 180.0 / core::f64::consts::PI;
        let lat = (lat_rad * rad_to_deg).clamp(-90.0, 90.0);
        let mut lon = lon_rad * rad_to_deg;
        lon = lon - 360.0 * (lon / 360.0).floor();
        if !(0.0..360.0).contains(&lon) {
            lon = 0.0;
        }
        let rlon = self.width as f64 / 360.0;
        let rlat = (self.height - 1) as f64 / 180.0;
        let mut fx = lon * rlon;
        let mut fy = (90.0 - lat) * rlat;
        let mut ix = fx.floor() as usize;
        let mut iy = fy.floor() as usize;
        fx -= ix as f64;
        fy -= iy as f64;
        if ix >= self.width {
            ix = 0;
            fx = 0.0;
        }
        if iy >= self.height - 1 {
            iy = self.height - 2;
            fy = 1.0;
        }
        let ix1 = if ix + 1 == self.width { 0 } else { ix + 1 };
        let (v00, v01) = (self.at(iy, ix), self.at(iy, ix1));
        let (v10, v11) = (self.at(iy + 1, ix), self.at(iy + 1, ix1));
        let a = (1.0 - fx) * v00 + fx * v01;
        let b = (1.0 - fx) * v10 + fx * v11;
        let c = (1.0 - fy) * a + fy * b;
        Ok(self.offset + self.scale * c)
    }
}

/// The geoid model a run uses.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Geoid {
    /// No grid configured: N = 0 everywhere (MSL is the ellipsoid).
    #[default]
    Zero,
    /// A gridded model.
    Grid(GeoidGrid),
}

impl Geoid {
    /// Geoid height N at geodetic latitude/longitude (radians), metres.
    pub fn height_m(&self, lat_rad: f64, lon_rad: f64) -> Result<f64, GeoError> {
        match self {
            Self::Zero => {
                if lat_rad.is_finite() && lon_rad.is_finite() {
                    Ok(0.0)
                } else {
                    Err(GeoError::NonFinite("geoid latitude or longitude"))
                }
            }
            Self::Grid(g) => g.height_m(lat_rad, lon_rad),
        }
    }
}
