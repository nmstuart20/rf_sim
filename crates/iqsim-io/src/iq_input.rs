//! Streaming validation and bounded random access for SigMF IQ-file inputs.
//!
//! Opening an input hashes its exact metadata bytes and makes one sequential
//! pass over the data. Readers verify each page whenever it is loaded; cached
//! pages retain the validated bytes until eviction. Sample normalization and
//! statistics use a fixed operation order independent of cache state.

use std::fs::{self, File};
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use serde_json::Value;
use sha2::{Digest, Sha256, Sha512};

use crate::IoError;

/// Sample datatype of an IQ input.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IqDatatype {
    /// Little-endian signed 16-bit I followed by signed 16-bit Q.
    Ci16Le,
    /// Little-endian IEEE 754 32-bit I followed by 32-bit Q.
    Cf32Le,
}

impl IqDatatype {
    /// SigMF `core:datatype` spelling.
    pub const fn name(self) -> &'static str {
        match self {
            Self::Ci16Le => "ci16_le",
            Self::Cf32Le => "cf32_le",
        }
    }

    /// Bytes per complex sample, including both I and Q.
    pub const fn bytes_per_sample(self) -> u64 {
        match self {
            Self::Ci16Le => 4,
            Self::Cf32Le => 8,
        }
    }
}

/// Bytes per verification page (1 MiB), a multiple of both sample sizes.
pub const PAGE_BYTES: u64 = 1 << 20;

/// Validated `.sigmf-meta` content.
#[derive(Clone, Debug, PartialEq)]
pub struct IqInputMeta {
    /// Sample representation declared by `core:datatype`.
    pub datatype: IqDatatype,
    /// Finite, positive `core:sample_rate` in complex samples per second.
    pub sample_rate_hz: f64,
    /// SigMF `core:version` string.
    pub version: String,
    /// Finite `core:frequency` of the first capture, when present.
    pub frequency_hz: Option<f64>,
    /// SHA-256 of the exact metadata bytes, including whitespace.
    pub meta_sha256: [u8; 32],
    /// Length of the metadata file in bytes.
    pub meta_byte_len: u64,
}

/// Result of the single streaming validation pass over `.sigmf-data`.
#[derive(Clone, Debug, PartialEq)]
pub struct IqInputSummary {
    /// Number of complex samples, always positive.
    pub sample_count: u64,
    /// Length of the data file in bytes.
    pub byte_len: u64,
    /// SHA-256 of all data bytes in file order.
    pub sha256: [u8; 32],
    /// SHA-256 of each verification page; the last page may be short.
    pub page_sha256: Vec<[u8; 32]>,
    /// Maximum normalized complex magnitude, computed in f64.
    pub peak_magnitude: f64,
    /// Sum of normalized `re * re + im * im` in f64, in file order.
    pub sum_power: f64,
}

/// An opened, validated IQ input.
#[derive(Debug)]
pub struct IqInput {
    meta: IqInputMeta,
    summary: Arc<IqInputSummary>,
    data_path: PathBuf,
}

impl IqInput {
    /// Validate a `.sigmf-meta` file and its `.sigmf-data` file.
    ///
    /// Reads data sequentially with at most one page of sample bytes in memory.
    /// Metadata and per-page hashes are retained separately from sample bytes.
    ///
    /// # Errors
    ///
    /// Rejects malformed metadata, empty or truncated data, non-finite cf32
    /// samples, invalid capture offsets, and mismatched `core:sha512`. File
    /// access errors retain their original I/O kind and associated path.
    pub fn open(meta_path: &Path, data_path: &Path) -> Result<Self, IoError> {
        let parsed = parse_meta(meta_path)?;
        let (summary, sha512) =
            scan_data(data_path, parsed.meta.datatype, parsed.sha512.is_some())?;
        if parsed.sha512 != sha512 {
            return Err(
                IoError::invalid("core:sha512 does not match the data file").at_path(meta_path)
            );
        }
        for (capture, start) in parsed.capture_starts {
            if start >= summary.sample_count {
                return Err(IoError::invalid(format!(
                    "captures[{capture}].core:sample_start must be less than sample_count ({})",
                    summary.sample_count
                ))
                .at_path(meta_path));
            }
        }
        Ok(Self {
            meta: parsed.meta,
            summary: Arc::new(summary),
            data_path: data_path.to_owned(),
        })
    }

    /// Validated metadata and its exact-byte identity.
    pub fn meta(&self) -> &IqInputMeta {
        &self.meta
    }

    /// Statistics and hashes from the streaming validation pass.
    pub fn summary(&self) -> &IqInputSummary {
        &self.summary
    }

    /// Path of the validated data file.
    pub fn data_path(&self) -> &Path {
        &self.data_path
    }

    /// Create an independent reader caching at most `max_pages` pages.
    ///
    /// The data file is opened lazily on the first nonempty read. Readers own
    /// their metadata and may outlive this input.
    ///
    /// # Panics
    ///
    /// Panics if `max_pages` is zero.
    pub fn reader(&self, max_pages: usize) -> IqFileReader {
        assert!(max_pages >= 1, "max_pages must be at least 1");
        IqFileReader {
            datatype: self.meta.datatype,
            summary: Arc::clone(&self.summary),
            data_path: self.data_path.clone(),
            file: None,
            max_pages,
            cache: Vec::new(),
        }
    }
}

struct ParsedMeta {
    meta: IqInputMeta,
    sha512: Option<[u8; 64]>,
    capture_starts: Vec<(usize, u64)>,
}

fn parse_meta(path: &Path) -> Result<ParsedMeta, IoError> {
    let bytes = fs::read(path).map_err(|error| IoError::io(error).at_path(path))?;
    parse_meta_bytes(&bytes).map_err(|error| error.at_path(path))
}

fn parse_meta_bytes(bytes: &[u8]) -> Result<ParsedMeta, IoError> {
    let value: Value = serde_json::from_slice(bytes)
        .map_err(|error| IoError::invalid(format!("invalid metadata JSON: {error}")))?;
    let root = value
        .as_object()
        .ok_or_else(|| IoError::invalid("metadata top-level must be a JSON object"))?;
    let global = root
        .get("global")
        .and_then(Value::as_object)
        .ok_or_else(|| IoError::invalid("global must be an object"))?;
    let captures = root
        .get("captures")
        .and_then(Value::as_array)
        .ok_or_else(|| IoError::invalid("captures must be an array"))?;
    if root.get("annotations").is_some_and(|v| !v.is_array()) {
        return Err(IoError::invalid("annotations must be an array"));
    }
    let datatype = match global.get("core:datatype").and_then(Value::as_str) {
        Some("ci16_le") => IqDatatype::Ci16Le,
        Some("cf32_le") => IqDatatype::Cf32Le,
        Some(other) => {
            return Err(IoError::invalid(format!(
                "unsupported core:datatype {other:?}"
            )));
        }
        None => return Err(IoError::invalid("core:datatype must be a string")),
    };
    let sample_rate_hz = finite_number(global.get("core:sample_rate"), "core:sample_rate")?;
    if sample_rate_hz <= 0.0 {
        return Err(IoError::invalid("core:sample_rate must be positive"));
    }
    let version = global
        .get("core:version")
        .and_then(Value::as_str)
        .ok_or_else(|| IoError::invalid("core:version must be a string"))?
        .to_owned();
    if let Some(channels) = global.get("core:num_channels")
        && channels.as_u64() != Some(1)
    {
        return Err(IoError::invalid("core:num_channels must be the integer 1"));
    }
    let sha512 = global
        .get("core:sha512")
        .map(|value| {
            let text = value
                .as_str()
                .filter(|text| text.len() == 128)
                .ok_or_else(|| IoError::invalid("core:sha512 must be a 128-hex-char string"))?;
            let mut hash = [0; 64];
            hex::decode_to_slice(text, &mut hash)
                .map_err(|_| IoError::invalid("core:sha512 must be a 128-hex-char string"))?;
            Ok(hash)
        })
        .transpose()?;
    let mut capture_starts = Vec::new();
    let mut previous_start = None;
    let mut frequency_hz = None;
    for (index, capture) in captures.iter().enumerate() {
        let capture = capture
            .as_object()
            .ok_or_else(|| IoError::invalid(format!("captures[{index}] must be an object")))?;
        if let Some(value) = capture.get("core:sample_start") {
            let start = value.as_u64().ok_or_else(|| {
                IoError::invalid(format!(
                    "captures[{index}].core:sample_start must be a non-negative integer"
                ))
            })?;
            if previous_start.is_some_and(|previous| start < previous) {
                return Err(IoError::invalid(format!(
                    "captures[{index}].core:sample_start must be non-decreasing"
                )));
            }
            capture_starts.push((index, start));
            previous_start = Some(start);
        }
        if let Some(value) = capture.get("core:frequency") {
            let frequency =
                finite_number(Some(value), &format!("captures[{index}].core:frequency"))?;
            if index == 0 {
                frequency_hz = Some(frequency);
            }
        }
    }
    Ok(ParsedMeta {
        meta: IqInputMeta {
            datatype,
            sample_rate_hz,
            version,
            frequency_hz,
            meta_sha256: Sha256::digest(bytes).into(),
            meta_byte_len: bytes.len() as u64,
        },
        sha512,
        capture_starts,
    })
}

fn finite_number(value: Option<&Value>, field: &str) -> Result<f64, IoError> {
    value
        .and_then(Value::as_f64)
        .filter(|value| value.is_finite())
        .ok_or_else(|| IoError::invalid(format!("{field} must be a finite number")))
}

fn scan_data(
    path: &Path,
    datatype: IqDatatype,
    need_sha512: bool,
) -> Result<(IqInputSummary, Option<[u8; 64]>), IoError> {
    let mut file = File::open(path).map_err(|error| IoError::io(error).at_path(path))?;
    let mut page = Vec::with_capacity(PAGE_BYTES as usize);
    let mut sha256 = Sha256::new();
    let mut sha512 = need_sha512.then(Sha512::new);
    let mut page_sha256 = Vec::new();
    let mut byte_len = 0_u64;
    let mut peak_power = 0.0_f64;
    let mut sum_power = 0.0;
    let sample_bytes = datatype.bytes_per_sample();
    loop {
        page.clear();
        file.by_ref()
            .take(PAGE_BYTES)
            .read_to_end(&mut page)
            .map_err(|error| IoError::io(error).at_path(path))?;
        if page.is_empty() {
            break;
        }
        if !(page.len() as u64).is_multiple_of(sample_bytes) {
            return Err(
                IoError::invalid("truncated data file: incomplete complex sample").at_path(path),
            );
        }
        let first = byte_len / sample_bytes;
        byte_len = byte_len
            .checked_add(page.len() as u64)
            .ok_or_else(|| IoError::invalid("data byte length overflow").at_path(path))?;
        for (index, bytes) in page.chunks_exact(sample_bytes as usize).enumerate() {
            let [re, im] = decode(datatype, bytes);
            if !re.is_finite() || !im.is_finite() {
                return Err(IoError::invalid(format!(
                    "non-finite cf32 value at sample index {}",
                    first + index as u64
                ))
                .at_path(path));
            }
            let re_power = re * re;
            let im_power = im * im;
            let power = re_power + im_power;
            peak_power = peak_power.max(power);
            sum_power += power;
        }
        sha256.update(&page);
        if let Some(hash) = &mut sha512 {
            hash.update(&page);
        }
        page_sha256.push(Sha256::digest(&page).into());
    }
    if byte_len == 0 {
        return Err(IoError::invalid("empty data file").at_path(path));
    }
    Ok((
        IqInputSummary {
            sample_count: byte_len / sample_bytes,
            byte_len,
            sha256: sha256.finalize().into(),
            page_sha256,
            peak_magnitude: peak_power.sqrt(),
            sum_power,
        },
        sha512.map(|hash| hash.finalize().into()),
    ))
}

fn decode(datatype: IqDatatype, bytes: &[u8]) -> [f64; 2] {
    match datatype {
        IqDatatype::Ci16Le => [
            f64::from(i16::from_le_bytes([bytes[0], bytes[1]])) / 32768.0,
            f64::from(i16::from_le_bytes([bytes[2], bytes[3]])) / 32768.0,
        ],
        IqDatatype::Cf32Le => [
            f64::from(f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
            f64::from(f32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]])),
        ],
    }
}

#[derive(Debug)]
struct CachedPage {
    index: u64,
    bytes: Vec<u8>,
}

/// Random-access reader over a data file through a bounded LRU page cache.
///
/// Every cache miss verifies the page against its hash at open. A cache hit
/// returns the already verified bytes even if the file has since changed.
#[derive(Debug)]
pub struct IqFileReader {
    datatype: IqDatatype,
    summary: Arc<IqInputSummary>,
    data_path: PathBuf,
    file: Option<File>,
    max_pages: usize,
    // Oldest use first, newest use last: no timestamps, ties, or counter wrap.
    cache: Vec<CachedPage>,
}

impl IqFileReader {
    /// Fill `out` with normalized complex samples `[first, first + out.len())`.
    ///
    /// ci16 components are divided by 32768 in f64; cf32 components are widened
    /// to f64. Both conversions are exact. An empty read at the end is allowed.
    ///
    /// # Errors
    ///
    /// Out-of-range indices or index overflow are `InvalidInput`. Changed pages
    /// and short reads are I/O errors with kind `InvalidData`. On an I/O error,
    /// earlier pages may already have filled part of `out`.
    pub fn read(&mut self, first: u64, out: &mut [[f64; 2]]) -> Result<(), IoError> {
        let count = u64::try_from(out.len()).map_err(|_| {
            IoError::invalid("sample range length overflow").at_path(&self.data_path)
        })?;
        let end = first
            .checked_add(count)
            .ok_or_else(|| IoError::invalid("sample range overflow").at_path(&self.data_path))?;
        if first > self.summary.sample_count || end > self.summary.sample_count {
            return Err(IoError::invalid("sample range is out of range").at_path(&self.data_path));
        }
        let datatype = self.datatype;
        let sample_bytes = datatype.bytes_per_sample() as usize;
        let samples_per_page = PAGE_BYTES / datatype.bytes_per_sample();
        let mut position = first;
        let mut remaining = out;
        while !remaining.is_empty() {
            let page_index = position / samples_per_page;
            let offset = (position % samples_per_page) as usize;
            let page = self.page(page_index)?;
            let count = remaining.len().min(page.len() / sample_bytes - offset);
            let (target, rest) = remaining.split_at_mut(count);
            for (sample, bytes) in target
                .iter_mut()
                .zip(page[offset * sample_bytes..].chunks_exact(sample_bytes))
            {
                *sample = decode(datatype, bytes);
            }
            position += count as u64;
            remaining = rest;
        }
        Ok(())
    }

    /// Number of pages currently cached, never greater than `max_pages`.
    pub fn cached_pages(&self) -> usize {
        self.cache.len()
    }

    fn page(&mut self, index: u64) -> Result<&[u8], IoError> {
        if let Some(position) = self.cache.iter().position(|page| page.index == index) {
            let page = self.cache.remove(position);
            self.cache.push(page);
        } else {
            if self.cache.len() == self.max_pages {
                self.cache.remove(0);
            }
            if self.file.is_none() {
                self.file = Some(
                    File::open(&self.data_path)
                        .map_err(|error| IoError::io(error).at_path(&self.data_path))?,
                );
            }
            let file = self.file.as_mut().expect("file opened above");
            let length = file
                .metadata()
                .map_err(|error| IoError::io(error).at_path(&self.data_path))?
                .len();
            if length != self.summary.byte_len {
                return Err(changed(&self.data_path));
            }
            let start = index * PAGE_BYTES;
            let len = PAGE_BYTES.min(self.summary.byte_len - start) as usize;
            let mut bytes = vec![0; len];
            file.seek(SeekFrom::Start(start))
                .map_err(|error| IoError::io(error).at_path(&self.data_path))?;
            file.read_exact(&mut bytes).map_err(|error| {
                if error.kind() == io::ErrorKind::UnexpectedEof {
                    changed(&self.data_path)
                } else {
                    IoError::io(error).at_path(&self.data_path)
                }
            })?;
            if <[u8; 32]>::from(Sha256::digest(&bytes)) != self.summary.page_sha256[index as usize]
            {
                return Err(changed(&self.data_path));
            }
            self.cache.push(CachedPage { index, bytes });
        }
        Ok(&self.cache.last().expect("page inserted above").bytes)
    }
}

fn changed(path: &Path) -> IoError {
    IoError::io(io::Error::new(
        io::ErrorKind::InvalidData,
        "IQ data file changed since it was opened",
    ))
    .at_path(path)
}
