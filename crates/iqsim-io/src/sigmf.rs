//! Sequential SigMF recordings with checked writes and separate finalization.
//!
//! Only consumer-visible core metadata is accepted. Truth and run manifests
//! belong in separate sidecars. Dropping a writer without finishing leaves no
//! metadata file and does not establish a successful recording.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use iqsim_core::SampleRate;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

use crate::IoError;

/// Little-endian format of one interleaved complex sample.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SampleFormat {
    /// Signed 16-bit I followed by signed 16-bit Q.
    Ci16Le,
    /// IEEE 754 32-bit I followed by IEEE 754 32-bit Q.
    Cf32Le,
}

impl SampleFormat {
    /// SigMF `core:datatype` spelling.
    pub const fn datatype(self) -> &'static str {
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

/// Consumer-visible recording metadata shared by every capture.
#[derive(Clone, Debug)]
pub struct SigmfGlobal {
    /// Exact rational rate in complex samples per second.
    pub sample_rate: SampleRate,
    /// Optional consumer-visible description.
    pub description: Option<String>,
    /// Recording software identification.
    pub recorder: String,
}

/// Consumer-visible RF context beginning at a complex sample index.
#[derive(Clone, Debug)]
pub struct SigmfCapture {
    /// First complex sample governed by this capture's context.
    pub sample_start: u64,
    /// RF reference frequency in hertz; must be finite.
    pub frequency_hz: f64,
    /// Optional UTC label from the receiver's reported timestamp clock.
    pub datetime: Option<String>,
}

/// Checked sample encoding over an arbitrary byte writer.
///
/// Invalid input writes nothing and can be corrected. An underlying write or
/// flush error permanently fails the sink, including when a partial write has
/// already occurred. The hash includes exactly the bytes accepted by `W`.
/// For a buffered `W`, acceptance does not imply durability; finalization must
/// flush and synchronize the underlying file.
#[derive(Debug)]
pub struct CheckedSampleSink<W: Write> {
    writer: W,
    format: SampleFormat,
    sample_count: u64,
    hash: Sha256,
    failed: Option<IoError>,
}

impl<W: Write> CheckedSampleSink<W> {
    /// Start an empty sink using the specified complex sample format.
    pub fn new(writer: W, format: SampleFormat) -> Self {
        Self {
            writer,
            format,
            sample_count: 0,
            hash: Sha256::new(),
            failed: None,
        }
    }

    /// Number of complex samples in fully successful write calls.
    ///
    /// After an I/O error the sink is unusable; its last partial call is not
    /// included in this count.
    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }

    /// SHA-256 of the exact bytes accepted so far, including any partial write.
    pub fn sha256(&self) -> [u8; 32] {
        self.hash.clone().finalize().into()
    }

    /// Encode and write interleaved signed 16-bit I,Q values in little endian.
    ///
    /// # Errors
    ///
    /// Rejects a wrong format, odd element count, count/length overflow, or a
    /// failed sink. Propagates every underlying write failure.
    pub fn write_ci16(&mut self, iq: &[i16]) -> Result<(), IoError> {
        let next_count = self.prepare(SampleFormat::Ci16Le, iq.len())?;
        self.write_encoded(iq, i16::to_le_bytes)?;
        self.sample_count = next_count;
        Ok(())
    }

    /// Encode and write interleaved finite 32-bit I,Q values in little endian.
    ///
    /// # Errors
    ///
    /// Rejects a wrong format, odd element count, any non-finite value,
    /// count/length overflow, or a failed sink. Validation precedes all writes.
    /// Propagates every underlying write failure.
    pub fn write_cf32(&mut self, iq: &[f32]) -> Result<(), IoError> {
        let next_count = self.prepare(SampleFormat::Cf32Le, iq.len())?;
        if iq.iter().any(|value| !value.is_finite()) {
            return Err(IoError::invalid("cf32 values must all be finite"));
        }
        self.write_encoded(iq, f32::to_le_bytes)?;
        self.sample_count = next_count;
        Ok(())
    }

    /// Flush the byte writer, permanently failing the sink if flushing fails.
    pub fn flush(&mut self) -> Result<(), IoError> {
        self.check_ready()?;
        match self.writer.flush() {
            Ok(()) => Ok(()),
            Err(error) => Err(self.fail(error)),
        }
    }

    fn check_ready(&self) -> Result<(), IoError> {
        match &self.failed {
            Some(error) => Err(error.clone()),
            None => Ok(()),
        }
    }

    fn fail(&mut self, error: io::Error) -> IoError {
        let error = IoError::io(error);
        self.failed = Some(error.clone());
        error
    }

    fn prepare(&self, format: SampleFormat, len: usize) -> Result<u64, IoError> {
        self.check_ready()?;
        if self.format != format {
            return Err(IoError::invalid("sample format does not match this sink"));
        }
        if !len.is_multiple_of(2) {
            return Err(IoError::invalid("interleaved I,Q length must be even"));
        }
        let additional = u64::try_from(len / 2)
            .map_err(|_| IoError::invalid("sample count does not fit in u64"))?;
        let next = self
            .sample_count
            .checked_add(additional)
            .ok_or_else(|| IoError::invalid("sample count overflow"))?;
        next.checked_mul(format.bytes_per_sample())
            .ok_or_else(|| IoError::invalid("recording byte length overflow"))?;
        Ok(next)
    }

    fn write_encoded<T: Copy, const N: usize>(
        &mut self,
        iq: &[T],
        encode: impl Fn(T) -> [u8; N],
    ) -> Result<(), IoError> {
        // Only called with N = 2 or 4. Bound temporary memory independently of
        // the caller's chunk length; never allocate a second whole IQ block.
        let mut bytes = [0_u8; 8_192];
        for chunk in iq.chunks(bytes.len() / N) {
            for (value, target) in chunk.iter().zip(bytes.as_chunks_mut::<N>().0.iter_mut()) {
                target.copy_from_slice(&encode(*value));
            }
            self.write_bytes(&bytes[..chunk.len() * N])?;
        }
        Ok(())
    }

    fn write_bytes(&mut self, mut bytes: &[u8]) -> Result<(), IoError> {
        while !bytes.is_empty() {
            match self.writer.write(bytes) {
                Ok(0) => {
                    return Err(self.fail(io::Error::new(
                        io::ErrorKind::WriteZero,
                        "sample writer accepted zero bytes",
                    )));
                }
                Ok(written) => {
                    self.hash.update(&bytes[..written]);
                    bytes = &bytes[written..];
                }
                Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
                Err(error) => return Err(self.fail(error)),
            }
        }
        Ok(())
    }
}

/// A recording whose metadata is created only by explicit successful finalization.
#[derive(Debug)]
pub struct SigmfWriter {
    sink: CheckedSampleSink<BufWriter<File>>,
    data_path: PathBuf,
    meta_path: PathBuf,
    global: SigmfGlobal,
    captures: Vec<SigmfCapture>,
}

impl SigmfWriter {
    /// Create `<stem>.sigmf-data` exclusively in an existing directory.
    ///
    /// The stem must be nonempty, must not start with a dot, and must not contain
    /// either path separator or NUL. At least one capture is required, starting
    /// at zero and ordered by nondecreasing `sample_start`, with finite RF
    /// frequencies. An existing data or metadata path is never overwritten.
    ///
    /// # Errors
    ///
    /// Returns validation errors or a path-qualified filesystem error.
    pub fn create(
        dir: &Path,
        stem: &str,
        format: SampleFormat,
        global: SigmfGlobal,
        captures: Vec<SigmfCapture>,
    ) -> Result<Self, IoError> {
        let (data_path, meta_path) = recording_paths(dir, stem)?;
        validate_captures(&captures).map_err(|error| error.at_path(&meta_path))?;
        // symlink_metadata also rejects dangling metadata symlinks as existing.
        match fs::symlink_metadata(&meta_path) {
            Ok(_) => {
                return Err(IoError::io(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "recording metadata already exists",
                ))
                .at_path(&meta_path));
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(IoError::io(error).at_path(&meta_path)),
        }
        let file =
            File::create_new(&data_path).map_err(|error| IoError::io(error).at_path(&data_path))?;
        Ok(Self {
            sink: CheckedSampleSink::new(BufWriter::new(file), format),
            data_path,
            meta_path,
            global,
            captures,
        })
    }

    /// Write interleaved signed 16-bit I,Q; see [`CheckedSampleSink::write_ci16`].
    pub fn write_ci16(&mut self, iq: &[i16]) -> Result<(), IoError> {
        self.sink
            .write_ci16(iq)
            .map_err(|error| error.at_path(&self.data_path))
    }

    /// Write finite interleaved 32-bit I,Q; see [`CheckedSampleSink::write_cf32`].
    pub fn write_cf32(&mut self, iq: &[f32]) -> Result<(), IoError> {
        self.sink
            .write_cf32(iq)
            .map_err(|error| error.at_path(&self.data_path))
    }

    /// Flush and sync data, verify its length, then create and sync metadata.
    ///
    /// The pretty JSON has sorted object keys, a trailing newline, and only
    /// the allowed SigMF core fields. The returned count and hash belong in
    /// the separate run manifest, not in consumer metadata.
    ///
    /// # Errors
    ///
    /// Fails for an earlier write error, flush/sync/metadata I/O failure,
    /// inconsistent data length, or a metadata path created since `create`.
    /// A metadata write/sync failure triggers removal of the incomplete file;
    /// any cleanup failure is included in the returned error.
    pub fn finish(mut self) -> Result<SigmfSummary, IoError> {
        self.sink
            .flush()
            .map_err(|error| error.at_path(&self.data_path))?;
        self.sink
            .writer
            .get_ref()
            .sync_all()
            .map_err(|error| IoError::io(error).at_path(&self.data_path))?;
        let byte_len = self
            .sink
            .sample_count
            .checked_mul(self.sink.format.bytes_per_sample())
            .ok_or_else(|| IoError::invalid("recording byte length overflow"))
            .map_err(|error| error.at_path(&self.data_path))?;
        let actual_len = fs::metadata(&self.data_path)
            .map_err(|error| IoError::io(error).at_path(&self.data_path))?
            .len();
        if actual_len != byte_len {
            return Err(IoError::Io {
                path: Some(self.data_path.clone()),
                kind: io::ErrorKind::InvalidData,
                message: format!("data length is {actual_len} bytes, expected {byte_len}"),
            });
        }
        let meta = metadata_bytes(self.sink.format, &self.global, &self.captures)
            .map_err(|error| error.at_path(&self.meta_path))?;
        let mut file = File::create_new(&self.meta_path)
            .map_err(|error| IoError::io(error).at_path(&self.meta_path))?;
        if let Err(error) = file.write_all(&meta).and_then(|()| file.sync_all()) {
            drop(file);
            let mut error = error;
            if let Err(cleanup) = fs::remove_file(&self.meta_path) {
                error = io::Error::new(
                    error.kind(),
                    format!("{error}; removing incomplete metadata also failed: {cleanup}"),
                );
            }
            return Err(IoError::io(error).at_path(&self.meta_path));
        }
        Ok(SigmfSummary {
            data_path: self.data_path,
            meta_path: self.meta_path,
            format: self.sink.format,
            sample_count: self.sink.sample_count,
            byte_len,
            sha256: self.sink.sha256(),
        })
    }
}

/// Finalized recording identity for a separate run manifest or verification.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SigmfSummary {
    /// Path to the synchronized little-endian sample bytes.
    pub data_path: PathBuf,
    /// Path to the synchronized consumer metadata JSON.
    pub meta_path: PathBuf,
    /// Complex sample encoding.
    pub format: SampleFormat,
    /// Total number of complex samples.
    pub sample_count: u64,
    /// Verified data file length in bytes.
    pub byte_len: u64,
    /// SHA-256 of the sample bytes only.
    pub sha256: [u8; 32],
}

impl SigmfSummary {
    /// Lowercase hexadecimal SHA-256 of the data file's bytes.
    pub fn sha256_hex(&self) -> String {
        hex::encode(self.sha256)
    }
}

/// Parsed consumer metadata and raw sample bytes read from a recording pair.
#[derive(Clone, Debug)]
pub struct Recording {
    /// Parsed `.sigmf-meta` JSON, without added fields or transformations.
    pub meta: Value,
    /// Exact `.sigmf-data` bytes.
    pub data: Vec<u8>,
}

/// Read a recording pair for inspection or subsequent replay validation.
///
/// This helper validates the stem and parses JSON. Manifest identities, hashes
/// and metadata semantics remain the replay caller's responsibility.
pub fn read_recording(dir: &Path, stem: &str) -> Result<Recording, IoError> {
    let (data_path, meta_path) = recording_paths(dir, stem)?;
    let bytes = fs::read(&meta_path).map_err(|error| IoError::io(error).at_path(&meta_path))?;
    let meta = serde_json::from_slice(&bytes).map_err(|error| IoError::Io {
        path: Some(meta_path),
        kind: io::ErrorKind::InvalidData,
        message: format!("invalid SigMF JSON: {error}"),
    })?;
    let data = fs::read(&data_path).map_err(|error| IoError::io(error).at_path(&data_path))?;
    Ok(Recording { meta, data })
}

/// Decode little-endian i16 elements, retaining interleaved I,Q ordering.
///
/// Returns an error if the byte count is not a multiple of two.
pub fn decode_ci16(bytes: &[u8]) -> Result<Vec<i16>, IoError> {
    if !bytes.len().is_multiple_of(2) {
        return Err(IoError::invalid(
            "ci16 element byte length must be a multiple of 2",
        ));
    }
    Ok(bytes
        .as_chunks::<2>()
        .0
        .iter()
        .map(|part| i16::from_le_bytes(*part))
        .collect())
}

/// Decode little-endian f32 elements, retaining interleaved I,Q ordering.
///
/// Returns an error if the byte count is not a multiple of four. All IEEE 754
/// bit patterns are preserved; finite-value validation belongs to the caller.
pub fn decode_cf32(bytes: &[u8]) -> Result<Vec<f32>, IoError> {
    if !bytes.len().is_multiple_of(4) {
        return Err(IoError::invalid(
            "cf32 element byte length must be a multiple of 4",
        ));
    }
    Ok(bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|part| f32::from_le_bytes(*part))
        .collect())
}

fn recording_paths(dir: &Path, stem: &str) -> Result<(PathBuf, PathBuf), IoError> {
    if stem.is_empty() || stem.starts_with('.') || stem.contains(['/', '\\', '\0']) {
        return Err(IoError::invalid(
            "recording stem must be nonempty, must not start with '.', and must not contain separators or NUL",
        )
        .at_path(dir));
    }
    Ok((
        dir.join(format!("{stem}.sigmf-data")),
        dir.join(format!("{stem}.sigmf-meta")),
    ))
}

fn validate_captures(captures: &[SigmfCapture]) -> Result<(), IoError> {
    if captures
        .first()
        .is_none_or(|capture| capture.sample_start != 0)
    {
        return Err(IoError::invalid(
            "captures must be nonempty and start at sample 0",
        ));
    }
    if captures
        .windows(2)
        .any(|pair| pair[0].sample_start > pair[1].sample_start)
    {
        return Err(IoError::invalid("captures must be sorted by sample_start"));
    }
    if captures
        .iter()
        .any(|capture| !capture.frequency_hz.is_finite())
    {
        return Err(IoError::invalid("capture frequency must be finite"));
    }
    Ok(())
}

fn metadata_bytes(
    format: SampleFormat,
    global: &SigmfGlobal,
    captures: &[SigmfCapture],
) -> Result<Vec<u8>, IoError> {
    // BTreeMaps also preserve sorted insertion if serde_json's preserve_order
    // feature is ever enabled through workspace feature unification.
    let mut global_fields = BTreeMap::from([
        ("core:datatype", json!(format.datatype())),
        (
            "core:sample_rate",
            json!(global.sample_rate.num() as f64 / global.sample_rate.den() as f64),
        ),
        ("core:version", json!("1.2.0")),
        ("core:num_channels", json!(1)),
        ("core:recorder", json!(global.recorder)),
    ]);
    if let Some(description) = &global.description {
        global_fields.insert("core:description", json!(description));
    }
    let capture_fields: Vec<_> = captures
        .iter()
        .map(|capture| {
            let mut fields = BTreeMap::from([
                ("core:sample_start", json!(capture.sample_start)),
                ("core:frequency", json!(capture.frequency_hz)),
            ]);
            if let Some(datetime) = &capture.datetime {
                fields.insert("core:datetime", json!(datetime));
            }
            fields
        })
        .collect();
    let json_error = |error: serde_json::Error| IoError::Io {
        path: None,
        kind: io::ErrorKind::InvalidData,
        message: format!("serializing SigMF metadata: {error}"),
    };
    let document = BTreeMap::from([
        ("annotations", json!([])),
        (
            "captures",
            serde_json::to_value(capture_fields).map_err(json_error)?,
        ),
        (
            "global",
            serde_json::to_value(global_fields).map_err(json_error)?,
        ),
    ]);
    let mut bytes = serde_json::to_vec_pretty(&document).map_err(json_error)?;
    bytes.push(b'\n');
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug)]
    struct TestWriter {
        bytes: Vec<u8>,
        fail_after: usize,
        max_write: usize,
        interrupt_once: bool,
        fail_flush: bool,
        calls: usize,
    }

    impl TestWriter {
        fn new(fail_after: usize) -> Self {
            Self {
                bytes: Vec::new(),
                fail_after,
                max_write: usize::MAX,
                interrupt_once: false,
                fail_flush: false,
                calls: 0,
            }
        }
    }

    impl Write for TestWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            self.calls += 1;
            if self.interrupt_once {
                self.interrupt_once = false;
                return Err(io::Error::from(io::ErrorKind::Interrupted));
            }
            if self.bytes.len() >= self.fail_after {
                return Err(io::Error::new(
                    io::ErrorKind::StorageFull,
                    "injected disk failure",
                ));
            }
            let len = bytes
                .len()
                .min(self.fail_after - self.bytes.len())
                .min(self.max_write);
            self.bytes.extend_from_slice(&bytes[..len]);
            Ok(len)
        }

        fn flush(&mut self) -> io::Result<()> {
            if self.fail_flush {
                Err(io::Error::other("injected flush failure"))
            } else {
                Ok(())
            }
        }
    }

    #[test]
    fn partial_write_failures_propagate_and_permanently_fail_the_sink() {
        let iq = [-32_768_i16, -1, 0, 32_767];
        let expected: Vec<_> = iq.iter().flat_map(|value| value.to_le_bytes()).collect();
        for fail_after in 0..expected.len() {
            let mut sink =
                CheckedSampleSink::new(TestWriter::new(fail_after), SampleFormat::Ci16Le);
            let error = sink.write_ci16(&iq).unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::StorageFull);
            assert!(error.to_string().contains("injected disk failure"));
            assert_eq!(sink.writer.bytes, expected[..fail_after]);
            assert_eq!(
                sink.sha256(),
                <[u8; 32]>::from(Sha256::digest(&expected[..fail_after]))
            );
            assert_eq!(sink.sample_count(), 0);
            let calls = sink.writer.calls;
            assert_eq!(sink.write_ci16(&iq).unwrap_err(), error);
            assert_eq!(sink.write_ci16(&[]).unwrap_err(), error);
            assert_eq!(sink.write_cf32(&[]).unwrap_err(), error);
            assert_eq!(sink.flush().unwrap_err(), error);
            assert_eq!(sink.writer.calls, calls);
        }
    }

    #[test]
    fn short_writes_and_interrupts_preserve_every_byte() {
        let mut writer = TestWriter::new(usize::MAX);
        writer.max_write = 3;
        writer.interrupt_once = true;
        let mut sink = CheckedSampleSink::new(writer, SampleFormat::Cf32Le);
        let iq = [1.0_f32, -0.0, f32::MIN_POSITIVE, f32::from_bits(1)];
        sink.write_cf32(&iq).unwrap();
        sink.flush().unwrap();
        let expected: Vec<_> = iq.iter().flat_map(|value| value.to_le_bytes()).collect();
        assert_eq!(sink.writer.bytes, expected);
        assert_eq!(sink.sample_count(), 2);
        assert_eq!(sink.sha256(), <[u8; 32]>::from(Sha256::digest(&expected)));
    }

    #[test]
    fn zero_writes_and_flush_errors_are_sticky() {
        let mut writer = TestWriter::new(usize::MAX);
        writer.max_write = 0;
        let mut sink = CheckedSampleSink::new(writer, SampleFormat::Ci16Le);
        assert_eq!(
            sink.write_ci16(&[1, 2]).unwrap_err().kind(),
            io::ErrorKind::WriteZero
        );
        assert!(sink.write_ci16(&[1, 2]).is_err());

        let mut writer = TestWriter::new(usize::MAX);
        writer.fail_flush = true;
        let mut sink = CheckedSampleSink::new(writer, SampleFormat::Ci16Le);
        sink.write_ci16(&[1, 2]).unwrap();
        let error = sink.flush().unwrap_err();
        assert!(error.to_string().contains("injected flush failure"));
        assert_eq!(sink.write_ci16(&[3, 4]).unwrap_err(), error);
        assert_eq!(sink.writer.bytes, [1, 0, 2, 0]);
    }

    #[test]
    fn validation_and_count_overflow_write_nothing() {
        let mut sink = CheckedSampleSink::new(Vec::new(), SampleFormat::Cf32Le);
        // The invalid value is beyond one internal byte-encoding chunk.
        let mut iq = vec![0.0; 8_194];
        iq[8_193] = f32::NAN;
        assert!(sink.write_cf32(&iq).is_err());
        assert!(sink.write_cf32(&[1.0]).is_err());
        assert!(sink.write_ci16(&[1, 2]).is_err());
        assert!(sink.writer.is_empty());
        sink.write_cf32(&[1.0, 2.0]).unwrap();
        assert_eq!(sink.sample_count(), 1);

        let mut sink = CheckedSampleSink::new(Vec::new(), SampleFormat::Ci16Le);
        sink.sample_count = u64::MAX;
        assert!(sink.write_ci16(&[1, 2]).is_err());
        sink.sample_count = u64::MAX / 4;
        assert!(sink.write_ci16(&[1, 2]).is_err());
        assert!(sink.writer.is_empty());
    }
}
