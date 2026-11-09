//! SigMF recording round trips, metadata isolation, failures, and golden bytes.

use std::fs;
use std::io::ErrorKind;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use iqsim_core::SampleRate;
use iqsim_io::sigmf::{
    SampleFormat, SigmfCapture, SigmfGlobal, SigmfSummary, SigmfWriter, decode_cf32, decode_ci16,
    read_recording,
};
use serde_json::Value;
use sha2::{Digest, Sha256};

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

struct TempDir(PathBuf);

impl TempDir {
    fn new() -> Self {
        loop {
            let count = NEXT_DIR.fetch_add(1, Ordering::Relaxed);
            let path =
                std::env::temp_dir().join(format!("iqsim-io-{}-{count}", std::process::id()));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(error) if error.kind() == ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("creating {}: {error}", path.display()),
            }
        }
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).expect("clean up recording test directory");
    }
}

fn global(description: Option<&str>) -> SigmfGlobal {
    SigmfGlobal {
        sample_rate: SampleRate::from_ratio(10_000_000, 3).unwrap(),
        description: description.map(str::to_owned),
        recorder: "iqsim-io test".to_owned(),
    }
}

fn capture(sample_start: u64, datetime: Option<&str>) -> SigmfCapture {
    SigmfCapture {
        sample_start,
        frequency_hz: 1_090_000_000.25,
        datetime: datetime.map(str::to_owned),
    }
}

fn assert_keys(object: &Value, expected: &[&str]) {
    let keys: Vec<_> = object
        .as_object()
        .unwrap()
        .keys()
        .map(String::as_str)
        .collect();
    assert_eq!(keys, expected);
}

fn assert_summary(summary: &SigmfSummary, bytes: &[u8], count: u64, format: SampleFormat) {
    assert_eq!(summary.format, format);
    assert_eq!(summary.sample_count, count);
    assert_eq!(summary.byte_len, count * format.bytes_per_sample());
    assert_eq!(summary.byte_len, bytes.len() as u64);
    assert_eq!(summary.sha256, <[u8; 32]>::from(Sha256::digest(bytes)));
    assert_eq!(summary.sha256_hex(), hex::encode(Sha256::digest(bytes)));
    assert_eq!(
        fs::metadata(&summary.data_path).unwrap().len(),
        summary.byte_len
    );
    assert!(summary.meta_path.is_file());
}

#[test]
fn round_trip_ci16_with_only_allowed_metadata_and_sorted_keys() {
    let dir = TempDir::new();
    let datetime = "2026-01-01T00:00:00.000000000000Z";
    let mut writer = SigmfWriter::create(
        &dir.0,
        "integer",
        SampleFormat::Ci16Le,
        global(Some("Consumer recording")),
        vec![capture(0, Some(datetime)), capture(3, None)],
    )
    .unwrap();
    let iq = [i16::MIN, -123, 0, 1, 0x1234, i16::MAX, -1, 17];
    writer.write_ci16(&iq[..2]).unwrap();
    writer.write_ci16(&iq[2..6]).unwrap();
    writer.write_ci16(&[]).unwrap();
    writer.write_ci16(&iq[6..]).unwrap();
    assert!(!dir.0.join("integer.sigmf-meta").exists());
    let summary = writer.finish().unwrap();
    let recording = read_recording(&dir.0, "integer").unwrap();
    let expected: Vec<_> = iq.iter().flat_map(|value| value.to_le_bytes()).collect();
    assert_eq!(recording.data, expected);
    assert_eq!(decode_ci16(&recording.data).unwrap(), iq);
    assert_summary(&summary, &recording.data, 4, SampleFormat::Ci16Le);
    assert_eq!(summary.data_path, dir.0.join("integer.sigmf-data"));
    assert_eq!(summary.meta_path, dir.0.join("integer.sigmf-meta"));

    let meta = &recording.meta;
    assert_keys(meta, &["annotations", "captures", "global"]);
    assert_eq!(meta["annotations"], serde_json::json!([]));
    assert_keys(
        &meta["global"],
        &[
            "core:datatype",
            "core:description",
            "core:num_channels",
            "core:recorder",
            "core:sample_rate",
            "core:version",
        ],
    );
    assert_eq!(meta["global"]["core:datatype"], "ci16_le");
    assert_eq!(meta["global"]["core:sample_rate"], 10_000_000_f64 / 3.0);
    assert_eq!(meta["global"]["core:version"], "1.2.0");
    assert_eq!(meta["global"]["core:num_channels"], 1);
    assert_eq!(meta["global"]["core:recorder"], "iqsim-io test");
    assert_eq!(meta["global"]["core:description"], "Consumer recording");
    assert_eq!(meta["captures"].as_array().unwrap().len(), 2);
    assert_keys(
        &meta["captures"][0],
        &["core:datetime", "core:frequency", "core:sample_start"],
    );
    assert_keys(
        &meta["captures"][1],
        &["core:frequency", "core:sample_start"],
    );
    assert_eq!(meta["captures"][0]["core:datetime"], datetime);
    assert_eq!(meta["captures"][0]["core:sample_start"], 0);
    assert_eq!(meta["captures"][1]["core:sample_start"], 3);
    for context in meta["captures"].as_array().unwrap() {
        assert_eq!(context["core:frequency"], 1_090_000_000.25);
    }
    let meta_bytes = fs::read(&summary.meta_path).unwrap();
    assert_eq!(meta_bytes.last(), Some(&b'\n'));
    let mut expected_meta = serde_json::to_vec_pretty(meta).unwrap();
    expected_meta.push(b'\n');
    assert_eq!(meta_bytes, expected_meta);
}

#[test]
fn round_trip_cf32_in_several_chunks_preserves_bits() {
    let dir = TempDir::new();
    let mut writer = SigmfWriter::create(
        &dir.0,
        "float",
        SampleFormat::Cf32Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    let iq = [0.0_f32, -0.0, f32::MAX, f32::MIN, f32::from_bits(1), -0.125];
    for chunk in iq.chunks(2) {
        writer.write_cf32(chunk).unwrap();
    }
    let summary = writer.finish().unwrap();
    let recording = read_recording(&dir.0, "float").unwrap();
    let expected: Vec<_> = iq.iter().flat_map(|value| value.to_le_bytes()).collect();
    assert_eq!(recording.data, expected);
    let decoded = decode_cf32(&recording.data).unwrap();
    assert_eq!(
        decoded
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        iq.map(f32::to_bits)
    );
    assert_summary(&summary, &recording.data, 3, SampleFormat::Cf32Le);
    assert_keys(&recording.meta, &["annotations", "captures", "global"]);
    assert_keys(
        &recording.meta["global"],
        &[
            "core:datatype",
            "core:num_channels",
            "core:recorder",
            "core:sample_rate",
            "core:version",
        ],
    );
    assert_eq!(recording.meta["global"]["core:datatype"], "cf32_le");
    assert_keys(
        &recording.meta["captures"][0],
        &["core:frequency", "core:sample_start"],
    );
    assert_eq!(SampleFormat::Ci16Le.bytes_per_sample(), 4);
    assert_eq!(SampleFormat::Cf32Le.bytes_per_sample(), 8);
}

#[test]
fn rejected_input_writes_nothing_and_valid_writes_can_follow() {
    let dir = TempDir::new();
    let mut writer = SigmfWriter::create(
        &dir.0,
        "integer",
        SampleFormat::Ci16Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    let error = writer.write_ci16(&[1, 2, 3]).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::InvalidInput);
    assert_eq!(
        error.path(),
        Some(dir.0.join("integer.sigmf-data").as_path())
    );
    assert!(writer.write_cf32(&[1.0, 2.0]).is_err());
    writer.write_ci16(&[7, 8]).unwrap();
    let summary = writer.finish().unwrap();
    assert_summary(
        &summary,
        &fs::read(&summary.data_path).unwrap(),
        1,
        SampleFormat::Ci16Le,
    );

    let mut writer = SigmfWriter::create(
        &dir.0,
        "float",
        SampleFormat::Cf32Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    for value in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        assert!(writer.write_cf32(&[1.0, 2.0, 3.0, value]).is_err());
    }
    assert!(writer.write_cf32(&[1.0]).is_err());
    assert!(writer.write_ci16(&[1, 2]).is_err());
    writer.write_cf32(&[0.5, -0.5]).unwrap();
    let summary = writer.finish().unwrap();
    let bytes = fs::read(&summary.data_path).unwrap();
    assert_eq!(decode_cf32(&bytes).unwrap(), [0.5, -0.5]);
    assert_summary(&summary, &bytes, 1, SampleFormat::Cf32Le);
}

#[test]
fn existing_data_or_metadata_is_never_overwritten() {
    let dir = TempDir::new();
    for suffix in ["sigmf-data", "sigmf-meta"] {
        let path = dir.0.join(format!("existing.{suffix}"));
        fs::write(&path, b"untouched").unwrap();
        let error = SigmfWriter::create(
            &dir.0,
            "existing",
            SampleFormat::Ci16Le,
            global(None),
            vec![capture(0, None)],
        )
        .unwrap_err();
        assert_eq!(error.kind(), ErrorKind::AlreadyExists);
        assert_eq!(error.path(), Some(path.as_path()));
        assert_eq!(fs::read(&path).unwrap(), b"untouched");
        assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 1);
        fs::remove_file(path).unwrap();
    }
}

#[test]
fn invalid_stems_and_captures_create_no_files() {
    let dir = TempDir::new();
    for stem in ["", ".", "..", ".hidden", "a/b", "a\\b", "/absolute", "a\0b"] {
        assert!(
            SigmfWriter::create(
                &dir.0,
                stem,
                SampleFormat::Ci16Le,
                global(None),
                vec![capture(0, None)],
            )
            .is_err(),
            "{stem:?}"
        );
        assert!(read_recording(&dir.0, stem).is_err());
    }
    let mut invalid = vec![
        vec![],
        vec![capture(1, None)],
        vec![capture(0, None), capture(2, None), capture(1, None)],
    ];
    for frequency_hz in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut context = capture(0, None);
        context.frequency_hz = frequency_hz;
        invalid.push(vec![context]);
    }
    for captures in invalid {
        assert!(
            SigmfWriter::create(
                &dir.0,
                "invalid",
                SampleFormat::Cf32Le,
                global(None),
                captures,
            )
            .is_err()
        );
    }
    assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 0);
}

#[test]
fn unfinished_recording_has_no_metadata() {
    let dir = TempDir::new();
    let mut writer = SigmfWriter::create(
        &dir.0,
        "unfinished",
        SampleFormat::Ci16Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    writer.write_ci16(&[1, 2]).unwrap();
    drop(writer);
    assert!(dir.0.join("unfinished.sigmf-data").exists());
    assert!(!dir.0.join("unfinished.sigmf-meta").exists());
    assert!(read_recording(&dir.0, "unfinished").is_err());
}

#[test]
fn finish_rejects_length_mismatch_without_metadata() {
    let dir = TempDir::new();
    let mut writer = SigmfWriter::create(
        &dir.0,
        "damaged",
        SampleFormat::Ci16Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    writer.write_ci16(&[1, 2]).unwrap();
    let path = dir.0.join("damaged.sigmf-data");
    fs::OpenOptions::new()
        .write(true)
        .open(&path)
        .unwrap()
        .set_len(99)
        .unwrap();
    let error = writer.finish().unwrap_err();
    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert_eq!(error.path(), Some(path.as_path()));
    assert!(error.to_string().contains("expected 4"));
    assert!(!dir.0.join("damaged.sigmf-meta").exists());
}

#[test]
fn metadata_created_before_finish_is_not_overwritten() {
    let dir = TempDir::new();
    let writer = SigmfWriter::create(
        &dir.0,
        "race",
        SampleFormat::Ci16Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    let path = dir.0.join("race.sigmf-meta");
    fs::write(&path, b"untouched").unwrap();
    let error = writer.finish().unwrap_err();
    assert_eq!(error.kind(), ErrorKind::AlreadyExists);
    assert_eq!(error.path(), Some(path.as_path()));
    assert_eq!(fs::read(path).unwrap(), b"untouched");
}

#[test]
fn empty_recording_can_be_finalized() {
    let dir = TempDir::new();
    let writer = SigmfWriter::create(
        &dir.0,
        "empty",
        SampleFormat::Cf32Le,
        global(None),
        vec![capture(0, None)],
    )
    .unwrap();
    let summary = writer.finish().unwrap();
    assert_summary(&summary, &[], 0, SampleFormat::Cf32Le);
}

#[test]
fn read_and_decode_failures_are_reported() {
    let dir = TempDir::new();
    let error = read_recording(&dir.0, "missing").unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert_eq!(
        error.path(),
        Some(dir.0.join("missing.sigmf-meta").as_path())
    );
    fs::write(dir.0.join("bad.sigmf-meta"), b"not JSON").unwrap();
    let error = read_recording(&dir.0, "bad").unwrap_err();
    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert!(error.to_string().contains("invalid SigMF JSON"));
    fs::write(dir.0.join("bad.sigmf-meta"), b"{}").unwrap();
    let error = read_recording(&dir.0, "bad").unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert_eq!(error.path(), Some(dir.0.join("bad.sigmf-data").as_path()));
    assert!(decode_ci16(&[0]).is_err());
    for length in [1, 2, 3, 5, 6, 7] {
        assert!(decode_cf32(&vec![0; length]).is_err());
    }
    assert!(decode_ci16(&[]).unwrap().is_empty());
    assert!(decode_cf32(&[]).unwrap().is_empty());
    // Helpers decode elements; an unpaired single element is still valid.
    assert_eq!(decode_ci16(&[0x34, 0x12]).unwrap(), [0x1234]);
    assert_eq!(decode_cf32(&1.0_f32.to_le_bytes()).unwrap(), [1.0]);
}

#[test]
fn fixed_recording_golden_bytes() {
    let dir = TempDir::new();
    let mut metadata = global(Some("Fixed SigMF test"));
    metadata.sample_rate = SampleRate::from_ratio(5, 2).unwrap();
    let mut writer = SigmfWriter::create(
        &dir.0,
        "golden",
        SampleFormat::Ci16Le,
        metadata,
        vec![
            capture(0, Some("2026-01-01T00:00:00.000000000000Z")),
            capture(2, None),
        ],
    )
    .unwrap();
    writer.write_ci16(&[i16::MIN, -1]).unwrap();
    writer.write_ci16(&[0, 1, 0x1234, i16::MAX]).unwrap();
    let summary = writer.finish().unwrap();
    iqsim_testkit::assert_golden("io/sigmf_meta", &fs::read(summary.meta_path).unwrap());
    iqsim_testkit::assert_golden("io/sigmf_data", &fs::read(summary.data_path).unwrap());
}
