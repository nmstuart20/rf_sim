//! SigMF emitter-input validation, streaming statistics, and verified paging.

#![forbid(unsafe_code)]

use std::fs::{self, OpenOptions};
use std::io::{ErrorKind, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use iqsim_io::IoError;
use iqsim_io::iq_input::{IqDatatype, IqInput, PAGE_BYTES};
use serde_json::{Value, json};
use sha2::{Digest, Sha256, Sha512};

static NEXT_DIR: AtomicU64 = AtomicU64::new(0);

struct Fixture {
    directory: PathBuf,
    meta: PathBuf,
    data: PathBuf,
}

impl Fixture {
    fn new(meta: &Value, data: &[u8]) -> Self {
        let directory = loop {
            let count = NEXT_DIR.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir()
                .join(format!("iqsim-io-iq-input-{}-{count}", std::process::id()));
            match fs::create_dir(&path) {
                Ok(()) => break path,
                Err(error) if error.kind() == ErrorKind::AlreadyExists => continue,
                Err(error) => panic!("creating {}: {error}", path.display()),
            }
        };
        let fixture = Self {
            meta: directory.join("input.sigmf-meta"),
            data: directory.join("input.sigmf-data"),
            directory,
        };
        fs::write(&fixture.meta, serde_json::to_vec_pretty(meta).unwrap()).unwrap();
        fs::write(&fixture.data, data).unwrap();
        fixture
    }

    fn open(&self) -> IqInput {
        IqInput::open(&self.meta, &self.data).unwrap()
    }

    fn invalid(&self, field: &str, metadata: bool) {
        let error = IqInput::open(&self.meta, &self.data).unwrap_err();
        assert!(matches!(error, IoError::InvalidInput { .. }), "{error}");
        assert!(error.to_string().contains(field), "{error}");
        assert_eq!(
            error.path(),
            Some(if metadata { &self.meta } else { &self.data }.as_path())
        );
    }

    fn overwrite_byte(&self, offset: u64) {
        let mut file = OpenOptions::new().write(true).open(&self.data).unwrap();
        file.seek(SeekFrom::Start(offset)).unwrap();
        file.write_all(&[0xff]).unwrap();
        file.sync_all().unwrap();
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.directory).expect("clean up IQ input fixture");
    }
}

fn meta(datatype: &str) -> Value {
    json!({
        "global": {
            "core:datatype": datatype,
            "core:sample_rate": 2_500_000.0,
            "core:version": "1.2.0",
            "core:num_channels": 1
        },
        "captures": [{"core:sample_start": 0, "core:frequency": 1_090_000_000.25}],
        "annotations": []
    })
}

fn cf32_bytes(samples: &[[f32; 2]]) -> Vec<u8> {
    samples
        .iter()
        .flatten()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

fn assert_changed(error: &IoError, fixture: &Fixture) {
    assert!(matches!(error, IoError::Io { .. }), "{error}");
    assert_eq!(error.kind(), ErrorKind::InvalidData);
    assert!(
        error.to_string().contains("changed since it was opened"),
        "{error}"
    );
    assert_eq!(error.path(), Some(fixture.data.as_path()));
}

#[test]
fn ci16_round_trip_hashes_and_exact_statistics() {
    let values = [[i16::MIN, 0], [0, 16384], [24576, i16::MIN]];
    let bytes: Vec<_> = values
        .iter()
        .flatten()
        .flat_map(|value| value.to_le_bytes())
        .collect();
    let fixture = Fixture::new(&meta("ci16_le"), &bytes);
    let input = fixture.open();
    let metadata = input.meta();
    assert_eq!(metadata.datatype, IqDatatype::Ci16Le);
    assert_eq!(metadata.datatype.name(), "ci16_le");
    assert_eq!(metadata.datatype.bytes_per_sample(), 4);
    assert_eq!(metadata.sample_rate_hz, 2_500_000.0);
    assert_eq!(metadata.version, "1.2.0");
    assert_eq!(metadata.frequency_hz, Some(1_090_000_000.25));
    let meta_bytes = fs::read(&fixture.meta).unwrap();
    assert_eq!(metadata.meta_byte_len, meta_bytes.len() as u64);
    assert_eq!(
        metadata.meta_sha256,
        <[u8; 32]>::from(Sha256::digest(&meta_bytes))
    );
    assert_eq!(input.data_path(), fixture.data);
    let summary = input.summary();
    assert_eq!(summary.sample_count, 3);
    assert_eq!(summary.byte_len, 12);
    assert_eq!(summary.sha256, <[u8; 32]>::from(Sha256::digest(&bytes)));
    assert_eq!(summary.page_sha256, vec![summary.sha256]);
    assert_eq!(summary.peak_magnitude, 1.25);
    assert_eq!(summary.sum_power, 2.8125);
    let mut reader = input.reader(1);
    assert_eq!(reader.cached_pages(), 0);
    let mut samples = [[0.0; 2]; 3];
    reader.read(0, &mut samples).unwrap();
    assert_eq!(samples, [[-1.0, 0.0], [0.0, 0.5], [0.75, -1.0]]);
    assert_eq!(reader.cached_pages(), 1);
    reader.read(1, &mut samples[..1]).unwrap();
    assert_eq!(samples[0], [0.0, 0.5]);
}

#[test]
fn cf32_round_trip_preserves_widened_values_and_signs() {
    let samples = [
        [0.75, -1.0],
        [0.0, -0.0],
        [f32::MAX, f32::MIN],
        [f32::from_bits(1), -f32::from_bits(1)],
        [0.125, -0.5],
    ];
    let bytes = cf32_bytes(&samples);
    let fixture = Fixture::new(&meta("cf32_le"), &bytes);
    let input = fixture.open();
    assert_eq!(input.meta().datatype, IqDatatype::Cf32Le);
    assert_eq!(input.meta().datatype.name(), "cf32_le");
    assert_eq!(input.meta().datatype.bytes_per_sample(), 8);
    assert_eq!(input.summary().sample_count, samples.len() as u64);
    assert_eq!(input.summary().byte_len, bytes.len() as u64);
    assert_eq!(
        input.summary().sha256,
        <[u8; 32]>::from(Sha256::digest(&bytes))
    );
    let mut reader = input.reader(1);
    drop(input);
    let mut actual = [[0.0; 2]; 5];
    reader.read(0, &mut actual).unwrap();
    for (actual, expected) in actual.iter().flatten().zip(samples.iter().flatten()) {
        assert_eq!(actual.to_bits(), f64::from(*expected).to_bits());
    }
}

#[test]
fn cf32_statistics_accumulate_power_in_file_order() {
    // The unit powers disappear when accumulated after 2^60. Reordering them
    // would produce a different sum, so this also fixes the operation order.
    let mut samples = vec![[1.0, 0.0]; 257];
    samples[0] = [1_073_741_824.0, 1.0];
    let fixture = Fixture::new(&meta("cf32_le"), &cf32_bytes(&samples));
    let input = fixture.open();
    assert_eq!(input.summary().sum_power, 1_152_921_504_606_846_976.0);
    assert_eq!(input.summary().peak_magnitude, 1_073_741_824.0);

    let fixture = Fixture::new(
        &meta("cf32_le"),
        &cf32_bytes(&[[0.75, -1.0], [-1.0, 0.0], [0.0, 0.5]]),
    );
    let input = fixture.open();
    assert_eq!(input.summary().peak_magnitude, 1.25);
    assert_eq!(input.summary().sum_power, 2.8125);
}

#[test]
fn unsupported_datatypes_are_named() {
    for datatype in ["ci16_be", "ri16_le", "cf64_le", "ci8", "CI16_LE", ""] {
        let fixture = Fixture::new(&meta(datatype), &[0; 4]);
        fixture.invalid("unsupported core:datatype", true);
    }
}

#[test]
fn metadata_structure_and_field_types_are_validated() {
    let cases = [
        ("/global", Some(json!([])), "global"),
        ("/global", None, "global"),
        ("/captures", Some(json!({})), "captures"),
        ("/captures", None, "captures"),
        ("/annotations", Some(json!(null)), "annotations"),
        ("/global/core:datatype", None, "core:datatype"),
        ("/global/core:datatype", Some(json!(16)), "core:datatype"),
        ("/global/core:sample_rate", None, "core:sample_rate"),
        (
            "/global/core:sample_rate",
            Some(json!(-1)),
            "core:sample_rate",
        ),
        (
            "/global/core:sample_rate",
            Some(json!(0)),
            "core:sample_rate",
        ),
        (
            "/global/core:sample_rate",
            Some(json!("NaN")),
            "core:sample_rate",
        ),
        (
            "/global/core:sample_rate",
            Some(json!(null)),
            "core:sample_rate",
        ),
        ("/global/core:version", None, "core:version"),
        ("/global/core:version", Some(json!(12)), "core:version"),
        (
            "/global/core:num_channels",
            Some(json!(2)),
            "core:num_channels",
        ),
        (
            "/global/core:num_channels",
            Some(json!(1.0)),
            "core:num_channels",
        ),
        (
            "/global/core:num_channels",
            Some(json!("1")),
            "core:num_channels",
        ),
        ("/captures/0", Some(json!([])), "captures[0]"),
        (
            "/captures/0/core:sample_start",
            Some(json!(-1)),
            "core:sample_start",
        ),
        (
            "/captures/0/core:sample_start",
            Some(json!(0.5)),
            "core:sample_start",
        ),
        (
            "/captures/0/core:sample_start",
            Some(json!("0")),
            "core:sample_start",
        ),
        (
            "/captures/0/core:frequency",
            Some(json!("Inf")),
            "core:frequency",
        ),
        (
            "/captures/0/core:frequency",
            Some(json!(null)),
            "core:frequency",
        ),
    ];
    for (pointer, replacement, field) in cases {
        let mut metadata = meta("ci16_le");
        if let Some(replacement) = replacement {
            *metadata.pointer_mut(pointer).unwrap() = replacement;
        } else {
            let (parent, key) = pointer.rsplit_once('/').unwrap();
            metadata
                .pointer_mut(parent)
                .unwrap()
                .as_object_mut()
                .unwrap()
                .remove(key);
        }
        let fixture = Fixture::new(&metadata, &[0; 4]);
        fixture.invalid(field, true);
    }
}

#[test]
fn malformed_json_and_nonobject_root_are_rejected() {
    let fixture = Fixture::new(&meta("ci16_le"), &[0; 4]);
    for bytes in [b"{".as_slice(), b"not json", b"{\"global\": NaN}"] {
        fs::write(&fixture.meta, bytes).unwrap();
        fixture.invalid("JSON", true);
    }
    fs::write(&fixture.meta, b"[]").unwrap();
    fixture.invalid("top-level", true);
}

#[test]
fn optional_fields_empty_captures_and_capture_order() {
    let mut metadata = meta("ci16_le");
    metadata.as_object_mut().unwrap().remove("annotations");
    metadata["global"]
        .as_object_mut()
        .unwrap()
        .remove("core:num_channels");
    metadata["captures"] = json!([]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    assert_eq!(fixture.open().meta().frequency_hz, None);

    metadata["captures"] = json!([{}, {"core:frequency": -12.0}]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    assert_eq!(fixture.open().meta().frequency_hz, None);

    metadata["captures"] = json!([
        {"core:sample_start": 0},
        {"core:sample_start": 0},
        {"core:sample_start": 2}
    ]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    fixture.open();

    metadata["captures"] = json!([{"core:sample_start": 2}, {}, {"core:sample_start": 1}]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    fixture.invalid("non-decreasing", true);

    metadata["captures"] = json!([{"core:sample_start": 3}]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    fixture.invalid("core:sample_start", true);

    metadata["captures"] = json!([{}, {"core:frequency": "NaN"}]);
    let fixture = Fixture::new(&metadata, &[0; 12]);
    fixture.invalid("captures[1].core:frequency", true);
}

#[test]
fn sha512_is_validated_against_exact_data_bytes() {
    let bytes = cf32_bytes(&[[0.75, -1.0], [-0.0, 0.5]]);
    let mut metadata = meta("cf32_le");
    metadata["global"]["core:sha512"] = json!(hex::encode_upper(Sha512::digest(&bytes)));
    let fixture = Fixture::new(&metadata, &bytes);
    fixture.open();

    for invalid in [
        json!("0".repeat(128)),
        json!("z".repeat(128)),
        json!("abc"),
        json!(512),
    ] {
        metadata["global"]["core:sha512"] = invalid;
        let fixture = Fixture::new(&metadata, &bytes);
        fixture.invalid("core:sha512", true);
    }
    metadata["global"]["core:sha512"] = json!("0".repeat(128));
    let fixture = Fixture::new(&metadata, &bytes);
    fixture.invalid("core:sha512 does not match the data file", true);
}

#[test]
fn empty_and_truncated_data_are_rejected_for_both_datatypes() {
    for datatype in [IqDatatype::Ci16Le, IqDatatype::Cf32Le] {
        let fixture = Fixture::new(&meta(datatype.name()), &[]);
        fixture.invalid("empty", false);
        for len in [
            1,
            datatype.bytes_per_sample() - 1,
            datatype.bytes_per_sample() + 1,
        ] {
            let fixture = Fixture::new(&meta(datatype.name()), &vec![0; len as usize]);
            fixture.invalid("truncated", false);
        }
    }
}

#[test]
fn cf32_nonfinite_values_report_first_bad_sample_index() {
    for bad in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        for sample in [[bad, 0.0], [0.0, bad]] {
            let fixture = Fixture::new(
                &meta("cf32_le"),
                &cf32_bytes(&[[0.0, 0.0], sample, [f32::NAN, f32::NAN]]),
            );
            fixture.invalid("sample index 1", false);
        }
    }
    let mut bytes = vec![0; PAGE_BYTES as usize];
    bytes.extend_from_slice(&cf32_bytes(&[[0.0, f32::NAN]]));
    let fixture = Fixture::new(&meta("cf32_le"), &bytes);
    fixture.invalid(&format!("sample index {}", PAGE_BYTES / 8), false);
}

fn page_sample(index: u64) -> [f32; 2] {
    [
        (index % 8192) as f32 / 1024.0 - 4.0,
        ((index * 3) % 4096) as f32 / 512.0 - 4.0,
    ]
}

#[test]
fn multipage_reads_and_evictions_are_independent_of_cache_capacity() {
    let per_page = PAGE_BYTES / 8;
    let count = 3 * per_page + 17;
    let bytes: Vec<_> = (0..count)
        .flat_map(|index| page_sample(index).into_iter().flat_map(f32::to_le_bytes))
        .collect();
    let mut metadata = meta("cf32_le");
    metadata["global"]["core:sha512"] = json!(hex::encode(Sha512::digest(&bytes)));
    let fixture = Fixture::new(&metadata, &bytes);
    let input = fixture.open();
    assert_eq!(input.summary().sample_count, count);
    assert_eq!(input.summary().byte_len, bytes.len() as u64);
    assert_eq!(
        input.summary().sha256,
        <[u8; 32]>::from(Sha256::digest(&bytes))
    );
    let hashes: Vec<[u8; 32]> = bytes
        .chunks(PAGE_BYTES as usize)
        .map(|page| Sha256::digest(page).into())
        .collect();
    assert_eq!(input.summary().page_sha256, hashes);
    assert_eq!(hashes.len(), 4);
    let mut one = input.reader(1);
    let mut two = input.reader(2);
    for (first, len) in [
        (per_page - 3, 9),
        (2 * per_page - 3, 9),
        (3 * per_page - 3, 20),
        (0, 5),
        (per_page - 1, per_page as usize + 3),
        (count - 1, 1),
        (0, count as usize),
    ] {
        let mut a = vec![[0.0; 2]; len];
        let mut b = vec![[0.0; 2]; len];
        one.read(first, &mut a).unwrap();
        two.read(first, &mut b).unwrap();
        assert_eq!(a, b);
        for (offset, actual) in a.iter().enumerate() {
            assert_eq!(*actual, page_sample(first + offset as u64).map(f64::from));
        }
        assert!(one.cached_pages() <= 1);
        assert!(two.cached_pages() <= 2);
    }
}

#[test]
fn changed_page_is_rejected_on_initial_load_and_after_eviction() {
    let fixture = Fixture::new(&meta("cf32_le"), &vec![0; 2 * PAGE_BYTES as usize + 8]);
    let input = fixture.open();
    let mut reader = input.reader(1);
    let mut out = [[0.0; 2]; 1];
    reader.read(0, &mut out).unwrap();
    fixture.overwrite_byte(0);
    // Cached bytes remain validated until eviction.
    reader.read(0, &mut out).unwrap();
    assert_eq!(out, [[0.0; 2]; 1]);
    reader.read(PAGE_BYTES / 8, &mut out).unwrap();
    assert_changed(&reader.read(0, &mut out).unwrap_err(), &fixture);
    assert!(reader.cached_pages() <= 1);
    assert_changed(&input.reader(2).read(0, &mut out).unwrap_err(), &fixture);
}

#[test]
fn cache_hits_update_lru_order() {
    let fixture = Fixture::new(&meta("ci16_le"), &vec![0; 3 * PAGE_BYTES as usize]);
    let input = fixture.open();
    let mut reader = input.reader(2);
    let mut out = [[0.0; 2]; 1];
    let per_page = PAGE_BYTES / 4;
    for first in [0, per_page, 0, 2 * per_page] {
        reader.read(first, &mut out).unwrap();
        assert!(reader.cached_pages() <= 2);
    }
    fixture.overwrite_byte(0);
    fixture.overwrite_byte(PAGE_BYTES);
    reader.read(0, &mut out).unwrap();
    assert_eq!(out, [[0.0; 2]; 1]);
    // Page zero was refreshed, so page one must have been evicted.
    assert_changed(&reader.read(per_page, &mut out).unwrap_err(), &fixture);
}

#[test]
fn changed_file_length_is_rejected() {
    for new_len in [7, 9] {
        let fixture = Fixture::new(&meta("cf32_le"), &[0; 8]);
        let input = fixture.open();
        OpenOptions::new()
            .write(true)
            .open(&fixture.data)
            .unwrap()
            .set_len(new_len)
            .unwrap();
        assert_changed(
            &input.reader(1).read(0, &mut [[0.0; 2]; 1]).unwrap_err(),
            &fixture,
        );
    }
}

#[test]
fn out_of_range_and_overflow_reads_are_rejected_without_loading_pages() {
    let fixture = Fixture::new(&meta("ci16_le"), &[0; 8]);
    let mut reader = fixture.open().reader(1);
    for (first, len) in [(2, 1), (1, 2), (3, 0), (u64::MAX, 1), (u64::MAX - 1, 3)] {
        let mut out = vec![[9.0; 2]; len];
        let error = reader.read(first, &mut out).unwrap_err();
        assert!(matches!(error, IoError::InvalidInput { .. }), "{error}");
        assert_eq!(error.path(), Some(fixture.data.as_path()));
        assert_eq!(out, vec![[9.0; 2]; len]);
        assert_eq!(reader.cached_pages(), 0);
    }
    reader.read(0, &mut []).unwrap();
    reader.read(2, &mut []).unwrap();
    assert_eq!(reader.cached_pages(), 0);
}

#[test]
#[should_panic(expected = "max_pages must be at least 1")]
fn zero_cache_capacity_is_rejected() {
    let fixture = Fixture::new(&meta("ci16_le"), &[0; 4]);
    fixture.open().reader(0);
}

#[test]
fn file_access_errors_retain_the_associated_path() {
    let fixture = Fixture::new(&meta("ci16_le"), &[0; 4]);
    let input = fixture.open();
    fs::remove_file(&fixture.data).unwrap();
    let error = IqInput::open(&fixture.meta, &fixture.data).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert_eq!(error.path(), Some(fixture.data.as_path()));
    let error = input.reader(1).read(0, &mut [[0.0; 2]; 1]).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert_eq!(error.path(), Some(fixture.data.as_path()));
    fs::remove_file(&fixture.meta).unwrap();
    let error = IqInput::open(&fixture.meta, &fixture.data).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert_eq!(error.path(), Some(fixture.meta.as_path()));
}
