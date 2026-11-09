//! Property tests for Philox keying, ln and Gaussian noise (prompt 3).

use iqsim_core::{
    KEY_DOMAIN, NoiseStream, PhiloxStream, SampleIndex, StreamId, ln, mix64, noise, philox4x32_10,
    stream_key,
};
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;

/// Cases per property when `PROPTEST_CASES` is not set.
const DEFAULT_CASES: u32 = 2048;

/// Proptest settings, as in `numeric_props.rs` (explicit persistence path, because
/// integration tests have no lib.rs beside them).
fn config() -> ProptestConfig {
    let mut config = ProptestConfig {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "proptest-regressions/noise_props.txt",
        ))),
        ..ProptestConfig::default()
    };
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = DEFAULT_CASES;
    }
    config
}

/// Multiplicative inverse of an odd u64 mod 2⁶⁴ (Newton iteration).
fn inv_odd(a: u64) -> u64 {
    let mut x = a; // correct to 3 bits
    for _ in 0..5 {
        x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
    }
    x
}

/// Inverse of `z ^= z >> s`.
fn unxorshift(z: u64, s: u32) -> u64 {
    let mut x = z;
    let mut k = s;
    while k < 64 {
        x = z ^ (x >> s);
        k += s;
    }
    x
}

/// Inverse of mix64, built from the inverse of each step.
fn unmix64(z: u64) -> u64 {
    let z = unxorshift(z, 31);
    let z = z.wrapping_mul(inv_odd(0x94D0_49BB_1331_11EB));
    let z = unxorshift(z, 27);
    let z = z.wrapping_mul(inv_odd(0xBF58_476D_1CE4_E5B9));
    unxorshift(z, 30)
}

/// f64 test oracle for ln (std ln is allowed in tests only as a reference).
#[allow(clippy::disallowed_methods)]
fn ln_ref(x: f32) -> f64 {
    (x as f64).ln()
}

/// One f32 ulp at |v| (v finite).
fn ulp_at(v: f64) -> f64 {
    let a = (v as f32).abs();
    if a == 0.0 {
        f32::from_bits(1) as f64
    } else {
        (f32::from_bits(a.to_bits() + 1) - a) as f64
    }
}

fn positive_finite_f32() -> impl Strategy<Value = f32> {
    (1u32..0x7F80_0000).prop_map(f32::from_bits)
}

proptest! {
    #![proptest_config(config())]

    /// mix64 is a bijection: an explicit inverse undoes it.
    #[test]
    fn mix64_is_invertible(z in any::<u64>()) {
        prop_assert_eq!(unmix64(mix64(z)), z);
        prop_assert_eq!(mix64(unmix64(z)), z);
    }

    /// Within one seed, the key determines (entity, stream): unmixing recovers them.
    #[test]
    fn stream_key_is_injective_within_a_seed(
        seed in any::<u64>(), entity in any::<u32>(), instance in any::<u16>(),
    ) {
        let id = StreamId::receiver_thermal_noise(instance);
        let k = stream_key(seed, entity, id);
        let k = (k[1] as u64) << 32 | k[0] as u64;
        let packed = unmix64(k) ^ mix64(seed ^ KEY_DOMAIN);
        prop_assert_eq!(packed, (entity as u64) << 32 | id.to_bits() as u64);
    }

    /// Distinct (entity, stream) pairs under one seed give distinct keys.
    #[test]
    fn distinct_pairs_distinct_keys(
        seed in any::<u64>(),
        a in (any::<u32>(), any::<u16>()),
        b in (any::<u32>(), any::<u16>()),
    ) {
        prop_assume!(a != b);
        let ka = stream_key(seed, a.0, StreamId::receiver_thermal_noise(a.1));
        let kb = stream_key(seed, b.0, StreamId::receiver_thermal_noise(b.1));
        prop_assert_ne!(ka, kb);
    }

    /// The stream counter layout and the free noise function agree with the raw block
    /// function, at any index.
    #[test]
    fn noise_is_a_pure_function_of_its_inputs(
        seed in any::<u64>(), entity in any::<u32>(), instance in any::<u16>(), n in any::<u64>(),
    ) {
        let id = StreamId::receiver_thermal_noise(instance);
        let key = stream_key(seed, entity, id);
        let block = philox4x32_10([n as u32, (n >> 32) as u32, 0, 0], key);
        prop_assert_eq!(PhiloxStream::new(seed, entity, id).block(n, 0), block);
        let idx = SampleIndex::new(n);
        let a = noise(seed, entity, id, idx);
        let b = NoiseStream::new(seed, entity, id).sample(idx);
        prop_assert_eq!(a.0.to_bits(), b.0.to_bits());
        prop_assert_eq!(a.1.to_bits(), b.1.to_bits());
    }

    /// Every noise sample is finite and inside the radius bound sqrt(65·ln 2) ≈ 6.71.
    #[test]
    fn noise_is_finite_and_bounded(
        seed in any::<u64>(), entity in any::<u32>(), instance in any::<u16>(), n in any::<u64>(),
    ) {
        let (i, q) = noise(seed, entity, StreamId::receiver_thermal_noise(instance), SampleIndex::new(n));
        prop_assert!(i.is_finite() && q.is_finite());
        let r2 = (i as f64) * (i as f64) + (q as f64) * (q as f64);
        prop_assert!(r2 < 6.72 * 6.72, "r² = {}", r2);
    }

    /// ln is within 1 ulp of the f64 reference for any positive finite f32.
    #[test]
    fn ln_within_one_ulp(x in positive_finite_f32()) {
        let want = ln_ref(x);
        let got = ln(x).unwrap() as f64;
        prop_assert!((got - want).abs() <= ulp_at(want), "x = {:e}: {} vs {}", x, got, want);
    }

    /// ln is monotone non-decreasing between neighbouring f32 values.
    #[test]
    fn ln_is_monotone(bits in 1u32..0x7F7F_FFFF) {
        let a = ln(f32::from_bits(bits)).unwrap();
        let b = ln(f32::from_bits(bits + 1)).unwrap();
        prop_assert!(a <= b, "{:#x}: {} > {}", bits, a, b);
    }
}

/// Monotonicity over every f32 in [0.5, 2), where the range reduction switches at √2.
#[test]
fn ln_monotone_around_one() {
    let lo = 0.5f32.to_bits();
    let hi = 2.0f32.to_bits();
    let mut prev = ln(f32::from_bits(lo)).unwrap();
    for bits in lo + 1..hi {
        let v = ln(f32::from_bits(bits)).unwrap();
        assert!(v >= prev, "{bits:#x}");
        prev = v;
    }
}
