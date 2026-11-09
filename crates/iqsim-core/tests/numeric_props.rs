//! Property and edge-case tests for the numeric core (prompt 2): 128-bit primitives,
//! block polynomials, carrier-offset phase and sin/cos.

use iqsim_core::wide::{U128, add_with_carry, mul_hi_u64, sub_with_borrow};
use iqsim_core::{
    AcqInstant, Delay, DelayPoly, MAX_BLOCK_LEN, OffsetFrequency, PS_PER_SECOND, Phase, PhasePoly,
    PolyOutput, TruthTime, binomials, sincos,
};
use proptest::prelude::*;
use proptest::test_runner::FileFailurePersistence;

/// Cases per property when `PROPTEST_CASES` is not set.
const DEFAULT_CASES: u32 = 2048;

/// Proptest settings: `PROPTEST_CASES` overrides the case count; other settings
/// (including `PROPTEST_RNG_SEED`) come from proptest's environment defaults.
fn config() -> ProptestConfig {
    let mut config = ProptestConfig {
        // Integration tests have no lib.rs beside them, so proptest's default
        // source-relative persistence silently saves nothing. Use an explicit path
        // (relative to the crate directory, where cargo runs tests).
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            "proptest-regressions/numeric_props.txt",
        ))),
        ..ProptestConfig::default()
    };
    if std::env::var_os("PROPTEST_CASES").is_none() {
        config.cases = DEFAULT_CASES;
    }
    config
}

const DAY_PS: i64 = 86_400 * PS_PER_SECOND;

/// SplitMix64: a fixed test-input generator (not simulation randomness).
struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

fn u128s() -> impl Strategy<Value = U128> {
    any::<u128>().prop_map(U128::from_u128)
}

/// Native-u128 reference for the jump formula.
fn jump_ref(d: [U128; 4], k: u64) -> u128 {
    let (c1, c2, c3) = binomials(k);
    d[0].to_u128()
        .wrapping_add(d[1].to_u128().wrapping_mul(u128::from(c1)))
        .wrapping_add(d[2].to_u128().wrapping_mul(u128::from(c2)))
        .wrapping_add(d[3].to_u128().wrapping_mul(u128::from(c3)))
}

/// Check every k of a block: stepping from 0 equals jumping, and jumping equals the
/// native reference.
fn assert_full_block<T: PolyOutput + PartialEq + std::fmt::Debug>(
    p: &iqsim_core::BlockPoly<T>,
    check_reference_every: u64,
) {
    let mut cur = p.cursor(0).unwrap();
    for k in 0..p.block_len() {
        let j = p.jump(k).unwrap();
        assert_eq!(cur.register(), j, "k = {k}");
        assert_eq!(cur.value(), p.eval(k).unwrap());
        if k % check_reference_every == 0 || k + 1 == p.block_len() {
            assert_eq!(j.to_u128(), jump_ref(p.registers(), k), "k = {k}");
        }
        cur.step();
    }
}

proptest! {
    #![proptest_config(config())]

    #[test]
    fn primitives_match_u128(a: u64, b: u64, c in 0u64..2) {
        let full = u128::from(a) * u128::from(b);
        prop_assert_eq!(mul_hi_u64(a, b), (full >> 64) as u64);
        let s = u128::from(a) + u128::from(b) + u128::from(c);
        prop_assert_eq!(add_with_carry(a, b, c), (s as u64, (s >> 64) as u64));
        let d = i128::from(a) - i128::from(b) - i128::from(c);
        prop_assert_eq!(sub_with_borrow(a, b, c), (d as u64, u64::from(d < 0)));
        prop_assert_eq!(U128::mul_wide_u64(a, b).to_u128(), full);
    }

    #[test]
    fn u128_ops_match_native(x in u128s(), y in u128s(), m: u64) {
        let (xn, yn) = (x.to_u128(), y.to_u128());
        prop_assert_eq!(x.wrapping_add(y).to_u128(), xn.wrapping_add(yn));
        prop_assert_eq!(x.wrapping_sub(y).to_u128(), xn.wrapping_sub(yn));
        prop_assert_eq!(x.wrapping_neg().to_u128(), xn.wrapping_neg());
        prop_assert_eq!(x.wrapping_mul_u64(m).to_u128(), xn.wrapping_mul(u128::from(m)));
        // Signed view: the same bits as i128 wrapping arithmetic.
        prop_assert_eq!(
            x.wrapping_mul_u64(m).to_i128(),
            x.to_i128().wrapping_mul(i128::from(m))
        );
        prop_assert_eq!(U128::from_i128(x.to_i128()), x);
    }

    #[test]
    fn jump_matches_reference(d in prop::array::uniform4(u128s()), k in 0..MAX_BLOCK_LEN) {
        let p = PhasePoly::from_registers(d, MAX_BLOCK_LEN).unwrap();
        prop_assert_eq!(p.jump(k).unwrap().to_u128(), jump_ref(d, k));
    }

    /// Stepping from any k0 reaches the same registers as jumping to k1, bit for bit,
    /// for arbitrary (wrapping) registers.
    #[test]
    fn stepping_matches_jumping(
        d in prop::array::uniform4(u128s()),
        k0 in 0..MAX_BLOCK_LEN,
        n in 0u64..3000,
    ) {
        let p = DelayPoly::from_registers(d, MAX_BLOCK_LEN).unwrap();
        let k1 = (k0 + n).min(MAX_BLOCK_LEN - 1);
        let mut cur = p.cursor(k0).unwrap();
        for _ in k0..k1 {
            cur.step();
        }
        prop_assert_eq!(cur, p.cursor(k1).unwrap());
        prop_assert_eq!(cur.register(), p.jump(k1).unwrap());
    }

    /// Exact floor semantics for the carrier-offset phase against an independent
    /// bit-serial big-integer reference.
    #[test]
    fn offset_phase_matches_bignum(f: i64, t: i128) {
        let got = OffsetFrequency::from_bits(f).phase_at(AcqInstant::from_q64_64(t));
        prop_assert_eq!(got.to_bits(), offset_phase_ref(f, t));
    }

    /// φ is additive in t to within the floor: φ(t1+t2) − φ(t1) − φ(t2) ∈ {0, 1}, with
    /// t1 up to ±104 days. Holds at any magnitude, so precision does not degrade with
    /// absolute time.
    #[test]
    fn offset_phase_additive_in_time(
        f: i64,
        t1 in -(104 * i128::from(DAY_PS)) << 64..(104 * i128::from(DAY_PS)) << 64,
        t2 in -(2 * i128::from(DAY_PS)) << 64..(2 * i128::from(DAY_PS)) << 64,
    ) {
        let f = OffsetFrequency::from_bits(f);
        let ph = |t: i128| f.phase_at(AcqInstant::from_q64_64(t)).to_bits();
        let e = ph(t1 + t2).wrapping_sub(ph(t1)).wrapping_sub(ph(t2));
        prop_assert!(e <= 1, "defect {e}");
    }

    /// φ is additive in Δf to within the floor.
    #[test]
    fn offset_phase_additive_in_frequency(
        f1 in i64::MIN / 2..i64::MAX / 2,
        f2 in i64::MIN / 2..i64::MAX / 2,
        t: i128,
    ) {
        let t = AcqInstant::from_q64_64(t);
        let ph = |f: i64| OffsetFrequency::from_bits(f).phase_at(t).to_bits();
        let e = ph(f1 + f2).wrapping_sub(ph(f1)).wrapping_sub(ph(f2));
        prop_assert!(e <= 1, "defect {e}");
    }

    /// Near 100 days, a whole-hertz offset has exactly the phase it has near t = 0,
    /// because 100 days is a whole number of seconds: nothing is lost to magnitude.
    #[test]
    fn offset_phase_exact_near_hundred_days(
        hz in -(1i64 << 31)..(1i64 << 31),
        delta in 0i128..(i128::from(PS_PER_SECOND) << 64),
    ) {
        let f = OffsetFrequency::from_hz(hz).unwrap();
        let big = i128::from(100 * DAY_PS) << 64;
        let near_zero = f.phase_at(AcqInstant::from_q64_64(delta));
        let near_100d = f.phase_at(AcqInstant::from_q64_64(big + delta));
        prop_assert_eq!(near_zero, near_100d);
    }
}

/// Independent reference for `phase_at`: schoolbook product in 32-bit digits, then
/// bit-serial long division by D = 2³²·10¹², then the floor for negative products.
fn offset_phase_ref(f: i64, t: i128) -> u64 {
    let a = f.unsigned_abs();
    let b = t.unsigned_abs();
    let ad = [a as u32, (a >> 32) as u32];
    let bd = [
        b as u32,
        (b >> 32) as u32,
        (b >> 64) as u32,
        (b >> 96) as u32,
    ];
    let mut prod = [0u64; 7]; // 32-bit digits, little-endian, with carry room
    for (i, &x) in ad.iter().enumerate() {
        let mut carry = 0u64;
        for (j, &y) in bd.iter().enumerate() {
            let cur = prod[i + j] + u64::from(x) * u64::from(y) + carry;
            prod[i + j] = cur & 0xFFFF_FFFF;
            carry = cur >> 32;
        }
        prod[i + bd.len()] += carry;
    }
    let d: u128 = (1u128 << 32) * 1_000_000_000_000;
    let mut rem: u128 = 0;
    let mut q: u64 = 0;
    for bit in (0..6 * 32).rev() {
        let digit = prod[bit / 32];
        rem = (rem << 1) | u128::from((digit >> (bit % 32)) & 1);
        if rem >= d {
            rem -= d;
            if bit < 64 {
                q |= 1 << bit;
            }
        }
    }
    let negative = (f < 0) != (t < 0);
    if negative {
        q.wrapping_add(u64::from(rem != 0)).wrapping_neg()
    } else {
        q
    }
}

#[test]
fn full_blocks_bit_for_bit() {
    let mut rng = SplitMix(0x5EED_0002);
    let mut regs = || {
        [
            U128::new(rng.next(), rng.next()),
            U128::new(rng.next(), rng.next()),
            U128::new(rng.next(), rng.next()),
            U128::new(rng.next(), rng.next()),
        ]
    };
    // 1 ms at 61.44 MS/s and at 250 MS/s, arbitrary registers.
    for len in [61_440, 250_000] {
        assert_full_block(&DelayPoly::from_registers(regs(), len).unwrap(), 997);
        assert_full_block(&PhasePoly::from_registers(regs(), len).unwrap(), 997);
    }
    // The longest allowed block.
    assert_full_block(
        &PhasePoly::from_registers(regs(), MAX_BLOCK_LEN).unwrap(),
        10_007,
    );
}

#[test]
fn phase_wraps_many_times_per_block() {
    // 0.37 cycles per sample plus curvature: the phase wraps ~23 000 times in 61 440
    // samples. Stepping and jumping must agree throughout.
    let p = PhasePoly::from_f64_differences([0.9, 0.37, 1e-6, -3e-12], 61_440).unwrap();
    assert_full_block(&p, 1);
    let mut wraps = 0;
    let mut prev = p.eval(0).unwrap();
    for k in 1..61_440 {
        let now = p.eval(k).unwrap();
        if now.to_bits() < prev.to_bits() {
            wraps += 1;
        }
        prev = now;
    }
    assert!(wraps > 20_000, "{wraps} wraps");
}

#[test]
fn negative_delay_rate() {
    // An approaching target: delay 30 000 samples, shrinking by 8.3e-7 samples per
    // sample (≈ 3 600 output LSB), slight deceleration. The output strictly decreases,
    // and stepping matches jumping.
    let p = DelayPoly::from_f64_differences([30_000.0, -8.3e-7, 2e-15, -1e-22], 250_000).unwrap();
    assert_full_block(&p, 1);
    let first = p.eval(0).unwrap();
    let last = p.eval(249_999).unwrap();
    assert!(last < first);
    let mut prev = first;
    for k in 1..250_000 {
        let now = p.eval(k).unwrap();
        assert!(now < prev, "k = {k}");
        prev = now;
    }
}

#[test]
fn delay_crossing_zero_and_negative() {
    // Delay output is signed: a polynomial through zero floors correctly on both sides.
    let p = DelayPoly::from_f64_differences([0.5, -0.25, 0.0, 0.0], 5).unwrap();
    let want = [0.5, 0.25, 0.0, -0.25, -0.5];
    for (k, w) in want.into_iter().enumerate() {
        assert_eq!(
            p.eval(k as u64).unwrap(),
            Delay::from_samples_f64(w).unwrap()
        );
    }
    assert_full_block(&p, 1);
}

/// Quantization of a realistic delay polynomial stays within 2 output LSB of the
/// f64 monomial (whose own evaluation error is < 10⁻¹¹ samples here).
#[test]
fn realistic_delay_quantization_error() {
    let len = 61_440u64;
    // 150 km at 61.44 MS/s ≈ 30 740 samples; 250 m/s closing; 3 g; small jerk.
    let a = 30_740.123_456_789;
    let b = -8.339e-7;
    let c = 3.9e-14;
    let dd = -1.0e-21;
    let p = DelayPoly::from_f64_differences([a, b + c + dd, 2.0 * c + 6.0 * dd, 6.0 * dd], len)
        .unwrap();
    let lsb = 1.0 / 4_294_967_296.0;
    let mut worst: f64 = 0.0;
    for k in 0..len {
        let x = k as f64;
        let exact = a + x * (b + x * (c + x * dd));
        let got = p.eval(k).unwrap().to_samples_f64();
        worst = worst.max((got - exact).abs() / lsb);
    }
    eprintln!(
        "worst delay error: {worst:.3} LSB = {:.3e} samples",
        worst * lsb
    );
    assert!(worst < 2.0, "{worst} LSB");
}

#[test]
fn offset_phase_closed_form_near_hundred_days() {
    // 1.5 Hz at 100 days + 1/3 s: 0.5 cycles exactly, by exact rational time.
    let f = OffsetFrequency::from_hz_f64(1.5).unwrap();
    let third_s = (i128::from(PS_PER_SECOND) << 64) / 3; // floor of 1/3 s in Q64.64 ps
    let t = (i128::from(100 * DAY_PS) << 64) + third_s;
    let got = f.phase_at(AcqInstant::from_q64_64(t)).to_bits();
    // 1.5 Hz · (1/3 s − ε) is just under 0.5 cycle: ⌊2⁶³ − tiny⌋ = 2⁶³ − 1.
    assert_eq!(got, (1u64 << 63) - 1);
    let t_ps = AcqInstant::from_truth_time(TruthTime::from_ps(100 * DAY_PS));
    assert_eq!(f.phase_at(t_ps), Phase::ZERO);
}

/// Max absolute error of `sincos` over all 2²⁴ values of the top 24 phase bits, with
/// fixed pseudo-random low bits. Contract: ≤ 3·10⁻⁷.
// Reference oracle only: std f64 sin/cos (~1e-16 accurate) measures the f32 result and
// never feeds simulation output.
#[allow(clippy::disallowed_methods)]
#[test]
fn sincos_error_sweep() {
    let mut rng = SplitMix(0x51C0_5000);
    let mut max_sin: f64 = 0.0;
    let mut max_cos: f64 = 0.0;
    let two_pi_over_2_64 = 2.0 * std::f64::consts::PI / 18_446_744_073_709_551_616.0;
    for top in 0u64..(1 << 24) {
        let bits = (top << 40) | (rng.next() >> 24);
        let (s, c) = sincos(Phase::from_bits(bits));
        assert!(s.abs() <= 1.0 && c.abs() <= 1.0);
        let x = (bits as f64) * two_pi_over_2_64;
        let (rs, rc) = x.sin_cos();
        max_sin = max_sin.max((f64::from(s) - rs).abs());
        max_cos = max_cos.max((f64::from(c) - rc).abs());
    }
    eprintln!("sincos max abs error: sin {max_sin:.3e}, cos {max_cos:.3e}");
    assert!(max_sin <= 3e-7 && max_cos <= 3e-7);
}
