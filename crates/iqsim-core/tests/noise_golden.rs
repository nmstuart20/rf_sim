//! Golden hashes for Philox noise and ln (prompt 3).
//!
//! A change to any output bit fails the test. Bless only on purpose
//! (`cargo xtask bless`) and explain it in docs/progress.md. The CUDA backend
//! (prompt 20) must reproduce these hashes.

use iqsim_core::{NoiseStream, SampleIndex, StreamId, ln, noise};
use iqsim_testkit::assert_golden;

const N: usize = 1_000_000;

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

fn push(bytes: &mut Vec<u8>, (i, q): (f32, f32)) {
    bytes.extend(i.to_bits().to_le_bytes());
    bytes.extend(q.to_bits().to_le_bytes());
}

/// 10⁶ noise samples: 500 000 consecutive samples of one stream, starting just below
/// 2³² so the counter's high word changes, then 500 000 samples at arbitrary
/// (seed, entity, channel, index).
#[test]
fn noise_golden() {
    let mut bytes = Vec::with_capacity(N * 8);
    let s = NoiseStream::new(0x5EED, 3, StreamId::receiver_thermal_noise(1));
    let start = (1u64 << 32) - 250_000;
    for n in start..start + (N as u64) / 2 {
        push(&mut bytes, s.sample(SampleIndex::new(n)));
    }
    let mut rng = SplitMix(0x0003_401E);
    for _ in 0..N / 2 {
        let seed = rng.next();
        let ids = rng.next();
        let n = rng.next();
        let id = StreamId::receiver_thermal_noise((ids >> 32) as u16);
        push(&mut bytes, noise(seed, ids as u32, id, SampleIndex::new(n)));
    }
    assert_golden("core/noise_1e6", &bytes);
}

/// ln of 10⁶ inputs: half uniform over the positive finite bit patterns (subnormals
/// included), half from the Box-Muller range [2⁻⁶⁵, 1).
#[test]
fn ln_golden() {
    let mut rng = SplitMix(0x0003_0106);
    let mut bytes = Vec::with_capacity(N * 4);
    for k in 0..N {
        let r = rng.next();
        let bits = if k % 2 == 0 {
            1 + (r % 0x7F7F_FFFF) as u32
        } else {
            // Exponent field 62..=126 (2⁻⁶⁵ ≤ x < 1), arbitrary mantissa.
            let e = 62 + (r >> 32) as u32 % 65;
            (e << 23) | (r as u32 & 0x007F_FFFF)
        };
        bytes.extend(ln(f32::from_bits(bits)).unwrap().to_bits().to_le_bytes());
    }
    assert_golden("core/ln_1e6", &bytes);
}
