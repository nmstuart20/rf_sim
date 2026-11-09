//! Philox4x32-10 counter-based RNG, stream keying and the stream-id registry.
//!
//! Every random draw in the simulator comes from here (AGENTS.md, contracts §1.5).
//! A counter-based generator has no state: output block `n` is a pure function of
//! `(key, counter)`. Any thread, block partition or device can jump straight to any
//! sample and produce the same bits.
//!
//! **Philox4x32-10** is implemented exactly as in Salmon et al., "Parallel random
//! numbers: as easy as 1, 2, 3" (SC'11), and the Random123 reference code. A round maps
//! the counter `(c0, c1, c2, c3)` with key `(k0, k1)` to
//!
//! ```text
//! (hi0, lo0) = M0 · c0      (32 × 32 → 64-bit product, high and low halves)
//! (hi1, lo1) = M1 · c2
//! out = (hi1 ⊕ c1 ⊕ k0, lo1, hi0 ⊕ c3 ⊕ k1, lo0)
//! ```
//!
//! with `M0 = 0xD2511F53` and `M1 = 0xCD9E8D57`. Ten rounds run in sequence. The key
//! is bumped by the Weyl constants `W0 = 0x9E3779B9`, `W1 = 0xBB67AE85` (mod 2³²)
//! between rounds, nine bumps in all. The known-answer tests use the Random123
//! `kat_vectors` file.
//!
//! **Stream key** (contracts §1.5). The 64-bit Philox key comes from
//! (scenario seed `u64`, entity id `u32`, stream id `u32`):
//!
//! ```text
//! mix(z)   = SplitMix64 finalizer: z ^= z>>30; z *= 0xBF58476D1CE4E5B9;
//!            z ^= z>>27; z *= 0x94D049BB133111EB; z ^= z>>31     (all mod 2⁶⁴)
//! s        = mix(seed ⊕ KEY_DOMAIN)
//! k        = mix(s ⊕ (entity << 32 | stream))
//! key      = (k0, k1) = (low 32 bits of k, high 32 bits of k)
//! ```
//!
//! `mix` is a bijection on u64 (xor-shifts and odd multipliers are invertible). For a
//! fixed seed, `(entity, stream) ↦ k` is therefore injective: two different
//! `(entity, stream)` pairs in one scenario can never share a key.
//!
//! **Counter** = (sample index `u64`, sub-index `u64`) as
//! `(c0, c1, c2, c3) = (index_lo, index_hi, sub_lo, sub_hi)`. The sample index is
//! taken at the rate where the draw happens. The sub-index lets one model take several
//! blocks per sample. Each model documents which sub-indices it uses.
//!
//! **Stream ids** come only from [`StreamId`]'s registry constructors. A stream id is
//! `kind << 16 | instance`, where [`StreamKind`] codes are fixed forever. Adding a
//! kind appends a new code. Renumbering or reusing a code would change output bits
//! and is a numeric-contract change.

/// Philox4x32 round multiplier for counter word 0.
const M0: u32 = 0xD251_1F53;
/// Philox4x32 round multiplier for counter word 2.
const M1: u32 = 0xCD9E_8D57;
/// Weyl key-bump constant for key word 0 (golden ratio).
const W0: u32 = 0x9E37_79B9;
/// Weyl key-bump constant for key word 1 (√3 − 1).
const W1: u32 = 0xBB67_AE85;

/// Domain-separation constant XORed into the seed before mixing (ASCII `"iqsim-p1"`).
pub const KEY_DOMAIN: u64 = 0x6971_7369_6D2D_7031;

/// High and low 32-bit halves of the 64-bit product `a · b`.
#[inline]
const fn mul_hi_lo(a: u32, b: u32) -> (u32, u32) {
    let p = (a as u64) * (b as u64);
    ((p >> 32) as u32, p as u32)
}

/// One Philox4x32 round.
#[inline]
const fn round(c: [u32; 4], k: [u32; 2]) -> [u32; 4] {
    let (hi0, lo0) = mul_hi_lo(M0, c[0]);
    let (hi1, lo1) = mul_hi_lo(M1, c[2]);
    [hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0]
}

/// Philox4x32 with `rounds` rounds (`rounds ≥ 1`). Only the 10-round form is public;
/// the 7-round form exists so the unit tests can check the round structure against
/// the 7-round known-answer vectors too.
#[inline]
const fn philox4x32_rounds(ctr: [u32; 4], key: [u32; 2], rounds: u32) -> [u32; 4] {
    let mut c = round(ctr, key);
    let mut k = key;
    let mut r = 1;
    while r < rounds {
        k = [k[0].wrapping_add(W0), k[1].wrapping_add(W1)];
        c = round(c, k);
        r += 1;
    }
    c
}

/// The Philox4x32-10 block function: 128-bit counter and 64-bit key to 128 random bits.
///
/// Word order matches Random123 (`ctr[0]`, `key[0]` and `out[0]` are the first words
/// printed in its `kat_vectors` file).
pub const fn philox4x32_10(ctr: [u32; 4], key: [u32; 2]) -> [u32; 4] {
    philox4x32_rounds(ctr, key, 10)
}

/// SplitMix64 finalizer, a bijection on u64 (module docs).
#[inline]
pub const fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// A kind of random process. Codes are part of the numeric contract: append new
/// kinds, never renumber or reuse one. Code 0 is reserved and never used.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Hash)]
#[repr(u16)]
pub enum StreamKind {
    /// Thermal noise of one receiver channel (instance = channel index).
    ReceiverThermalNoise = 1,
}

impl StreamKind {
    /// The kind's fixed 16-bit code.
    pub const fn code(self) -> u16 {
        self as u16
    }

    /// The kind with code `code`, if it is registered.
    pub const fn from_code(code: u16) -> Option<Self> {
        match code {
            1 => Some(Self::ReceiverThermalNoise),
            _ => None,
        }
    }
}

/// A registered stream id: `kind << 16 | instance`.
///
/// There is no constructor from raw bits, so every stream id names a registered
/// [`StreamKind`]. Two models never share a stream as long as each model uses its own
/// kind.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, Hash)]
pub struct StreamId(u32);

impl StreamId {
    /// The stream `(kind, instance)`.
    pub const fn new(kind: StreamKind, instance: u16) -> Self {
        Self(((kind.code() as u32) << 16) | instance as u32)
    }

    /// Thermal-noise stream of receiver channel `channel`. Key it with the receiver's
    /// entity id (contracts §5).
    pub const fn receiver_thermal_noise(channel: u16) -> Self {
        Self::new(StreamKind::ReceiverThermalNoise, channel)
    }

    /// The 32-bit id that enters the key mix.
    pub const fn to_bits(self) -> u32 {
        self.0
    }

    /// The stream's kind.
    pub const fn kind(self) -> StreamKind {
        match StreamKind::from_code((self.0 >> 16) as u16) {
            Some(kind) => kind,
            // Unreachable: StreamId is only built from a registered kind.
            None => panic!("StreamId holds an unregistered kind"),
        }
    }

    /// The instance within the kind (for example the channel index).
    pub const fn instance(self) -> u16 {
        self.0 as u16
    }
}

/// The Philox key for (scenario seed, entity id, stream id), by the mix in the module
/// docs.
pub const fn stream_key(seed: u64, entity: u32, stream: StreamId) -> [u32; 2] {
    let s = mix64(seed ^ KEY_DOMAIN);
    let k = mix64(s ^ (((entity as u64) << 32) | stream.to_bits() as u64));
    [k as u32, (k >> 32) as u32]
}

/// One keyed Philox4x32-10 stream. It holds only the key and is cheap to copy.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PhiloxStream {
    key: [u32; 2],
}

impl PhiloxStream {
    /// The stream for (scenario seed, entity id, stream id).
    pub const fn new(seed: u64, entity: u32, stream: StreamId) -> Self {
        Self {
            key: stream_key(seed, entity, stream),
        }
    }

    /// The stream's Philox key.
    pub const fn key(self) -> [u32; 2] {
        self.key
    }

    /// The 128-bit block at counter (sample index, sub-index).
    pub const fn block(self, index: u64, sub_index: u64) -> [u32; 4] {
        let ctr = [
            index as u32,
            (index >> 32) as u32,
            sub_index as u32,
            (sub_index >> 32) as u32,
        ];
        philox4x32_10(ctr, self.key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The philox4x32 lines of the Random123 `kat_vectors` file, verbatim (provenance
    /// in the file header).
    const KAT: &str = include_str!("../tests/data/random123_kat_philox4x32.txt");

    fn hex(s: &str) -> u32 {
        u32::from_str_radix(s, 16).unwrap()
    }

    /// (rounds, ctr, key, expected) of one `philox4x32` line.
    type Kat = (u32, [u32; 4], [u32; 2], [u32; 4]);

    /// Every `philox4x32` line of the file.
    fn kat_vectors() -> Vec<Kat> {
        KAT.lines()
            .filter(|l| l.starts_with("philox4x32 "))
            .map(|l| {
                let f: Vec<&str> = l.split_whitespace().collect();
                assert_eq!(f.len(), 12, "malformed KAT line: {l}");
                (
                    f[1].parse().unwrap(),
                    [hex(f[2]), hex(f[3]), hex(f[4]), hex(f[5])],
                    [hex(f[6]), hex(f[7])],
                    [hex(f[8]), hex(f[9]), hex(f[10]), hex(f[11])],
                )
            })
            .collect()
    }

    #[test]
    fn random123_known_answers() {
        let v = kat_vectors();
        assert_eq!(
            v.len(),
            6,
            "expected three 7-round and three 10-round vectors"
        );
        let mut ten = 0;
        for (rounds, ctr, key, want) in v {
            assert_eq!(philox4x32_rounds(ctr, key, rounds), want, "R={rounds}");
            if rounds == 10 {
                assert_eq!(philox4x32_10(ctr, key), want);
                ten += 1;
            }
        }
        assert_eq!(ten, 3);
    }

    #[test]
    fn stream_block_counter_layout() {
        let s = PhiloxStream::new(1, 2, StreamId::receiver_thermal_noise(3));
        let idx = 0x0123_4567_89AB_CDEF;
        let sub = 0xFEDC_BA98_7654_3210;
        let want = philox4x32_10(
            [0x89AB_CDEF, 0x0123_4567, 0x7654_3210, 0xFEDC_BA98],
            s.key(),
        );
        assert_eq!(s.block(idx, sub), want);
    }

    #[test]
    fn key_mix_matches_definition() {
        let id = StreamId::receiver_thermal_noise(5);
        assert_eq!(id.to_bits(), 0x0001_0005);
        let s = mix64(42 ^ KEY_DOMAIN);
        let k = mix64(s ^ (7u64 << 32 | 0x0001_0005));
        assert_eq!(stream_key(42, 7, id), [k as u32, (k >> 32) as u32]);
    }

    #[test]
    fn mix64_known_value() {
        // SplitMix64 with state 0: the first output is mix(0x9E3779B97F4A7C15), whose
        // published value is 0xE220A8397B1DCDAF.
        assert_eq!(mix64(0x9E37_79B9_7F4A_7C15), 0xE220_A839_7B1D_CDAF);
    }

    #[test]
    fn registry_round_trip() {
        let id = StreamId::receiver_thermal_noise(0xBEEF);
        assert_eq!(id.kind(), StreamKind::ReceiverThermalNoise);
        assert_eq!(id.instance(), 0xBEEF);
        assert_eq!(StreamKind::from_code(0), None);
        assert_eq!(
            StreamKind::from_code(StreamKind::ReceiverThermalNoise.code()),
            Some(StreamKind::ReceiverThermalNoise)
        );
    }
}
