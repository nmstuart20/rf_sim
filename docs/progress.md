# Progress log

Append one entry per prompt. Newest at the bottom. Each entry covers what was built, deviations from the design and why, new dependencies, and open questions.

---

## Prompt 0: Repository bootstrap (2026-11-08)

**Built**

- Cargo workspace `iqsim` (edition 2024, resolver 3). Toolchain pinned to Rust 1.99.1 in `rust-toolchain.toml`. `overflow-checks = true` in dev and release.
- 19 crates under `crates/`, each with a crate-level doc giving its responsibility and allowed dependencies. They hold no simulation logic yet. `iqsim-cli` builds the `sim` binary: `--version` and `--help` work, and every subcommand exits 2 as not yet implemented.
- Internal dependency DAG declared in the manifests now, so it is fixed from day one. See `docs/contracts.md` §9.
- `#![forbid(unsafe_code)]` in every crate except `iqsim-render-cuda` and `iqsim-hil`, which use `#![deny(unsafe_code)]` until their FFI modules exist.
- Invariants enforced by tooling:
  - `clippy.toml` denies the std f32/f64 transcendentals and `mul_add` (use libm, no FMA), and `HashMap`/`HashSet` (deterministic iteration).
  - Workspace lints warn on `missing_docs`. CI turns warnings into errors.
- xtask (std only):
  - `ci`: fmt check, clippy `-D warnings`, tests, rustdoc `-D warnings`, deps check.
  - `bless`: runs tests with `IQSIM_BLESS=1`; extra arguments go to cargo test.
  - `gpu-test`: placeholder.
  - `deps`: forbids control → runtime/hil and runtime → control/hil, and allows only cli to depend on cuda/hil.
- `iqsim-testkit::assert_golden(name, bytes)`: compares SHA-256 against `tests/golden/<name>.sha256`. Mismatch, missing and malformed files fail with a clear message and the bless hint. `IQSIM_BLESS=1` rewrites the file. Includes a self-check golden `tests/golden/testkit/self_check.sha256`.
- `AGENTS.md`: invariants, dependency policy, commands, Definition of Done. `CLAUDE.md` only imports `@AGENTS.md`, so Claude Code and other agents read the same file.
- `docs/contracts.md`, contract version 0 (proposed): numeric rules, amplitude normalization and overflow policy, ClockMap/timestamp/LO models, the commit boundary, tap responses and state ownership, run manifest, transport and delivery guarantees, resource budgeting, and the list of decisions phase 1 must settle.

**Deviations**

- `AGENTS.md` replaces `CLAUDE.md` because several agents (Claude Code, Codex) work on this repo. `agent_prompts.md` references were updated to match.
- Added the crate `iqsim-contracts`, which the original crate list does not have. It holds the versioned protocol, handler-trait, IQ-framing and manifest types, with no I/O. Control adapters and the runtime both depend on it, so no runtime ↔ control cycle can form. This satisfies the prompt's dependency-graph requirement.
- Added the crate `iqsim-testkit` to hold the golden helper, used only as a dev-dependency.
- Added `cargo xtask deps` and included it in `ci`.
- `NUMERIC_CONTRACT_VERSION` is planned for `iqsim-core` (contracts.md §9), not `iqsim-contracts`.

**Dependencies**

- `sha2` 0.10.9 and `hex` 0.4.3, MIT/Apache-2.0, used only by iqsim-testkit. They are the de-facto standard RustCrypto SHA-256 and hex encoding, with 9 small transitive packages and no build scripts that use the network. SHA-256 is also the manifest hash, so sha2 will be reused there.

**Open questions**

- `libm` is not added yet. Prompt 1 or 2 adds it when the first control-plane math lands.
- The pinned toolchain must be present in the air-gapped RHEL 9 image (rustup offline install, or a vendored toolchain tarball). This belongs in phase 6 packaging.
- Every *proposed* value in contracts.md §10 still needs settling in phase 1.

---

## Prompt 1: Time and rate types (2026-11-08)

**Built** (all in `iqsim-core`)

- `TruthTime` and `TruthDuration` (`time.rs`): i64 picoseconds with checked add, subtract, multiply and negate. Overflow returns `TimeError::TimeOverflow`.
- `SampleRate` (`rate.rs`): reduced rational `num/den` S/s. `from_hz` and `from_ratio` reject zero and anything outside the settled bounds (den ≤ 2³², num ≤ 2⁴⁸, ≤ 1 GS/s, checked after reduction).
- `SampleIndex` and `ChannelTimebase` (`timebase.rs`): exact index → time (`ExactTime { floor, rem, rem_den }`), floor time → index (`index_at_or_before`) and a ceiling companion (`index_at_or_after`). The u128 overflow argument is in the module docs and in contracts §1.2.
- `AcqInstant` (`clock.rs`): `ps + frac·2⁻⁶⁴` ps, a signed Q64.64 value, with exact `to_q64_64`/`from_q64_64`.
- `ClockMap` trait: `instant`, `index_at_or_after`, `index_at_or_before`, `id`. `IdealClockMap` floors the exact map onto the 2⁻⁶⁴ ps grid. `ClockMapId` is SHA-256 of a versioned canonical encoding.
- A separate `TimestampClock` trait with `IdealTimestampClock`. It returns `ReportedTime`, a different type from `AcqInstant`, so a reported timestamp cannot drive geometry or phase.
- `Block` (half-open `[start, end)`), `BlockSpec` (grid aligned to multiples of the block length from index 0) and the `Blocks` iterator (`tile` clips to a span, and `iter_from` runs to the end of the u64 index space). It yields `(block number, piece)`.
- `TimeError`: one error enum for the crate's time, rate, index and block failures.
- Tests: 18 unit tests and 11 integration tests in `crates/iqsim-core/tests/timebase_props.rs`.
  - Properties (2048 cases each): exact index → time reconstruction; round trips through both inverses of the timebase and the clock map; correct bracketing by the inverse for arbitrary sub-ps instants; strict monotonicity of instants and non-decreasing time → index; overflow reported exactly when the true time does not fit in i64 (arbitrary u64 index, epoch and accepted rate); block tiling with no gaps, overlaps or empty pieces.
  - Edge cases: 61.44 MS/s (16276 + 1/24 ps period, exact every 24th sample, 1 ms = 61 440 samples, 100 days, last representable sample) and 250 MS/s (exact 4000 ps steps, negative epoch, neighbour selection, 1 ms blocks). Also the corner rates of the accepted region at 100 days.

- An independent review (Codex CLI) found no correctness bugs in the conversions. It raised four minor issues, all fixed:
  - `ExactTime` fields are now private, so callers cannot build one that panics in `to_acq_instant`.
  - `SampleRate` now orders by numeric value instead of the derived field order, which put 1/2 S/s above 1 S/s.
  - Two test oracles could panic at the slowest accepted rate (1/2³² S/s), where sample 1 is already beyond the i64 time range. They now handle that, and `slowest_rate_overflow_cases` covers it.

- Proptest policy, decided with the user: seeds stay random so runs keep finding new cases, and every `proptest-regressions/` file proptest writes is committed so a failure is never lost (AGENTS.md, Commands). The case count is 2048 by default and `PROPTEST_CASES` overrides it.

**Contracts changed** (`docs/contracts.md`, still version 0 because §10 is not fully settled)

- §10 decision 1 settled: rate bounds den ≤ 2³², num ≤ 2⁴⁸, ≤ 1 GS/s. The no-overflow proof is in §1.2.
- §10 decision 2 settled: `AcqInstant` uses a fixed 2⁻⁶⁴ ps fraction (floor). `index_to_time` also returns the exact rational remainder. The inverse map returns only an index; callers get the sub-sample offset by exact Q64.64 subtraction (§3.1).
- Tolerances: carrier phase ≤ 0.01 rad and delay ≤ 0.001 samples are settled. Time representation uses < 10⁻²⁰ rad and < 2⁻⁶⁴ ps of them, so prompt 2 owns the quantization share. Added an acquisition-instant representation error row (< 2⁻⁶⁴ ps).
- §1.1 and §3.1 now describe the representations and the `IdealClockMap` identity encoding.

**Deviations**

- The prompt asks for floor time → index. `ChannelTimebase` and `ClockMap` also provide the ceiling (`index_at_or_after`), which the contracts §3.1 trait sketch already names. A time before the epoch is an error for the floor and maps to 0 for the ceiling.
- `BlockSpec::block_len()` instead of `len()`, because clippy's `len_without_is_empty` lint applies to the name `len`.
- Block length has no upper bound yet. Prompt 2 pins it with the block polynomials (C(k,3) must fit in u64).
- `IdealTimestampClock` is the identity. Offset and drift models arrive in prompt 15.

**Dependencies**

- `sha2` (already a workspace dependency, used by testkit) is now a normal dependency of `iqsim-core`, for `ClockMapId`. Contracts §3.1 asks for a content hash, and SHA-256 is the project-wide hash. It adds no new packages to the lockfile.
- `proptest` 1.12.0, MIT/Apache-2.0, dev-dependency of `iqsim-core` only, declared with `default-features = false, features = ["std"]` to drop fork and timeout support. It is the standard property-testing crate, and writing shrinking and strategies in-house is not worth it. It adds 12 packages: bitflags, chacha20, core_detect, cpufeatures 0.3, getrandom 0.4, num-traits, r-efi, rand 0.10, rand_core 0.10, rand_xorshift, regex-syntax and unarray. None download anything at build time. It is test-only and never reaches a simulator binary. Proptest's internal RNG drives test-case generation only; simulation randomness still comes only from Philox.

**Open questions**

- `SampleRate` rejects rates above 1 GS/s. A wider-band future emitter source would need the bounds and proof revisited.

---

## Prompt 2: Fixed-point delay, carrier phase, polynomials and sin/cos (2026-11-08)

**Built** (all in `iqsim-core`, written so the CUDA backend can port it line by line)

- `wide` (public module): `U128 {hi, lo}` with wrapping add, sub, neg, `× u64` and `u64 × u64 → U128`. It is built only from `mul_hi_u64` (four 32 × 32 products), wrapping 64-bit multiply, `add_with_carry` and `sub_with_borrow`. It never uses native u128 arithmetic. Native conversions exist for host construction and tests.
- `fixed`:
  - `Delay` (Q32.32 samples in i64, range-checked, no arithmetic operators) and `Phase` (u64 cycle fraction, wrapping add/sub).
  - Ties-to-even float → integer helpers: `f64_to_i64_ties_even`, `f32_to_i64_ties_even`, `require_finite_*`.
  - Exact f64 → 128-bit fixed conversions: checked for delay, mod one cycle via the exact `x − trunc(x)` for phase.
- `poly`:
  - `BlockPoly<T>` with `DelayPoly`/`PhasePoly` aliases. Newton-form D0..D3 128-bit registers (Q32.96 for delay, 0.128 cycle for phase).
  - `jump(k)`, `eval(k)`, and `cursor(k)` → `PolyCursor::step()` (three wrapping adds).
  - Exact u64 `binomials(k)`. `MAX_BLOCK_LEN = 4 801 281`.
- `trig`:
  - `sincos(Phase) -> (f32, f32)`: 10-bit table index, 24-bit remainder, cubic sin δ / quadratic cos δ correction, every f32 op explicit.
  - `SIN_TABLE` (1024 × f32) is built at compile time by an integer-only Taylor routine in Q2.62 with explicit ties-to-even f32 rounding and quarter-wave symmetry. No float library is used.
- `offset`: `OffsetFrequency` (i64, 2⁻³² Hz) and `phase_at(AcqInstant)`. It computes the exact ⌊F·T/(2⁴⁴·5¹²)⌋ mod 2⁶⁴ with a 192-bit product, a shift, and division by 5¹² over 32-bit digits, all in u64 arithmetic.
- `NumericError` (non-finite, out of range, block length, index outside block).
- Tests: 23 new unit tests, 15 in `tests/numeric_props.rs` and 4 goldens in `tests/numeric_golden.rs`.
  - Primitives and `U128` against native u128/i128 (proptest).
  - Jump against a native reference at random k, and step against jump from random k0 (arbitrary wrapping registers).
  - Every k of full blocks: 61 440 and 250 000 samples, delay and phase, plus one block of `MAX_BLOCK_LEN` samples.
  - A phase polynomial that wraps more than 20 000 times per block.
  - A negative delay rate: strictly decreasing, crossing zero.
  - Realistic delay-quantization error, measured 1.0 LSB worst case.
  - Ties-to-even and range boundaries, and non-finite rejection on every constructor.
  - Offset phase:
    - against an independent bit-serial bignum reference
    - additivity in t (|t| up to 104 days) and in Δf
    - φ(100 days + δ) = φ(δ) exactly for whole-Hz offsets
    - closed-form values, the floor of tiny negative phases, and extreme inputs
  - sincos:
    - table correctly rounded (against f64), symmetric, and no entry near a rounding tie
    - exact quadrant points
    - a 2²⁴-point sweep of the top 24 bits with random low bits. **Max abs error: sin 5.98·10⁻⁸, cos 5.97·10⁻⁸** (target 3·10⁻⁷).
  - Goldens, each 10⁶ evaluations except the table: `core/sincos_table`, `core/sincos_1e6`, `core/block_poly_1e6` (jump and step, delay and phase) and `core/offset_phase_1e6`. All four were blessed for the first time; no existing golden changed.

**Contracts changed** (`docs/contracts.md`, still version 0: §10 decisions 8–10 remain)

- §10 decisions 3–7 settled: maximum block length and quantization bound; sincos design; subnormals preserved; offset-frequency representation; normalization (R₀ = 50 Ω RMS, S = 40), headroom limit L = 2⁶² (one guard bit), null floor and V_ref definition.
- Amplitude tolerance, decided with the user. Above the −40 dB floor the bound is 0.01 dB. Below it, the absolute bound is **10⁻⁵·V_ref**, not the proposed 10⁻³·V_ref, which would have been 100× looser just below the floor than just above it. The bound is now continuous.
- §1.1 now names where intentional wrapping lives and where it is forbidden. §1.4 pins the float → integer and integer → float rules, the subnormal policy, and CPU/CUDA compiler flags (CUDA `--fmad=false -ftz=false -prec-div=true -prec-sqrt=true`, no fast intrinsics). §1.5 and §1.6 describe the implemented sincos and offset-phase algorithms.

**Deviations**

- The headroom guard is L = 2⁶² (the contract text proposed "2⁶⁰ of slack" while also saying "about 3 guard bits", which disagree). One full guard bit is simple to state and still allows 2²² full-scale paths.
- `MAX_BLOCK_LEN` is enforced by `BlockPoly`, not by `BlockSpec`. `BlockSpec` keeps tiling arbitrary lengths (prompt 1 tests use 2⁴⁰ and 2⁶²). The render configuration must check `block_len ≤ MAX_BLOCK_LEN` when it builds work lists.
- `OffsetFrequency::phase_at` returns the exact floor. Nothing in prompt 2 consumes it yet; the LO model (prompt 15) and receiver chain will.
- Fixed a prompt 1 test-infrastructure bug. proptest's default failure persistence finds no `lib.rs` beside `tests/*.rs`, so it printed "failed to find lib.rs" and saved nothing. That silently broke the AGENTS.md rule to commit regression files. Both property files now use `FileFailurePersistence::Direct("proptest-regressions/<file>.txt")`, which I verified writes the file. AGENTS.md notes this for future integration tests.

**Dependencies**

- `iqsim-testkit` is now a dev-dependency of `iqsim-core`, for the goldens. It is an internal crate with no iqsim dependencies, so no cycle forms. There are no new external crates. `libm` is still not needed: the table builder is integer-only and nothing in this prompt computes a control-plane transcendental.

**Open questions**

- The delay polynomial's range check over a block (true delay inside Q32.32) belongs to the geometry fitter (phase 3). Evaluation cannot detect wrap-around.
- The CUDA backend must embed `SIN_TABLE` bytes exactly. The `core/sincos_table` golden pins them.

**Independent review (Codex CLI)**

- The first run, at `xhigh` reasoning, stalled after about 1.5 hours with no output for 80 minutes, and I stopped it. Before stalling it had confirmed the 128-bit primitives, the Newton jump/step identity and the offset-phase division, with all 71 core tests passing. It also raised two issues, both fixed:
  - The quantization bound was stated as "within 1.5 LSB". The exact worst case at `MAX_BLOCK_LEN` is 1.5 + 2⁻²¹ LSB, so the bound now reads < 2 LSB: 4.7·10⁻¹⁰ samples of delay and 6.9·10⁻¹⁹ rad of phase. Both are still negligible against the budgets.
  - A true delay just inside −2³¹ samples could be pushed past the range edge by quantization and wrap to nearly `i64::MAX`. I added `DELAY_RANGE_MARGIN_LSB = 2`. Contracts §1.3 now requires the geometry fitter to keep the true delay that far inside the Q32.32 range.
- A second, narrower run at `high` reasoning covered `fixed.rs` and `trig.rs`. It found no implementation bugs (ties, range checks, wrapping, signed zero, table construction, CUDA/FMA portability). It found one wrong comment: a power-of-two-scaled product *can* be subnormal (2⁻¹⁰⁷⁴·2³²). The scaling is still exact, so only the wording changed.

---

## Prompt 3: Philox and Gaussian noise (2026-11-08)

**Built** (all in `iqsim-core`, written so the CUDA backend can port it line by line)

- `philox`:
  - `philox4x32_10(ctr, key)`: Philox4x32-10 exactly as in Random123 (Salmon et al., SC'11). Passes all six philox4x32 known-answer vectors (three 7-round, three 10-round) from Random123's `tests/kat_vectors`. The vectors are checked in verbatim at `crates/iqsim-core/tests/data/random123_kat_philox4x32.txt`, with the upstream commit (`9545ff6`), the SHA-256 of the full file and the BSD-3 notice. A private round-count parameter exists only so the 7-round vectors can test the round structure.
  - `stream_key(seed, entity, stream)`: `k = mix(mix(seed ⊕ KEY_DOMAIN) ⊕ (entity << 32 | stream))`, where `mix64` is the SplitMix64 finalizer. Because `mix64` is a bijection, the key is injective in (entity, stream) for a fixed seed.
  - `PhiloxStream::block(index, sub_index)`: counter `(index_lo, index_hi, sub_lo, sub_hi)`.
  - Stream-id registry: `StreamKind` (u16 codes, append-only, 0 reserved) and `StreamId = kind << 16 | instance`, which can only be built from a registered kind. The only kind so far is `ReceiverThermalNoise` (code 1, instance = channel).
- `log`: `ln(f32) -> Result<f32, NumericError>`, the fdlibm decomposition with Taylor coefficients and a 12-bit/f32 split of ln 2. It uses integer range reduction and handles subnormals. There are no std or libm calls, and every operation is explicit. **Measured error ≤ 0.86 ulp**. On the Box-Muller range [2⁻⁶⁵, 1] the absolute error is ≤ 1.932·10⁻⁶ (asserted ≤ 2·10⁻⁶).
- `noise`: `noise(seed, entity, stream, SampleIndex) -> (f32, f32)` and `NoiseStream::sample` (the key is mixed once per stream). Box-Muller over one Philox block per sample, with **unit complex power** (`E[I²] = E[Q²] = ½`, `E|z|² = 1`):
  - The radius uniform is built exactly from 64 bits (truncated to 24 significant bits, x = 0 → 2⁻⁶⁵).
  - `r = sqrt(−ln u)`, so `r ≤ 6.71`.
  - The angle is word 2 as a u32 phase into `sincos`. Word 3 is unused.
- Tests: 11 new unit tests, 8 property tests (`tests/noise_props.rs`), 2 statistics tests (`tests/noise_stats.rs`) and 2 goldens (`tests/noise_golden.rs`).
  - Known-answer vectors; the counter layout; the key-mix definition; the published SplitMix64 first output; registry round trip.
  - ln: the split constants (e·LN2_HI exact for every exponent), the rejected inputs, every power of two, an exhaustive sweep of 13 binades plus subnormals, ≤ 1 ulp over random positive finite f32, and monotonicity (exhaustive over [0.5, 2), random elsewhere).
  - Radius uniform: equals the exact truncated value, is monotone, and stays in (0, 1). Noise is finite and within the radius bound at random inputs. The free function matches the stream form and the raw block function.
  - Key mix: an explicit inverse of `mix64` round-trips, and unmixing a key recovers (entity, stream).
  - Statistics, with thresholds fixed in advance (6 standard errors; Anderson-Darling 1 % value 3.857). Over 10⁷ samples:
    - E|z|² = 0.999931 and E[IQ] = 1.4·10⁻⁴
    - I: variance 0.499866, skewness 1.96·10⁻³, excess kurtosis −1.3·10⁻⁴
    - Q: variance 0.500064, skewness −6.8·10⁻⁴, excess kurtosis 1.05·10⁻³
    - A²: I 1.153, Q 0.282, |z|² vs Exp(1) 0.390
  - Cross-correlation over 10⁶ samples (complex, lags −2..2, plus power correlation) for adjacent channels, entities and seeds, and for the same stream at lags 1..5. The largest |ρ| was 1.7·10⁻³ against a bound of 4.2·10⁻³.
  - Goldens `core/noise_1e6` (500 000 consecutive samples across the 2³² counter-word boundary plus 500 000 random tuples) and `core/ln_1e6`. Both are blessed for the first time; no existing golden changed.

**Contracts changed** (`docs/contracts.md`, still version 0: §10 decisions 9 and 10 remain)

- §10 decision 8 settled. §1.5 now pins the ln algorithm and error, the Philox constants and known answers, the key mix, the counter layout, the stream-id registry table and the noise transform and variance. §5's noise bullet names the registry constructor.

**Deviations**

- The prompt says counter = (sample index, sub-index). The sub-index is a full u64 in counter words 2–3, so no counter word is left reserved.
- The output is unit *complex* power (each component variance ½), not unit variance per component. A noise power P then scales by √(P·R₀) directly, matching §2's RMS convention. Box-Muller therefore uses `sqrt(−ln u)` instead of `sqrt(−2 ln u)`.
- The radius uniform takes 64 bits (two words), not 32, so the tail reaches 9.5 σ per component instead of 6.7 σ. That costs one Philox block per complex sample, with one word unused. If noise generation turns out to limit CPU real-time margin (plan, phase 5), a two-samples-per-block layout is the fallback, and it would be a contract change.
- `ln` returns `Result` at the public boundary (zero, negative and non-finite inputs are errors, §1.4). Noise uses an unchecked internal form whose input is positive by construction.
- `noise_stats` takes about 23 s in the debug test profile (10⁷ samples plus three sorts and Anderson-Darling). It is not `#[ignore]`d, because the prompt makes it part of the required tests.

**Dependencies**

- `libm` 0.2.16, MIT, **dev-dependency of `iqsim-core` only**, declared in `[workspace.dependencies]`. It provides `erfc`, `expm1` and `log` in f64 for the Anderson-Darling reference CDFs. The normal tails need an accurate erfc, which std lacks, and writing one in-house for a test oracle is not worth it. It has no dependencies. Its `build.rs` only emits cfg flags and never touches the network. It is already the project's planned control-plane math crate (AGENTS.md), so later prompts can promote it to a normal dependency without new vetting. Simulation code does not use it here.

**Open questions**

- The registry has one kind. Later prompts append theirs (SSR PRF jitter in prompt 12, the clock/GPSDO and phase-noise processes in prompt 15) and document their entity and sub-index use in the §1.5 table.
- Entity ids are raw u32 for now. The scenario resolver (phase 3+) must assign them deterministically, because the key depends on them.

**Independent review (Codex CLI, `high` reasoning)**

- It confirmed that Philox matches Random123's `philox.h` (round function, word order, key bumps) and found no correctness bugs in key mixing, ln, the radius uniform, Box-Muller or CUDA portability. It raised two minor issues, both fixed:
  - The ln absolute-error figure was rounded down: the worst case is 1.9313·10⁻⁶ at `0x1FB2E007`. The docs now say ≤ 1.932·10⁻⁶ (asserted ≤ 2·10⁻⁶).
  - The lagged cross-correlation normalized by n while summing n − 2 products, which attenuated ρ slightly. It now divides by the true product count.

---

## Prompt 4: Units, scenario schema skeleton and `sim lint` (2026-11-08)

**Built**

- `iqsim-units`: quantities parsed from strings with case-sensitive units, each with `FromStr`, serde `Deserialize` and a schemars `JsonSchema` (string with pattern, description and example). The full vocabulary is in contracts §11.
  - An exact decimal parser (`decimal.rs`, u128 mantissa and base-10 exponent) feeds every type, so there is no float round trip.
  - Exact types: `Frequency` (i64 µHz), `SampleRate` (wraps the core rational: `"61.44 MS/s"` is exactly 61 440 000/1), `Duration` (wraps `TruthDuration`, i64 ps) and `DataRate` (u64 bit/s). A value finer than the resolution is an error, not a rounding.
  - Float types: `Power` (W; dBm/dBW via `libm::exp10`), `Gain`, `Distance`, `Speed` and `Angle`. Decimal unit factors are applied exactly before one correctly rounded conversion.
  - `Angle` takes deg, °, rad and DMS in both `12°34'56.7"` and `12d34m56.7s` forms. `Latitude` and `Longitude` take numbers, decimal or DMS strings, and N/S/E/W letters.
  - `Altitude` needs an explicit reference (`MSL`, `AGL`, `ellipsoid`/`HAE`), and `FLnnn` is a standard pressure altitude.
  - `Position` is `{ lat, lon, alt }` or `{ mgrs, alt }`.
  - A bare number where a unit is required gets a message such as `frequency needs a unit: got the bare number 1090; write it as a string such as "1090 MHz" (units: Hz, kHz, MHz, GHz)`.
- **MGRS and UTM, in house** (`mgrs.rs`, `utm.rs`):
  - AA lettering, band-based selection of the 2000 km northing cycle, the Svalbard non-zones and 0–5 digit precision.
  - The UTM inverse is Karney's Krüger series to n⁶ with Newton iteration, through `libm`.
  - A position is the south-west corner of the square, and `center()` is also exposed.
- **Reference tests** (`tests/data/`, regenerated by `gen_oracle.py`):
  - 60 MGRS references decoded by GeographicLib's published GeoConvert 2.5.2 service. They are hand-picked edge cases (32V, 31X–37X, bands C and X, both sides of the equator, the antimeridian, 0–5 digit precision) plus random points, and their strings were encoded by the NGA GEOTRANS-based Python `mgrs` package.
  - 400 random UTM points inverted by PROJ 9.8.1.
  - Worst differences: **5·10⁻¹² °** against GeoConvert and **5·10⁻¹³ °** against PROJ, both at the files' print precision (under 1 µm). The asserted tolerance is 10⁻¹⁰ ° (11 µm). The UTM centre of each square also matches GeoConvert's easting and northing exactly.
- `iqsim-scenario`: strict serde types for the top level of a scenario.
  - `schema_version`, `metadata`, `epoch` (UTC, POSIX seconds + ns, in-house civil-date conversion) and `seed`.
  - `run_mode`, `library.paths` and `platforms` (position, `trajectory = { kind = "static" }`, mounts).
  - `emitters` (frequency, power, behavior and waveform library references).
  - `receivers` with channels: centre frequency, sample rate, bandwidth, and separate `clock` and `lo` placeholders (`{ model = "ideal" }`).
  - `antenna_ports` (B_out bandwidth), `sinks` (kind, resources) and `resources` (NIC/disk/USB/PCIe, optional `capacity`).
  - Every table has `deny_unknown_fields`, and names carry `toml::Spanned` source spans.
- Validation (`validate.rs`):
  - one entity namespace, with mount and channel names scoped to their platform and receiver
  - typed references, including mounts on the referenced platform
  - positive frequencies, and bandwidth ≤ complex sample rate (compared exactly)
  - at least one channel per receiver
  - diagnostics in source order, with a `note:` at a duplicate's first definition
- `Diagnostic::render` prints `file:line:col: error: message`, the source line and a caret under the span. Columns count characters.
- JSON Schema: `cargo xtask schema` runs `iqsim-scenario`'s `write_schema` example into `schema/scenario.schema.json`. xtask stays std-only. A test fails while the checked-in file is stale. The null alternatives schemars emits for `Option` are stripped (TOML has no null), and `schema_version` is pinned with `const`.
- `sim lint <FILE>`: exit 0 and `FILE: ok`, exit 1 with diagnostics and an error count on stderr (or for an unreadable file), exit 2 for usage errors.
- `examples/minimal.toml`: one static CW emitter, one receiver/channel on an MGRS-placed site, a SigMF sink and a disk resource. It lints clean (CLI test and scenario test).
- Negative fixtures in `crates/iqsim-cli/tests/lint/`: 18 files, each `minimal.toml` with one mistake and an `# expect: LINE:COL: fragment` header. The line and column were computed by the generator script, independently of the Rust code.
  - The mistakes: missing unit, bare number, wrong unit kind, unit case, unknown field, bad MGRS column, polar MGRS, altitude without reference, non-UTC epoch, dangling platform, dangling sink, wrong-kind reference, missing mount, duplicate name, bandwidth > rate, bad run mode, schema version, DMS minutes ≥ 60.
  - The test checks the first diagnostic's location and text, and the caret column.
- Cross-check, not part of CI: Python `jsonschema` (Draft 2020-12) accepts the example against the checked-in schema and rejects the 8 syntactic fixtures. The other 10 are semantic errors that only `sim lint` can catch.
- Tests: 17 unit tests and 14 integration tests in units (8 properties at 2048 cases), 6 unit tests and 8 integration tests in scenario, and 6 CLI tests.

**Contracts changed**

- New §11, "Scenario quantities and files": quantity text rules, representations, altitude references, DMS, MGRS conventions and accuracy, scenario file rules, lint checks and diagnostics format. No numeric-contract values changed, and the version stays 0.

**Deviations**

- `AltitudeRef` has a fourth reference, `PressureStd`, for flight levels. The prompt lists FL as a unit beside ellipsoid/MSL/AGL references, but FL350 is a pressure altitude, not a height above any of those, so it cannot take one. The Mode C chain (plan, "Mode C altitude") will convert it.
- Added quantity types the prompt does not list: `SampleRate` (receiver channels need one) and `DataRate` (resource capacities, contracts §8).
- MGRS positions use the square's south-west corner. GeographicLib uses the centre, so the tests compare `center()` against it.
- Polar UPS references are rejected with a clear message, not implemented.
- The scenario has a required `schema_version = 1` (plan: every interface is versioned).
- Clock and LO placeholders are per channel, following the prompt's wording. A receiver whose channels share one oscillator can be modeled when the clock models arrive (prompt 15).
- toml only reports the first syntax or type error. Semantic errors are all reported.

**Dependencies** (all declared in `[workspace.dependencies]`)

- `serde` 1.0.229 (derive), MIT/Apache-2.0: the standard serialization framework, and toml and schemars require it. Used by units and scenario.
- `toml` 1.2.1, MIT/Apache-2.0, with `default-features = false, features = ["std", "serde", "parse"]`: TOML 1.1 parser with byte spans for errors and `Spanned` values, which file:line:column diagnostics need. Writing a spec-compliant TOML parser in house is not worth it. It is used only by iqsim-scenario, plus iqsim-units as a dev-dependency.
- `schemars` 1.2.3 (derive), MIT: JSON Schema from the Rust types, which the prompt requires. Used by units and scenario.
- `serde_json` 1.0.154, MIT/Apache-2.0: already a schemars dependency. Used to write the schema and in schema transforms.
- `libm` (already vetted in prompt 3) is now a **normal** dependency of iqsim-units, for dB conversion and the UTM inverse.
- New packages built: serde, serde_core, serde_derive, toml, toml_parser, toml_datetime, serde_spanned, winnow, schemars, schemars_derive, serde_derive_internals, serde_json, zmij, itoa, memchr, dyn-clone, ref-cast, ref-cast-impl, proc-macro2, quote, syn and unicode-ident. Their licenses are MIT, Apache-2.0, Unlicense and Unicode-3.0 (unicode-ident's tables).
- Build scripts (serde, serde_core, serde_json, zmij, ref-cast, proc-macro2, quote) only emit `cargo:rustc-cfg` flags. None touch the network.
- `Cargo.lock` also lists four optional toml dependencies that are never built: indexmap, hashbrown, equivalent (`preserve_order`) and toml_writer (`display`). Cargo locks optional dependencies regardless of features, so `cargo vendor` will include their sources.

**Open questions**

- §10 decision 10 (error-budget sampling and refinement) is due before phase 3, and decision 9 (event-time grid for C) before phase 4–5. Neither is in prompt 4's scope, so the numeric contract stays at version 0 at the end of phase 1.
- Library resolution: `library.paths` and the `LibraryRef` names (antennas, behaviors, waveforms) are parsed but not resolved. `sim lint`'s "references" and "unattached antennas" checks for library items wait for the library prompt.
- UPS (polar MGRS) support, if a scenario ever needs it.
- Resource headroom allowances (contracts §8) have no scenario field yet. The dry-run prompt should add one.
- Seeds are limited to 2⁶³−1 by TOML's signed integers.

**Independent review (Codex CLI, `high` reasoning, read-only)**

- It confirmed the column sets, the even-zone row offset, digit scaling, the Krüger coefficients, the Newton iteration, longitude normalization, the integer conversions and sample-rate reduction, epoch conversion, typed references, the bandwidth comparison and character columns. It found no prohibited transcendentals, HashMap/HashSet, FMA, unsafe or nondeterminism.
- **Two MGRS bugs, both fixed, with regression tests:**
  - Band validity depended on digit precision. `31NEJ0000090000` was rejected while `31NEJ` was accepted. The band test now uses the 100 km block, as GeographicLib does.
  - `31MEA0000000000` was accepted at a southern northing of 10 000 km, a block that lies north of the equator. Southern blocks must now end at or before the false northing.
- **Five minor issues, all fixed:**
  - `SampleRate` displayed `5/2 S/s`, which does not parse. It now prints exact decimals.
  - A flight level displayed as `10668 m FL`. It now prints `FLnnn`.
  - The DMS sum overflowed on `1d0m1e-100s`, and a 38-digit length in ft overflowed the exact product. Both now fall back to f64 arithmetic instead of reporting a false "too large" error.
  - The schema accepted `mgrs` together with only `lat`. Each `oneOf` branch now excludes the other form, and Python `jsonschema` confirms the rejection.
- Added a property test: every decimal-denominator rate's `Display` parses back to the same rate.

---

## Prompt 5: Work list and CPU render backend (2026-11-08)

**Built**

- `iqsim-core`:
  - `amplitude`: R₀ = 50 Ω RMS (`rms_volts`), `TargetScale` (V_fs → the f32 factor `lsb_per_volt = f32(2⁴⁰/V_fs)`, which defines the scale), `Headroom` (Σ ceil(peak LSB) ≤ 2⁶², error otherwise), `FULL_SCALE_PATHS = 2²²` and `LOSSLESS_CULL_LSB = ¼`.
  - `response`: `TapKind`, `ResponseSpec` (integer f_p, f_s, A, τ_g), `TapResponse` (Gaussian-transition linear-phase lowpass: |H|, phase, H(f), truncated closed-form h(t), support) and `TapResponseId` (SHA-256 of kind, definition, rate and centre).
- `iqsim-worklist`: `#[repr(C)]`, padding-free types with sizes pinned by a test.
  - `PolyRegisters`, `AmplitudePoly`, `PathSegment`, `EmissionRef`, `CwEmission`, `PulseTrain`, `Pulse`, `EdgeTableHeader` and `TargetDesc`.
  - `EdgeTableArena` and `BlockWorkList`, with `validate()` (every index, code, length, coordinate bound, pulse order, and the plateau shared by a pulse's two tables) and `content_hash()`.
- `iqsim-waveform::edge`: deterministic band-limited edge tables for linear ramps or steps through any tap response at any offset, with adaptive resolution and a measured interpolation-error bound. Keyed by `EdgeTableKey`.
- `iqsim-render`: the `Backend` trait and `CpuBackend` (`std::thread::scope`, contiguous chunks or a fixed chunk length dealt round-robin). The jump-then-step kernel uses a fixed f32 operation order, ties-even conversion and checked i64 adds. `edge_value` is the table lookup.
- `iqsim-geo`: `Ecef` and `SPEED_OF_LIGHT_M_S` (minimal).
- `iqsim-runtime::static_scene`: the temporary static-geometry builder. Friis with isotropic antennas, τ = R/c, the exact `phase_at` block-start phase plus −f_c·τ, the waveform coordinate from `ClockMap` instants, the CW response fold, culling, the peak-overlap sweep for pulse trains and the headroom proof. It builds one work list per target, so antenna-port and virtual-receiver lists are always separate.
- Tests (including review regressions): 22 new unit tests (core 9, worklist 5, waveform 5, render 1, geo 1, runtime 1) and 11 integration tests in `crates/iqsim-runtime/tests/static_render.rs` (16 MS/s channel at 1090 MHz, H_rx 6/7.5 MHz 80 dB τ_g = 1 µs, B_out 7/8 MHz 60 dB τ_g = 0, 30 km link). Measured:
  - **Pulse TOA** (energy centroid, group delay removed, six fractional offsets, both taps): worst 2.5·10⁻⁶ samples (bound 10⁻³).
  - **Pulse shape against an independent Simpson quadrature** (no tables, no truncation): 3.4·10⁻⁶·V at +3 MHz and 4.3·10⁻⁶·V for a 1030 MHz pulse 60 MHz out of band (bound 10⁻⁵·V). The 1030 MHz residue at −59 dB is the pulse's real spectral splatter inside H_rx, not an alias.
  - **CW level** within 1.1·10⁻⁷ relative of independent Friis·|H(f)| at five offsets (passband and transition) on both taps. Frequency is within 1 Hz, and phase (−f_c·τ − f·τ_g) within 10⁻⁴ rad.
  - **Out of band:** 1030 MHz and 1099 MHz CW are culled. With culling off, the would-be images at −4 and −7 MHz are ≤ 10⁻⁴ of the unfiltered level.
  - **Lossless culling:** identical output with and without culling a 0.03 LSB path.
  - **Determinism:** identical across 1, 2 and 8 threads and chunk lengths 777 and 1, for two blocks and both taps, with a pulse train straddling the block boundary. New golden `render/static_scene`, blessed for the first time. No existing golden changed.
  - **Headroom:** a 10⁻²⁰ V full scale is rejected. A hand-made list that bypasses the proof fails with `Overflow` (checked add) or `Contribution` (above i64), and never wraps.
  - **Tap isolation (basic):** the antenna-port work-list hash and output are unchanged across two different virtual-receiver responses and full scales. The virtual output does change.
- Docs: contracts §2 (scale factor, peak overlap, culling), new §5.1 (response family), §5.2 (work-list layout and kernel) and §5.3 (edge tables, keys, cache preparation), and §9. AGENTS.md has invariants for normalization/headroom and band-limiting before sampling.

**Deviations**

- **Response family.** A Gaussian-transition lowpass (sinc × Gaussian impulse response) instead of an unspecified filter. "Stopband attenuation A" means ≥ A dB at and beyond f_s, decreasing monotonically, not a flat floor. A raised cosine was rejected because its 1/t³ tails need ~25× longer edge tables for the same accuracy. A measured or arbitrary H (FIR taps, a receiver profile's own curve) would be a new response kind with its own identity tag.
- `fs ≤ F_s/2` is required, so the anti-alias condition holds before folding. The scenario's single `bandwidth` field is not mapped to (f_p, f_s, A, τ_g) yet. That happens when the scenario is wired to the renderer (prompt 6).
- **Edge shapes:** linear ramp (and step) only. Other shapes become new `EdgeShape` variants and keys.
- **Interpolation error** is measured at every cell midpoint against a 2× table, not proven between points (stated in contracts §5.3).
- **The delay polynomial** is in the work list for truth and diagnostics. The kernel reads time only from the waveform coordinate, and phase from the phase polynomial.
- **Placement:** the response and amplitude code went into `iqsim-core` (waveform, receiver and the builder all need it, and core is their only common dependency), and `libm` became a normal dependency of core for that control-plane math. The edge-table arena lives in `iqsim-worklist` and the builder in `iqsim-waveform`; the runtime copies tables into the arena. The static builder lives in `iqsim-runtime`. No dependency edges were added.
- **The static builder** ignores Earth rotation (no Sagnac) and light-time iteration, and accepts only clock maps that are linear over a block (midpoint residual ≤ 10⁻⁶ sample). Phase 3 replaces it.
- **Pulse offsets** convert from ps to Q32.32 target samples by floor, with the whole-sample quotient split from the remainder so no intermediate overflows. The result is within 2⁻³² samples.
- **Work-list validation** accepts only linear, non-decreasing waveform coordinates for now (contracts §5.2). Phase 3 must add an extremum check for cubic coordinates.
- The static builder rejects links with |f_c·τ| > 2⁴⁰ cycles, where f64 can no longer resolve the propagation phase to 7.7·10⁻⁴ rad. The phase 3 fitter needs a modular propagation phase for longer links.
- A CW emission is always on. Keyed CW must be rendered as pulse edges so that it stays band-limited.

**Dependencies**

- No new external crates. `libm` (vetted in prompt 3) is now a normal dependency of `iqsim-core` and `iqsim-waveform`. `sha2` (already used by core) is a normal dependency of `iqsim-worklist` for `content_hash`. `iqsim-testkit` and `libm` are dev-dependencies of `iqsim-runtime`.

**Open questions**

- With the exact `OffsetFrequency::from_uhz`, the residual frequency rounding (≤ ½ · 2⁻³² Hz) still gives up to ~6 mrad of phase over 100 days for offsets not representable exactly. It is within the 0.01 rad budget, but the phase 3 fitter should count it.
- The edge-table cache has no eviction yet. The runtime (phase 5) should drop tables whose keys no block after C references.
- The CPU kernel filters active pulses per chunk with a linear scan. That is fine for Mode A/C. Dense pulse trains may want a binary search (phase 7 performance work).
- Noise shaping by a "digital equivalent of H_rx" (prompt 6) can sample `TapResponse::response` directly.

**Independent review (Codex CLI, `high` reasoning, read-only, about 25 minutes)**

- It confirmed no padding, prohibited-transcendental, HashMap or FMA issues, and the 8 integration tests passed. It found eight issues, all fixed, each with a regression test:
  - **P1, edge-table plateau.** At 1 MS/s with a 512 MHz offset, one Gauss–Legendre rule per 15.6 ns cell undersampled the carrier, and the base-grid plateau came out 0.355 instead of ≈ 0. Refinement then copied that value into every table. Fixed with panels of at most half a cycle of |f_off| + 2f_c, the analytic H(f_off) as the shared plateau, and a quadrature convergence check against it (`EdgeError::Quadrature`).
  - **P1, pulse extent.** A rise that outlasts the fall table (rise 10 µs, fall start 1 µs, fall 0) left a non-zero tail after the extent the kernel and builder assumed, so small chunks dropped it. The extent is now the later of both table ends (`pulse_extent`), and a test renders the pulse with chunk lengths 1, 7, 64, 777 and the whole block.
  - **µHz → f64 → `OffsetFrequency`** drifted 0.0375 rad over 100 days at 7 MHz + 1 µHz. Added the exact integer `OffsetFrequency::from_uhz`.
  - **σ solve** ignored the far Gaussian tail: f_p = 1 Hz, f_s = 1 kHz, A = 20 dB gave |H(f_p)| = 0.80 instead of ≥ 0.9. It now solves against the full droop expression.
  - **`ps_to_coord`** needed 134-bit intermediates for in-range values and could panic on `abs(i64::MIN)`. The quotient and remainder are now split before the shift.
  - **Headroom** allowed `a·lsb_per_volt` to overflow f32 when a tiny pulse amplitude hid a strong path. The envelope bound is now taken as at least 1, and |pulse amplitude| is capped at 10⁶.
  - **Validation** checked only endpoints, so a curved coordinate could wrap in the middle. Coordinates must now be linear and non-decreasing until phase 3.
  - **Propagation phase** −f_c·τ in f64 loses the fractional cycle at absurd ranges (2.1 rad at 10¹⁵ m). Such links are now rejected (|f_c·τ| ≤ 2⁴⁰ cycles).
- It also noted that the determinism test's non-zero check could be met by the CW paths alone, and that both out-of-band CW tones quantized to zero even with culling off. The test now asserts that pulse trains are present in both blocks. A new case renders an 8.3 MHz tone through B_out whose −7.7 MHz image is non-zero: measured 3.8·10⁻⁷ of the unfiltered level, against a 60 dB (10⁻³) stopband.
- No golden hash changed.

---

## Prompt 6: Receiver chain, SigMF output and first validation checks (2026-11-08)

**Built**

- `iqsim-receiver` (contracts §5.4, §5.5):
  - `chain`: `ChainConfig` with one `Option` per stage, and `ReceiverChannelState`, the single owner of a channel's state. Stages in order: accumulator → volts; front-end gain (its noise figure sets the noise level); kTBF noise, a correlation of counter-indexed Philox noise with the sampled H_rx scaled to N₀ = k·T₀·F·g²; LO offset at each sample's acquisition instant with the exact modular phase (phase noise is a stub); IQ imbalance and DC; decimation (sampled §5.1 lowpass, carried history); AGC (off by default); ADC quantization and clipping, left-justified in sc16. All f32, every operation separate.
  - Blocks must arrive in order with contiguous indices; anything else is an error that leaves the state unchanged. Output is identical for any block split, down to 1-sample blocks with every stage on.
  - `port`: `PortConverter`, integer-only antenna-port conversion (sc16 `round_ties_even(acc/2²⁵)`, clipping is an error; cf32 `f32(acc)·2⁻⁴⁰`).
- `iqsim-io` (written by Codex CLI to my spec, then reviewed): `sigmf::SigmfWriter` (ci16_le and cf32_le) over a generic `CheckedSampleSink<W: Write>` that checks every write, counts samples and hashes the exact bytes. It never overwrites, and it writes `.sigmf-meta` only on a successful `finish` (flush, `sync_all`, length check). The meta holds only SigMF core fields. Also `read_recording`, `decode_ci16`/`decode_cf32`, and `utc::format_datetime` (integer, picosecond, proleptic Gregorian).
- `iqsim-contracts::manifest`: `RunManifest` (schema version 1), `BuildIdentity`, `ScenarioRecord`, `AssetKind`/`AssetRecord`/`AssetRegistry` (register with conflict detection, verify through a caller-supplied reader, no I/O), `OutputRecord`, `TruthRecord`, `EventLogRecord` (unused until prompt 16), and canonical JSON. `deterministic_sha256` excludes `backend` and `nondeterministic`.
- `iqsim-core::NUMERIC_CONTRACT_VERSION = 0`. `iqsim-geo`: `Ecef::from_geodetic` (WGS84, libm).
- `iqsim-scenario`: `[render]` (duration, block), `response` details for H_rx/B_out, the receiver-chain stages, `fixed_offset` clock and LO models, ADC and antenna-port full scales, and sink `format`. New lint checks for these. Schema regenerated.
- `iqsim-runtime`:
  - `plan`: `resolve` → `ResolvedScenario` (every default filled in, serializable, hashed into the manifest) → `RunPlan`. Receiver entity ids come from a hash of the name.
  - `offline`: `run_offline` (static work list → CPU render → chain or port → SigMF; truth record; manifest) and `replay`. `build.rs` records the rustc release, target, profile and a fingerprint of the workspace sources.
- `sim run --offline <FILE> --out <DIR> [--threads N]` and `sim replay <MANIFEST> [--threads N]`. Run directory: `consumer/<sink>/<stream>.sigmf-*`, `run/manifest.json`, `run/truth.json`.
- `examples/minimal.toml` now renders 10 ms: the CW beacon, a channel with a 30 dB / 5 dB NF front end and a −30 dBm ADC, and a new antenna port. `sim run` then `sim replay` succeeds.
- Tests: 8 chain tests and 2 port unit tests, 2 manifest unit tests, 1 geo test, 11 + 8 io tests (Codex), 9 validation tests in `crates/iqsim-runtime/tests/validation.rs`, and updated CLI and lint tests. Measured:
  - **Delay:** pulse TOA through gain and the 16-bit ADC, worst 7.6·10⁻⁵ samples (bound 10⁻³).
  - **Power and in-band level:** six CW tones from DC to the transition band on both taps against Friis·|H(f)|²: worst 0.0004 dB (bound 0.1 dB).
  - **Noise floor:** total power −0.037 dB against k·T₀·F·G·∫|H_rx|²df (B_n = 13.70 MHz), Welch passband PSD −0.036 dB (bound 0.2 dB).
  - **Out of band:** a 1030 MHz, 60 dBm emitter into the 1090 MHz, 16 MS/s channel leaves nothing at the +4 MHz fold on either tap (the path is culled), and a 1091 MHz reference keeps its level within 0.1 dB.
  - **Tap isolation:** H_rx stopband, gain, NF, noise off, seed, clock offset, LO offset, ADC bits, IQ, decimation and AGC each change the virtual-receiver hash and never the antenna-port hash. The port's own full scale and group delay do change it.
  - **Golden:** `runtime/minimal_rx_ch0`, `runtime/minimal_port` and `runtime/minimal_truth`, plus `io/sigmf_data` and `io/sigmf_meta` (Codex). All are new; no existing golden changed. The test also checks that the manifest hashes match the files and that consumer metadata mentions no truth.
  - **Truth mapping:** with decimation by 2 and a 13 ns clock offset, truth gives 8 000 output samples at 8 MS/s, the decimation delay, and the shifted acquisition instant and reported label. The antenna port keeps its own mapping.
  - **Replay:** identical across 4 and 1 threads. Rejected: an edited or removed scenario, a different numeric contract version, lockfile or source fingerprint, a tampered resolved scenario (with and without a matching rehash), a changed or missing asset, an event log, unknown manifest fields, and different output or truth hashes.

**Contracts changed** (`docs/contracts.md`, still version 0)

- New §5.4 (receiver chain, state ownership, handoff) and §5.5 (antenna-port output). §2 output conversion settled. §1.5 entity ids. §6 implemented fields, replay order, run directory, and the canonical JSON decision ("iqsim canonical JSON v1", RFC 8785 except float formatting; prompt 9/10 had been the proposed decision point). §7 SigMF fields. §9 placements. §10 power, noise and out-of-band checks marked tested. §11 new scenario fields and what `sim run` supports.

**Deviations**

- **Bandwidth rule.** Prompt 4 allowed bandwidth = sample rate. A response needs a transition band, so lint now requires bandwidth < stopband ≤ rate (stopband defaults to the rate). Two prompt 4 tests were updated to the new rule; the negative fixtures still pass unchanged.
- **Scope of `sim run`.** Library items are not resolved yet, so emitters must use the built-in `cw` / `always_on`. MSL and AGL heights are treated as ellipsoidal (no geoid or terrain), flight levels are rejected, mounts are ignored (isotropic antennas at the platform origin), and only SigMF sinks are written. Each is a resolution error, not a silent default.
- **Clock and LO configuration.** The prompt asks to vary the clock configuration in the tap-isolation test, but perturbed clock maps belong to prompt 15. I added the smallest real models: `fixed_offset` clock (an ideal map with a later epoch) and `fixed_offset` LO (constant frequency error).
- **Noise needs no carried state.** The design lists filter history among the stateful processes. Noise shaping is instead a correlation over counter-indexed white noise, so each block recomputes its L − 1 sample lookahead. That is deterministic and split-invariant without state. Decimation and AGC do carry state.
- **Decimation** starts from zero history, so its first (L − 1)/2 outputs are a start-up transient, and it adds (L − 1)/2 input samples of delay. Truth reports that delay; consumer metadata does not.
- **The manifest records no git commit.** Instead, a source fingerprint over every crate's `Cargo.toml`, `build.rs` and `src/` (plus the root manifest and toolchain file) covers uncommitted edits too, without git at build time.
- **Effective full scale.** The antenna port's sc16 code 32768 is 2⁴⁰ accumulator LSB, which is `2⁴⁰/lsb_per_volt`, within 2⁻²⁴ relative of the declared V_fs, because the f32 factor defines the scale (§2). The truth record reports it.
- **`sim replay` checks only.** It re-renders into hashing sinks and writes nothing.
- **Out of band.** In this pipeline the 1030 MHz CW is culled (|H| ≈ 0), so the alias check is trivially met. The non-culled image level was measured in prompt 5 (render level).
- `serde_json` now has `float_roundtrip` enabled workspace-wide. Without it a manifest's resolved scenario re-hashed differently after a write/read round trip; replay caught this on the first run.
- `examples/minimal.toml`: the beacon moved from 15 m to 3 km from the site, so the example is a realistic link that stays inside the antenna port's full scale.

**Independent review (Codex CLI, `high` reasoning, read-only)**

- Codex also wrote `iqsim-io` to a spec I gave it (workspace-write, limited to that crate). I reviewed it before wiring it in and changed nothing.
- The review found four issues. All are fixed, each with a test:
  - **P1, build identity.** A code change with the same version and lockfile kept the same identity, so replay would render with different code. Added `source_sha256` (above), and a replay test for it.
  - **P1, truth time mapping.** Truth reported render-rate counts and delays without saying which index space they were in, ignored decimation, and did not change with a clock offset. Truth now has a per-stream mapping (rates, counts, index factor, decimation delay, first-sample acquisition instant and reported label, effective full scale). `delay_samples` was renamed to `delay_render_samples`. The `runtime/minimal_truth` golden was re-blessed for this; it is new in this prompt, and both sample goldens are unchanged.
  - **P2, port scale.** Conversion assumed the declared V_fs/2⁴⁰ while the f32 factor defines the scale, so a value just below full scale could report a false clip. Defined and documented the effective full scale (`TargetScale::full_scale_effective_volts`).
  - **P2, ADC overdrive.** A finite but huge value (`dc_i = 1e35`) overflowed the f32 ADC multiply and failed the run instead of clipping. The input is now checked for finiteness first, and ±∞ after scaling clamps to a clip. My own pass had already fixed the related i64-range case (`i64::MAX` input).

**Dependencies**

- No new external crates. `sha2` is also a build-dependency of `iqsim-runtime`, for the source fingerprint. `serde_json` gains the `float_roundtrip` feature, which adds no packages. Workspace dependencies added to more crates: contracts (serde, serde_json, sha2, hex), io (sha2, hex, serde_json), receiver (libm; dev libm), geo (libm), and runtime (libm, serde, serde_json, sha2, hex). No internal dependency edges were added; `cargo xtask deps` passes.
- `iqsim-runtime/build.rs` only runs the configured `rustc -vV` and reads `TARGET`/`PROFILE`. It never touches the network.

**Open questions**

- The truth record is per link (geometry and level), not per emission. Parquet truth and consumer-knowledge exports arrive with the exports prompt.
- Consumer metadata should label the decimation group delay and start-up transient once the receiver-profile prompt pins how a consumer learns about them. Truth already records the delay.
- AGC is a simple multiplicative loop. A receiver-profile AGC (attack/decay in time constants, dB-domain) may replace it.
- `ResolvedScenario` is the first canonical "resolved scenario". Library and IQ-file prompts must add their items to it and register assets so replay covers them.

---

## Prompt 6a: IQ-file emitters and CPU resampling (2026-11-09)

**Built**

- **Signal model and resampler (contracts §5.6).** A sampled source is reconstructed by φ, a zero-phase member of the §5.1 family (new tap kind `SourceReconstruction`, code 3). Each target sees `z = x ∗ K`, where `K = φ ∗ h` shifted to the carrier offset (spectrum Φ(f)·H(f + Δf), including H's group delay). The carrier phase polynomial does the frequency translation at the target rate, so a low native rate never translates its own samples.
  - **Stage 1** (`iqsim-waveform::continuous`, control plane, f64, fixed order): K comes from its frequency integral (Gauss–Legendre, convergence by doubling), is truncated to supp φ + supp h, and is evaluated on U = 2^j points per native sample. U is the smallest that keeps a 6-tap Lagrange interpolator within 2·10⁻⁶ at the highest frequency z can contain. Windows are pure functions of the grid index, so block splits agree bit for bit.
  - **Stage 2** (`iqsim-render::farrow_value`, f32): 6-tap quintic Lagrange Farrow with a fixed coefficient table and operation order, recorded for the CUDA port.
- `iqsim-worklist`: `EMISSION_CONTINUOUS`, `ContinuousEmission`, `SourceWindowHeader` and `SourceWindowArena` (a new `BlockWorkList::continuous` and `windows`). `validate()` checks native units and that every Farrow tap lies inside the window. `content_hash` covers the windows.
- `iqsim-runtime::static_scene`: `StaticWaveform::Continuous` (`ContinuousSpec`, `SharedSource`). Kernels are cached by `ResampleKey`. Per block:
  - the native coordinate is exact at k = 0 (Q64.64 instant × rational rate in i128, τ rounded once);
  - links that share a source, kernel and start share one window;
  - coordinates are rebased to the window origin;
  - the level factor applied in stage 1 (f64), the `stop` exhaustion check over the samples the run keeps, the window-size limit, lossless pre-filter culling, and the headroom bound from the window peak times the Farrow Lebesgue constant.
- `iqsim-io::iq_input` (written by Codex CLI to my spec, reviewed, unchanged): SigMF metadata validation, one streaming pass (SHA-256, per-page hashes, peak, Σ|x|², optional `core:sha512`), and an LRU page-cache reader that re-verifies each page on load.
- `iqsim-contracts`: `AssetRegistry::register_digest` / `verify_digests`, so large files are hashed and verified without loading them.
- `iqsim-scenario`: `Emitter.iq_file` (`IqFileSpec`, `EndPolicySpec`, `NormalizationSpec`, `SourceResponseSpec`) and lint rules. Schema regenerated.
- `iqsim-runtime::plan` / `offline`:
  - `resolve(scenario, base_dir)` opens the recording and records exact rate, N, start, policy, normalization factor, φ and both file hashes in the resolved scenario. CW emitters serialize as before.
  - `RunPlan::new` reopens the recording and requires the same hashes.
  - The manifest registers `<emitter>.data` and `<emitter>.meta` as `iq_input` assets. Replay verifies them by streamed digest before anything else.
  - Truth gains an `iq_file` entry per IQ-file link (CW link records are unchanged).
- `examples/iq_file.toml` with `examples/iq/chirp_2msps.sigmf-{meta,data}` (10 ms ci16 chirp + tone, regenerated and compared by a test). `sim lint`, `sim run --offline` and `sim replay` work on it, and a flipped data byte makes replay fail.
- Tests: waveform 5, render 1, worklist 1, scenario 2, Codex's io 17, and 16 in `crates/iqsim-runtime/tests/iq_file.rs` (4 of them review regressions). Measured on both taps (16 MS/s channel at 1090 MHz, 2 MS/s source at +5 MHz):
  - **Frequency and level:** four tones from −0.65 to +0.7 MHz: level within 3.5·10⁻⁷ dB of Friis·|Φ|·|H| (bound 0.1 dB), frequency within 3·10⁻¹⁰ Hz (bound 1 Hz).
  - **No native-rate images:** every 2 MHz image inside the band ≤ 1.6·10⁻⁸ of the tone (bound 10⁻⁴).
  - **Fractional delay:** a chirp at two start offsets against an independent direct sum whose kernel is a time-domain Simpson convolution of the test's own sinc·Gaussian responses: ≤ 2.2·10⁻⁷ of peak (bound 2·10⁻⁵), envelope delay ≤ 7.8·10⁻⁷ samples (bound 10⁻³).
  - **Out of band:** a 1030 MHz IQ emitter is culled, and is exactly zero with culling off. A straddling source at +8 MHz: the 7.3 MHz component matches |H| within 5·10⁻⁷ dB, and the would-be −7.3 MHz alias is ≤ 2·10⁻¹⁷ (rx) and 5·10⁻¹⁵ (port) of the unfiltered level.
  - **Port power through SigMF**, ci16 and cf32, `full_scale` and `rms`: within 10⁻⁴ dB of the link budget.
  - **Determinism:** identical across 1/2/8 threads and chunk lengths 1 and 777; blocks of 16 000, 4000 and 1777 samples are bit-identical sample for sample (the test allows 10⁻⁶); `loop` repeats bit for bit every period; `zero_pad` is exactly empty after the tail; `stop` fails the block (and the run) that outlasts the file.
  - **Replay and inputs:** a changed IQ byte, a changed meta file and a missing data file are rejected before rendering. Truncated, `ci16_be`, NaN, mismatched-rate and missing recordings are rejected at resolution.
  - **Golden:** new `runtime/iq_file_rx_ch0`, `runtime/iq_file_port`, `runtime/iq_file_truth`, blessed for the first time. No existing golden changed.

**Decisions and deviations**

- **End policy `stop`** is strict, decided with the user: the recording must cover the run, and a block whose retarded coordinate passes the last sample fails. `zero_pad` is silence afterwards, and `loop` repeats with the last sample followed directly by the first. Before `start`, every policy is silent. Offline runs render the last block whole, but the check covers only the samples the run keeps.
- **Stage 1 runs on the control plane in f64**, not f32 on the render path. Only the Farrow interpolation is a render kernel. Stage 1 computes at U·F_n, which is wasteful when the native rate far exceeds the target band; a decimating stage 1 is phase 7 performance work.
- **"Filter history"** is recomputed: every window reads its own halo from the immutable source, instead of a carried state handed off between blocks. This is deterministic and split-invariant by construction, like the receiver noise in prompt 6.
- **Level normalization** is explicit and required: `full_scale` or `rms` (Σ|x|²/N over the whole file in f64, file order). The factor is applied in stage 1 in f64, not in the f32 amplitude. The headroom bound uses the actual window peak, not a declared crest factor.
- **Native rate** comes from `core:sample_rate` via its shortest round-trip decimal, parsed exactly. A scenario `sample_rate` is an assertion, not an override.
- **φ's tap kind:** φ is not a tap, but `TapKind` gained a third code so the response machinery and identity hashing could be reused unchanged.
- **Recordings are read three times** at run start (resolution, plan and, on replay, verification). This is simple and correct. Passing the opened input through is a later optimization.

**Independent review (Codex CLI, `high` reasoning, read-only, about 20 minutes)**

- Codex also wrote `iqsim-io::iq_input` to my spec (workspace-write, limited to that crate). I reviewed it before wiring it in and changed nothing.
- The review reproduced five issues with scratch programs. All are fixed, each with a regression test. No golden changed.
  - **P1, RMS silence.** A 10⁻³⁰ recording with an RMS factor of 10³⁰ was culled: the window peak was computed from f32 squares, which underflowed, and the factor sat in the f32 amplitude. The factor is now applied in stage 1 (f64), and the window peak is computed in f64 and rounded up to an f32 bound.
  - **P1, wrapping coordinate.** `validate()` checked wrapped endpoints only, so a huge linear slope could wrap mod 2¹²⁸ back into the window and then index outside it in the Farrow kernel (a panic). It now requires `D0 + (len − 1)·D1` to fit in i128. This also closes the same gap for pulse coordinates.
  - **P2, end policy missing from window grouping.** Two links on one source with `zero_pad` and `loop` shared the first link's policy. The grouping key now includes the timeline (length, policy) and the level factor.
  - **P2, window size checked too late.** An oversized window (e.g. a 16 MS/s source in 4.8 M-sample blocks, 77 M points) was allocated and filtered before the 2²⁶ limit rejected it. It is now `WindowTooLong` before any allocation.
  - **P2, `stop` and a final partial block.** The check covered the discarded tail of the last block, so a 1 ms run in 20 ms blocks rejected a 10 ms recording. `StaticSceneBuilder::with_sample_limit`, set by the offline runner, limits the check to the kept samples.

**Dependencies**

- No new external crates. `sha2::Sha512` (already vendored) is used by `iqsim-io` for `core:sha512`. No internal dependency edges were added: the runtime adapts the io reader to the waveform trait. `cargo xtask deps` passes.

**Open questions**

- Multipath (phase 8) will put several paths on one window. Grouping by (source, kernel, start) already supports it, but the window then covers the union of their delays.
- Moving geometry (phase 3) makes the coordinate cubic. `validate()` must then bound the coordinate's interior (as for pulses), and the window must cover the extremum.
- For live or generated continuous sources, stage-1 lookahead (supp K + 3 cells) must enter the commit boundary's filter lookahead (§4).
- A recording whose content reaches native Nyquist is attenuated by φ's transition (0.8–1.0 × rate by default). That is the model, not an error, but users replaying full-band captures may want a wider φ.
