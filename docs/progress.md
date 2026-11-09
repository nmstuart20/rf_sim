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
