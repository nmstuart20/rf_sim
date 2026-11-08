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
