# AGENTS.md — iqsim

Instructions for every coding agent (Claude Code, Codex, or others) working in this repository. `CLAUDE.md` only imports this file; edit this one.

## Project

`iqsim` is a real-time, deterministic simulator that turns RF scenarios into IQ samples at simulated receivers. The same build, resolved scenario, seed, event log and verified input assets always produce bit-identical output on any backend (CPU, CUDA, multi-GPU), identified by a run manifest. Output streams to software consumers (framed ZMQ, VITA-49, SigMF, preview) or to SDRs for open-loop HIL. Truth stays out of every consumer-visible stream.

The binary is `sim` (crate `iqsim-cli`). The reference scenario is three Mode A/C aircraft, four SSR interrogators and one receiver feeding a geolocation consumer.

## Start of every session

Read these first, in order:

1. `docs/progress.md`: what earlier sessions built, deviations and open questions. The newest entry is at the bottom.
2. `docs/contracts.md`: numeric, clock, commit-boundary, tap, manifest and transport contracts. Code must match it. Change it deliberately, and record the change in progress.md.
3. The sections of `project_plan.md` (the canonical design) that your prompt cites. `agent_prompts.md` holds the sequenced prompts.

If the design, contracts and your prompt disagree, stop and ask rather than picking one silently.

## Layout

```
crates/iqsim-core          time, fixed point, phase, sin/cos, Philox RNG, ClockMap, commit-boundary time
crates/iqsim-units         unit-checked quantities parsed from strings
crates/iqsim-contracts     versioned shared interface types: envelope, commands, acks, IQ framing, manifest
crates/iqsim-scenario      scenario file types (TOML, serde, JSON Schema)
crates/iqsim-geo           frames, trajectories, attitude, light time
crates/iqsim-antenna       patterns, polarization, mounts, blockage
crates/iqsim-terrain       DTED/GeoTIFF, line of sight, diffraction, reflection
crates/iqsim-behavior      discrete-event engine, interrogators, transponders
crates/iqsim-waveform      waveform generators, edge tables, IQ-file emitter
crates/iqsim-worklist      flat per-tap, per-block work lists (the only render input)
crates/iqsim-render        render backend trait + CPU golden-reference backend
crates/iqsim-render-cuda   CUDA backend (cudarc host side, CUDA C kernels)
crates/iqsim-receiver      virtual receiver chain and antenna-port finalizer
crates/iqsim-io            sinks (SigMF, framed ZMQ, VITA-49, preview) and exports
crates/iqsim-control       transport adapters (ZMQ) over handler traits from contracts
crates/iqsim-runtime       pipeline, pacing, commit boundary, event log, replay
crates/iqsim-hil           UHD TX/RX sinks and sources, calibration
crates/iqsim-testkit       test support: golden-hash helper (dev-dependency only)
crates/iqsim-cli           the `sim` binary
xtask/                     cargo xtask ci | bless | gpu-test | deps
tests/golden/              checked-in SHA-256 golden files
schema/                    generated JSON Schemas (checked in; `cargo xtask schema`)
examples/                  example scenarios (must lint clean)
docs/                      contracts.md, progress.md
```

The internal dependency graph is fixed in the crate manifests and checked by `cargo xtask deps`. `iqsim-control` never depends on `iqsim-runtime`, and the runtime never depends on control or HIL. Shared protocol and manifest types live in `iqsim-contracts`, which has no handlers and does no I/O. Adding an edge is a design change: record it in progress.md.

## Invariants

These are not negotiable without a design change recorded in `docs/contracts.md` and `docs/progress.md`.

- **No floating point for time, delay or carrier phase on the render path.** Truth time is i64 picoseconds. Sample indices are u64. Delay is Q32.32 in i64. Phase is a u64 fraction of a cycle. Polynomials are 128-bit Newton form.
- **Integer (i64) accumulation** with documented physical units, scale and checked headroom. Overflow is an error, never a wrap; only phase and polynomial arithmetic wrap on purpose. Emitter partitions return raw i64 sums. Receiver effects and final conversion happen once, after the global reduction.
- **Amplitude normalization** (contracts §2). Samples are complex-envelope volts at R₀ = 50 Ω, RMS (`P = |v|²/R₀`). Each render target has its own full scale V_fs, so the antenna port never depends on virtual-ADC settings. One accumulator LSB is exactly `1 / f32(2⁴⁰/V_fs)`. Every block passes the headroom proof Σ ceil(peak contribution) ≤ 2⁶² (2²² full-scale paths) before it renders. Culling below ¼ LSB is lossless; anything larger is a declared approximation.
- **Band-limit before sampling.** Every render target declares a tap response (H_rx or B_out, contracts §5.1). Narrowband paths fold H(f_off) into amplitude and phase. Pulses use edge tables keyed by response identity, shape and offset. Never render an unfiltered emitter into a target and filter afterwards.
- **No std f64/f32 transcendentals** (sin, cos, tan, exp, ln, log, atan2, powf, hypot, …) in control-plane code. Use the `libm` crate so results match across machines and libc versions. The render path uses the in-house sin/cos and ln from iqsim-core. clippy (`clippy.toml`) enforces this.
- **No fast-math, no FMA contraction (`mul_add`), no nightly features**, no `target-cpu=native`. Write multiplies and adds as separate operations. The toolchain is pinned in `rust-toolchain.toml`.
- **`#![forbid(unsafe_code)]` in every crate except `iqsim-render-cuda` and `iqsim-hil`.** Those two `deny` it at crate level and allow it only in narrowly scoped, documented FFI modules.
- **Deterministic iteration.** Outputs never depend on HashMap/HashSet order. Use BTreeMap, BTreeSet or sorted Vecs (clippy denies HashMap and HashSet). Any `#[allow]` needs a comment explaining why order cannot leak.
- **Every random draw comes from the Philox streams in iqsim-core**, keyed by (scenario seed, entity id, stream id) with the sample index as counter. No `rand::thread_rng`, no time-based seeds.
- **Sampling instants, reported timestamps and LO errors are distinct models.** Geometry and phase are evaluated at the actual acquisition instants from a `ClockMap`. The timestamp clock only labels. Use the `ClockMap` interface from the first render implementation, even while it is ideal.
- **Tap isolation.** Antenna-port rendering uses B_out. Virtual-receiver rendering uses H_rx. No virtual-receiver effect (H_rx, gain, noise, clock or LO error, ADC) reaches an antenna-port output. Cached-table keys include the tap and response identity.
- **Commit boundary.** Commands are admitted atomically at the shared exclusive commit boundary C. C covers behavior lookahead and every committed dependency. Effective time = max(requested, C). Replay uses the logged effective times and ingress sequence and never recomputes them.
- **Replay validates** the run manifest and every external input hash before rendering. Truth exports and manifests stay separate from consumer-visible metadata.
- **Delivery.** Strict delivery uses acknowledged framed streams. Raw PUB is best effort, and VITA-49/UDP only allows consumer-side loss detection. Budget aggregate traffic on every shared resource (NIC, disk, USB, PCIe).
- **Determinism scope.** Sample hashes and canonical truth hashes are deterministic. Telemetry, wall-clock diagnostics and container metadata are excluded from them.

## Dependency policy

- Use few crates, all widely used and well maintained, that vendor cleanly (`cargo vendor`) for an air-gapped RHEL 9 build. The build never touches the network: no build scripts that download, and no git dependencies.
- Justify each new dependency in your summary and progress entry: what it does, why in-house is worse, its license, and its transitive footprint. Declare versions in `[workspace.dependencies]`.
- Prefer std. Keep heavy dependencies (zmq, cudarc, uhd, parquet) inside the one crate that needs them.
- Commit `Cargo.lock`. CI uses `--locked`.

## Commands

```
cargo xtask ci        # fmt check, clippy -D warnings, all tests, rustdoc -D warnings, dependency-graph check
cargo xtask bless     # regenerate golden hashes (IQSIM_BLESS=1); extra args go to cargo test
cargo xtask gpu-test  # CUDA tests (placeholder until prompt 20)
cargo xtask deps      # dependency-graph rules only
cargo xtask schema    # regenerate schema/scenario.schema.json (a test fails while it is stale)
cargo run -p iqsim-cli -- --help
```

Property tests use proptest. They run 2048 cases per property by default; set `PROPTEST_CASES=<n>` for a longer run and `PROPTEST_RNG_SEED=<n>` to reproduce a run. Seeds are random by default. When a property fails, proptest saves the failing case under a `proptest-regressions/` directory and replays it first on every later run. Commit those files with the fix, and never delete them to make a test pass. Integration tests under `tests/` have no `lib.rs` beside them, so proptest's default persistence silently saves nothing there. Set `failure_persistence` to `FileFailurePersistence::Direct("proptest-regressions/<file>.txt")`, as `crates/iqsim-core/tests/*_props.rs` does.

Golden tests call `iqsim_testkit::assert_golden(name, bytes)`, which compares SHA-256 against `tests/golden/<name>.sha256`. Golden hashes change only through `cargo xtask bless`, and only on purpose: explain every blessed change in progress.md.

## Working rules

- One prompt per session. Plan first, then implement. Do not reach ahead into later prompts' scope.
- Several agents may work in this repo. Keep changes scoped to your task. Don't reformat or rewrite unrelated files, and don't revert changes you didn't make.
- Commit only when the human asks.

## Definition of Done

A task is done only when all of these hold:

1. `cargo xtask ci` passes.
2. New behavior has unit tests. Where `project_plan.md` names a validation check (Validation strategy table), there is a test for that check.
3. Every public item is documented (rustdoc builds with `-D warnings`).
4. `docs/contracts.md` is updated if any contract changed.
5. `docs/progress.md` has a new entry for this prompt: what was built, deviations from the design and why, new dependencies with justification, and open questions.
