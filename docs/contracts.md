# iqsim contracts

Shared contracts every crate and backend must honor. They come from `project_plan.md` ("Time, numerics and determinism", "Receiver model", "Compute backends and scaling", "Interfaces", "Validation strategy").

**Status:** numeric contract version **0 (proposed)**. Values marked *proposed* are starting points. Phase 1 (prompts 1–4) fixes them before any render code lands, and each change is recorded here and in `docs/progress.md`. The version becomes 1 when §10 is fully settled. After that, any change that can alter output bits increments it, and the version is written into every run manifest.

Contents:

1. Numeric representation and operation rules
2. Amplitude normalization, accumulation and overflow policy
3. Clocks: ClockMap, timestamp clock, LO model
4. Exclusive commit boundary
5. Tap points: responses and state ownership
6. Run manifest and replay
7. Transport choice, framing and delivery guarantees
8. Resource budgeting
9. Crate placement and dependency graph
10. Proposed tolerances and decisions phase 1 must settle

---

## 1. Numeric representation and operation rules

### 1.1 Integer quantities

| Quantity | Type | Unit / scale | Range | Overflow |
| --- | --- | --- | --- | --- |
| Truth time `TruthTime` | i64 | picoseconds from scenario epoch | ±2⁶³ ps ≈ ±106.7 days | checked; error |
| Sub-ps fraction (`AcqInstant.frac`) | u64 | 2⁻⁶⁴ ps, floor of the exact instant (§3.1) | [0, 1) ps | n/a |
| Exact index-time remainder (`ExactTime.rem`) | u64 over the reduced rate numerator | fraction of 1 ps | [0, 1) ps | n/a |
| Sample index `SampleIndex` | u64 | samples, per receiver channel, from channel epoch | 2⁶⁴ | checked; error |
| Sample rate `SampleRate` | (u64 num, u64 den), reduced, nonzero | samples per second = num/den | den ≤ 2³², num ≤ 2⁴⁸, num/den ≤ 10⁹ | rejected at construction |
| Path delay (output) | i64 Q32.32 | output samples | ±2³¹ samples (2.1 s at 1 GS/s) | range-checked at fit time |
| Carrier phase (output) | u64 | 2⁻⁶⁴ cycle | one cycle | **wraps by design** |
| Delay/phase polynomial coefficients | 128-bit, Newton forward differences D0..D3 | output scale + 64 extra fractional bits | n/a | **wraps by design** (§1.3) |
| Accumulator | i64 per I and Q | per-target LSB voltage (§2) | ±(2⁶³−1) | **forbidden**; error |

Wrapping is allowed only for carrier phase and the Newton-form polynomial registers, where arithmetic mod 2¹²⁸ (or mod 2⁶⁴ for phase) is the definition. Every other integer operation on time, index, delay output or accumulators is checked, and overflow returns an error. It never wraps or saturates silently. Release builds keep `overflow-checks = true`.

### 1.2 Time and index conversions

- Ideal map: `t(n) = channel_epoch + n · den / num` seconds, evaluated in ps as `n · den · 10¹² / num` with u128 intermediates. The quotient is the integer ps, and the remainder `(n · den · 10¹²) mod num` over `num` is the exact sub-ps fraction.
- Time → index: floor (`ChannelTimebase::index_at_or_before`: largest n with exact t(n) ≤ t; a time before the epoch is `BeforeEpoch`). A ceiling companion (`index_at_or_after`) returns the smallest n with t(n) ≥ t, and 0 for times before the epoch.
- Index → time (`ChannelTimebase::index_to_time`): returns `ExactTime { floor, rem, rem_den = num }`, exact for every index. It is an integer ps time when `rem = 0`. Converting to an `AcqInstant` floors `rem / num` onto the 2⁻⁶⁴ ps grid. The fraction is never dropped on the geometry path.
- **Rate bounds (settled, prompt 1):** after reduction, `den ≤ 2³²`, `num ≤ 2⁴⁸` and `num/den ≤ 10⁹` (1 GS/s). `SampleRate::from_ratio` rejects anything else, and zero numerators or denominators.
- **No-overflow proof (prompt 1).** Let `M = den · 10¹² < 2⁷²`.
  - Index → time: `n · den < 2⁹⁶` always. The product with 10¹² is checked. For runs ≤ 100 days at ≤ 1 GS/s, `n ≤ 8.64·10¹⁵ < 2⁵³`, so it is `< 2¹²⁵` and never overflows. If it does overflow (far larger n), the product is ≥ 2¹²⁸ and the quotient ≥ 2¹²⁸ / 2⁴⁸ = 2⁸⁰, which is not a valid i64 ps time anyway. So the check never rejects a representable result, and the error is `TimeOverflow` either way. The epoch is added in i128 and range-checked into i64.
  - Time → index: `d = t − epoch` is computed in i128, `|d| < 2⁶⁴`. `d · num < 2¹¹²`. The quotient is ≤ d / 1000 < 2⁵⁴, so it always fits in u64.
  - Remainder → fraction: `rem < num ≤ 2⁴⁸`, so `rem · 2⁶⁴ < 2¹¹²`.
  - A property test checks that `index_to_time` fails exactly when the true time does not fit in i64, for arbitrary u64 indices, epochs and accepted rates.

### 1.3 Block polynomials (delay and phase)

- Order ≤ 3 in local sample index k ∈ [0, block_len). Stored as Newton forward differences with 64 fractional bits beyond the output format: delay coefficients are Q32.96, phase coefficients are 0.128 cycle fractions.
- Jump: `p(k) = D0 + C(k,1)·D1 + C(k,2)·D2 + C(k,3)·D3` mod 2¹²⁸. Step: three 128-bit adds per sample. Jump and step agree bit for bit at every k.
- Output extraction keeps the top 64 bits with an arithmetic shift (floor toward −∞). It does not round.
- Maximum block length: `C(k,3)` must fit in u64, so k < ≈ 4.8·10⁶. Prompt 2 pins the exact bound.
- The 128-bit arithmetic is written with explicit 64-bit primitives (`mul_hi_u64`, wrapping mul, add with carry) in one module, so the CUDA C port mirrors it line by line. It is tested against an i128 reference.
- The geometry engine checks that the delay polynomial's true range over the block stays inside the Q32.32 range. Wrap-around in the delay output is a fit error.

### 1.4 Floating point

Floating point is allowed only where its rounding is fixed by construction.

- **Render path:** f32 only, for amplitude polynomials, sin/cos results, waveform samples and per-path complex products. Never f64 (T1000 FP64 is 1:32). Never time, delay or carrier phase.
- **Control plane:** f64 for geometry, fitting, link budgets and attitude. It produces integer or f32 coefficients for the work list. Transcendentals come from the `libm` crate, never std (clippy enforces this). Results must be identical on x86-64 (RHEL 9 target) and aarch64 (development). Only IEEE basic operations plus libm are allowed. No x87, no `target-cpu=native`.
- **Operation order:** the source order is the order. No reassociation, no fast-math, no FMA contraction. `mul_add` is banned. CUDA kernels build with `--fmad=false`.
- **Rounding mode:** IEEE 754 round-to-nearest-even for all arithmetic.
- **Division and sqrt:** correctly rounded IEEE on both backends. CUDA uses `-prec-div=true -prec-sqrt=true`.
- **Subnormals:** *proposed*: preserved (IEEE gradual underflow) on both backends. CUDA uses `-ftz=false`. Prompt 2 confirms the cost or switches both backends to a documented flush rule.
- **Float → integer:** ties-to-even, never Rust `as` (which truncates and saturates). Use `round_ties_even()` followed by a checked range conversion, and `__float2ll_rn` on CUDA. Out-of-range is an error.
- **Non-finite values:** NaN and ±Inf are rejected at every input validation and at every float → integer boundary. A non-finite value never reaches an accumulator.
- **Compiler settings:** the toolchain is pinned (`rust-toolchain.toml`). No nightly features. CUDA fatbins are built offline for sm_75 and sm_86 with the flags above, and their hashes go in the run manifest.

### 1.5 sin/cos, ln, Gaussian

- `sincos(phase: u64) -> (f32, f32)`: the top phase bits index a table built by a documented deterministic routine, and a fixed-order polynomial corrects the remainder. Every operation is written explicitly. Target max absolute error ≤ 3·10⁻⁷ (prompt 2). The table's bytes have a golden hash.
- `ln` for Box-Muller: in-house range reduction plus a fixed polynomial, with no std or libm call, so the GPU port matches. Its error is documented in prompt 3.
- Noise: Philox4x32-10. Key = deterministic mix of (scenario seed: u64, entity id: u32, stream id: u32). Counter = (sample index, sub-index). The Gaussian draw is Box-Muller, with the angle passed as a u32 phase to `sincos`. Stream ids are allocated in a single registry in iqsim-core, so two models never share a stream.

### 1.6 Carrier-offset phase

Offset phase `φ(t) = frac(Δf · t)` comes from exact time with modular arithmetic. Δf is an exact fixed-point or rational frequency (*proposed*: i64 in units of 2⁻³² Hz). t is integer ps plus a sub-ps fraction. The product is reduced mod 1 cycle in 128-bit integer arithmetic, then scaled to 2⁻⁶⁴ cycle. Nothing ever computes a large floating-point `f · t` and discards its integer cycles. Precision near 100 days must match precision near t = 0 (prompt 2 test).

## 2. Amplitude normalization, accumulation and overflow policy

**Physical normalization (*proposed*).** A sample value is the complex-envelope voltage at a reference impedance R₀ = 50 Ω, using the RMS convention `P = |v|² / R₀`. A CW tone of power P has `|v| = √(P·R₀)`. Path amplitude polynomials carry path loss, antenna gains, polarization mismatch and the selected tap response. They are in volts at the render target's reference plane.

**Per-target integer scale.** Each render target (one tap of one channel) declares:

- reference plane: antenna port of the receiver, or the virtual ADC input
- full-scale voltage `V_fs`
- LSB voltage `q = V_fs / 2^S`, with *proposed* S = 40

Contribution component value = `round_ties_even(component_volts / q)` as i64. The antenna-port scale is independent of any virtual ADC full scale.

**Headroom proof.** Before rendering, and again whenever the path set or amplitude bounds change, the control plane computes, per target and block:

```
Σ_paths  ceil(max|a_p| · max|s_e| / q)  ≤  2⁶³ − 1 − guard
```

`max|a_p|` comes from the block's amplitude polynomial bound. `max|s_e|` is the waveform's declared peak (pulses: edge-table peak; continuous: declared crest factor). *Proposed* guard: 2⁶⁰ of slack, about 3 guard bits. With S = 40 and contributions bounded by full scale, about 2²² full-scale paths fit per sample. A configuration that fails the proof is rejected (`sim lint` and run start). Within the bound, integer addition is associative, so partitioning and thread order cannot change the sum. The CPU backend also uses checked adds, so a proof bug fails loudly instead of wrapping.

**Reduction and finalization.** Emitter partitions return raw i64 I/Q sums. The host, or one designated device, reduces them exactly. Receiver effects and output conversion run once on the complete sum (§5).

**Output conversion.** sc16 or cf32 per sink. i64 → sc16 is `round_ties_even` of the scaled value with clip detection. Clipping is a reported receiver effect on the virtual ADC. On the antenna port, clipping is an error.

**Amplitude tolerance and null floor (*proposed*, phase 1 decides).** Each link declares a reference voltage `V_ref`, its nominal main-beam amplitude.

- Where `|a| ≥ floor`, the error must be ≤ 0.01 dB relative. *Proposed* floor = `V_ref · 10⁻²`, i.e. −40 dB.
- Below the floor, the absolute complex-gain error must be ≤ `10⁻³ · V_ref`.
- Phase is tested only above the floor.

## 3. Clocks: ClockMap, timestamp clock, LO model

Three distinct models. No code may use one as a stand-in for another.

### 3.1 ClockMap (sampling oscillator)

A ClockMap maps a sample index to the actual truth-time acquisition instant. Geometry, propagation delay and carrier phase are evaluated at these instants.

```text
trait ClockMap {
    /// Actual acquisition instant of sample n. Strictly increasing in n.
    fn instant(&self, n: SampleIndex) -> Result<AcqInstant>;
    /// Smallest n with instant(n) >= t (and the companion "largest n with instant(n) <= t").
    fn index_at_or_after(&self, t: AcqInstant) -> Result<SampleIndex>;
    /// Stable identity (content hash of the model + parameters) for cache keys and manifests.
    fn id(&self) -> ClockMapId;
}
AcqInstant = { ps: TruthTime, frac: u64 }   // ps + frac·2⁻⁶⁴ ps, i.e. signed Q64.64 ps (settled, prompt 1)
```

- **Ideal map:** the exact rational map of §1.2, floored onto the 2⁻⁶⁴ ps grid (`IdealClockMap`). This is the v1 implementation, and the first render code uses the trait, not a shortcut. Its identity is SHA-256 of `"iqsim.clockmap.ideal.v1\0" ‖ epoch_ps (i64 LE) ‖ num (u64 LE) ‖ den (u64 LE)`.
- **AcqInstant representation (settled, prompt 1):** `ps + frac · 2⁻⁶⁴` picoseconds, floor of the exact instant. A fixed binary grid works the same for the ideal and perturbed maps, and it is a plain 128-bit integer for the GPU port. The quantization error is < 2⁻⁶⁴ ps ≈ 5.4·10⁻²⁰ ps. At 12 GHz that is < 5·10⁻²¹ rad of carrier phase, so it uses none of the §10 phase or delay budget. Rounding keeps the map strictly increasing, because samples are ≥ 1000 ps apart.
- **Inverse:** `index_at_or_after(t)` returns the smallest n with `instant(n) ≥ t`, and `index_at_or_before(t)` the largest n with `instant(n) ≤ t` (`None` before sample 0). Both compare against the floored instants, so `index_at_or_after(instant(n)) = index_at_or_before(instant(n)) = n` exactly. The inverse returns only the index. A caller that needs the sub-sample offset computes `t − instant(n)` exactly in Q64.64. When no representable sample is at or after t, `index_at_or_after` is `TimeOverflow`.
- **Perturbed maps (prompt 15):** frequency offset, drift, aging, optional GPSDO discipline, and phase noise driven by Philox. They are deterministic and strictly monotone, so the inverse is well defined (integer search over a monotone function). Cumulative random processes are integrated on the control plane in block order and frozen into per-block coefficients. Evaluating a block never depends on partitioning or on which device renders it.
- **Coherent channels:** the M channels of one receiver share one oscillator. Each channel's index → instant goes through the shared oscillator time at shared physical times, so channels with different rates stay coherent.
- **Antenna-port tap:** uses the declared output-device timebase, which is the ideal map at the device's nominal rate. It bypasses every virtual clock impairment.

### 3.2 Timestamp clock

Maps an acquisition instant to the time reported to the consumer (stream headers, VITA-49 timestamps, SigMF `core:datetime`). It has its own offset and drift. It relabels only and never moves sampling instants. Truth exports record the truth instant, the reported timestamp and the sample index together.

### 3.3 LO model

Controls mixing phase and frequency error, i.e. the `e^{j2π(f_c − f_LO,r)t}` term and phase noise, evaluated with §1.6 modular arithmetic at actual acquisition instants. Coupling to the sampling oscillator (common reference) is explicit in the receiver profile. It is never implied.

## 4. Exclusive commit boundary

- **Horizons.** Render progress R (per channel, in sample indices) and the commit boundary C (truth time) are distinct. Telemetry reports R, C and command latency separately.
- **Meaning of C.** C is exclusive. Every truth-time dependency before C is immutable. C covers everything already consumed by behavior (including lookahead), geometry, work-list construction, receiver processing and filter lookahead.
- **Reservation.** A stage reserves its dependency interval `[start, end)` before advancing. Advancing C rounds fractional endpoints upward to the event-time grid (*proposed* grid: 1 ps).
- **Serialization.** Command admission and stage advancement are serialized. In v1 this is one commit actor, with no rollback.
- **Admission.** A command with requested time `t_req` gets effective time `t_eff = max(t_req, C)`. It is appended to the durable event log, fsync'd, and only then acknowledged. The ack carries `t_eff` and the ingress sequence. All of this completes before any stage advances across `t_eff`. Every command pending at a boundary is ordered before the advance.
- **Ordering.** Behavior events are ordered by `(truth time, entity id, sequence)`. *Proposed*: an operator command at time t precedes entity events at t, and commands are ordered among themselves by ingress sequence.
- **Latency.** Behavior lookahead counts toward the same latency budget as pipeline depth. *Proposed* defaults: block 1 ms, depth 50–100 ms.
- **In-flight signals.** A signal already emitted keeps its emission-time parameters. A command changes future emissions only.
- **Log records.** Start, stop, accepted, rejected (with reason), duplicate (returns the original outcome) and pending-at-stop each have an explicit record. Every record holds: requested time, effective time, ingress sequence, client/session id, request id and outcome. Stop may drain committed output and reports the final output boundary per channel.
- **Replay.** Replay applies logged effective times and ingress order directly. It never recomputes them from replay speed or wall clock.

## 5. Tap points: responses and state ownership

Every receiver channel has two taps.

| | Antenna-port tap | Virtual-receiver tap |
| --- | --- | --- |
| Purpose | Drive real hardware (HIL TX) | Feed software consumers |
| Pre-sampling response | B_out(f), output band selection for the DAC/output channel | H_rx(f), analog front end + anti-alias |
| Timebase | output-device timebase (ideal) | receiver ClockMap (§3.1) |
| Receiver effects | none | chain 1–7: gain/NF, kTBF noise, LO offset + phase noise, IQ imbalance + DC, decimation, AGC, ADC quantization/clip |
| Scale | its own `V_fs` (calibrated output) | virtual ADC full scale |

- **Response identity.** `TapResponseId` = content hash of (tap kind, response definition, output rate, center frequency). Work lists, edge tables and narrowband response samples are keyed by it. Changing any virtual-receiver setting cannot change an antenna-port key. This is the tap-isolation check.
- **Separate scenes.** Antenna-port and virtual-receiver work lists are separate whenever their response, clock map, rate, center or scale differ. Sharing is an optimization, allowed only when all of those are identical.
- **Response application.** Narrowband emitters: the response at the emitter offset folds into the amplitude and phase polynomials, and paths below the in-band threshold are culled. Analytic pulses: complex band-limited edge tables (edge ⊛ response, shifted to the offset), built deterministically before their effective boundary. Wideband emitters: filtered inside the resampler before the fractional delay.
- **Stateless kernels.** Render kernels see only the work list and immutable tables. They own no state across blocks.
- **Stateful processes.** Decimation filter history, AGC, phase noise, clock processes and resampler source history live in one owner per channel (`ReceiverChannelState`, or the waveform source on the control plane). They are handed off in block order: block b+1 starts from the state after block b. Parallel or multi-GPU rendering never resets or duplicates that state. Finalization of a channel happens on exactly one owner.
- **Noise.** Thermal noise is added once, after the complete path sum. Its Philox stream is keyed by (seed, receiver entity, channel noise stream), with counter = sample index at the rate where noise is injected.
- **Truth.** Truth and consumer-knowledge exports leave from the control plane. They never enter IQ samples, IQ headers or context packets.

## 6. Run manifest and replay

A versioned document written alongside every run. It is never part of consumer-visible output.

**Fields (initial).**

- `manifest_schema_version`, `numeric_contract_version`.
- **build:** package version, git commit and dirty flag, rustc version, target triple, profile, enabled features, `Cargo.lock` SHA-256, CUDA toolkit version and fatbin SHA-256s when the GPU is used.
- **scenario:** the resolved configuration in canonical serialization with its SHA-256, plus the source file paths and their SHA-256s. Canonical form: *proposed* RFC 8785 JSON canonicalization, decided when the manifest is implemented (prompt 9 or 10).
- `seed`, run mode (offline or real time), block length and pipeline depth.
- **backend and partitioning:** recorded for diagnostics. By contract they cannot change sample or truth hashes.
- **event log:** path, SHA-256, record count, first and last ingress sequence, final C.
- **assets:** one entry per external input: `{kind, logical name, resolved path, byte length, SHA-256, format/version}`. Kinds: library item, terrain tile, geoid grid, antenna pattern, IQ input, calibration table, live recording, live timing/gap record. Live recordings and their timing/gap records become assets when finalized.
- **outputs:** per stream: stream id, tap, sample format, first and last sample index, sample count, SHA-256 of the sample bytes, finalization status. Also the canonical truth-record SHA-256.
- **nondeterministic section** (excluded from every deterministic hash): wall-clock start and end, host, telemetry summary, performance counters, container-file metadata.

**Replay rules.** `sim replay <manifest>` verifies the build, numeric contract version, resolved scenario, event log and every asset hash before rendering anything. Missing or changed content is rejected. Different builds are not replay. They are compared against stored golden outputs with stated tolerances. Hash algorithm everywhere: SHA-256.

## 7. Transport choice, framing and delivery guarantees

Handlers and framing are independent of transport adapters (§9). The adapters are thin, so an alternative transport can be added later without touching handlers.

**Decision for v1 and alternatives** (from the design's transport table):

| Use | v1 choice | Guarantee | Alternatives kept open |
| --- | --- | --- | --- |
| Commands and status | ZMQ ROUTER (sim) / DEALER (clients) | acknowledged, deduplicated, logged | HTTP/WebSocket (phase 12, same handlers); gRPC if a consumer requires it |
| Operator telemetry | ZMQ PUB | best effort, sequence-numbered, snapshot resync | — |
| Strict software IQ | framed ZMQ ROUTER/DEALER over TCP or IPC | end-to-end acknowledged or the run fails | plain framed TCP, shared memory, DDS: revisit for a concrete integration or measured bottleneck |
| Receiver-interop IQ | VITA-49.2 over UDP | paced production; the consumer detects discontinuities | — |
| Preview | raw ZMQ PUB, sc16/fc32 (GNU Radio ZMQ SUB compatible) | best effort, no header, losses uncounted | — |
| Recording | SigMF files | checked writes, finalized count and hash | — |

This is a design choice to validate with the real consumer and target host, not a throughput claim.

**Control envelope.** JSON: `{schema_version, type, client_id, session_id, request_id, payload}`. The JSON Schema is generated from the Rust types and published.

- The ROUTER uses `ROUTER_MANDATORY`, bounded queues (HWM) and send deadlines. Delivery failure is reported.
- Ack = after durable log append. It returns `t_eff`, the ingress sequence and the outcome.
- Retries are deduplicated by `(client_id/session_id, request_id)`. A retry returns the original outcome and never schedules the action again.
- Clients can query command and run status.
- Telemetry carries sequence numbers. A status/snapshot request resynchronizes a client.

**Strict framed IQ.**

- Session per configured receiver-channel stream, with a dedicated queue per consumer. Fan-out is explicit.
- Handshake: consumer READY.
- Binary message header (versioned): stream/session id, block sequence, first sample index, sample count, datatype and byte order, receiver-reported first-sample timestamp (§3.2), context version. Context carries RF reference frequency, sample rate, bandwidth and gain. No truth fields and no manifest.
- Consumers send cumulative acknowledgements of contiguous blocks. The window is bounded in outstanding bytes and blocks. End-of-stream is an acknowledged final boundary.
- Queueing a message is not delivery.
- Real time: a disconnect, a missing ack, or backpressure past the configured deadline fails the sink and the run. Reconnect never silently resumes a stream as successful.
- Offline: the sender may wait up to its configured operation timeout.
- A reference receiver and a GNU Radio adapter are supplied.

**VITA-49 over UDP.** The interoperability profile is pinned when implemented (phase 5): VITA 49.2 signal data packets with stream id, 4-bit packet count, integer and fractional timestamps (units pinned with the profile, *proposed* UTC seconds + picoseconds), and payload sample count. Context packets carry RF reference frequency, sample rate, gain and bandwidth at a pinned cadence and on every change. MTU is pinned (*proposed* 1500 bytes, jumbo optional). UDP send success proves nothing about delivery. The supplied receiver checks continuity from packet count, timestamp and payload count together, handling count wrap. It reports loss, reordering and duplicates, and fails strict consumption on any gap.

**SigMF.** Sequential writes with every error checked. The sample count and SHA-256 are finalized at close. Truth annotations go to a separate sidecar file, never into the consumer's `.sigmf-meta`. In real time, recording shares the sink deadline and resource budget. Offline recording is unpaced.

**Failure classes** (kept distinct in status, logs and exit codes):

1. Producer underrun: hard error.
2. Strict-sink deadline or delivery failure: hard error.
3. Preview loss: allowed and never blocks required sinks.
4. UDP loss detected by the consumer: reported by the consumer.

## 8. Resource budgeting

- Every sink and device declares resource identities (NIC, disk, USB controller, PCIe link/root) with usable capacities.
- `sim lint` and `sim dry-run` sum concurrent traffic along every shared resource: payload, packet/framing/ack overhead, fan-out copies, and i64 emitter-partition transfers (16 B per complex sample per partition).
- Each resource gets a configured headroom allowance.
- Rough payload limits: sc16 about 30 MS/s per 1 GbE and 300 MS/s per 10 GbE. fc32 halves both.
- Two 20 MS/s sc16 sinks on one 1 GbE link need 1.28 Gb/s before overhead, so they fail the real-time check even though each fits alone.
- Unknown capacity is reported as *unverified*, never as passing.

## 9. Crate placement and dependency graph

Shared contracts sit low in the graph, so no runtime ↔ control cycle can form:

```text
iqsim-core ─┬─ iqsim-units ─┬─ iqsim-contracts ─┬─ iqsim-control      (adapters → handler traits)
            │               │                   ├─ iqsim-io ── iqsim-hil
            │               │                   ├─ iqsim-receiver
            │               │                   ├─ iqsim-scenario ── iqsim-behavior
            │               │                   └─ iqsim-runtime     (implements handler traits)
            │               ├─ iqsim-geo ── iqsim-antenna, iqsim-terrain
            │               └─ iqsim-waveform
            └─ iqsim-worklist ── iqsim-render ── iqsim-render-cuda
iqsim-cli = runtime + control + hil + scenario + contracts (wiring only)
```

- **iqsim-core:** TruthTime, SampleIndex, SampleRate, `ClockMap` trait and ideal map, commit-boundary time types, fixed point, sin/cos, Philox, `NUMERIC_CONTRACT_VERSION`.
- **iqsim-contracts:** pure data and pure encode/decode, with no sockets, files or simulation logic. It holds:
  - the control envelope, commands, acks and status
  - the handler traits (`CommandHandler`, `StatusProvider`)
  - telemetry records
  - the IQ frame header and context, and the VITA-49 profile constants
  - the run-manifest schema
- **iqsim-control:** ZMQ adapters that decode envelopes and call handler traits.
- **iqsim-io:** sinks that use contracts framing.
- **iqsim-runtime:** implements the handlers and owns the commit actor, pipeline and replay. It never depends on control or HIL.
- **iqsim-cli:** wires runtime to adapters.
- `cargo xtask deps` enforces the forbidden edges.

## 10. Proposed tolerances and decisions phase 1 must settle

**Proposed tolerances**

| Item | Proposed | Settled by |
| --- | --- | --- |
| Carrier phase error (above null floor) | ≤ 0.01 rad (**settled**). Time representation contributes < 10⁻²⁰ rad, so the full budget goes to fitting and polynomial quantization. | prompt 1; prompt 2 shows the quantization share |
| Delay error | ≤ 0.001 samples (**settled**). Time representation contributes < 2⁻⁶⁴ ps. | prompt 1; prompt 2 shows the quantization share |
| Acquisition-instant representation error | < 2⁻⁶⁴ ps, floor (**settled**) | prompt 1 |
| Amplitude error (above null floor) | ≤ 0.01 dB | prompt 2 |
| Null floor | V_ref · 10⁻² (−40 dB) | prompt 2, before prompt 5 |
| Absolute complex-gain error below floor | ≤ 10⁻³ · V_ref | prompt 2 |
| sincos max abs error | ≤ 3·10⁻⁷ | prompt 2 |
| Power check | ±0.1 dB vs closed-form link budget | phase 2 |
| Noise floor check | ±0.2 dB vs kTBF·NF | phase 2 |
| Diffraction check | ±0.1 dB vs independent P.526 | phase 8 |
| Out-of-band in-band level | ±0.1 dB vs selected H_rx/B_out | phase 2 |
| Determinism | identical hashes (CPU, GPU models, split/unsplit) | always |

**Decisions phase 1 must settle** (record the outcome here, then bump to contract version 1):

1. **Settled (prompt 1).** `SampleRate` bounds: den ≤ 2³², num ≤ 2⁴⁸, rate ≤ 1 GS/s, enforced at construction after reduction. The u128 no-overflow proof is in §1.2.
2. **Settled (prompt 1).** `AcqInstant` is a fixed 2⁻⁶⁴ ps fraction (floor). `ChannelTimebase::index_to_time` also returns the exact rational remainder. The inverse map returns an index only, and the caller gets the offset by exact Q64.64 subtraction (§3.1).
3. Exact maximum block length and coefficient-quantization error bounds for delay and phase (prompt 2).
4. sin/cos table size, polynomial order and table-construction routine (prompt 2).
5. Subnormal policy (preserve vs documented flush) on CPU and CUDA (prompt 2).
6. Carrier-offset frequency representation for modular phase (prompt 2).
7. Physical normalization (R₀, RMS convention), integer scale S and guard size; null floor and V_ref definition per link (prompt 2, before prompt 5).
8. Philox key-mixing function and stream-id registry (prompt 3).
9. Event-time grid for rounding C, and command vs entity-event ordering at equal times (before phase 4–5).
10. Error-budget sampling and refinement policy, including derivative bounds vs adaptive refinement and maximum depth (before phase 3).

Later phases pin these with the implementation: the manifest canonical encoding, the VITA-49 profile (timestamp units, MTU, context cadence) and the strict-sink window defaults.
