# fs-gum-gpu — Phase G2 RESULTS: the wgpu f64 backend on lavapipe

Machine-checked results of the GPU rung of ROADMAP_v4_SCALE.md (Phase G2),
built on the environment facts pinned by scale_survey_v4.json DIMENSION 3.
Everything below is reproduced by `cargo run --release --bin gum_gpu_gates`
(total wall time ≈ 2.4 s on this 4-core container; budget was 15 min) and
`cargo test --release` (tier-1 naga validation, no GPU needed).

## 1. Environment proof

Reproduced the survey's headless check first (minimal adapter enumeration +
trivial f64 compute dispatch) before building anything:

```
adapter: name="llvmpipe (LLVM 20.1.2, 256 bits)" backend=Vulkan type=Cpu
         driver="llvmpipe" driver_info="Mesa 25.2.8-0ubuntu0.24.04.2 (LLVM 20.1.2)"
SHADER_F64 = true
```

* Vulkan ICD: `/usr/share/vulkan/icd.d/lvp_icd.json` (mesa-vulkan-drivers
  25.2.8-0ubuntu0.24.04.2, already installed in this container; the survey's
  CI recipe is `apt-get update && apt-get install mesa-vulkan-drivers`).  No
  `VK_ICD_FILENAMES` or other env vars needed — loader discovery is automatic.
* Crates through the proxy: wgpu 30.0.0, naga 30.0.0, pollster 1.0.1 all
  resolved and built from crates.io with zero failures (majors pinned in
  Cargo.toml; wgpu upgrades are certified-surface changes).
* f64 WGSL: naga's extension (`f64` type, `lf` literals), validating only
  with `Capabilities::FLOAT64` = wgpu `Features::SHADER_F64` — both the
  positive and the negative probe are unit tests (`tests/naga_validate.rs`).
* One survey claim did NOT generalise and is corrected here: lavapipe f64
  execution is *not* unconditionally bit-identical to same-order CPU f64.
  The trivial-dispatch probe measured 1–2 ulp differences on mul-add chains
  (llvmpipe's LLVM backend FMA-contracts, expression-dependently).  This is
  exactly the survey's own "real drivers may FMA-contract" risk, just already
  live on the rig — it confirms the tolerance-band class for every CPU-vs-GPU
  gate below.  Same-device replay is unaffected (§4).

## 2. What was built

`GpuEngine` (src/lib.rs): headless Vulkan device init preferring the
lavapipe adapter class, SHADER_F64 required at device creation, the field
device-resident as **4 per-component f64 storage buffers** in fs-gum-field's
exact padded SoA layout (P = N + 2·GHOST, offset `(ip·P + jp)·P + kp`),
ghost rind included; upload/download of raw padded data and of `Field3`
real cells; guard bookkeeping (dev/devf → w6s/cstw/e0c/rotfac, verbatim
`fs_gum_statics::engine` pass-2 dressing) on the CPU, shipped to the device
as a 16-slot f64 params buffer.

### Kernel inventory (all WGSL f64, workgroup size 64 fixed — part of the contract)

| kernel | dispatch domain | what it does |
|---|---|---|
| `fill_ghosts` | P³ | fixed-vacuum ghost rind (pure stores) |
| `forward` | N³ (1 thread/interior cell) | central4 derivs + nn/dd + 6-term PAIRS det; per-workgroup partials for (s_e2, s_e4, s_det, s_det², s_e0, s_i) via workgroup-shared reduction in fixed ascending lane order; ~KB partials buffer |
| `flux` | N³ | pointwise (2+4) flux P[3][4] + sextic/guard gq[4] (`point_flux` math, identical expression order), 16 f64/point to the intermediate buffer |
| `gather` | N³ | gather-form adjoint of the central4 stencil (transposed coefficients, ghost flux dropped), cell density terms (E0 + Routhian), h³ weight and tangent projection fused; writes the Grad-layout gradient |
| `renormalize` | N³ | pointwise unit-norm fix (IEEE sqrt) |
| `axpy_v` | 4N³ | ANF velocity kick v -= (dt/h³) g |
| `step_q` | N³ | ANF field step q += dt·v (real cells) |
| `project_v` | N³ | ANF velocity tangent re-projection v -= (v·q) q |

Final reduction of the forward workgroup partials happens on the **CPU in
fixed ascending workgroup order** (readback 81 KB at N = 48).  ANF snapshot
save/restore are device-side buffer copies; the arrest logic (one objective
readback per iteration) is the only per-iteration CPU work.

The scatter→gather rewrite is the survey's "one GPU-hostile construct"
resolved: per output cell the add order is FIXED (gq(self); axis x, y, z
with neighbour offsets −1, +1, −2, +2; cell terms; h³; projection).  It is a
deliberate reordering of the CPU engine's eval-point-ascending scatter, so
CPU-vs-GPU gradient agreement is a band gate by design, never a bit gate.

## 3. Determinism / golden policy (the Phase G2 contract, exact)

1. **Same-device replay MUST be bit-identical and is gated** (gate K2):
   every kernel writes a disjoint output location per invocation and the
   only cross-invocation reduction (forward) is in fixed lane order, so two
   dispatches of the same device state must produce byte-equal buffers.
2. **CPU-vs-GPU is the tolerance-band/metric class** (the ROADMAP_v4
   sanctioned fallback for backend boundaries, stated per gate): the shader
   compiler may FMA-contract (observed on llvmpipe at 1–2 ulp/op), and the
   workgroup reduction + gather adjoint reorder the CPU accumulation.
   Sector sums and gradient fields are gated within the stated bands below;
   bit agreement is recorded as informative only, never required.
3. Cross-device identity is NOT promised: a real GPU must reproduce ITS OWN
   bits under replay (rule 1) and land inside the same bands (rule 2); the
   adapter/driver identity is printed by gate E1 for the evidence record.

### Band derivation (reduction-order analysis)

* **Forward sums (K1, band 1e-12 relative)**: n = N³ = 110,592 terms per
  sum at N = 48.  The a-priori Higham bound on the *difference* between the
  CPU sequential order and the GPU blocked order (64-lane blocks, then 1,728
  workgroup partials sequentially) is γ_n·Σ|x|/|Σx| class: with the
  well-conditioned sums here (E2/E4/det²/E0/I positive-term; det
  single-signed on the hedgehog to ~5%), that is ≈ n·eps ≈ 1.2e-11
  worst-case, dominated by the CPU's own sequential order (the GPU's blocked
  reduction is the *tighter* order, γ ≈ (64+1728)·eps ≈ 2e-13).  Measured:
  5.7e-14 worst.  Band 1e-12 = measured × ~20 headroom, one decade under the
  coarse a-priori class; a future driver exceeding it escalates to the
  documented 1e-11 bound, not to silent tolerance creep.
* **Gradient field (K3/K3b, bands 1e-11 L2 / 1e-10 Linf relative)**: each
  output entry is a fixed ~57-term resummation (12 adjoint-stencil terms ×
  3 axes + gq + cell terms + projection) — no n-growth, per-entry class
  √57·eps ≈ 1e-15.  With both one-sided guards forced active the MU = 5000 /
  MU_F = 400 dressed weights amplify cancellation between the ±w6-scaled
  PAIRS terms; measured 9.7e-13 L2 (guarded) vs 4.2e-15 (bare).  Bands sit
  ×10–×2,000 above measured.
* **ANF endpoints (K4, bands 1e-9 relative energies/I, 1e-10 abs deg)**:
  conditional on the gated identical accept/arrest sequence (K4b), state
  drift is per-iteration 1e-14-class relative with at-most-linear growth
  over 50 iterations → 1e-12 class expected, 1.9e-14 measured.  The band is
  deliberately loose (×5e4) because it must absorb a *different-FMA* real
  driver that still makes the same 50 binary descent decisions; a driver
  that flips a decision fails K4b, which is the real gate.

## 4. Gate table (all PASS; `gum_gpu_gates`, this container, 2026-07-12)

| gate | class | band | measured |
|---|---|---|---|
| E1 device: Vulkan adapter + SHADER_F64 + 8 pipelines | hard | — | llvmpipe (LLVM 20.1.2, 256 bits), Mesa 25.2.8, init 0.1 s |
| E2 fill_ghosts: corrupted rind → vacuum, real cells untouched | **bit** (pure stores) | exact | bit-exact over 4·P³ = 562,432 values |
| K1 forward sums vs `measure()`, ε=0.05 hedgehog N=48 central4 | band | rel ≤ 1e-12 each of E2,E4,E6,E0,I,deg | worst 5.7e-14 (E0); E2 5.3e-14, E4 3.3e-14, E6 2.2e-16, I 4.1e-14, deg 2.8e-14 |
| K2 same-device replay | **bit** | exact | forward partials (10,368 f64) identical; gradient buffer (442,368 f64) identical |
| K3 gradient vs CPU analytic, bare E_static N=48 | band | L2 ≤ 1e-11, Linf ≤ 1e-10 rel | L2 4.2e-15, Linf 5.4e-15 |
| K3b gradient with anchor+wall guards ACTIVE (epen 5.1e0, efpen 5.0e-1) | band | L2 ≤ 1e-11, Linf ≤ 1e-10 rel | L2 9.7e-13, Linf 1.1e-12 |
| K4a ANF 50 it N=32: GPU objective monotone non-increasing | hard | — | 3.2116716 → 3.1952312 over 51 states |
| K4b identical descent decisions vs CPU `anf` | hard | equality | arrests 1 = 1, iters 50 = 50, status iter_cap = iter_cap |
| K4c objective checkpoints (it 0,10,…,50) | band | rel ≤ 1e-9 | worst 4.0e-14 |
| K4d endpoint obj/estat/E2/E4/E6/E0/I + deg | band | rel ≤ 1e-9, deg abs ≤ 1e-10 | worst rel 1.9e-14 (E2), deg 7.8e-16 |

Tier-1 static gates (`cargo test --release`, no GPU): all 8 WGSL kernels
parse + validate under naga `Capabilities::FLOAT64` at N ∈ {8, 32, 48}, and
the forward kernel is rejected under `Capabilities::empty()` (the f64
extension is genuinely capability-gated).

Timings (informative only — lavapipe is a correctness rig): forward
dispatch N=48 ≈ 0.01–0.06 s, flux+gather ≈ 0.05–0.15 s; ANF 50 it N=32:
GPU 1.0 s vs CPU 0.7 s.  Throughput on a software rasterizer is irrelevant
and unrecorded as evidence.

## 5. What a real GPU changes

Nothing in this crate but the adapter that `Instance::request_adapter`
returns (and the E1 evidence line).  The kernels, buffers, dispatch
topology, fixed reduction orders, CPU combine and the gate suite are
device-independent.  Specifically:

* the WGSL requires only SHADER_F64 (Vulkan `shaderFloat64`) — any
  Vulkan-class discrete GPU with f64 (and lavapipe) qualifies; Metal and
  the browser are out (no f64), per the survey.
* expected differences on real hardware: different FMA-contraction /
  scheduling → different bits *within the same bands*; K2 (its own replay)
  and K4b (decision equality) are the gates that must still hold exactly.
* scale ceilings under wgpu's default limits (128 MiB/binding): the
  per-component field buffers hold to N ≈ 251, the 4N³ gradient buffer to
  N = 161, and the 16-f64/point flux intermediate to N = 101.  Past that,
  raise limits at device request (adapter-permitting) or split the flux
  buffer per axis — a mechanical change outside the certified math.
* consumer-GPU f64 ALU ratios (1:64) temper throughput expectations for the
  sextic-heavy gradient pass (survey risk entry) — a correctness inheritance
  is claimed here, not a speedup number.

## 6. Files

* `Cargo.toml` — pinned wgpu/naga 30.0.0, out-of-workspace sandbox pattern.
* `src/shaders.rs` — WGSL f64 kernel generators (grid extents baked, all
  run-time constants through the params buffer so their bits equal the CPU
  engine's).
* `src/lib.rs` — `GpuEngine` (device init, buffers, pipelines, upload /
  download, forward combine + guard dressing, gradient, ANF device ops)
  and `anf_gpu` (the device-resident arrested-Newton-flow driver).
* `src/bin/gum_gpu_gates.rs` — the gate suite above.
* `tests/naga_validate.rs` — tier-1 static WGSL validation (no GPU).

## 7. Epistemic notice (binding)

Everything here is within-model numerical engineering on a speculative
theory's functional.  A passing gate certifies the GPU port of a
discretisation, its adjoint, and the backend plumbing on a software Vulkan
rig — it says nothing about nature, and nothing yet about real-GPU drivers
beyond the design's claim that they inherit the code and the gate suite
unchanged.  lavapipe CI green ≠ hardware certified: the first run on metal
must re-pass this suite (same bands, own-replay bit rule) before any
real-GPU result is used as evidence.
