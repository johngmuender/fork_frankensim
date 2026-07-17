# fs-gum-gpu — K3 RESULTS: full-pipeline WGSL f64 verification battery (llvmpipe)

Workstream K3 (ROADMAP_v7_HARDENING.md, the last scheduled item): the wgpu
f64 on-ramp verified for the FULL physics pipeline — sector measure
(E2, E4, E6, E0, I) and the E_static gradient — WGSL f64 on llvmpipe,
tolerance-band against the fs-gum-statics CPU implementation at
N = 24, 32, 48, 64.  Reproduced by
`cargo run --release --bin gum_gpu_k3_gates` (12.5 s wall on this 4-core
container; writes `k3_results.json` next to Cargo.toml) and
`cargo test --release` (tier-1 naga validation, no GPU).
**All 20 gates PASS** (this container, 2026-07-17).

## 0. Scope and relation to Phase G2

Phase G2 (RESULTS.md) already ported the full kernel set (forward sweep,
flux + gather adjoint, renormalize, ANF device ops) and gated it at a
single grid (N = 48 sums/gradient, N = 32 ANF).  K3 is the readiness
battery on top of the UNCHANGED kernels and device setup (reused verbatim —
no shader or engine change was needed; the only source additions are the
new gate binary `src/bin/gum_gpu_k3_gates.rs` and two extra grid sizes in
the tier-1 naga test):

* sector + gradient GPU-vs-CPU at N = 24, 32, 48 **and 64** (llvmpipe
  wall-time allowed it: the whole battery is 12.5 s);
* gradient agreement reported as **RMS and max-abs relative difference
  over all sites/components** (K3 task metric; G2 reported L2/Linf —
  L2 rel ≡ RMS rel, Linf rel ≡ max-abs rel, so the numbers are the same
  metric class);
* a **directional-derivative check on BOTH backends** (FD of E_static
  vs ⟨grad, dir⟩, each backend's own forward path);
* **golden FNV-1a-64 fingerprints** of every GPU output buffer
  (within-backend bit goldens);
* an arrested-descent smoke stepping with the **GPU gradient**, energy
  monotone under the **CPU evaluator** (G2's K4 gated monotonicity under
  the GPU evaluator; K3 closes the cross-evaluator loop);
* the per-N wall-time record with the llvmpipe honesty statement.

Regression state after the extension: the Phase G2 suite
(`gum_gpu_gates`) re-run green (ALL GATES PASS, 1.4 s warm), and the
tier-1 naga test now validates all 8 kernels at N ∈ {8, 24, 32, 48, 64}
(2 tests pass, no GPU).  The shipped `k3_results.json` is run 1 (cold
timings, matching §7); re-runs legitimately rewrite timing fields only.

Protocol constants: LBOX = 4.5, Scheme::Central4, t = T_FROZEN =
0.008276434949296802, hedgehog seeds from the frozen Step-1 radial solve
(radial_solve(T_FROZEN, 4000, 6.0)); E_static objective
`Opts::estatic` (obj = t·E2 + t·E4 + E6 + E0 = E_static); all seeds fixed
(DD base 11, DD dirs 7000..7002, AD bump 424242).

## 1. Environment proof

```
adapter: name="llvmpipe (LLVM 20.1.2, 256 bits)" backend=Vulkan type=Cpu
         driver="llvmpipe" driver_info="Mesa 25.2.8-0ubuntu0.24.04.2 (LLVM 20.1.2)"
SHADER_F64 = true   (required at device creation; init fails without it)
```

Vulkan ICD `/usr/share/vulkan/icd.d/lvp_icd.json`; wgpu = 30.0.0,
naga = 30.0.0 (pinned).  Device setup is `GpuEngine::new` from Phase G2,
reused verbatim.

## 2. Sector measure GPU-vs-CPU, per N (gates S24–S64: PASS)

GPU: central4 forward sweep, one thread per interior cell, per-workgroup
(WG = 64) partials accumulated by lane 0 in fixed ascending lane order,
final combine on the CPU in fixed ascending workgroup order
(tile-canonical reduction — the within-backend determinism device).
CPU: `fs_gum_field::measure()` (plain sequential sweep).  Same
interior/wall handling: 2-layer fixed-vacuum ghost rind on both backends
(`fill_ghosts` kernel ≡ `fill_ghosts_vacuum`), stencil never reads past
the rind.

Measured relative differences, eps = 0.05 hedgehog (deg is the sixth
tracked quantity; E2 = (1/4π)∫|D_i q|² etc. per the frozen conventions):

| N  | E2      | E4      | E6      | E0      | I       | deg     | worst   |
|----|---------|---------|---------|---------|---------|---------|---------|
| 24 | 1.1e-15 | 1.3e-15 | 3.6e-15 | 1.1e-14 | 5.1e-15 | 5.2e-15 | 1.1e-14 |
| 32 | 8.7e-15 | 2.2e-16 | 5.6e-15 | 6.2e-15 | 1.0e-14 | 2.0e-15 | 1.0e-14 |
| 48 | 5.3e-14 | 3.3e-14 | 2.2e-16 | 5.7e-14 | 4.1e-14 | 2.8e-14 | 5.7e-14 |
| 64 | 1.1e-13 | 4.2e-14 | 3.6e-15 | 3.3e-14 | 8.3e-14 | 4.2e-14 | 1.1e-13 |

**Band: rel ≤ 1e-12 per sector (PASS at every N).**  Justification from
the measured trend: the difference is reduction-order + FMA-contraction
noise; it grows roughly like the sequential-order γ_n term of the CPU sum
(n = N³ terms: measured worst 1.1e-14 → 1.1e-13 from N = 24 → 64,
consistent with ~n·eps·(Σ|x|/|Σx|) growth dominated by the CPU's own
sequential order — the GPU's blocked order is the tighter one, γ ≈
(64 + n/64)·eps).  The 1e-12 band keeps ≥ ×9 headroom at N = 64 and one
decade under the coarse a-priori bound ≈ 6e-12 at N = 64; extrapolating
the same growth, the band holds to N ≈ 190-class grids before a
principled re-derivation (not a silent widening) would be required.

## 3. E_static gradient GPU-vs-CPU, per N (gates G24–G64: PASS)

GPU: gather-form adjoint — pointwise flux pass (P[3][4], gq[4], exact
`point_flux` expression order) then a gather kernel where **each thread
computes its own site's gradient** (transposed central4 stencil, ghost
flux dropped, cell density terms, h³ weight, tangent projection fused; no
scatter races by construction).  t = T_FROZEN weighting via
`Opts::estatic`.  CPU: the analytic scatter adjoint of
`fs_gum_statics::eval`.  All 4N³ sites/components compared:

| N  | entries   | RMS rel  | max-abs rel |
|----|-----------|----------|-------------|
| 24 | 55,296    | 1.73e-15 | 3.16e-15    |
| 32 | 131,072   | 2.47e-15 | 4.54e-15    |
| 48 | 442,368   | 4.22e-15 | 5.35e-15    |
| 64 | 1,048,576 | 7.66e-15 | 1.07e-14    |

(RMS rel = ‖g_gpu − g_cpu‖₂/‖g_cpu‖₂; max-abs rel = max|Δ|/max|g_cpu|.)

**Bands: RMS ≤ 1e-11, max-abs ≤ 1e-10 (PASS at every N).**  Each entry is
a fixed-length (~57-term) resummation — per-entry √57·eps ≈ 1e-15 class
with no n-growth, matching the flat measured profile; the mild rise with
N is the max over more entries plus FMA-contraction variance.  Bands are
the G2 bands (≥ ×1,000 headroom here), sized to absorb a real driver's
different contraction pattern; guard-dressed weights (the harder case,
×230 amplification) were gated in G2 K3b and are not re-run here.

## 4. Directional-derivative check (gate DD: PASS)

N = 32, base = hedgehog + seeded bump (seed 11, amp 0.05, renormalized);
3 seeded random tangent directions (seeds 7000–7002: `bump_field` →
tangent-project → L2-normalize); central FD of E_static at eps = 1e-5
(the statics G-A protocol), FD evaluated by each backend's OWN forward
path; ⟨g, u⟩ in fixed CPU order from each backend's own gradient:

| dir seed | FD-vs-⟨g,u⟩ res (CPU) | FD-vs-⟨g,u⟩ res (GPU) | cross ⟨g_gpu,u⟩ vs ⟨g_cpu,u⟩ |
|----------|----------------------|----------------------|------------------------------|
| 7000     | 2.10e-7              | 2.10e-7              | 2.66e-15                     |
| 7001     | 7.58e-7              | 8.40e-9              | 3.79e-15                     |
| 7002     | 2.74e-6              | 3.67e-7              | 1.74e-14                     |

**Bands: FD residual ≤ 1e-5 on both backends; cross-backend ⟨g,u⟩ rel ≤
1e-9 (PASS).**  The residuals are FD-truncation/roundoff dominated
(O(eps²) + O(machine/eps) at eps = 1e-5), identical class on both
backends — the analytic gradients themselves agree to 1.7e-14 worst in
the directional projection.  That the CPU FD residual is sometimes larger
than the GPU one (7002) is FD noise, not a backend asymmetry.

## 5. GPU determinism + golden fingerprints (gates D, F: PASS)

* **D24–D64**: two dispatches of the same device state are bit-identical
  on every output buffer — forward partials (1,296 / 3,072 / 10,368 /
  24,576 f64) and gradient (55,296 / 131,072 / 442,368 / 1,048,576 f64)
  compared bitwise.  This is the payoff of the tile-canonical reduction
  (fixed lane order in-workgroup, fixed workgroup order in the CPU
  combine) and race-free gather.
* **F24–F64**: FNV-1a-64 fingerprints over the LE bytes, the
  within-backend goldens:

| N  | forward partials     | gradient             |
|----|----------------------|----------------------|
| 24 | `3f75d11ca2d4a539`   | `abd80c50faa4e9c2`   |
| 32 | `bae2844c642b0785`   | `a8002235e1de3caf`   |
| 48 | `3381707ad6d6b68d`   | `7d9c07e1f904f901`   |
| 64 | `67f84dcfd8469265`   | `c69b6954ae0e791d`   |

* **Cross-run reproduction (measured, this container)**: a full second
  run of the gates binary passed all 20 gates, reproduced every
  fingerprint above, and produced a `k3_results.json` identical in every
  physics field (compared with timing fields stripped; only timings
  differed).  Fingerprints are goldens **for this stack** (llvmpipe /
  Mesa 25.2.8 / LLVM 20.1.2): a driver or Mesa upgrade may legitimately
  move the bits — then D (own-replay) must still hold and S/G/DD must
  still land in band; the fingerprints are re-baselined, not "fixed".

## 6. Arrested-descent smoke with the GPU gradient (gates AD-a, AD-b: PASS)

N = 24, seed = hedgehog + bump (seed 424242, amp 0.02), guards referenced
to the N = 24 hedgehog (`with_guards(deg_h, floor_gap_h − FLOOR_BAND)`),
40 iterations, dt0 = 0.01, dt_max = 0.05.  The descent loop steps with
the **GPU gradient** and arrests on the **GPU objective**; every ACCEPTED
state is downloaded and re-evaluated by the CPU `eval`:

* 38 accepted steps, 2 arrests, no stall;
* **CPU-evaluated objective monotone non-increasing over all 39 accepted
  states** (zero slack): 3.781355477 → 3.235577496 (net −5.458e-1);
* GPU objective also monotone; worst GPU-vs-CPU objective relative
  difference over the accepted trajectory 1.16e-13.

## 7. Wall-time honesty (record, NOT a performance claim)

llvmpipe is **software Vulkan executing on the same 4 CPU cores** as the
CPU reference — these timings demonstrate correctness/readiness cost
only; no GPU speed claim is made or implied (a software rasterizer losing
to a native sequential loop is the expected sign, and it does):

| N  | CPU measure | CPU eval+grad | GPU init | GPU forward+readback | GPU flux+gather+readback |
|----|-------------|---------------|----------|----------------------|--------------------------|
| 24 | 0.001 s     | 0.003 s       | 8.60 s¹  | 0.231 s              | 0.141 s                  |
| 32 | 0.002 s     | 0.007 s       | 0.137 s  | 0.037 s              | 0.107 s                  |
| 48 | 0.007 s     | 0.024 s       | 0.068 s  | 0.114 s              | 0.533 s                  |
| 64 | 0.015 s     | 0.061 s       | 0.063 s  | 0.226 s              | 0.744 s                  |

¹ first-engine cold start: llvmpipe/LLVM pipeline compilation of the 8
f64 kernels; subsequent engines 0.06–0.14 s, and Mesa's on-disk shader
cache makes a repeat run of the whole binary finish in 1.5 s (run 2:
N = 24 init 0.078 s, forward 0.003 s).  AD smoke (40 GPU iterations +
39 CPU re-evals + downloads, N = 24): 0.5 s.  Whole battery: 12.5 s
cold, 1.5 s warm.

## 8. Real-GPU migration path (extends RESULTS.md §5)

Nothing in the certified math changes on hardware — only
`Instance::request_adapter` returns a different adapter.  Concretely:

* **SHADER_F64 availability**: required at device creation (hard E1
  fail otherwise).  Vulkan discrete/integrated GPUs with
  `shaderFloat64` qualify (most NVIDIA/AMD discrete parts; Intel
  integrated often NOT); Metal and WebGPU-in-browser are out — no f64.
  Consumer parts run f64 at 1:32–1:64 ALU ratio: correctness inherits,
  throughput must be re-measured there, never extrapolated from here.
* **Workgroup sizing**: WG = 64 is part of the bit contract (partials
  layout + lane order).  Retuning WG for occupancy on hardware changes
  the partials buffer and the within-backend bits — allowed, but it
  re-baselines the F fingerprints and is a certified-surface change
  (re-run this suite).  The CPU-side fixed-order combine is unaffected.
* **Tolerance bands**: real drivers FMA-contract/schedule differently →
  different bits *within the same bands* (S: 1e-12; G: 1e-11 RMS /
  1e-10 max-abs; DD: 1e-5 / 1e-9 cross).  The bands were sized from
  measured llvmpipe headroom (×9 worst-case at N = 64 sectors, ×1,000
  gradient) precisely so a different contraction pattern fits; a
  hardware run exceeding them is a finding, not a tolerance to widen.
* **Own-replay rule**: D gates (bit-identical replay on the SAME device)
  and the AD/K4 decision-equality class must hold exactly on hardware;
  cross-device bit identity is explicitly NOT promised.
* **Scale ceilings** (wgpu default 128 MiB/binding): field buffers to
  N ≈ 251, gradient buffer to N = 161, flux intermediate to N = 101;
  past that raise device limits or split the flux buffer per axis — a
  mechanical change outside the certified math.
* First run on metal: re-run `gum_gpu_k3_gates` + `gum_gpu_gates`
  unchanged; record the new adapter line and fingerprints.

## 9. Defects / limitations (printed as measured)

1. Cold-start pipeline compile is 8.6 s on llvmpipe (first engine only) —
   recorded in §7; harmless here, but a real-GPU CI budget should expect
   a shader-compile warmup.
2. GPU-vs-CPU timing is uniformly a loss on this rig (expected: software
   Vulkan vs native loop on the same 4 cores) — stated in §7; no
   performance readiness is claimed, only correctness readiness.
3. The fingerprints (§5) are stack-goldens, not device-class goldens; a
   Mesa/LLVM upgrade moves them legitimately (re-baseline protocol
   stated).
4. AD-a gates monotonicity with zero slack across backends; it passed
   here (accepted decreases ≫ the 1e-13 cross-backend objective noise),
   but a descent tuned to near-tie acceptances could flip a comparison
   at the noise floor — the gate would then correctly surface it as a
   finding to adjudicate, not silently absorb it.
5. Guard-dressed gradient weights are exercised in the AD smoke
   (guards referenced to the hedgehog) and were band-gated directly in
   G2 K3b at N = 48; K3 adds no separate per-N guarded-gradient gate.
6. N = 64 is the largest grid in the battery (12.5 s total wall);
   nothing was truncated for time.

## 10. Epistemic notice (binding)

Everything here is within-model numerical engineering on a speculative
theory's functional, on a software Vulkan rig.  A passing gate certifies
the port, the adjoint, the deterministic reduction design, and the
backend plumbing — never anything about nature.  llvmpipe green ≠
hardware certified: the first run on metal must re-pass this suite (same
bands, own-replay bit rule) before any real-GPU number is used as
evidence.
