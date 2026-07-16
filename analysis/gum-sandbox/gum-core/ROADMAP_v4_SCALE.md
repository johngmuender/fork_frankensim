# GUM Physics Core — Roadmap v4: Scale (CPU → GPU)

Successor to ROADMAP_v3 (executed in full). Closes the scale gap of
GAP_ANALYSIS_v4_SCALE.md. Design rules carried forward: machine-checked
gates against frozen referees; deterministic-by-construction wherever
possible; deliberate, documented golden bumps when a canonical order
changes; tolerance-band/metric goldens only at backend boundaries (GPU,
and ANF endpoints across backends); standalone crates under gum-core/.

## Phase G1 — `fs-gum-kern`: the deterministic parallel kernel layer ▶ EXECUTING
The one parallel primitive: a domain-decomposed sweep (stencil-read →
per-point physics → per-point flux records) with (a) per-tile partial
reductions in a thread-count-INDEPENDENT tile decomposition (i-slabs,
TILE_I = 4, part of the bit contract) combined by a fixed pairwise tree,
and (b) the adjoint inverted from scatter to GATHER. Backends: serial
(the new canonical order — one deliberate golden bump, documented and
re-frozen) and threaded (std::thread::scope, tile dispenser; bit-identical
to serial BY CONSTRUCTION — gate asserts hash equality at 1/2/3/4
threads). Vendors fs-exec's pairwise_fold/Compensated (~200 LOC,
swap-back documented — the B3/F3 precedent). Gates: equivalence to the
old engine at machine-eps-class once, then re-frozen bit goldens;
measured speedup table (target ≥ 3× at N=96 on 4 threads); fs-gum-statics
descent gates re-run through the new layer.

## Phase G2 — wgpu f64 backend on lavapipe ▶ EXECUTING
The GPU rung, TESTED HERE via software Vulkan: WGSL f64 compute kernels
(SHADER_F64 / naga FLOAT64) for the forward sums + flux + gather sweep,
device-resident field, per-component storage buffers, fixed-workgroup
reduction with a deterministic combine on readback. Goldens: the
tolerance-band/metric class at kernel level (sector sums vs CPU within
stated bands) and at ANF endpoints (deg/E_static/κ/halo bands) — the
sanctioned fallback, stated per gate. Lavapipe throughput is irrelevant
(it is a correctness rig); the deliverable is a working, gated backend
that a real GPU inherits unchanged.

## Phase G3 — the G*-vs-𝔠₀ precision push ▶ EXECUTING
Single-core, 2-D spectral, sub-hour: drive the halo-saturated objective
G* = min(E_static − I/16π) to ≤ 1e-7-class precision (basis/grid
ladders, Richardson) and adjudicate whether G* = 𝔠₀ = 128√42/105π
exactly (currently 1.2×10⁻⁴ apart, drifting down). Either outcome is a
result: a new exact identity for the F-R5 record, or a pinned gap.

## Phase G4 — the I.2 two-knot bond-equation flagship ▶ EXECUTING
The corpus protocol never executed anywhere: two B=1 knots at
separations, product-ansatz seed, guarded relaxation on 144×96×96
(h = 0.09375 frozen), 3 relative orientations; extract the interaction
energy vs separation, fit the screened-dipole ℬ_eff form, close the
bond equation and compare x₀ = 1.90 ± 0.05 (𝔟 = 42) vs the corpus's
measured 1.92 ± 0.08. Runs as a 4-process task farm on the EXISTING
serial engine (measured 3.96× multiprocess scaling) — per-run
bit-identity intact; the fitted x₀ is a tolerance-band observable by
nature. Reduced protocol if budget forces it (fewer separations,
documented), full protocol as the follow-on.

## Phase G5 — the N=192 F-R5 endpoint (after G1)
The threaded kernel's proof workload: the halo-descent referee at 8×
the campaign's cells, refining the F-R5 endpoint and hardening the
DISCHARGE_PACKAGE. Scheduled once G1's speedup lands.

## Out of envelope (unchanged)
Blue-fog condensate/disclination boxes and tower-depth studies at scale
(real-GPU class — the G2 backend is their on-ramp); gravitation;
nucleation.

## Execution record
- Phase G-A (measured scaling survey, 4 surveyors): ✅ docs + JSON.
- **G1 ✅** fs-gum-kern: 13/13 gates; serial↔1/2/3/4-thread BIT-IDENTICAL
  (one BLAKE3 hash); 3.61× at N=96; golden bump 2.8e-14 documented;
  corner adjoint inversion cost zero bits; SIMD honestly unclaimed.
- **G2 ✅** fs-gum-gpu: f64 WGSL on llvmpipe; all gates green; identical
  descent decisions vs CPU; tolerance-band doctrine confirmed by
  measured 1–2 ulp FMA contraction.
- **G3 ✅** G* = 16√2/9 EXACTLY (oblate compacton closed form); the 𝔠₀
  identity disproven (24√21 ≈ 35π accident); saturated closure exact:
  𝔠_paper = 64√2/(9π) at κ = 1/√2.
- **G4 ✅** fs-gum-twoknot: App I.2 executed for the first time anywhere
  — 23 guarded relaxations, 160×96×96, 3 orientations × 7 separations +
  far anchor, 49.5 min wall (3.7× farm scaling); anisotropic engine
  bitwise-gated (20/20) against fs-gum-statics. **Bond loop lands in
  band**: 𝔟 = 42 closed-loop x₀ = 1.88/1.97 vs predicted 1.90 ± 0.05;
  direct x₀ = 2.0 ± 0.15 vs corpus-measured 1.92 ± 0.08 (0.5σ);
  attractive well at x = 1.89 (2.1σ); channel sign structure confirmed;
  tail mass ±3%. NOT recoverable: the 𝔟_eff absolute normalization
  (~7 orders — App-A convention at ε = 0.05 with overlapping cores;
  same spec-underdetermination class as the Tier-2b 𝔭 amplitude).
  Caveats: iter_cap runs, 2.1σ well depth, product-ansatz bias — the
  longer-cap rerun is the natural G5 companion workload.
- **G5a ✅** N=192 F-R5 referee on fs-gum-kern (fs-gum-kern/
  n192_RESULTS.md): state-space path grid-robust; clock ladder resolved
  (state function, floor ≈ 1.245–1.2525, ratio > 1.24 on every B=1
  state); control cleaner at N=192; 2T/4T bit-identity at 192³.
  Corrections printed: iteration milestones are solver-clock; the
  fixed-h box axis is a NULL (4A's big-box plank was h-confounded);
  the N=192 crossing was subsequently OBSERVED directly at
  it 2070–2080 (G5c, 46.3 min quiet run): kappa through 1/sqrt(8pi) with
  I at the crossing 29.7–29.8 at EVERY grid in the ladder — the F-R5
  referee is now fully observational at h = 0.046875.
- **G5b ✅** two-knot refinement (fs-gum-twoknot/twoknot_REFINE.md):
  cap systematic measured (G4's well depth was under-converged 22× its
  σ; location right); well now **5.0σ** (Aitken, −6.53e-3 ± 1.3e-3,
  two-point feature); **x₀ = 1.95 ± 0.07 vs corpus 1.92 ± 0.08
  (0.3σ)**; 𝔟 = 42 crossing honest-negative at 1.81 (1.8σ);
  𝔟_eff normalization still unrecoverable (𝔭-class).
  DISCHARGE_PACKAGE updated with the G5 refinements + corrections.

## Program status: ROADMAP v4 EXECUTED (G-A, G1–G4, G5a–G5b)
Next (unscheduled): explicit
k-contiguous SIMD (the documented unclaimed step); real-GPU execution
of fs-gum-gpu; blue-fog boxes on that on-ramp; the corpus's response.
