# H3.2 — T-H5 adjudicated by computation: what honest deformation from the rigid rung actually selects

**Phase:** H3.2 (ROADMAP_v5_THEORY Phase H3; numeric adjudication of T-H5 =
`T2_closure_pillars.md` §d.3 DEFECT-CANDIDATE).  **Date:** 2026-07-16.
**Engines:** `fs-gum-field` (radial solver, Field3), `fs-gum-statics`
(guarded ANF protocol, halo referee), `fs-gum-kern` (deterministic tiled
threaded layer) — none modified.
**Solver:** `h32_solve/` (this directory), bin `h32_selection`
(usage: `h32_selection <outdir> [m_desc] [ncyc] [threads]`, run at 25 400 4);
analysis layer `h32_analyze.py` (appends the `analysis` block to the JSON).
**Compute:** 577 s (3 flows × 10,000 guarded descent iterations, 4 threads).
**Output:** `h32_results.json` (full per-cycle trajectories + analysis).
**Determinism:** in-run bitwise replay of free cycles 1–3 asserted PASS;
cross-execution + cross-thread-count (4 vs 2) reruns byte-identical
(modulo the recorded `threads` arg and runtime), prefix-identical to the
production run.  See §6 for the honest cross-host boundary.

**Epistemic frame (binding).**  Within-model numerical engineering on a
speculative theory's functional.  Everything below is a measured fact about
the campaign's certified functional and the corpus's argument structure —
never about nature.  These are adjudication-relevant facts; the
adjudication itself is the coordinator's.

---

## 1. The claim under test, and the computation that tests it

Corpus IV.D: *"the rigid spinning knot is superradiant and is not a
solution.  Self-consistency forces centrifugal deformation **until κ < 1**:
the in-gap character of matter is a derived necessity"* — with the
convergence endpoint exhibited as the ⟨r1⟩ benchmark (κ = 0.802 ± 0.018 at
ε = 0.05).  T2 §d graded the superradiant **rejection** SOUND
(convention-robust) and the **selection narrative** DEFECT-CANDIDATE
(asserted, not derived; criterion applied asymmetrically).  T3.1 predicts
the honest flow passes κ = 1 without stopping and heads for the saturation
locus κ_paper = 1/√2.  H2.6 measured the descent-relaxed clock *root*
drifting to κ_paper ≈ 0.70 with budget.  What was still missing: the
narrative executed **as a dynamics** — deformation from the rigid rung
itself, with the corpus's own self-consistency in the loop.

**The flow (quasi-static selection dynamics).**  On the frozen N = 48 gate
protocol (LBOX = 4.5, ε = 0.05 ⇔ t = 0.0082764349, corner scheme, frozen
guards MU = 5000 at deg_ref = 0.945803 / MU_F = 400 at fgap_ref = −0.125856;
static anchor = hedgehog + bump(424242, K = 6, amp 0.02) + 900 guarded ANF
iterations ⇒ I_s = 22.60711, E_s = 2.9865333), alternate:

- **(a)** 25 guarded ANF descent iterations of the fixed-L Routhian
  R(q; L) = E_static + L²/2I — *energy-decreasing deformation at fixed
  charge and spin*, monotone non-increasing within each block by
  construction;
- **(b)** the corpus's own clock self-consistency, L ← √((2/3)·I·E_static)
  (Section-K form of E_tot = 𝔠ω at j = ½; forces E_rot/E = ¼ identically at
  each update),

for 400 cycles = 10,000 descent iterations per run.  **Start = the rigid
rung**: the relaxed spherical static solution spun at the rigid-rung clock
rate, L_rigid = √(7/6)·I_s/(2√π) = 6.888323 via the validated unit map
κ_paper = 2√π·L/I — the seed measures κ_paper = 1.08012345 = √(7/6) exactly
by construction.  Instrumented per cycle: κ_paper (block and post-update),
E_rot/E, halo fraction, I-radial profile (fractions beyond 1.0/1.25/1.5 R*,
I-rms radius), degree, bare-Routhian descent rate dR/diter, guard penalties,
arrests.

**Controls (the F-R5-class rescue test):** the identical flow with far-field
growth re-seeded after every descent block — every cell beyond r_c of the
current (1−q0) centroid reset to the static-anchor field, then renormalised
(a deliberately crude in-gap projection), at two radii: **1.5 R*** (the
campaign's halo diagnostic radius — far-field suppression as F-R5's referee
defines it) and **1.25 R*** (compact support just above the corpus's own
in-family support, V^{1/3}R* = 1.12 R* at V = √2; the tightest reading under
which the corpus's claimed deformations remain admissible).  These test
whether the narrative is rescuable by an unstated constraint.

## 2. The free flow: the measured κ(t) trajectory

| level crossed (downward) | meaning | cycle (≈ iteration) | dκ/dcyc before / at / after | plateau? (min/med rate near crossing) |
|---|---|---|---|---|
| **1.000** | the corpus's "until κ < 1"; its own measured onset 1.000 ± 0.004 | **16.3 (~410)** | −3.25 / −2.82 / −2.45 ×10⁻³ | **none** (0.97) |
| **0.802** | the ⟨r1⟩ benchmark = the narrative's exhibited endpoint | **130.3 (~3,260)** | −1.18 / −1.09 / −1.02 ×10⁻³ | **none** (0.95) |
| **0.7638** | √(7/12), the corpus's deep-BPS κ_phys | **169.3 (~4,230)** | −0.93 / −0.87 / −0.81 ×10⁻³ | **none** (0.96) |
| **0.7071** | 1/√2, the F-R5/T3.1 saturation locus | **249.2 (~6,230)** | −0.61 / −0.57 / −0.54 ×10⁻³ | **none** (0.96) |

**Nothing in the dynamics distinguishes κ = 1.**  At the crossing:
dR/diter = −2.23×10⁻⁴ (smoothly interpolating the monotone −2.9×10⁻⁴ →
−0.4×10⁻⁴ decay of the whole run), halo fraction 0.00066 growing at a
negligible 5.7×10⁻⁶/cycle, d ln I/dcyc = 0.57 %, zero arrests, degree
guarded flat (0.9402).  The local |dκ/dcyc| minimum near the crossing is
97 % of the neighbourhood median — the rate curve has no feature there at
the 3 % level.  The corpus narrative's entire "deformation until κ < 1"
phase is **the first ~4 % of the honest flow**; the remaining 96 % of the
budget is spent continuing the same descent *in-gap*, through the
benchmark, through the deep-BPS value, through the saturation locus.

**Endpoint (10,000 iterations):** κ_paper = 0.6443, still descending at
−3.0×10⁻⁴/cycle (dR/diter = −3.6×10⁻⁵ < 0, zero arrests) — no stationary
point reached at budget.  I/I_s = 2.73, E_static +2.4 %, 𝔠_paper = 2.85.
Rate-vs-κ fit dκ/dcyc = −λ(κ − κ∞): tail half gives λ = 4.5×10⁻³/cycle,
**κ∞ = 0.580** (rms 7.5×10⁻⁶; full-trajectory fit 0.623) — the attractor
estimate sits at/below the saturation locus, 0.13–0.22 *below* the
benchmark.  The raw endpoint alone is already 8.8σ below the benchmark in
the corpus's own error bar.

**Where the growth lives (measured):** the I-fraction beyond 1.5 R* stays
≤ 0.007 the whole run — at this box the F-R5 *far*-halo channel is barely
populated; the inertia growth is **near-field spreading** (fraction beyond
1.0 R*: 0.13 → 0.40; I-rms radius 0.79 → 0.95 R*).  The descent finds its
in-gap deformation dollars just outside the compacton support, inside the
halo diagnostic radius.

**The clock never resists.**  E_rot/E = ¼ exactly at every update
(identity), block-tracking lag ≤ 8.6×10⁻³ — the corpus's self-consistency
condition is satisfied continuously *along the whole descent* and supplies
no stabilisation whatsoever (the H2.6 clock-blindness fact, now measured at
flow level).

**Coincidence-class observation (flagged, not charged):** the flow's
(κ, 𝔠_paper) transit passes (0.802, 2.252) at cycle 130 and
(0.764, 2.372) at cycle 169 — the ⟨r1⟩ benchmark tuple (0.802 ± 0.018,
2.37 ± 0.09) is consistent, within its own bars, with a *snapshot of this
transit* (F-R5's budget-stamped-saddle reading), not with a terminus of it.

## 3. The controls: is the narrative rescuable by an unstated constraint?

| run | projection bite (tail dI cut/cycle) | κ = 1 | 0.802 | 0.7071 | endpoint (cyc 400) | tail drift /cyc | tail-fit κ∞ |
|---|---|---|---|---|---|---|---|
| free | — | 16.3 | 130.3 | 249.2 | 0.6443 | −3.0×10⁻⁴ | 0.580 |
| control 1.5 R* (far-field/F-R5 radius) | 2.8×10⁻⁴ (inert) | 16.3 | 130.8 | 251.5 | 0.6470 | −2.9×10⁻⁴ | 0.587 |
| control 1.25 R* (compact support) | 2.3×10⁻³ (binds; 153 arrests) | 16.8 | **142.3, no pause** | 365.6 | 0.7008 | −1.75×10⁻⁴ | **0.679** |

- **Far-field suppression (1.5 R*) is inert**: the projection cuts almost
  nothing because the growth is near-field (§2); the trajectory is the free
  flow to within 0.4 %.  The narrative is NOT rescued by suppressing what
  F-R5's referee calls the halo.
- **Compact support (1.25 R*) binds hard and still fails**: it cuts
  inertia every cycle (the descent re-grows mass against the cut — the
  "stationarity" it eventually approaches is constraint-enforced, with
  dR/diter still negative at −7×10⁻⁶), *slows* the flow (saturation-locus
  crossing delayed 249 → 366), but does **not stop it at the benchmark**:
  0.802 is crossed at −9.4×10⁻⁴/cycle with no plateau (rate ratio 0.95),
  and the constrained flow converges toward **κ∞ ≈ 0.68 — the
  saturation-class locus (1/√2 − 4 %), not the benchmark**.  Its tail mean
  0.7025 sits 5.5σ below 0.802 on the corpus's own error bar and is still
  drifting down.

So even granting the corpus an unstated constraint of the F-R5 lock class —
at either the referee's radius or the tightest compact-support reading —
the "stop in-gap at the benchmark" clause stays false; the constrained
flow's destination is the same saturation-class point T3.1 derives.

## 4. Connection to T3.1 and the H2.6 drift

Theorem T3.1 (C2 ≡ C3): the stability-constrained closure coincides exactly
with the saturated closure — constrained optimum ON the threshold boundary,
KKT-marginal, **no intermediate stabilized tuple between the rung and
saturation**.  The flow measurements realize this dynamically: no run
decelerates toward *anything* in (0.71, 1.08); the only attractor-like
behaviour anywhere is the saturation-class region (tail κ∞ = 0.58–0.68
across runs, with the constrained flow converging to 1/√2 − 4 % and h26's
budget-matched cold-start clock root at 0.6997 ≈ 1/√2 − 1 %).  The free
lattice flow *overshoots* 1/√2 slowly through the near-field channel
(measured, §2): the exact marginality at 1/√2 is the idealized *far*-halo
ledger's statement, and the ε = 0.05 / N = 48 near-field channel keeps
descending below it — consistent with F-R5's "no budget-independent root",
and a-fortiori fatal to "stops in-gap": the honest infimum direction lies
at or *below* saturation, never between the onset and the benchmark.

## 5. Adjudication-relevant summary (facts; coordinator adjudicates)

1. **The rejection half of IV.D survives untouched** (T2 §d.2 SOUND): the
   rigid-rung seed κ = √(7/6) is above every operative threshold and the
   honest flow leaves it immediately and permanently.
2. **"Deformation until κ < 1" is true only as the first ~4 % of the
   flow; "…and stops in-gap" has no honest reading that survives.**  The
   κ = 1 crossing carries no dynamical signature at the 3 % level in any
   instrumented observable (rate, dR/diter, halo, I-growth, arrests,
   degree); descent continues through the entire in-gap band.
3. **The benchmark κ ≈ 0.802 is not selected by the corpus's own
   dynamics**: crossed without pause at cycle 130/400; free endpoint 8.8σ
   below it; even the tightest constrained flow ends 5.5σ below it and
   still descending.  Within-model, the ⟨r1⟩ tuple is consistent with a
   budget-stamped snapshot of the transit (§2, flagged coincidence-class).
4. **The selection endpoint is the saturation-class locus, per T3.1**:
   every attractor estimate lands at/below 1/√2 (0.58–0.68 fitted; 0.6997
   at the h26-matched budget class), never in the corpus's claimed regime.
5. **The F-R5-class rescue fails**: far-field suppression at the halo
   radius is inert (the growth is near-field — itself a new measured
   fact); compact-support suppression slows the flow but converges it to
   ~1/√2 − 4 %, not to the benchmark.
6. **Self-consistency has no stabilizing content**: the clock identity
   E_rot/E = ¼ holds along the entire descent (lag ≤ 0.9 %); it re-arms
   the descent (L grows with I) rather than arresting it.

**Net for T-H5:** the DEFECT-CANDIDATE grade of the *convergence/selection*
narrative is confirmed by direct dynamics on the corpus's own functional
and own clock, with the discharge path (unstated-constraint rescue) also
measured and closed for the two natural constraint readings.  The repair
that survives is exactly T3.1's: the narrative's honest endpoint is the
saturated closure, not the in-gap benchmark.

## 6. Determinism and honest boundaries

- **Determinism:** frozen LCG stream, fs-gum-kern bit contract.  In-run
  bitwise replay (free cycles 1–3) asserted PASS inside the production
  binary.  Two fresh executions at 4 and 2 threads: JSON byte-identical
  except the recorded `threads` arg and runtime; both prefix-identical to
  the production trajectories.  Same binary + args + host ⇒ bit-identical
  JSON.
- **Cross-host boundary (disclosed):** a container restart re-executed the
  production run mid-flight; the interrupted instance (different container
  placement) agreed with the recorded one to all printed digits through
  ~cycle 300 and diverged at ulp-amplification level in the tails (free
  endpoint 0.665 vs 0.644; compact-support 0.755 vs 0.700) — the
  arrest-cascade ulp amplification the fs-gum-kern docs themselves flag for
  cross-backend comparison.  Every crossing, plateau-absence, and endpoint
  fact above holds identically in both instances; the recorded artifact is
  the single complete execution (replay PASS).
- **Budget:** the free flow is not stationary at 10,000 iterations; κ∞
  values are fits (extrapolation), all level-crossing and no-plateau facts
  are raw measurement.  The below-1/√2 drift is a lattice/near-field-channel
  fact at ε = 0.05, N = 48, LBOX = 4.5 — not a claim about the continuum
  infimum.
- **Protocol dependence:** one quasi-static granularity (25-iteration
  blocks) was run; the block size only sets the clock-tracking lag
  (≤ 0.9 %).  The independent H2.6 protocol (cold-start roots, 350/700
  budgets) gives the same direction and the same endpoint class.
- The crude projection is a *sufficient* test of the two natural constraint
  readings, not of every conceivable lock; a corpus discharge would have to
  exhibit its constraint explicitly (axi3 §5 route (a)), at which point it
  can be run through this same harness.

## Files

- `h32_solve/Cargo.toml`, `h32_solve/src/bin/h32_selection.rs` — the
  selection-flow binary (free + two controls + replay; deterministic).
- `h32_results.json` — protocol, rigid-rung seed, full per-cycle
  trajectories of all three runs, references, determinism flags, and the
  appended `analysis` block (crossings, window stats, plateau tests, rate
  fits, verdict facts).
- `h32_analyze.py` — the deterministic analysis layer (recomputes and
  re-appends `analysis` from the raw trajectories; idempotent).
- `h32_RESULTS.md` — this memo.
