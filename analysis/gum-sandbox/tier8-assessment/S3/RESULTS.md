# S3 — O-1 pilot: validated forward flow-map covering (F-T8-S3)

**Goal.** Pilot the extension of the Q1 validated-enclosure technology from
single backward paths to a FORWARD box covering — the only identified
technology for obligation O-1 of Conjecture C (Q2/ANALYTIC_RECON.md §5:
early-time flow control on [0, T0] × supp ρ₀, "the dominant obligation").
A pilot is NOT an O-1 discharge; the deliverable is the feasibility datum:
fraction of a pre-registered subdomain certified crossed-or-localized by T0,
per-box validated-step cost, width/wrapping growth, and the honest
extrapolated cost of a full O-1 covering.

**Campaign:** Tier 8 Phase S (ROADMAP_v10_ALTERNATIVES.md, S3).
**Date:** 2026-07-18. **Code:** `s3_cover.py` (phases g1 / bench / cover /
ladder / truegrowth / gates / fig; per-box checkpoints `_s3_box_*.json`,
resumable `_s3_ck_*.json`; production under nohup, logs `s3_produce.log`,
`s3_truegrowth.log`). **Raw output:** `s3_results.json`, `_s3_g1.json`,
`_s3_bench.json`, `_s3_truegrowth.json`. **Figure:** `s3_fig.png`.
**Prior art (mandatory reads, honored):** tier7-program/Q1/q1_enclose.py +
RESULTS.md (the interval engine, REUSED — all interval semantics imported
unchanged), Q2/ANALYTIC_RECON.md §5 (O-1), tier6-foundations/L5′ RESULTS
(testbed; τ_max(n_y = 1) = 5.130950).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**S3-G1, S3-G2, S3-G4 PASS; S3-G3 PASS with a hard finding:** the forward
engine is VALIDATED (the forward replay of Q1 path 134's certified backward
tube re-encloses the Q1 anchor: containment TRUE, end widths 7.3e-6, gate
G1); the pre-registered 24-box covering plus one 2×2 subdivision round was
executed to terminal status with zero unaccounted subdomain (96/96 leaves,
leaf areas sum exactly to the band area) — and the honest result is
**0 % certified: every leaf box terminated in the blow-up class before
t ≈ 1.05** (95 f-disc-zero, 1 Picard-inflation z-exit), a diagnostic width
ladder found NO initial width down to 3.1e-4 that survives to T0 = 5.14
(deaths at t = 1.05 → 1.95), and the float-layer variational diagnostic
shows why: the TRUE flow-map amplification over [1, 5.14] is only 17–279
(median 127), while the enclosure growth exceeds it by a WRAPPING factor
that compounds from 2.2× to 55× before t = 2 — the pilot engine's
first-order wrapping, not the flow, is the bottleneck. Extrapolated full
O-1 cost with THIS engine: 1.2e13–9.4e16 boxes per |n_y| family
(3.7e11–2.8e15 core-hours) — **infeasible**; the wrapping-free adaptive
lower bound from the true amplification is ~3.5e6 boxes (~1.0e5
core-hours/family) — expensive but conceivable. **The pilot's feasibility
datum: O-1 is not discharged by scaling this engine up; it needs
higher-order (Taylor-model/CAPD-grade) enclosures plus adaptive
subdivision.** The kill-bin theorem remains open; no register moves.

## 1. Setup, subdomain, and engine adaptation

**Certified object (identical to O1/P2/Q1).** The exact Nx = 2048
band-limited field with the analytic z-factor; guidance velocity
v_x = Im(f′/f) − n_y π cot(πz), v_z = n_y Re(f′/f); n_y = +1. All interval
semantics — outward-rounded float64 scalar intervals, complex disc
arithmetic with Gargantini–Henrici formulas, mpmath.iv prec-120
transcendentals, binary-exponentiation mode phases, Lohner QR with rigorous
2×2 adjugate inverse, verified Picard a-priori boxes — are **imported
unchanged from q1_enclose.py**; the trust base is Q1's, verbatim.

**Forward adaptation (printed).** (i) Time direction: forward step
s0 → s1 = s0 + h via the exact identity X(s1) = X0 + h v(X0,s0) +
∫ (s1−u) D du, (s1−u) ∈ [0,h], so X(s1) ∈ X0 + h v + (h²/2) hull(D(B,I_s)),
with S = I + h J(B,I_s) in the mean-value form (Q1: mirrored signs); the
h²/2 factor is applied as two outward-rounded scalings (0.5·h exact), no
pre-rounded h² constant. (ii) Fat initial data: the Lohner frame starts at
A = I with q = box − center an interval RECTANGLE (Q1: a point).
(iii) Picard pad is h-scaled (0.35 w(hv) + 4e-3 h + 1e-14, ×2 inflation on
failure, ≤ 8 attempts) instead of Q1's tube-width-scaled pad, which for fat
boxes would inflate B by ~55 % of the box size; the pad is heuristic, the
Picard containment CHECK — the rigorous part — is unchanged. Warm-started
velocity guess across steps. Final partial step lands exactly on fl(5.14)
(Sterbenz-exact hh; outward time intervals).

**Subdomain (pre-registered) and coordinate clarification (spec item).**
The Q2 T_ball band is w ∈ [−1.2, −0.6], z ∈ [0.05, 0.95] with
w = (x − x0)/t − k0 (x0 = −5, k0 = 2). On the single slice t0 = 1.0 this
maps EXACTLY to an (x, z) rectangle: x = x0 + (k0 + w)t0 = w − 3, i.e.
**x ∈ [−4.2, −3.6]** — no awkwardness, rectangle to rectangle. Grid:
6 × 4 = 24 boxes of 0.1 × 0.225. Horizon: t ∈ [1.0, 5.14]
(T0 = 5.14 ≳ measured T_max(1) = 5.1310). Detector d = 1.

**Classification (pre-registered + printed clarifications).**
(i) CROSSED: lower(X_x) > d at some s ≤ T0; (ii) LOCALIZED:
upper(X_x) < d at T0, enclosure printed; (iii) blow-up class: width > 1
(pre-registered cap) — **plus, grouped as blow-up-equivalents
(clarification):** FAIL_FDISC (the f-disc reaches 0 so v is no longer
enclosable — in practice this fires FIRST, at widths 0.23–0.48, i.e. it is
a strictly earlier terminal of the same divergence; the width-1 cap was
never reached by any box), FAIL_WALL (a z-interval leaves (0,1)),
FAIL_PICARD (no verified a-priori box after 8 inflations). A box reaching
T0 with width ≤ 1 straddling d would be UNRESOLVED, counted with FAIL
(clarification; the case never occurred). Failure-class parents subdivided
once 2×2, children re-run, then stop (as specced).

## 2. S3-G1 — engine validation by forward replay (PASS)

Q1 path 134 (n_y = +1, t0 = 6.545, z0 = 0.4125, h = 2⁻¹⁴, 1359 steps,
CERTIFIED) ended at s_stop = 6.46205322265625 with enclosure
x ∈ [1.001003452765, 1.001006184883], z ∈ [0.399170402049, 0.399173212674]
(widths 2.7e-6/2.8e-6). The true backward trajectory through the anchor
(d, z0) = (1.0, 0.4125) lies in that box, so a valid forward enclosure of
the box over the same 1359 time steps (bitwise-identical time grid
s = t0 − (K−k)h) MUST contain the anchor. Measured: end box
x ∈ [0.99999634, 1.00000366], z ∈ [0.41249793, 0.41250207] — **anchor
contained** (widths 7.32e-6 / 4.15e-6, growth 2.6× over depth 0.083, 7 s).
Containment, not distance, as gated.

## 3. Benchmark and production protocol

Bench (box p31, 0.1 × 0.225, h = 2⁻¹²): 5.0 ms/validated step (two field
packs per step: thin nmax=1 + fat nmax=3 over the Picard box), full-horizon
worst case 16,958 steps ≈ 86 s/box. **h = 2⁻¹² fixed before production**
(growth is dominated by J-interval widths, which are h-independent per unit
time; halving h doubles cost with no width benefit — the h² local term is
~1e-8, negligible against box widths). No scale-down of the pre-registered
M = 24 was needed; the whole production (g1 + bench + 24 parents + 96
children + 8 ladder rungs + gates + figure) took 92 s wall plus 203 s for
the float diagnostic, far under budget. All runs under nohup with atomic
per-box JSON checkpoints (`_s3_ck_*` every 4000 steps; every terminal state
in `_s3_box_*.json`); no restart was needed.

## 4. S3-G2 — the covering, executed honestly (PASS; 0 % certified)

All 24 parents terminated FAIL_FDISC within 4–46 steps (t ≤ 1.0112),
enclosure widths having grown from 0.1/0.225 to 0.23–0.48. All 24 were
subdivided (one 2×2 round, as specced) → 96 children (0.05 × 0.1125):
**95 FAIL_FDISC + 1 FAIL_WALL**, deaths t ∈ [1.0017, 1.0486], median 133
validated steps, median 0.7 s/box, total covering wall cost 51 s.
Growth factors at termination (cutoff-limited by the f-disc death, not a
free growth measurement): min 1.92, median 2.33, max 4.68. The one
FAIL_WALL (child p22_c00, final z ∈ [0.46, 0.69]) is a Picard-INFLATION
z-exit at an already-blown box — the ×2⁷ inflated a-priori pad pushed the
B-box z-interval outside (0,1) — not wall proximity; blow-up-equivalent
either way.

**Accounting: zero unaccounted subdomain.** Every parent is either terminal
or replaced by its 4 children; 96/96 leaves have recorded terminal
statuses; leaf areas sum to the band area 0.54 exactly (measured area
fraction 1.0000000000000004, float summation).

**Fractions (area): CROSSED 0 %, LOCALIZED 0 %, FAIL 100 %.** An honest
FAIL is a finding: at the pre-registered widths the Q1 engine cannot carry
a fat box even 0.05 time units forward in the early-epoch field.

## 5. Width ladder (diagnostic, added; printed as such)

Because no covering box certified, the pre-registered covering alone cannot
locate the feasible-width scale that S3-G3's extrapolation needs. Added
diagnostic (coordinator-owned addition, printed): single boxes at the band
center (−3.9, 0.5), aspect fixed at the grid's 2.25, widths halving from
0.05 to 3.125e-4 (the last two rungs extend the planned six):

| wx | death t | steps | G_enc at death |
|---|---|---|---|
| 0.05 | 1.049 | 201 | 2.8 |
| 0.02 | 1.115 | 473 | 6.0 |
| 0.01 | 1.190 | 780 | 10.1 |
| 0.005 | 1.289 | 1182 | 33.3 |
| 0.0025 | 1.410 | 1681 | 41.9 |
| 0.00125 | 1.563 | 2305 | 81.4 |
| 6.25e-4 | 1.747 | 3059 | 255.9 |
| 3.125e-4 | 1.949 | 3886 | 409.3 |

All FAIL_FDISC. **No certifying width w\* exists down to 3.1e-4.** The
per-halving death-time gain grows (0.066 → 0.202), i.e. the effective
enclosure growth rate λ_enc falls from 13.8/unit (t ≈ 1.08) to 3.4/unit
(t ≈ 1.85) — still ≫ the true flow's rate (§6).

## 6. Wrapping vs true flow growth (float layer, NOT rigorous; labeled)

DOP853 (rtol 1e-10) trajectories + the variational system V′ = JV give the
true flow-map amplification A_true(1→t) = ‖DΦ‖₂ at 9 points spanning the
band. **A_true(1 → 5.14) ∈ [17.3, 279] (median 127)** — the true flow is
only modestly expanding (the large shear J_xz = π²csc²(πz) ≈ 10–17 pairs
with small J_zx of opposite sign: near-elliptic local rotation, which
Lohner QR should track). The enclosure growth divided by the band-center
A_true at the matching epoch — the WRAPPING factor — compounds:

| wx | t_death | G_enc | A_true | wrapping |
|---|---|---|---|---|
| 0.05 | 1.049 | 2.8 | 1.27 | 2.2× |
| 0.01 | 1.190 | 10.1 | 2.35 | 4.3× |
| 0.00125 | 1.563 | 81.4 | 5.79 | 14.0× |
| 3.125e-4 | 1.949 | 409.3 | 7.40 | 55.3× |

The overhead is engine wrapping (first-order interval J over the fat
Picard box; the J-interval widths scale with box size and re-fatten the
box every step), not physics. This is precisely the "wrapping-effect
growth" datum the roadmap pre-registered for S3-G3.

## 7. S3-G3 — feasibility datum and full-O-1 extrapolation (PASS)

**Measured pilot data:** fraction (i)/(ii)/FAIL = 0/0/100 % (area); median
133 validated steps and 0.7 s per (early-dying) box; full-horizon box cost
86 s (bench, 16,958 steps at 5.0 ms); growth factors §4–§5; wrapping §6.

**Full O-1 needs** [0, T0] × supp ρ₀: x-span 7.0 ([−8.5, −1.5]) ×
z ∈ [0.05, 0.95] (the near-wall strips remain the O-2/O1 flux-bounded
exception class — NOT covered), horizon [0, 5.14] (1.24× the pilot's), per
certified |n_y| family (0.75 and 1.0; mirrors by symmetry; a continuum-n_y
statement would add an interval dimension on top). Cost model: 107 s per
surviving box (bench × 1.24), N = ⌈7/wx⌉·⌈0.9/(2.25 wx)⌉.

- **With THIS engine (uniform covering), bracketed by two projections of
  the ladder death-rate sequence:** model A (rate frozen at the last
  measured 3.4/unit — pessimistic): w\* ≈ 5.5e-9 → **9.4e16 boxes,
  2.8e15 core-hours** per family. Model B (rate ∝ 1/t from t = 1.95 —
  optimistic, matches dispersive smoothing): w\* ≈ 4.8e-7 → **1.2e13
  boxes, 3.7e11 core-hours** per family. **Verdict: infeasible in any
  bracketing — more compute does not discharge O-1 with this engine.**
- **Wrapping-free lower bound (from A_true; what ANY interval technology
  must pay):** initial widths ~0.25/A_true_max ≈ 9.0e-4 →
  **~3.5e6 boxes ≈ 1.0e5 core-hours per family** for a uniform covering;
  adaptive splitting pays the amplification locally (leaf count ~ the
  transported-volume integral) and would reduce this further. Expensive
  but conceivable.
- **The datum's conclusion:** the gap between 1e13–1e17 (this engine) and
  ~3.5e6 (wrapping-free bound) is entirely wrapping overhead. A full O-1
  covering requires higher-order enclosures (Taylor models / CAPD-COSY
  grade, which Q1's G5 already named as the next rung) plus adaptive
  subdivision — a technology upgrade, not a scale-up. That is what this
  pilot was pre-registered to decide.

## 8. Gates

| id | spec (pre-registered) | measured | verdict |
|----|------|----------|---------|
| S3-G1 | forward-integrate the reverse of one Q1 certified backward path; endpoint enclosure must CONTAIN the Q1 anchor (containment, not distance) | Path 134 forward replay, 1359 steps, bitwise-matching time grid: end box x ∈ [0.99999634, 1.00000366], z ∈ [0.41249793, 0.41250207] ∋ (1.0, 0.4125); widths 7.3e-6/4.2e-6; 7 s | **PASS** |
| S3-G2 | all 24 boxes (+ ≤ 1 round of 2×2 subdivisions) run to T0 or terminal blow-up, all checkpointed, zero uncertified gaps in the accounting | 24/24 parents terminal (FAIL_FDISC, t ≤ 1.0112) → all subdivided → 96/96 children terminal (95 FAIL_FDISC, 1 FAIL_WALL, t ≤ 1.0486); every leaf has a recorded terminal status + checkpoint file; leaf areas sum to the band area exactly; certified fraction 0 % — an honest FAIL-class outcome for every leaf, fully accounted | **PASS** |
| S3-G3 | feasibility datum printed: fraction (i)/(ii)/FAIL, median validated-step cost, width-growth factor distribution, extrapolated full-O-1 cost, honest | 0 %/0 %/100 %; median 133 steps, 0.7 s/box (early-death), 86 s full-horizon bench; growth factors 1.92/2.33/4.68 (cutoff-limited) + ladder curve 2.8→409; wrapping factors 2.2×→55× vs A_true ∈ [17, 279]; extrapolation bracketed 1.2e13–9.4e16 boxes/family (this engine: infeasible) vs ~3.5e6 boxes (wrapping-free bound): O-1 needs Taylor-model-grade technology + adaptivity | **PASS** |
| S3-G4 | scope honesty: pilot ≠ discharge; kill-bin theorem open; register unchanged | Stated throughout: this pilot discharges NOTHING of O-1; it prices the technology. Conjecture C obligations O-1…O-5 and the T4-W5 analytic theorem remain open; no grade language rises; within-model | **PASS** |

## 9. Key numbers

- G1 replay: growth 2.6× over 0.083 backward-depth units; end widths
  7.32e-6 / 4.15e-6; anchor contained.
- Engine cost: 5.0 ms/validated step (h = 2⁻¹²); full horizon
  16,958 steps ≈ 86 s/box; whole production 92 s + 203 s diagnostic.
- Covering: 24 parents + 96 children, 100 % blow-up class; parent deaths
  t ∈ [1.0010, 1.0112]; child deaths t ∈ [1.0017, 1.0486]; median
  133 steps / 0.7 s; total 51 s.
- Effective blow-up mode: f-disc zero at enclosure widths 0.23–0.48
  (the pre-registered width-1 cap never fired).
- Ladder deaths: t = 1.049 (wx 0.05) → 1.949 (wx 3.125e-4); λ_enc
  13.8 → 3.4 /unit; no w\* found.
- True amplification A_true(1→5.14): min 17.3, median 127, max 279
  (9 points; float layer). Wrapping factor 2.2× → 55×, compounding.
- Extrapolation (one family): this engine 1.2e13–9.4e16 boxes
  (3.7e11–2.8e15 core-h); wrapping-free bound ~3.5e6 boxes (~1.0e5
  core-h).

## 10. Honest scope and caveats

- **Pilot, not discharge.** Nothing here certifies any part of the
  subdomain as crossed-or-localized; O-1 stands exactly as open as Q2 left
  it. The pilot's value is the priced verdict on the technology path.
- **Spec deviations, all printed:** (a) FDISC/WALL/PICARD grouped as
  blow-up-equivalents (the pre-registered width-1 cap is unreachable
  because the f-disc dies first — a strictly earlier terminal of the same
  divergence); (b) the Picard pad is h-scaled for fat boxes (heuristic
  layer only; the verified-containment semantics are Q1's); (c) the width
  ladder and float-layer true-growth diagnostic are ADDED
  (coordinator-owned) because the covering alone left G3's extrapolation
  unquantifiable; the ladder got two extra rungs; (d) UNRESOLVED category
  defined for T0-straddling boxes (never occurred).
- The forward integrator inherits Q1's step semantics verbatim (mean-value
  S = I + hJ over the verified a-priori box, intersected with the direct
  interval evaluation); no change to the Q1 trust base (IEEE-754 float64,
  math.nextafter, mpmath.iv prec-120, exact model data).
- The true-growth layer (§6) is FLOAT diagnostics (DOP853 + variational
  matrix), not an enclosure; it informs the extrapolation only. The
  wrapping factors compare box growth against the band-CENTER trajectory's
  A_true; off-center trajectories differ by O(1) factors.
- The extrapolation models A/B are honest brackets, not measurements: A
  freezes the last measured enclosure rate (pessimistic), B decays it as
  1/t (optimistic). Both land ≥ 7 orders of magnitude above the
  wrapping-free bound; no conclusion depends on which bracket is right.
- The full-O-1 estimate excludes the near-wall strips z ∉ [0.05, 0.95]
  (O-2's flux-bounded exception class), assumes per-family coverage at
  discrete n_y ∈ {0.75, 1.0} (+ mirrors by the re-verified z → 1−z
  symmetry), and starts at t = 0 where the pilot starts at t0 = 1.0; the
  early [0,1] epoch is, if anything, harder (larger |f″/f| scales).
- The 1 FAIL_WALL is a Picard-inflation artifact at an already-blown box
  (final z ∈ [0.46, 0.69], nowhere near the wall); it is classified in the
  same failure class it would otherwise reach one step later.
- Growth factors of §4 are death-cutoff-limited (max 4.68 because boxes
  die at width ~0.25–0.5 from 0.1-scale starts); the ladder (§5) is the
  uncensored growth measurement.

Within-model; nothing here bears on nature. The kill-bin theorem — and
with it the analytic T4-W5 proof — remains the terminal open item.
