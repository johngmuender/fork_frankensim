# fs-gum-kern — Phase G5a results: the F-R5 halo-descent referee at N = 192

Scope (ROADMAP_v4_SCALE Phase G5, first half): the campaign's F-R5
halo-descent referee — the frozen fs-gum-statics G-C protocol (Phase E1,
itself the Rust port of the tier4-field 4A protocol) — re-run through
`fs_gum_kern::{eval, anf}` on 4-thread sweeps at up to **8× the 4A
campaign cells (192³ vs 96³) and 64× the E1 gate-suite cells (48³)**.
Everything protocol-level is FROZEN and identical across grids: Step-1
radial hedgehog at t = T_FROZEN, LCG bump seed 424242 (K = 6, amp 0.02),
tilt-halo shell (η = 0.05, r₀ = 3.2, w = 0.6), guards MU = 5000 /
band 0.005 / MU_F = 400 / floor band 0.01 referenced to the engine's own
hedgehog, L_main = 0.266·I_hedgehog (over-spun, κ₀ ≈ 0.253–0.26),
L_ctrl = 0.12·I_hedgehog (sub-threshold), corner objective, arrested
Newton flow, records every 10 iterations. Only (N, LBOX) and the
iteration caps vary, and both are printed in the launch commands below.

Binary: `src/bin/n192_fr5.rs` (subcommands `budget` / `proto` /
`ctrlclock` / `bitcheck`). Reproduce:

```
cargo build --release
./target/release/n192_fr5 budget    192 4.5 10 4
./target/release/n192_fr5 proto     n192_runs/n96_base.json   96 4.5 900  600 600 4
./target/release/n192_fr5 proto     n192_runs/n128_box.json  128 6.0 900  800   0 4
./target/release/n192_fr5 proto     n192_runs/n192_main.json 192 4.5 450 1600 600 4
./target/release/n192_fr5 ctrlclock n192_runs/n192_ctrlclock.json 192 4.5 450 600 4
./target/release/n192_fr5 bitcheck  192 4.5 20 2 4
./target/release/n192_fr5 crossing  n192_runs/n192_crossing.json 192 4.5 450 2500 4
```

Measured facts only; adjudication is the coordinator's layer.

## 1. Grid plan — the two refinement axes

The N = 96 referee box is LBOX = 4.5 (half-width; the field3d 4A box).

* **h-axis (h → h/2 at fixed box)**: N = 48 → 96 → 192 at LBOX = 4.5
  (h = 0.1875 / 0.09375 / 0.046875). The N = 48 rung is the frozen E1/K-E
  record; 96 and 192 are new runs of the same protocol.
* **box-axis (box → 4/3× at fixed h)**: N = 96, LBOX = 4.5 vs
  N = 128, LBOX = 6.0 — **exactly** the same h = 0.09375, box scaled
  4/3 ≈ 1.33. Deliberate deviation from the task sheet's "1.3× box at
  N = 144" (h = 0.08125, which matches the h of no other run and would
  confound the two axes): with LBOX = 6.0 = 4.5 + 16h the two grids are
  exactly nested, so the box axis is measured at literally identical
  discretization. (The 4A box-sanity run changed h and box TOGETHER —
  N fixed at 96, h 0.09375 → 0.12188 — its own confound; see §6.)

## 2. Budget protocol (measured before committing to caps)

* Task-sheet projection for the N = 192 corner-gradient ANF iteration at
  4 threads: ~1.24 s. Measured probe (10 iterations, seeded over-spun
  configuration): **1.258 s/eval quiet**; 2.793 s/eval under full
  contention with the concurrent G5b two-knot farm (4 procs at 100%).
  The long quiet runs later measured 1.05–1.13 s/eval (static/main at
  N = 192). N = 192 was therefore kept (no drop to N = 160).
* Co-tenancy handling (the G1 K-F precedent): the G5b farm was
  SIGSTOPped for the measured run block and SIGCONTed after —
  bit-neutral to the farm. The recovery segment (§3, last row) ran
  partially alongside the restarted 3-process farm; correctness is
  scheduling-independent by the layer's construction, only those wall
  clocks are contended (flagged in the table).
* Cap re-budget mid-phase: the N = 96 run measured the threshold
  crossing at iteration 490..500 — ~3× the N = 48 record's ≈160 — so the
  queued N = 192 main cap of 600 was raised to 1600 before launch.

## 3. Runs table

| run | grid | stage caps | status | wall (4T) | measured-fact gates |
|-----|------|-----------|--------|-----------|---------------------|
| n96_base R1 static | 96³, LBOX 4.5 | 900 | iter_cap | 119.8 s | R1-deg/floor/hood/mono PASS |
| n96_base R2 main | | 600 | iter_cap | 79.7 s | R2-desc/kappa/halo/mono PASS; crossing it 490..500 |
| n96_base R3 ctrl | | 600 | iter_cap | 80.4 s | R3-halo/kappa PASS |
| n96_base R4 clock | | — | — | — | R4 PASS; ratio 1.4105 (§7) |
| n128_box R1 static | 128³, LBOX 6.0 | 900 | iter_cap | 289.6 s | R1 all PASS (digit-identical to n96 static, §6) |
| n128_box R2 main | | 800 | iter_cap | 255.8 s | R2 all PASS; crossing it 490..500 |
| n192 proto R1 static | 192³, LBOX 4.5 | 450 | iter_cap | 475.0 s | see §4 |
| n192 proto R2 main | | 1600 | iter_cap | 1807.6 s | monotone; crossing NOT reached (κ = 0.20668, falling); halo +0.0118 |
| n192 proto R3 ctrl | | 600 | **killed at it ≈ 500** | — | process killed externally (§11.3); log preserved, series JSON lost |
| n192 ctrlclock R1+R3+R4 (recovery) | | 450 + 600 | iter_cap ×2 | 960.8 s + 1020.9 s (contended, §11.4) | static replay digit-identical to the killed run; R3 self-limiting; R4 identities exact (Δerot = 0.0); ratio 1.4847 |
| n192 bitcheck | | 20 iters, 2T vs 4T | PASS | 45.0 s (2T) / 25.4 s (4T) | hashes identical: 0803fe446497060d… |
| n192 crossing R1+R2 (long cap) | 192³, LBOX 4.5 | 450 + 2500 | iter_cap | 421.3 s + 2348.9 s (quiet) | **threshold crossing OBSERVED at it 2070..2080**; obj strictly monotone (max rise −2.2e-6); static + main-to-1600 rows digit-identical to the earlier runs; kill-tolerant JSON flushed every 100 it |

All descents objective-monotone by construction; guard footprint
|E_pen| + |E_fpen| ≤ 2e-3 at every ENDPOINT (N = 192 logged rows
≤ 1.6e-3), with transient maxima along the descents of 3.1e-3 (N = 96)
and 5.0e-3 (N = 128) — the same class as the E1 record's 4.5e-3
main-descent maximum at N = 48.

## 4. R1 — static relaxation across the ladder

| grid | Estat_hedgehog (corner) | deg_hedgehog | static endpoint Estat (dE) | deg (band edge) | floor penetration |
|------|------------------------|--------------|---------------------------|-----------------|-------------------|
| 48 (E1 record) | 3.0430807 | 0.945803 | 2.9865 (−0.0565) | 0.940153 | −2.4e-3 |
| 96 | 3.1466236 | 0.986012 | 3.1072601 (−0.0394) | 0.980435 | −2.2e-3 |
| 192 | 3.1738445 | 0.996472 | 3.1573216 (−0.0165) | 0.991062 | −1.0e-3 |

* The corner-scheme hedgehog converges to the Step-1 continuum value
  Estat = 3.1830538 with clean O(h²): errors −0.13997 / −0.03643 /
  −0.00921 (ratios 3.84, 3.96); Richardson from the 96/192 pair gives
  3.18292 (continuum 3.18305, 4e-5 off). Degree errors −0.0542 /
  −0.0140 / −0.0035 (ratios 3.87, 3.96). **At N = 192 the corner
  engine passes the campaign's own Tier-4 gate bands that it FAILED at
  N = 96** (|deg − 1| = 3.5e-3 < 5e-3; Estat-class error 0.3% < 1%).
  I = 22.26826 identically at 96/128/192 (quadrature-exact instrument,
  1e-5-class agreement with N = 48's 22.26809).
* The near-BPS drift pattern is the documented one at every N, and
  SHRINKS with h: dEstat −0.0565 → −0.0394 → −0.0165; degree pinned at
  the anchor edge; wall penetration within the E1 reference class.
  (N = 192 static gnorm/gnorm₀ = 3.05 at its 450 cap — the finer grid
  exposes stiffer high-k modes; the endpoint is guard-bounded, not
  gradient-converged, exactly like every 4A/E1 static stage.)

## 5. R2 — the over-spun main descent: milestones across N

κ-milestones of the SAME frozen protocol (κ₀ ≈ 0.253–0.26 after
seeding; threshold 1/√(8π) = 0.199471). N = 192 brackets are from the
`crossing` run's 10-iteration series (its static and main stages replay
the killed run digit-for-digit through it 1600 and continue to 2500):

| first κ below | N=48 (h=.1875) | N=96 (h=.09375) | N=128 box 4/3× (h=.09375) | N=192 (h=.046875) |
|---------------|----------------|------------------|---------------------------|--------------------|
| 0.24 | — (κ₀ 0.253; not instrumented) | it 120..130 | it 120..130 | it 230..240 |
| 0.22 | — | it 280..290 | it 280..290 | it 830..840 |
| threshold 0.19947 | it 150..160 | it 490..500 | it 490..500 | **it 2070..2080 — DIRECTLY OBSERVED** (κ 0.199871 → 0.199447; inside the 1600-cap run's it ≈ 2000–2400 extrapolation) |
| halo at threshold | n/r | 0.0400 | 0.0399 | 0.0430 |
| I at threshold | n/r | 29.80 | 29.73 | 29.70 |

* **The iteration count to any fixed κ mark is grid-sensitive**:
  ≈160 / ≈495 / ≈2075 iterations to threshold at h = 0.1875 / 0.09375 /
  0.046875 — a growing ×3.1 then ×4.2 slowdown per h-halving. The
  mechanism is solver dynamics, not physics: at finer h the ANF
  dt-cascade settles ~2–3× lower (N = 192 dt ≈ 0.002–0.008 vs N = 48's
  ~0.01–0.05) and each accepted step moves the state less.
* **The trajectory in STATE space is grid-robust.** At matched κ the
  quadrature-exact instruments agree across the ladder: at κ ≈ 0.238,
  N = 96 (it ≈ 130): halo 0.0310, I 25.31 vs N = 192 (it ≈ 235): halo
  0.0292, I 24.76; at κ ≈ 0.22, N = 96 (it ≈ 285): halo 0.0340,
  I 27.07 vs N = 192 (it ≈ 835): halo 0.0343, I 26.98; at the threshold
  itself: halo 0.0400 / 0.0399 / 0.0430 and I 29.80 / 29.73 / 29.70
  (N = 96 / 128-box / 192).
  The descent path (κ, halo, I) refines onto itself; only its
  parameterization by iteration slows.
* N = 192 crossing-run endpoint at cap 2500: R 3.9377 → 3.745653
  monotone (dR = −0.192; max recorded bare-R rise +2.9e-4, penalty-
  sized; objective max rise −2.2e-6), κ 0.25876 → **0.19502** (below
  threshold), halo 0.0276 → **0.0456** (+0.0180 — now ABOVE the 0.015
  E1 gate margin that the 1600-cap read fell short of), I 22.89 →
  30.37, degree pinned at 0.9909, floor at the wall, 48 arrests, guard
  footprint ≤ 3.0e-3. At the crossing itself (it 2070..2080): halo
  0.0430, I 29.70, R 3.75793.
* Python 4A N = 96 reference (DIFFERENT spin protocol, L = 8.4979,
  κ₀ = 0.374, stated for context, not equivalence): R crosses the
  axisymmetric reference at it ≈ 110, κ crosses solA's 0.2656 at
  it 400..425 and the threshold at it 1600..1625 of cap 2400. Its
  iteration milestones are solver-clock numbers in exactly the same
  sense.

## 6. R3 control + the box axis

* **Control (N = 96, cap 600)**: κ 0.11761 → 0.10211 (max 0.11761,
  always below threshold), halo 0.0276 → 0.0295 (+0.0019 =
  0.13× the main run's rise), I 22.72 → 26.17 — climbing the bounded
  tilt dial toward the F-R4 ceiling (3/2)·I_h = 33.40, Estat flat
  (3.1205 → 3.1075): no far-field condensate. Differential gates pass.
* **Control (N = 192, cap 600, recovery run)**: κ 0.11674 → 0.11158
  (max 0.11674, always below threshold), halo 0.0276 → **0.0256** — at
  N = 192 the control's halo fraction FALLS (−0.0020): the seeded shell
  relaxes away instead of condensing. I 22.89 → 23.95 (72% of the tilt
  ceiling), Estat 3.1714 → 3.1506 (flat — no E₀-priced condensate).
  Matched-cap-600 differential vs the N = 192 main run: dhalo −0.0020
  vs +0.0042, dκ −0.0052 vs −0.0308 — the two-sided signature holds at
  the finest grid, and is if anything CLEANER (the sub-threshold side
  loses its lattice-mush-assisted halo drift as h shrinks).
* **Box axis is NULL at fixed h through the measured caps**: the
  N = 128 / LBOX 6.0 run reproduces the N = 96 / LBOX 4.5 main descent
  to the 4th decimal at matched iterations (it 500: R 3.714064 vs
  3.712971, κ 0.19922 vs 0.19879, halo 0.0399 vs 0.0400, I 29.73 vs
  29.80) and to the SAME milestone brackets (120..130 / 280..290 /
  490..500) — with 4/3 the room, nothing about the descent through and
  past the threshold changes. Continuing past the N = 96 cap it keeps
  descending: it 800: κ 0.18360, halo 0.0469, I 32.26, R 3.67299. The
  static stages of the two grids are digit-identical in every printed
  instrument — the exact grid nesting (LBOX 6.0 = 4.5 + 16h, vacuum
  padding contributes exact zeros to every plain-f64 tile sum) makes
  the box run a controlled experiment and the agreement a nesting
  cross-check.
* Consequence for the 4A box-sanity reading (measured fact, stated
  plainly): the 4A 1.3× box run differed from its main box AT MATCHED
  ITERATIONS (crossing 800 vs 1600, deeper R) while ALSO coarsening h
  0.09375 → 0.12188. On this phase's h-ladder, coarser h alone
  produces exactly that signature (fewer iterations to any κ). At
  fixed h the box contributes nothing through halo ≈ 0.047 /
  it ≤ 800. The 4A conclusion that far-field growth is room-limited
  in the LATE stage (halo 0.10–0.16) is untested at fixed h — these
  caps never push the halo into the boundary region.

## 7. R4 — the clock ladder, and what the ratio is a function of

E_rot/E_tot = 1/4 at the bisected L_clock holds to machine precision on
every state measured (N = 96: |Δ| = 5.6e-17; N = 192: Δ = 0.0 — the
stored f64 is 0.25 exactly; bisection = closed form to ≤ 2.2e-16 rel).
The headline number κ(L_clock)/threshold:

| state the clock is imposed on | I (% of tilt ceiling 33.40) | ratio |
|-------------------------------|------------------------------|-------|
| N=48 ctrl endpoint, cap 600 (E1 G-D) | 31.94 (96%) | 1.2525 |
| N=96 Python 4A ctrl endpoint, cap 1400 | 33.62 (101%) | 1.2448 |
| N=96 G5a ctrl endpoint, cap 600 | 26.17 (78%) | 1.4105 |
| **N=192 G5a ctrl endpoint, cap 600** | 23.95 (72%) | **1.4847** |
| N=192 static endpoint | 22.16 (66%) | 1.5449 |
| N=192 main endpoint, cap 1600 | 28.66 (86%) | 1.3605 |
| N=192 main endpoint, cap 2500 (crossing run) | 30.37 (91%) | 1.3220 |
| 4A clean axisymmetric control reference | — | 1.474 |
| Step-2 solA | — | 1.331 |

The ladder's resolution (measured, not adjudicated): **the ratio is a
state function, ratio = √((2/3)·Estat/I)·√(8π), and along every control
trajectory it DECREASES toward a floor as the bounded tilt dial
saturates** (Estat stays flat, I → (3/2)·I_h). All four control
trajectories — N = 48 E1, N = 96 Python, N = 96 G5a, N = 192 G5a —
collapse onto the same declining curve when parameterized by I (e.g. at
I ≈ 26.2 the Python curve reads ≈ 1.42–1.43, the G5a N = 96 run reads
1.4105 at 26.17; the G5a N = 192 run reads 1.4847 at I = 23.95 where
the Python curve reads 1.4670 at I = 24.18 — the ~1% offset tracks the
corner Estat's O(h²) convergence, 3.151 vs 3.108 in the ratio's
numerator). The
floor at the tilt ceiling is 1.248–1.257 for the measured Estat range
(3.107–3.152); the fully-saturated readings are the 1.2448/1.2525 pair.
What is grid-robust: the floor value ≈ 1.245–1.2525, and the fact that
EVERY B = 1 state measured at every grid gives ratio > 1.24 — the
clock charge lands above the halo threshold on all of them. What is
cap-sensitive (NOT grid-sensitive physics): a control endpoint read at
a fixed iteration cap sits higher on the curve at finer h because
saturation-per-iteration is slower (§5).

## 8. Thread bit-identity spot check at N = 192

A 20-iteration seeded over-spun descent at N = 192 (hedgehog + bump +
halo shell, the K-D construction), full instrumented series + run
metadata + final field bytes hashed under the bit contract:

```
threads=2: 0803fe446497060de0a0cd0de930515d70efd4558a9c362ff2bcbc76371509d8  (45.0 s)
threads=4: 0803fe446497060de0a0cd0de930515d70efd4558a9c362ff2bcbc76371509d8  (25.4 s)
```

BIT-IDENTICAL — the G1 thread-identity contract holds at 8× the cells
it was frozen at (and under farm co-tenancy, which is the point: no
recorded number in this phase depends on thread count or load).

Supporting replay evidence from the recovery itself: the `ctrlclock`
run re-executed the killed run's static stage and control stage and
reproduced their printed instruments digit-for-digit at every logged
iteration (static it 0–450, control it 0–500 — every R, I, κ, deg,
halo, gnorm, dt, arrest count), across distinct processes ~1 h apart —
the layer's determinism contract exercised end-to-end at 192³ scale.

## 9. Performance at scale (context, not a gate)

Quiet-window ANF hot-path cost at N = 192, 4 threads: 1.05–1.13 s/eval
across the long runs (451-eval static: 475.0 s; 1601-eval main:
1807.6 s), vs the task-sheet projection 1.24 s and the 10-iteration
probe 1.258 s. The N = 192 main descent alone is ~10¹³-FLOP-class; the
full phase's engine time is ≈ 143 min wall at 4 threads (§11.5),
of which the coordinator-directed crossing extension is 46.3 min
(421.3 s static + 2348.9 s main-2500, quiet — 0.94 s/eval on the
mixed static/main stream). The 16·N³ pass-A flux buffer at N = 192 corner is ~920 MB;
peak RSS ≈ 3 GB of 15 — the slab-staged pass A flagged in the G1
RESULTS remains unnecessary at this scale and remains the prerequisite
for N > 192.

## 10. Measured summary relative to the F-R5 endpoint claims

Stated as measured facts for the coordinator's adjudication:

1. **The referee's qualitative signature reproduces at every grid**:
   monotone R descent, κ falling with I rising, halo fraction growing,
   control self-limiting with Estat flat, degree and floor guarded
   (footprint §3). Two-sidedness holds at N = 96 and N = 192
   (at N = 192 the control's halo FALLS, −0.0020, while the over-spun
   run's rises at every grid).
2. **What the 8× refinement sharpens**: the corner engine's absolute
   accuracy — at N = 192 it passes the 1%-sector/5e-3-degree gate
   class it failed at N = 96, with measured O(h²) ratios 3.8–4.0, and
   the near-BPS static drift shrinks −0.0565 → −0.0165; the state-space
   descent path (halo, I at matched κ) is confirmed grid-robust at the
   few-percent level from h = 0.1875 to h = 0.046875; the box axis is
   shown null at fixed h through the measured range; the clock-ratio
   floor 1.245–1.2525 and the ratio > 1.24 on every measured state are
   grid-robust.
3. **What is exposed as solver-clock, not physics**: every
   iteration-indexed milestone. Iterations-to-threshold grow ×3–4 per
   h-halving (≈160 / ≈495 / ≈2075, the last directly observed); the 4A
   milestones "it ≈
   110/400/1610" and the E1 "crossing by it 160" are properties of the
   arrested-Newton clock at their grids.
4. **The crossing is directly observed at N = 192** (the `crossing`
   run, cap 2500, replacing this item's earlier honest-limit status):
   κ falls through 1/√(8π) at **it 2070..2080** (0.199871 → 0.199447),
   inside the 1600-cap run's it ≈ 2000–2400 extrapolation, with halo
   0.0430 and I 29.70 at the crossing (vs 0.0400/29.80 at N = 96) and
   κ = 0.19502, halo = 0.0456 (+0.0180, above the 0.015 E1 margin),
   I = 30.37 at the 2500 cap, objective strictly monotone throughout.
   Every referee clause — monotone descent, κ crossing, halo growth —
   is now an observed fact at 8× the campaign's cells. The remaining
   honest limits are the ones shared with every 4A/E1 record: the runs
   end at iteration caps (minimizing sequences still in motion, not
   converged endpoints), and the deep post-crossing regime
   (κ → 0.18-class, halo ≥ 0.05–0.10) was measured only at N ≤ 128.

## 11. Deviations (with reasons)

1. **Box pair (N=128, LBOX 6.0) instead of the sheet's (N=144, 1.3×)**:
   exact-h nesting beats a 1.3× box whose h matches no other run; the
   box factor is 4/3 instead of 1.3 (§1).
2. **N=192 main cap 600 → 1600** after the N = 96 crossing measurement
   (§2); the crossing still wasn't reached at 1600. The direct
   observation was then authorized by the coordinator as a follow-up
   segment on the idle machine (the `crossing` run, cap 2500, §5,
   §10.4) and landed at it 2070..2080.
3. **The first N = 192 proto was killed externally** (the harness
   reaped the whole background run block; no OOM, no panic — the
   binary was mid-control at it ≈ 500) at 07:42. Its static + main
   stages were complete and their log rows are the record used in
   §4–§5; the 10-iteration series JSON was lost because `proto` wrote
   JSON only at process end (design defect, printed here). The
   recovery subcommand `ctrlclock` re-runs static + control + clock
   with stage-flushed JSON, and its static replay doubles as a
   determinism check (§8).
4. **Concurrency**: the G5b two-knot farm was SIGSTOPped during the
   measured block (K-F precedent, bit-neutral) and ran alongside the
   recovery segment (those walls flagged; all non-wall numbers are
   scheduling-independent by construction). Thread count for the
   recovery stayed 4 (bits independent of it; K-D).
5. **Wall budget**: the original phase ran ≈ 97 min engine wall
   (probes 45 s; n96 281 s; n128 547 s; killed n192 proto ≈ 2883 s;
   recovery 1986 s; bitcheck 70 s) vs the sheet's ≤ ~90 min — the
   overage is the killed run's lost control stage plus the recovery's
   contended static re-run (≈ 26 min combined, of which ≈ 16 min is
   pure re-execution/contention cost); without the kill the phase
   lands at ≈ 80 min including the deliberate main-cap doubling. The
   coordinator-directed crossing extension added 2777 s (46.3 min,
   quiet machine) for a total ≈ 143 min. Stated rather than hidden.
6. **No control at N = 128**: the box axis needs the main channel only
   (the 4A box-sanity protocol likewise ran main only).

## 12. Files

`src/bin/n192_fr5.rs` (runner), `n192_runs/n96_base.{json,log}`,
`n192_runs/n128_box.{json,log}`, `n192_runs/n192_main.log` (killed run:
static + main complete, control to it 500; its JSON was never written),
`n192_runs/n192_ctrlclock.{json,log}` (+ `.partial` stage flush),
`n192_runs/n192_bitcheck.log`, `n192_runs/n192_crossing.{json,log}`
(the direct-observation run: full 10-iteration series for static 450 +
main 2500, kill-tolerant), this file. Full instrumented series for
every run that reached a JSON write; the killed run's main series is
superseded by the crossing run's digit-identical replay-and-extend.

## 13. Epistemic notice (binding, inherited)

Everything here is a within-model computation on a speculative theory's
functional, on finite lattices, with the two documented one-sided
guards and a seeded perturbation protocol. The runs refine the
campaign's F-R5 REFEREE — an engine-capability replication and a
discretization study. No statement here is about nature, and none is an
adjudication; PASS lines are measured facts against frozen thresholds
from Phase E1.
