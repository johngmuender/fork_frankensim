# H2.2 — NR-D2 P4 reconciliation test: two mismatched textures on the certified engines
## SDiff-reachability of the "reconciliation field" (the P4 GAP of T1 §2.3), measured

**Phase:** H2.2 (numeric follow-on item 3 of `T1_NR_gates.md` §6).
**Date:** 2026-07-16.  **Engines:** `fs-gum-field` (measurement path, radial
solver), `fs-gum-statics` (guarded ANF protocol) — neither modified; seeding
pattern from `fs-gum-twoknot`.  **Solver:** `h22_solve/` (this directory),
deterministic (no RNG anywhere; analytic seeds; sequential f64 loops; same
binary + args ⇒ bit-identical JSON).  **Compute:** ~24 min total (five runs +
seed scan), grid 96³ at the frozen h = 0.09375 (certified LBOX = 4.5 box).
**Inputs:** `corpus/NR-D2-Collision-Problem-x-Preview-v0_1.md` §2 (P4),
`theory-audit/T1_NR_gates.md` §2.3 + §5.3, `tier2-closure/axi_FR4_NOTE.md`.

**Epistemic frame (binding).** Within-model numerical engineering on a
speculative theory's functional.  Everything below is a measured fact about
the campaign's own certified functional and about the corpus's argument
STRUCTURE as transferred to it — never about nature, and never corpus-level
adjudication (that is the coordinator's).

---

## 1. Construction and mapping (the spec-reconstruction risk, printed in full)

NR-D2 §2's setup: *"Two overlapping double-twist textures of different axes
mismatch at O(1) in the director; the naive cost is K/p² × p³ = f²p per
crossing — but the mismatch **is** SDiff-reachable (a smooth interpolating
reorientation lives on the soft manifold), so the medium screens it at the
ε-lifted cost V_int ≈ ε·f²p."*  The engine has no cholesteric/Frank-elastic
sector, so the geometry is reconstructed through the hedgehog-class
machinery (the same interpretation class the twoknot campaign used for
App I.2):

| corpus object | engine realisation |
|---|---|
| texture patch (pitch p) | B = 1 hedgehog knot (size R* = 2^{5/6}) |
| two overlapping patches | product ansatz q = q₁(x+s)·(a q₂(x−s) a\*), d = 2s = 1.6875, x = d/R* = 0.947 (deep overlap) |
| mismatched local anisotropy axes, O(1) | knot 2 isorotated θ = π about the separation axis ("repulse" channel).  By hedgehog equivariance q(Rx) = a q(x) a\* this is exactly a relative rotation of the second texture's internal frame — the most SDiff-flattering mismatch (its reconciliation looks most like a spatial rearrangement) |
| reconciled configuration | identical geometry, identity relative rotation ("align" channel) |
| the soft manifold | the compositional SDiff orbit Q → Q∘T, det DT ≡ 1 (F-R4's reading; the corpus never defines its manifold more precisely — that missing text is itself part of the P4 gap) |
| ε-lifting | the t-weighted kinetic sectors E2, E4: at the frozen dial ε = t(E2+E4)/(E6+E0) = 0.05, these are the ONLY sectors that vary on the SDiff orbit (F-R4 lemma); E6, E0 are frozen |
| naive wall f²p | dE_naive = E_seed(mismatch) − E_seed(aligned) — the engine's own realised O(1)-mismatch cost |
| screening scale ε·f²p | S = 0.05·\|dE_naive\| = **0.1315** engine units (unit conversion: f²p ↔ \|dE_naive\| = 2.6296, ε ↔ 0.05, stated once, used throughout) |

**Mapping caveats, up front:** (i) the corpus medium is Frank-elastic chiral
with a pitch; the engine is the quaternionic (2+4+6+0) sector — the transfer
tests the *mechanism class* (value-function freezing vs ε-lifting on SDiff),
not the corpus's magnitudes; (ii) the corpus's operative ε band is
[10⁻⁵, 6×10⁻⁴] — per task the engine dial ε = 0.05 is mapped instead
(measured at the relaxed single knot: 0.0472); (iii) chirality/pitch
structure of double-twist is not represented — the O(1) axis mismatch is.

**Protocol.** Cubic `Field3` n = 96, half = 4.5 (frozen h = 0.09375); frozen
radial profile (N = 4000, rmax = 6); corner-scheme guarded ANF exactly per
the statics contract.  Guard deviations from the twoknot convention, both
documented in-source: (a) the Bogomolny wall is referenced to the
**continuum per-unit-degree B = 2 bound** (fgap_ref = 32√2/15 − 0.01)
instead of the seed's own gap, because the overlapping seeds carry ~2.4–5.0
units of E6+E0 excess and a seed-referenced wall would freeze the value
sectors even in the FULL relaxation — i.e. would manufacture the very
result under test; the wall keeps its one purpose (blocking the near-BPS
lattice collapse below the continuum bound); (b) deg_ref = min over the
{align, channel} seed degrees (align run: 1.95399, repulse run: 1.95069 — a
0.0033 asymmetry; endpoint anchor penalties 2.8×10⁻⁴ on both sides bound
its effect).

**The five runs** (+ seed scan), all deterministic single processes:
`single` (800 it ANF), `pair align` and `pair repulse` (1200 it guarded ANF
from the two seeds), `sdiff repulse` (SDiff-projected descent, below),
`twist repulse` (exact volume-preserving twist-family scan, 21 points),
`seedscan` m = 6..20.

**SDiff-projected relaxation — the documented choice:** projection of the
descent step onto divergence-free displacement flows, applied
*compositionally with analytic reseeding* (the F-R4 family-A pullback
pattern: the state is always q_seed(X(x)) through the composed departure
map — never a field resampled from itself, so cumulative interpolation
diffusion cannot fake a value-sector change).  Per step: advective steepest
direction w_c = Σ_a g_a ∂_c q_a from the engine's exact guarded gradient;
Leray projection u = w − ∇φ with the CONSISTENT wide Laplacian
(div_central∘grad_central; the first attempt used the compact 7-point
operator and the descent exploited the projector inconsistency for a −0.48
value leak — a measured lesson, kept in the log); warm-started CG to 10⁻⁸;
map composition X ← X(x − τu/|u|∞) with trilinear interpolation of the
smooth displacement; arrest on objective increase OR on a frozen-quantity
budget violation (|Δ(E6+E0)| ≤ 0.05, |Δdeg| ≤ 0.01 — exact SDiff keeps both
constant, so the budget enforces the constraint without biasing the
objective).  Map distortion measured directly: max|det∇X − 1| = 7.8×10⁻²
(localised), mean 2.5×10⁻⁵.

**Exact twist family:** T_α(x, y, z) = (x, R(α·χ(x))(y, z)), χ = smoothstep
on [−s, s]: det DT = 1 **exactly**.  This is the corpus's "smooth
interpolating reorientation" executed literally as an SDiff element — at
α = π it locally undoes the π mismatch at knot 2 (the rotation axis passes
through both centres).

---

## 2. The cost table (engine units, corner scheme, bare E_static)

| quantity | value | note |
|---|---:|---|
| E_ref = E(align pair, relaxed 1200 it) | **6.4005** | the reconciled reference (iter_cap; tail flatness 2.8×10⁻³) |
| 2 × E(single, relaxed) | 6.2176 | far/uniform reference |
| E_seed(align) | 11.4159 | E6+E0 = 11.0237 |
| E_seed(repulse = mismatch) | 8.7863 | E6+E0 = 8.4624 |
| **dE_naive** (mismatch − aligned, seeds) | **−2.6296** | **SIGN-INVERTED vs the corpus's +f²p**; value part −2.5613 (97.4%), kinetic part −0.0683 (2.6%) |
| **cost_full** = E(repulse relaxed) − E_ref | **−0.0755** | achieved by moving E6+E0 by −2.449 (99.5% of the total −2.461 reduction) — an off-orbit path by construction |
| **cost_sdiff_projected** = E(sdiff endpoint) − E_ref | **+2.3374** | = 17.8 × S |
| **cost_twist_best** (α = 0, the family's minimum) | **+2.3858** | the twist family never improves on the identity |
| screening scale S = 0.05·\|dE_naive\| | **0.1315** | the corpus's ε·f²p, mapped |
| rigorous SDiff value floor: (E6+E0)_seed(mis) − E_ref | **+2.0619** | = **15.7 × S**; every SDiff image of the mismatch seed has E ≥ this (t(E2+E4) ≥ 0 + frozen value sectors) |

**Seed scan (sign structure, m = 6..20, x = 0.63..2.10):**
E_seed(repulse) < E_seed(align) at EVERY overlapping separation, the
difference peaking at −2.95 (x ≈ 1.16) and → 0⁻ as the overlap vanishes;
the value-sector share of the difference is 94–98% throughout.  **There is
no separation at which the corpus's positive naive wall exists in this
sector.**

## 3. The SDiff-projection result

* **Projected descent: nothing beyond numerical leak.**  Total reduction
  achieved before stalling against the frozen-quantity budget:
  0.0484 — of which the measured E6+E0 leak is 0.0500 (clamped at budget;
  12 budget hits, then τ-stall at iteration 55).  Reduction beyond the
  leak: **−0.0016 ≤ 0** — the kinetic sectors actually rose slightly.  The
  full ANF from the same seed releases **2.4613** (50.8× more), of which
  2.4490 (99.5%) is value-sector.  SDiff moves capture **0%** of the
  legitimate relaxation within the measured error bar.
* **Exact twist family (zero-leak exhibit): the interpolating reorientation
  does not reconcile.**  Along α: 0 → π, E0 is flat to **1.8×10⁻¹⁵**
  (machine; the value density is exactly axisymmetric about the twist axis
  by hedgehog equivariance, so the exact map preserves it cell-by-cell),
  E6 flat to 0.39% (pure quadrature of the twist gradient; deg flat to
  1.7×10⁻³), while E_static **rises monotonically** (+0.0033 at α = π) even
  as the pointwise distance to the reconciled seed drops from 6.67 to 3.06.
  The map moves the configuration "toward" reconciliation spatially and
  still cannot lower the energy: the mismatch is stored in the value
  DISTRIBUTION, which no pullback can change — the F-R4 frozen-value lemma,
  confirmed numerically on this engine's functional in its guarded descent
  form.

## 4. The floor test (the ε-lifted cost as a LOWER bound)

For SDiff-restricted reconciliation the floor holds with a large margin:
any SDiff image of the mismatch seed satisfies E − E_ref ≥ 2.0619 =
**15.7 × S** (rigorous, from the frozen value sectors alone); the measured
projected-descent endpoint sits at 17.8 × S and the twist-family minimum at
18.1 × S.  **cost_sdiff ≥ S: PASS (×17.8).  cost_twist ≥ S: PASS (×18.1).**
The unrestricted relaxation lands BELOW the reference (cost_full = −0.0755,
so cost_full ≥ S: FAIL) — but it gets there only by changing the value
distribution (−2.449 of E6+E0), i.e. by leaving the orbit entirely; it is
not a counterexample to the floor, it is the demonstration that the
corpus's reconciliation is **not an SDiff move**.

## 5. What this does to NR-D2's P4 gap (measured facts; coordinator adjudicates)

1. **The reachability parenthesis fails in the engine's conventions.**  P4
   asserts the mismatch "is SDiff-reachable … on the soft manifold."  On
   the campaign's own functional, with the soft manifold read as the
   compositional SDiff orbit (the only reading the record defines, per
   F-R4), reconciliation-class moves change nothing: the projected descent
   achieves zero beyond its measured leak, and the literal interpolating
   reorientation (exact SDiff element) raises the energy while E6+E0 stay
   pinned.  The F-R4 transfer that T1 §2.3 flagged as "structural, not
   literal" is now a measured within-engine fact.
2. **The ε-lift-as-floor direction of F-R4 holds numerically** (15.7–17.8×
   margin here): if the murk sector's soft manifold is SDiff-like, no
   zero-cost (or sub-ε) reconciliation channel exists — escape (γ) stays
   closed, exactly as the audit argued, now with numbers.
3. **Direction-safety of the kill is preserved and sharpened.**  Both
   dispositions still land on KILL: reachability failing reverts the cost
   toward the unscreened scale (the audit's ~1/ε harder); and even granting
   the soft manifold, the measured floor is ≥ 15× the ε-lifted scale, not
   below it.
4. **A new mapping-level fact the corpus text does not anticipate:** in
   this sector the "naive wall" P3 relies on has the WRONG SIGN at every
   overlap — the aligned pair is the costlier configuration (coherent
   sextic overlap), and 94–98% of the mismatch-vs-aligned difference lives
   in exactly the sectors that are frozen on the SDiff orbit.  Within the
   engine, the two-texture overlap energetics is owned by value-function
   densities in BOTH configurations, so an ε-scaled screening estimate has
   no sector to act on.  This is a transfer caveat as much as a finding:
   it may reflect the Skyrme sextic rather than Frank elasticity — but it
   means the engine offers NO support for the corpus's f²p → ε·f²p
   structure at any separation tested.

## 6. Caveats (complete list)

1. **Sector transfer** (§1): mechanism-class test, not a Frank-elastic
   replication; no pitch/chirality; magnitudes not corpus-mapped beyond
   the stated ε and \|dE_naive\| conversions.
2. **Unconverged tails:** both pair runs end at iter_cap = 1200 (tail
   flatness 2.8×10⁻³ align / 7.5×10⁻⁴ repulse; gnorm/gn0 2.0×10⁻² /
   9.2×10⁻³).  cost_full = −0.0755 is therefore resolved only to a few
   ×10⁻² — "degenerate to slightly negative" is the safe reading; its
   sign plays no role in §3–§4 (which compare against gaps of 2.0–2.4).
   The single-knot reference is similarly slow (its seed is already
   near-optimal); the 2×single row is context only.
3. **Guard deviations** from the twoknot per-seed convention (wall
   re-referenced to the continuum B = 2 bound; deg_ref common-min with a
   0.0033 asymmetry) — rationale in §1; endpoint penalties ≤ 1.0×10⁻³
   (wall) and ≤ 2.8×10⁻⁴ (anchor) bound their effect.
4. **SDiff-descent numerics:** the discrete flow is volume-preserving only
   approximately; the enforced budget (|Δ(E6+E0)| ≤ 0.05 ≈ 0.6% of the
   value sectors, |Δdeg| ≤ 0.01) is the error bar of the zero-capture
   claim, and the descent stalled against it — a stronger statement than
   convergence: no legal direction remained.  det∇X distortion max
   7.8×10⁻², mean 2.5×10⁻⁵.  The twist exhibit is leak-free (exact map)
   and independently carries the same verdict.
5. **One mismatch channel run to endpoint** (θ = π about the separation
   axis — the maximal, most SDiff-flattering mismatch); the attract
   channel (π about ẑ) was not campaigned.  The seed scan brackets the
   repulse channel across all separations.
6. **d_eff diagnostic** drifts (1.69 seeded → 2.20/2.31 relaxed): ANF has
   no separation constraint; the relaxed states are the free local minima
   of each basin, which only LOWERS both endpoints symmetrically —
   direction-safe for every comparison above.  Boundary tails ≤ 7×10⁻⁸.

## 7. Files

* `h22_solve/` — tiny crate (Cargo.toml + `src/bin/h22_solve.rs`), builds
  against fs-gum-field/fs-gum-statics unmodified; modes: single, pair,
  sdiff, twist, seedscan.
* `h22_runs/` — per-run JSONs + logs (bit-deterministic given binary+args):
  `h22_single.json`, `h22_pair_{align,repulse}_m09.json`,
  `h22_sdiff_repulse_m09.json`, `h22_twist_repulse_m09.json`,
  `h22_seedscan_m06_m20.json`, `*.log`.
* `h22_analyze.py` — assembles the table above; writes `h22_results.json`.
* `h22_results.json` — the verdict-relevant numbers, machine-readable.

— end of memo —
