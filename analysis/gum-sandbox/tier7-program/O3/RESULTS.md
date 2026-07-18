# O3 — F-R14 rescue decision matrix [F-T7-O3]

**Phase:** Tier 7 program, Phase O addendum (ROADMAP_v9_PROGRAM.md item O3).
**Date:** 2026-07-18.  **Code:** `o3_matrix.py` (exact rational arithmetic
for every matrix cell, float scan for the window figure; deterministic,
21/21 internal checks PASS, ~2 s).  **Outputs:** `o3_results.json`,
`o3_fig.png`.
**Sources (inputs, independently re-implemented and cross-checked here):**
`theory-audit/h33_RESULTS.md` Part B + §8 and `h33_pq.py` (master form,
factor dictionary, anchors, rescue window); `theory-audit/h25_RESULTS.md`
§§3, 5 (sector N_s; knot rows K2/K3); `corpus3/01-GUM-Omega-Paper-v4.1-ext.md`
(F-R14 dispute box at VI.C; App. E.6 rewrite; the "(6.3) margin 4.0,
boundary-exact" arithmetic).

**Epistemic frame (binding).  Within-model; nothing here bears on nature.**
This workstream is **DECISION SUPPORT ONLY**: the corpus's one unrepaired
structural defect (F-R14) is turned into a fully computed matrix — every
textually available reading priced, **no reading recommended**.  The
constitutive choice is the authors' (charter language).  The F-R14 regrade
language is **unchanged** by anything below.

---

## 1. What is implemented

- **Master form** (h33, re-implemented in exact `Fraction` arithmetic):
  𝔞₁(p,q) = −(5p+q)/(12(2p+q)) per scalar mode; positive iff
  (5p+q)(2p+q) < 0.
- **Factor dictionary** (exponents of (1+h), additive): √g₃ → +3/2; each
  derivative-index contraction → −1; field pair covector/frame/vector →
  −1/0/+1; group-trace kinetic → 0 (Killing form, no spatial index).
- **Structural invariants** (re-verified exactly): B2 photon doublet
  **p = q − 1** under every variance reading and every density weight;
  Q̃-sectors **q = 3/2 group-protected**.
- **Knot band-edge (Dirac-type)**: 𝔞₁ = 4(1/6 − 1/4) = **−1/3 per field
  before the loop sign** (h25 K2); supertrace-signed loop (h25 K3,
  Γ_F = −Tr ln) = +2 minimal-scalar equivalents = **+1/3 per field**.
  Every matrix row is shown **both ways**.
- **P-acoustic family**: density weight w shifts (p,q) by (+3w,+3w) off the
  frame baseline; rescue window w ∈ (−5/18, −2/9), width 1/18;
  distinguished point w = −1/4 on q = −3p.
- **Sector content for Σ_s w_s = 1/16πG** (w_s = N_s Λ_s²𝔞₁): N_s = 2, 2, 1
  for B2±, B3, B4 (printed via the corpus's own operator, h25 §3 table);
  **N_knot is NOT printed anywhere in the corpus → parametrized (default 3),
  with every threshold reported as a function of N_knot**; no per-sector Λ_s
  is printed → equal Λ assumed, all weights quoted in units of Λ².

## 2. Gate G1 — the five h33 anchors, reproduced exactly

| anchor | (p,q) computed | 𝔞₁ computed | h33 printed | match |
|---|---|---|---|---|
| covector, B2 | (−1/2, 1/2) | −1/3 | −1/3 | EXACT |
| frame/soldered, B2 | (1/2, 3/2) | −2/15 | −2/15 | EXACT |
| vector, B2 | (3/2, 5/2) | −5/33 | −5/33 | EXACT |
| cone-only | (0, 1) | −1/12 | −1/12 | EXACT |
| P-acoustic w = −1/4 | (−1/4, 3/4), on q = −3p | +1/6 | +1/6 | EXACT |

## 3. THE DECISION MATRIX (centerpiece — every cell computed, none asserted)

Units: weights /Λ² (equal Λ_s); Σ_bos = 2𝔞₁^(B2) + 3𝔞₁^(B3/B4);
totals at the default N_knot = 3.  "Flip thr." = the exact N_knot the
supertrace-signed knot count must **exceed** for the G-sign to flip
positive (equal-Λ).  (6.3) margin arithmetic in all cells:
log₁₀(10⁻¹⁵/10⁻¹⁹) = **4.0 orders, boundary-exact** — corpus language
"IF the bound held" retained.

| reading | 𝔞₁^(B2) | 𝔞₁^(B3/B4) | Σ_bos | knots (K2 naive / K3 supertrace) | Σw_s (naive / super) | G-sign (naive / super; flip thr.) | hull status (naive / super) | (6.3) margin status | priced cost |
|---|---|---|---|---|---|---|---|---|---|
| **frame / soldered** (favored) | −2/15 | −2/15 | −2/3 | −1 / +1 | **−5/3 / +1/3** | − / + (thr. N_knot > 2) | formally convex (all −) / **MIXED — c_GW² can exit** | naive: ratio arithmetic retained (4.0 orders IF bound held) but **1/16πG < 0 voids sector upstream**; super: demoted to order-of-magnitude | favored: F10′ objectivity + KBKK soldering, zero new postulates; outcome cost: inverted G (naive) or hull exit (super) |
| **covector** (unsoldered) | −1/3 | +1/6 | −1/6 | −1 / +1 | **−7/6 / +5/6** | − / + (thr. N_knot > 1/2) | **MIXED / MIXED** — c_GW² can exit either way | demoted to order-of-magnitude both ways (4.0 orders IF bound held) | violates F10′ frame-relative discipline (the corpus's cured M-1); abandons KBKK soldering; unmodeled soldering-gradient connection terms (h33 §9) make its +1/6 the least stable cell |
| **vector** | −5/33 | −1/6 | −53/66 | −1 / +1 | **−119/66 / +13/66** | − / + (thr. N_knot > 53/22) | formally convex (all −) / **MIXED** | as frame row | no printed text selects the vector variance; same F10′ tension without the +1/6 payoff |
| **cone-only** (VI.B displayed g_s) | −1/12 | −1/12 | −5/12 | −1 / +1 | **−17/12 / +7/12** | − / + (thr. N_knot > 5/4) | formally convex (all −) / **MIXED** | as frame row | reads the displayed time-time g_s literally; ignores that spatial covariantization is forced by VI.A (strain-only self-destructs, F-H33-3) |
| **P-acoustic w = −1/4** (+ frame) | +1/6 | +1/6 | +5/6 | −1 / +1 | **−1/6 / +11/6** | − / + | MIXED / **CONVEX: all w_s > 0, (6.2) holds as inequality** | naive: still fails (knots −); **super: (6.3) restored as inequality, margin 4.0 orders boundary-exact** | NEW [CJ]-class quantization-measure postulate, not an erratum — three printed costs (§4) + tuning fraction **1/4** (§4) |

**The one fully sound cell is P-acoustic + supertrace-signed knots** — i.e.
exactly h33's "two independent repairs is the minimum bill for one appendix
line."  Every other cell fails computationally in one of two ways: inverted
Newton constant (all-negative weights: Σw_s < 0 while the ratio (6.2)
formally survives), or hull exit (mixed signs: c_GW² not a convex
combination, the slaving bound (6.2)/(6.3) fails as an inequality).

## 4. P-acoustic pricing (Gate G3)

**Window:** w ∈ (−5/18, −2/9), width **1/18** — verified as the exact wedge
of the frame baseline + 3w: the w = −2/9 edge is the 5p+q = 0 zero
(𝔞₁ → 0), the w = −5/18 edge is the **2p+q = 0 pole** (𝔞₁ diverges — h25's
"the (6.1) ansatz fails structurally" ray is the window's own edge).

**Tuning fraction (quantified):** both natural conventions are excluded and
**equidistant** — d(w=0) = d(w=−1/2) = 2/9 from the window as a set (= 1/4
from its center, which is exactly the distinguished point w = −1/4):

- primary (set-distance): (1/18)/(2/9) = **1/4 = 0.25**
- center-distance variant: (1/18)/(1/4) = **2/9 ≈ 0.222**

**Window scan (61 interior points; figure `o3_fig.png`):** 𝔞₁(w) > 0
throughout (min 0.0027 at the zero edge, diverging toward the pole edge);
per-sector weights, Σw_s and hull status computed at every point:

- **supertrace bookkeeping:** Σw_s > 0 and all-positive **at every point of
  the window** — convex, G > 0.
- **naive bookkeeping:** Σw_s > 0 only on the pole-side subwindow
  w ∈ (−5/18, −25/99) — exact crossover **w = −25/99** (where
  5𝔞₁ = N_knot/3 = 1/5), subwindow width **5/198 = 5/11 of the window** —
  and even there the hull is **MIXED** (knots negative): P-acoustic alone
  never yields a sound cell.  (Scan fraction 0.459 vs exact 5/11 = 0.4545.)

**Three printed costs** (carried verbatim from h33 §8 / the E.6 rewrite):
(1) reopens Sec. III's flat-measure [DF] derivations (Fisher/Madelung,
Wallstrom, Born rates) on curved backgrounds; (2) ρ₀, J ∝ 1/det g₃ against
VI.A's own KBKK wedge-insertion picture and the F5/F8 bookkeeping;
(3) unforced: nothing printed selects the window (tuning by VI.D's own
standard) — now quantified by the tuning fraction above.

## 5. The supertrace erratum, both ways (Gate G3)

Knot row shown both ways in every reading (§3).  **Does the G-sign
conclusion flip?**  Computed answer: **naive (K2, formula as printed):
G-sign is negative in all five rows** — including P-acoustic at w = −1/4.
**Supertrace-signed (K3): at the default N_knot = 3 (equal Λ) the G-sign
flips to positive in every row**; the exact flip thresholds are
N_knot > 2 (frame), > 1/2 (covector), > 53/22 (vector), > 5/4 (cone-only).
**But the flip never buys soundness by itself:** in the four geometric
readings the supertrace sign creates mixed signs (bosonic −, knots +) —
the failure mode *changes* from inverted-G to hull-exit; it does not go
away.  The G-sign conclusion is therefore bookkeeping- **and
multiplicity-dependent** (N_knot, Λ_s are unprinted), while the
**convexity failure outside the P-acoustic + supertrace cell is not** —
that disjunction is the robust content, matching h33's exhaustive
disjunction exactly.

## 6. Gates

| gate | requirement | measured | verdict |
|---|---|---|---|
| G1 | reproduce h33's five anchors EXACTLY | −1/3, −2/15, −5/33, −1/12, +1/6 — all exact (Fraction equality, §2) | **PASS** |
| G2 | matrix fully computed, no cell asserted | 5 rows × 2 bookkeepings, ~104 computed cells; every status derived from computed sign patterns | **PASS** |
| G3 | P-acoustic pricing quantified; supertrace both ways | tuning fraction 1/4 (set) / 2/9 (center); window scan 61 pts; flip thresholds exact per reading; naive-G>0 subwindow 5/198 | **PASS** |
| G4 | honest bottom line printed | §7 below; F-R14 regrade language unchanged; no reading recommended | **PASS** |

## 7. The honest bottom line (Gate G4)

**DECISION SUPPORT ONLY.  The constitutive choice is the authors'
(charter language).  No reading is recommended here; each is priced.**
The F-R14 regrade language is unchanged by this workstream: Theorem VI.1
remains OPEN-DEFECT (unconditional as printed); (6.3) remains demoted to
an order-of-magnitude estimate (the verified arithmetic — margin 4.0
orders, boundary-exact — holds only IF the bound held); P-M3 remains
premise-suspended pending the E.6 discharge.  Adoption of P-acoustic plus
the knot supertrace erratum — the only computationally sound cell — would
be a new [CJ]-class constitutive postulate plus an erratum against the
printed weight formula: **two independent repairs is the minimum bill**,
priced above at tuning fraction 1/4 and three structural costs.

## 8. Key numbers

| quantity | value |
|---|---|
| anchors reproduced | −1/3, −2/15, −5/33, −1/12, +1/6 (all exact) |
| Σ_bos /Λ² by reading | −2/3 (frame), −1/6 (covector), −53/66 (vector), −5/12 (cone-only), +5/6 (P-acoustic w=−1/4) |
| knot row /Λ² (N_knot = 3) | −1 (naive K2) / +1 (supertrace K3) |
| Σw_s /Λ² naive / super | frame −5/3 / +1/3; covector −7/6 / +5/6; vector −119/66 / +13/66; cone-only −17/12 / +7/12; P-acoustic −1/6 / +11/6 |
| G-flip thresholds (N_knot >) | 2, 1/2, 53/22, 5/4 (frame, covector, vector, cone-only) |
| window; width | (−5/18, −2/9); 1/18 |
| tuning fraction | **1/4** (set-distance; both conventions equidistant at 2/9) / 2/9 (center) |
| naive-G>0 subwindow | (−5/18, −25/99), width 5/198 = 5/11 of window; hull still MIXED there |
| (6.3) margin | 4.0 orders, boundary-exact (IF the bound held) |
| sound cells in the matrix | exactly 1 of 10: P-acoustic + supertrace |

## 9. Method (one paragraph)

Exact rational arithmetic throughout the matrix: the factor dictionary
(exponents additive in log(1+h)) reproduces h33's sector (p,q) per reading;
the master form is evaluated as a `Fraction`; sector sums use the h25
operator-table multiplicities (2, 2, 1) and a parametrized N_knot; hull and
G-sign statuses are derived from the computed sign patterns (all-positive →
convex; all-negative → formally convex with Σ < 0; mixed → hull exit), not
asserted.  The window scan evaluates 𝔞₁(w) = −(4+18w)/(30+108w) (frame
baseline + 3w, closed form) at 61 interior points; the crossover −25/99 and
subwindow width 5/198 are solved exactly and confirmed by the scan.

## 10. Caveats (honest limits)

1. **N_knot and Λ_s are unprinted.** All Σw_s totals and G-flip thresholds
   assume equal Λ_s and default N_knot = 3; the thresholds are exact in
   N_knot, but any Λ-hierarchy reweights them.  The hull/convexity statuses
   per cell do **not** depend on N_s, Λ_s (they depend only on sign
   patterns) — they are the parametrization-robust column.
2. **The knot ±1/3 is the minimal/Lichnerowicz anchor** (h25 K2/K3,
   torsion-free import).  Its response across the four geometric readings
   and across the P-acoustic window was not derived in h25/h33 and is
   carried as constant here, flagged; h25 states K2's w_knot < 0 holds in
   every geometric reading, which is what the naive column uses.
3. **Linear order in h, conformal probe class** — inherited from h33; the
   covector +1/6 additionally carries h33 §9's unmodeled connection-term
   instability (noted in its cost cell).
4. The (6.3) margin column re-states the corpus's own verified arithmetic
   (4.0 orders) under each reading's computed bound status; it introduces
   no new bound.
5. This memo prices; it does not adjudicate.  Promotion or adoption is the
   coordinator's / authors' call.

## Deliverables

- `o3_matrix.py` — the computation (21/21 checks).
- `o3_results.json` — machine-readable: anchors, full matrix, scan arrays,
  pricing, gates, check list.
- `o3_fig.png` — weights vs w across the P-acoustic window (panel A: the
  window in context with both natural conventions outside; panel B:
  per-sector weights and both totals, crossover marked).
- `RESULTS.md` — this memo.  Filed as **F-T7-O3**.
