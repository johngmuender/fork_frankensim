# J4 — multi-ε saturated re-scan: the (4.8′) exponent measured (analysis layer)

**Phase:** ROADMAP_v6 J4 · 2026-07-16 · theory-audit numerics subagent (Fable-class).
**Charge:** J2 defect D3 — the repaired law 𝔠_sat(ε) = 3.2011·(1 + 1.152ε) is a LINEAR
class fitted on TWO measured saturated points. Densify the small-ε saturated branch,
measure the exponent, and re-propagate (4.9′) per surviving reading.
**Computations:** `j4_epslaw.py` (deterministic, i1 machinery reused unmodified) →
`j4_results.json` (stage-flushed). **WITHIN-MODEL ONLY; the coordinator adjudicates.**

## §0 Provenance: which two points J2 used, and how i1's curve relates

**J2's linear fit used exactly two points, both from `tier2-closure/axi4_locked_results.json`
(ring-saturated engine):** `saturated.ring_t0` 𝔠 = 3.203459 (t → 0; exact 𝔠_sat = 3.2011247,
engine bias +7.3×10⁻⁴ = J2's D4) and `saturated.ring_eps` 𝔠 = 3.385470 (ε = 0.05) →
s_exact = 1.1518, s_engine = 1.1359 (G0b/G0c reproduce both).
**i1's curve is a different family, same objective:** the 12-mode axi3-basis saturated
closure G_t on the ε-dial hedgehog base (i1_r5test.py `saturated_at`), measured at
ε ∈ {0.02, …, 2.0}. At ε = 0.05: i1 gives 3.382843 — 7.8×10⁻⁴ rel BELOW axi4's ring value,
i.e. the tighter of the two upper bounds on the true G_t*; its exact-base slope reading at
0.05 is 1.135 vs axi4's 1.152 (≈1.5% family-amplitude spread, same class as J2's D4).
Consequence for J4: **within-family consistency (i1 machinery reused unmodified at every
new ε) decides the EXPONENT; the axi4-vs-i1 family spread is an AMPLITUDE band** carried
through §3 as separate readings. Both families are variational upper bounds; the exact
ε → 0 anchor 64√2/9π is the only point known exactly.

## §1 The densified small-ε dataset (MEASURED)

STATUS: PENDING (stage 1)

## §2 Exponent fits (window, classes, exclusions)

STATUS: PENDING (stage 2)

## §3 (4.9′) propagation per reading (J2/J3 formulas verbatim)

STATUS: PENDING (stage 3)

## §4 Cross-checks and gates

STATUS: PENDING (stage 4)

## §5 Defects (printed honestly)

STATUS: PENDING

## §6 Verdict (analysis layer)

STATUS: PENDING
