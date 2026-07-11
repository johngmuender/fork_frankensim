# Tier 2a — The Reduced Closure Loop: REPORT & ADJUDICATION

**Status: COMPLETE (reduced model). Tier 2b (full 2-D field solve) remains
future work — and this tier now quantifies exactly what it must supply.**

## What was run
`closure_reduced.py`: the paper's own reduced structure for the isorotating
knot (Sec. IV.C/IV.J; App. G.5 oblate branch as corrected by AUD-15
F-A15-7) — shape moduli entering only through `g(λ)·V`, `E_field` on the BPS
locus, spin + clock constraints — plus an ε-sector penalty `ε·c_p/λ^p`
(both `p = 1`, the form consistent with the corpus's geometric exponent,
and `p = 4/3`, the naive slab-gradient estimate). **Exactly one constant
(c_p) calibrated from exactly one datum (λ* = 0.42 at ε = 0.05); every
other quantity predicted.**

## Results

**Anchors (must-pass):**
- ε→0 limit reproduces the exact family to 6 decimals: 𝔠 = 2.514753,
  V = √2, κ = √(7/12), E_rot/E = ¼. PASS.
- Clock residual 0 to machine precision; **E_rot/E = ¼ automatic at every ε**
  — confirming the paper's Cor. IV.2 claim that ¼ is pure kinematics of the
  spin+clock constraints, independent of the potential's composition.

**Predictions vs the ⟨r1⟩ benchmark (ε = 0.05; p = 1 form):**

| Quantity | Reduced model | ⟨r1⟩ | Pull |
|----------|--------------|------|------|
| g* | 1.2855 | 1.31 ± 0.04 | −0.6σ ✓ |
| 𝔠 | 2.4478 | 2.37 ± 0.09 | +0.9σ ✓ |
| κ | 0.8454 | 0.802 ± 0.018 | **+2.4σ** |
| V | 1.4510 | 1.409 ± 0.010 | **+4.2σ** |
| E_rot/E | 0.250000 | 0.2500 ± 0.0002 | 0 ✓ |

(p = 4/3 form: 𝔠 = 2.4162 (+0.5σ), κ = 0.8402 (+2.1σ), V = 1.4412 (+3.2σ).)

**ε-scan law** (`epsscan_p*.csv`): fitted deficit exponents 0.749 (p = 1)
and 0.589 (p = 4/3) over ε ∈ [10⁻³, 10⁻¹], vs the corpus's 2/3 with
acknowledged O(ε^{1/3} ln ε) corrections. The p = 1 form is the one
consistent with 2/3-plus-logs (a positive log correction raises the
finite-window fit, as observed); the p = 4/3 form gives the analytically
expected 3/5 — my own scaling analysis confirmed in the numerics. Local
c_g at ε = 0.05 (exponent pinned to 2/3): 0.195 (p = 1) vs corpus
0.42 ± 0.04 — **~2× low**.

## Adjudication

1. **The loop's coarse structure replicates**: with one calibrated constant,
   the reduced model lands g* and 𝔠 inside the benchmark's 1σ band, exactly
   reproduces the ε→0 endpoint and the ¼ invariant, and produces the
   ε^{2/3}-class scan law from pure spheroid geometry. The paper's claim
   that the closure's physics is "shape entering only through g(λ)V"
   is *structurally* correct at the few-percent level in 𝔠.
2. **The tuple's fine structure does NOT replicate in the reduced model**:
   κ (+2.4σ), V (+4.2σ), and c_g (~2× low) all miss. The pattern is
   informative: my V-independent penalty *raises* V above √2, while the
   benchmark measured V = 1.409 slightly *below* √2, and the corpus quotes
   a benchmark binding fraction δ = 0.29 ± 0.05 (the ε-sector *binds* and
   couples to dilation). **Conclusion: the ⟨r1⟩ tuple's fine structure
   genuinely requires the field-level ε-sector (V-coupled, binding), not
   any one-constant shape penalty.** This is now a quantified statement of
   what the full Tier-2b solve must supply — and equally, a quantified
   falsification target: if Tier 2b's field solve also failed to pull
   κ → 0.802 and V → 1.409, the benchmark tuple would be in trouble.
3. **FINDING F-R2 (minor):** AUD-15 §V15.2 quotes "g(0.37–0.47) spans
   1.28–1.32" for its own corrected oblate formula; direct evaluation gives
   **[1.257, 1.315]** (g(0.47) = 1.2567, not ≈1.28). The verification's
   conclusion survives (overlap with g* = 1.31 ± 0.04 holds) but the quoted
   span is inconsistent with the formula it cites — same class as F-R1:
   a print-level arithmetic slip, conclusion unaffected.
4. **Not a falsification of ⟨r1⟩**: the reduced model is *my* approximation;
   the corpus's benchmark came from a field-level solve. The misses bound
   the reduced-model gap, not the corpus's correctness.

## What Tier 2b must now do (sharpened)
Solve the axisymmetric field problem with the actual ε-sector densities
(W₂ + W_χ + W₄ evaluated on the deformed compacton, with the log-divergent
edge boundary layer resolved) and check: κ → 0.802 ± 0.018,
V → 1.409 ± 0.010 (below √2), c_g → 0.42 ± 0.04, binding δ → 0.29 ± 0.05,
over-spin onset κ = 1.000 ± 0.004. The reduced model's misses are the
sensitivity map: these five numbers carry the field-level information.

## Epistemic notice
Within-model replication of a speculative theory's internal consistency;
nothing here bears on nature.
