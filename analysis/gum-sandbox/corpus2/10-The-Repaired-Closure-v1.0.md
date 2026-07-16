# 10 — THE REPAIRED CLOSURE v1.0 (Theorem T3.1)
## The constructive companion, filed as a corpus-style theorem document: the one well-posed ħ-closure the corpus's own equations select, in closed form · Ω v3.0-ext "Appendix M" candidate · 2026-07-16
### Cited wherever text derived from 𝔠₀ = 2.5147 (charter rule 6) · Source memo: theory-audit/T3_repaired_closure.md · Computations: theory-audit/t3_closures.py (deterministic; every closed form ≥ 9 digits; the new exact identities 12–14 digits) · Engine record: tier2-closure/{axi, axi3, axi4_locked, gstar}_RESULTS.md, tier4-field/TIER4_ADJUDICATION.md

**Standing notice (inherited verbatim in force).** Within-model only. "Repair" means exactly: *make the printed mathematics select its printed numbers.* Nothing here bears on nature; no stake is adjudicated by it; adoption is the authors'.

---

## §0 What is being repaired, and the identity that organizes it

The closure's kinematic skeleton survives every external test and is configuration-space-free: E_rot/E = j/2 from L = j𝔠 and the clock alone; the selection algebra w·j(1−j) = 1; the M.0 loop κ = √2ê/𝔠; the ±½ error-signal exponents. What broke is the **specification of the configuration space**: the (V, g)-family presupposes an energy-flat inertia dial that its own SDiff invariance forbids (**F-R4**: 𝕀 = ∫2(q₁²+q₂²) is a value function — invariant to 10⁻¹⁵ under any measure-preserving pullback), and the enlarged space contains a tilt/ring halo channel that makes every printed tuple a saddle (**F-R5**: the fixed-L Routhian decreases under halo growth whenever κ > 1/√2; rung 0.935, benchmark 0.802, deep-BPS 0.764 all sit above threshold).

A well-posed closure is therefore (C-spin) + (C-clock) **plus a declared configuration class 𝒦 and/or stability constraint** on which the fixed-L minimization is bounded and attained. Organizing everything, one new exact identity from the clock fixed point alone (T3.0; verified 10⁻⁹, and to 12 digits in-family where it reproduces κ² = 7/8g):

> **κ² = ê_tot / (2 𝔦_tot) — for any closure solution, on any 𝒦.**

This is the honest generalization of the S4′ stake, and it prices every candidate at a glance.

## §1 Theorem T3.1 (constrained ⟺ saturated; the closure is unique)

**Statement.** On the unrestricted degree-1 class, impose the only self-consistent stability side condition κ ≤ κ_crit = 1/√2 (the infimum of the halo-channel menu κ_crit(w) = 1/√(2⟨sin²θ⟩_w); any weaker channel's constraint leaves the functional unbounded below through equatorial rings — measured, axi4 §5). Then **the stability-constrained closure coincides exactly with the saturated closure**: the unconstrained fixed-L minimum already sits ON the constraint boundary — the KKT multiplier is zero, the constraint is marginally active — and **no intermediate tuple exists between the scaling rung and saturation** for the corpus to land on.

**Proof sketch** (six lines in full at t3_closures.py header; small-η algebra printed in axi3_RESULTS.md §3). In the binding channel (⟨s²⟩ = 1) a marginal far halo adds statics dE₀ = dI/(16π) — the halo ledger is *linear*: E_static = G + I/(16π), G halo-blind. The fixed-L objective R(I) = G + I/(16π) + L²/(2I) is strictly convex in the halo-fed inertia with interior stationarity at 1/(16π) = L²/(2I²) ⟺ **κ = 1/√2 exactly**. Substituting I* and the clock L² = (2/3)·I·E_static gives, for any core, L = √(8π)·G and **𝔠 = (4/π)·G**; minimizing over cores gives **G → G\* = 16√2/9** — G3's sharp Bogomolny-type bound, *attained* by the closed-form **oblate BPS compacton** g(cos(F/2)sinθ) = ρ³/(6√2) with a pinned marginal equatorial halo (both defining identities hold on the G3 tuple to 10⁻¹⁵). The general-channel form 𝔠 = (4/π)√⟨s²⟩·G_w is engine-validated (reproduces axi4's printed uniform-channel pair to 3×10⁻⁷). ∎

## §2 The tuple (all closed forms; the corpus's (4.10) shown for the record)

| quantity | closed form | value | corpus (4.10) |
|---|---|---|---|
| 𝔠_phys | **64√2/(9π)** | 3.201125 | 128√42/(105π) = 2.514754 |
| κ_phys | **1/√2** | 0.707107 | √(7/12) = 0.763763 |
| binding depth 1−κ | 1 − 1/√2 | **29.29%** | 23.61% |
| ω_th = 2κω₀ | **√2·ω₀** | 1.414214 | √(7/3)ω₀ = 1.527525 |
| E_rot/E | ¼ | exact | exact — **survives** |
| g_core · g_tot | **35/24 · 35/12** | 1.4583 · 2.9167 | (g → 3/2 claimed, unreachable) |
| ê_core/𝔦_core · ê_tot/𝔦_tot | **5/4 · 1** | — | 7/4 (spherical backbone only) |
| ê_tot = 𝔦_tot (saturation identity) | **64/(9π)** | 2.263537 | — |
| support-volume ratio | **3/2 exactly** | 1.500000 | V = √2 (restricted family only) |
| mass-law factor E_tot/(Λm̃) | √2ê₀ × 5/(3√2) = 64/(9π) | **+17.85%** | √2ê₀ |
| Λ√J recalibration | ×6√21/35 | 0.785584 | 1 |
| interaction range | κ·ħ/Mc | **0.7071 ƛ_C** (−7.4%) | 0.7638 ƛ_C |
| **S4′ successor stake** | **κ²g_tot = 35/24** (κ²g_core = 35/48; cleanest: ê_tot = 𝔦_tot) | 1.4583 | κ²g = 7/8 |

New exact identities of record (each 12–14 digits in t3_closures.py): (1) **saturation ⟺ ê_tot = 𝔦_tot** — via T3.0, κ = 1/√2 is precisely where the total energy and inertia numbers coincide; (2) **g_core = 35/24** — 97.2% of the old supremum, but with ê_core/𝔦_core = 5/4: *the "7" that generates √(7/8g), √(7/12), and √(7/3) does not survive off the spherical backbone* (the √42 in 𝔠₀ inherits the 7 of 105 = 3·5·7; the saturated geometry replaces it by 1, whence the pure √2 of G\*); (3) **the oblate compacton's support is exactly 3/2× the spherical compacton's** (∫₀^π g(sinθ)/sin²θ dθ = 2; quad-verified 1.3×10⁻¹⁵, not yet symbolically proven — flagged).

**Where the corpus's own measured numbers reappear, lawfully:** its over-spin onset 1.000 ± 0.004 is the polar-channel threshold (⟨sin²θ⟩ = ½ ⟹ exactly 1); its κ₀ = √(7/12) is the sin²-channel threshold (exact 6/7 moment identity); its benchmark 2.37 ± 0.09 sits at the locked uniform-saturated closure 2.37096 (identity numeric; the reading — *the ⟨r1⟩ numbers are thresholds/saturation values of the halo family its solver was blind to* — remains a flagged hypothesis only the frozen code can decide).

## §3 The constrained alternative (retained, per charter rule 6)

**The scaling rung (C1):** 𝔠 = 2√(ê₀𝔦₀) = 256/(15√7π) = **2.053288**, κ = √(7/8), g = 1, V = √2 — well-posed on the declared hedgehog(+dilation) class, replicated to 7 digits. Its honest in-family extension (the corpus's own G.5 spheroidal family, costed with P(λ) ≥ 1): λ\* = 0.8173, g\* = **1.0781**, 𝔠\* = **2.1687**, κ\* = 0.9009 — *the honest g-dial buys 5.6% of 𝔠, not 22.5%*, and **κ²g = 7/8 is an exact identity on the entire locked-transport family** (12 digits) — family-robust, therefore value-free: it cannot discriminate g and cannot adjudicate the closure (the S4′ retarget follows). Both C1 solutions are halo-unstable in every enlarged space; C1 is internally consistent and externally indefensible — printed as the *constrained alternative*, not the repair.

## §4 The caveats (the price of adoption, printed at full size)

1. **Marginal, not gapped, stability.** The halo direction is exactly flat at the solution — a degenerate manifold (core unique; halo pinned only in its ∫η² moments). The zero mode must be sold as a modulus; radiative silence at threshold is marginal, and the consistency ring's third link is repaired only as "no strict descent," not "gapped."
2. **The ε-law flips sign.** Saturated 𝔠 *rises* with ε: 3.2011 → 3.3855 at ε = 0.05 (+5.76%, ≈ +1.15ε class). ⟨r6⟩'s falling 𝔠₀[1 − 0.42ε^{2/3}] law and its geometric-exponent narrative die with the old endpoint; the B-U1′ ceiling's 3/2-exponent derivation dies with them (repaired linear law ⟹ ceiling ≈ 2.6×10⁻³; the admission window [1.5×10⁻⁶, ~3×10⁻³] survives numerically).
3. **The convention question travels with it.** κ = 1/√2 is stated in the closure algebra's own normalization (κ_paper = ω/(√2μ), pinned by the validated unit map); App I.4's control-run sentence must be re-worded to name its channel (s² = |cos|, exact onset 1.00000). A convention error in the halo algebra/unit map would move §§1–2 wholesale — the three attack points are listed (axi3 §5(b)) and every entering number is reproducible from REPRODUCE.md.
4. **The benchmark is re-labelled, not erased.** ⟨r1⟩ becomes the *restricted-branch code-validation point* — a budget-stamped snapshot of the honest flow's transit (F-R15: the (κ, 𝔠) flow passes within the tuple's own bars at cycles 130–170).

## §5 What adopting T3.1 changes (the Group-1 stake/constant register)

- **Ω (4.10) and every echo of it** (Course 12.1/22.3; Primer Ch 10 + NUMBERS TO HOLD; TE Unit 5; audit-note concordance; 00-README stakes line; Session Map §6): 𝔠_phys = 3.2011, κ_phys = 0.7071, binding 29.29%, ω_th = √2ω₀, support ratio 3/2, ħ = 3.201√(ΛJ) = (4/π)G\*·√(ΛJ).
- **The S4′ stake, retargeted** (per S-30, a dated re-stake delta exactly as F-A15-2's): **κ²g_tot = 35/24**; kill re-banded; 0.764 → 0.7071; AUD-15's WS10-S4 ε-extraction assignment moot-for-discrimination.
- **Bench targets** (Ω IX.D/IX.F; Course 23.5; Primer Ch 17): S3 unchanged; V = √2 re-scoped to the restricted branch with support-ratio 3/2 as the saturated diagnostic (⟨r1⟩'s measured 1.409 sits 9σ from 3/2 — the bench must know which number it tests); P8-7 range κ = 0.7071; re-entry protocol channel-resolved via κ_crit(w).
- **The mass law**: total normalization E_tot = 64/(9π)Λm̃ (+17.85%), absorbed by recalibrating Λm̃ (Λ√J ×0.7856); species-universal.
- **Theorem IV.4 + bracket (4.6)**: demoted — supremum of an inertia diagnostic, closure-inert; both bracket endpoints unstable saddles; replaced by the threshold/saturation structure. **App G.5**: struck as mechanism; replaced by the threshold theorem κ_crit(w), Theorem T3.1, G\* = 16√2/9 + the oblate compacton, and the new exact identities. **App B.6/IV.J**: ε-law superseded (sign-flipped). **IV.D**: final table row "saturated (well-posed limit): 3.201, attained"; 2.515 demoted to open-endpoint supremum of non-solutions (C4).
- **Downstream κ/ω_th inheritors** (WS-A8 envelope and well parameter ×0.857; WS-A6 tuple anchor; WS-A census range 2.73×10⁻¹³ m; MIP-2 row 8; the K–K̄ lines below √2ω₀): per their deltas.

## §6 What T3.1 does NOT change (printed at equal size)

**The family sector is 𝔠₀-blind — S1 unaffected either way.** P-O1's argument formalized and verified (T3 §5.2): the mass law's geometric factor is species-universal, so **M_k/M_j = (ê_k m̃_k)/(ê_j m̃_j)** — 𝒢 and Λ cancel; every Tier-5 quantity (Σ(p) structure, ⟨r8⟩ pulls, ⟨r11⟩ belts, the ⟨r10⟩/S1 pipeline) is built from such ratios. **The kinematic pillars survive untouched and externally replicated:** E_rot/E = ¼ exact; the spin-selection theorem and T-B1's algebra (family-scoped); the M.0 loop; the ±½ exponents and −sinΔφ torque; the superradiant rejection of the rigid rung; the clock identities; Corollary IV.3′ (κ = √2ê/𝔠, g-free); the g(λ) *geometry* (F-A15-7's oblate branch, reproduced to 6+ digits — right geometry, wrong energetics). The spectrum, Born statistics, electrodynamics identities, and every sector outside §IV/App G/I.4's variational argument: no contact.

## §7 ⟦H4 resolved⟧ — the box this appendix carried, patched

> **⟦H4 resolved⟧.** The closure's admissibility premise (j ∈ ½ℤ) routes through QUANT/App F.5's consistency part, which consumed the then-unproven minimal hypothesis **MH**. The named computation has executed: **∂(gen) = ±2·[rot]** (two independent routes; 28/28 machine checks). **Completion III is REFUTED** — the fatal outcome this box armed against did not occur: the rotation class survives with order exactly 2, χ(rot) = −1 exists at odd w, **j ∈ ½ℤ stays admissible and this appendix, and the closure it repairs, stand solvable together.** Completion II (ℤ₄) is refuted at odd w. Completion I holds **in corrected form**: ∂ ≠ 0, and the printed "ℤ₂ × A" is corrected to **π₁(𝒞_strat) = ℤ₂ × A/⟨2a_rot⟩** (the winding factor torsioned; the base π₁ absent at w ≠ 0 — F-R16, minor). Residuals, at full size: the N2/F.3 exchange lift stays the GAP; curved-line moduli and even w ≥ 2 stay scoped opens; confidence HIGH in-model, MEDIUM-HIGH on model identification. [H4: theory-audit/h4_completion.md, F-R16]

> **⟦I4 reported⟧.** The residuals have executed: the exchange lift **FAILS** (F-R17 — F.3's own proviso constructs the geon obstruction, σ(exchange) = −1 ≠ +1 = σ(rotation); repair χ_exch = χ(σ)·χ_rot with χ(σ) a new free ℤ₂; "ARE fermions" stays overreach, now with the mechanism exhibited); curved-line moduli close **theorem-grade**; χ(rot) = −1 extends to all w ≠ 0. The closure's admissibility premise is untouched either way — **j ∈ ½ℤ stands on strictly larger ground.** [I4: theory-audit/i4_topology.md, F-R17]

## §7′ ⟦Phase I reported⟧ — archive-side evidence and propagation

**The ⟨r5⟩ corroboration test (I1).** The corpus's own frozen quartic-regime
closure ⟨r5⟩ (𝔠_q = 3.1 ± 0.2) was confronted with both branches computed on
the certified machinery across ε = 0.02–2.0. Verdict: **saddle-excluding and
repair-consistent** — under every reading in which a quartic-regime 𝔠_q is
ε-stable, the archive's number selects the saturated (repaired) branch: the
branches separate by ≈7σ of ⟨r5⟩'s own error bar, and the t-independent
anchor 64√2/9π = 3.2011 sits at **+0.51σ**; no reading reaches it from the
printed (N2-saddle) chain. Two new rigorous floors travel with the test
(𝔠_sat ≥ 64√2/9π at every ε; a Faddeev-type quartic floor), and every
measured value respects them. Normalization caveat printed in the memo.
Filed as archive-side evidence FOR the repair — the archive appears to
contain the repair's fingerprint. [theory-audit/i1_RESULTS.md, i1_results.json]

**Propagation (I2).** Adopting T3.1 end-to-end: every data confrontation
SURVIVES (S1 verified untouched; the ⟨r10⟩ bridge moves ≤ 0.12σ); the bench
sharpens — S4′ becomes discriminating at last (κ²g_tot = 35/24; ⟨r1⟩'s
0.843 ± 0.046 sits 13.4σ from it and 0.7σ from the restricted branch's 7/8:
branch declaration mandatory), and the four-rung over-spin ladder is exact.
Honest costs, printed at full size: the ε-anchor is forfeited; the Q-6′
censorship margin thins ×4.3. Obligations filed: the repaired (4.8′)/(4.9′)
ε-scan; the Majoron battery re-run at g ≈ 7.8×10⁻¹⁰. [theory-audit/i2_RESULTS.md]

## §8 The steelman record (why 2.515 is not adoptable instead)

Every recovery route was constructed and tested (C4): the G\* value-coincidence **disproven exactly** (gap 5.96×10⁻⁴ pinned; the near-miss is 24√21 ≈ 35π); the uniform-channel reference is a non-binding channel accident (twice over); the zero-cost spheroidal reading reproduces 𝔠₀ *exactly* — as the sup of a family of non-solutions, at an open endpoint, on a premise (zero cost) the F-R4 lemma refutes; κ₀ = √(7/12) is recoverable only as a threshold, not a solution; the g = 3/2 algebra does not compose with the true geometry (gives 2.4796, neither number). **One reading stays NOT-AUDITABLE-FROM-TEXT:** the frozen ⟨r1⟩ code may contain an explicit constraint realizing the idealization as an honest penalty — discharge route (a), standing open on the Watch-Mode memo v2's new gate register. Absent that, within the model, 2.515 has no derivation from the printed equations.

---

*Closing sentence, in the corpus's own idiom: the centerpiece was asked to produce its number and could not; asked instead what number its equations do produce, it answered in closed form — **ħ = (64√2/9π)·Λ√J, κ = 1/√2, marginal at its own threshold** — and the repair is filed the way the corpus files everything: with its caveats at full size, its alternative retained, its pending box armed and since patched in its favor (⟦H4 resolved⟧, §7), and its adoption left to the authors.*
