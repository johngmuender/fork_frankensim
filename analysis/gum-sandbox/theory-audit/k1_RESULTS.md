# K1 — the /6 vs /(15π/8) discriminator (J1's residual split): VERIFIED NULL

**ANALYSIS LAYER.** Every claim below is WITHIN-MODEL only: does any printed
archive quantity separate the two surviving edge-referenced amplitude
conventions. Measured facts only; adjudication is the coordinator's.

Reproduce: `python3 k1_discriminator.py` (deterministic, no RNG, ~2 s;
reads only frozen records; writes `k1_results.json` stage-flushed;
re-run verified byte-identical). Console: `k1_run.log`. 17/17 gates PASS.

The candidates (J1, c₂ = 0 benchmark, 2N): 𝔭 = A·e^{−μR*}/6 = **0.843379**
(+0.11σ) vs A·e^{−μR*}/(15π/8) = **0.859059** (+0.64σ). Exact split
48/(15π) = 1.018592: **1.86% in 𝔭, 3.75% in 𝔟** (𝔟 ∝ 𝔭⁻²).

## (a) Enumeration — every archive quantity vs amplitude-convention sensitivity

Sensitivity exponent n (Q ∝ 𝔭ⁿ); power = discrimination in σ of the archive's
own printed precision (instrument error folded where it applies):

| archive quantity | source | n | printed precision | power |
|---|---|---|---|---|
| 𝔭 = 0.84 ± 0.03 | VII.C / I.1 / map §47 / Course 15 | 1 | ±3.6% + 2 dp | **0.52σ** |
| 𝔟_eff = 42 ± 6 | VII.D / I.2 / A.4 | −2 | ±14.3% (+ C_d ±19–41%, ×54 unconverged) | **0.37σ max** |
| x₀(bond eq, 𝔟=42) = 1.90 ± 0.05 | VII.D | via 𝔟 | ±2.6% | **0.075σ** (B3) |
| x₀(measured) = 1.92 ± 0.08 | VII.D | 0 exact | — | 0 |
| 𝔟_∞ ≈ 4.5 | VII.D / Course 16 | −2 | **no band printed** | unusable (B2) |
| B.4 window edges [−m̃²/2, m̃²/6) | B.4 | 0 exact | exact rationals | 0 (B4) |
| C₆ = 64Λm̃/15π | (7.3) / map §47 | 0 | exact | 0 (B5) |
| μ² = m̃²/(2a_ψ) | B.5 | 0 | formula | 0 |
| I.1 gates (virial 3×10⁻⁴, ∫b=1.000, ε=0.05) | I.1 | 0 | thresholds | 0 |
| dipole coeff p ∝ 𝔭/μ² | (7.4)/(15.6) | 1 | **no coefficient printed** | unusable |
| x₀^lock = 2.42 ± 0.12, shift 0.50 ± 0.14 | VII.D (p-2) | 0 (ln-shift; 𝔭² cancels) | — | 0 (+⟨r6⟩ superseded, B.6) |
| C_d (pair-law amplitude) | (7.5) | — | **never printed numerically** | enters only via 𝔟_eff |

The archive prints exactly TWO amplitude-bearing calibrands (𝔭, 𝔟_eff) plus
one unbanded echo (𝔟_∞) — i3's enumeration inputs re-confirmed. Everything
else is shape/energy/potential-space and convention-free by construction.

## (b) Executed discriminators — all fall short of 2σ

**B1 — the 𝔭 band.** Split 0.01568 vs ±0.03 → **0.52σ**. 2σ would need
±0.0078; the archive prints ±0.03. In c₂ terms (slope d𝔭/dc₂ = −0.407 at
benchmark) the switch ≡ Δc₂ = 0.039: 0.84-crossings at c₂ = +0.0083 (/6) vs
+0.0460 (/(15π/8)) — both in-window; the corpus never prints its central c₂.
*Conditional printed-digit channel:* "0.84" (2 dp) has half-ulp window
[0.835, 0.845): 𝔭/6 = 0.8434 rounds to 0.84; 𝔭/(15π/8) = 0.8591 would print
0.86 — a 6.6σ-equivalent exclusion of /(15π/8) IF three unprinted
assumptions hold: (i) benchmark ≡ c₂ = 0; (ii) printed 0.84 = pipeline
central value rounded, not a c₂-scan band center (J1: ±0.03 ↔ c₂ half-width
~0.07, so the scan-summary reading is live); (iii) cross-code amplitude
agreement A_ours = A_frozen ≪ 0.5% (unverifiable — that is the original
OPEN itself). Not archive-backed ⇒ does NOT meet the unconditional bar.

**B2 — the 𝔟_eff channel.** 𝔟 splits by ×1.0375. Max |Δpull| over all
4 C_d estimators × the full named-factor family: **0.37σ** (G4-seed-central4,
F_b = 8π; the in-band G5b/8π² pair gives 0.13σ: 57.0/+1.19 vs 55.0/+1.06).
2σ needs total σ_𝔟 ≤ 1.9%; archive band alone is 14.3%. Two structural
blockers on top: (i) the F_b identity is unpinned (i3 caveat 1) and a
factor relabel absorbs the candidate switch to 0.06σ of the band (pair
(3, π): |Δln − 2ln(48/15π)| = 0.0093); (ii) C_d is unconverged (×54
cap-400→1200, ±19% at best). 𝔟_∞ ≈ 4.5 adds nothing: no printed error, the
far zone is unmeasured by the instrument (archived extrap running-𝔟 at the
largest x are sign-indefinite), and the F_b degeneracy voids even a
generous ±ulp reading.

**B3 — the bond-equation x₀ loop (1.90 ± 0.05).** Rescaling the archived
running-𝔟_eff(x) arrays (seed/relaxed/extrap × F_b ∈ {2, 4π², 8π²}) between
the candidates moves the 42-crossing by at most **Δx₀ = 0.0038 = 0.075σ**
(the crossing sits in the steep rise; i3 caveat 5 quantified). x₀(measured)
= 1.92 ± 0.08 is amplitude-free exactly (well location of E_int(x)).

**B4 — the B.4 window edges re-derive, convention-free (new closure).**
Lower edge: positivity — 𝒱(π) = 2m̃² + 4c₂ ≥ 0 ⇔ c₂ ≥ −m̃²/2 (numeric: min 𝒱
= 0 at the edge, < 0 just below). Upper edge: the quartic Taylor
coefficient of 𝒱 about the vacuum, a₄ = c₂/4 − m̃²/24, **changes sign at
c₂ = m̃²/6 exactly** (numeric root 0.166666667; matches B.4's
"c₂-corrections to the amplitude-mode coupling" language). Both edges are
functionals of 𝒱(f; m̃, c₂) alone — the dimensionless amplitude enters only
the exterior tail readout (VII.C) — so they re-derive bit-identically under
either candidate: discrimination power exactly 0.

**B5 — the amplitude-free anchors re-verified.** C₆ Haar quadrature =
64/(15π) to 1.6×10⁻¹³ (the 15π here is the energy-bound constant, NOT the
tail normalization — numerological proximity to 15π/8 carries no
constraint); μ = m̃/√(2a_ψ) convention-free at ≤ 10⁻⁴ across the window
(J1, all 9 points).

**B6 — the two-knot amplitude cannot be re-fit from the archive.** The
frozen fs-gum-twoknot records hold scalar diagnostics + iteration series
only (largest array 41 entries vs 1.47M cells; boundary_tail is a scalar):
no field snapshot exists, so a direct 3-D tail-amplitude fit is impossible
from the records, and a fresh 3-D re-run would add OUR precision, not
ARCHIVE precision — the bottleneck is the printed ±0.03/±6 bands.

## (c) Gate table

17/17 PASS (`k1_run.log`, booleans in `k1_results.json/gates`): B1 input
integrity vs j1_run.log; B1 power < 2σ; B1 rounding in/out; B2 max Δpull
< 2σ; B2 F_b-relabel absorption; B2 𝔟_∞ blockers; B3 power < 2σ; B4 lower
edge exact + sharp; B4 upper edge = 1/6 to 10⁻⁹ + analytic a₄ match; B4
amplitude-free (structural); B5 C₆ + μ; B6 no-snapshot; C verified-null.

## Verdict (analysis layer — measured facts; adjudication is the coordinator's)

**VERIFIED NULL — the /6 vs /(15π/8) split is INDISTINGUISHABLE WITHIN THE
ARCHIVE.** The archive prints two amplitude-bearing calibrands; the
strongest unconditional discriminator is the 𝔭 band itself at **0.52σ**
(𝔟_eff: 0.37σ max; x₀ loop: 0.075σ; every other quantity has sensitivity
exactly 0 or no printed digits). Nothing reaches the 2σ bar and no
combination can (the channels share the same two bands). The item closes
as "indistinguishable within the archive" — J1's edge-referenced CLASS
recovery stands as the terminal within-model statement. One conditional
pointer is printed, not claimed: the second printed digit of "0.84" is
consistent with /6 (0.8434) and not with /(15π/8) (0.8591 → "0.86"),
a 6.6σ-equivalent split that rests on three unprinted assumptions (B1).

## Defects (printed per protocol)

1. **j1_results.json on disk is a partial stage-flush** (4/9 grid points,
   no summary block) inconsistent with j1_run.log/j1_RESULTS.md — K1
   restricted grid arithmetic to archived points (both 0.84-crossings lie
   in the surviving c₂ ≥ 0 range) and gated benchmark values against the
   log. The J1 record should be regenerated (53 s re-run).
2. **Prior split figure overstated:** J1/ROADMAP_v7 carried "/6 vs /(15π/8)
   0.7σ apart"; the measured split is 0.52σ (0.0157/0.03). The null verdict
   is *stronger* than advertised; no prior conclusion changes.
3. The rounding-channel σ-equivalent (6.6σ) treats the half-ulp as a
   uniform error (σ = 0.005/√3); it is a conditional illustration, not a
   calibrated statistic, and is excluded from the verdict gate.
4. The B4 upper-edge re-derivation identifies c₂ = m̃²/6 with the quartic
   a₄ sign change; the corpus states the window without derivation
   ("positivity + saturability"), so the identification is our reading —
   exact numerically, but the corpus's own reason is not printed.
5. B2's Δpull maximum occurs at out-of-band rows (pull ≈ +3); among
   jointly in-band solutions the split is smaller still (0.13σ). Either
   way < 2σ; the maximum is quoted to be conservative toward the null.
6. Everything is WITHIN-MODEL bookkeeping of a frozen pipeline; nothing
   validates physics/nature.

## Files

```
theory-audit/
  k1_discriminator.py   enumeration + all channels + gates (~2 s, deterministic)
  k1_results.json       stage-flushed records (catalog, B1-B6, gates, verdict)
  k1_run.log            console transcript incl. full gate table
  k1_RESULTS.md         this file
```
