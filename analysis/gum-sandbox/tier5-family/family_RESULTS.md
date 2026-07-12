# Tier 5a — Family / Neutrino / Color Sector Audit (RESULTS)

**Campaign:** GUM replication, Tier 5a. **Scope:** every recomputable number downstream of the printed Σ(p) frustration-integral values in Ω-paper §VII.H–J and Appendix K, Course Ch. 18, and the WS-ν/WS-K/WS-N worksheets, plus a print-level consistency audit of the ⟨r7⟩→E-F1→⟨r8⟩ correction chain (Tier-0-gauntlet style).
**Method:** `family_audit.py` (Python 3, deterministic, no RNG; run it to reproduce the full table). Machine-readable output: `family_results.json`.
**Division of labor honored:** the ⟨r10⟩ bridge arithmetic (24.36 ± 0.90; 0.0468 eV; [0.019, 0.115]; the [0.058, 0.11] window; the covariance sign = **F-R3**) was replicated in Tier 0 (RESULTS.md Groups D & I) and is **cited, not redone**; only App K's new printed number (the qR\* parenthetical) is examined here.

---

## Epistemic notice

All verdicts are **within-model**: they test whether the corpus's printed arithmetic is internally consistent and whether its data-side numbers reproduce from PDG inputs. The "≈2σ-equivalent" grade is the **corpus's own** self-grading metric (RT-8), audited here only for arithmetic reproducibility. Nothing in this document bears on whether GUM describes nature. Print-defect candidates are **flagged, not adjudicated** — adjudication belongs to the coordinator.

## Data-side inputs (documented)

PDG/CODATA: m_e = 0.51099895 MeV, m_μ = 105.6583755 MeV, m_τ = 1776.86 MeV (2022/23; 1776.93 for 2024 also checked). Quarks (MS-bar, MeV): u 2.16, d 4.67, s 93.4, c 1270, b 4180; top both ways: 172 500 (direct/MC) and 162 500 (MS-bar m_t(m_t)). M_Pl: 1.2209×10¹⁹ GeV (full) and 2.4353×10¹⁸ (reduced). T_rec = 0.26 eV.

---

## Check table

| # | Check | Printed | Recomputed | Pull / residual | Verdict |
|---|-------|---------|------------|-----------------|---------|
| C1a | ln(m_τ/m_μ) | 2.822 / 2.8222 | 2.82239 (2.82243 w/ PDG24) | — | **PASS** |
| C1b | ln(m_μ/m_e) | 5.332 / 5.3316 | 5.33160 | — | **PASS** |
| C1c | data ratio | 1.889 / 1.8893 | 1.88903 | — | **PASS** |
| C2a | ½A vs ln(m_τ/m_μ) | 2.80 ± 0.45 (0.05σ) | pull +0.050σ | 0.05σ | **PASS** |
| C2b | ½(A+B) vs ln(m_μ/m_e) | 5.65 ± 0.75 (0.4σ) | σ = 0.750 exact (indep.); pull +0.425σ | 0.4σ | **PASS** |
| C2c | "≈2σ-equivalent" self-grade ("prior weight ≈2–4%") | 2–4% | product-of-pull-probabilities 1.3% (≈2.5σ); joint χ²₂ closeness 8.7% (≈1.7σ) | — | **CONVENTION-LIMITED** (self-grade bracketed, plausible) |
| C3a | chain A: 2.9 → ×2 → +14% ⇒ 5.6? | A = 5.6 ± 0.9 | best chain 5.8·1.14 = **6.612**; residual +1.01 (+18%, 1.12σ of printed band); required net correction **−3.4%** vs printed +14 ± 5% | +1.12σ | **PRINT-DEFECT-CANDIDATE (T5-F1)** |
| C3b | chain B: 3.2 → ×2 → (−9+7−4)% ⇒ 5.7? | B = 5.7 ± 1.2 | mult. 5.982 (+0.24σ); additive 6.016 (+0.26σ); required net −10.9% | +0.24σ | **MARGINAL** (within corrections' own ±) |
| C3c | ⟨r8⟩ errors under chain propagation | ±0.9, ±1.2 | chain forces ≥1.37 (A), ≥1.68 (B) — printed errors 30–35% **smaller** | — | **PRINT-DEFECT-CANDIDATE (T5-F1)** |
| C3d | independent-recomputation reading | narrative "with corrections … the ⟨r8⟩ values" | implied ⟨r8⟩ leading: A 2.456 (0.74σ from r7's 2.9), B 3.049 (0.17σ from 3.2) → numbers consistent **only** as re-evaluation | — | consistent-as-recomputation; narrative false as written |
| C4a | r7 ratio (A+B)/A | 2.10 ± 0.35, in N-F1 band [1.6, 2.2] | 2.1034; in-band ✓ (0.097 under edge); indep. σ = 0.385 → printed 0.35 needs ρ ≈ +0.18 (unpublished) | +0.61σ vs data | **PASS** (central/band); error convention-limited |
| C4b | r8 ratio | 2.02 ± 0.28 vs 1.889 (0.5σ) | 2.0179 ± 0.270; pull +0.46σ | 0.5σ ✓ | **PASS** |
| C5a | ΔN_eff = (4/7)(10.75/106.75)^{4/3} | 0.0268; "+0.027" | 0.026772 | — | **PASS** |
| C5b | m_ν ~ ħqc; pitch | "~4 µm" | 1/q = ħc/0.0468 eV = 4.22 µm | — | **PASS** (only checkable number) |
| C5c | B-ν1 = (T·M_Pl·m_ν²)^{1/4} | 1.4 MeV | 1.63 MeV (full M_Pl) / 1.09 (reduced) @ T = 0.26 eV; exact 1.4 needs T = 0.14 eV (full) or 0.71 eV (reduced) | — | **CONVENTION-LIMITED** (WS-ν-P4: "with thermal factors", unprinted) |
| C5d | Majoron g-band endpoints | [8e-10, 1.3e-8] (v2.0); [1.7e-9, 1.3e-8] (v2.0.1) | 0.047/60 = 7.8e-10 ✓; 0.0503/29 = 1.73e-9 ✓; upper 1.3e-8 needs m ≈ 0.050–0.052 (live window), not 0.047 | — | **PASS-with-ambiguity** (m convention mixed across endpoints) |
| C6a | up-type data pair | (1.30, 9.82) | (1.2983, 9.8228) with **m_t = 172.5 GeV direct**; MS-bar top gives (1.31, 9.70) | — | **PASS** (pins convention) |
| C6b | down-type data pair | (0.79, 7.60) | (0.7881, 7.6024); R<1 ⇒ negative effective B — data-side sign flip real | — | **PASS** |
| C6c | theory R from K.4 parts | 1.27 ± 0.09; 0.82 ± 0.16 | rule recovered: R = 1 + B_geo/A − B_tube/A → 1.27 ± 0.089 and 0.82 ± 0.161 (exact incl. error digits) | — | **PASS** |
| C6d | "four numbers at 0.2–0.3σ" | 0.2–0.3σ | +0.31, +0.22, −0.20, −0.17 σ | ✓ | **PASS** (only under direct-top convention; MS-bar top → 0.49σ) |
| C7a | sinθ_c = √(1−1/𝔪), 𝔪 = 1.9 ± 0.4 | 0.69 ± 0.11 | 0.6882; symmetric ±0.080; downside excursion 0.111 → printed band = conservative one-sided | — | **PASS**; NOT the Cabibbo angle (see below) |
| C7b | f = sinθ_c·√ε·𝔪_Sk window & ε-floor | [4, 29] MeV; ε_e ≥ 1.5e-6 | f(1e-5) = 3.70 → "4"; f(6e-4) = 28.7 → "29"; floor (1.4/(0.69·1.7 GeV))² = 1.43e-6 → "1.5e-6" — all with 𝔪_Sk = **1.7 GeV** | — | **PASS** ([4,60] in body superseded in-corpus by F-A15-3) |
| C7c | ε_dress = (f_q/𝔪_Sk)² ≈ 1e-2, "twenty-fold above ceiling" | ~1e-2; 20× | [0.0068, 0.0138]; 11–23× above 6e-4 | — | **PASS** |
| C8a | termination: p=3 at ~100 eV; ε-thresholds | ~100 eV; 8e-6 / 8.4e-6; 4e-5 @1σ-low | ½(A+2B) = 8.50; 104.0 eV; 8.38e-6; 1σ-low (indep. σ = 1.28) → 4.6e-5 vs "4e-5" (one-digit rounding) | — | **PASS** |
| C8b | undressed-ν kill: ε_ν ~ 1e9 ε_e | ~10⁹ | (0.511 MeV/0.05 eV)^{4/3} = 2.2×10⁹ | — | **PASS** (order claim) |
| C8c | bond loop x₀ = 1.90 ± 0.05 vs 1.92 ± 0.08 | closed loop ✓ | pull 0.21σ; **bond equation itself never printed** | 0.21σ | **PASS** (agreement only) |
| C8d | lock shift ln2/0.96 = 0.72 vs 0.50 ± 0.14; pull "1.3σ"→"1.4σ" | 1.3σ (VII.D) / 1.4σ (F-A15-5) | 0.7220 ✓; pull 1.54σ (quad) / 1.59σ (±0.14) | 1.5–1.6σ | **PRINT-DEFECT-CANDIDATE (T5-F3, minor)** |
| C8e | K.2 "qR\* ∼ 10⁻¹⁷" vs own inversion | 10⁻¹⁷ | e^{−24.36} = 2.63×10⁻¹¹; implied R\* = 1.11×10⁻¹⁶ m = τ reduced Compton wavelength exactly | 6 orders / 14.8 ln-units | **PRINT-DEFECT-CANDIDATE (T5-F2)** |
| C8f | bridge graze & band | 0.3σ→0.08σ; ±0.90/ρ | ln(0.0503/0.0468)/0.90 = 0.080 ✓ (one-liner) | — | **CITE-ONLY** (Tier-0 D/I, F-R3; in-corpus C-ν1/F-A15-5) |

Tally: 19 PASS-class · 2 convention-limited · 1 marginal · 4 print-defect-candidate rows (3 distinct candidates) · 1 cite-only.

---

## Finding candidates (flagged for the coordinator; NOT adjudicated)

### T5-F1 (moderate) — the ⟨r7⟩→⟨r8⟩ correction chain does not close for A, and the printed ⟨r8⟩ errors exclude *every* chain reading

The corpus narrates (Ω §VII.H, Course 18.4, App K.1): leading A = 2.9 ± 0.6, B = 3.2 ± 0.9 → exact ×2 (E-F1) → "+14 ± 5% cone stiffening (A); −9 ± 4% rim relief, +7 ± 3% back-reaction, −4 ± 2% inter-face (B)" → **A = 5.6 ± 0.9, B = 5.7 ± 1.2**.

- **A does not close under any of 8 tested readings** (multiplicative / additive, pre- vs post-doubling — these commute, swapped labels, flipped signs, un-doubled corrections, fraction-of-corrected convention). The doubled leading value (5.8) already *exceeds* the printed 5.6 before the +14% is applied; the naive chain gives **6.61 vs 5.6** (residual +18%, 1.12× the printed ±0.9). Closure requires a net **−3.4%** correction — 15.3 percentage points from the printed +14%, more than 3× the printed ±5% correction uncertainty.
- **B nearly closes** (5.98 vs 5.7, 0.24σ; within its own correction uncertainties) — B is not the defect.
- **Errors:** any chain reading forces σ_A ≥ 2×0.6×1.14 = 1.37 and σ_B ≥ 1.68 (the ×2 is exact; corrections only add uncertainty). Printed ±0.9 / ±1.2 are 30–35% *smaller*. Errors cannot shrink under exact rescaling plus uncertain corrections.
- **The only self-consistent reading:** ⟨r8⟩ is an *independent, higher-precision re-evaluation* whose leading values shifted to A_lead ≈ 2.46, B_lead ≈ 3.05 (each within ⟨r7⟩'s stated errors: 0.74σ and 0.17σ — so no cross-record contradiction). Under that reading the printed numbers are coherent, but the narrative sentence "with beyond-leading corrections (+14% …), the ⟨r8⟩ values: A = 5.6, B = 5.7" is **arithmetically false as written**. Same genus as Tier-0's F-R1: a print-level chain that does not reproduce its own endpoint. Note the defect is direction-neutral for the physics claim (5.6 lands *closer* to data than the chain's 6.61 would — ½×6.61 = 3.31 vs 2.822 would be a 1.1σ miss instead of 0.05σ, which is exactly why the chain's non-closure deserves a pinned explanation).

### T5-F2 (minor) — App K.2's "qR\* ∼ 10⁻¹⁷" contradicts its own inversion in the same paragraph

ln(1/qR\*) = 24.36 ⟹ qR\* = e^{−24.36} = 2.63×10⁻¹¹, six orders (14.8 ln-units) from the printed "∼10⁻¹⁷". The inversion side is load-bearing and verified (WS-ν-P2, Tier-0 Group D); its implied R\* = 1.11×10⁻¹⁶ m equals the τ reduced Compton wavelength exactly (consistent with the declared base m̃_τ, since m₃ = m_τ·qR\*). A literal 10⁻¹⁷ would need R\* ≈ 4×10⁻²³ m — matching neither that scale nor the substrate scale a ≲ 10⁻²⁶ m. No downstream number moves (core-term q-blindness needs only qR\* ≪ 1). New App-K item, distinct from F-R3.

### T5-F3 (minor) — the locked-bond pull annotation is unreproducible from printed errors

Predicted shift ln2/0.96 = 0.722 ("0.72" ✓); measured 2.42 − 1.92 = 0.50 ± 0.14 ✓; but the pull from printed errors is 1.54σ (quadrature 0.12⊕0.08) or 1.59σ (±0.14) — vs the body's "1.3σ" and even the corpus's own harmonized "1.4σ" (F-A15-5). Reproducing 1.4σ needs σ_shift ≈ 0.16, larger than any printed combination — presumably an uncertainty on the 0.96 denominator, whose source is never printed. Direction: flatters agreement by ~0.15σ. (F-A15-5 caught the drift but harmonized to a still-unreproducible value.)

**Not new findings (already in-corpus or Tier-0):** the f ∈ [4, 60] MeV body value (superseded by F-A15-3's [4, 29], which we verify closes exactly); the 0.3σ→0.08σ floor-graze annotation (F-A15-5/C-ν1); the ±0.90 band's correlation sign (Tier-0 F-R3).

---

## Positive reconstructions worth recording

1. **⟨r11⟩ quark belts (C6) reconstruct completely** — the sector's cleanest replication. Data pairs (1.30, 9.82) and (0.79, 7.60) reproduce exactly from PDG MS-bar masses with the **direct top mass 172.5 GeV** (MS-bar top fails the R pull); spacing convention = lepton convention (heaviest pair = ½A). The unprinted composition rule is uniquely pinned by the printed numbers: **R = 1 + B_geo/A − B_tube/A** (tube term subtracts; anti-aligned tube 0.51 > geo 0.33 forces the negative effective B — the "derived sign flip" arithmetic is real). Errors reproduce to all printed digits under independence (√(0.08²+0.04²) = 0.089 → 0.09; √(0.08²+0.14²) = 0.161 → 0.16). All four pulls in [0.17, 0.32]σ — "0.2–0.3σ" is fair.
2. **The soft-sector chain closes with one key** (C7): sinθ_c = √(1−1/1.9) = 0.688 → 0.69 ✓; f-window [4, 29] MeV and the ε-floor 1.5×10⁻⁶ both reproduce with 𝔪_Sk = **1.7 GeV**; ε_dress ≈ 10⁻² and the "twenty-fold" margin check.
3. **C7 meaning check (the mission's "Cabibbo 0.225??"):** K.5's "sinθ_c = 0.69 ± 0.11" is the **condensate tilt** (blue-fog cone angle, Thm H-3), *not* the Cabibbo angle. WS-K charter RK-7/H-K1 names the symbol collision as a hazard ("the collision is symbolic, not numeric; s_fog = 0.69 vs sinθ_C = 0.225") and seal q-h forbids any Cabibbo-from-tilt identification. There is no data comparison to check; K.5's retention of the θ_c glyph (against its own H-K1 lock) is a reader hazard, not an arithmetic defect. Adjacent hazard: K.5's "𝔪 = 1.9 ± 0.4" is the dimensionless margin, while the f-formula's 𝔪_Sk = 1.7 GeV is the Skyrme mass floor — two quantities, one glyph family, one appendix line.
4. **Termination arithmetic** (C8a): ~100 eV (104.0), the 8×10⁻⁶ / 8.4×10⁻⁶ threshold (8.38×10⁻⁶), and F-A15-3's 1σ-low 4×10⁻⁵ (4.6×10⁻⁵, one-digit rounding) all reproduce — "exactly three charged families" is genuine downstream arithmetic of (A, B) and the ε-window.
5. **The ≈2σ self-grade** (C2c) is bracketed by the two natural statistics — 1.3% (joint product of two-sided landing probabilities, ≈2.5σ) and 8.7% (χ²₂ closeness, ≈1.7σ) — so the corpus's "2–4%" is achievable but its exact statistic is unprinted. The magnitude-×1.9 pre-E-F1 record is graded separately by the corpus (post-hoc demotion), not part of the 2σ arithmetic.

## What is NOT recoverable from the text

- **The Σ(p) frustration integrals themselves** — the belt integrands, the pancake background, the cone/rim/back-reaction/inter-face correction integrals. Only their quoted outputs (2.9/3.2; 5.6/5.7; the four percentages) are printed. Everything in this audit is downstream of those outputs.
- **The ⟨r8⟩ A–B covariance** — never published (contrast ⟨r10⟩'s ρ = +0.45). The r7 ratio band ±0.35 (needs ρ ≈ +0.18), and the exact 1σ-low termination threshold, are reproducible only up to this unpublished correlation.
- **FQ belt geometry** (K.4): λ\* = 0.25 ± 0.03 and the collapse function c(λ\*); the A_geo/A_tube split (only A_class = 9.4/7.9 outputs are printed; the implied split A_geo ≈ 8.65, A_tube ≈ 0.75 is not independently derivable); quark base scales m̃_t, m̃_b (inputs, per the corpus's own "honest residue").
- **The bond equation** (VII.D/I.2) in closed form, and the source of the 0.96 in ln2/0.96; only loop-consistency (0.21σ) is checkable.
- **γ_soft, χ_soft** (K.3): E ~ γ_soft/q ⟹ m_ν ~ ħqc is dimensionally checkable only via the pitch (4.22 µm ✓); no coefficient is printed.
- **B-ν1's "thermal factors"** (K.5/WS-ν-P4): the raw (T·M_Pl·m_ν²)^{1/4} brackets 1.4 MeV (1.09–1.63 across M_Pl conventions at T_rec = 0.26 eV) but the exact 1.4 needs an unprinted O(1) factor.
- **σ = πf_q²ln κ_q**: κ_q never printed; f_q = 0.14–0.20 GeV ⟺ ln κ_q = 1.5–3.1 — consistent, underdetermined by construction (corpus's own label: "inversion, not derivation").

## Files

- `family_audit.py` — all checks, printed table, reading enumeration for C3 (deterministic).
- `family_results.json` — machine-readable checks, C3 reading table, finding candidates, data inputs.
- `family_RESULTS.md` — this document.
