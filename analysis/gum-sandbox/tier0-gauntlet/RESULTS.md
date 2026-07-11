# Tier 0 — The Constants Gauntlet: RESULTS

**Status: COMPLETE. 37/37 checks PASS.**
Run: `cargo run --release` (frankensim pinned nightly; std-only dependency
closure `fs-ivl → fs-math`, `fs-evidence → fs-obs`, `fs-package → fs-blake3 +
fs-crosswalk`, `fs-checker`). Full output: `RESULTS.txt`.

## What this is
The first tier of the GUM replication program (per
`../GUM_SIMULATION_ASSESSMENT_v2_FABLE.md`): every exact constant and identity
of GUM-Ω v2.0.1 Sections IV/VII — the AUD-15 §V15.1–V15.9 verification sweep —
**re-executed independently, outside the corpus**, on frankensim's certified
arithmetic. This is the external re-derivation the corpus's own governance
(SWP Branch S) mandates and could not perform on itself.

## Method
- **fs-ivl outward-rounded intervals** for every computed quantity (π enclosed
  by `[f64::PI, next_up]`; quadratures by interval Riemann sums with rigorous
  O(1/n) enclosures at n = 200,000).
- **Exact integer/rational cross-multiplication** wherever the corpus claims
  exactness (7/4, 7/3, 8/35, 2/15, the 𝔠₀² square identity, w·j(1−j)=1, the
  Yukawa-trace polynomial cancellation).
- **fs-evidence `Certified<f64>`** — every passing enclosure re-validated by
  the fail-closed certifier.
- **fs-package + fs-checker** — all 37 claims bundled as Verified claims with
  content-addressed BLAKE3 certificates (`hash_domain("gum-tier0:cert:v1",
  id|statement|lo_bits|hi_bits)`), and checked three ways:
  | Mode | Capability | Result | Meaning |
  |------|-----------|--------|---------|
  | 1 | deny-all | `passed = false` | unauthenticated Verified claims refused — **the anti-laundering rule demonstrated** |
  | 2 | recompute-and-compare certificate verifier | `passed = true` | every certificate re-derived from the claim's own typed content |
  | 3 | tampered root | `passed = false` | tamper detection |

## Headline verifications (37 checks; selection)

| ID | Claim | Enclosure (midpoint) | Corpus | Verdict |
|----|-------|---------------------|--------|---------|
| A3 | ê₀/𝔦₀ = 7/4 **exact** | 1.75 (+ integer identity) | 7/4 | PASS |
| A4 | 𝔠₀ = 128√42/105π **two routes + exact square identity** | 2.514753626 | 2.5147 | PASS |
| A5 | κ = √(7/12); ω_th = √(7/3) | 0.763762616; 1.527525232 | 0.7638; 1.5275 | PASS |
| A6 | rigid rung √(7/6) **> 1 certified** (superradiant) | 1.080123450 | 1.0801 | PASS |
| A8 | j-family: V = √2, ¼ exact; j=1 pole; j=3/2 → V²<0 | — | — | PASS |
| A9 | 𝔠(0.05) band ∩ ⟨r1⟩ 2.37±0.09 | [2.358, 2.385] | 2.37 | PASS |
| B1+B2 | **two computations, one number**: Haar quadrature ∩ compacton direct integral, both = 64/15π | width 6×10⁻⁵ | 1.3581 | PASS |
| B3 | hedgehog degree K = 1 (π₃(S³)) | 1.0 ± 1×10⁻⁵ | 1 | PASS |
| D1–D3 | lepton logs 2.8222/5.3325/1.8894; pulls 0.05σ/0.42σ/0.47σ | — | 0.05/0.4/0.5σ | PASS |
| D5 | ⟨r10⟩ bridge: 24.357 → m₃ = 0.0468 eV | [0.04694, 0.04695] | 0.0468 | PASS |
| D6 | minimal-NO Σ = 0.0590 eV (stake floor 0.058) | 0.05896 | 0.059 | PASS |
| D7 | ΔN_eff = 0.02677 | — | 0.0268 | PASS |
| E1 | g_s=2 chain: ½qωƛ² ≡ μ_B (identity + CODATA) | 9.2740101×10⁻²⁴ | 9.274×10⁻²⁴ | PASS |
| E2 | Unruh T(9.81) = 3.978×10⁻²⁰ K | — | 4×10⁻²⁰ | PASS |

## Findings of note (the honest ledger)
1. **No arithmetic defect found in the corpus.** Every quoted constant in the
   swept set is reproducible within quoting precision; the exactness claims
   (7/4, ¼, √2, 7/3, 8/35, the 𝔠₀ square identity) are exact as claimed.
2. **One harness defect (mine), caught and fixed:** the first run's A5b "FAIL"
   was a quote-precision bug in this gauntlet's check (a point enclosure vs a
   4-decimal corpus quote), not a corpus error — filed here in the corpus's own
   spirit of printing one's own defects.
3. **The checker caught a real schema violation:** the first packaging attempt
   used a provenance string with characters outside the identity-field charset;
   fs-checker refused with `invalid-identity`. Fail-closed behavior confirmed
   against its own author.
4. **D14 edge note:** the phason window's lower edge computes to 3.71 MeV; the
   corpus rounds to "4 MeV". Same class as the corpus's own F-A15-5 pull
   harmonization — a rounding, not an error.
5. **Kill-content:** none fired. Tier 0's prediction ("none expected — the
   chain was spot-verified") held.

## Epistemic notice (binding)
Every PASS verifies a **within-model** claim of a speculative theory.
Reproducible arithmetic is not evidence about nature; GUM's contact with
nature remains its staked experiments (S1, S2, bench, A3 gate).

## Addendum — Group F: NR-D2 murk-kill recomputation (42/42 total)
The archive's only fired kill (T-D2) is explicitly **replication-gated** —
"the ink dries at replication." This is the first external recomputation of
its arithmetic spine (its two *physics* gates — SDiff-reachability of the
reconciliation cost, and hard-wall applicability — are argument-level and not
discharged by arithmetic):

| ID | Check | Result |
|----|-------|--------|
| F1 | overlap coefficient f·p = 7.6×10⁷ at f = 3.75 MeV, p = 4 µm | PASS (7.60×10⁷) |
| F2a | transparency bound from the corpus's KE quote: ε ≤ 2.6×10⁻¹² | PASS (corpus "3×10⁻¹²") |
| **F2b** | **FINDING F-R1**: at the *stated* normalization (v = 10⁻³c, m = 16f), KE = 8×10⁻⁶f — the printed 2×10⁻⁴f is **×25 too large** (corresponds to v = 5×10⁻³c). Literal bound: ε ≤ 1.05×10⁻¹³ | PASS — **conservative-direction slip; the kill fires ~25× harder** |
| F3 | σ/m = π(2p)²/m = 4.7×10¹⁸ cm²/g at mid-band; Bullet ≲ 1 → over by 18+ orders | PASS |
| F4 | classicality λ_dB/p ≈ 1.3×10⁻⁶ ≪ 1 | PASS |

**F-R1 is the campaign's first genuine corpus finding**: an internal
inconsistency in NR-D2 §1's printed KE coefficient. It does not reopen the
kill — every reading tightens it — but it is exactly the class of defect
(cf. the corpus's own C-ν1, "too harsh on itself") that the replication gate
exists to catch. Golden root updated by the added claims (bit-deterministic).

## Addendum 2 — Groups G & H: the remaining arithmetic archive gates (50/50 total)
**Group G — WS-A4 precision-perimeter demand tables:** all five sealed-law
rows (`|c| ≤ δ/ε_max^p`, p = 1, 2, 3, 4, 2/3) reproduce from one consistent
window δ = 1.0×10⁻¹² at ε_max = 6×10⁻⁴ (1.67e-9 / 2.78e-6 / 4.63e-3 / 7.72 /
1.41e-10 vs the corpus's 1.7e-9 / 2.8e-6 / 4.6e-3 / 7.7 / 1.4e-10); the
compositeness scale Λ* = 1.97×10¹⁰ GeV; the chirally-protected margin
~1.5×10¹⁵ and the **thin linear corner at 38.6** (corpus: "∼40 — thin,
printed rather than hidden") both confirmed; the eEDM coefficient 1.06×10⁻¹⁹.

**Group H — NR-D1b core-energy floor:** c_h band [8π/3, 8π] = [8.378, 25.13];
the floor clears the eternal threshold 5.6 by ×1.496 ("×1.5"); S = π·8² =
201.1 vs survival ≈92 gives exponent margin 109 > 100 ("τ ≥ e^{100+} t_U") —
and the corpus's G→8 rounding before squaring is conservative-direction (the
honest floor gives S = 220, margin 128). **W-eternal replicated: the fate is
decided by the one-sided bound, as the instrument claims.** The replication
gate's two *physics* inputs (c_h in the chiral medium; the r_min ≈ 1/f floor)
remain argument-level, not discharged by arithmetic.

## Addendum 3 — Group I: the ⟨r10⟩ error band (52/52 total; FINDING F-R3)
The bridge band's components reproduce exactly (s_A = 0.845, s_κ = 0.732,
quadrature 1.118 — matching WS-nu-P2's own check-lines). But under standard
propagation the **printed ρ = +0.45 with same-sign partials gives ±1.344**;
anti-correlating coupling gives ±0.833; the corpus's printed **±0.90
corresponds to an effective ρ = −0.356, not the printed +0.45**. The phrase
"partial correlation cancellation" (App. K.2) implies the anti-correlating
reading, in which case ±0.90 is mildly conservative; under the literal
reading the band is ~1.5× understated. **The S1 stake window is unaffected**
(its floors are physical: the oscillation floor and the Σ bound), but the
covariance sign convention should be pinned — filed as F-R3, the campaign's
third finding.
