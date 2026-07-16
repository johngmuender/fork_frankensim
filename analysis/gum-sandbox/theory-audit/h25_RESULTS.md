# H2.5 — per-sector Seeley–DeWitt 𝔞₁: the App. E.6 positivity computation (G2 / F-T4-4)

**Campaign:** replication Phase H2, follow-on №5 of the theory audit (`T4_gravity_oneworld.md` §G2;
ledger item F-T4-4). **Date:** 2026-07-16. **Code:** `h25_sdw.py` (python3 + sympy/numpy,
deterministic, seed 20260716, runtime 3.0 s, 27/27 internal checks PASS).
**Raw output:** `h25_results.json`. **Status: measured/derived facts below; promotion is the
coordinator's call.**

## 0. Verdict in one line

**E.6's premise 𝔞₁^(s) > 0 is NOT a heat-kernel fact and is NOT derivable from the corpus's own
operators: the per-sector R-coefficient is 𝔞₁ = −(5p+q)/(12(2p+q)) per mode, a function of a
constitutive response ray (p, q) the corpus never fixes — NEGATIVE on every same-sign response
ray including the corpus's own displayed perturbation form (cone-only, 𝔞₁ = −1/12 per mode),
positive only on an 8.5%-measure opposite-sign wedge; and the knot band-edge sector, fermionic
by the corpus's own Theorem c″, has 𝔞₁ = −1/3 per Dirac field, so w_s < 0 as printed in every
geometric reading. DEFECT-CANDIDATE (conditional on reading), GAP confirmed and sharpened.**

## 1. The question (restated)

Omega paper VI.B–VI.C + App. E.6: Γ_ind = Σ_s (N_s/32π²)∫√(−g_s)[2Λ_s⁴𝔞₀ + Λ_s²𝔞₁R[g_s] + …],
summed over **B2±, B3, B4, and knot band-edges**, g_s = ḡ + 2δ_s c² n⊗n, and (6.2):
c_GW² = Σ w_s c_s²/Σ w_s with **w_s = N_s Λ_s² 𝔞₁^(s) > 0**, justified in E.6 by the single
phrase "heat-kernel positivity". The T4 audit graded this a GAP: positivity is a theorem about
𝔞₀, not 𝔞₁ (𝔞₁-integrand = tr(R/6 − E), field-content-dependent). This run computes 𝔞₁^(s)
per sector from the corpus's own quadratic action.

## 2. Conventions (stated)

- Euclidean signature (+,+,+,+), d = 4, static backgrounds; one-loop determinants over ℝ⁴ with
  flat functional measure (the medium's fields carry the flat measure — this is what forces the
  divergence-form operator below).
- Laplace-type reduction per Gilkey/Vassilevich: D = −(g^{μν}∂_μ∂_ν + a^σ∂_σ + b) ⟹
  D = −(Δ_g + E) with the standard ω_δ/E formulas; Tr e^{−tD} = (4πt)^{−d/2}Σ tⁿ∫√g tr aₙ,
  **a₁-integrand = tr(E + R/6)**. Curvature sign: R(S³, ρ) = +6/ρ² (verified symbolically).
- **𝔞₁^(s) in the corpus's (6.1) sense = the coefficient of R** in tr(E + R/6): the sector-gap
  part of E (−m_s²) renormalizes the Λ⁴/potential sector, not the R term. (The uncharitable
  literal reading — 𝔞₁ includes −m_s² — is reported in §6: it is negative outright for every
  gapped sector on weakly curved backgrounds.)
- Bundle traces: after polarization decoupling every micropolar sector is N_s identical
  scalars (tr 1 = N_s); Dirac option is rank 4 with Lichnerowicz E = −R/4 (torsion-free import).
- Baseline: Case A (objective mass), χ₃ = 0 — the corpus's vacuum. With χ₃ ≠ 0 the extra
  parity-odd first-order term is absorbed into the Gilkey connection (still Laplace-type per
  polarization); corrections to 𝔞₁ enter at O(χ₃²) and do not move the sign map.

## 3. Sector operators (derived symbolically from the corpus's W₂, validated against tier-1)

The 6×6 K(k) was rebuilt in sympy exactly as `tier1-spectrum/spectrum.py` assembles it; every
closed form below was cross-validated against the direct 6×6 eigensolve on the full k-grid
(max rel dev 1.4×10⁻⁹, and 10⁻¹³-level agreement with the campaign's published eigensolve
numbers at k = 1, e.g. B3: 5.181847289748).

| sector | N_s | inertia w_t | stiffness w_x | gap M | dispersion (exact/derived) | in E.6 sum |
|---|---|---|---|---|---|---|
| B1 long. acoustic | 1 | ρ₀ | λ+2μ | 0 | ω² = c_L²k² exact, c_L² = (λ+2μ)/ρ₀ | **no** (elliptic constraint sector) |
| B2± light doublet (photon) | 2 | ρ₀ | μ | 0 | ω² = (μ/ρ₀)k² + O(k⁴/A_c), A_c = 2μ_c + m_V²/2 | yes |
| B3 heavy KG doublet | 2 | J | Jc_ψ² | (4μ_c+m_V²)/J | ω² = ω₀² + c_ψ²k² + O(k⁴); **c_ψ² = (β+γ)/(2J) + A_c/(2ρ₀)** (= 5.75 at benchmark; new closed form, consistent with tier-1's measured 5.181847… at k = 1 to 10⁻¹⁴) | yes |
| B4 long. twist | 1 | J | α+β | (4μ_c+m_V²)/J | ω² = ((α+β)k² + 4μ_c + m_V²)/J **exactly KG** (no k⁴ term) | yes |
| knot band-edges | per species | 1 | c² | ω₀² | ω² = ω₀² + c²k² + O(k⁴a²) (corpus (7.1)); **fermionic by Theorem c″** | yes |

Reducibility, stated honestly:
- The **full 6×6 operator is NOT Laplace-type** (its leading symbol carries distinct cones
  c_L², μ/ρ₀, c_ψ², (α+β)/J — no single metric ×𝟙). The per-sector decomposition of (6.1) is
  therefore *forced*, and is legitimate only after (i) branch projection (removing the
  first-order gyroscopic u–φ mixing, |K_{uφ}| = A_c k, exactly as the corpus's own "eliminating
  the heavy relative field") and (ii) truncation at second order in derivatives. Declared
  truncation error: O(k²/gap) ~ O(R/m̃²) — subleading on weakly curved backgrounds, so it does
  not move the sign conclusions below. **No sector is NOT-REDUCIBLE**; the audit point is what
  the reduction *yields*.
- Admissibility scan (20 000 draws over the Eringen-positive moduli region, λ and α allowed
  negative within 3λ+2μ > 0, 3α+β > 0): every sector cone and gap is positive on the whole
  admissible range (mins ~10⁻⁶ at the range edge, all > 0) — the acoustic metrics exist
  everywhere; **the sign of 𝔞₁ is NOT decided by moduli positivity** (it drops out of the
  R-coefficient entirely, see §4).

## 4. The computation: 𝔞₁ for the corpus's variational sector operator

Every summed sector reduces, at second order in derivatives, to
**D = −w_t(x)∂_τ² − ∂_i w_x(x)∂_i + M(x)** — the divergence form is forced by the corpus's own
energy functional (K = ΣcR†R with position-dependent coefficients) and the flat functional
measure. The Gilkey metric is then *not a choice*: g^{μν} = diag(w_t, w_x δ_ij).

How the defect background ḡ enters the moduli is nowhere specified in the corpus — that is
exactly T4's "no computation of the effective ξ_s is anywhere in the text". Parametrizing the
response at linear order, w_x ∝ 1 + p h(x), w_t ∝ 1 + q h(x), M ∝ 1 + r h(x) (h an arbitrary
static profile), the symbolic reduction gives (all anchors PASS: S³ convention; minimal
Laplacian ⟹ E = 0; conformal operator ⟹ a₁ = 0; unimodular ray recovers +1/6):

- R[g] = (q + 2p)Δh + O(h²);  a₁-integrand = −((5p+q)/12)Δh − M₀rh + O(h²)
- the gap modulation (r) separates cleanly into the potential/Λ⁴ sector — it never touches R;
- **master result (closed form, per scalar mode):**

  **𝔞₁(p, q) = −(5p + q) / (12(2p + q))**

  moduli-independent; only the response ray matters. Positive **iff (5p+q)(2p+q) < 0** — an
  opposite-sign (stiffness up ⟺ inertia down) wedge of angular measure **8.5%** of all rays
  (numeric ray scan 0.0847 vs analytic 0.0848).

Named rays (identical for B1, B2±, B3, B4 and the bosonic-KG knot reading — same functional form):

| response reading | (p, q) | 𝔞₁ per mode | sign |
|---|---|---|---|
| **cone-only — the corpus's own displayed form g_s = ḡ + 2δ_s c²n⊗n (time-time only)** | (0, 1) | **−1/12** | **−** |
| stiffness-only | (1, 0) | −5/24 | − |
| conformal (p = q) | (1, 1) | −1/6 | − |
| GR weak-field (q = −p) | (1, −1) | −1/3 | − |
| unimodular / exact-acoustic (q = −3p; √g = const — the Unruh case) | (1, −3) | **+1/6** | + |
| minimal covariantization imposed by hand (unforced; not the medium's operator) | — | +1/6 · N_s | + |

Two structural pathologies beyond the sign:
- **Pole ray q = −2p:** R vanishes at linear order while the a₁ term does not — there the
  induced term is not of the form 𝔞₁R at all; (6.1)'s ansatz fails structurally.
- **Non-universality:** for a general defect background ḡ the ratio p/q varies in space unless
  the constitutive response is rigidly tied, so "𝔞₁^(s) = a per-sector constant" already
  presupposes a fixed response ray — an assumption nowhere in the text.

## 5. The knot band-edge sector (the fermionic flank)

Corpus (7.1) writes a bosonic KG form; corpus Theorem c″ makes knots fermions. The three
bookkeeping readings:

| reading | 𝔞₁^(knot) | w_s = N_sΛ_s²𝔞₁ as printed |
|---|---|---|
| K1: bosonic KG taken literally | master formula (minimal: +N/6; cone-only: −N/12) | sign per §4 |
| K2: fermionic, formula as printed (N_s > 0 mode count, no statistics sign) | 4(1/6 − 1/4) = **−1/3 per Dirac field** (Lichnerowicz E = −R/4) | **w_knot < 0 in EVERY geometric reading** |
| K3: fermionic, statistics-signed loop (Γ_F = −Tr ln) | contribution = **+2 minimal-scalar equivalents** per Dirac field | positive under the minimal reading, but then the printed formula w_s = N_sΛ_s²𝔞₁ is wrong as written (needs N_s → supertrace-signed count) |

## 6. Sign map and what it does to Theorem VI.1's convexity core

Assembled across readings (bosonic sectors B2±/B3/B4 share one sign; knots separately):

| (geometry reading, statistics bookkeeping) | B2±, B3, B4 | knots | Σ_s w_s | convexity of (6.2) |
|---|---|---|---|---|
| minimal + statistics-signed (double repair) | + | + | + | **SOUND** — but both repairs are unforced and underived; minimal coupling is *false* for the medium's own variational operators except on the unimodular ray |
| minimal + corpus-literal | + | − | ? | **mixed signs ⟹ c_GW² can exit the hull; the slaving bound fails as an inequality** |
| cone-only (corpus's displayed form) + statistics-signed | − | + | ? | **mixed signs ⟹ fails** |
| cone-only + corpus-literal | − | − | < 0 | ratio (6.2) survives formally (common sign divides out) but induced 1/16πG = Σw_s < 0: **inverted-sign Newton constant**, contradicting VI.B's G ~ c³/(ħN_effΛ_UV²) |
| literal-𝔞₁ reading (mass included) | 𝔞₁ = −m_s² < 0 flat-background for every gapped sector (benchmark −21 for B3/B4) | − | — | fails outright; confirms the R-coefficient reading is the only charitable one |

**Verdict: DEFECT-CANDIDATE (conditional) — GAP confirmed and sharpened.** The E.6 one-liner
("heat-kernel positivity") is false as a theorem about 𝔞₁, and the computation shows the
positivity premise fails in the corpus's own displayed perturbation form and in every
same-sign constitutive reading; it holds only on (a) the 8.5% opposite-sign wedge containing
the unimodular ray q = −3p, together with (b) a statistics-sign erratum for the fermionic knot
sector. The negative-weight region, stated precisely: **every response ray with
(5p+q)(2p+q) > 0** (includes cone-only, stiffness-only, conformal, GR-weak-field), **plus the
knot sector as printed in every reading.** Theorem VI.1's algebra (audited SOUND given P1, P2
in T4) is untouched; its premise P2 is now *quantified*: the convexity core survives only if
the corpus derives (i) the unimodular/minimal constitutive response for all four micropolar
sectors and (ii) the supertrace-signed weight for knots — neither derivation exists in the
text or archive. Until then the UHECR-slaved bound |Δc_GW/c| ≤ δ_UV stands only as an
order-of-magnitude estimate (T4's fallback), not an inequality.

## 7. Honest limits of this run

- Linear order in the background modulation: sufficient for the sign at weak curvature (the
  regime of every (6.2) application); O(h²) terms cannot rescue a wrong-sign linear term.
- The (p,q) response is *genuinely underdetermined by the corpus*: this run does not prove
  every reading fails — it proves the assertion "𝔞₁^(s) > 0" is reading-dependent, exhibits
  the failing readings (including the textually closest one), and locates the unique rescuing
  wedge. A corpus derivation landing exactly on q = −3p would discharge the geometric half of
  the gap (and is the obvious repair target); nothing rescues the printed fermionic bookkeeping
  except an erratum.
- Heavy-field-elimination corrections to B2 (the tier-1 O(k⁴/A_c) terms) contribute to the
  R-coefficient only at O(R/m̃²) — truncated, declared, sign-irrelevant.
- The Dirac reduction for knots imports Lichnerowicz on the torsion-free acoustic background;
  torsion corrections (the sector's own EC claim) would shift E = −R/4 by contact terms but
  cannot change the statistics-sign structure of the bookkeeping conflict.

## 8. Deliverables

- `h25_sdw.py` — symbolic derivation + numeric scans (single file, deterministic, 3 s;
  27/27 internal checks).
- `h25_results.json` — machine-readable: dispersions, master formula, ray table, scans,
  sector table, verdict, full check list.
- `h25_RESULTS.md` — this memo.
