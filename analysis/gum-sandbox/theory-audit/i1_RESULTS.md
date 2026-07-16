# Phase I1 — the ⟨r5⟩ corroboration test (ROADMAP_v6 headline question)

Reproduce with `python3 i1_r5test.py` (deterministic, no RNG; 756.6 s wall at
full resolution, `--quick` 340.7 s; writes `i1_results.json`; every printed
table in `i1_run.log`). Companion files: `i1_r5test.py`,
`i1_results.json`. Sections 1–7 are **measured facts and exact derivations**;
Section 8 is **analysis, clearly marked** — the corroboration claim itself is
adjudicated by the coordinator.

## 0. Question

The corpus's frozen quartic-regime closure release ⟨r5⟩ reports
**𝔠_q = 3.1 ± 0.2** (paper Sec. VII.G, used at 6σ in the N-ν1 shallow-knot
neutrino kill; App. I.8 manifest: "r5 (quartic-regime closure 𝔠_q = 3.1 ±
0.2; **scaling audit**)"). The campaign's repaired (saturated) closure gives
exactly **64√2/(9π) = 3.2011247** at ε → 0 (G3/T3.1), and the F-R5 halo
threshold κ_ours = 1/√(8π) is **t-independent**, so the saturated branch
exists at every ε. Was ⟨r5⟩ measuring the saturated branch?

Corpus quartic-regime definition (Sec. VII.C′): ε ≳ 1, knot Skyrme-stabilized
at core R_Sk = √(a₄/a₂), E2+E4 stabilization — i.e. the same a₂ = a₄ = t dial
continued to t(E2+E4) ≳ (E6+E0). Both that continuation (reading a) and the
literal pure-quartic limit a₆ = 0 (reading b) are computed.

## 1. Formulation (everything in the campaign's frozen conventions)

Conventions: a₆ = a₀ = m = 1, a₂ = a₄ = t; step-1/2 sector definitions;
compacton-anchored unit map (c_paper = c_ours/√(2π³), κ_paper = 2√π κ_ours),
zero calibrated constants. Engines REUSED unmodified: `radial_solve`
(step-1 hedgehog), `axi_solve.sector_pieces` + families A/B (step 2),
`axi3_solve` GenBasis/Objective/Nelder–Mead/finalize_G (step 3),
`gstar_solve.qfun_min` (G3 exact minimiser, gate only).

**The ε-dial** (campaign convention, axi_solve line 8): ε =
t(E2+E4)/(E6+E0) measured on the solved step-1 hedgehog profile. The
kinetic *fraction* of the static energy is ε/(1+ε): at ε ≥ 1 it passes
0.5 — the "0.5+" regime; the dial itself is unbounded, so ε = 2 is
reachable (t = 0.5835).

**Saturated closure at t > 0** (axi3 §4, algebra unchanged): a soft
delocalised halo has vanishing E2/E4 cost at any finite t (cost ∝ t k²η²
→ 0 with wavelength), while dE0 = (1/8π)∫η² dV and dI = 2∫η² dV exactly;
the halo ledger is therefore t-blind, saturation pins κ_ours = 1/√(8π)
(κ_paper = 1/√2) at every ε, and the closure factorises with the
t-sectors carried by the core:

```
G_t = t(E2+E4) + E6 + E0 − I/(16π)   minimised over cores (+ analytic d)
c_sat(ε) = (4/π) G_t*,   κ_paper = 1/√2 exactly,  E_stat_tot = (3/2)G_t*.
```

**Two rigorous analytic floors** (this step; both respected by every
measured point):

1. G_t ≥ G₀ ≥ 16√2/9 (drop t(E2+E4) ≥ 0; G3's sharp bound) ⟹
   **c_sat(ε) ≥ 64√2/(9π) = 3.2011 at every ε** — the saturated branch can
   never sit below its ε → 0 anchor.
2. Faddeev-type bound in these conventions: pointwise Σλᵢ² ≥ 3(λ₁λ₂λ₃)^{2/3},
   Σλᵢ²λⱼ² ≥ 3(λ₁λ₂λ₃)^{4/3}, AM-GM ⟹ E2+E4 ≥ 3π|B| ⟹
   **c_sat(ε) ≥ 3.2011 + 12 t(ε)**.

**Direct (saddle) branch** — the corpus's method class (restricted family,
interior Newton-type fixed point): family A (hedgehog + dilation, the
scaling-rung continuation, exact d-scalings) and family B (spheroidal
shell, the paper's G.5 ansatz; λ-grid 0.15–1.45), fixed-L Routhian minimum
+ clock fixed point, wide clock bracket.

## 2. Gates (all PASS)

| gate | what it checks | measured |
|---|---|---|
| G-A | frozen engine on G3's exact minimiser, N = 448 | G − 16√2/9 = −1.125e-5 (gstar table value −1.125e-5; identical), B = 0.99998 |
| G-B | dial(0.05) reproduces the frozen t | t = 0.00827897 vs 0.00827643 (rel +3.1e-4, inside the dial tolerance) |
| G-C1 | famA direct at ε = 0.05 vs frozen step-2 solA | 2.15826 vs 2.15825 (rel +5e-6) |
| G-C2 | famB direct at ε = 0.05 vs frozen step-2 solB | 2.27884 vs 2.27886 (rel −1e-5) |
| G-C3 | saturated 12-mode at ε = 0.05 vs frozen axi3 R2 | 3.38284 vs 3.38115 (rel +5.0e-4; axi3 used N = 720 and 2.7× the NM budget — both values are family upper bounds, ours sits above as it must) |
| G-C4 | saturated hedgehog-core at ε = 0.05 vs axi3 | 3.44034 vs 3.44018 (rel +5e-5) |
| GQ1 | generalised radial assembler vs `radial_solve.assemble` at a₆ = 1 | dE = 0, dgrad = 1.3e-15, dHess = 4.5e-13 |
| GQ2 | quartic-model closed-form bound vs quadrature | 6.079304072186 vs 6.079304072184 (12 digits) |

## 3. Measured: the ε-dial ladder

| ε (dial) | t | achieved | kinetic fraction | E_stat_rad | μ_true | dial conv (N×2 / rmax+2) |
|---|---|---|---|---|---|---|
| 0.02 | 0.0030105 | 0.02000 | 0.0196 | 3.08276 | 12.89 | −1.4e-6 / +5.5e-7 |
| 0.05 | 0.0082785 | 0.05000 | 0.0476 | 3.18309 | 7.77 | −5.2e-7 / +2.1e-7 |
| 0.1 | 0.0179505 | 0.09999 | 0.0909 | 3.35289 | 5.28 | −2.3e-7 / +9.0e-8 |
| 0.2 | 0.0393006 | 0.19997 | 0.1666 | 3.69970 | 3.57 | −7.0e-8 / +2.4e-8 |
| 0.5 | 0.1128992 | 0.49998 | 0.3333 | 4.78534 | 2.10 | +1.6e-7 / −4.8e-8 |
| 1.0 | 0.2550655 | 0.99994 | 0.5000 | 6.71280 | 1.40 | +3.3e-7 / −9.2e-8 |
| 1.5 | 0.4133485 | 1.49957 | 0.5999 | 8.75580 | 1.10 | +3.7e-7 / −1.0e-7 |
| 2.0 | 0.5835226 | 1.99956 | 0.6666 | 10.89407 | 0.93 | +3.6e-7 / −1.2e-7 |

(radial grids N = 4000, rmax = 6–12 by ε; Newton to |grad| ≤ 1e-12; E_stat_rad
= t(E2+E4)+E6+E0 on the hedgehog solution; ε → 0 anchor E_stat_rad → 32√2/15
= 3.01699 exactly.)

## 4. Measured: the two branches, raw (frozen compacton unit map)

⟨r5⟩ band: 3.1 ± 0.2. Pull = (value − 3.1)/0.2.

| ε | c_dirA | pull | c_dirB | pull | **c_sat** (12-mode) | pull | sat hh-core | analytic floor | grid conv (A/B/sat) |
|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 2.09516 | −5.02 | 2.21254 | −4.44 | **3.27396** | +0.87 | 3.32299 | 3.2373 | 4e-5/5e-5/4e-5 |
| 0.05 | 2.15826 | −4.71 | 2.27884 | −4.11 | **3.38284** | +1.41 | 3.44034 | 3.3005 | 6e-5/7e-5/6e-5 |
| 0.1 | 2.26358 | −4.18 | 2.38975 | −3.55 | **3.56404** | +2.32 | 3.64024 | 3.4165 | 1e-4/1e-4/1e-4 |
| 0.2 | 2.47489 | −3.13 | 2.61280 | −2.44 | **3.94853** | +4.24 | 4.05227 | 3.6727 | 2e-4/2e-4/1e-4 |
| 0.5 | 3.11885 | +0.09 | 3.29453 | +0.97 | **5.21676** | +10.58 | 5.36259 | 4.5559 | 3e-4/3e-4/3e-4 |
| 1.0 | 4.23536 | +5.68 | 4.47982 | +6.90 | **7.55466** | +22.27 | 7.73005 | 6.2619 | 5e-4/5e-4/5e-4 |
| 1.5 | 5.40803 | +11.54 | 5.72679 | +13.13 | **10.08335** | +34.92 | 10.26852 | 8.1613 | 7e-4/8e-4/7e-4 |
| 2.0 | 6.63461 | +17.67 | 7.03225 | +19.66 | **12.75624** | +48.28 | 12.94273 | 10.2034 | 1e-3/1e-3/8e-4 |

Measured facts of record:

- **Both branches RISE with ε** (the T3 §5.1 sign flip, now measured through
  the whole quartic regime): the corpus's printed falling law
  𝔠₀[1 − 0.42ε^{2/3}] is wrong in sign for every branch of the repaired
  model. Small-ε linear coefficients: saturated +1.14ε (T3's +1.15ε class),
  famA/famB ≈ +1.0ε (Step 2's −1.03ε^{1.0} deficit class, sign-flipped
  vs (4.8)).
- The saturated curve respects both analytic floors at every point and
  **never comes back down to 3.1** (rigorously ≥ 3.2011 for all ε).
- The direct branch crosses 3.1 at **ε ≈ 0.491 (famA) / 0.414 (famB)** —
  kinetic fraction ≈ 0.29–0.33, i.e. *not* the quartic regime under any
  reading of "ε ≳ 1".
- Direct-branch κ_paper (A/B): 0.938/0.903 at ε = 0.02 rising to
  1.146/1.101 at ε = 2 — **above the t-independent threshold 1/√2 at every
  ε**: every direct quartic-regime solution remains a halo-unstable saddle
  of the fixed-L Routhian, exactly as at ε = 0.05 (F-R5). A restricted
  Newton–Krylov (the corpus's method class) would sit on this branch.
- Saturated-core diagnostics: d* ∈ [1.31, 1.54] (interior; d-window 0.35–2.8
  never binding), degree B = 0.99940–0.99996, core sector ratio at the G-optimum
  slightly below the dial (0.83 at dial 1.0), total kinetic fraction at
  saturation 0.365/0.435/0.479 at ε = 1/1.5/2. The soft `amax` shrink-cap
  contact seen by axi3 recurs at ε ≥ 0.05 (halo-like skirt directions are
  G-flat; axi3 measured the wall costs ≤ 4e-4 in G).

## 5. Measured: normalisation readings ("scaling audit" enumeration)

The I.8 manifest says r5 contained a *scaling audit*. A raw c(ε) varies by
×1.7 across ε ∈ [1, 2], so a quartic-regime number quoted to ±0.2 must have
had its t-scaling removed (or be pinned at a specific undisclosed ε).
Enumerated conventions, both branches (pulls vs 3.1 ± 0.2 in parentheses):

- **N0** raw c(ε) (table above);
- **N1** c(ε)/(1+ε) — divide by the dial's kinetic-growth factor;
- **N2** c(ε)/[E_stat_rad(ε)/(32√2/15)] — divide by the measured
  hedgehog static-energy (mass) growth, the model's own ε → 0-anchored
  scale factor (non-circular: the normaliser comes from the step-1 mass
  law, not from either closure branch).

| ε | E-growth | sat N0 | sat N1 | **sat N2** | dirB N0 | dirB N1 | dirB N2 |
|---|---|---|---|---|---|---|---|
| 0.02 | 1.0218 | 3.274 (+0.9σ) | 3.210 (+0.5σ) | **3.204 (+0.5σ)** | 2.213 (−4.4σ) | 2.169 (−4.7σ) | 2.165 (−4.7σ) |
| 0.05 | 1.0551 | 3.383 (+1.4σ) | 3.222 (+0.6σ) | **3.206 (+0.5σ)** | 2.279 (−4.1σ) | 2.170 (−4.6σ) | 2.160 (−4.7σ) |
| 0.1 | 1.1113 | 3.564 (+2.3σ) | 3.240 (+0.7σ) | **3.207 (+0.5σ)** | 2.390 (−3.6σ) | 2.173 (−4.6σ) | 2.150 (−4.7σ) |
| 0.2 | 1.2263 | 3.949 (+4.2σ) | 3.290 (+1.0σ) | **3.220 (+0.6σ)** | 2.613 (−2.4σ) | 2.177 (−4.6σ) | 2.131 (−4.8σ) |
| 0.5 | 1.5861 | 5.217 (+10.6σ) | 3.478 (+1.9σ) | **3.289 (+0.9σ)** | 3.295 (+1.0σ) | 2.196 (−4.5σ) | 2.077 (−5.1σ) |
| 1.0 | 2.2250 | 7.555 (+22.3σ) | 3.777 (+3.4σ) | **3.395 (+1.5σ)** | 4.480 (+6.9σ) | 2.240 (−4.3σ) | 2.013 (−5.4σ) |
| 1.5 | 2.9022 | 10.083 (+34.9σ) | 4.033 (+4.7σ) | **3.474 (+1.9σ)** | 5.727 (+13.1σ) | 2.291 (−4.0σ) | 1.973 (−5.6σ) |
| 2.0 | 3.6109 | 12.756 (+48.3σ) | 4.252 (+5.8σ) | **3.533 (+2.2σ)** | 7.032 (+19.7σ) | 2.344 (−3.8σ) | 1.947 (−5.8σ) |

Measured facts:

- **N2 is the only enumerated convention that makes 𝔠_q an ε-stable O(1)
  number** on either branch: across the whole quartic regime (ε ∈ [1, 2])
  the N2-saturated value moves by 4% (3.395 → 3.533) and the N2-saddle
  value by 3% (2.013 → 1.947); N0 moves by ×1.7, N1 by 13%.
- Under N2, the two branches are separated by **≈ 1.4–1.6 in c, i.e. ≈ 7σ
  in ⟨r5⟩ units**: the datum 3.1 ± 0.2 is diagnostic between them.
- The N2-saturated curve's ε → 0 anchor is the exact repaired constant
  64√2/(9π) = 3.2011 (pull **+0.51σ**); the N2-saddle anchor is the scaling
  rung 2.0533 (pull −5.23σ).

## 6. Measured: reading (b) — the pure quartic limit a₆ = 0

With a₆ = 0 the model E = t(E2+E4) + E0 has t as pure scale (x → √t x gives
E, G ∝ t^{3/2}): a canonical normalisation exists only at a₂ = a₄ = 1
(t = 1), and *any* 𝔠_q extracted from this reading is convention-relative.
Measured at t = 1 (own hedgehog radial solve, N-ladder 3000/6000, rmax
12/16, ΔE ≤ 1.5e-6; virial ≤ 7e-6; engine cross-check ≤ 3e-3):

| quantity | value | floor (analytic, rigorous) |
|---|---|---|
| hedgehog E_static (radial) | 12.70600 | — |
| direct closure (famA): c_paper | 8.2657 (κ_paper = 1.2015) | — |
| saturated G_q* (12-mode) | 12.2509 → c = **15.598** (hh-core 15.660; conv 2.2e-3) | G_q ≥ 3π ⟹ c ≥ 12 exactly |
| Faddeev-type | E2+E4 ≥ 3π = 9.42478 (measured at solution: 12.484) | — |
| quartic+potential BPS-submodel bound | E4+E0 ≥ 16·2^{1/4}Γ(7/4)Γ(3/2)/Γ(13/4) = **6.07930** (closed form; quadrature check 12 digits) | derivation: E4-density ≥ 3(2π²b)^{4/3}, 3:1 AM-GM with the potential, degree pushforward |
| saturated quartic floor | G_q ≥ (2^{−1/4}/π)∫_{S³}[(1−cosψ)² + sin²ψcos²Θ]^{1/4} dΩ = **5.47137** (quadrature) | not sharp (isotropy + split not co-attainable) |

**Fact: reading (b) puts both branches at c ≥ 8.3 (direct) / ≥ 12 (saturated,
rigorous) in the frozen unit map — nowhere near 3.1.** If the corpus's
quartic regime meant a₆ → 0 literally, its 𝔠_q can only be compared after
an unknowable renormalisation (the frozen ⟨r5⟩ code is absent from the
archive).

## 7. Convergence / determinism ladder (every quoted number)

- **Dial**: N = 4000 → 8000 and rmax → rmax+2 shift the achieved ε by
  ≤ 1.4e-6 relative (worst ε = 0.02); Newton residual ≤ 1.6e-12.
- **Direct branch**: engine N = 320 vs 224: |Δc|/c from 4e-5 (ε = 0.02) to
  1.0e-3 (ε = 2). Frozen-reference agreement at ε = 0.05 (N = 512 archive
  values): famA +5e-6, famB −1e-5. Clock residual ≤ 1e-14 (bisection).
- **Saturated branch**: final N = 480 vs coarse 320: |Δc|/c from 4e-5
  (ε = 0.02) to 8e-4 (ε = 2). Mode ladder: hedgehog core (0 modes) sits
  +1.5–2.5% above the 12-mode optimum at every ε; axi3's full-budget
  12-mode value at ε = 0.05 sits 5.0e-4 *below* ours — every family value
  is an upper bound on the true G_t*, so the family-truncation systematic
  is downward-open but bounded by the hh-core/12-mode gap (≤ 2.5%),
  35× smaller than the branch separation in Section 5. The exact ε → 0
  anchor 16√2/9 is engine-verified independently (gate G-A).
- **Quartic model**: radial ladder ΔE ≤ 1.5e-6 across N and rmax; closed-form
  bound vs quadrature 12 digits; saturated conv 2.2e-3.
- **Determinism**: no RNG anywhere (fixed NM seeds, golden sections,
  bisections); `--quick` (half grids, one-third budgets) reproduces every
  ε-curve number to ≤ 8e-4 relative and the reading-(b) saturated value to
  3.5e-3.
- All measured saturated values respect the rigorous floors of §1
  (`floor_ok` true at every ε in `i1_results.json`).

## 8. VERDICT (analysis layer — marked as such)

**Reference pulls of ⟨r5⟩'s 𝔠_q = 3.1 ± 0.2:**

| candidate value | source | pull |
|---|---|---|
| **64√2/(9π) = 3.20112** | saturated branch, exact t-independent anchor (repair) | **+0.51σ** |
| N2-saturated at ε = 1 / 1.5 / 2 | this step, quartic regime, scaling-normalised | +1.5σ / +1.9σ / +2.2σ |
| corpus's own 𝔠₀ = 2.51475 | paper (4.10) | −2.93σ |
| corpus's own ε-law at ε = 1: 𝔠₀[1−0.42] = 1.459 | paper (4.8) extrapolated | −8.2σ |
| scaling rung 2.05329 | C1 (= N2-saddle anchor) | −5.23σ |
| honest G.5 spheroidal endpoint 2.16873 | T3 C1′ | −4.66σ |
| N2-saddle at ε = 1 / 1.5 / 2 | this step | −5.4σ / −5.6σ / −5.8σ |
| raw saddle at ε ≳ 1 | this step | ≥ +5.7σ |
| raw saturated at ε ≳ 1 | this step | ≥ +22σ |
| reading (b) direct / saturated (8.27 / 15.60, floor 12) | this step | +26σ / +62σ |

**Which branch does 3.1 ± 0.2 match?**

1. **Under the literal raw reading (N0)** there is no quartic-regime match
   at all: at ε ≳ 1 the saddle branch is ≥ +5.7σ away and the saturated
   branch ≥ +22σ. The only raw crossing of 3.1 is the *saddle* branch at
   ε ≈ 0.41–0.49 — kinetic fraction ≈ 0.3, not "ε ≳ 1" under any reading
   of Sec. VII.C′, and a transient curve-crossing rather than a regime
   value. A raw reading also cannot explain the ±0.2 precision (raw c
   varies ×1.7 across the quartic regime).
2. **Under the scaling-audit reading (N2 — the only enumerated convention
   in which a quartic-regime 𝔠_q is an ε-stable number, and the manifest
   says r5 did a scaling audit)** the datum is decisive: **3.1 ± 0.2
   matches the saturated branch (+0.5σ at the anchor, +1.5σ to +2.2σ
   across ε ∈ [1, 2]) and excludes the saddle branch — the corpus's own
   printed method class — at 5.4–5.8σ.** The two branches are ≈ 7σ apart
   in r5 units; the number picks one, and it is not the printed one.
3. **Independent sign fingerprint**: 𝔠_q = 3.1 sits *above* the corpus's
   own deep-BPS 𝔠₀ = 2.515 (+2.9σ), while the corpus's printed ε-law (4.8)
   *falls* with ε (it would put the quartic regime at ≈ 1.5, −8σ). A
   quartic-regime closure larger than the ε → 0 constant requires a
   *rising* ε-law — which, measured here, both repaired branches have and
   only the saturated branch combines with a value near 3.1–3.2. Within
   the corpus's own printed constants there is no way to generate 3.1;
   within the repaired model there is exactly one: the saturated branch
   read at (or scaled back to) its t-independent anchor 64√2/(9π).
4. **What ⟨r5⟩ cannot have been (within this model)**: a raw quartic-regime
   measurement on either branch (wrong by ≥ 5.7σ); the corpus's printed
   ε-law applied at ε ≳ 1 (−8σ); the pure a₆ = 0 model in the frozen unit
   map (c ≥ 12 rigorously on the saturated side, 8.3 direct); an N1-type
   (1+ε) normalisation (neither branch: +3.4σ/−4.3σ at ε = 1).

**Verdict statement (for the coordinator).** The measured structure is
*corroboration-consistent and saddle-excluding*: under every reading in
which a quartic-regime 𝔠_q = 3.1 ± 0.2 is a meaningful ε-stable quantity,
the archive's frozen number selects the SATURATED branch — the repair's
branch, whose exact anchor 64√2/(9π) = 3.2011 it matches at 0.51σ — and
excludes the corpus's own printed/restricted branch at ≈ 5σ, its printed
𝔠₀ at 2.9σ, and its printed ε-law at 8σ. What blocks an outright
"corroboration proven": the ⟨r5⟩ normalisation is NOT-AUDITABLE-FROM-TEXT
(the frozen r5 code is absent from the archive — the same discharge route
(a) as T3 §4); a literal raw reading matches neither branch in the quartic
regime and would instead pin a fresh discrepancy (the second of the two
outcomes the roadmap named); and the one deflationary alternative — the
saddle branch crossing 3.1 raw at ε ≈ 0.45 — would require the corpus's
"quartic regime" to mean a kinetic fraction of 0.3, contradicting its own
Sec. VII.C′ definition. On the arithmetic, the archive's ⟨r5⟩ number
carries the saturated closure's fingerprint and cannot be reproduced from
the corpus's own printed closure chain.

## 9. Caveats (the price list, printed honestly)

1. **The N2 normaliser is this audit's construction**, chosen as the
   model's own mass-growth factor (step-1 hedgehog E_static, ε → 0-anchored
   at 32√2/15); the corpus's actual audit convention is unknown. The
   branch-selection conclusion is robust across N1/N2 (saddle excluded at
   ≥ 3.8σ in both; saturated preferred in N2, neither in N1), but the
   *size* of the saturated-branch pull at ε ≳ 1 (1.5–2.2σ) is
   convention-dependent.
2. All saturated values at ε > 0 are **12-mode family upper bounds** on
   G_t* (systematic downward-open, bounded ≤ 2.5% by the mode ladder;
   exact only at the ε → 0 anchor). Tightening them would *lower* c_sat
   toward 3.1-compatibility at small ε and change nothing at ε ≳ 1.
3. The saturated factorisation at t > 0 (halo t-blindness) is axi3's
   algebra, verified there at ε = 0.05 to 1e-16 (internal identity) and
   by the engine halo probe; it is inherited, not re-derived, here.
4. The direct branch is measured on the two *restricted* step-2 families;
   the corpus's frozen solver family is unknown. axi3's enlarged-basis
   direct value at ε = 0.05 (3.1871, wall-pinned, non-converged,
   family-size-dependent 2.61 → 2.94 under cap growth) shows a
   sufficiently large *unstable-basis* direct computation can transiently
   pass near 3.1 while running toward saturation — coincidence-class, but
   it is the one reading under which a *direct* corpus code lands 3.1
   with a saturated-branch cause (the halo channel pulling it upward).
5. Everything is within-model (the campaign's frozen conventions and unit
   map); nothing here bears on nature or on the corpus's numerical
   methods beyond the compositional definitions in force since Step 1.
   The corroboration *claim* — what ⟨r5⟩'s authors measured — stays with
   the coordinator.
