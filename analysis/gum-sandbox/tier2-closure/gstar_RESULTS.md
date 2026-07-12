# Phase G3 — G* vs 𝔠₀: the saturated constant in closed form

Reproduce with `python3 gstar_solve.py` (deterministic, no RNG; 31.2 s wall;
writes `gstar_results.json`, `gstar_fig.png`; every printed table in
`gstar_run.log`). Companion files: `gstar_solve.py`, `gstar_results.json`,
`gstar_fig.png`, `gstar_run.log`. Unlike the earlier tier-2 reports, this
step was commissioned to DO the adjudicating mathematics, so verdict
language appears below and is marked as such; every number is measured or
proven, and the proof is verified against the frozen campaign engine.

## 1. Question and answer

**Question** (from axi3_RESULTS.md §5b, axi3_ADJUDICATION.md,
axi4_locked_RESULTS.md §8, TIER4_ADJUDICATION.md §2): the halo-saturated
objective `G* = min (E_static − I/(16π))` at t → 0 converged downward across
families to 2.51442–2.51446 (axi3) / 2.515991 (axi4 locked; same
normalisation — axi4 line 179 calls `axi3_solve.G_from` unchanged), within
~1.2×10⁻⁴ relative of the corpus's deep-BPS constant
𝔠₀ = 128√42/(105π). Is G* = 𝔠₀ exactly?

**Answer (proven, then engine-verified): no.**

```
G*  =  16√2 / 9  =  2.514157444218836…       (EXACT, closed form)
𝔠₀  =  128√42/(105π)  =  2.514753626327608…
𝔠₀ − G*  =  5.961821088×10⁻⁴   (2.3707×10⁻⁴ relative)   — PINNED GAP
𝔠₀ / G*  =  24√21/(35π)  =  1.000237129982…
```

The identity fails *exactly*, not merely at a numerical precision: G* is
algebraic while 𝔠₀ = (algebraic)/π is transcendental, so G* = 𝔠₀ is
impossible. The proximity that motivated the question is the numerical
accident 24√21 = 109.98182 ≈ 35π = 109.95574 (0.0237 %) — equivalently,
the identity would require π = 24√21/35 = 3.14234, wrong in the fourth
digit.

(Bookkeeping note: the task sheet quoted 𝔠₀ = 2.5147350803; the exact value
of 128√42/(105π) is 2.5147536263… — the quoted digits appear transposed.
Nothing below depends on this.)

## 2. Derivation (in the campaign's frozen conventions; nothing refitted)

At t → 0 the saturated objective is `G = E6 + E0 − I/(16π)` with
`E6 = π³∫b² dV`, `E0 = (1/4π)∫(1−q0) dV`, `I = ∫2(q1²+q2²) dV` (axi_solve
sector constants, unchanged). Since E0 and I are derivative-free and E6 is
quadratic in the topological density, G is a sextic-plus-potential
(BPS-Skyrme-type) functional,

```
G = ∫ [ π³ b² + W(q) ] dV,
W(q) = (1/4π)(1−q0) − (1/8π)(q1²+q2²)
     = (1/8π) [ (1−cos ψ)² + sin²ψ cos²Θ ]  ≥ 0,
```

(q0 = cos ψ, (q1,q2,q3) = sin ψ(sinΘ cos φ, sinΘ sin φ, cosΘ)). The
inertia discount −I/(16π) turns the potential into a **Θ-dependent**
effective potential — that is the entire content of the problem.

**(1) Sharp lower bound (rigorous; any degree-1 field, no symmetry
assumed).** Pointwise AM-GM `π³b² + W ≥ 2√(π³W)|b|`, plus the fact that
b dV is 1/(2π²) times the pullback of the target volume form, give for
|deg| = 1:

```
G ≥ 2√π³·(1/2π²)·∫_{S³} √W dΩ
  = (1/√2) ∫∫ sin²ψ sinΘ √[(1−cosψ)² + sin²ψ cos²Θ] dψ dΘ
  = 16√2/15 + 8√2·∫₀^{π/2} sin⁵x cos x ln cot(x/2) dx      (x = ψ/2)
  = 16√2/15 + 8√2·(4/45)                                    (by parts)
  = 16√2/9.
```

axi3's bound π/√2 = 2.22144 is this bound with the cos²Θ term discarded;
keeping it sharpens 2.221 → 2.514157.

**(2) Attainment (explicit closed-form minimiser).** Equality needs
`π³b² = W` pointwise with single-signed b. Within the direction-locked
ansatz Θ = θ this is a per-θ ODE that integrates in closed form. With
`g(x) := arcsin x − x√(1−x²)` (monotone [0,1] → [0,π/2]) the solution is
the implicit surface

```
g( cos(F/2)·sinθ ) = ρ³/(6√2),        ρ = r sinθ  (cylindrical radius),
```

i.e. `F = 2 arccos(min(1, r·x/ρ))` with `x = g⁻¹(min(ρ³/(6√2), π/2))`.
Substituting back, the BPS condition reduces to the algebraic identity
`√[(1−cosF)² + sin²F cos²θ] = 2 sin(F/2)√(1−x²)` at `x = cos(F/2) sinθ` —
exact, so the field saturates the bound identically. Hence

```
inf G  =  G*  =  16√2/9,     attained.
```

Properties of the minimiser (each verified in §4):

- on the symmetry axis it is EXACTLY the campaign's compacton
  `F = 2 arccos(r/2^{5/6})`; support radius R(0) = 2^{5/6} = 1.781797;
- the support is an **oblate compacton**: edge `g(sinθ) = ρ³/(6√2)`,
  equatorial radius R(π/2) = (3√2π)^{1/3} = 2.370984 (matches the
  oblate-waisted core axi3's figure recorded);
- it is **direction-locked** (Θ = θ exactly): the true minimiser needs no
  tilt. axi3's equatorward tilt (ΔΘ ≈ 0.39) was the 18-mode family's way
  of imitating the per-θ profile-shape freedom its radius-warp
  parametrisation lacked; axi4's locked family had the lock right but not
  the shape freedom (its base profile is the spherical compacton, warped
  only in radius);
- G is SDiff-invariant (E6 quadratic in a density; E0, I derivative-free),
  so the minimiser is a whole SDiff orbit; the field above is its locked
  axisymmetric representative. At t → 0⁺ the degeneracy is lifted at
  O(t log t) in G by t(E2+E4) — irrelevant to the t → 0 value.

**(3) Sector values and the saturated tuple, all in closed form** (via the
pushforward dV = dΩ/(2√(πW)) on the minimiser):

| quantity | closed form | value |
|---|---|---|
| E6 | 8√2/9 | 1.257078722109418 |
| E0 | 4√2/3 | 1.885618083164127 |
| I (core) | 64√2π/9 | 31.593834226903937 |
| E0 − I/16π | 8√2/9 (= E6: d-stationary at d = 1) | 1.257078722109418 |
| **G*** | **16√2/9** | **2.514157444218836** |
| L | √(8π)G* = 64√π/9 | 12.604116273105891 |
| c_paper | (4/π)G* = 64√2/(9π) | 3.201124679669711 |
| κ_paper | 1/√2 (exact, unchanged) | 0.707106781186548 |
| I_tot | 8πG* = 128√2π/9 | 63.187668453807870 |
| halo_D | 32√2π/9 | 15.796917113451968 |
| E0_halo | 4√2/9 | 0.628539361054709 |
| E_stat_tot | 24√2/9 = (3/2)G* (identity exact) | 3.771236166328254 |

The axi3 c_paper bracket [2.82843, 3.20146] therefore closes to the point
value **c_paper(t→0) = 64√2/(9π) = 3.2011247**.

## 3. Numerical verification (all gates PASS)

| gate | what it checks | measured |
|---|---|---|
| A1 | sharp bound: 2-D target-space Gauss quadrature vs 16√2/9 | +3.1e-15 |
| A2 | sharp bound: 1-D reduced (closed-form inner integral) | +4.9e-15 |
| A3 | sharp bound: closed-form split 16√2/15 + 32√2/45 | 0.0 |
| A4 | attainer E0 by pushforward quadrature vs 4√2/3 | −1.3e-11 |
| A5 | attainer I by pushforward quadrature vs 64√2π/9 | −6.6e-10 |
| B1 | Newton inversion residual of g⁻¹ (worst of 5002 hard points) | 8.9e-16 |
| B2 | axis profile vs exact compacton 2 arccos(r/2^{5/6}) | 6.0e-14 |
| B3 | analytic F_r (implicit differentiation) vs finite differences | 3.6e-10 rel |
| C | pointwise BPS identity π³b² = W on quadrature nodes (F > 0.05) | ≤ 1.5e-11 |
| D | engine degree at N = 896 | B = 0.999996 |
| D | d-stationarity from engine sectors, worst of ladder | d*−1 = −2.4e-5 (N=160) → −4.5e-7 (N=896) |

**Spectral quadrature of the minimiser** (Gauss–Legendre in μ and in the
edge-regularising variable r = R(θ)(1−v²), analytic F_r; independent of the
closed-form sector evaluation): G = 2.514157444219 at every ladder rung
n = 48², 96², 192², 320², max |G − 16√2/9| = 9.8e-14; E6, E0, I match their
closed forms to ≤ 1.4e-12.

**Frozen-engine ladder** (`axi_solve.sector_pieces`, the exact code that
produced every reported G* in axi3/axi4, REUSED unmodified; fixed domain
1.06·R_eq; d = 1):

| N | G (engine) | G − 16√2/9 | degree B |
|---|---|---|---|
| 160 | 2.514020500 | −1.369e-04 | 0.999824 |
| 224 | 2.514095592 | −6.185e-05 | 0.999911 |
| 320 | 2.514132139 | −2.531e-05 | 0.999962 |
| 448 | 2.514146193 | −1.125e-05 | 0.999981 |
| 640 | 2.514152982 | −4.463e-06 | 0.999992 |
| 720 | 2.514154170 | −3.275e-06 | 0.999993 |
| 896 | 2.514155654 | −1.791e-06 | 0.999996 |

Richardson: free-exponent fit over the last five rungs gives
G∞ = 2.514158041, p = 2.34, fit residual 6.8e-8; the last three pairwise
p = 2 extrapolants are 2.514159504 / 2.514158642 / 2.514158359. Honest
**error budget = 1.5e-6** (spread of all extrapolants);
|G∞ − 16√2/9| = 6.0e-7 — **CONFIRMED within budget**. The convention chain
(sector constants, unit map, degree normalisation) is thereby verified
end-to-end against the same engine that produced the family values; the
task's ≤1e-7 precision target is met by exactness of the closed form, with
the independent engine confirmation at 6e-7 ± 1.5e-6.

**Apples-to-apples**: axi3's 18-mode extended optimum re-evaluated on the
same engine at N = 720 gives G = 2.5144238 (reproduces its JSON value to
all printed digits) — sitting **+2.663e-4 above** the exact minimum; axi3
joint +2.936e-4; axi4's locked ring-G* (radius-warp family, K = 5)
+1.834e-3. Every family value is above 16√2/9, none below — as the bound
requires.

**Stationarity/minimality probes** (engine N = 320, six perturbation
families through the minimiser: contour warps l = 0, 2 at two radii, tilt
modes at two radii, plus dilation): every central-difference dG scales as
h² between h = 0.05 and h = 0.025 (pure truncation; measured ratios
1.8–5.2, vertex offsets |dG/d²G| ≤ 2.3e-3 mode-amplitude), every
d²G/dp² > 0 (1.07–19.06), every perturbed G exceeds G(minimiser), and the
analytic dilation optimum is d* = 0.9999950. Consistent with an interior
global minimum; flat SDiff directions exist by invariance but do not change
G*.

## 4. Adjudication (this step's assigned layer)

1. **The identity G* = 𝔠₀ FAILS.** G* = 16√2/9 exactly (rigorous sharp
   Bogomolny bound + explicit closed-form attainer, both verified above);
   the gap 𝔠₀ − G* = 128√42/(105π) − 16√2/9 = 5.9618e-4 is pinned in
   closed form. It is 408× the engine-verification budget and infinitely
   many σ analytically (algebraic vs transcendental). The precision
   question is closed with a **disproof**.
2. **The downward family trend is explained and terminates below 𝔠₀.**
   axi3/axi4's G* values were converging toward 16√2/9 = 2.5141574, i.e.
   *through and past* 𝔠₀ = 2.5147536, not to it. The families crossed 𝔠₀'s
   value because 𝔠₀ happens to lie 6e-4 above the true constant, not
   because 𝔠₀ is the limit. TIER4_ADJUDICATION's open item ("whether the
   halo-saturated G* and the paper's 𝔠₀ coincide exactly in some limit is
   open") is now closed: they do not.
3. **The √42 signature does not decompose.** 𝔠₀ = 2√(3ê₀𝔦₀/2) carries
   √42 = √2·√21 from the compacton Haar moments ê₀ = 64/(15π),
   𝔦₀ = 256/(105π) and the g = 3/2 dial; the saturated constant carries a
   pure √2 (all π's cancel between the sextic normalisation and the
   potential). No closed-form identity connects them; the ratio is
   24√21/(35π). The two constants agree to 2.4e-4 by coincidence of
   sin²-moment arithmetic — the same class of coincidence the Tier-4
   adjudication already flagged for its threshold identities, and this one
   fails where those held exactly.
4. **Within-model status of the saturated closure is upgraded**: its tuple
   is now fully closed-form — κ_paper = 1/√2, c_paper = 64√2/(9π) =
   3.2011247, G* = 16√2/9 — replacing the bracket [2.82843, 3.20146]. The
   minimiser is the explicit oblate BPS compacton
   g(cos(F/2)sinθ) = ρ³/(6√2), direction-locked, reducing to the
   campaign's spherical compacton on the axis.

## 5. Epistemic notice

Everything above is a *within-model* computation about the defined
saturated-closure objective in the campaign's frozen conventions
(a6 = a0 = m = 1, compacton-anchored unit map, zero calibrated constants).
The proof is elementary (AM-GM + degree pushforward + an integrable ODE)
and its verification uses the campaign's own frozen engine; no statement
here bears on nature, on the corpus's numerical methods, or on readings of
the corpus's text beyond the compositional definitions in force since
Step 1. The disproof concerns the specific conjecture G* = 𝔠₀; what that
means for the corpus's use of 𝔠₀ stays with the coordinator's layer.
