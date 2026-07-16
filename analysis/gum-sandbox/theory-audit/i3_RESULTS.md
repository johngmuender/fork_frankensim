# I3 — the joint convention solve (𝔭 × 𝔟_eff)

Reproduce: `python3 i3_conventions.py` (deterministic, ~1 s; reads only the
frozen campaign records, writes `i3_results.json`).  Companion files:
`i3_conventions.py`, `i3_results.json`.

## The two open items and why they are one problem

Two spec-underdeterminations survived single-constraint recovery:

* **𝔭 = 0.84 ± 0.03** (App B.4/I.1 radial tail amplitude).  The Tier-2b
  Step-1 pilot solved the ε = 0.05 profile exactly (μ, R\*, form all PASS)
  but none of five natural amplitude conventions (0.365–5.06) landed in
  band; nearest `A e^{−μR*}/2π = 0.805` (~1.2σ low).
* **𝔟_eff = 42 ± 6** (App I.2 bond constant).  G4/G5b closed the
  bond-equation SHAPE in band (x₀ = 1.95 ± 0.07 vs 1.92 ± 0.08) but the
  absolute normalization under the App-A reading `𝔟 = C_dμ/(2π𝔭²)`,
  𝔭 = 0.84, came out **~7 orders of magnitude high** (1.4×10⁷–7.5×10⁸).

They are coupled: 𝔟's formula contains 𝔭².  A single convention assignment
must satisfy both simultaneously, plus everything already pinned
(μ = 7.7725 = m/√(2a₂), R\* = 2^{5/6}, the compacton unit map, the in-band
bond location x₀).

## Inputs (all read from the frozen records)

| quantity | value | source |
|---|---|---|
| A (k₁-normalized tail amplitude) | 5.23331×10⁶ | radial_results.json (2N) |
| μ (frozen referee) | 7.77254686 | radial_results.json |
| R\* = 2^{5/6} | 1.78179744 | radial_results.json |
| μR\* | 13.8491 | derived — **the only large number available** |
| e^{−μR\*} | 9.6696×10⁻⁷ | derived — the edge envelope |
| f_fit(R\*) | 0.39177 | radial_results.json |
| C_d — G4 relaxed-anchored corner | 7.865×10⁶ ± 41% | twoknot_results.json (superseded: G5b showed cap-400 under-convergence) |
| C_d — G4 seed corner | 1.7211×10⁸ ± 0.5% | twoknot_results.json |
| C_d — G4 seed central4 | 1.4756×10⁸ ± 0.5% | twoknot_results.json |
| C_d — G5b relaxed cap-1200 | 4.2608×10⁸ ± 19% | twoknot_refine_results.json |

## Enumeration

    𝔭_cand = base · e^{−γ_p μR*} · μ^{a} R*^{b} · F_p^{±1},   base ∈ {A, f_fit(R*)}
    𝔟_cand = C_d · μ^{c} R*^{e} · e^{−γ_b μR*} / (F_b^{±1} · 𝔭_cand²)

with γ ∈ {0,1(,2)}, |a|,|b|,|c|,|e| ≤ 2 or 1, C_d over the four estimators,
and F drawn from the named factor family {m·πᵏ : m ∈ {1, √2, 3/2, 2, 2√2,
3, 4, 6, 8, 16/15, 4/3, 15/8, 16}, k ∈ {0,1,2}} (median multiplicative
gap 1.067×).  Every assignment is scored jointly: pull_𝔭 against
0.84 ± 0.03, pull_𝔟 against 42 ± 6 with the estimator's C_d error folded in.
15,322 assignments kept at (|pull_𝔭| ≤ 2, |pull_𝔟| ≤ 3).

## Result 1 — the discrete pin: the 7-order gap is exactly ONE edge envelope

Residual factor F_req that the convention family must supply at the baseline
`𝔟 = C_d μ e^{−γ_b μR*}/(42 𝔭²)`, 𝔭 = 0.84:

| C_d estimator | γ_b = 0 | **γ_b = 1** | γ_b = 2 |
|---|---|---|---|
| G4 relaxed corner | 2.06×10⁶ | **1.9947 (≈ 2 to 0.3%)** | 1.9×10⁻⁶ |
| G4 seed corner | 4.51×10⁷ | **43.65** | 4.2×10⁻⁵ |
| G4 seed central4 | 3.87×10⁷ | **37.42** | 3.6×10⁻⁵ |
| G5b relaxed 1200 | 1.12×10⁸ | **108.1** | 1.0×10⁻⁴ |

γ_b = 0 (no envelope) needs 10⁶–10⁸; γ_b = 2 (edge-to-edge screening)
needs 10⁻⁶–10⁻⁴ — both unbridgeable by any 2π/4π/√2-class factor.
**γ_b = 1 needs 2–108: only ONE factor of e^{−μR\*} = 9.67×10⁻⁷ brings the
bond amplitude into the named-factor range, for every estimator.**  This is
convention-family-independent (it is a choice among values 10⁶ apart) and
is the joint solve's hard finding: the frozen pipeline's calibration-zone
𝔟 normalization references screening from the knot **edge** (one envelope:
e^{−μ(d−R\*)}-class, edge-to-center), exactly as the frozen 𝔭 references
the tail amplitude **at the edge** (Result 2).  The "~7 orders unrecovered"
of G4/G5b = e^{+μR\*} × O(2–110), nothing else.

## Result 2 — 𝔭's convention class is pinned to the same edge reference

Of the full enumeration, in-band 𝔭 candidates at complexity ≤ 1 (a bare
amplitude over a single named factor) are **exclusively of the
edge-evaluated class** `A e^{−μR*}/F_p`:

| candidate | value | pull |
|---|---|---|
| A e^{−μR\*} / 6 | 0.8434 | +0.11 |
| A e^{−μR\*} / (15π/8) | 0.8591 | +0.64 |
| A e^{−μR\*} / 2π (the pilot's nearest) | 0.8054 | −1.15 |

No raw-A class (A/μ = 6.7×10⁵ is the paper's formula read literally in
pilot units), no f(R\*) class, lands in band without ≥ 3 convention
elements.  The Tier-2b "CONVENTION-LIMITED, nearest 1.2σ low" verdict is
upgraded: 𝔭 is recoverable to +0.1σ as the edge-evaluated tail amplitude
over a 2π-class factor; the residual 6-vs-2π ambiguity (a 4.7% effect) is
below the discrimination power of the ±0.03 band combined with the pilot's
c₂ = 0-only solve (App I.1's "𝔭-grid over the (B.4) window" — the c₂
window — was never scanned; caveat 4).

## Result 3 — strict App-A closure still fails; the joint family closes

The STRICT App-A formula `𝔟 = ℬ/(2π𝔭²)` with ℬ = C_dμe^{−μR\*} (γ_b = 1
now forced) and no further factor:

| 𝔭 | estimator | 𝔟 | pull |
|---|---|---|---|
| 0.8054 (/2π) | G4 relaxed | 14.50 | −3.2 |
| 0.8054 | G4 seed corner | 317.4 | +44 |
| 0.8054 | G5b relaxed | 785.7 | +4.9 |
| 0.8434 (/6) | G4 relaxed | 13.23 | −3.5 |
| 0.8434 | G4 seed corner | 289.4 | +40 |
| 0.8434 | G5b relaxed | 716.5 | +4.9 |

None closes: the literal 2π reading remains open even after the envelope is
supplied.  Allowing ONE more named factor closes it — the joint table
(economy-ranked; cx = total convention complexity; x0(B) = the reading-B
bond-equation re-check, seed/extrapolated arrays, corpus 1.90 ± 0.05):

| 𝔭 convention | 𝔭 | pull | 𝔟 convention | 𝔟 | pull | cx | x₀(B) s/e |
|---|---|---|---|---|---|---|---|
| A e^{−μR\*}/6 | 0.843 | +0.11 | C_d[G4-relaxed] μ e^{−μR\*}/(2𝔭²) | 41.55 | −0.02 | 2 | 1.98/1.82 |
| A e^{−μR\*}/6 | 0.843 | +0.11 | C_d[G4-relaxed] μ e^{−μR\*}/((15/8)𝔭²) | 44.32 | +0.12 | 2 | 1.98/1.81 |
| A e^{−μR\*}/6 | 0.843 | +0.11 | C_d[G4-seed-c4] μ e^{−μR\*}/(4π²𝔭²) | 39.49 | −0.42 | 2 | none/1.92 |
| A e^{−μR\*}/(15π/8) | 0.859 | +0.64 | C_d[G4-relaxed] μ e^{−μR\*}/((15/8)𝔭²) | 42.72 | +0.04 | 2 | 1.98/1.82 |
| A e^{−μR\*}/(15π/8) | 0.859 | +0.64 | C_d[G4-seed] μ e^{−μR\*}/(4π²𝔭²) | 44.40 | +0.40 | 2 | none/1.92 |
| A e^{−μR\*}/6 | 0.843 | +0.11 | C_d[G4-seed] e^{−μR\*}/(2π𝔭²) | 37.24 | −0.79 | 2 | none/1.93 |
| A e^{−μR\*}/(15π/8) | 0.859 | +0.64 | C_d[G5b] μ e^{−μR\*}/(8π²𝔭²) | 54.95 | +1.06 | 2 | none/1.97 |

Most economical in-band (≤ 1.5σ both) assignment per estimator:

| estimator | 𝔟 convention | 𝔟 | pull_𝔟 | max pull |
|---|---|---|---|---|
| G4 relaxed corner | μ e^{−μR\*}/(2𝔭²), 𝔭 = A e^{−μR\*}/6 | 41.55 | −0.02 | 0.11 |
| G4 seed central4 | μ e^{−μR\*}/(4π²𝔭²), 𝔭 = /6 | 39.49 | −0.42 | 0.42 |
| G4 seed corner | μ e^{−μR\*}/(4π²𝔭²), 𝔭 = /(15π/8) | 44.40 | +0.40 | 0.64 |
| G5b relaxed 1200 | μ e^{−μR\*}/(8π²𝔭²), 𝔭 = /(15π/8) | 54.95 | +1.06 | 1.06 |

**Every estimator admits a ≤ 1.5σ joint solution with exactly the same
structure: one edge envelope + one named O(1)–O(4π²) factor** — but 4,495
of 15,322 enumerated assignments land inside 1.5σ (1,724 inside 1.0σ):
the factor identity is not unique (caveat 1).

## Result 4 — validation against the already-pinned quantities

* μ = 7.7725, R\* = 2^{5/6}, the compacton unit map, the direct
  x₀ = 1.95 ± 0.07, and every gate of Tier-2b/G4/G5b are **amplitude-
  convention-free** (shape/energy measurements) — untouched by any
  assignment above, by construction.
* **Reading-B re-check** (running 𝔟_eff(x) crossing +42, re-scaled to each
  candidate convention, on the archived seed / cap-extrapolated arrays):
  the crossings stay at 1.81–1.98 vs the corpus's 1.90 ± 0.05 (within
  ~1.7σ); for the seed-C_d conventions the extrapolated arrays cross at
  **1.92** — on the corpus number.  "none" entries mean the re-scaled
  plateau sits just below 42, i.e. under the reconstructed convention
  **𝔟_eff = 42 is the pair-law plateau amplitude itself**, not an onset
  threshold — which is the natural reading of the corpus text ("𝔟_eff = 42
  vs 𝔟_∞ ≈ 4.5, the calibration-zone factor") and retro-explains why G4's
  reading B "worked" independently of the 7-order mismatch (any threshold
  crosses in the steep rise when the plateau is 10⁶× too high).

## Verdict

**Partial reconstruction — the discrete structure of the convention is
pinned; the residual O(1) factor is degenerate.  Neither a unique
convention nor a pinned incompatibility.**

What is discharged:

1. **𝔟's "~7 orders unrecovered" is closed structurally.**  The gap is
   exactly one edge envelope e^{+μR\*} = 1.03×10⁶ times a named factor of
   2–110 (inside the estimators' own ×54 systematic).  The frozen
   pipeline's calibration-zone normalization is edge-referenced; it is not
   evidence of inconsistent physics between B.4 and I.2 — on the contrary,
   both residuals resolve under the SAME convention element.
2. **𝔭's convention class is reconstructed**: the edge-evaluated tail
   amplitude A e^{−μR\*} over a single 2π-class factor (best 0.843, +0.1σ),
   upgrading Tier-2b's convention-limited 1.2σ-low verdict to compatible.

What is NOT discharged: the exact factor pair (6 vs 2π on 𝔭; 2 vs 4π² vs
8π² on 𝔟).  The C_d estimator systematic (7.9×10⁶ → 4.3×10⁸ between
cap-400-relaxed, seed, and cap-1200-relaxed) spans the entire named-factor
family, so the residual factor is unmeasurable with the current two-knot
instrument; the suggestive exact hits (G4-relaxed needs 1.9947 ≈ 2 to
0.3%) ride on the estimator G5b superseded and are flagged, not claimed.

## Honest caveats (binding)

1. **Look-elsewhere / degeneracy.**  The named-factor family has median
   granularity 1.067×; the 𝔟 band (±14%) is ~4 grid steps wide and the 𝔭
   band (±3.6%) ~1 step, so in-band hits are guaranteed by construction —
   4,495 assignments land ≤ 1.5σ.  Only the DISCRETE findings (γ_b = 1;
   the A e^{−μR\*} class for 𝔭) carry evidential weight: they are choices
   among alternatives 10⁵–10⁶ apart, not 1.07× apart.
2. **C_d is not a converged measurement.**  All two-knot runs ended at
   iter_cap; the relaxed amplitude moved ×54 from cap 400 to cap 1200 and
   was still rising; nu pinned at fit bounds; the core-repulsion
   exponential class is itself an analysis assumption (G4 caveat 5).  The
   residual-factor column inherits all of this.
3. **The ε = 0.05 overlap caveat (G4) stands.**  μR\* = 13.85: at bond
   separations the cores overlap and the linear-zone multipole picture is
   marginal; 𝔟_eff at the benchmark dial is a calibration-zone quantity
   (the corpus itself distinguishes 𝔟_eff = 42 from 𝔟_∞ ≈ 4.5), so the
   edge-envelope finding describes the frozen pipeline's bookkeeping, not
   a far-zone physical constant.
4. **The B.4 𝔭-grid was never scanned.**  App I.1 specifies "𝔭-grid over
   the (B.4) window" — the c₂ ∈ [−m̃²/2, m̃²/6) potential window; the
   pilot solved c₂ = 0 only.  A c₂-dependence of even a few % would
   discriminate /6 from /2π (1.26σ apart); that solve, not more
   convention enumeration, is the follow-on with discrimination power.
5. **Reading B is a weak third constraint**: the running-𝔟 crossing moves
   < 0.15 in x while the threshold normalization moves 10⁶× (the rise is
   that steep), so its 1.81–1.98 agreement band cannot break the factor
   degeneracy either.
6. Everything here is WITHIN-MODEL reverse-engineering of a frozen
   pipeline's bookkeeping from two calibration scalars and archived run
   records; nothing validates the underlying theory.

## Files

```
theory-audit/
  i3_conventions.py   the enumeration + scoring + reading-B re-check (~1 s)
  i3_results.json     full tables (envelope requirement, p candidates,
                      canonical readings, 24 top joint rows + per-estimator
                      best, reading-B re-checks, factor family)
  i3_RESULTS.md       this file
```
