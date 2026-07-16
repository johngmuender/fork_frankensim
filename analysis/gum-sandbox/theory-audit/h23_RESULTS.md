# H2.3 — NR-D1b Frank-energy anisotropy scan (the c_h band, decided)

**Campaign:** replication Phase H2, follow-on №4 of the theory audit (`T1_NR_gates.md` §3.3(a)/§6).
**Date:** 2026-07-16. **Code:** `h23_frank.py` (python3 + sympy/numpy/scipy, seed 20260716, runtime 163 s).
**Raw output:** `h23_results.json`. **Status: measured/derived facts below; §§5–6 marked as analysis;
promotion is the coordinator's call.**

## 0. Verdict in one line

**The corpus's band [8π/3, 8π] is recovered EXACTLY with no anisotropy at all — it is the
{hyperbolic, radial} pair of one-constant coefficients under the standard Frank convention
(E_hyp = (8π/3)Kr, E_rad = 8πKr) — so the audit's premise that the band "encodes unstated
anisotropy" is itself wrong: its degeneracy-at-8π computation is correct only for the corpus's
*printed* functional ∫K(∇n̂)², and the two conventions differ by the saddle-splay divergence.
The unstated assumption that closes the corpus's band is a convention, not an anisotropy.
One new finding cuts against the corpus: the hyperbolic profile is not stationary — relaxing
ψ(θ) lowers the one-constant coefficient to 7.7024 (2.4517π, ~8% below 8π/3), so 8π/3 is an
ansatz value, not a texture floor; the eternal clearance shrinks from ×1.496 to ×1.375 and
W-eternal still lands.**

## 1. The question (restated)

NR-D1b §1(i) prints E₁ = c_h·K·r_geo with c_h ∈ [8π/3, 8π] [IM: standard defect elasticity],
and uses c_h,min = 8π/3 for the floor 𝒢_c ≥ 8.4. The audit (T1 §3.3(a)) computed that radial
(n = r̂) and hyperbolic (n = (−x,−y,z)/r) hedgehogs are degenerate at 8π under one-constant
elasticity, concluded the band's lower edge encodes unstated anisotropy/escape structure, and
filed the band as unrecoverable-from-text. This run derives the three-constant closed forms and
decides what reaches 8π/3.

## 2. Closed forms (sympy, cross-checked by scipy quadrature to ≤ 5×10⁻¹⁶ rel.)

Standard Frank convention f = ½[K₁(∇·n)² + K₂(n·∇×n)² + K₃|n×(∇×n)|²]; textures scale-invariant,
so over a shell r ∈ [r_min, R] every term integrates to (coefficient)·K_i·(R − r_min):

| texture | splay | twist | bend | one-constant total | corpus's printed ∫K(∇n̂)² |
|---|---|---|---|---|---|
| radial n = r̂ | **8π K₁** | 0 | 0 | **8π K** | 8π K |
| hyperbolic n = (−x,−y,z)/r | **(8π/5) K₁** | **0 (exactly)** | **(16π/15) K₃** | **(8π/3) K** | 8π K |

- **E_hyp = (8π/15)(3K₁ + 2K₃)·r — twist-free.** K₂ never enters either texture (the scan is
  genuinely 2D in (K₁, K₃)).
- **One-constant degeneracy: REFUTED for the standard Frank form** (8π vs 8π/3, ratio exactly 3),
  **CONFIRMED for the corpus's printed unhalved Dirichlet form** ∫K(∇n̂)² (8π both — the audit's
  computation). The two statements are reconciled by the saddle-splay identity
  |∇n|² = S² + T² + B² + ∇·[(n·∇)n − n(∇·n)], verified pointwise-symbolically for both textures
  (residual ≡ 0); the divergence term is r-linear for scale-invariant textures and differs
  between them (+4π radial, −4π/3 hyperbolic, per unit r, halved convention), so it cannot be
  discarded as "just a surface term" in an r-linear energy.
- The azimuth-offset family n = (sinθcos(φ+φ₀), sinθsin(φ+φ₀), cosθ) interpolates the band
  continuously: **E(φ₀) = (8π/15)[K₁(7cos²φ₀+6cosφ₀+2) + 5K₂sin²φ₀ + K₃(2−cosφ₀−cos2φ₀)]**,
  = (8π/3)(2+cosφ₀) at one constant — from 8π (radial, φ₀=0) to 8π/3 (hyperbolic, φ₀=π).
  The degree −1 anti-hedgehog gives (8π/15)(3K₁+2K₃) for every φ₀ — the hyperbolic value.

## 3. What reaches 8π/3 (the scan; measured)

- **8π/3 is attained with NO anisotropy: it is the hyperbolic hedgehog's exact one-constant
  coefficient.** General locus (fixed-profile textures): (8π/15)(3k₁+2k₃) = 8π/3 ⟺
  **3k₁ + 2k₃ = 5** (k_i = K_i/K), which passes through isotropy k₁=k₃=1.
- The defect adopts the cheaper texture: hyperbolic iff k₃ < 6k₁ (radial-pure-splay wins only in
  strongly bend-stiff media).
- Floor-guarantee region (fixed profiles): c_min ≥ 8π/3 ⟺ k₁ ≥ 1/3 AND 3k₁+2k₃ ≥ 5.
- Eternal-threshold region: c_min ≥ 5.6 ⟺ k₁ ≥ 0.2228 AND 3k₁+2k₃ ≥ 3.3423 (uniform softening:
  k ≥ 0.6685).

## 4. New finding: 8π/3 is not a floor — the relaxed hyperbolic sits 8% lower (measured)

Minimizing over the profile ψ(θ) (ψ(0)=0, ψ(π)=π, i.e. fixed degree 1; 8→16→24 sine modes,
Gauss quadrature 400→1200 points, 12 random restarts, L-BFGS gradients ≤ 1.4×10⁻⁵):

- **min E(one-constant) = 7.7023807 = 2.451744π**, at φ₀ = π, stable to 8 digits under
  mode/quadrature refinement; profile monotone, max|ψ−θ| = 0.192 rad; split: splay 3.6009 +
  bend 4.1015, twist ≤ 3×10⁻¹⁶.
- Full-texture-space check (512² angular grid, evaluator validated to ≤ 5×10⁻⁶ rel. on both
  closed forms): 240 random smooth perturbations around the *ansatz* hyperbolic find descent
  (min ΔE = −5.5×10⁻³ — the ansatz is not stationary); the same 240-perturbation battery around
  the *relaxed* minimizer all raise the energy (min ΔE = +8.4×10⁻² > 0) — locally minimal.
- Relaxed three-constant scan (12×12 in (k₁,k₃) at k₂=1, warm-started): relaxation buys 0–57%
  over the fixed-profile ansatz (ratio 1.00 at bend-stiff corners where pure-splay radial is
  genuinely stationary; 0.43 at k₁=3, k₃=0.25 where the texture converts splay to cheap bend);
  at isotropy the ratio is 0.9194.
- Scope, printed: the relaxed value is the minimum over the m=1 hedgehog family plus a local
  full-space test; a *global* lower bound over all degree-1 textures is not proven here (no
  sharp topological lower bound for the Frank form is in the literature this run can lean on).

## 5. Reconciliation of audit vs corpus (analysis)

1. **The audit's GAP-(a) parenthetical is itself defective as printed.** Its step
   |∇(Mr̂)|² = |∇r̂|² is true but proves degeneracy of the *Dirichlet* energy only; the corpus's
   band matches the *standard Frank* convention, under which the textures split 8π vs 8π/3
   exactly. The audit graded the band "unrecoverable from the text"; it is recoverable in one
   line once the convention is fixed — [IM: standard defect elasticity] checks out against the
   standard convention.
2. **But the corpus's own printed functional contradicts its own band.** NR-D1b §1 prints
   E = ∫K(∇n̂)², under which c_h = 8π for *both* textures, no K₁/K₃ exists, and 8π/3 is
   unreachable. The corpus cannot have both the printed functional and the printed band.
3. **Direction analysis: every reading lands W-eternal at one constant.** Printed-functional
   reading: floor 8π ≈ 25.1 (clearance ×4.5). Standard-Frank ansatz reading: floor 8π/3 ≈ 8.38
   (×1.496, the corpus's number). Standard-Frank honest (relaxed-texture) reading: floor 7.702
   (×1.375). The corpus's choice of the band minimum was direction-safe but is *not* the
   texture floor under the convention that produces its own band.

## 6. Implication for the floor's c_h,min input (analysis; answers the audit's GAP)

- **The assumption that closes the band, pinned:** *the near-core energy is the standard
  three-term Frank form (the ½-normalized splay/twist/bend split, retaining the saddle-splay
  difference between textures), with the defect texture ranging over the hedgehog class
  {hyperbolic … radial}.* Nothing about anisotropy is needed; the band is the texture range at
  K₁=K₂=K₃=K. No anisotropy assumption *can* close it under the printed functional (no elastic
  constants survive in ∫K(∇n̂)²) — so the repair is the convention correction, full stop.
- **The honest floor under that convention is 7.702, not 8.378** (𝒢_c ≥ 7.7 at f·r_min = 1):
  the corpus's F-D1b-1 line "𝒢_c ≥ 8" overstates the true texture floor by 8.8%; W-eternal
  still lands with clearance ×1.375.
- **Anisotropy exposure, priced:** the floor stays above the eternal threshold 5.6 for
  3k₁+2k₃ ≥ 3.342 (fixed profile) / uniform softening k ≥ 0.727 (relaxed). A chiral medium
  whose splay/bend constants are soft relative to K = f² by ≳ 27–33% kills the eternal verdict.
  Note this is the *same* ~33% margin the audit found for the unflagged binding clause
  (both restate threshold/floor = 5.6/8.38-class): elastic softening, profile relaxation
  (8.8% of it already spent), and binding corrections all draw on **one** shared margin budget,
  which the corpus prices as if each were separately clear.

## 7. Deliverables

- `h23_frank.py` — derivation + scan (deterministic, seed 20260716; ~3 min single-threaded).
- `h23_results.json` — closed forms, quadrature checks, both scans, relaxed-minimum
  convergence/profile, perturbation tests.
- This memo.

— end of memo —
