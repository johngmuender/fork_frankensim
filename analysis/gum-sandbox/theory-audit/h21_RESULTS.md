# H2.1 — NR-K3a second moment: the decisive Monte-Carlo (T-H2)

**Campaign:** replication Phase H2, follow-on №6 of the theory audit (`T1_NR_gates.md` §4.2/§6).
**Date:** 2026-07-16. **Code:** `h21_mc.py` (python3 + numpy/scipy, seed 20260716, runtime 117 s).
**Raw output:** `h21_results.json`. **Figure:** `h21_fig.png`. **Status: measured facts below;
promotion of T-H2 is the coordinator's call.**

## 0. Verdict in one line

**The audit's second-moment estimate is CONFIRMED by direct Monte-Carlo: across the corpus's
own scale window, R = ⟨A²⟩/⟨A⟩² = 10¹⁴–10¹⁷ (grid corners 7.9×10¹² – 5.0×10¹⁸), never within
twelve orders of the corpus's implicit R ≈ 1; the audit's analytic 1/φ is accurate to better
than a factor 2.5 everywhere in its regime; and the motional-narrowing rescue requires the web
to translate ~0.3 of its own spacing inside one event's coherence time — a web speed of
~10⁸–10⁹ c. NR-K3a Step 2's "fluctuations negligible" clause fails as printed.**

## 1. The question (restated)

NR-K2b/K3a Step 2 (corpus): the 0νββ insertion amplitude is governed by the volume mean m̄ of
the line-supported pitch density — "line-support and bulk-support are amplitude-identical in
the mean"; the fluctuation term is "(fm/μm)-class, negligible". The audit (T-H2): the
observable is the *rate*, quadratic in the amplitude; for sparse line support the variance
dominates the squared mean by ~1/φ (φ = tube volume fraction) ~ 10¹⁴–2.5×10¹⁶. This run
decides between them by simulation.

## 2. Model and setup (measured facts about the code, not the corpus)

- **Web:** static network of straight tubes, radius `a`, uniform pitch amplitude m̂ inside,
  0 outside; volume mean m̄ = φ·m̂ held fixed (R is normalization-independent). Two
  geometries: (i) isotropic **Poisson line process**, length density ρ_L = 1/L_web²
  (audit convention, φ = πa²/L²); (ii) **jitter-perturbed cubic lattice** of lines along
  x/y/z at the same ρ_L (jitter ±0.25 of the lattice spacing).
- **Kernel:** the corpus prints no 0νββ kernel. Per the audit's reading (double-propagator,
  range ~1/q̄ ~ 2 fm) we use the normalized w(r) = exp(−r/r0)/(4π r0 r²), **r0 = 1 and 2 fm
  variants** — an acknowledged modeling choice, priced in §7.
- **Amplitude:** A(x0) = ∫ m(y) w(x0−y) d³y at random insertion points x0;
  R = ⟨A²⟩/⟨A⟩².
- **Scales, read from the corpus:** L_web = the μm spacing (NR-K3a §2); the heliknoton's
  pitch from m₃ = 0.047 eV gives λ = ħc/m₃ ≈ 4.2 μm (1σ 1.7–10 μm, Ω-paper VII.I; NR-D2
  uses p = 4 μm) → **L scan {1, 1.7, 4.2, 10} μm** (1 μm = the audit's conservative floor).
  Tube radius 1/f-class: f ∈ [4, 29] MeV (AUD-15 operative band) → 1/f ∈ [6.8, 49] fm;
  the audit used 52.6 fm → **a scan {2, 6.8, 20, 52.6, 200} fm** (2 = thin regime,
  200 = generously thick).
- **Method, three legs, each cross-checked:**
  1. one-tube response g(d) by deterministic quadrature, validated against brute-force
     kernel-sampling MC;
  2. exact Poisson second moment via **Campbell's theorem** (⟨A⟩ = m̂ρ_L∫2πd·g dd,
     Var = m̂²ρ_L∫2πd·g² dd) — exact at any scale separation, no MC noise;
  3. **direct network MC**: (a) at moderate separation (φ ~ 10⁻³–10⁻², reachable by uniform
     sampling of x0 in explicit Poisson and lattice networks) — validates leg 2 and the
     Poisson↔lattice equivalence; (b) at physical scales, direct MC with importance sampling
     on the nearest-tube distance (uniform sampling at φ ~ 10⁻¹⁴ would need >10¹⁴ points —
     that impossibility is itself the phenomenon under test).

## 3. Verification chain (all PASS)

**g(d) quadrature vs brute-force kernel MC** (4×10⁶ samples/point, 7 points spanning thin,
thick, inside/outside): all pulls within ±1.9σ; e.g. a=2, r0=1, d=4: quad 0.010598 vs MC
0.010664 ± 5.1×10⁻⁵. (The a=52.6, d=20 row is g = 1 with zero MC variance — its "pull" in
the JSON is an artifact of SE = 0; quadrature error there is 6×10⁻⁵ relative.)

**Mean check (the corpus's Step-2 mean statement):** ⟨A⟩/m̄ = 1.0000 ± 0.0004 on every grid
row (quadrature) and 0.996–1.003 (physical-scale MC). **The corpus's first-moment claim is
verified — line- and bulk-support ARE amplitude-identical in the mean.** The defect is
entirely in the second moment.

**Direct network MC vs Campbell, moderate separation** (r0 = 1; 4 network realizations ×
7.5×10⁵ points each; blocks give the error bar):

| a | φ | R Campbell | R MC Poisson | R MC jittered lattice |
|---|---|---|---|---|
| 0.5 | 1.5×10⁻³ | 140.4 | 137.3 ± 3.6 | 141.0 ± 0.3 |
| 2 | 3.0×10⁻³ | 212.3 | 207.4 ± 9.7 | 212.5 ± 2.3 |
| 10 | 1.0×10⁻² | 93.1 | 93.2 ± 3.8 | 91.9 ± 0.6 |

Poisson and perturbed-lattice geometries agree with each other and with the exact Campbell
value within error bars → the Campbell numbers at physical scales are trustworthy for both
geometry classes, and R is not an artifact of the Poisson idealization.

**Direct MC at physical scales vs Campbell** (2×10⁶ importance-sampled environments/point,
jackknife errors):

| a [fm] | r0 [fm] | L [μm] | R direct-MC | R Campbell |
|---|---|---|---|---|
| 52.6 | 1 | 1.0 | (1.133 ± 0.002)×10¹⁴ | 1.133×10¹⁴ |
| 52.6 | 2 | 1.0 | (1.115 ± 0.003)×10¹⁴ | 1.116×10¹⁴ |
| 52.6 | 2 | 4.2 | (1.967 ± 0.005)×10¹⁵ | 1.968×10¹⁵ |
| 6.8 | 2 | 1.0 | (5.319 ± 0.008)×10¹⁵ | 5.314×10¹⁵ |
| 2.0 | 1 | 4.2 | (8.868 ± 0.017)×10¹⁷ | 8.900×10¹⁷ |
| 200 | 2 | 10.0 | (7.925 ± 0.018)×10¹⁴ | 7.894×10¹⁴ |

MC error bars are 0.2–0.3% — the verdict is unambiguous not just in log10 but at the
percent level.

## 4. The R table (main result; Campbell-exact, MC-anchored per §3)

R = ⟨A²⟩/⟨A⟩², kernel r0 = 2 fm (r0 = 1 fm in parentheses), across the (a, L_web) grid:

| a \ L_web | 1 μm | 1.7 μm | **4.2 μm (corpus central)** | 10 μm |
|---|---|---|---|---|
| 2 fm | 3.2×10¹⁶ (5.0×10¹⁶) | 9.4×10¹⁶ (1.5×10¹⁷) | 5.7×10¹⁷ (8.9×10¹⁷) | 3.2×10¹⁸ (5.0×10¹⁸) |
| **6.8 fm (1/f, f=29 MeV)** | 5.3×10¹⁵ (6.1×10¹⁵) | 1.5×10¹⁶ (1.8×10¹⁶) | **9.4×10¹⁶** (1.1×10¹⁷) | 5.3×10¹⁷ (6.1×10¹⁷) |
| 20 fm | 7.3×10¹⁴ (7.6×10¹⁴) | 2.1×10¹⁵ (2.2×10¹⁵) | 1.3×10¹⁶ (1.3×10¹⁶) | 7.3×10¹⁶ (7.6×10¹⁶) |
| **52.6 fm (audit's 1/f)** | **1.1×10¹⁴** (1.1×10¹⁴) | 3.2×10¹⁴ (3.3×10¹⁴) | **2.0×10¹⁵** (2.0×10¹⁵) | 1.1×10¹⁶ (1.1×10¹⁶) |
| 200 fm | 7.9×10¹² (7.9×10¹²) | 2.3×10¹³ (2.3×10¹³) | 1.4×10¹⁴ (1.4×10¹⁴) | 7.9×10¹⁴ (7.9×10¹⁴) |

- **Minimum over the entire scanned window: R = 7.9×10¹²** (a = 200 fm — 4–30× thicker than
  the corpus's 1/f-class tube — at the audit's conservative L = 1 μm).
- **Corpus-central cells: R ≈ 10¹⁵–10¹⁷** (a = 1/f-class, L = 4.2 μm).
- The corpus's Step-2 claim corresponds to R ≈ 1. **The measured gap is 13–18 orders of
  magnitude; the "negligible" fluctuation term is the dominant term by that margin.**

## 5. Analytic vs MC (the audit's estimate, priced)

- **Thick regime (a ≫ r0):** audit's R ≈ 1/φ. Measured R·φ = 0.97–0.99 at a = 52.6 fm,
  0.92–0.96 at 20 fm, 0.77–0.88 at 6.8 fm, 0.99 at 200 fm — **the 1/φ law is correct to
  ≤ 25% in its stated regime**, and never off by more than a factor 2.5 anywhere
  (worst: R·φ = 0.41 at a = 2 fm, r0 = 2 fm, where a < r0 and 1/φ is not the right law).
- **Thin regime (a ≪ r0):** audit's R ≈ L²/(π²r0²) (2.5×10¹⁶ at L = 1 μm, r0 = 2 fm).
  Measured at a = 2 fm: R = 3.2×10¹⁶ — the formula is right in order; the exact thin-tube
  answer carries the expected log(r0/a)-class correction (measured R exceeds it by ×1.3
  and falls below the r0 = 1 fm version by ×2, i.e. the audit's thin formula is an
  envelope, good to a factor ~2).
- The audit's quoted band "10¹⁴–2.5×10¹⁶" is exactly reproduced at its own reference scales
  (a = 52.6 fm, L = 1 μm → 1.1×10¹⁴; thin limit → 3×10¹⁶). At the corpus's *central* spacing
  (4.2 μm) the enhancement is larger: **10¹⁵–10¹⁸.**
- Physical meaning, confirmed by the simulation: a φ-fraction of insertion points sit on
  tubes and carry A ≈ m̂ = m̄/φ; the rest carry ≈ 0. Rate ∝ ⟨A²⟩ = m̄²/φ, not m̄².

## 6. The motional-narrowing rescue, priced (the discharge condition)

The corpus text gives the web no dynamics; a rescue would need the web to move enough
*within a single event's coherence time* that the kernel averages over tube positions.
Model: rigid random transverse displacement, 2D Gaussian rms δ per axis, time-averaged
kernel g_eff = g ⊛ Gauss(δ); R(δ) recomputed exactly (mass check ⟨A⟩ = m̄ preserved to
≤ 4×10⁻⁴ at every δ).

Measured crossings (insensitive to a and r0, as the asymptotic R−1 ≈ L²/(4πδ²) predicts):

| target | required δ/L_web |
|---|---|
| R = 10 | 0.094–0.095 |
| R = 2 | 0.29–0.30 |

**The discharge condition: δ ≈ 0.3 L_web per coherence time.** With τ_coh ≈ r0/c
(= 6.7×10⁻²⁴ s at r0 = 2 fm) this is a web transport speed of

- L = 4.2 μm, r0 = 2 fm: **v ≈ 6.1×10⁸ c** (R = 2); 2.0×10⁸ c (R = 10)
- L = 1 μm, r0 = 2 fm: v ≈ 1.5×10⁸ c
- L = 4.2 μm, r0 = 1 fm: v ≈ 1.2×10⁹ c

i.e. **superluminal by 8–9 orders of magnitude.** Equivalently: to rescue the mean-field
answer the web must traverse ~μm distances in ~10⁻²³ s. Slow web motion (any v ≪ 10⁸ c)
leaves each decay event sampling a frozen snapshot, and averaging |A|² over snapshots
reproduces the static R. This is a *pricing* of the rescue class the audit named, not a
proof that no other rescue exists (see caveats).

## 7. Honest model caveats (all four flagged before running; none removed by the run)

1. **Straight static tubes, uniform m̂.** Real network curvature/junctions enter only at
   μm scale ≫ kernel range and cannot change R by more than O(1); tested Poisson vs
   perturbed-lattice — geometry class shifts R by < 3%. Non-uniform |m| along tubes could
   only *increase* the second moment (Cauchy–Schwarz direction).
2. **Insertion linear in local m(y)** — the corpus's own model (m̃(k) superposition,
   NR-K3a §2). A saturation/nonlinearity of the insertion on tubes is the one in-model
   rescue this MC cannot price; it is nowhere in the corpus text.
3. **Kernel choice.** The corpus never prints the kernel; we used exp(−r/r0)/r² with
   r0 = 1, 2 fm (standard 0νββ virtuality ~100 MeV). R depends on the kernel only through
   its range where a ≲ r0; in the corpus-central thick regime (a = 1/f ≫ r0) R = 1/φ is
   kernel-blind — the r0 = 1 vs 2 fm columns differ by < 1% there. A hypothetical kernel
   with range ≳ L_web (~μm, i.e. q̄ ~ 0.05 eV) would erase R, but that contradicts the
   fm-scale virtuality the corpus itself uses ("the fm-scale virtual exchange", NR-K3a §2).
4. **A smooth bulk component of m(y)** carrying an O(1) fraction of m̄ would restore
   R → O(1); the corpus's Theorem K-2 asserts the opposite ("the mixing density m(z) is
   line-supported, by the same theorem that licensed it"), so this rescue would cost the
   corpus its own Step-1 premise.
5. What this run does **not** adjudicate: whether the resulting rate *enhancement* (m̄²/φ)
   or a proper re-derivation of the four-cell table is phenomenologically fatal — that is
   the coordinator's promotion question (T-H2's "OPPOSITE direction from branch (b)" point
   stands: both original fork branches are unsupported).

## 8. Deliverables

- `h21_mc.py` — the simulation (deterministic, seed 20260716; ~2 min single-threaded).
- `h21_results.json` — full grids, MC error bars, validation pulls, narrowing curves.
- `h21_fig.py` / `h21_fig.png` — panel A: R vs a across L with MC points over
  Campbell-exact curves and the audit's 1/φ dashed; panel B: R(δ/L) rescue curve.
- This memo.

— end of memo —
