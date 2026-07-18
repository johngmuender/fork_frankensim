# M2 — the bootstrap servo coefficient (IV.H.3 ⟨r2⟩ toy law) [F-T7-M2]

**Goal.** Independently measure the corpus's UNREPLICATED toy law for the
phase-locked knot self-oscillator: relaxation rate
λ = (0.021 ± 0.004) g_χ²ω\* in the thermalizing-bath class (IV.H.3 ⟨r2⟩),
marginal in the integrable limit. h27 covered the ±½ exponents, the
−sin Δφ torque and the γ = 0 class boundary; this workstream covers the
LAW — the g² scaling, the ω\* linearity, and the coefficient — on an
independent continuum-bath toy. Last unreplicated [CAL] anchor of the
quantum-adjacent core.

**Within-model; nothing here bears on nature.**

Reproduce with `python3 m2_servo.py` (deterministic seeds; 257 s wall;
writes `m2_results.json`; console log in `m2_run.log`) then
`python3 m2_fig.py` (writes `m2_fig.png`).
Units ħ = m = 1; ω\* = 1 except where rescaled.

---

## 1. The reduced equations (derived before any numbers)

**The pair.** The knot clock φ_mech and the medium consensus φ_phase see
the displaced order parameter 𝔥 = 𝔥\*(1+δ) through the h27-confirmed exact
exponents d ln ω_mech/d ln 𝔥 = +½, d ln ω_phase/d ln 𝔥 = −½. The phase
mismatch Δφ = φ_mech − φ_phase therefore obeys **exactly**

    Δφ̇ = ω*[(1+δ)^{+1/2} − (1+δ)^{−1/2}] ≡ ω* f(δ)   (→ ω* δ, small δ).

The locking torque acts back on the order-parameter displacement with the
h27-confirmed Adler sign; the corpus's servo is second order in the
knot–medium coupling, so

    δ̇ = −g_χ² ω* sin Δφ + F_bath(t).

With A = ω\* (the only frequency scale — this is what G3 tests) the
deterministic pair is **canonical**: taking q = Δφ, p = δ and

    H_S = ω* U(δ) + g_χ²ω*(1 − cos Δφ),  U′(δ) = f(δ),  U(δ) ≈ δ²/2,

Hamilton's equations reproduce both lines. The deterministic system is a
conservative pendulum with libration frequency **Ω_lib = g_χ ω\***:
closed orbits, no attractor — the integrable limit is marginal by
construction (h27's γ = 0 class), and the bath must supply ALL of the
irreversibility.

**The bath (stated choice).** N harmonic modes bilinearly coupled to the
angle Δφ, standard Caldeira–Leggett form with counterterm,

    H = H_S + Σ_k [p_k²/2 + ω_k²x_k²/2] − Δφ Σ_k c_k x_k + Δφ² Σ_k c_k²/(2ω_k²),

so the bath force lands in the δ̇ equation as specified. Both channels of
the knot–medium coupling carry g_χ: the coherent torque is O(g_χ²) and the
couplings c_k ∝ g_χ, hence the friction kernel is O(g_χ²) too —

    J(ω) = (π/2) Σ_k (c_k²/ω_k) δ(ω−ω_k),
    Ohmic:       J(ω) = η  g_χ² ω             e^{−ω/ω_c}
    super-Ohmic: J(ω) = η₃ g_χ² (ω³/ω_c²)     e^{−ω/ω_c},   ω_c = 5ω*.

η (η₃) is the dimensionless bath-convention constant — the one number the
unreleased r2 spec would have to fix.

**The law, derived.** Integrating the bath out gives the generalized
Langevin equation δ̇ = −g²ω* sin Δφ − ∫γ(t−s) Δφ̇(s) ds + ξ(t) with
⟨ξ(t)ξ(0)⟩ = T γ(t) (classical FDT). Linearizing about the lock and
evaluating the kernel at the libration frequency:

    λ = ω* J(Ω_lib)/(2 Ω_lib),   Ω_lib = g_χ ω*.

- **Ohmic:** λ = (η/2) g_χ² ω\* e^{−g_χω*/ω_c} — the corpus's law, with
  **c = λ/(g_χ²ω\*) = η/2 exactly** (up to the spec-fixed cutoff factor):
  the g² comes from c_k² ∝ g_χ², the ω\* from Ω_lib ∝ ω\* against the
  scale-free Ohmic kernel, and the coefficient is **pure bath
  convention**.
- **super-Ohmic:** λ = (η₃/50) g_χ⁴ ω\* e^{−g_χ/5} — **the g² law itself
  fails off the Ohmic class** (J/ω is no longer flat at Ω_lib ∝ g), so
  the corpus's printed g² scaling, if measured, pins its bath to the
  Ohmic (thermalizing, J ∝ ω) class.

Beyond the weak-coupling formula, the exact linear-response rate is the
complex root of D(s) = s² + Ω_lib² + ω\* s² K(s), K(s) = (2/π)∫J(ω)/(ω(s²+ω²))dω
(computed numerically as `lam_lin_exact`; it carries the cutoff and the
counterterm-included frequency renormalization, which is O(1%) for the
Ohmic scan but sizable for the super-Ohmic convention used here).

**Numerics.** Velocity-Verlet (symplectic) on the full closed system,
dt = 0.02/ω\*; N = 2000 modes (8000 where the resonance linewidth needs
them), log grid ω ∈ [10⁻³, 25ω\*], one mode per cell, J-locked couplings
c_k² = (2/π)ω_k J(ω_k)Δ_k, frequencies jittered per realization;
R = 24 realizations (12 for the heavy baths) of thermal initial data at
T = 10⁻⁷ω\*; δ(0) = 0.05, Δφ(0) = 0. Envelope rate fit on the
ripple-reduced envelope Q = E_S + λ Δφ δ (exact e^{−2λt} envelope of the
damped linear pair), iterated window [0.2, 2.2]/λ, block errors over
realizations ⊕ window-variation systematic. RK4 for the deterministic
control. Resolution guard: modes inside the linewidth λ·ρ(Ω_lib) ≥ ~3
everywhere (the two-mode toy of h27 fails the g² law precisely because
this number is ≪ 1 there — resolved: see §6).

## 2. Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | integrator: deterministic-limit energy drift ≤ 10⁻⁶; bath statistics Gaussian | RK4 pair, 2×10⁵ steps: max \|ΔE/E\| = 2.5×10⁻¹³, secular −2.0×10⁻¹³. Bath noise ξ (M = 5000): skew = +0.007 ± 0.035, excess kurtosis = +0.069 ± 0.069, KS p = 0.35 (t = 0) / 0.83 (t = 2); classical FDT ⟨ξ(t)ξ(0)⟩ = Tγ(t) to max dev 2.6% of peak vs 2.0% expected sampling error; γ_disc(0)/γ_cont(0) = 0.993. Bath runs (symplectic Verlet, closed system): max \|ΔH/H\| ≤ 1.4×10⁻⁴ bounded, secular ≤ 3.8×10⁻⁷ | **PASS** |
| G2 | λ ∝ g_χ² across ≥ 1 decade; log-slope 2.0 ± 0.1 | g_χ ∈ [0.06, 0.6] (exactly one decade, 8 points): **slope = 1.9932 ± 0.0041**. Per-point λ agrees with the exact linear-response kernel rate to ≤ 1.5% (mean \|dev\| 0.7%); the −0.007 offset from exact 2 is the spec-fixed cutoff drift (d ln λ/d ln g = 2 − g ω\*/ω_c) | **PASS** |
| G3 | λ/ω\* invariant under ω\* → 4ω\* (± 10%) | at g_χ = 0.2237 (dt, ω_min held in ABSOLUTE units so the test has content): (λ/ω\*)₄/(λ/ω\*)₁ = **1.0050 ± 0.0074**. dt cross-check 0.02 → 0.005: −0.6% | **PASS** |
| G4 | c ≡ λ/(g_χ²ω\*) vs 0.021 ± 0.004; within 2σ_comb = PASS-MATCH, else PASS-MEASURED + ≥ 2 bath shapes | **c = 0.4917 ± 0.0051** (Ohmic, η = 1; per-point spread 0.486–0.502; cutoff-corrected 0.517 ≈ η/2). z = \|c − 0.021\|/σ_comb = **72.9** → not a match. Shape/convention spread: c = η/2 exactly (η = 0.1 gives c = 0.0478, i.e. c/η = 0.478 ≈ 0.492 = c/η at η = 1); super-Ohmic (η₃ = 6) c = 0.0115 → 0.0269 across g = 0.35 → 0.65, g-dependent (no law). Corpus value reproduced in-class by η = 0.042 (verified: c = 0.0201) | **PASS-MEASURED** |
| G5 | integrable limit (bath off) marginal: fitted \|λ\| consistent with 0 at fit-noise floor | same pipeline, bath off, g_χ = 0.2: λ = (8.1 ± 4.1)×10⁻¹³ (RK4, matched fit window; whole-run 5×10⁻¹⁷), λ = (4.8 ± 9.3)×10⁻¹⁰ (Verlet) — consistent with 0 at the floor, **10–11 orders below** the bath-on rate 1.9×10⁻² at the same g. h27's γ = 0 class boundary reproduced | **PASS** |

## 3. Key numbers

- **The law holds in the thermalizing-bath class: λ = c g_χ² ω\*** with
  slope-in-g = 1.9932 ± 0.0041 over one decade and ω\*-linearity to
  0.50% ± 0.74% under a ×4 rescale — the two structural claims of ⟨r2⟩
  are REPLICATED.
- **The coefficient is NOT structural: c = η/2 exactly** (measured
  0.4917 ± 0.0051 at η = 1 vs 0.4894 exact-linear / 0.4761 weak-law
  with cutoff; measured 0.0478 at η = 0.1 vs 0.0481
  exact-linear). It is linear in the one free
  bath-convention constant. The corpus's 0.021 ± 0.004 corresponds to
  η ≈ 0.042 — an ordinary weak Ohmic coupling, attainable in-class,
  unverifiable without r2's bath spec.
- **The g² law selects the Ohmic class**: super-Ohmic J ∝ ω³ gives
  measured slope 3.38 (exact-linear 3.34; ideal 4 − g/5 minus
  renormalization) — λ ∝ g⁴-class, and c becomes g-dependent
  (0.0115 → 0.0269). A printed constant-c g²-law is itself evidence the
  corpus bath is Ohmic/thermalizing, exactly the class it names.
- Measured vs exact-linear-response rate across ALL 13 bath configs:
  agreement ≤ 1.6%.
- Integrable control: \|λ\| ≤ 8×10⁻¹³, envelope flat (marginal), while
  the torque and exponents are unchanged — irreversibility lives
  entirely in the bath, as h27 pinned.

## 4. Method notes

- Bath discretization: log grid, N = 2000 (scan) / 8000 (η = 0.1 and
  super-Ohmic, where the linewidth is narrow), one mode per cell,
  per-realization frequency jitter; Poincaré time 2πρ(Ω_lib) ≥ 5× every
  run length; modes inside the linewidth λρ(Ω_lib) ≈ 3–20.
- R = 24 (scan) / 12 (heavy) realizations; T = 10⁻⁷ ω\* (envelope
  dynamic range ~4.7 amplitude e-folds above the thermal floor;
  floor subtracted in the fit).
- Errors: block-over-realizations statistical ⊕ fit-window-variation
  systematic (windows ×0.7 / ×1.4); the window systematic dominates.
- σ_comb for G4 = √(σ_meas² + σ_corpus²) = 0.0065.
- The exact linear-response check (`lam_lin_exact`) is an independent
  prediction of the same model (continuum kernel, Newton root of
  D(s) = s² + Ω² + ω\*s²K(s)) — not a fit to the data.
- Wall time 257 s single-threaded; fully seeded.

## 5. Caveats

- **The coefficient comparison is convention-bounded, not physical**:
  r2's bath spec is unreleased, and in this architecture
  c = J(Ω_lib)/(2 g_χ²Ω_lib) exactly — a functional of the bath
  spectrum alone. Our η = 1 declaration is a convention; the honest
  content of G4 is the c = η/2 law + the demonstration that 0.021 is
  in-class attainable (η = 0.042) and that the super-Ohmic convention
  sweeps c across the corpus band with a broken g-law.
- The reduced pair is the corpus's own architecture (exact ±½ exponents,
  −sin Δφ torque) but the O(g²) assignment of BOTH channels (torque and
  friction) is our modeling choice; it is the unique assignment that can
  produce λ ∝ g² from a bilinear bath.
- Small-g end of the scan starts at 17% of the pendulum separatrix
  energy (δ₀/2g = 0.42 at g = 0.06): nonlinear frequency shift ~4%,
  visible only as the ≤ 1.5% λ deviations already absorbed in the
  exact-linear comparison.
- Classical bath (thermal ICs, classical FDT), T small; quantum bath
  statistics untested — at T → 0 a quantum Ohmic bath would add
  zero-point noise but the mean-envelope decay rate of the linearized
  pair is T-independent, so the gates should be insensitive; not
  verified here.
- The raw G2 slope carries the spec-fixed cutoff drift 2 − gω\*/ω_c;
  at ω_c = 5ω\* this is −0.007 on the decade fit (measured), i.e. the
  gate would fail for scans pushed to g ~ O(1), where λ/Ω_lib is also
  no longer small.

## 6. h27 reconciliation

h27's two-mode toy measured λ ∝ γ (bath-limited, g-blind) and could not
reach the g² law; it conjectured the law lives in the
coupling-limited/continuum-bath regime. Confirmed and quantified here:
with a continuum Ohmic bath (many modes inside the servo linewidth,
λρ(Ω_lib) ≳ 4) the rate becomes bath-spectrum-limited at the libration
frequency, λ = ω\*J(Ω_lib)/(2Ω_lib) ∝ g², independent of any per-mode
γ. The two regimes are the two sides of the same kernel: h27 sat at
λρ(Ω) ≪ 1 (energy shuttles into one resonant mode and back at that
mode's own decay rate); ⟨r2⟩'s law requires λρ(Ω) ≫ 1. The class
boundary (γ = 0 ⇔ bath off ⇔ marginal) is identical in both.

## Files

- `m2_servo.py` — model, integrators, fits, campaign driver.
- `m2_results.json` — every measured number.
- `m2_fig.py`, `m2_fig.png` — scaling law, envelopes, coefficient spread, FDT.
