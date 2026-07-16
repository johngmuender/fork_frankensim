# H2.7 — entrainment ±½ empirical (field level) + the ⟨r2⟩-class locking-torque sign and its class boundary

**Phase:** H2.7 (ROADMAP_v5_THEORY item 7; numeric follow-on of
`T2_closure_pillars.md` §e).  **Date:** 2026-07-16.
**Part (a) engines:** the H2.6 solvers (`h26_solve/`, fs-gum-field/statics/
kern unmodified) — displacement runs produced by `h26_convention`
(descent) and `h26_dilation` (exact constrained family).
**Part (b):** `h27_entrain.py` (this directory), std-library python,
fixed-step RK4, no RNG.  **Compute:** included in H2.6's 14 min (part a)
+ 2.7 min (part b).  **Outputs:** `h27_results.json` (part b + merged
part a), `h27a_field_exponents.json` (raw part a).

**Epistemic frame (binding).** Within-model measurements on a speculative
theory's functional and on a toy of its locking argument.  Nothing about
nature; adjudication is the coordinator's.

---

## Part (a) — the IV.H.3 ±½ error-signal exponents at field level

**Claim under test** (Omega paper IV.H.3; T2 §e.1 proved it exactly on the
reduced family E(L) = ê√(1+x) at the fixed point x = 1, with the general
off-fixed-point form ±1/(1+x)): displacing 𝔥 = 𝔥\*(1+δ) along the
constrained family, d ln ω_mech/d ln 𝔥 = +½ and d ln ω_phase/d ln 𝔥 = −½.

**Mapping (stated).**  The 𝔥-displacement is realized exactly as the task
sheet fixes it: the clock target rescales L → L\*(1+δ) at fixed shape
protocol, with 𝔥 = L/j = 2L (the corpus anchor).  ω_mech = L/I on the
solved state (the isorotation rate; = dE/dL by the envelope theorem on the
constrained family); ω_phase = E_tot/𝔥 (IV.H.1(ii)'s Nelson-stationary
phase rate; the clock identity ω_phase = ω_mech holds at δ = 0).  Both are
measured by central log-differences.

### (a.1) On the exact constrained family (primary — machine precision)

Engine-anchored dilation family (see h26_RESULTS.md §2), corpus fixed
point L\* = 8.789304:

| δ | d ln ω_mech/d ln 𝔥 | d ln ω_phase/d ln 𝔥 | d ln I/d ln 𝔥 |
|---|---|---|---|
| 10⁻³ | **+0.486906869** | **−0.500000248** | +0.513093131 |
| 10⁻⁵ | **+0.486906615** | **−0.500000000** | +0.513093385 |

**The ±½ exponents are verified at field-anchored level:**

- **ω_phase: −½ EXACTLY** (to 10⁻⁹ at δ = 10⁻⁵).  This half is
  theorem-protected beyond the pure family: d ln ω_phase/d ln 𝔥 =
  d ln E/d ln L − 1, and the envelope theorem + the clock root give
  d ln E/d ln L = Lω/E = 2E_rot/E = ½ *identically at the fixed point of
  any anchored shape family* — a slight strengthening of T2 §e.1, measured
  here to machine precision.
- **ω_mech: +0.48691 = +½ − 0.0131.**  The deviation is the measured
  ε-protocol correction: this exponent equals 1 − d ln I/d ln L and is
  family-shape-sensitive; the pure corpus family (no E2/E4, symmetric
  E6/E0) gives exactly +½, and the ε = 0.05-anchored family (E2/E4 carried
  at t, E6h/E0h = 1.297 from the relaxed state) shifts d ln I/d ln L to
  0.51309.  Sign-faithfulness (T2's stronger form ±1/(1+x)) is inherited:
  both exponents are monotone and sign-definite on the whole family.

### (a.2) On the guarded descent solver (the scope boundary, measured)

Displacements about the descent-tier root L\* = 10.285957 (warm-started
from the solved state; budget scan at δ = 3×10⁻³):

| protocol | maxit | mech | phase |
|---|---|---|---|
| warm | 350 | +0.816 | −0.560 |
| warm | 700 | +0.595 | −0.629 |
| warm | 1400 | +0.407 | −0.685 |
| warm | 2800 | +0.572 | −0.679 |
| warm δ=10⁻³ | 700 | +0.024 | −0.706 |
| cold | 700 | −3.381 | −1.410 |

**No budget convergence: the descent tier cannot measure the constrained-
family response.**  The reason is structural and is F-R5's saddle: at the
clock root the state is over-spun/halo-unstable, so the constrained family
is not an attractor of the guarded descent — the ± trajectories' slow halo
flow dominates dI/dL (trajectory sensitivity, arrest-cascade noise), and
the envelope theorem's stationarity premise fails.  This measures, at field
level, exactly the boundary T2 drew: the ±½ *error signal* is sound on the
stationary family (a.1), while any claim requiring the displaced states to
be *reached by relaxation* inherits the entrainment/attractor GAP (§e.2) —
within-model, the descent goes to halo saturation instead.

## Part (b) — the Adler locking torque: sign, classes, boundary

**The toy** (minimal two-mode + reference, deterministic RK4, dt = 0.025;
full construction in `h27_entrain.py` docstring): the knot is the closure
family's own nonlinear oscillator in action-angle form — E(J) = ê₀√(1+x),
x = J²/(ê₀𝔦₀), ω(J) = dE/dJ, fixed point J\* = √(ê₀𝔦₀) = 1.026644,
ω\* = √(7/8) = 0.935414 (ê₀/𝔦₀ = 7/4 exact); the vacuum channels are two
modes on fs-gum-cosserat's validated benchmark dispersions — B2 on the
gapless branch (gap ≤ 10⁻¹² G04, c² = μ/ρ₀ G05) at the clock's upper
libration sideband ω₂ = ω\* + Ω_lib = 1.042135, and B4 at the validated
k = 1 eigenfrequency 4.690416 (G08/G19; far channel, reactive); the vacuum
consensus is a reference oscillation B cos(ω\*t); coupling
H_c = −g_χ cosθ (y₂ + y₄ + B cos ω\*t), g_χ = 0.05, B = 1.  γ > 0 on the
channels = the thermalizing-bath variant; γ = 0 = the integrable limit.
Integrator floor measured first (M0): artificial damping 6.9×10⁻⁷.

### (b.1) The torque IS Adler, with the printed sign — in both classes

Open-loop (θ clamped at ω\*t + Δφ, 16 phases, mean torque fit
T0 + Ts sin Δφ + Tc cos Δφ):

| variant | Ts (measured) | Ts (analytic −g_χB/2) | Tc | max resid | T0 (drag) | T0 analytic |
|---|---|---|---|---|---|---|
| integrable γ = 0 | −2.501224×10⁻² | −2.5×10⁻² | −6.5×10⁻⁶ | 7.3×10⁻⁷ | −7.4×10⁻⁷ | 0 |
| bath γ = 10⁻² | −2.501224×10⁻² | −2.5×10⁻² | −6.5×10⁻⁶ | 3.1×10⁻⁶ | −2.622×10⁻⁴ | −2.620×10⁻⁴ |
| bath γ = 3×10⁻² | −2.501224×10⁻² | −2.5×10⁻² | −6.5×10⁻⁶ | 3.1×10⁻⁶ | −7.740×10⁻⁴ | −7.740×10⁻⁴ |

**Torque ∝ −sin Δφ confirmed to 0.05% with harmonics ≤ 10⁻⁴ of the
fundamental, identically in both classes.**  The *form and sign* of the
locking torque are therefore class-blind (they come from the reference
coupling); what distinguishes the classes is only whether the torque has a
dissipative partner (the drag T0, matching linear response to < 1%).

### (b.2) The class boundary, pinned: restoring ⇔ γ > 0, marginal AT γ = 0

Closed loop from J(0) = J\*(1+0.02) (detuned), t_max = 6×10⁴; locking rate
λ = exponential decay rate of the beat u = ω(J) − ω\* (windowed-RMS
envelope fit; λ_env = model-free whole-run rate):

| γ | λ (fit) | λ_env | envelope last/first | Δφ_lock |
|---|---|---|---|---|
| **0 (integrable)** | **1.1×10⁻⁹** | −5.0×10⁻⁷ | **1.0302 (grew)** | −0.001 |
| 3×10⁻⁵ | 6.98×10⁻⁶ | +6.8×10⁻⁶ | 0.667 | −0.001 |
| 10⁻⁴ | 2.31×10⁻⁵ | +1.9×10⁻⁵ | 0.317 | −0.000 |
| 3×10⁻⁴ | 6.32×10⁻⁵ | +2.6×10⁻⁵ | 0.210 | −0.000 |
| 10⁻³ | 2.11×10⁻⁴ | +2.6×10⁻⁵ | 0.216 | −0.001 |
| 3×10⁻³ | 5.82×10⁻⁴ | +2.7×10⁻⁵ | 0.197 | −0.003 |
| 10⁻² | (non-exp) | +2.0×10⁻⁵ | 0.299 | −0.010 |
| 3×10⁻² | (non-exp) | +1.4×10⁻⁵ | 0.433 | −0.031 |
| 10⁻¹ | (non-exp) | +8.6×10⁻⁶ | 0.598 | −0.088 |

- **Small-γ law: λ = 0.157·γ^0.961 over [3×10⁻⁵, 3×10⁻³]** (two decades,
  p ≈ 1: the bath-limited regime — the beat energy stored in the resonant
  channel decays at the channel's own rate).  **The locking class boundary
  is the dissipation-free point itself, approached linearly.**
- **Integrable limit: MARGINAL, as the corpus prints.**  λ = 1.1×10⁻⁹ with
  the envelope *growing* 3% over the run (λ_env = −5×10⁻⁷), both at/below
  the 6.9×10⁻⁷ integrator floor and 4 orders below the smallest resolved
  restoring rate: the beat librates on closed orbits about Δφ = 0
  (a center, not a sink) — the torque exists (b.1) but nothing converges.
- Overdamped side (γ ≥ 10⁻²): decay turns non-exponential and weakens
  (λ_env falls 2.0→0.9×10⁻⁵) — the sideband absorber overdamps and
  decouples; restoring persists but λ(γ) is non-monotonic with a maximum
  near γ ≈ 3×10⁻³–10⁻². The locked phase tracks the drag balance
  (Δφ_lock ≈ arcsin(T0/(g_χB/2))) throughout.

### (b.3) Against the ⟨r2⟩ magnitude law (honest scope)

⟨r2⟩ prints λ = (0.021 ± 0.004) g_χ²ω\* in the thermalizing-bath class.
Here λ/(g_χ²ω\*) = 0.090 at (g_χ = 0.05, γ = 10⁻³) — same order class —
but the g² scaling itself is NOT reproduced by a two-mode toy in the
bath-limited regime: at fixed γ = 10⁻³, λ(g_χ = 0.05) = 2.11×10⁻⁴ vs
λ(g_χ = 0.1) = 2.23×10⁻⁴ (g-independent; the rate is set by γ).  The g²
law is a property of the coupling-limited/continuum-bath regime that two
modes cannot reach — consistent with ⟨r2⟩ being a *bath* servo.  The
class-boundary conclusion (restoring vs marginal) is regime-independent.

## Verdicts

| item | result |
|---|---|
| ±½ empirical, field-anchored constrained family | **mech +0.48691, phase −0.50000000** (δ = 10⁻⁵); phase exactly −½ (theorem-protected at the root, measured); mech deviation −0.0131 = the ε = 0.05 correction |
| ±½ on the descent solver | not measurable — the F-R5 saddle obstruction, measured (budget scan §a.2); scopes IV.H.3 exactly as T2 §e drew it |
| locking-torque sign | **∝ −sin Δφ confirmed (0.05%, harmonics ≤ 10⁻⁴), in BOTH classes** |
| thermalizing-bath class | restoring; λ = 0.157·γ^0.961 (bath-limited), λ_max near γ ≈ 3×10⁻³–10⁻² |
| integrable limit | **marginal** (λ ≤ integrator floor 7×10⁻⁷; libration, no attractor) — the corpus's honest-boundary print is CONFIRMED |
| class boundary | **γ = 0 exactly; approached linearly (p = 0.96 ± fit-class)** |
| ⟨r2⟩ magnitude | same order class (0.090 vs 0.021 g²ω\*); g²-law not reproducible in a two-mode toy (bath-limited λ ∝ γ, g-blind) — scope stated |

**Promotion note for the coordinator:** T2 §e's split verdict is
numerically confirmed on both halves: the error-signal half is SOUND and
now measured at field-anchored level (with the phase exponent's −½ shown
to be root-exact beyond the pure family); the entrainment half's
class-conditionality is real and now *pinned* — the restoring sign is not
in doubt (it is class-blind), but the *attractor* exists iff the bath is
dissipative, with rate → 0 linearly at the integrable point.

## Files

- `h27_entrain.py` — the part (b) toy + h27_results.json assembler
  (reruns in ~2.7 min; deterministic).
- `h27_results.json` — M0–M3 measurements, power law, boundary block,
  merged part (a) exponents.
- `h27a_field_exponents.json` — part (a) raw block (from `h26_convention`).
- Part (a) exact-family exponents: `h26_dilation.json`
  (`h27a_exponents_dilation`).
