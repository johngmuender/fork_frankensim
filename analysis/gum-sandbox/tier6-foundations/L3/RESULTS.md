# L3 — Lock-2 conservation lemma: the leakage-exponent ladder (F-T6-5-EXEC)

**Goal.** Toy-grade discharge of named debt **T4-W3**: GUM VIII.D Lock 2 prints the
stiff-sector leakage exponent (c/c_L)^5 at finite c_L, which presupposes that the
monopole and dipole channels of the elliptic-constraint source vanish — a step the
corpus never proves. Here we show by exact computation that **source conservation
supplies exactly that step**: the radiated-power exponent in 1/c_L climbs the ladder
1 → 3 → 5 as the source acquires charge conservation and then momentum conservation.

*Within-model; nothing here bears on nature.*

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G0 | flux consistency: P(R1) vs P(R2) agree ≤ 1%, every case, every c_L | max rel. diff = 2.50e-4 | **PASS** |
| G1 | monopole slope = 1 ± 0.1 | 1.000000 | **PASS** |
| G2 | dipole slope = 3 ± 0.15 | 2.999844 | **PASS** |
| G3 | quadrupole slope = 5 ± 0.25 | 4.999869 | **PASS** |
| G4 | conservation lemma stated and proof-sketched | see below | **PASS** |

## Key numbers

- ω = 1, array half-size a = 0.05, so k a = ωa/c_L ∈ [1.58e-3, 0.05] over the sweep.
- c_L sweep: 12 points, logspace over 1.5 decades, c_L ∈ [1, 31.62].
- Far-field radii: R1 = 2000, R2 = 4000. Spec floor 50·c_L,max/ω = 1581; R1 exceeds it,
  giving kR = ωR/c_L ≥ 63.2 at every c_L (deep wave zone for the largest c_L used).
- Fitted slopes of ln P vs ln(1/c_L): **1.000000 / 2.999844 / 4.999869**
  (max fit residuals 1.8e-15 / 1.9e-4 / 1.0e-4 in ln P).
- Max R1-vs-R2 flux discrepancy: monopole 0 (exact), dipole 1.9e-4, quadrupole 2.5e-4.

## Method

Exact retarded solution of (1/c_L²)∂²p/∂t² − ∇²p = s for arrays of point sources with
harmonic amplitude cos(ωt): p(x,t) = (1/4π) Σᵢ qᵢ cos(ω(t − rᵢ/c_L))/rᵢ, rᵢ = |x−xᵢ|.
No PDE grid. ∂p/∂t is a sum of sinusoids at frequency ω, so the time average
⟨(∂p/∂t)²⟩ = ½|A|² is computed **analytically** from the complex amplitude
A = (ω/4π) Σᵢ qᵢ e^{−iωrᵢ/c_L}/rᵢ (common phase e^{−ikR} factored out for accuracy).
The far-field intensity proxy I = ⟨(∂p/∂t)²⟩/c_L is integrated over the sphere at radius
R by 64-node Gauss–Legendre quadrature in cosθ (sources on the z-axis ⇒ φ symmetry).
R-independence of P (two radii, ≤1%) is the flux-consistency gate G0.

Arrays:
- **monopole** (non-conserved source): q = +1 at origin — total "charge" oscillates.
- **dipole** (charge-conserving, momentum-violating): q = ±1 at z = ±a — Σq = 0,
  dipole moment oscillates.
- **linear quadrupole** (conserved-source class): q = +1 at z = +a, q = −2 at 0,
  q = +1 at z = −a — Σq = 0 **and** Σq·z = 0.

Slope of ln P vs ln(1/c_L) fitted by least squares over the full sweep.

## G4 — The conservation lemma (stated and proof-sketched)

**Lemma.** Let the constraint source be s(x,t) = ∂ρ/∂t where ρ is a locally conserved
density (∂ρ/∂t + ∇·j = 0 with total ρ and total momentum ∫x ρ d³x conserved, i.e.
d/dt ∫ρ = 0 and d²/dt² ∫xρ = 0 up to internal stresses). Then the l = 0 and l = 1
radiation channels of the retarded solution cancel identically, and the leading
leakage power is quadrupolar: P ∝ (ωa/c_L)⁴ · (monopole normalization) i.e.
**P ∝ (c/c_L)⁵** in the corpus's normalization.

**Proof sketch (retarded multipole expansion).** In the far zone the retarded field is
p ≈ (1/4πR) Σᵢ qᵢ(t − rᵢ/c_L). Expanding the retardation about the array center,
rᵢ ≈ R − n̂·xᵢ, gives

  p ≈ (1/4πR) Σₗ (1/l!) (∂/∂t)ˡ (n̂·xᵢ/c_L)ˡ-weighted moments
    = (1/4πR) [ Q(t_r) + (1/c_L) d/dt (n̂·D)(t_r) + (1/2c_L²) d²/dt² (n̂n̂:Q₂)(t_r) + … ],

with Q = Σqᵢ (monopole), D = Σqᵢxᵢ (dipole), Q₂ = Σqᵢxᵢxᵢ (quadrupole). Each
successive order carries one extra factor ωa/c_L in the far-field **amplitude**, hence
(ωa/c_L)² per order in **power**. The radiated power ladder, including the 1/c_L flux
factor I = ⟨p_t²⟩/c_L, is therefore

  P ∝ (1/c_L) · (ωa/c_L)^{2l}  →  exponents 1, 3, 5 for l = 0, 1, 2.

Conservation kills the bottom rungs: if s = ∂ρ/∂t with ∫ρ conserved, then Q̇ = 0 and
the l = 0 amplitude ∝ Q̈ vanishes identically; if additionally total momentum is
conserved, D̈ = d/dt ∫j d³x = 0 and the l = 1 amplitude vanishes identically. The first
surviving channel is l = 2, giving P ∝ (1/c_L)⁵ ≡ (c/c_L)⁵. The exact computation above
realizes each rung: the monopole array has Q ≠ 0 (slope 1.000000); the dipole array has
Q = 0 but D ≠ 0 (slope 2.999844); the quadrupole array has Q = 0 and D = 0 (slope
4.999869). ∎

## Connection to the GUM corpus

GUM VIII.D Lock 2 prints leakage (c/c_L)^5 ≤ 1e-20; the T4 audit flagged the exponent
as presupposing the vanishing of the monopole ((c/c_L)^1) and dipole ((c/c_L)^3)
channels, with honest fallback (c/c_L)^3. This workstream isolates the missing step as
a **conservation lemma** and discharges it at toy grade: **IF** the constraint source is
built from conserved knot mass/momentum densities (the residual named assumption, now
isolated), the printed exponent 5 is restored. Filed as **F-T6-5-EXEC**.

Fallback arithmetic: with c_L/c ≥ 1e4 (the corpus's Bell-baseline floor), exponents
1 / 3 / 5 give leakage **1e-4 / 1e-12 / 1e-20**. Without the conservation lemma the
defensible bound is the dipole rung, 1e-12; with it, the printed 1e-20 stands.

## Caveats

- Toy grade: scalar wave equation with prescribed point-source arrays, not the actual
  GUM elliptic-constraint sector; the lemma transfers only insofar as the Lock-2 source
  is genuinely of the form s = ∂ρ/∂t with conserved ρ and momentum — that remains a
  **named assumption** of the corpus, isolated but not proven here.
- The measured slopes carry small O((ka)²) corrections (max fit residual ~2e-4 in ln P),
  visible as the tiny deviations of the dipole/quadrupole slopes from 3 and 5.
- The intensity proxy I = ⟨p_t²⟩/c_L equals the true flux only in the wave zone;
  residual O(1/(kR)) contamination is bounded by the G0 two-radius check (≤2.5e-4).
- Fixed ω = 1 throughout; the ladder is in 1/c_L at fixed source geometry, matching the
  corpus's stiff-limit setup.
