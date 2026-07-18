# Q2 — analytic reconnaissance of the terminal T4-W5 proof (F-T7-Q2)

**Goal.** A scoped analytic memo on the terminal open item: the fully analytic
proof of the L7b kill-bin zeros (t in [6.4, 6.8] and [6.8, 7.2], |n_y| >= 0.75),
whose required shape N3 fixed (transport/flow-map, not line-local) and O1/P2/Q1
demonstrated computationally. Deliverables per the pre-registered spec:
(i) one PROVEN lemma — the late-time ballistic asymptotic for this field class,
with an explicit, rigorous, computable error bound, numerically verified at
>= 100 points against the engine's field; (ii) the transport/no-slow-lane
mechanism written as a precise conjecture with its proof obligations
enumerated; (iii) an explicit statement of what the lemma alone already
yields, with the honest delta to the kill-bin theorem.

**Campaign:** Tier 7 Phase Q (ROADMAP_v9_PROGRAM.md addendum, Q2).
**Date:** 2026-07-18. **Code:** `q2_check.py` (runtime ~2 s).
**Raw output:** `q2_results.json`.
**Prior art (mandatory reads, honored):** tier7-program/N3/RESULTS.md (the
refuted line criterion + transport diagnosis), O1/RESULTS.md (the
backward-reachability theorem + measured crossing epochs),
tier6-foundations/L5prime/l5prime_run.py (initial data and engine),
P2/p2_certify.py (exact band-limited evaluator pattern, reused verbatim).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**ALL THREE GATES PASS, with the honest delta stated plainly:** the late-time
ballistic law is PROVEN for the continuum field of the compactly supported
L5' datum as an EXACT identity plus a moment-bounded error — a(x,t) =
(x − x0)/t + E with |E| <= Dbeta(w,t)/t = O(1/t^2) explicit, and the same
bound for the spin-drift generator b — verified at 122 (x,t) points against
the engine's exact 2048-mode field with ZERO violations (min headroom
1.8e-3, engine-vs-continuum delta median 1.6e-3); the lemma already
certifies, by pen-and-paper constants plus 1-D quadratures, that the line
criterion's hypothesis is false at every admissible t (continuum
a(d,t) > 0 on t in [4, 200], kill-window bounds a in [0.776, 0.999]) and
that the z-advection has the mechanism's sign (b in [0.104, 0.228] > 0 in
the kill window at the line); the no-slow-lane transport conjecture is
stated with five proof obligations O-1..O-5; but the lemma's own rigorous
no-crossing reach is T_ball(n_y = 1) = 67.3 (T_ball(0.75) = 139.7) under
confinement conditions — a factor ~9–19 in time ABOVE the kill bins
[6.4, 7.2]. **The lemma controls the far ballistic tail, not the kill
window.** The kill-bin theorem still requires the early-time flow-map
control (obligation O-1), which no asymptotic argument can supply.

## 1. Setup, and which part of v the lemma controls

The bench state is exactly factorized (L5prime, engine-inherited):

    phi(x, z, t) = f(x, t) sin(pi z) e^{-i pi^2 t/2},

and the guidance velocity of member n_y decomposes exactly (L7b/N3):

    v_x(x, z, t) = a(x, t) + n_y S(z),   S(z) = -pi cot(pi z),
    v_z(x, z, t) = n_y b(x, t),
    a = Im(f'/f),   b = Re(f'/f).

S(z) is a known, exact, t-independent function — it is NOT part of any
error term. **The lemma proven below controls a and b — i.e., the entire
axial part of v_x and the entire z-advection v_z.** The z-dependence of v
is exactly known (ground mode); nothing about z is asymptotic.

The x-datum is f0(y) = A0(y − x0) e^{i k0 y} with x0 = −5, k0 = 2, and

    A0(u) = e^{-u^2/2} W(|u|),  W = 1 on [0, 2.5], S2-smoothstep down to 0
    at 3.5, supp A0 = [-3.5, 3.5],  A0 real, even, nonnegative, C^2.

**Spectrum check (spec item).** The datum is NOT strictly band-limited and
its spectrum is NOT strictly positive: Ahat0 (the transform of the
envelope) is entire; the boosted spectrum is centered at k0 = 2 with
Gaussian-class width 1, so it has mass at k < 0 (relative envelope
amplitude e^{-2} ≈ 0.13 at k = 0), and the C2 truncation adds |k|^-4-class
tails on all of k-space. Honest handling chosen (spec's option (a),
adapted): the lemma is proven EXACTLY for the **continuum evolution of the
compactly supported datum** — compact support in x is what makes every
constant finite and computable, and no band-limited proxy is needed on the
proof side. The gap to the engine's field (band-limited periodic
interpolant, the campaign's certified object since O1/Q1) is then bounded
NUMERICALLY at every verification point (Section 3), and the exact
finite-mode-sum decomposition is given too (Lemma A'), with the honest
reason no t-uniform ballistic bound can exist for the periodic field.

## 2. The proven lemma

### Lemma A (ballistic law with explicit error, continuum field)

Let f(·, t) be the free evolution (hbar = m = 1) of f0 above. Define

    xi = x − x0,   vtil = xi/t   (ray velocity from the source center),
    w = vtil − k0  (lag/lead label),
    Ahat0(w) = INT A0(u) e^{-i w u} du     (real, even; entire),
    beta(w)  = Ahat0'(w) / Ahat0(w),
    m2 = INT u^2 A0(u) du = 2.424780,   m3 = INT |u|^3 A0(u) du = 3.734202.

**Claim.** For every (x, t) with t > 0 satisfying the admissibility
condition Ahat0(w) > m2/(2t):

    (B1)  | a(x,t) − (x − x0)/t | <= C_lem(x,t),
    (B2)  | b(x,t) − beta(w)/t  | <= C_lem(x,t),
    (B3)  | a(x,t) − x/t |        <= ( |x0| + Dbeta(w,t) ) / t,

    C_lem(x,t) = Dbeta(w,t)/t,
    Dbeta(w,t) = [ m3/(2t) + |beta(w)| m2/(2t) ] / [ Ahat0(w) − m2/(2t) ].

At fixed w, Dbeta = O(1/t), so C_lem = O(1/t^2): the ballistic law
a → (x − x0)/t holds with a rate, and the spin-drift generator satisfies
b → beta(w)/t with the same rate. All constants are explicit 1-D integrals
of the known datum.

### Proof (exact Fresnel decomposition; no asymptotics folklore)

**Step 1 (exact factorization).** The free propagator gives, substituting
y = x0 + u and expanding (xi − u)^2/2t = xi^2/2t − xi u/t + u^2/2t:

    f(x,t) = (2 pi i t)^{-1/2} INT e^{i(x−y)^2/2t} f0(y) dy
           = (2 pi i t)^{-1/2} e^{i(k0 x0 + xi^2/2t)} Phi_t(w),

    Phi_t(w) := INT A0(u) e^{-i w u} e^{i u^2/2t} du .

This is an identity for every t > 0 — the factor e^{i u^2/2t} is the exact
Fresnel remainder of the stationary-phase decomposition, not an expansion.
(Numerically confirmed: direct Fresnel quadrature of f, f' vs the
decomposition agrees to 4.4e-16 at 8 points; `identity_check`.)

**Step 2 (exact logarithmic derivative).** x enters only through xi^2/2t
and w (with dw/dx = 1/t), so wherever Phi_t(w) ≠ 0:

    f'/f = i xi/t + Phi_t'(w) / ( t Phi_t(w) ),

    =>  a = xi/t + Im[Phi_t'/Phi_t]/t ,      b = Re[Phi_t'/Phi_t]/t .

**Step 3 (moment bounds).** On supp A0, |e^{i u^2/2t} − 1| =
2|sin(u^2/4t)| <= u^2/2t, hence

    |Phi_t(w)  − Ahat0(w) | <= m2/(2t),
    |Phi_t'(w) − Ahat0'(w)| <= m3/(2t)      (Phi_t' has the extra −iu).

**Step 4 (ratio transfer).** Since A0 is real and even, Ahat0 is real
(imag residue 4.4e-16 measured) and Ahat0' real; so beta(w) is real. With
Ahat0(w) > m2/(2t) we get |Phi_t| >= Ahat0 − m2/(2t) > 0 and

    Phi_t'/Phi_t − beta = ( (Phi_t'−Ahat0') Ahat0 − Ahat0' (Phi_t−Ahat0) )
                          / ( Phi_t Ahat0 ),

    |Phi_t'/Phi_t − beta| <= [ (m3/2t) Ahat0 + |Ahat0'| (m2/2t) ]
                             / [ Ahat0 (Ahat0 − m2/2t) ]  =  Dbeta(w,t).

beta real then gives BOTH |Im(Phi'/Phi)| <= Dbeta and
|Re(Phi'/Phi) − beta| <= Dbeta; dividing by t yields (B1) and (B2).
(B3) follows from x/t = xi/t + x0/t. ∎

**Corollary A.1 (certified sign of the z-advection).** b >=
(beta(w) − Dbeta(w,t))/t. In the lag cone w < 0 one has beta(w) > 0
(Ahat0 even, peaked, decreasing in |w| on the admissible window), so the
sign of the spin drift is CERTIFIED wherever Dbeta < beta — e.g. at every
kill-window detector point (Section 4 numbers).

**Corollary A.2 (ray rigidity).** For the axial member (v_x = a), along
any trajectory confined to an admissible cone on [T0, t], the ray label
u(s) = (X(s) − x0)/s obeys |u(t) − u(T0)| <= sup_cone Dbeta ·
(1/T0 − 1/t): late-time axial trajectories are asymptotically straight
rays through the source center, with an explicit rate.

### Lemma A' (exact decomposition for the engine's finite mode sum, and
### why no t-uniform ballistic bound exists for it)

The engine's field is the finite sum f_N(x,t) = SUM_k Ctil_k e^{i kappa_k
x − i kappa_k^2 t/2} (2048 modes, spacing pi/30). The same completion of
the square (i kappa x − i kappa^2 t/2 = i x^2/2t − i t (kappa − x/t)^2 /2)
gives the EXACT discrete-Fresnel decomposition

    f_N = e^{i x^2/2t} H_t(v),  H_t(v) = SUM_k Ctil_k e^{-i t (kappa_k − v)^2/2},
    f_N'/f_N = i x/t + H_t'(v) / ( t H_t(v) ),   v = x/t,

so a_N = x/t + Im[H_t'/H_t]/t exactly — the "Fresnel decomposition of the
finite mode sum" the spec asked for. But H_t is a discrete Gauss/Fresnel
sum: by Poisson summation f_N(x,t) = SUM_j f_BL(x + 60 j, t) is the
60-periodization of the band-limited interpolant's free evolution, and as
t grows the spreading images ENTER the window — the ballistic law for the
periodic field must fail at large t, provably in the measurement:
a_N(1, 15.9) = −0.460 while the continuum value is +0.385 with certified
lower bound +0.352 (late-layer table below). Hence Lemma A is stated for
the continuum field, and its control of the ENGINE field is exactly
co-extensive with the wrap-controlled window — where it is verified next.
(This also resolves an N3 footnote: N3's late-time A(t) sign flip at
t ≈ 15.9, and hence the "axial support endpoint 15.87", are properties of
the periodized model field; the continuum line velocity is provably
positive there. The campaign's certified object is the discretized field
— O1/Q1 doctrine — so no prior result moves; the lemma quantifies the
continuum-vs-model layer that doctrine names.)

## 3. Numerical verification of the lemma (gate G1)

Constants (piecewise Gauss–Legendre, 200 nodes/piece on the four analytic
pieces of A0; node-doubling convergence <= 4.2e-14 on everything; mpmath
dps-30 cross-check of Phi_t max diff 3.0e-14):

    m0 = 2.498621, m2 = 2.424780, m3 = 3.734202;
    Ahat0 first zeros at |w| = 3.715, 4.551 (verification band |w| <= 1.8
    is far inside the first zero; min Ahat0 on band = 0.490).

**122 verification points** (121 distinct; 9 t-values in [4.8, 10] x 13
w-values in [−1.7, 1.5], plus 5 exact detector points (x = 1, t in
[6.4, 7.2])), each checked for all of (B1), (B2), (B3) against BOTH:

- the **continuum field** (independent high-accuracy quadrature of
  Phi_t, Phi_t'), and
- the **engine's field** (P2's exact evaluator: direct 2048-mode sum of
  the Nx = 2048 FFT spectrum on [−15, 45]; the object O1/P2/Q1 certify).

Result: **ZERO violations in all 5 x 122 = 610 checks.** Headroom
(bound minus actual error): engine-a min 1.8e-3 / median 2.6e-2; engine-b
min 5.2e-3; continuum-a min 7.2e-3. Engine-vs-continuum delta (the
band-limitation + periodization layer): a median 1.6e-3, max 3.0e-2 (the
max at (x, t) = (26.5, 9), a low-density lead-cone point, rho_rel = 0.11);
at the detector points the delta is 4e-4..1.9e-3. All points are
well-conditioned (rho_rel in [0.057, 0.99]; engine floor/clamp never
bind). The tightest headroom sits at t = 9–10 near w ≈ 0 where C_lem is
smallest (~0.009) and the wrap onset contributes ~6e-3 — consistent with
the late-layer table: the engine delta grows to 3.6e-2 at t = 12 and
0.84 at t = 15.9. The lemma's verified engine-window is t <= 10 here;
the kill window t <= 7.2 sits deep inside it.

Late-layer diagnostic (x = d = 1):

| t | a_continuum | certified lower bound | a_engine | delta |
|---|---|---|---|---|
| 8.0 | +0.7700 | 0.6970 | +0.7572 | 1.3e-2 |
| 10.0 | +0.6150 | 0.5563 | +0.6031 | 1.2e-2 |
| 12.0 | +0.5114 | 0.4637 | +0.4755 | 3.6e-2 |
| 14.0 | +0.4375 | 0.3982 | +0.0257 | 4.1e-1 |
| 15.9 | +0.3846 | 0.3515 | −0.4598 | 8.4e-1 |

## 4. What the lemma certifies in the kill window (continuum layer)

Rigorous two-sided bounds at the detector line x = d = 1 (all from
(B1)/(B2) with the computed constants; N3's measured engine line field
A(t) lies inside the a-bounds at every sampled t, all 9 rows):

| t | a(1,t) in | b(1,t) in | A(t) measured (N3) |
|---|---|---|---|
| 6.4 | [0.876, 0.999] | [0.1039, 0.2278] | 0.9614 |
| 6.8 | [0.823, 0.942] | [0.1051, 0.2242] | 0.9043 |
| 7.2 | [0.776, 0.891] | [0.1054, 0.2200] | 0.8574 |

Consequences, both now ANALYTIC (pen-and-paper constants + 1-D
quadratures, no evolution engine involved):

1. **N3's refutation is confirmed analytically.** a(d, t) > 0 throughout
   the kill window — indeed the scan of the certified lower bound
   vtil − C_lem at x = d is positive for ALL admissible t in [4.0, 200]
   (min 0.0297, at t = 200). The line criterion's hypothesis
   sup_z v_x(d, z, t) < 0 is provably false at every late time for the
   continuum field: v_x(d, 1/2, t) = a(d, t) > 0. No line-local
   certificate exists — now a theorem for the continuum layer, not just
   a measurement.
2. **The transport mechanism's sign structure is certified.** b(d, t) >=
   0.104 > 0 in the kill window: for n_y > 0 the z-advection v_z = n_y b
   pushes UP toward the wall z = 1 where S(z) → +infinity (fast lane);
   for n_y < 0, down toward z = 0 where n_y S(z) → +infinity likewise.
   The Das–Dürr lifting that empties the not-yet-crossed set has the
   right sign, rigorously, exactly where the kill bins live.

## 5. The conjecture (gate G2): no-slow-lane transport, precisely

**Near-wall spin drift, derived (spec item).** Because the z-factor is
exactly the ground mode, rho = 2|f|^2 sin^2(pi z) and

    v_z = n_y Re(f'/f) = n_y (d rho/dx) / (2 rho)   EXACTLY,

independent of z — the "b_z ~ +n_y (d rho/dx)/(2 rho)-class" near-wall
asymptotic the spec anticipated is exact and z-UNIFORM in this bench;
there is no wall correction to derive. The wall structure enters v only
through the known S(z) in v_x. Late-time form: v_z = n_y beta(w)/t +
O(1/t^2) with beta = (ln Ahat0)'(w) > 0 in the lag cone (proved,
Corollary A.1): lagging trajectories (w < 0) are advected toward the
n_y-favored wall at rate ~ |n_y| |w|/t.

**Conjecture C (finite emptying time / no slow lane).** For the bench
field model and every |n_y| >= n* = 0.75 there exists T_max(n_y) <
infinity such that every trajectory starting in supp rho_0 (up to
explicitly flux-bounded near-wall strips of the O1 class) either first
crosses x = d by T_max(n_y) or never crosses; the first-crossing support
is compact. Empirically T_max(1) ≈ 5.14, T_max(0.75) ≈ 6.02.

**Structure of the intended proof and its obligations** (a competent
analyst could pick these up; each is stated with what would suffice):

- **O-1 (early-time flow control on [0, T0]).** The lemma is vacuous
  below t_adm(w) = m2/(2 Ahat0(w)) (≈ 1.0–2.5 on the lag band) and its
  constants are useful only for t >~ 4–5, while the measured cutoffs
  (5.1/6.0) sit mostly INSIDE [0, T0]. Sufficient: a validated-enclosure
  covering of the flow map on [0, T0] x supp rho_0 (Q1's interval
  technology extended from single backward paths to a Lipschitz covering
  — N3's diagnosis item 4), certifying where the not-yet-crossed set IS
  at T0. No asymptotic argument can replace this piece. **This is the
  dominant obligation.**
- **O-2 (lag-cone localization for t >= T0).** Show the not-yet-crossed
  set stays in an admissible lag cone {w(x,t) in [w_-, w_+], x < d} with
  w_- above the first Ahat0 zero (−3.715). Sufficient: a barrier argument
  along w-level sets using (B1) (dw/dt = (a − xi/t + n_y S(z))/t, with
  |a − xi/t| <= Dbeta/t), plus a rho-weight bound for the far-lag tail
  (equivariance + transported-measure smallness, the O1/N3 flux-bounding
  pattern made rigorous).
- **O-3 (quantitative monotone z-lift).** On the lag cone, b >=
  (beta(w) − Dbeta(w, t))/t > 0 (proved for the continuum; needs
  engine-layer transfer, see O-5). Monotone functional along
  not-yet-crossed trajectories: **V(t) = sigma Z(t), sigma = sign(n_y)**
  — nondecreasing while the trajectory lags (dV/dt = |n_y| b >= 0,
  z-uniformly), giving the explicit clock

      sigma Z(t) >= sigma Z(T0) + |n_y| beta_eff ln(t/T0),
      beta_eff = min over the lag band of (beta − Dbeta)|_{t=T0},

  so z reaches any level z-bar < 1 by T0 exp((z-bar − z0)/(|n_y|
  beta_eff)). Compounding invariant-region structure for the (z, x − d)
  system: once past mid-channel on the favored side, dw/dt >=
  (|n_y| |S(z)| − Dbeta/t)/t > 0 — the wall push RAISES the ray label,
  so z-lift and x-catch-up reinforce; candidate Lyapunov functional for
  the pair: Lambda = w + c1 sigma z (obligation: exhibit c1 > 0 and
  phi > 0 with dLambda/dt >= phi(Lambda)/t on the not-yet-crossed set).
- **O-4 (wall-lane crossing).** Near the favored wall (1 − sigma-side
  distance eps), S contributes |n_y| pi cot(pi eps) ~ |n_y|/eps to v_x
  while |a| is bounded on the cone by |vtil| + Dbeta/t; sufficient: an
  explicit eps-bar and left barrier x_inf (from the (B1) lower envelope
  d xi/dt >= xi/t − Dbeta/t) such that every trajectory entering the
  lane crosses within Delta-t <= (d − x_inf) eps-bar / |n_y|, and cannot
  leave the lane while lagging (v_z sign, O-3).
- **O-5 (assembly, uniformity, engine transfer).** T_max(n_y) = T0 +
  T0 exp((z-bar − z_min)/(|n_y| beta_eff)) + Delta-t(n_y); monotone
  decreasing in |n_y|, divergent as n_y → 0 — consistent with the axial
  family's unbounded first-crossing support and with the measured
  ordering tau_max(1) < tau_max(0.75). All continuum inequalities must
  finally be transferred to the certified band-limited periodic field on
  the wrap-free window: either a rigorous aliasing/image bound in the
  Lemma-A' decomposition, or Q1-grade interval certification of the
  finitely many field inequalities used (the deltas measured here,
  median 1.6e-3, say the transfer is cheap in principle).

**Heuristic consistency (flagged, NOT evidence, no tuning):** the O-3
clock with raw beta(−1.1) = 1.10, T0 = 2, Delta-z = 0.9 gives T_max
ballparks 4.53 (n_y = 1) and 5.95 (0.75) vs measured 5.13 / 6.02 — the
mechanism's magnitude and its |n_y| ordering come out right, which is
what a mechanism-level conjecture should reproduce before anyone invests
in O-1.

## 6. What the lemma alone already yields (gate G3) — and the honest delta

**Y-1 (ray rigidity / no-early-crossing, axial).** Corollary A.2 pins any
cone-confined axial trajectory to its T0 ray label up to sup Dbeta / T0;
in particular a first crossing of d at time t forces |6/t − u(T0)| <=
that rate — crossing times are determined by the ballistic label, with an
explicit error.

**Y-2 (line-criterion refutation, analytic).** Section 4: continuum
a(d, t) > 0 for all admissible t in [4, 200] — the field-sign route to
the kill bins is closed by proof, not just by measurement (N3's
diagnosis, elevated one layer).

**Y-3 (the no-slow-lane far-tail bound, T_ball).** Rigorous statement,
from Corollary A.1 + the O-3 clock, continuum layer: fix the lag band
w in [−1.2, −0.6] and the mid-channel band z in [0.05, 0.95]. For
t >= T0 = 7.5, beta_eff = min(beta − Dbeta)|_{T0} = 0.410 on this band.
**No trajectory can remain in (lag band) x (z band) beyond**

    T_ball(n_y) = T0 exp( 0.9 / (|n_y| beta_eff) ):
    T_ball(1) = 67.3,   T_ball(0.75) = 139.7

(optimized over T0 in [3, 12]; the z-band exit is through the favored
wall side, by monotonicity of sigma z). Under conditions C =
{cone-and-band confinement since T0; engine-layer transfer per O-5},
there are consequently no first crossings at the detector after
T_ball(n_y) from such trajectories: they have left the slow lane.

**The honest delta.** The kill bins are [6.4, 7.2]. T_ball is 67.3 resp.
139.7 — a factor **9.3–19.4 in time beyond the kill window** — and even
its hypotheses (confinement since T0 = 7.5) begin after the kill window
ends. **The kill bins fall OUTSIDE the lemma's reach: the lemma controls
the far ballistic tail, not the kill window.** What the kill-bin theorem
needs and the lemma cannot give: (a) O-1, control of the flow on
[0, ~5–6] where the actual emptying happens — a finite-time,
non-asymptotic statement over ALL of supp rho_0, for which the only
identified technology is the Q1 validated-enclosure covering; (b) O-2's
localization; (c) O-5's engine transfer. The lemma supplies the
mechanism's sign, rate, and late-time closure — the *shape* of the
terminal proof — but not its kill-window core. Anyone claiming the
kill-bin zeros from ballistic asymptotics alone would be wrong by an
order of magnitude in time.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | lemma proven (checkable derivation, explicit error), numerically verified at >= 100 points, error within the derived bound, ZERO violations | Lemma A proven exactly (Section 2; identity confirmed to 4.4e-16 against direct Fresnel quadrature; mpmath dps-30 cross-check 3.0e-14; all constants explicit: m2 = 2.424780, m3 = 3.734202, Ahat0 computable, quadrature converged to 4.2e-14). Verified at 122 points (121 distinct, incl. 5 kill-window detector points) x 5 checks (B1/B2 continuum, B1/B2 engine, B3 engine) = 610 checks, **0 violations**; min headroom 1.8e-3 (engine a), 7.2e-3 (continuum a); engine-vs-continuum delta median 1.6e-3 / max 3.0e-2, reported per point | **PASS** |
| G2 | conjecture + proof obligations printed precisely (pick-up-able) | Conjecture C stated with quantifiers and the exception class; obligations O-1..O-5 each with a sufficient condition; near-wall spin drift derived (exact, z-uniform: v_z = n_y d(rho)/dx/(2 rho)); monotone functional V = sigma z and Lyapunov candidate Lambda = w + c1 sigma z named; mechanism sign certified in the kill window (b >= 0.104); heuristic magnitude check 4.53/5.95 vs measured 5.13/6.02, flagged non-evidence | **PASS** |
| G3 | the honest delta between lemma and needed theorem stated | T_ball(1) = 67.3, T_ball(0.75) = 139.7 (conditional, continuum layer) vs kill bins [6.4, 7.2]: outside the lemma's reach by 9.3–19.4x in time, with the dominant missing piece named (O-1 early-time flow control, [0, ~5–6], non-asymptotic) and the engine-layer transfer quantified (deltas measured; rigorous transfer = O-5) | **PASS** |

## Key numbers

- Lemma constants: m0 = 2.498621, m2 = 2.424780, m3 = 3.734202 (GL
  convergence <= 4.2e-14); Ahat0(0)/sqrt(2 pi) = 0.9968 (taper
  correction 0.32%); first Ahat0 zeros |w| = 3.715, 4.551; min Ahat0 on
  the verification band [−1.8, 1.6] = 0.490.
- Verification: 122 points, t in [4.8, 10], w in [−1.7, 1.5]; 0/610
  violations; headroom (engine a) min 1.8e-3 at (22, 10), median 2.6e-2;
  rho_rel in [0.057, 0.99]; engine-vs-continuum delta a: median 1.6e-3,
  max 3.0e-2; detector points delta 4.5e-4..1.9e-3.
- Kill window (continuum, rigorous): a(1, t) in [0.776, 0.999],
  b(1, t) in [0.104, 0.228] across t in [6.4, 7.2]; N3's measured A(t)
  inside the a-bounds at all 9 sampled t. Line-positivity scan: certified
  a(1, t) lower bound > 0 for all admissible t in [4.0, 200] (min
  0.0297).
- Late layer (x = 1): engine-vs-continuum delta 1.3e-2 (t = 8) →
  3.6e-2 (t = 12) → 0.84 (t = 15.9), where a_engine = −0.460 vs
  certified continuum lower bound +0.352 — the periodization wall.
- T_ball: beta on lag band [−1.2, −0.6] in [0.588, 1.207]; best
  T0 = 7.5 gives beta_eff = 0.410, T_ball(1) = 67.3, T_ball(0.75) =
  139.7; kill bins inside reach: **false**.
- Heuristic clock (flagged): T_max est 4.53 / 5.95 vs measured
  5.1310 / 6.0174.

## Method

All constants and transforms are 1-D piecewise Gauss–Legendre quadratures
(200 nodes per analytic piece of A0; pieces split at the taper junctions
±2.5 and at 0 for the |u|^3 kink), with node-doubling convergence checks
and an mpmath dps-30 cross-check of Phi_t. The decomposition identity is
verified against direct Fresnel quadrature of f and f'. The engine field
is evaluated by P2's exact pattern: direct summation of the 2048-mode
FFT spectrum of the grid datum on [−15, 45] (the object O1/P2/Q1
certify), no clamps binding at any verification point (rho_rel monitored,
min 0.057 vs floor 1e-14). N3's measured line field A(t) is read from
n3_raw.npz unmodified. T_ball scans T0 in [3, 12] x the lag band with
Dbeta evaluated at t = T0 (Dbeta decreasing in t makes beta_eff valid for
all t >= T0). hbar = m = 1. Runtime ~2 s; no long runs (nohup not
required).

## Honest scope

The proven object is a statement about the CONTINUUM free evolution of
the exact compactly supported datum, with every constant computable; its
application to the campaign's certified object (the band-limited periodic
engine field) is verified pointwise on the wrap-controlled window with
measured deltas, not proven — making that transfer rigorous is obligation
O-5, and the Lemma-A' decomposition shows why no t-uniform version can
exist for the periodic field. The conjecture section is a proof PLAN with
named sufficient conditions, not a proof; the heuristic T_max numbers are
flagged and play no evidentiary role. The kill-bin zeros remain, after
this workstream, exactly what O1/P2/Q1 left them: a field/flow statement
at proof-grade-modulo-discretization (Q1: partially interval-certified),
with the analytic T4-W5 proof open and its remaining core now isolated to
obligation O-1 (early-time flow-map control) plus O-2/O-5. Within-model;
nothing here bears on nature.

## Caveats

- Quadrature values (m2, m3, Ahat0, Phi_t) carry ~1e-13 numerical error;
  every gate margin is >= 10 orders above that. They are converged
  numerics, not interval enclosures; a fully formal version of Lemma A's
  CONSTANTS would enclose four 1-D integrals of an explicit elementary
  integrand — trivial for the Q1 machinery, not done here.
- The 122-point grid contains one duplicated point ((1.0, 7.2) appears
  as both a grid and a detector point); 121 distinct points, still >= 100.
- The verification window t <= 10 is where the engine field tracks the
  continuum (deltas <= 3e-2); the lemma is NOT verified (and provably
  fails, Lemma A') for the engine field at t >~ 12–14 — this is the
  periodization layer, quantified in the late-layer table, not a defect
  of the lemma.
- beta_eff and T_ball use the continuum bounds on a fixed band; different
  band choices move T_ball by factors of order one — none brings it
  within an order of magnitude of the kill bins (the T0 scan shows the
  minimum over T0 explicitly).
- The min engine-a headroom (1.8e-3 at t = 10) is small because C_lem is
  smallest exactly where the wrap onset begins to bite; at the kill
  window the headroom is >= 2.2e-2 with deltas <= 1.9e-3.
- Conjecture C's exception class (never-crossing zero-weight strips) is
  inherited from N3/O1's covering diagnostics; its precise measure-zero
  formulation is part of obligation O-2.
