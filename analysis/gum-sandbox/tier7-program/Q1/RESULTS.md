# Q1 — validated interval enclosures for the pre-crossing certificates (F-T7-Q1)

**Goal.** Execute the formal rung P2's G5 named: turn pre-crossing
certificates of kill-bin backward paths into MACHINE-CHECKED INEQUALITIES
by rigorous interval enclosure of the backward guidance flow over the exact
band-limited field — scoped (per the Phase-Q addendum) to the SHORT backward
segment from the detector point (d, z0, t0) until the enclosure tube's X_x
lower bound exceeds d + delta (delta = 1e-3), which is all the O1 theorem
consumes.

**Campaign:** Tier 7 Phase Q (ROADMAP_v9_PROGRAM.md addendum, Q1).
**Date:** 2026-07-18. **Code:** `q1_enclose.py` (phases test / g1 / certify /
gates / fig; per-path checkpoint files; production launched under nohup).
**Raw output:** `q1_results.json` + per-path `_q1_path_<idx>.json`,
`_q1_g1.json`, logs `q1_g1.log`, `q1_certify.log`. **Figure:** `q1_fig.png`.
**Prior art (mandatory reads, honored):** tier7-program/P2/RESULTS.md +
p2_certify.py (exact spectral representation, path data, worst-margin
record), tier7-program/O1/RESULTS.md (the theorem being certified);
P2/_selection.npz + _dop.npz + _ladder.json + p2_results.json inspected for
path selection.

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**ALL FIVE GATES PASS — including the stretch case:** the interval field
encloses all 16,000 high-precision reference evaluations at 10,000 random
points with 0 violations (G1); the P2 WORST-MARGIN path (selection index
232, n_y = +0.75, t0 = 6.405, z0 = 0.6825, P2 margin 0.6206790034) is
CERTIFIED by a validated enclosure 3.555 time units deep — 233,002 rigorous
steps at h = 2^-16 ending at s_stop = 2.8496716 with certified bound
X_x >= 1.0010189 > d + delta and tube width 7.84e-4 <= 1e-3 at the
certificate point (G2); six additional paths spanning both kill bins and
both |n_y| families, including the small-margin path 395 (P2 margin 1.305,
depth 2.73, final width 4.2e-5), are certified (G3); there are ZERO
failures to account — no tube ever wrapped, no Picard step ever needed
inflation retries, no width cap was hit (G4); and the trust statement is
printed (G5). Combined with O1's backward-reachability theorem, the
statement "these seven kill-bin detector points carry no first-crossing
arrival" is now a machine-checked inequality GIVEN the trust base below.

## The certificate semantics

O1's theorem: if the unique backward trajectory through detector point
(d, z0) at t0 satisfies max_{s<t0} X_x(s) >= d, the point carries no
first-crossing (arrival) density. Q1 certifies, for each path, the
rigorous statement

    X_x(s_stop)  >=  x_lo  >  d + delta,   s_stop < t0,  delta = 1e-3,

where [x_lo, x_hi] is a validated enclosure of the true backward
trajectory at s_stop (existence + uniqueness from Picard-Lindelof; the
tube keeps |f|^2 >= 0.022 and z in [0.016, 0.69] — far inside the
smooth region). The engine's regularization provably NEVER binds on any
tube (rho >= 3.9e12 x the engine floor, velocity components >= 237x below
the +-500 clamp, z inside [zclip, 1-zclip]), so the enclosed flow is
simultaneously that of the ideal velocity field and of the engine's
regularized law — preserving O1/P2 semantics exactly.

## Gates

| id | spec (Phase-Q addendum) | measured | verdict |
|----|------|----------|---------|
| G1 | interval field encloses >= 1e4 random exact-point evaluations, 0 violations; outward-rounding sanity >= 100 cases | 16,000 containment checks at 10,000 random (x,z,t,n_y) points — f and f' at all 10,000, velocity components a, b at 2,000 — against gmpy2/MPFR 200-bit full 2048-mode reference values (reference error <= 1e-55, far below enclosure widths): **0 violations**. Outward-rounding sanity: 120/120 cases (60 fat-box interior-disc containments, 40 thin-point width-positivity + reference containments, 20 scalar-layer cross-checks vs mpmath.iv). Max thin-field enclosure width 1.07e-11 | **PASS** |
| G2 | P2 worst-margin path certified: rigorous lower bound X_x > d + delta, widths <= 1e-3 at the certificate point | Path idx 232 (O1 flat 32136; n_y = +0.75, t0 = 6.405, z0 = 0.6825; P2 margin_mp 0.620679003382...): **CERTIFIED**, s_stop = 2.8496716309, certified X_x >= 1.0010188761 (delta achieved 1.019e-3 > 1e-3), tube width at certificate **7.845e-4 <= 1e-3**, certification depth 3.5553 (the full pre-crossing excursion: the tube dips to x ~ -1.26 and returns), 233,002 steps at h = 2^-16, 1,225 s, 0 Picard retries | **PASS** |
| G3 | >= 5 additional paths certified, spread over both kill bins and both \|n_y\| families, incl. >= 1 small-margin path | **6/6 additional paths CERTIFIED**: 134 (+1.00, bin1, depth 0.083, w 2.8e-6), 194 (+1.00, bin2, 0.135, 5.7e-6), 372 (+0.75, bin1, 0.286, 1.3e-5), 375 (+0.75, bin2, 1.017, 1.4e-4), 212 (+1.00, bin1, wall-adjacent z0 = 0.0175 immediate-exit, w 4.9e-10), 395 (+0.75, bin1, **small-margin** P2 margin 1.305, depth 2.725, 178,578 steps, w 4.2e-5). Families {0.75, 1.0} both covered, bins [6.4,6.8] and [6.8,7.2] both covered | **PASS** |
| G4 | failure accounting: every wrap/blow-up reported with widths | **Zero failures**: 7/7 attempted paths certified; Picard containment succeeded on first try at every one of ~436k total steps (0 retries); max tube width ever reached = the 7.84e-4 at 232's certificate point; width cap 2e-2 never approached | **PASS** |
| G5 | trust statement + scope | Printed below and in q1_results.json | **PASS** |

## Key numbers

- **Worst-margin path (G2).** s_stop = 2.8496716309 sits 1.93e-4 below
  P2's refined crossing epoch s_pc = 2.8498644664 — the enclosure stops at
  the first step whose LOWER bound clears d + delta, i.e. just past the
  crossing; the two epochs agree as they must. Certified clearance
  x_lo - d = 1.019e-3; tube width at certificate 7.84e-4 (gate headroom
  1.27x — the thinnest margin in this workstream, reported as such); the
  width history is monotone (max over tube = width at stop) and the final
  growth spurt (5e-5 -> 7.8e-4 over s in [3.5, 2.85]) is the genuine
  variational amplification of the small-margin approach, matching the
  float-recon transition-matrix integral (~46) within a factor ~2.
- **Small-margin path (G3).** idx 395 is the shallowest genuinely
  small-margin path in the P2 subset (margin 1.305; all margins < 1.5 have
  depth >= 2.7): certified 2.725 deep at h = 2^-16, final width 4.2e-5 —
  24x under the gate.
- **Depth vs width scaling.** Certified widths grow with depth exactly as
  the remainder model predicts: 2.8e-6 (depth 0.08, h = 2^-14) -> 1.4e-4
  (1.02, 2^-14) -> 4.2e-5 (2.73, 2^-16) -> 7.8e-4 (3.56, 2^-16).
- **Field-enclosure tightness.** Fat-box disc radii match the analytic
  dependency-width model sum_k |C_k| kappa^n (kappa w_x + kappa^2 h/4)
  to within 5% (measured during design; the tight sqrt-based magnitude
  bound matters — a |Re|+|Im| bound compounds sqrt(2) per squaring in the
  21-level phase-power chain and inflates top-mode radii ~500x, which was
  caught and fixed before production).
- **Engine-law coincidence (rigorous).** Over all 7 tubes: min |f|^2 lower
  bound 0.0222 vs engine floor upper bound 5.7e-15 (headroom >= 3.9e12);
  max |f'/f| upper bound 2.11 vs clamp 500 (headroom >= 237x); z ranges
  within [0.0175, 0.6825] vs clip 1e-9. No clip, floor, or z-clip can bind
  anywhere on any tube.
- **Runtime.** g1: 517 s. certify: 1,225 s wall for all 7 paths (4
  workers; the worst path is the critical path at 233k steps, ~5.3 ms per
  validated step incl. both field evaluations). Per-path checkpoints every
  20k steps (never needed — no restart occurred).

## Method (the validated construction, as printed by the script)

**Field.** The certified object is P2's exact band-limited field: the full
2048-mode sum f(x,t) = sum_k C_k e^{i kappa_k (x - xmin)} e^{-i kappa_k^2
t/2} (float64 C_k declared exact model data; NO truncation — the C2-taper
spectrum's tails are too heavy for a useful truncation bound, checked and
rejected: dropping even 1,848 of 2,048 modes leaves sum|C| tail 1.3e-2).
Complex quantities are midpoint-radius (disc) intervals: center complex
float64 + radius upper bound; Gargantini-Henrici product/quotient
formulas; radius arithmetic padded outward by relative 2^-46 (= 128 ulp)
+ 5e-324; center magnitudes bounded by sqrt(re^2+im^2)(1+2^-40) (IEEE-754
sqrt is correctly rounded). Discs are rotation-tight — no box wrapping —
so mode phases can be built by binary exponentiation (11 + 21 doubling
levels) of two base rotations W = e^{i b1}, U = e^{-i b2}: the only
transcendental evaluations per field call are cos/sin of b1, b2 midpoints
and the z-factor cot/csc, all via mpmath.iv at prec 120 converted outward
to float64. The per-mode sums for f, f', f'', f''' run vectorized in numpy
float64 (IEEE-754 correctly-rounded +,-,*,/) with radius arrays carried
alongside; scalar real intervals use float64 pairs with one math.nextafter
outward step per endpoint per operation.

**Integrator.** Backward step s0 -> s1 = s0 - h via the exact identity
X(s1) = X0 - h v(X0,s0) + int (u - s1) D(X(u),u) du, D = dt_v + (v.grad)v
(needs only f...f''': dt_vx = Re(G3)/2, dt_vz = -(n_y/2) Im(G3), G3 =
f'''/f - (f'/f)(f''/f); grad v from G2 = f''/f - (f'/f)^2 and the analytic
csc^2 term), giving X(s1) in X0 - h v(X0,s0) + (h^2/2) hull(D(B, I_s))
for any a-priori box B containing the step segment. B comes from a
verified Picard inclusion Y0 + [0,h].(-v(B,I_s)) subset B (checked in
interval arithmetic, inflated on failure — never needed). The
X0-dependence propagates in Lohner QR mean-value form: X0 - c in A q,
X(s1) in [c - h v(c,s0) + R2] + (I - h J(B,I_s)) A q, with the interval
Jacobian J (componentwise mean-value theorem), per-step QR re-factoring
(rigorous 2x2 adjugate interval inverse), and intersection with the
direct box evaluation. Stop when lower(X_x) > d + delta. Step sizes
(exact binary floats, so s_k = t0 - k h is one rounding): 2^-14 for the
four interior extras, 2^-16 for the two deep paths, 2^-19 for the
wall-adjacent immediate exit — fixed from the float reconnaissance
(amplification integrals, Jacobian scales) BEFORE the production run.
hbar = m = 1.

## Honest scope + trust statement (G5)

TRUST BASE: (i) IEEE-754 float64 semantics of numpy/CPython arithmetic
(+,-,*,/,sqrt correctly rounded to nearest); (ii) math.nextafter
correctness; (iii) mpmath.iv outward-rounding correctness for pi, sin,
cos at prec 120 (used only for the two base rotations and the z-factor
per evaluation); (iv) gmpy2/MPFR correct rounding — used ONLY in G1
reference checks, never in a certificate; (v) the float64 spectral
coefficients C_k and the path data (t0, z0, n_y, d) are EXACT model data.
GIVEN this base, every certificate above is a machine-checked inequality:
the true backward trajectory of the exact Nx = 2048 band-limited field
provably pre-crosses the detector line. What this adds to P2: the
certified paths' pre-crossing statements no longer rest on converged-but-
unenclosed numerics — each carries a rigorous enclosure with printed
widths. What CAPD/COSY-grade work would still add: a formally verified
interval kernel replacing (i)-(iii), higher-order Taylor models to push
widths down orders of magnitude, and certification of ALL 64,000 kill-bin
paths rather than 7 witnesses. The certified object remains the
DISCRETIZED field (its distance to the continuum is the separately
quantified layer, as in O1/P2). The analytic T4-W5 proof (no computation
at all) remains the terminal open item.
Within-model; nothing here bears on nature.

## Caveats

- The worst-margin path's certificate width (7.84e-4) passes the 1e-3
  gate with only 1.27x headroom. The h-roster was fixed before production
  (no tuning); h = 2^-17 would roughly halve the width at double the
  cost, recorded here for any successor workstream.
- Seven paths are certified — witnesses, not the population: the
  remaining ~63,993 kill-bin points keep their O1/P2 classification at
  precision-certification grade only. The seven were selected (from the
  float recon, before any interval run) to cover both kill bins, both
  |n_y| families, the wall-adjacent immediate-exit class, the shallowest
  genuinely small-margin path, and the global worst margin.
- delta = 1e-3 matches P2's s_pc convention (crossing of d + 1e-3); the
  certificates therefore witness clearance ABOVE d + 1e-3, strictly
  stronger than X_x >= d.
- The G1 reference values are 200-bit MPFR point sums — "exact-point
  evaluations" per the roadmap wording; float64 numpy sums (rounding
  ~1e-15) would NOT generally lie inside the ~1e-11 thin enclosures and
  were not the registered comparison.
- Mirror families (-0.75, -1.0) inherit by the z -> 1-z symmetry
  re-verified in P2 (3.5e-8 level, classification-exact); not re-certified
  here.
- The two deep certificates take the tube through x ~ -1.26 (232) /
  x ~ 0.005 (395) — regions where |f|^2 stays >= 0.022 (monitored
  rigorously); no certificate approaches the wall lane except 212, whose
  10-step tube stays at z = 0.0175 with csc^2-scale Jacobians handled at
  h = 2^-19.
- Checkpoint/resume machinery was in place (20k-step atomic checkpoints
  per path) but no restart occurred during production; the g1 and certify
  phases ran concurrently under nohup.
