# Tier-2b STEP 1 — 1-D radial (hedgehog) field-level solve at eps = 0.05

Reproduce with `python3 radial_solve.py` (deterministic, ~1 s wall time; writes
`radial_profile.png`, `radial_results.json`, `radial_run.log` holds the printed
table). Companion files: `radial_solve.py`, `radial_profile.png`,
`radial_results.json`, `radial_run.log`.

## Setup and conventions

Hedgehog `P = exp(i f(r) x_hat.sigma)`, `f(0) = pi`, `f(inf) = 0`. Total energy
`E = a2 E2 + a4 E4 + a6 E6 + a0 E0` with the radial sector integrals as
specified (E2 sigma-model, E4 Skyrme, E6 sextic, E0 potential with c2 = 0).
Fixed `a6 = a0 = 1`, `m = 1`; dial `a2 = a4 = t` tuned so that at the solved
profile `ratio = t (E2+E4)/(E6+E0) = 0.05`.

**Solver** (documented in the script header): midpoint discretization of the
energy on a graded grid (N = 4000 cells, refined to 2N = 8000; node density
boosted near r = 0, near the former compacton edge r ~ R\*, and through the
tail window), `R_max = 6`; hard Dirichlet `f = pi, 0` at the endpoints.
Minimization: preconditioned energy-monotone gradient flow into the basin,
then damped (Levenberg) Newton with the *exact* tridiagonal Hessian of the
discretized energy (Thomas solves), run to `max|grad| < 1e-12` (achieved
~4e-13 at N, ~7e-13 at 2N; last accepted Newton step changed E by < 1e-16,
i.e. the < 1e-12-per-step gate is met at machine floor). Outer fixed-point
iteration on t converged in 3 solves:
t = 0.008 -> 0.0082510 -> **t = a2 = a4 = 0.00827643494930** (ratio 0.0499847).

The eps = 0 reference is the BPS compacton `f0 = 2 arccos(r/R*)` with
`R* = 2^(5/6) = 1.781797...` in this normalization (derived below). At
eps = 0.05 the edge is smoothed into a Yukawa (massive-dipole) tail; the graded
grid puts ~40% of its nodes within +-3 widths of the edge region.
`f` is monotone decreasing at both resolutions; `max f` beyond `0.97 R_max` is
2.5e-15 < 1e-8 (tail-containment check passed, no R_max extension needed).

## Energy accounting (values at 2N = 8000; delta = |N - 2N| refinement shift)

| quantity | value (2N) | refinement delta |
|---|---|---|
| a2 = a4 = t | 0.00827643494930 | (held fixed across N, 2N) |
| E2 | 12.1174455 | 1.4e-05 |
| E4 | 6.1911247 | 1.7e-06 |
| E6 | 1.5239360 | 4.5e-07 |
| E0 | 1.5075871 | 1.4e-06 |
| E_total | 3.1830528 | 1.1e-06 |
| eps ratio t(E2+E4)/(E6+E0) | **0.0499847** | 2.6e-08 |
| degree K | 1.00000007 | 2.2e-07 |
| virial residual (rel.) | 5.96e-07 | 1.8e-06 (N value 2.38e-06) |
| max EL residual (grad/weight) | 1.9e-09 | — |
| mu_fit | 7.772564 | 6.2e-05 |
| A_fit | 5.2333e+06 | 7.5e+02 (1.4e-04 rel.) |
| BPS deficit (E6+E0)/bound - 1 | 4.8174e-03 | 3.2e-07 |

All refinement deltas are far below the gate tolerances — the honest error bar
on every gated quantity is at the 1e-6 level or better.

## Gate 1 — Degree

`K = -(2/pi) INT sin^2 f f' dr` (signed; f is monotone decreasing, verified),
midpoint quadrature on the solved profile:

**K = 1.0000003 (N = 4000), 1.0000001 (2N = 8000)** — i.e. `K = 1.0000` to
better than 4 decimals, converging toward 1 under refinement as O(h^2).

**Verdict: PASS.**

## Gate 2 — Derrick virial

Derivation (own, as instructed). Under `f_lambda(r) = f(r/lambda)`:
`E2 -> lambda E2` (the `f'^2 r^2` and `2 sin^2 f` pieces each pick up one net
power of lambda), `E4 -> lambda^{-1} E4`, `E6 -> lambda^{-3} E6`
(`sin^4 f (f'^2/lambda^2)/(lambda^2 r^2) * lambda dr`), `E0 -> lambda^3 E0`.
Stationarity `d/dlambda|_1 = 0` of `lambda a2 E2 + lambda^{-1} a4 E4 +
lambda^{-3} a6 E6 + lambda^3 a0 E0` gives, **with the sector weights
included**:

```
a2 E2 - a4 E4 - 3 a6 E6 + 3 a0 E0 = 0
```

The paper's unweighted form `E2 - E4 - 3 E6 + 3 E0 = 0` only coincides with
this when all weights are equal; at eps = 0.05 the weighted identity is the
correct one (the unweighted combination evaluates to 1.42 here — clearly not a
solution property — so the paper's statement must be read in its own
sector-weighting convention).

Measured relative residual `|t E2 - t E4 - 3 E6 + 3 E0| / E_total`:

**2.38e-06 (N), 5.96e-07 (2N)** — two to three orders below the corpus
acceptance of 3e-4, and shrinking ~4x under refinement (pure O(h^2)
quadrature error, as expected for the exact discrete minimizer).

**Verdict: PASS.**

## Gate 3 — Tail (massive l = 1 Yukawa)

Linearization (own derivation): for small f the quadratic energy is
`INT [a2 (f'^2 r^2 + 2 f^2) + (a0 m^2/2) f^2 r^2] dr` (note
`1 - cos f = f^2/2 + O(f^4)`), whose EL equation is

```
a2 (f'' + 2 f'/r - 2 f/r^2) = (a0 m^2 / 2) f   =>   mu = m sqrt(a0 / (2 a2))
```

i.e. `f ~ A k_1(mu r) = A (1/(mu r) + 1/(mu r)^2) e^{-mu r}`.
**Convention caveat:** the task sheet's quoted `mu_analytic = m sqrt(a0/a2)`
omits the 1/2 from the potential expansion; the data unambiguously
adjudicates (below).

Fit of `f = A (1/(mu r)^2 + 1/(mu r)) e^{-mu r}` over the window
`f in [1e-6, 1e-2]` (r in [2.224, 3.353], 809 points at N / 1618 at 2N;
log-space least squares, golden-section in mu):

| | fitted | derived m/sqrt(2 a2) | task-sheet m/sqrt(a2) |
|---|---|---|---|
| mu | **7.772564 (2N)** | 7.772547 | 10.992041 |

Agreement with the derived mu is 2.2e-6 relative (and improves 8x under
refinement); rms log-residual of the fit is 2.0e-6 over four decades of f —
the tail is an essentially perfect massive dipole. The sqrt(2)-larger
task-sheet value is excluded.

Amplitude: **A = 5.2333e+06** (refinement shift 1.4e-04 relative). A is large
because k_1 normalization hides `e^{-mu R*} ~ 1e-6`; natural dimensionless
combinations (m = 1 makes everything dimensionless already):

| combination | value |
|---|---|
| A | 5.2333e+06 |
| A mu | 4.0676e+07 |
| A mu^2 | 3.1616e+08 |
| A / mu | 6.7331e+05 |
| A / mu^2 | 8.6626e+04 |
| A e^{-mu R*} | 5.0603 |
| A e^{-mu R*} / mu | 0.6510 |
| A e^{-mu R*} / (mu R*) | 0.3654 |
| A e^{-mu R*} / (2 pi) | 0.8054 |
| fit value at R\* (= A k_1(mu R\*)) | 0.3918 |

Comparison to the corpus's `p = 0.84 +- 0.03`: **none of the natural
combinations lands inside the quoted band.** The closest is
`A e^{-mu R*}/(2 pi) = 0.805` (just below the band's lower edge, ~1.2 corpus
sigma), with `A e^{-mu R*}/mu = 0.651` next. Since the corpus's normalization
convention for p is not recoverable from its text, this cannot be scored as a
numerical pass or fail; what our solve pins down without ambiguity is mu and
the full profile, and the fact that no order-unity convention we constructed
reproduces 0.84 exactly is itself an informative (negative) datum for
reverse-engineering the corpus convention.

**Verdict: PASS on the functional form and on mu (against the correctly
derived linearization); CONVENTION-LIMITED on the amplitude normalization.**

## Gate 4 — Energy accounting and BPS bound

Bogomolny derivation (own, this normalization): with
`1 - cos f = 2 sin^2(f/2)`,

```
E6 + E0 = INT [ (sqrt(a6) sin^2 f f'/r  +-  sqrt(a0) m sqrt(2) sin(f/2) r)^2
               -+ 2 sqrt(a6 a0) sqrt(2) m sin^2 f sin(f/2) f' ] dr
        >= 2 sqrt(a6 a0) sqrt(2) m | INT_0^pi sin^2 f sin(f/2) df |
         = 2 sqrt(2) m * (16/15)          [ INT_0^pi sin^2 f sin(f/2) df = 16/15 ]
         = 32 sqrt(2) / 15 = 3.0169889...   (a6 = a0 = m = 1)
```

(The task sheet's `64/(15 pi) * (2m) * ...` hint evaluates to a different
constant and evidently carries the corpus's angular/4pi bookkeeping; the bound
above is the one consistent with the radial integrals as specified here.)
Saturation (`sin^2 f f'/r = -sqrt(2) m sin(f/2) r`) integrates in closed form
to `cos(f/2) = r/R*`, i.e. exactly the compacton `f0 = 2 arccos(r/R*)` with
`R* = (4 sqrt(2))^{1/3} = 2^{5/6} = 1.7818` — this is the R\* used in the
figure.

Measured: `E6 + E0 = 3.0315231`, bound `= 3.0169889`, so

**BPS deficit (E6+E0)/bound - 1 = 4.82e-03 (delta 3.2e-07)** — the eps = 0.05
solution sits only ~0.5% off-BPS in its BPS sector, while the suppressed
sector carries `t (E2+E4) = 0.15151` (4.76% of E_total = 3.1830528). Achieved
dial: **eps ratio = 0.049985** (target 0.05 +- 0.001).

**Verdict: PASS** (accounting closes: `t(E2+E4) + E6 + E0 = E_total` to
1e-12 by construction; ratio within tolerance; deficit small and positive as
required by the bound).

## Figure

`radial_profile.png`: f(r) at eps = 0.05 (blue) vs the eps = 0 compacton
`2 arccos(r/R*)` (red dashed), with `R* = 2^(5/6)` from this normalization;
the compacton's square-root edge cusp at R\* is smoothed into the Yukawa tail.
Log-scale inset: the tail over 9 decades with the fitted `A k_1(mu r)`
(violet dotted) through the shaded fit window `f in [1e-6, 1e-2]`.

## Summary of verdicts

| gate | result | verdict |
|---|---|---|
| 1. Degree | K = 1.0000001 (2N) | **PASS** |
| 2. Virial (weighted, own derivation) | rel. residual 6.0e-07 vs gate 3e-4 | **PASS** |
| 3. Tail form + mu | mu_fit/mu_derived - 1 = 2.2e-6; rms 2e-6 | **PASS** |
| 3b. Tail amplitude vs corpus p = 0.84(3) | best candidate 0.805 (A e^{-mu R*}/2pi) | **CONVENTION-LIMITED** |
| 4. Energy/BPS | ratio 0.049985; deficit +0.48% | **PASS** |
