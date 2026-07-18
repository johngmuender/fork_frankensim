#!/usr/bin/env python3
"""
S1 -- P-acoustic adoption pilot: pricing cost (1) of the F-R14 sound cell.
[Phase S, Tier 8; ROADMAP_v10_ALTERNATIVES.md workstream S1]

Prices, symbolically and exactly, O3's printed cost (1) of adopting the
P-acoustic quantization-measure weight w = -1/4 (the only computationally
sound cell of the F-R14 decision matrix, tier7-program/O3/RESULTS.md sec. 3):

    "(1) reopens Sec. III's flat-measure [DF] derivations (Fisher/Madelung,
     Wallstrom, Born rates) on curved backgrounds"

by asking whether the Sec. III.B Fisher -> Bohm quantum-potential identity
(corpus3/01-GUM-Omega-Paper-v4.3-ext.md, eqs. (3.2)-(3.4), Theorem III.1)
survives on a curved 3-metric g3 = (1+h) delta with quantization measure
(det g3)^w, at linear order in h, for general symbolic w.

CAMPAIGN-UNITS CONVENTION (as printed in the corpus, Sec. III.B):
    E_Q[rho] = (hbar^2/8m) I_F[rho],  I_F = int |grad rho|^2 / rho d^3x   (3.2)
    delta E_Q/delta rho = (hbar^2/8m)[|grad rho|^2/rho^2 - 2 lap rho/rho] (3.3)
                        = -(hbar^2/2m) lap(sqrt rho)/sqrt rho = U_Q       (3.4)
i.e. Fisher stiffness hbar^2/8m, Bohm coefficient hbar^2/2m (Theorem III.1).
Measure-weight convention as in o3_matrix.py: det g3 = (1+h)^3, weight
(det g3)^w = (1+h)^{3w} (the "+3w shift" bookkeeping).

FORMALIZATION NOTE (coordinator-owned, printed): the corpus prints NO curved
Fisher functional; any curved extension requires a formalization choice for
where the weight (det g3)^w enters.  Both natural readings are computed:

  SCHEME A (covariant kinetic term, weighted quantization measure):
    E_A[rho] = (hbar^2/8m) int g3^{ij} d_i rho d_j rho / rho  (det g3)^w d^3x,
    chemical potential conjugate to rho under the quantization measure
    d mu_w = (det g3)^w d^3x  (so w = 1/2 is the Riemannian measure).
  SCHEME B (the printed flat functional applied verbatim to the weighted
    density):  rho = (det g3)^w sigma, sigma scalar;
    E_B[sigma] = (hbar^2/8m) int |grad rho|^2/rho d^3x  (flat contraction,
    coordinate measure, exactly as printed in (3.2)); chemical potential
    conjugate to sigma under d^3x.

  Target in both schemes: the covariant Bohm form
    U_cov = -(hbar^2/2m) Lap_g sqrt(f)/sqrt(f),
    Lap_g = (det g3)^{-1/2} d_i ( (det g3)^{1/2} g3^{ij} d_j . ),
  with f the density the scheme varies (rho for A, sigma for B).
  The CORRECTION is  C = (chemical potential) - U_cov, expanded to O(h)
  (h -> eps*h(x); |grad h|^2 terms are O(eps^2) and drop automatically).

Deterministic, symbolic (sympy), nothing fitted; all checks are exact
symbolic identities (simplify(...) == 0) plus one numeric O(eps^2)
scaling spot-check.  Runtime ~ tens of seconds.

Outputs: s1_results.json, s1_fig.png (correction coefficients vs w).
"""

import json
import time
import sympy as sp

T0 = time.time()
checks = []


def check(name, passed, detail=""):
    checks.append({"name": name, "pass": bool(passed), "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}  {detail}")
    return bool(passed)


# ---------------------------------------------------------------- setup
x, y, z = sp.symbols("x y z", real=True)
X = (x, y, z)
hbar, m, w, eps = sp.symbols("hbar m w epsilon", positive=True, real=True)
# w is real, not necessarily positive:
w = sp.Symbol("w", real=True)

rho = sp.Function("rho", positive=True)(x, y, z)
sig = sp.Function("sigma", positive=True)(x, y, z)
hf = sp.Function("h", real=True)(x, y, z)

pref = hbar**2 / (8 * m)          # printed Fisher stiffness (3.2)
bohm = hbar**2 / (2 * m)          # printed Bohm coefficient (3.4)


def grad(f):
    return [sp.diff(f, xi) for xi in X]


def dot(a, b):
    return sum(ai * bi for ai, bi in zip(a, b))


def lap(f):
    return sum(sp.diff(f, xi, 2) for xi in X)


def EL(L, f):
    """Euler-Lagrange delta/delta f of int L d^3x, L = L(f, grad f)."""
    return sp.diff(L, f) - sum(sp.diff(sp.diff(L, sp.diff(f, xi)), xi)
                               for xi in X)


def is_zero(expr):
    return sp.simplify(sp.together(sp.expand(expr))) == 0


# =================================================================
# (a) S1-G1: the flat-measure identity, exact
# =================================================================
print("\n== (a) flat-measure Fisher -> Bohm identity (Sec. III.B, "
      "(3.2)-(3.4)) ==")

L_flat = dot(grad(rho), grad(rho)) / rho              # I_F integrand
dEQ = pref * EL(L_flat, rho)                          # delta E_Q / delta rho
form_33 = pref * (dot(grad(rho), grad(rho)) / rho**2 - 2 * lap(rho) / rho)
U_Q = -bohm * lap(sp.sqrt(rho)) / sp.sqrt(rho)        # Bohm quantum potential

check("S1-G1a  delta E_Q/delta rho equals printed intermediate form (3.3)",
      is_zero(dEQ - form_33), "coefficient hbar^2/8m")
check("S1-G1b  (3.3) == -(hbar^2/2m) lap(sqrt rho)/sqrt rho  (3.4), exact",
      is_zero(dEQ - U_Q), "Theorem III.1 identity, exact symbolic zero")

# =================================================================
# (b) curved 3-metric g3 = (1+h) delta, weight (det g3)^w, general w
# =================================================================
print("\n== (b) curved background, general symbolic w ==")
Om = 1 + eps * hf                                     # conformal factor (1+h)
detg = Om**3                                          # det g3 (o3 convention)
sqrtg = Om**sp.Rational(3, 2)


def lap_g(f):
    """Laplace-Beltrami on g3 = Om*delta."""
    return (1 / sqrtg) * sum(
        sp.diff(sqrtg * (1 / Om) * sp.diff(f, xi), xi) for xi in X)


def lin(expr):
    """Expand to linear order in eps (h and its derivatives first order)."""
    return sp.expand(sp.series(sp.together(expr), eps, 0, 2).removeO())


# ---------------- SCHEME A: covariant kinetic, weighted measure ----------
L_A = detg**w * (1 / Om) * dot(grad(rho), grad(rho)) / rho
U_A = sp.together(detg**(-w) * pref * EL(L_A, rho))   # conj. to rho, d mu_w
U_cov_rho = -bohm * lap_g(sp.sqrt(rho)) / sp.sqrt(rho)
C_A_exact = sp.simplify(sp.together(U_A - U_cov_rho))
C_A_lin = sp.simplify(lin(C_A_exact))

# predicted structure: (hbar^2/8m) (3 - 6w) eps grad h . grad rho / rho
gradh_gradrho = dot(grad(hf), grad(rho)) / rho
coefA = sp.Rational(3, 1) - 6 * w
C_A_pred = pref * coefA * eps * gradh_gradrho
check("S1-G2a  scheme A linear correction == (hbar^2/8m)(3-6w) "
      "grad h.grad rho/rho",
      is_zero(C_A_lin - C_A_pred), f"coefficient (3 - 6w) x hbar^2/8m")

# constant h: correction vanishes exactly, all orders, all w
h0 = sp.Symbol("h0", positive=True)
C_A_const = C_A_exact.subs(hf, h0)   # derivatives of const vanish on subs
check("S1-G2b  scheme A, constant h: correction == 0 EXACTLY (all orders, "
      "general w)", is_zero(C_A_const),
      "measure factor cancels; only the covariant g^ij contraction remains")

# cancellation: solve (3 - 6w) = 0
solA = sp.solve(sp.Eq(coefA, 0), w)
check("S1-G3a  scheme A cancellation solved exactly: w = 1/2 (Riemannian "
      "measure), unique", solA == [sp.Rational(1, 2)], f"solve -> {solA}")
check("S1-G3b  scheme A at w=1/2: correction == 0 EXACTLY (all orders in h)",
      is_zero(C_A_exact.subs(w, sp.Rational(1, 2))),
      "w=1/2 is the fully covariant functional")

# rescue-window exclusion (O3: w in (-5/18, -2/9))
wlo, whi = sp.Rational(-5, 18), sp.Rational(-2, 9)
in_window = bool(wlo < sp.Rational(1, 2)) and bool(sp.Rational(1, 2) < whi)
check("S1-G3c  cancelling w = 1/2 is OUTSIDE the F-R14 rescue window "
      "(-5/18, -2/9)", not in_window,
      "no w both rescues F-R14 and preserves the identity (scheme A)")

# ---------------- SCHEME B: printed flat functional, weighted density ----
rho_B = detg**w * sig
L_B = dot(grad(rho_B), grad(rho_B)) / rho_B           # (3.2) verbatim
U_B = pref * EL(sp.together(sp.expand(L_B)), sig)     # conj. to sigma, d^3x
U_cov_sig = -bohm * lap_g(sp.sqrt(sig)) / sp.sqrt(sig)
C_B_lin = sp.simplify(lin(sp.together(U_B - U_cov_sig)))

# predicted structure: three independent h-structures
core = dot(grad(sig), grad(sig)) / sig**2 - 2 * lap(sig) / sig  # flat U_Q core
gradh_gradsig = dot(grad(hf), grad(sig)) / sig
cB1, cB2, cB3 = 3 * w + 1, 1 - 6 * w, -6 * w
C_B_pred = pref * eps * (cB1 * hf * core + cB2 * gradh_gradsig
                         + cB3 * lap(hf))
check("S1-G2c  scheme B linear correction == (hbar^2/8m)[(3w+1) h*core + "
      "(1-6w) grad h.grad sig/sig - 6w lap h]",
      is_zero(C_B_lin - C_B_pred),
      "core = |grad sig|^2/sig^2 - 2 lap sig/sig")

# constant h: only the multiplicative (3w+1) h * core piece survives
C_B_const = sp.simplify(lin(sp.together(U_B - U_cov_sig)).subs(
    [(sp.diff(hf, xi, 2), 0) for xi in X]
    + [(sp.Derivative(hf, xi, yi), 0) for xi in X for yi in X]
    + [(sp.diff(hf, xi), 0) for xi in X]))
check("S1-G2d  scheme B, constant h: correction == (3w+1) h x (flat U_Q), "
      "purely multiplicative",
      is_zero(C_B_const - pref * eps * cB1 * hf * core),
      "absorbable into the hbar^2/8m stiffness iff position-independent")

# cancellation: all three coefficients simultaneously?
solB = sp.solve([sp.Eq(cB1, 0), sp.Eq(cB2, 0), sp.Eq(cB3, 0)], w)
check("S1-G3d  scheme B: NO w cancels all three structures "
      "(zeros at w=-1/3, 1/6, 0 are distinct)",
      solB == [] and sp.solve(sp.Eq(cB1, 0), w) == [sp.Rational(-1, 3)]
      and sp.solve(sp.Eq(cB2, 0), w) == [sp.Rational(1, 6)]
      and sp.solve(sp.Eq(cB3, 0), w) == [0],
      f"simultaneous solve -> {solB} (empty)")

# =================================================================
# (c) evaluation at the printed w values
# =================================================================
print("\n== (c) evaluation at w = 0, -1/4 (P-acoustic), -1/2 ==")
w_vals = {"0": sp.Integer(0), "-1/4": sp.Rational(-1, 4),
          "-1/2": sp.Rational(-1, 2)}
tableA, tableB = {}, {}
for name, wv in w_vals.items():
    a = coefA.subs(w, wv)
    b = (cB1.subs(w, wv), cB2.subs(w, wv), cB3.subs(w, wv))
    tableA[name] = str(a)
    tableB[name] = [str(t) for t in b]
    print(f"  w = {name:>4}: scheme A coeff (x hbar^2/8m on grad h.grad "
          f"rho/rho) = {a}")
    print(f"            scheme B coeffs (x hbar^2/8m on [h*core, "
          f"grad h.grad sig/sig, lap h]) = {b}")

check("S1-G3e  w = -1/4 (P-acoustic): scheme A coefficient = 9/2 "
      "(i.e. 9 hbar^2/16m)", coefA.subs(w, sp.Rational(-1, 4))
      == sp.Rational(9, 2), "C_A = (9 hbar^2/16 m) grad h . grad rho / rho")
check("S1-G3f  w = -1/4 (P-acoustic): scheme B coefficients = "
      "(1/4, 5/2, 3/2)",
      (cB1.subs(w, sp.Rational(-1, 4)), cB2.subs(w, sp.Rational(-1, 4)),
       cB3.subs(w, sp.Rational(-1, 4)))
      == (sp.Rational(1, 4), sp.Rational(5, 2), sp.Rational(3, 2)), "")
check("S1-G3g  neither printed convention (w=0, -1/2) nor P-acoustic "
      "(w=-1/4) cancels either scheme",
      all(coefA.subs(w, wv) != 0 for wv in w_vals.values())
      and all(any(c.subs(w, wv) != 0 for c in (cB1, cB2, cB3))
              for wv in w_vals.values()), "")

# =================================================================
# numeric spot-check: linear-order formula is O(eps^2)-accurate
# =================================================================
print("\n== numeric O(eps^2) spot-check (scheme A) ==")
subs_f = {rho: sp.exp(-(x**2 + 2 * y**2 + 3 * z**2) / 2) * (1 + x * y / 10),
          hf: sp.sin(x + 2 * y - z) / 20}
pt = {x: sp.Rational(1, 3), y: sp.Rational(-1, 5), z: sp.Rational(1, 7),
      hbar: 1, m: 1, w: sp.Rational(-1, 4)}
resid = sp.together(C_A_exact - C_A_pred)  # should be O(eps^2)
rf = sp.lambdify(eps, resid.subs(subs_f).doit().subs(pt), "math")
r1, r2 = abs(rf(1e-3)), abs(rf(2e-3))
ratio = r2 / r1 if r1 > 0 else float("inf")
check("numeric: residual (exact - linear) scales as eps^2 "
      "(ratio at eps doubling in [3.5, 4.5])", 3.5 < ratio < 4.5,
      f"|r(2e-3)|/|r(1e-3)| = {ratio:.3f}")

# =================================================================
# figure: correction coefficients vs w
# =================================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

wgrid = np.linspace(-0.6, 0.6, 601)
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
C1, C2, C3, C4 = "#2a78d6", "#008300", "#e87ba4", "#eda100"

fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=150)
fig.patch.set_facecolor(SURF)
ax.set_facecolor(SURF)
ax.axhline(0, color=INK2, lw=0.8, alpha=0.5)
ax.axvspan(-5 / 18, -2 / 9, color=C4, alpha=0.18, lw=0,
           label="F-R14 rescue window (-5/18, -2/9)")
ax.plot(wgrid, 3 - 6 * wgrid, color=C1, lw=2,
        label=r"scheme A: $\nabla h\!\cdot\!\nabla\rho/\rho$  (3-6w)")
ax.plot(wgrid, 3 * wgrid + 1, color=C2, lw=2,
        label=r"scheme B: $h\times$core  (3w+1)")
ax.plot(wgrid, 1 - 6 * wgrid, color=C3, lw=2,
        label=r"scheme B: $\nabla h\!\cdot\!\nabla\sigma/\sigma$  (1-6w)")
ax.plot(wgrid, -6 * wgrid, color=C4, lw=2,
        label=r"scheme B: $\nabla^2 h$  (-6w)")
for wv, lbl in [(0.5, "w=1/2\n(A cancels)"), (-0.25, "w=-1/4\n(P-acoustic)"),
                (0.0, "w=0"), (-0.5, "w=-1/2")]:
    ax.axvline(wv, color=INK2, lw=0.7, ls=":", alpha=0.6)
    ax.text(wv, 6.55, lbl, ha="center", va="top", fontsize=7.5, color=INK2)
ax.plot([0.5], [0], marker="o", ms=8, mfc=C1, mec=SURF, mew=2, zorder=5)
ax.plot([-0.25], [4.5], marker="o", ms=8, mfc=C1, mec=SURF, mew=2, zorder=5)
ax.annotate("9/2", (-0.25, 4.5), textcoords="offset points", xytext=(8, 6),
            fontsize=9, color=INK)
ax.set_xlabel("quantization-measure weight w   [$\\rho_0, J \\propto (\\det g_3)^w$]",
              color=INK)
ax.set_ylabel("correction coefficient  ($\\times\\, \\hbar^2/8m$, linear in h)",
              color=INK)
ax.set_title("S1: h-linear corrections to the Fisher$\\to$Bohm identity vs "
             "measure weight w\n(no w cancels scheme B; scheme A cancels only "
             "at w=1/2, outside the rescue window)", fontsize=10, color=INK)
ax.tick_params(colors=INK2)
for s in ax.spines.values():
    s.set_color(INK2); s.set_alpha(0.4)
ax.grid(True, color=INK2, alpha=0.12, lw=0.6)
ax.legend(fontsize=8, loc="lower left", framealpha=0.9)
ax.set_xlim(-0.6, 0.6)
ax.set_ylim(-4.7, 6.7)
fig.tight_layout()
import os
OUT = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(OUT, "s1_fig.png"), facecolor=SURF)
print("figure written: s1_fig.png")

# =================================================================
# results JSON
# =================================================================
npass = sum(c["pass"] for c in checks)
results = {
    "workstream": "S1 (Tier 8, Phase S) -- P-acoustic adoption pilot: "
                  "pricing O3 cost (1)",
    "date": "2026-07-18",
    "conventions": {
        "campaign_units": "E_Q = (hbar^2/8m) I_F (corpus (3.2)); "
                          "U_Q = -(hbar^2/2m) lap(sqrt rho)/sqrt rho ((3.4), "
                          "Theorem III.1)",
        "metric": "g3 = (1+h) delta_ij; det g3 = (1+h)^3; weight (det g3)^w "
                  "= (1+h)^{3w} (o3_matrix.py bookkeeping)",
        "order": "linear in h (h -> eps h(x); |grad h|^2 = O(eps^2) dropped "
                 "by series, printed truncation O(eps^2))",
    },
    "formalization_amendment": "Coordinator-owned: corpus prints no curved "
        "Fisher functional; both natural formalizations computed (scheme A: "
        "covariant kinetic + weighted measure; scheme B: printed flat "
        "functional verbatim on the weighted density). Spec's single "
        "'curved-measure variation' is thus delivered as a two-scheme "
        "robustness pair; conclusions agree.",
    "flat_identity": {
        "functional": "E_Q = (hbar^2/8m) int |grad rho|^2/rho d^3x",
        "variation": "(hbar^2/8m)[|grad rho|^2/rho^2 - 2 lap rho/rho]",
        "bohm_form": "-(hbar^2/2m) lap(sqrt rho)/sqrt rho",
        "identity_exact": True,
    },
    "scheme_A": {
        "definition": "E = (hbar^2/8m) int g^{ij} d_i rho d_j rho /rho "
                      "(det g)^w d^3x; chemical potential wrt d mu_w = "
                      "(det g)^w d^3x; target U_cov = -(hbar^2/2m) "
                      "Lap_g sqrt(rho)/sqrt(rho)",
        "correction_linear": "(hbar^2/8m) (3 - 6w) grad h . grad rho / rho",
        "constant_h": "zero exactly, all orders, all w",
        "cancellation_w": "1/2 (unique; exact to all orders; = Riemannian "
                          "measure)",
        "cancellation_in_rescue_window": False,
        "at_w": {"0": "3", "-1/4": "9/2", "-1/2": "6",
                 "units": "x hbar^2/8m on grad h.grad rho/rho"},
        "at_w_minus_quarter": "C = (9 hbar^2/16 m) grad h . grad rho / rho",
    },
    "scheme_B": {
        "definition": "rho = (det g)^w sigma; E = (hbar^2/8m) int "
                      "|grad rho|^2/rho d^3x ((3.2) verbatim); chemical "
                      "potential wrt sigma, d^3x; target U_cov in sigma",
        "correction_linear": "(hbar^2/8m)[ (3w+1) h (|grad sig|^2/sig^2 - "
                             "2 lap sig/sig) + (1-6w) grad h.grad sig/sig "
                             "- 6w lap h ]",
        "constant_h": "(3w+1) h x (flat U_Q): multiplicative rescaling, "
                      "absorbable into stiffness; zero only at w=-1/3",
        "cancellation_w": None,
        "individual_zeros": {"h*core": "-1/3", "grad-grad": "1/6",
                             "lap h": "0"},
        "at_w": {"0": ["1", "1", "0"], "-1/4": ["1/4", "5/2", "3/2"],
                 "-1/2": ["-1/2", "4", "3"],
                 "units": "x hbar^2/8m on [h*core, grad h.grad sig/sig, "
                          "lap h]"},
    },
    "pricing_bottom_line": {
        "cost_1_priced": True,
        "statement": "No w cancels the h-gradient corrections in either "
            "formalization at any printed value {0, -1/4, -1/2}; the unique "
            "scheme-A cancellation w = 1/2 lies outside the F-R14 rescue "
            "window (-5/18, -2/9), so no weight both rescues F-R14 and "
            "preserves the Fisher->Bohm identity. At w = -1/4 the correction "
            "is a genuine new force term ~ (9 hbar^2/16 m) grad h.grad "
            "rho/rho (scheme A). Effect on the [DF] label: Theorem III.1 "
            "remains [DF] on flat backgrounds (h = 0, and exactly for "
            "constant h in scheme A); on curved backgrounds at w = -1/4 the "
            "chain does NOT close as printed -- the curved extension would "
            "need a new derivation and is not [DF] as it stands.",
        "no_adoption_recommended": True,
        "fr14_regrade_language": "unchanged",
    },
    "checks": checks,
    "checks_passed": f"{npass}/{len(checks)}",
    "runtime_s": round(time.time() - T0, 2),
}
with open(os.path.join(OUT, "s1_results.json"), "w") as f:
    json.dump(results, f, indent=1)

print(f"\n{npass}/{len(checks)} checks PASS; runtime "
      f"{results['runtime_s']} s; s1_results.json written")
if npass != len(checks):
    raise SystemExit(1)
