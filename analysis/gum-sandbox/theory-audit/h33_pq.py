#!/usr/bin/env python3
"""
H3.3 Part B (T-H10): derive the E.6 response exponents (p, q) from the
corpus's OWN action (W2, Sec. II.C) + the defect-geometry map (Sec. VI.A),
and feed them into the h2.5 master form

    a1(p, q) = -(5p + q) / (12 (2p + q))     [per scalar mode]

(w_x ~ 1 + p h, w_t ~ 1 + q h under a weak static defect background h).

Textual constraints extracted (01-GUM-Omega-Paper-v2.0.1.md):
  * II.C: the moduli lambda, mu, mu_c, alpha, beta, gamma, J, rho0 are
    CONSTANTS (numbers of the medium; no defect-response functions printed).
  * VI.A [IM: KBKK]: disclination density IS the curvature of the material
    (Riemann-Cartan) connection => a weak static disclination background is a
    non-flat material metric g3_ij; probe class: conformal g3 = (1+h) delta.
  * VI.B: "slowly varying defect-geometric background"; displayed sector
    perturbation g_s = gbar + 2 delta_s c^2 n(x)n -- CONE-ONLY (time-time).
  * II.B (F10'): objectivity -- energy depends on rotations only through the
    frame-RELATIVE texture (the M-1 lesson) => the algebra index of the
    micro-rotation is a material-FRAME index (soldered reading favored).
  * h2.5: the medium's fields carry the FLAT functional measure (forced by
    the corpus's Sec. III construction); divergence-form operators.

What the corpus does NOT print: how W2's spatial contractions covariantize
(index variance of u and of the soldering), any moduli dressing, any measure
change on defected backgrounds. This script enumerates ALL textually
consistent readings and shows the sign of a1 per sector per reading.

Machine checks:
  M1  independent sympy re-derivation of the master form (Gilkey/Vassilevich
      on D = -w_t dtau^2 - d_i w_x d_i + M, 4D Euclidean): replicates h2.5.
  M2  anchors: unimodular ray q = -3p has E = O(h^2) and a1 = +1/6; minimal
      covariant Laplacian anchor; cone-only (0,1) -> -1/12; mass response r
      never enters the R-coefficient.
  M3  sector (p,q) from exponent counting on conformal g3 (sympy series).
  M4  the displacement invariant p = q - 1 (kinetic vs gradient differ by
      exactly one derivative-index contraction) for every uniform variance
      and every field density weight.
  M5  wedge closure: positivity wedge (5p+q)(2p+q) < 0 requires q in
      (2/3, 5/6) when p = q-1 -- no half-integer reachable => a1(B2) < 0 in
      EVERY reading; Q-sector wedge p in (-3/4, -3/10) contains only the
      unsoldered covector reading.
  M6  rescue window: density weight w (flat-measure field = (det g3)^w x
      field) must lie in (-5/18, -2/9); w = -1/4 = exact-acoustic/unimodular
      measure; natural weights 0, -1/2 excluded.

Deterministic; python3 + sympy + numpy. Appends to h33_summary.json.
"""

import itertools
import json
import pathlib
from fractions import Fraction

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
SUMMARY = HERE / "h33_summary.json"

checks = []


def check(name, passed, detail=""):
    checks.append({"name": name, "pass": bool(passed), "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}  {detail}")


# ======================================================================
# M1. Independent re-derivation of the master form.
#     D = -w_t d_tau^2 - d_i w_x d_i + M,  w_t = 1 + q e h(x),
#     w_x = 1 + p e h(x),  M = m0^2 (1 + rr e h(x));  flat L^2 measure.
#     Gilkey leading symbol: g^{munu} = diag(w_t, w_x, w_x, w_x).
#     a1-integrand = tr(E + R/6);  a1 := R-coefficient => master form.
# ======================================================================
tau, X, Y, Z = sp.symbols("tau x y z", real=True)
coords = [tau, X, Y, Z]
e = sp.Symbol("epsilon")
p, q, rr = sp.symbols("p q r", real=True)
m0 = sp.Symbol("m0", positive=True)
h = sp.Function("h")(X)

wt = 1 + q * e * h
wx = 1 + p * e * h
Mass = m0**2 * (1 + rr * e * h)

g_up = sp.diag(wt, wx, wx, wx)
g_dn = g_up.inv()


def lin(expr):
    """Expand to first order in e."""
    return sp.expand(sp.series(sp.expand(expr), e, 0, 2).removeO())


# Christoffels of g_dn (exact, then linearized)
Gamma = [[[0] * 4 for _ in range(4)] for _ in range(4)]
for l in range(4):
    for m in range(4):
        for n in range(4):
            s = 0
            for k in range(4):
                s += g_up[l, k] * (sp.diff(g_dn[k, m], coords[n])
                                   + sp.diff(g_dn[k, n], coords[m])
                                   - sp.diff(g_dn[m, n], coords[k])) / 2
            Gamma[l][m][n] = lin(s)

# Ricci: R^l_{m l n} contraction; R = g^{mn} Ric_{mn}
Ric = sp.zeros(4, 4)
for m in range(4):
    for n in range(4):
        s = 0
        for l in range(4):
            s += sp.diff(Gamma[l][m][n], coords[l]) \
                - sp.diff(Gamma[l][m][l], coords[n])
            for k in range(4):
                s += Gamma[l][l][k] * Gamma[k][m][n] \
                    - Gamma[l][n][k] * Gamma[k][m][l]
        Ric[m, n] = lin(s)
Rscal = lin(sum(g_up[m, n] * Ric[m, n] for m in range(4) for n in range(4)))

hpp = sp.diff(h, X, 2)
R_coeff = sp.simplify(Rscal.coeff(e).coeff(hpp))
# expected: R = -(q + 2p) h'' e + O(e^2) in THIS sign convention set by the
# check below (sphere-positive scalar curvature); h2.5 quotes |(q+2p)|.
check("M1a R[g_Gilkey] proportional to (q+2p) h'' at O(h) [sympy]",
      sp.simplify(sp.Abs(R_coeff) - sp.Abs(q + 2 * p)) == 0
      or sp.simplify(R_coeff**2 - (q + 2 * p)**2) == 0,
      f"R-coeff = {R_coeff}")

# sign-convention anchor: 4D conformal sphere-like check via round S^2 x R^2
# is heavy; instead anchor on the known conformal result: for
# g_mn = (1 + e H) delta (i.e. g^mn = (1 - e H) delta -> wt = wx = 1 - e H):
# geometer convention (sphere positive) gives R = -3 e Lap H + O(e^2).
H2 = sp.Function("H")(X)
Rconf = Rscal.subs([(q, -1), (p, -1)]).subs(h, H2).doit()
Rconf_coeff = sp.simplify(lin(Rconf).coeff(e).coeff(sp.diff(H2, X, 2)))
check("M1b sign convention: conformal g_mn=(1+eH)d => R = -3 e Lap H",
      Rconf_coeff == -3, f"coeff = {Rconf_coeff}")

# Gilkey connection and E for D = -(g^{mn} d_m d_n + a^s d_s + b):
# a^x = w_x', others 0; b = -Mass.
a_up = [0, sp.diff(wx, X), 0, 0]
b = -Mass

omega = [0] * 4
for d in range(4):
    s = 0
    for n in range(4):
        contr = sum(g_up[mm, ss] * Gamma[n][mm][ss]
                    for mm in range(4) for ss in range(4))
        s += g_dn[n, d] * (a_up[n] + contr) / 2
    omega[d] = lin(s)

E = b
for m in range(4):
    for n in range(4):
        term = sp.diff(omega[n], coords[m]) + omega[m] * omega[n]
        for l in range(4):
            term -= omega[l] * Gamma[l][m][n]
        E -= g_up[m, n] * term
E = lin(E)

integrand = lin(E + Rscal / 6)
int_e = sp.expand(integrand.coeff(e))
c_hpp = sp.simplify(int_e.coeff(hpp))
c_h = sp.simplify(int_e.coeff(h))
check("M1c a1-integrand h''-coefficient = -(5p+q)/12 up to R-convention sign",
      sp.simplify(c_hpp**2 - ((5 * p + q) / 12)**2) == 0,
      f"coeff = {c_hpp}")
check("M1d mass response r enters only the potential sector (h-term ~ m0^2)",
      sp.simplify(c_h + m0**2 * rr) == 0, f"h-coeff = {c_h}")

# master a1 := (h''-coefficient of integrand) / (h''-coefficient of R)
a1_master = sp.simplify(c_hpp / R_coeff)
master_expected = -(5 * p + q) / (12 * (2 * p + q))
check("M1e MASTER FORM a1(p,q) = -(5p+q)/(12(2p+q)) [replicates h2.5]",
      sp.simplify(a1_master - master_expected) == 0,
      f"a1 = {sp.simplify(a1_master)}")

# ======================================================================
# M2. Anchors.
# ======================================================================
check("M2a unimodular ray q=-3p: a1 = +1/6",
      sp.simplify(a1_master.subs(q, -3 * p) - sp.Rational(1, 6)) == 0)
E_unimod = sp.simplify(lin(E.subs([(q, -3 * p), (rr, 0), (m0, 0)])).coeff(e))
check("M2b unimodular ray: E = O(h^2) (flat measure = covariant measure)",
      E_unimod == 0, f"E|_lin = {E_unimod}")
check("M2c cone-only (p,q)=(0,1): a1 = -1/12",
      sp.simplify(a1_master.subs([(p, 0), (q, 1)])
                  + sp.Rational(1, 12)) == 0)
check("M2d conformal ray p=q: a1 = -1/6",
      sp.simplify(a1_master.subs(q, p) + sp.Rational(1, 6)) == 0)


def a1_of(pv, qv):
    pv, qv = sp.Rational(pv), sp.Rational(qv)
    return sp.Rational(-(5 * pv + qv), 12 * (2 * pv + qv))


# ======================================================================
# M3. Sector (p, q) from the covariantized action, conformal material
#     metric g3_ij = (1+h) delta_ij (weak disclination background, KBKK).
#     Factor dictionary (exact, sympy series-verified below):
#       sqrt(g3)                    -> (1+h)^{3/2}   (covariant measure)
#       each derivative-index pair g3^{ij}  -> (1+h)^{-1}
#       each field-index pair: covector g3^{ij} -> (1+h)^{-1}
#                              frame  delta_ab  -> (1+h)^{0}
#                              vector g3_{ij}   -> (1+h)^{+1}
#       group-trace kinetic (Q~^{-1} d_t Q~): intrinsic Killing form, NO
#       spatial index => only sqrt(g3).
#     Sector structures (from II.C / II.E, cf. h2.5 sector table):
#       B2 (photon doublet): kinetic rho0 u' u' [field pair];
#                            stiffness mu e_(ij) e_(ij) [deriv pair+field pair]
#       B3/B4 (Q~ sectors):  kinetic J Tr[(Q~^-1 dt Q~)^2] [group trace];
#                            stiffness (beta+gamma)/(alpha+beta) Gamma Gamma
#                            [deriv pair + soldered/algebra pair]
#       knot band-edge (bosonic (7.1) reading): inherits Q~ structure.
# ======================================================================
hs = sp.Symbol("h_s")          # scalar background amplitude (series var)
one = sp.Integer(1)
sqrtg = (1 + hs)**sp.Rational(3, 2)
dpair = (1 + hs)**(-1)         # derivative-index contraction
fpair = {"covector": (1 + hs)**(-1), "frame": one, "vector": (1 + hs)}


def exponent(expr):
    """extract a in expr = (1+h)^a via series: a = d(log expr)/dh at 0."""
    val = sp.simplify(sp.diff(sp.log(expr), hs).subs(hs, 0))
    return sp.nsimplify(val)


sectors = {}
for reading in ("covector", "frame", "vector"):
    fp = fpair[reading]
    # B2: displacement field
    q_B2 = exponent(sqrtg * fp)               # kinetic: sqrt(g) * field pair
    p_B2 = exponent(sqrtg * dpair * fp)       # stiffness: + one deriv pair
    # B3/B4/knot(bosonic): group kinetic, soldered stiffness
    q_Q = exponent(sqrtg)                     # group trace: measure only
    p_Q = exponent(sqrtg * dpair * fp)        # wryness contraction
    sectors[reading] = {
        "B2": (p_B2, q_B2, a1_of(p_B2, q_B2)),
        "B3/B4/knot_bosonic": (p_Q, q_Q, a1_of(p_Q, q_Q)),
    }

exp_expect = {
    "covector": {"B2": (-sp.Rational(1, 2), sp.Rational(1, 2)),
                 "B3/B4/knot_bosonic": (-sp.Rational(1, 2),
                                        sp.Rational(3, 2))},
    "frame": {"B2": (sp.Rational(1, 2), sp.Rational(3, 2)),
              "B3/B4/knot_bosonic": (sp.Rational(1, 2), sp.Rational(3, 2))},
    "vector": {"B2": (sp.Rational(3, 2), sp.Rational(5, 2)),
               "B3/B4/knot_bosonic": (sp.Rational(3, 2), sp.Rational(3, 2))},
}
ok = all(sectors[rd][sec][:2] == exp_expect[rd][sec]
         for rd in sectors for sec in sectors[rd])
check("M3a sector (p,q) exponents per reading [sympy series]", ok,
      str({rd: {k: (str(v[0]), str(v[1])) for k, v in sectors[rd].items()}
           for rd in sectors}))

a1_expect = {
    ("covector", "B2"): sp.Rational(-1, 3),
    ("covector", "B3/B4/knot_bosonic"): sp.Rational(1, 6),
    ("frame", "B2"): sp.Rational(-2, 15),
    ("frame", "B3/B4/knot_bosonic"): sp.Rational(-2, 15),
    ("vector", "B2"): sp.Rational(-5, 33),
    ("vector", "B3/B4/knot_bosonic"): sp.Rational(-1, 6),
}
ok = all(sectors[rd][sec][2] == a1_expect[(rd, sec)]
         for rd in sectors for sec in sectors[rd])
check("M3b a1 per (reading, sector): {-1/3,+1/6,-2/15,-2/15,-5/33,-1/6}", ok,
      str({f"{rd}:{sec}": str(sectors[rd][sec][2])
           for rd in sectors for sec in sectors[rd]}))

# cone-only reading (VI.B displayed form): (p,q) = (0,1)
a1_cone = a1_of(0, 1)
check("M3c cone-only (corpus displayed g_s): a1 = -1/12 < 0",
      a1_cone == sp.Rational(-1, 12) and a1_cone < 0)

# ======================================================================
# M4. Displacement invariant p = q - 1 under all uniform variance readings
#     AND all field density-weight redefinitions (weight w shifts both by 3w
#     since det g3 = (1+h)^3).
# ======================================================================
w = sp.Symbol("w", real=True)
ok = True
for reading in ("covector", "frame", "vector"):
    fp = fpair[reading]
    dens = (1 + hs)**(3 * w)   # (det g3)^w per field, squared in quadratic
    qk = exponent(sqrtg * fp * dens**2)
    pk = exponent(sqrtg * dpair * fp * dens**2)
    ok = ok and sp.simplify(pk - (qk - 1)) == 0
check("M4 invariant p = q - 1 for B2 (all variances, all density weights)",
      ok)

# ======================================================================
# M5. Wedge closure.
#     positivity wedge: (5p+q)(2p+q) < 0.
#     (i) p = q-1: wedge <=> q in (2/3, 5/6); reachable q in {1/2,3/2,5/2}
#         (+3w shifts move p,q together, preserving p=q-1): NONE inside
#         at natural weights (w=0).  => a1(B2) < 0 in every reading.
#     (ii) Q-sectors: q = 3/2 fixed; wedge <=> p in (-3/4, -3/10);
#         reachable p in {-1/2, 1/2, 3/2}: only p = -1/2 (unsoldered
#         covector) inside.
# ======================================================================
qq = sp.Symbol("qq", real=True)
sol = sp.solve_univariate_inequality(
    (5 * (qq - 1) + qq) * (2 * (qq - 1) + qq) < 0, qq, relational=False)
check("M5a B2 wedge (p=q-1) <=> q in (2/3, 5/6) [sympy]",
      sol == sp.Interval.open(sp.Rational(2, 3), sp.Rational(5, 6)),
      str(sol))
reachable_q = [sp.Rational(1, 2), sp.Rational(3, 2), sp.Rational(5, 2)]
check("M5b no reachable B2 ray in the wedge: a1(B2) < 0 in EVERY reading",
      all(not sol.contains(qv) for qv in reachable_q)
      and all(a1_of(qv - 1, qv) < 0 for qv in reachable_q),
      f"a1(B2) = {[str(a1_of(qv - 1, qv)) for qv in reachable_q]}")

pp = sp.Symbol("pp", real=True)
solQ = sp.solve_univariate_inequality(
    (5 * pp + sp.Rational(3, 2)) * (2 * pp + sp.Rational(3, 2)) < 0, pp,
    relational=False)
check("M5c Q-sector wedge (q=3/2) <=> p in (-3/4, -3/10) [sympy]",
      solQ == sp.Interval.open(sp.Rational(-3, 4), sp.Rational(-3, 10)),
      str(solQ))
reachable_p = [sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(3, 2)]
inwedge = [pv for pv in reachable_p if solQ.contains(pv)]
check("M5d only the unsoldered covector reading (p=-1/2) rescues Q-sectors",
      inwedge == [sp.Rational(-1, 2)], f"in-wedge: {inwedge}")

# ======================================================================
# M6. The rescue window (minimal constitutive/measure postulate).
#     Flat-measure field -> (det g3)^w-weighted field shifts (p,q) by
#     (3w, 3w).  B2 frame baseline (1/2, 3/2): wedge <=> w in (-5/18, -2/9).
#     w = -1/4 (exact-acoustic / unimodular measure) lies inside and lands
#     exactly on q = -3p with a1 = +1/6.  Natural weights 0 and -1/2 are
#     excluded.  Window width = 1/18.
# ======================================================================
ww = sp.Symbol("ww", real=True)
pB2 = sp.Rational(1, 2) + 3 * ww
qB2 = sp.Rational(3, 2) + 3 * ww
solw = sp.solve_univariate_inequality(
    (5 * pB2 + qB2) * (2 * pB2 + qB2) < 0, ww, relational=False)
check("M6a rescue window: w in (-5/18, -2/9), width 1/18",
      solw == sp.Interval.open(sp.Rational(-5, 18), sp.Rational(-2, 9)),
      str(solw))
check("M6b w = -1/4 (unimodular/exact-acoustic measure) inside; on q=-3p",
      solw.contains(sp.Rational(-1, 4))
      and sp.simplify(qB2.subs(ww, sp.Rational(-1, 4))
                      + 3 * pB2.subs(ww, sp.Rational(-1, 4))) == 0
      and a1_of(pB2.subs(ww, sp.Rational(-1, 4)),
                qB2.subs(ww, sp.Rational(-1, 4))) == sp.Rational(1, 6))
check("M6c natural weights excluded: w = 0 (corpus flat measure), w = -1/2 "
      "(half-density) both outside",
      (not solw.contains(0)) and (not solw.contains(sp.Rational(-1, 2))))
# equivalent moduli-dressing form: rho0, J ~ (det g3)^{c/3} with c = -3
# (inertia lab-anchored, not material-anchored) reaches q = -3p with c_mu = 0
c_rho = sp.Symbol("c_rho")
qd = sp.Rational(3, 2) + c_rho
pd = sp.Rational(1, 2)
sol_c = sp.solve(sp.Eq(qd, -3 * pd), c_rho)
check("M6d moduli-dressing form: rho0, J propto (det g3)^{-1} (c_rho = -3)",
      sol_c == [-3])

# ======================================================================
# Assemble sign map and verdict
# ======================================================================
n_pass = sum(1 for c in checks if c["pass"])
print(f"\ninternal checks: {n_pass}/{len(checks)} PASS")

sign_map = {
    "cone_only_VIB_displayed": {
        "B2": "-1/12", "B3/B4": "-1/12", "knots_bosonic": "-1/12",
        "consequence": "all negative => Sum w_s < 0 => induced Newton "
                       "constant inverted (contradicts VI.B G>0); c_GW "
                       "ratio formally survives (common sign)"},
    "KBKK_frame_soldered (favored by F10' objectivity)": {
        "B2": "-2/15", "B3/B4": "-2/15", "knots_bosonic": "-2/15",
        "knots_fermionic_as_printed": "-1/3 per Dirac (h2.5 K2)",
        "consequence": "all negative => G_ind < 0 (anti-gravity), "
                       "contradicts VI.B; with statistics-signed knots "
                       "(h2.5 K3): mixed signs => hull exit"},
    "KBKK_covector_unsoldered": {
        "B2": "-1/3", "B3/B4": "+1/6", "knots_bosonic": "+1/6",
        "consequence": "mixed signs => c_GW^2 can exit the hull; the "
                       "slaving bound (6.2) fails as an inequality; also "
                       "violates F10' frame-relative discipline"},
    "KBKK_vector": {
        "B2": "-5/33", "B3/B4": "-1/6", "knots_bosonic": "-1/6",
        "consequence": "all negative => G_ind < 0"},
    "strain_only_no_covariantization": {
        "consequence": "(p,q) = (0,0): W2 is exactly quadratic with constant "
                       "coefficients => NO induced R-term at all => the "
                       "Sakharov term (6.1) is empty and Theorem VI.1 has "
                       "no kinetic weights: reading self-destructs, forcing "
                       "the covariant reading"},
}

section = {
    "task": "T-H10: pin (p,q) from W2 + VI.A",
    "verdict_one_line": (
        "the corpus's own action PINS the photon-doublet sign: the "
        "displacement kinetic term forces p = q - 1 for B2, the positivity "
        "wedge needs q in (2/3, 5/6), and every textually consistent "
        "reading puts q on the half-integer lattice {1/2, 3/2, 5/2} -- so "
        "a1(B2) < 0 UNCONDITIONALLY. E.6's premise w_s > 0 for all s is "
        "therefore FALSE on the corpus's own action under every "
        "textually-permitted reading: T-H10's conditional defect resolves "
        "unconditionally AGAINST App. E.6 / Theorem VI.1 as printed."),
    "master_form": "a1(p,q) = -(5p+q)/(12(2p+q)) [independently re-derived, "
                   "M1; replicates h2.5]",
    "sector_pq": {rd: {sec: [str(v[0]), str(v[1]), str(v[2])]
                       for sec, v in sectors[rd].items()}
                  for rd in sectors},
    "B2_invariant": "p = q - 1 for all uniform variances and all density "
                    "weights (M4); wedge needs q in (2/3,5/6) (M5a): "
                    "unreachable (M5b)",
    "Q_sector": "q = 3/2 protected by the group-trace kinetic term; only "
                "the unsoldered covector reading (p = -1/2, a1 = +1/6) "
                "reaches the wedge, and it violates F10' objectivity "
                "(the corpus's own M-1 lesson) and cannot rescue B2 anyway",
    "sign_map": sign_map,
    "minimal_postulate": {
        "statement": "P-acoustic: the one-loop functional measure on a "
                     "defected background is the acoustic-covariant one -- "
                     "equivalently quantize (det g3)^{-1/8}-weighted fields "
                     "(w = -1/4), or dress the inertial moduli rho0, J "
                     "propto (det g3)^{-1}; every bosonic sector then sits "
                     "on the unimodular ray q = -3p with a1 = +1/6 > 0",
        "window": "any rescue must place the measure weight w inside "
                  "(-5/18, -2/9), a width-1/18 window; natural choices "
                  "w = 0 (the corpus's printed flat measure) and w = -1/2 "
                  "are outside",
        "still_needed_separately": "the knot-sector supertrace erratum "
                                   "(h2.5 K3): (p,q) cannot fix the "
                                   "fermionic statistics bookkeeping",
        "costs": [
            "contradicts the flat functional measure the Sec. III quantum "
            "chain (Fisher/Madelung, Born rates) is built on: the measure "
            "becomes defect-dependent, reopening the Wallstrom discharge "
            "and III.D rate derivations on curved backgrounds "
            "(derivational cost, not numerical: all campaign tiers are "
            "flat-background)",
            "in the moduli-dressing form, rho0, J propto 1/det g3 says "
            "defect-inserted material carries no extra inertia -- tension "
            "with the KBKK wedge-insertion picture (VI.A) and the F5/F8 "
            "homogeneous-lattice bookkeeping (II.C)",
            "unforced: no printed principle selects the window except the "
            "desired E.6 conclusion (tuning against VI.D's own "
            "'without fine-tuning' standard)",
            "does NOT break Tier-1 dispersion or any campaign numeric "
            "(all flat-background, blind to (p,q)); it is Unruh-ALIGNED "
            "(the exact-acoustic ray is the Unruh case) -- the honest "
            "trade is E.6 saved at the price of the Sec. III measure "
            "premise",
        ],
    },
    "checks": checks,
    "n_pass": n_pass,
    "n_total": len(checks),
}

summary = json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {}
summary.setdefault("phase", "H3.3")
summary["part_B_pq"] = section
summary["overall"] = {
    "T-H8": "printed convexity proof FALSE (strict counterexample); "
            "conclusion survives via history-path repair (Lemma II.1'); "
            "pairs-only theorem VIII'.1 unaffected under the repair",
    "T-H10": "the action pins a1(B2) < 0 in every textually-permitted "
             "reading: E.6 positivity premise unconditionally FALSE as "
             "printed; rescue requires the P-acoustic measure postulate "
             "(w in (-5/18,-2/9)) plus the knot supertrace erratum",
}
SUMMARY.write_text(json.dumps(summary, indent=1, sort_keys=False))
print(f"wrote {SUMMARY}")
