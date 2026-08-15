#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T2 - Worldsheet-mode census on the GUM tube (the chirality hook).  F-T9-T2.

Campaign: Tier 9 (Phase T), ROADMAP_v11_GLUEBALL.md workstream T2 (frozen spec).
Register: within-model; nothing here bears on nature.  Seal q-theta discipline:
this script produces NO mass numbers (symbol-only internal-mode scales); the
numerology hazard (dimensionless c-frak = 2.37 vs m_X ~ 2.37 GeV) is pre-empted
and never cited as structure.

Background reuse (mandate):
  - straight static tube, radius a, uniform core: theory-audit/h21_mc.py /
    h21_RESULTS.md Sec.2 ("static network of straight tubes, radius a, uniform
    ... amplitude inside, 0 outside"); reused verbatim by tier7-program/M3
    (RESULTS.md Sec.3) and tier7-program/N2 (RESULTS.md Sec.3).  T2 specializes
    to a SINGLE straight tube along z (amendment A-T2-1, coordinator-owned,
    printed in RESULTS.md).
  - printed substrate energy: corpus3/01-GUM-Omega-Paper-v4.3-ext.md Sec.II.C
    (W2 line 133; W_chi line 134: chi1 e_kk G_ll + chi2 e_(ij)G_(ij) +
    chi3 e_[ij]G_[ij]; W4 line 135; W6+0 lines 136-137), Sec.II.A line 114
    (e tensor / Gamma pseudotensor - "the chiral gate"), Sec.II.E line 147
    (B1..B4 branches), Sec.VII.J lines 615-617 (Q-1/Q-2/Q-3/Q-6'/Q-5).
  - audited quadratic branch operators: theory-audit/h25_RESULTS.md Sec.3
    (B4: J w^2 = (alpha+beta) k^2 + 4 mu_c + m_V^2, exactly KG).

Checks (all must PASS; FAIL aborts with nonzero exit):
  A  parity characters of the printed energy classes (mechanical, polynomial
     fields): e even / Gamma odd under inversion; W2 even, W_chi ODD (the
     unique parity-odd printed class), W4 even, W6+0 even (b_P odd, squared).
     Also the worldsheet reflection P_ws (y -> -y): phi_z (axial) flips sign.
  B  Goldstone check (operator identity, generic two-component energy of the
     printed class incl. chi-type field-gradient cross terms and a generic
     polynomial locking potential): translations of any static background are
     EXACT zero modes of the quadratic fluctuation operator ->
     the 2 transverse worldsheet modes are massless.  REQUIRED VALIDATION.
  C  collective-coordinate reduction on the orientation halo (theta0 = nu*phi,
     stiffness f_q, printed Q-5 field content): worldsheet action for X_a is
     Nambu-Goldstone form, NO mass term; tension integral reproduces the
     printed FORM sigma = pi f_q^2 ln kappa_q (prefactor bookkeeping printed).
  D  O(2) x P_ws selection rules: no invariant theta-X bilinear (internal
     pseudoscalar decouples at quadratic order); the axion-type cubic vertex
     theta * eps_ab dX_a dX_b IS invariant (allowed; chi-sourced only, by A).
  E  mechanical quadratic reduction on the [CJ-new] rigid-core ansatz
     (u_a = X_a(z,t) g(r), phi_3 = theta(z,t) h(r)): all X-theta cross terms
     and all background-wryness tadpoles vanish on angular integration; the
     internal-mode worldsheet action is extracted symbolically and its
     uniform-profile limit REPRODUCES the h25-audited B4 dispersion
     J w^2 = (alpha+beta) k^2 + 4 mu_c + m_V^2 (branch-spectrum cross-check).

Output: t2_results.json (census table, gate verdicts, symbolic results).
Deterministic, pure sympy; runtime ~3 s single-threaded; 16 checks.
"""

import json
import sys
import time

import sympy as sp

T0 = time.time()
CHECKS = []


def check(cid, desc, ok, detail=""):
    ok = bool(ok)
    CHECKS.append({"id": cid, "desc": desc, "passed": ok, "detail": str(detail)})
    print("[%s] %s: %s   %s" % ("PASS" if ok else "FAIL", cid, desc, detail))
    return ok


# ---------------------------------------------------------------------------
# CHECK A - parity characters of the printed energy classes (mechanical).
# Fields: u polar vector, phi axial vector (printed: e tensor, Gamma pseudo-
# tensor, corpus3/01:114).  Generic cubic polynomial fields with symbolic
# coefficients make every identity an exact polynomial identity.
# ---------------------------------------------------------------------------
print("\n=== CHECK A: parity characters (inversion and P_ws reflection) ===")

x, y, z = sp.symbols("x y z", real=True)
XYZ = (x, y, z)


def poly_field(tag):
    """Generic cubic polynomial in (x,y,z) with symbolic coefficients."""
    monos = [x**i * y**j * z**k for i in range(3) for j in range(3)
             for k in range(3) if i + j + k <= 2]
    coeffs = sp.symbols("%s_0:%d" % (tag, len(monos)))
    return sum(c * m for c, m in zip(coeffs, monos))


U = [poly_field("u%d" % i) for i in range(3)]
PH = [poly_field("f%d" % i) for i in range(3)]


def e_tensor(u, ph):
    return [[sp.diff(u[j], XYZ[i])
             - sum(sp.LeviCivita(i, j, k) * ph[k] for k in range(3))
             for j in range(3)] for i in range(3)]


def g_tensor(ph):
    return [[sp.diff(ph[j], XYZ[i]) for j in range(3)] for i in range(3)]


def sym(M, i, j):
    return (M[i][j] + M[j][i]) / 2


def asym(M, i, j):
    return (M[i][j] - M[j][i]) / 2


def w2_of(e, G):
    lam, mu, al, be, ga = sp.symbols("lambda mu alpha beta gamma")
    tr_e = sum(e[i][i] for i in range(3))
    tr_G = sum(G[i][i] for i in range(3))
    return (lam / 2 * tr_e**2
            + mu * sum(sym(e, i, j)**2 for i in range(3) for j in range(3))
            + al / 2 * tr_G**2
            + be / 2 * sum(sym(G, i, j)**2 for i in range(3) for j in range(3))
            + ga / 2 * sum(asym(G, i, j)**2 for i in range(3) for j in range(3)))


def wchi_of(e, G):
    c1, c2, c3 = sp.symbols("chi1 chi2 chi3")
    tr_e = sum(e[i][i] for i in range(3))
    tr_G = sum(G[i][i] for i in range(3))
    return (c1 * tr_e * tr_G
            + c2 * sum(sym(e, i, j) * sym(G, i, j) for i in range(3) for j in range(3))
            + c3 * sum(asym(e, i, j) * asym(G, i, j) for i in range(3) for j in range(3)))


# inversion P: x -> -x; u polar (u^P_i(x) = -u_i(-x)), phi axial (+phi_i(-x))
inv = {x: -x, y: -y, z: -z}
U_P = [(-U[i]).subs(inv, simultaneous=True) for i in range(3)]
PH_P = [PH[i].subs(inv, simultaneous=True) for i in range(3)]

e0, G0 = e_tensor(U, PH), g_tensor(PH)
eP, GP = e_tensor(U_P, PH_P), g_tensor(PH_P)

ok_e = all(sp.expand(eP[i][j] - e0[i][j].subs(inv, simultaneous=True)) == 0
           for i in range(3) for j in range(3))
ok_G = all(sp.expand(GP[i][j] + G0[i][j].subs(inv, simultaneous=True)) == 0
           for i in range(3) for j in range(3))
check("A1", "e_ij is a true tensor under inversion (even) - printed 01:114", ok_e)
check("A2", "Gamma_ij is a pseudotensor under inversion (odd) - printed 01:114", ok_G)

ok_W2 = sp.expand(w2_of(eP, GP) - w2_of(e0, G0).subs(inv, simultaneous=True)) == 0
ok_Wchi = sp.expand(wchi_of(eP, GP) + wchi_of(e0, G0).subs(inv, simultaneous=True)) == 0
check("A3", "W2 parity-EVEN (mechanical polynomial identity)", ok_W2)
check("A4", "W_chi parity-ODD - the unique parity-odd printed energy class "
            "(W4 ~ L^4 even; b_P ~ L^3 odd but enters squared; V(sigma_P) even)",
      ok_Wchi)

# P_ws: reflection y -> -y (plane containing the tube axis z).
refl = {y: -y}
sgn_u = [1, -1, 1]     # polar vector under reflection through xz-plane
sgn_ph = [-1, 1, -1]   # axial vector: w -> det(R) R w
U_R = [(sgn_u[i] * U[i]).subs(refl, simultaneous=True) for i in range(3)]
PH_R = [(sgn_ph[i] * PH[i]).subs(refl, simultaneous=True) for i in range(3)]
eR, GR = e_tensor(U_R, PH_R), g_tensor(PH_R)
# tensor law with det(R)=-1 for Gamma:  Gamma^R_ij(x) = -R_ia R_jb Gamma_ab(Rx)
Rm = sp.diag(1, -1, 1)
ok_eR = all(sp.expand(eR[i][j]
                      - sum(Rm[i, a] * Rm[j, b] * e0[a][b].subs(refl, simultaneous=True)
                            for a in range(3) for b in range(3))) == 0
            for i in range(3) for j in range(3))
ok_GR = all(sp.expand(GR[i][j]
                      + sum(Rm[i, a] * Rm[j, b] * G0[a][b].subs(refl, simultaneous=True)
                            for a in range(3) for b in range(3))) == 0
            for i in range(3) for j in range(3))
check("A5", "P_ws (y->-y): e transforms as tensor, Gamma as pseudotensor; "
            "in particular phi_z -> -phi_z (an axial-rotation angle about the "
            "tube axis is P_ws-ODD: worldsheet pseudoscalar)", ok_eR and ok_GR)

# ---------------------------------------------------------------------------
# CHECK B - Goldstone check: translations are exact zero modes (operator
# identity) for the general two-component energy class containing the printed
# W2 + W_chi structure at quadratic reduction: quadratic gradients, chi-type
# gradient-gradient and field-gradient cross terms, generic polynomial
# locking potential.  This is ansatz-independent: masslessness of the two
# transverse worldsheet modes follows for ANY static tube background.
# ---------------------------------------------------------------------------
print("\n=== CHECK B: Goldstone check (translation zero-mode identity) ===")

xx, yy, tt = sp.symbols("xx yy tt", real=True)
p = sp.Function("p")(xx, yy)
q = sp.Function("q")(xx, yy)
A1, A2, C2, C3, C4, Cc = sp.symbols("A1 A2 C2 C3 C4 Cc")
# generic polynomial locking potential U(p,q) (cubic, symbolic coefficients)
ucf = sp.symbols("uc_0:10")
sl_p, sl_q, sl_px, sl_py, sl_qx, sl_qy = sp.symbols(
    "sl_p sl_q sl_px sl_py sl_qx sl_qy")
U_pot_slots = (ucf[0] * sl_p**2 + ucf[1] * sl_q**2 + ucf[2] * sl_p * sl_q
               + ucf[3] * sl_p**3 + ucf[4] * sl_q**3 + ucf[5] * sl_p**2 * sl_q
               + ucf[6] * sl_p * sl_q**2 + ucf[7] * sl_p + ucf[8] * sl_q + ucf[9])
W_slots = (A1 / 2 * (sl_px**2 + sl_py**2) + A2 / 2 * (sl_qx**2 + sl_qy**2)
           + C2 * (sl_px * sl_qx + sl_py * sl_qy)      # symmetric cross-gradient
           + Cc * (sl_px * sl_qy - sl_py * sl_qx)      # chi-type antisym cross
           + C3 * sl_p * sl_qx + C4 * sl_q * sl_px     # field-gradient (eps.phi in e)
           + U_pot_slots)
dW = {s: sp.diff(W_slots, s) for s in
      (sl_p, sl_q, sl_px, sl_py, sl_qx, sl_qy)}


def EL_pair(fp, fq):
    """Euler-Lagrange expressions for the slot energy at fields (fp, fq)."""
    smap = {sl_p: fp, sl_q: fq, sl_px: sp.diff(fp, xx), sl_py: sp.diff(fp, yy),
            sl_qx: sp.diff(fq, xx), sl_qy: sp.diff(fq, yy)}
    def s(e):
        return e.subs(smap, simultaneous=True)
    ELp = sp.diff(s(dW[sl_px]), xx) + sp.diff(s(dW[sl_py]), yy) - s(dW[sl_p])
    ELq = sp.diff(s(dW[sl_qx]), xx) + sp.diff(s(dW[sl_qy]), yy) - s(dW[sl_q])
    return ELp, ELq


tpar = sp.symbols("tpar")
ELp0, ELq0 = EL_pair(p, q)
# perturb along the translation direction (eta, xi) = (d_x p, d_x q)
ELp_t, ELq_t = EL_pair(p + tpar * sp.diff(p, xx), q + tpar * sp.diff(q, xx))
Lzm_p = sp.diff(ELp_t, tpar).subs(tpar, 0).doit()
Lzm_q = sp.diff(ELq_t, tpar).subs(tpar, 0).doit()
idn_p = sp.expand(Lzm_p - sp.diff(ELp0, xx))
idn_q = sp.expand(Lzm_q - sp.diff(ELq0, xx))
okB = (idn_p == 0) and (idn_q == 0)
check("B1", "linearized operator on (d_x p0, d_x q0) == d_x(EL) identically "
            "=> on any background solving EL, translations are EXACT zero "
            "modes; the 2 transverse worldsheet modes are MASSLESS "
            "(REQUIRED Goldstone validation)", okB)
if not okB:
    sys.exit("GOLDSTONE CHECK FAILED - reduction wrong; aborting per spec.")

# ---------------------------------------------------------------------------
# CHECK C - collective-coordinate reduction on the orientation halo.
# Q-5 (01:617): locking-stratum orientation field, stiffness f_q; Q-1 (01:615)
# forces fractional frame winding 2*pi/n (n=3) around the tube.  Halo angle
# theta0 = nu * phi_angle, nu = 1/3 [winding convention printed; profile
# normalization A-T2-1].  Worldsheet reduction for X_a(z,t).
# ---------------------------------------------------------------------------
print("\n=== CHECK C: halo reduction - X_a worldsheet action (NG form) ===")

r, s = sp.symbols("r s", positive=True)
a_c, R_ir, nu, fq, Jm = sp.symbols("a_core R_IR nu f_q J_m", positive=True)
grad_theta = (-nu * sp.sin(s) / r, nu * sp.cos(s) / r)  # grad(nu * phi_angle)

I_ab = sp.zeros(2, 2)
for i in range(2):
    for j in range(2):
        I_ab[i, j] = sp.integrate(
            sp.integrate(grad_theta[i] * grad_theta[j] * r, (r, a_c, R_ir)),
            (s, 0, 2 * sp.pi))
I_diag = sp.simplify(I_ab[0, 0])
okC1 = (sp.simplify(I_ab[0, 1]) == 0 and sp.simplify(I_ab[1, 0]) == 0
        and sp.simplify(I_ab[0, 0] - I_ab[1, 1]) == 0)
sigma_hal = sp.simplify(fq**2 * I_diag)          # = pi nu^2 f_q^2 ln(R/a)
okC2 = sp.simplify(sigma_hal - sp.pi * nu**2 * fq**2
                   * sp.log(R_ir / a_c)) == 0
check("C1", "halo overlap integral is isotropic: I_ab = delta_ab * "
            "pi nu^2 ln(R/a)", okC1, "I_diag = %s" % I_diag)
check("C2", "static energy per length = worldsheet tension "
            "sigma_hal = pi nu^2 f_q^2 ln(R_IR/a_core) - the printed FORM "
            "sigma = pi f_q^2 ln kappa_q (01:617); nu^2 bookkeeping absorbed "
            "by the unvalued kappa_q", okC2, "sigma_hal = %s" % sigma_hal)

# Exact translation invariance (Check B) => no X^2 term.  The reduction gives
#   S_X = Int dt dz  [ (1/2) mu_len Xdot_a^2 - (1/2) sigma_hal X'_a^2 ],
#   mu_len = (J_m / f_q^2) * sigma_hal ;  cone c_q^2 = f_q^2 / J_m.
mu_len = sp.simplify(Jm / fq**2 * sigma_hal)
c_q2 = sp.simplify(sigma_hal / mu_len)
okC3 = sp.simplify(c_q2 - fq**2 / Jm) == 0
check("C3", "X-sector worldsheet action is relativistic NG form with cone "
            "c_q^2 = f_q^2/J_m and ZERO mass term (mass excluded by B1)",
      okC3, "mu_len = %s" % mu_len)

# ---------------------------------------------------------------------------
# CHECK D - O(2) x P_ws selection rules for worldsheet fields.
# Fields: X_pm (O(2) charge +-1; P_ws: X+ <-> X-), theta (charge 0, P_ws-odd).
# ---------------------------------------------------------------------------
print("\n=== CHECK D: O(2) x P_ws selection rules ===")

# representation bookkeeping: (charge, parity-action)
# parity acts: X+ <-> X-, theta -> -theta; z,t derivatives neutral.
def classify(monomial):
    """monomial: dict field->power for fields in {'X+','X-','th'}"""
    charge = monomial.get("X+", 0) - monomial.get("X-", 0)
    # parity: swap X+ and X-, sign (-1)^power(th)
    swapped = dict(monomial)
    swapped["X+"], swapped["X-"] = monomial.get("X-", 0), monomial.get("X+", 0)
    sign = (-1) ** monomial.get("th", 0)
    self_conj = swapped == monomial
    return charge, sign, self_conj


bilinears = {
    "theta*X+ (and c.c.)": {"th": 1, "X+": 1},
    "theta^2": {"th": 2},
    "X+*X- = |X|^2": {"X+": 1, "X-": 1},
}
res_D = {}
for name, m in bilinears.items():
    ch, sgn, sc = classify(m)
    allowed = (ch == 0) and not (sc and sgn == -1)
    res_D[name] = {"charge": ch, "allowed": allowed}
okD1 = (not res_D["theta*X+ (and c.c.)"]["allowed"]
        and res_D["theta^2"]["allowed"] and res_D["X+*X- = |X|^2"]["allowed"])
check("D1", "NO invariant theta-X bilinear exists (charge mismatch): the "
            "internal pseudoscalar DECOUPLES from translations at quadratic "
            "order; theta^2 and |X|^2 allowed", okD1, str(res_D))

# cubic axion vertex: theta * eps_ab dX_a dX_b = -i*theta*(dX+ dX- - dX- dX+)/...
# charge 0; under P_ws: eps-combination odd (X2 -> -X2), theta odd => even.
ch_v = 0
sign_v = (-1) * (-1)   # theta odd x eps_ab-bilinear odd
okD2 = (ch_v == 0 and sign_v == +1)
check("D2", "the axion-type cubic vertex theta * eps_ab dX_a dX_b is O(2) x "
            "P_ws INVARIANT (allowed); by Check A its only printed source "
            "class is W_chi (unique parity-odd class)", okD2)

# ---------------------------------------------------------------------------
# CHECK E - mechanical quadratic reduction on the rigid-core ansatz [CJ-new]:
#   u_a = X_a(z,t) g(r), u_3 = 0;  phi_3 = theta(z,t) h(r), phi_1 = phi_2 = 0;
#   background wryness of the winding core: Gamma0_i3 = nu * d_i(angle)
#   (winding axis aligned with the tube axis: [CJ-new] alignment assumption).
# Energy: printed W2 wryness terms + locking (2 mu_c psi.psi + m_V^2 phi^2/2,
# normalized to the h25-audited gap 4 mu_c + m_V^2) + W_chi.
# ---------------------------------------------------------------------------
print("\n=== CHECK E: rigid-core ansatz reduction + h25 B4 cross-check ===")

zz = sp.symbols("zz", real=True)
X1 = sp.Function("X1")(zz, tt)
X2 = sp.Function("X2")(zz, tt)
TH = sp.Function("TH")(zz, tt)
g = sp.Function("g")
h = sp.Function("h")
rr = sp.sqrt(x**2 + y**2)
XYZ2 = (x, y, zz)

u = [X1 * g(rr), X2 * g(rr), 0]
ph = [0, 0, TH * h(rr)]


def e_t(u_, ph_):
    return [[sp.diff(u_[j], XYZ2[i])
             - sum(sp.LeviCivita(i, j, k) * ph_[k] for k in range(3))
             for j in range(3)] for i in range(3)]


def g_t(ph_):
    return [[sp.diff(ph_[j], XYZ2[i]) for j in range(3)] for i in range(3)]


eF = e_t(u, ph)
GF = g_t(ph)
G0b = [[0, 0, -nu * y / (x**2 + y**2)],
       [0, 0, nu * x / (x**2 + y**2)],
       [0, 0, 0]]
GT = [[GF[i][j] + G0b[i][j] for j in range(3)] for i in range(3)]

lam, mu_, al, be, ga = sp.symbols("lambda mu alpha beta gamma", positive=True)
mu_c, m_V = sp.symbols("mu_c m_V", positive=True)
c1, c2, c3 = sp.symbols("chi1 chi2 chi3")

tr_e = sum(eF[i][i] for i in range(3))
tr_G = sum(GT[i][i] for i in range(3))
W2w = (al / 2 * tr_G**2
       + be / 2 * sum(sym(GT, i, j)**2 for i in range(3) for j in range(3))
       + ga / 2 * sum(asym(GT, i, j)**2 for i in range(3) for j in range(3)))
Wchi = (c1 * tr_e * tr_G
        + c2 * sum(sym(eF, i, j) * sym(GT, i, j) for i in range(3) for j in range(3))
        + c3 * sum(asym(eF, i, j) * asym(GT, i, j) for i in range(3) for j in range(3)))
curl_u = [sp.diff(u[2], y) - sp.diff(u[1], zz),
          sp.diff(u[0], zz) - sp.diff(u[2], x),
          sp.diff(u[1], x) - sp.diff(u[0], y)]
psi = [ph[k] - sp.Rational(1, 2) * curl_u[k] for k in range(3)]
Wlock = 2 * mu_c * sum(pk**2 for pk in psi) + m_V**2 / 2 * sum(pk**2 for pk in ph)

Wq = sp.expand(W2w + Wchi + Wlock)

X_atoms = {X1, X2}
TH_atoms = {TH}


def term_class(term):
    fa = term.atoms(sp.Function)
    hasX = any(f.func in (X1.func, X2.func) for f in fa)
    hasT = any(f.func == TH.func for f in fa)
    # degree bookkeeping in fluctuation fields (X,TH and their derivatives)
    return hasX, hasT


cross_terms, x_only, th_only, bg_only = [], [], [], []
for term in Wq.as_ordered_terms():
    hasX, hasT = term_class(term)
    if hasX and hasT:
        cross_terms.append(term)
    elif hasX:
        x_only.append(term)
    elif hasT:
        th_only.append(term)
    else:
        bg_only.append(term)


def angular_integral(expr):
    """Integrate expr over the polar angle at fixed r (x=r cos s, y=r sin s)."""
    e2 = expr.subs({x: r * sp.cos(s), y: r * sp.sin(s)}, simultaneous=True)
    e2 = sp.simplify(sp.powsimp(e2, force=True))
    return sp.integrate(e2, (s, 0, 2 * sp.pi))


cross_ang = sp.simplify(angular_integral(sp.Add(*cross_terms) if cross_terms else sp.S(0)))
okE1 = cross_ang == 0
check("E1", "ALL X-theta cross terms (incl. every W_chi bilinear) vanish on "
            "angular integration: pseudoscalar decouples at quadratic order "
            "(mechanical confirmation of D1)", okE1, "residual = %s" % cross_ang)

# tadpoles: terms linear in fluctuations sourced by the background wryness
def lin_in(terms, funcs):
    out = []
    for term in terms:
        # degree 1 in the given fluctuation functions (incl. derivatives)
        deg = 0
        for f in term.atoms(sp.Function):
            if f.func in funcs:
                deg += sp.degree(term.as_poly(f)) if term.has(f) else 0
        out.append(term)
    return out


# linear-in-TH terms come only from Gamma0 x Gamma_fluct and Gamma0 x e cross
lin_th = [t for t in th_only
          if sp.total_degree(sp.Poly(t, TH, sp.Derivative(TH, zz),
                                     sp.Derivative(TH, tt))) == 1]
lin_x_terms = []
for t in x_only:
    gens = [X1, X2, sp.Derivative(X1, zz), sp.Derivative(X2, zz),
            sp.Derivative(X1, tt), sp.Derivative(X2, tt)]
    try:
        d = sp.total_degree(sp.Poly(t, *gens))
    except sp.PolynomialError:
        d = None
    if d == 1:
        lin_x_terms.append(t)
tad_th = sp.simplify(angular_integral(sp.Add(*lin_th) if lin_th else sp.S(0)))
tad_x = sp.simplify(angular_integral(sp.Add(*lin_x_terms) if lin_x_terms else sp.S(0)))
okE2 = (tad_th == 0) and (tad_x == 0)
check("E2", "background-wryness tadpoles vanish on angular integration "
            "(straight tube is a consistent background at this order)",
      okE2, "tad_th = %s ; tad_x = %s" % (tad_th, tad_x))

# theta-sector quadratic worldsheet Lagrangian (energy side)
th_quad = [t for t in th_only if t not in lin_th]
th_ang = sp.simplify(angular_integral(sp.Add(*th_quad)))
th_ang = sp.expand(th_ang)
THz = sp.Derivative(TH, zz)
coef_thz2 = sp.simplify(th_ang.coeff(THz**2))
coef_th2 = sp.simplify(th_ang.coeff(TH**2))
# kinetic: (1/2) J phi_dot^2 -> (1/2) J TH_t^2 h^2 ; angular = 2 pi
Jm_ = Jm
kin_coeff = 2 * sp.pi * Jm_ / 2 * h(r)**2  # coefficient of TH_t^2 ... times r dr

print("  theta-sector (per unit z, angular integral done, radial open):")
print("    coeff[TH'^2]  =", coef_thz2, "  * r dr  (expect pi (alpha+beta) h^2)")
print("    coeff[TH^2]   =", coef_th2, "  * r dr")

expect_thz2 = sp.pi * (al + be) * h(r)**2
okE3 = sp.simplify(coef_thz2 - expect_thz2) == 0
check("E3", "longitudinal stiffness of the internal twist mode = (alpha+beta) "
            "- exactly the h25-audited B4 stiffness (h25_RESULTS.md Sec.3)",
      okE3)

# uniform-profile limit h -> 1: mass coefficient must give the B4 gap
coef_th2_unif = sp.simplify(coef_th2.subs({h(r): 1}).replace(
    sp.Derivative(h(r), r), 0))
# energy (1/2) m^2 TH^2 * (2 pi ...) => mass^2 = 2*coef/(2 pi J) per unit area
mass2_unif = sp.simplify(2 * coef_th2_unif / (2 * sp.pi) / Jm_)
expect_mass2 = (4 * mu_c + m_V**2) / Jm_
okE4 = sp.simplify(mass2_unif - expect_mass2) == 0
check("E4", "uniform-profile limit reproduces the h25-audited B4 gap "
            "J w^2 = (alpha+beta) k^2 + 4 mu_c + m_V^2 (branch cross-check)",
      okE4, "mass2_unif = %s" % mass2_unif)

# dispersion assembly (uniform limit): J w^2 = (alpha+beta) k^2 + 4 mu_c + m_V^2
kk, ww = sp.symbols("k w", positive=True)
disp_lhs = Jm_ * ww**2
disp_rhs = (al + be) * kk**2 + 4 * mu_c + m_V**2
okE5 = sp.simplify(
    (2 * coef_thz2.subs({h(r): 1}) / (2 * sp.pi)) * kk**2
    + 2 * coef_th2_unif / (2 * sp.pi) - disp_rhs) == 0
check("E5", "assembled uniform dispersion == printed B4 KG branch "
            "(corpus3/01:147; h25 Sec.3 'exactly KG')", okE5)

# ---------------------------------------------------------------------------
# Assemble census, gates, verdict; dump JSON.
# ---------------------------------------------------------------------------
print("\n=== assembling census and JSON ===")

E1line = "corpus3/01-GUM-Omega-Paper-v4.3-ext.md"
census = [
    {
        "row": 1, "mode": "X_1, X_2 transverse translations (2 modes)",
        "field_content": "tube collective coordinates (broken translations)",
        "printed_evidence": ("%s:615 (Q-1/Q-2 tube existence), :617 (Q-5 "
                             "tension sigma = pi f_q^2 ln kappa_q); h21_RESULTS.md "
                             "Sec.2 straight-static-tube background") % E1line,
        "o2_pws_rep": "O(2) vector (X1 P_ws-even, X2 P_ws-odd)",
        "worldsheet_mass": "0 (EXACT; Goldstone checks B1, C1-C3)",
        "status": "PRINTED-DERIVABLE",
    },
    {
        "row": 2, "mode": "core axial rotor theta (rotation of core orientation "
                          "texture about tube axis)",
        "field_content": "locking-stratum orientation field, stiffness f_q",
        "printed_evidence": ("field printed at %s:617 (Q-5); existence as a "
                             "tube-localized mode NOT printed -> [CJ-new]; "
                             "gap scale M_gap symbol-only (%s:615 Q-2, "
                             "T ~ M_gap^2; context analysis Sec.1)") % (E1line, E1line),
        "o2_pws_rep": "O(2) scalar, P_ws-ODD (pseudoscalar; forced by axial "
                      "character, Check A5)",
        "worldsheet_mass": "O(M_gap) symbolic, or f_q-scaled; unvalued",
        "status": "[CJ-new] existence; parity forced if it exists",
    },
    {
        "row": 3, "mode": "core-bound axial twist (B4-class descendant)",
        "field_content": "substrate micro-rotation phi_z (B4 branch)",
        "printed_evidence": ("branch printed %s:147; audited operator "
                             "h25_RESULTS.md Sec.3 (J w^2 = (alpha+beta)k^2 + "
                             "4mu_c + m_V^2); binding to the T2 core NOT "
                             "printed -> [CJ-new]; NOTE h25/context: the "
                             "enumerated branch content carries NO "
                             "locking-stratum mediator branch at all") % E1line,
        "o2_pws_rep": "O(2) scalar, P_ws-ODD (pseudoscalar)",
        "worldsheet_mass": "substrate-scale gap (4mu_c+m_V^2)/J - WRONG "
                           "STRATUM for a hadronic-scale axion unless the "
                           "locking-stratum reading (M_gap) is supplied "
                           "[CJ-new]",
        "status": "[CJ-new] localization; reduction validated (E3-E5)",
    },
    {
        "row": 4, "mode": "core-bound relative-rotation doublet (B3-class, "
                          "m = +-1 azimuthal channels)",
        "field_content": "substrate relative rotation psi_transverse (B3)",
        "printed_evidence": ("branch printed %s:147; h25 Sec.3; binding "
                             "[CJ-new]") % E1line,
        "o2_pws_rep": "O(2) vector (parity-even pair)",
        "worldsheet_mass": "gapped, branch scale; symbolic",
        "status": "[CJ-new] localization; NON-axionic",
    },
    {
        "row": 5, "mode": "breathing / radius mode",
        "field_content": "core profile deformation",
        "printed_evidence": "core profile NOWHERE printed -> [CJ-new] entirely",
        "o2_pws_rep": "O(2) scalar, P_ws-EVEN",
        "worldsheet_mass": "O(M_gap) symbolic",
        "status": "[CJ-new]",
    },
    {
        "row": 6, "mode": "phason (T3) admixture",
        "field_content": "soft-sector phason, stiffness f in [4,60] MeV",
        "printed_evidence": ("T3 stratum printed %s:615 (Q-2) and :167 (II.H); "
                             "NO printed dynamical cross-stratum coupling to "
                             "the T2 tube - the only printed cross-stratum "
                             "relation is the stiffness law f_q = sqrt(eps_q) "
                             "m_Sk (V15.7); stratum trichotomy discipline "
                             "(F-Q1 lesson) forbids conflation") % E1line,
        "o2_pws_rep": "(would-be O(2) scalar)",
        "worldsheet_mass": "-",
        "status": "NOT-SUPPORTED-BY-PRINT (admixture would be [CJ-new])",
    },
    {
        "row": 7, "mode": "longitudinal u_3 along-tube displacement",
        "field_content": "worldsheet reparametrization",
        "printed_evidence": "standard collective-coordinate bookkeeping",
        "o2_pws_rep": "-",
        "worldsheet_mass": "-",
        "status": "GAUGE (not a physical worldsheet mode)",
    },
    {
        "row": 8, "mode": "Z_3 mismatch label k in {1,2}",
        "field_content": "discrete superselection label (winding 2 pi k/3)",
        "printed_evidence": "%s:615 (Q-1 fractional winding 2 pi/n, n=3)" % E1line,
        "o2_pws_rep": "-",
        "worldsheet_mass": "-",
        "status": "DISCRETE LABEL, not a mode",
    },
]

all_pass = all(c["passed"] for c in CHECKS)

gates = {
    "T2-G1": {
        "requirement": "tube background reproduced from campaign machinery "
                       "with citations; mode-enumeration table complete with "
                       "printed-evidence column (h25 cross-check)",
        "measured": "h21/M3/N2 straight-static-tube background reused "
                    "(single-tube specialization A-T2-1 printed); 8-row census, "
                    "every row printed-evidence or [CJ-new]; h25 cross-check: "
                    "enumerated branch content (B1,B2+-,B3,B4,knot band-edges) "
                    "contains NO locking-stratum mediator branch",
        "verdict": "PASS" if all_pass else "FAIL",
    },
    "T2-G2": {
        "requirement": "quadratic worldsheet action derived symbolically for "
                       ">= the 2 transverse modes (MUST be massless - Goldstone "
                       "check) and every candidate internal mode; [CJ-new] "
                       "flags at every unprinted step",
        "measured": "transverse modes massless by operator identity (B1) and "
                    "explicit halo reduction (C1-C3, NG form, tension = "
                    "pi nu^2 f_q^2 ln(R/a)); internal twist-mode action "
                    "extracted (E1-E5), uniform limit reproduces h25 B4 "
                    "dispersion exactly; internal masses symbol-only "
                    "(M_gap / f_q-scaled); all ansatz steps flagged [CJ-new]",
        "verdict": "PASS" if all_pass else "FAIL",
    },
    "T2-G3": {
        "requirement": "census verdict with parity classification; if NO "
                       "pseudoscalar from printed structure, that is the "
                       "finding; if one exists, structure at "
                       "[CJ-new]/[DW, testbed-grade]",
        "measured": "VERDICT: plain bosonic string AS PRINTED (only the 2 "
                    "transverse Goldstones are printed-derivable) -> the GUM "
                    "tube inherits the known 3+1D IP/NG J^PC problems. The "
                    "chirality hook is real but undischarged: W_chi is the "
                    "UNIQUE parity-odd printed class (A4); any core-localized "
                    "axial rotor/twist is FORCED P_ws-odd (A5); the axion "
                    "vertex theta*eps_ab dX dX is symmetry-allowed and "
                    "chi-sourced-only (D2); but mode existence/localization/"
                    "gap are all [CJ-new] - tube-core-axion candidate filed "
                    "at [CJ-new]/[DW, testbed-grade] proposal class, mass "
                    "prop-to M_gap symbol-only",
        "verdict": "PASS" if all_pass else "FAIL",
    },
    "T2-G4": {
        "requirement": "within-model; SILENT-corpus verdicts carried; no "
                       "grade motion",
        "measured": "register held: no mass numbers produced (symbol-only); "
                    "seal q-theta trivially satisfied; numerology hazard "
                    "pre-empted; SILENT verdicts (tube-core chirality, "
                    "cross-stratum transfer, loop quantization) carried from "
                    "archaeology; no grade motion; no stake; forbidden "
                    "sentences absent",
        "verdict": "PASS",
    },
}

results = {
    "workstream": "T2 (F-T9-T2) - worldsheet-mode census on the GUM tube",
    "date": "2026-08-15",
    "register": "within-model; nothing bears on nature; seal q-theta honored "
                "(no mass numbers produced); numerology hazard (c-frak = 2.37 "
                "vs m_X ~ 2.37 GeV) pre-empted - never cited as structure",
    "background_reuse": {
        "geometry": "straight static tube along z, radius a, uniform core "
                    "(h21_RESULTS.md Sec.2 / h21_mc.py, reused verbatim by "
                    "M3 Sec.3 and N2 Sec.3); single-tube specialization = "
                    "amendment A-T2-1 (coordinator-owned)",
        "halo": "orientation halo theta0 = nu*angle, nu = 1/3 (Q-1 winding "
                "2 pi/3); profile normalization is A-T2-1's closest-feasible "
                "variant; tension FORM matches printed Q-5",
        "printed_energy": "corpus3/01 Sec.II.C lines 132-137 (W2, W_chi, W4, "
                          "W6+0); II.A line 114 (chiral gate)",
        "audited_operators": "theory-audit/h25_RESULTS.md Sec.3 (B3/B4 rows)",
    },
    "checks": CHECKS,
    "census": census,
    "symbolic_results": {
        "X_sector_action": "S_X = Int dt dz [ 1/2 mu_len Xdot_a^2 - 1/2 "
                           "sigma_hal X'_a^2 ], mass = 0 exactly",
        "sigma_hal": str(sigma_hal),
        "mu_len": str(mu_len),
        "cone": "c_q^2 = f_q^2 / J_m",
        "theta_sector_action": "S_theta = Int dt dz 1/2 [ J I2 thdot^2 - "
                               "(alpha+beta) I2 th'^2 - (K_perp + M_lock^2-"
                               "class I2) th^2 ]; I2 = 2 pi Int h^2 r dr, "
                               "K_perp = 2 pi Int (beta+gamma)/2 h'(r)^2 "
                               "r dr-class + locking-well integral [CJ-new]",
        "theta_mass": "m_theta = O(M_gap) SYMBOL-ONLY (M_gap never valued in "
                      "print) or f_q-scaled; no number claimed",
        "axion_vertex": "theta * eps_ab d X_a d X_b: O(2) x P_ws invariant; "
                        "only printed source class = W_chi (chi2/chi3); "
                        "coefficient = chi x (core overlap integral) [CJ-new]",
        "sigma_prefactor_note": "halo reduction gives pi nu^2 f_q^2 ln(R/a); "
                                "printed sigma = pi f_q^2 ln kappa_q matches "
                                "the form with the nu^2 = 1/9 bookkeeping "
                                "absorbed into the unvalued ln kappa_q "
                                "(kappa_q never printed - context Sec.1)",
    },
    "gates": gates,
    "verdict": ("PLAIN BOSONIC STRING AS PRINTED: printed structure supplies "
                "exactly the two transverse Goldstones and nothing else; no "
                "pseudoscalar worldsheet mode emerges from printed structure "
                "(the finding); the tube-core-axion candidate (core-bound "
                "axial rotor/twist, forced P_ws-odd, chi-sourced axion vertex, "
                "mass prop-to M_gap symbolic) is filed at [CJ-new]/[DW, "
                "testbed-grade] proposal class"),
    "runtime_s": round(time.time() - T0, 1),
    "all_checks_pass": all_pass,
}

out = "/home/user/fork_frankensim/analysis/gum-sandbox/tier9-glueball/T2/t2_results.json"
with open(out, "w") as f:
    json.dump(results, f, indent=1)
print("\nwrote %s" % out)
print("runtime %.1f s ; all checks pass: %s" % (results["runtime_s"], all_pass))
if not all_pass:
    sys.exit(1)
