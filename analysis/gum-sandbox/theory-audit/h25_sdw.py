#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h25_sdw.py — Phase H2.5 (ROADMAP_v5_THEORY item 5): per-sector Seeley–DeWitt a1
for the App. E.6 positivity claim behind Theorem VI.1 (slaved graviton cone).

THE CLAIM UNDER AUDIT (Omega paper v2.0.1, Sec. VI.B/VI.C + App. E.6):
  Gamma_ind = Sum_s (N_s/32 pi^2) Int sqrt(-g_s) [2 Lambda_s^4 a0
              + Lambda_s^2 a1 R[g_s] + a2 ln(Lambda_s) R^2],   (6.1)
  summed over sectors s in {B2+-, B3, B4, knot band-edges},
  g_s = gbar + 2 delta_s c^2 n (x) n, and
  c_GW^2 = Sum_s w_s c_s^2 / Sum_s w_s,  w_s = N_s Lambda_s^2 a1^(s) > 0.  (6.2)
  E.6's entire justification of a1^(s) > 0 is the phrase "heat-kernel positivity",
  which is a theorem about a0, not a1 (audit T4, G2: GAP).  This run computes
  a1^(s) per sector from the corpus's own quadratic action.

CONVENTIONS (stated, fixed for the whole run):
  * Euclidean signature (+,+,+,+), d = 4, static backgrounds, coordinates
    (tau, x, y, z); one-loop determinants over R^4 with flat functional measure.
  * Laplace-type reduction (Gilkey / Vassilevich "Heat kernel expansion: user's
    manual", eqs. (2.1)-(2.15)): for D = -(g^{mu nu} d_mu d_nu + a^s d_s + b),
       omega_d = (1/2) g_{d nu} (a^nu + g^{mu s} Gamma^nu_{mu s}),
       E = b - g^{mu nu}(d_mu omega_nu + omega_mu omega_nu - omega_s Gamma^s_{mu nu}),
    so that D = -(Delta_g + E) and
       Tr exp(-tD) = (4 pi t)^{-d/2} Sum_n t^n Int sqrt(g) tr a_n,
       a1-integrand = tr(E + R/6).
  * Curvature sign convention: R(S^3, radius rho) = +6/rho^2 (verified below at
    machine-exact symbolic level; spheres have R > 0).
  * "a1^(s)" in the corpus's (6.1) sense = the coefficient of R in tr(E + R/6):
    the non-derivative part of E (the sector mass -m_s^2 and its modulation)
    renormalizes the Lambda^4 / potential sector, NOT the R term.  The
    uncharitable "literal" reading (a1 includes -m_s^2) is also reported.
  * Bundle traces: after polarization decoupling every micropolar sector is
    N_s identical scalars (tr 1 = N_s); the Dirac option for knot band-edges is
    rank 4 (tr 1 = 4) with Lichnerowicz E = -R/4 (torsion-free import).

WHAT IS COMPUTED:
  PART 1 — symbolic branch operators. The corpus's W2 (Omega II.C, Case A
    objective mass; chi3 = 0 baseline) is assembled exactly as in
    tier1-spectrum/spectrum.py (K += 2c R^dag R), and every branch dispersion is
    derived in closed form: B1 (excluded from the E.6 sum, kept for reference),
    B2+- light doublet, B3 heavy (Klein-Gordon) doublet, B4 longitudinal twist.
    Validated against the campaign's measured tier-1 numbers.
  PART 2 — Laplace-type reduction of the generic sector operator
       D = - w_t(x) d_tau^2 - d_i w_x(x) d_i + M(x),
    (w_t = sector inertia, w_x = sector stiffness, M = sector gap), which is the
    variational (divergence-form) operator every branch reduces to at second
    order in derivatives.  Background modulation ("how a defect enters the
    medium") is parametrized by response exponents (p, q, r):
       w_x = w_x0 (1 + eps p h(x)),  w_t = w_t0 (1 + eps q h(x)),
       M = M0 (1 + eps r h(x)),   h arbitrary static profile, linear order.
    The Gilkey metric is forced: g^{mu nu} = diag(w_t, w_x delta_ij); E and
    R[g] are computed symbolically and the R-coefficient
       a1_ratio(p, q) = (coefficient of Delta h in tr(E + R/6)) /
                        (coefficient of Delta h in R)
    is extracted in closed form.  Anchors: (i) S^3 curvature; (ii) minimal
    Laplacian => E = 0; (iii) conformal operator => a1 = 0; (iv) the
    unimodular ray q = -3p (= the exact-acoustic / Unruh case) => a1 = +1/6.
  PART 3 — the knot band-edge sector under the corpus's Theorem c''
    (knots are fermions): the three bookkeeping readings.
  PART 4 — deterministic numeric scans: (i) Eringen-admissible moduli scan
    (positivity of every sector cone and gap — the acoustic metrics exist);
    (ii) the (p,q) response-ray sign map; (iii) tier-1 cross-validation of the
    closed-form dispersions by direct eigensolve.

Deterministic: no randomness beyond numpy PCG64 seed 20260716.
Outputs: h25_results.json (+ stdout log).  Runtime target: << 10 min.
"""

import json
import time

import numpy as np
import sympy as sp

T0 = time.time()
SEED = 20260716
RESULTS = {"phase": "H2.5", "date": "2026-07-16", "seed": SEED}
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return bool(ok)


# ======================================================================
# PART 1 — symbolic branch operators from the corpus's own W2
# ======================================================================
print("\n=== PART 1: symbolic branch dispersions (corpus W2, Case A, chi3=0) ===")

kk = sp.symbols("k", positive=True)
lam, mu, mu_c, alp, bet, gam, Jm, rho0 = sp.symbols(
    "lambda mu mu_c alpha beta gamma J rho_0", positive=True)
mV = sp.symbols("m_V", nonnegative=True)

EPS = {}
for (i, j, l), s in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                     ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)]:
    EPS[(i, j, l)] = s


def eps3(i, j, l):
    return EPS.get((i, j, l), 0)


# operators exactly as in tier1-spectrum/spectrum.py (k || z, state (u, phi) in C^6)
Re_ = [[[sp.Integer(0)] * 6 for _ in range(3)] for _ in range(3)]
RG = [[[sp.Integer(0)] * 6 for _ in range(3)] for _ in range(3)]
for j in range(3):
    Re_[2][j][j] += sp.I * kk          # d_z u_j -> ik u_j
    RG[2][j][3 + j] += sp.I * kk       # d_z phi_j -> ik phi_j
for i in range(3):
    for j in range(3):
        for l in range(3):
            Re_[i][j][3 + l] -= eps3(i, j, l)

Rpsi = sp.zeros(3, 6)
Rpsi[0, 3] = 1
Rpsi[0, 1] = sp.I * kk / 2             # psi_x = phi_x + (ik/2) u_y
Rpsi[1, 4] = 1
Rpsi[1, 0] = -sp.I * kk / 2            # psi_y = phi_y - (ik/2) u_x
Rpsi[2, 5] = 1                         # psi_z = phi_z


def tensor_parts(R):
    tr = sp.Matrix(1, 6, lambda _, c: R[0][0][c] + R[1][1][c] + R[2][2][c])
    sym = sp.Matrix(9, 6, lambda r, c: (R[r // 3][r % 3][c] + R[r % 3][r // 3][c]) / 2)
    asym = sp.Matrix(9, 6, lambda r, c: (R[r // 3][r % 3][c] - R[r % 3][r // 3][c]) / 2)
    return tr, sym, asym


etr, esym, easym = tensor_parts(Re_)
Gtr, Gsym, Gasym = tensor_parts(RG)

K6 = (lam * etr.H * etr + 2 * mu * esym.H * esym + 2 * mu_c * easym.H * easym
      + alp * Gtr.H * Gtr + bet * Gsym.H * Gsym + gam * Gasym.H * Gasym
      + mV**2 * Rpsi.H * Rpsi)              # Case A (objective mass), chi3 = 0
K6 = sp.expand(K6)
check("K6 Hermitian (symbolic)", sp.simplify(K6 - K6.H) == sp.zeros(6, 6))

# --- exact decoupling: u_z (B1) | phi_z (B4) | transverse (B2 +- , B3)
w2_B1 = sp.simplify(K6[2, 2] / rho0)
w2_B4 = sp.simplify(K6[5, 5] / Jm)
check("B1: omega^2 = ((lambda+2mu)/rho0) k^2 exactly",
      sp.simplify(w2_B1 - (lam + 2 * mu) / rho0 * kk**2) == 0, w2_B1)
check("B4: omega^2 = ((alpha+beta)k^2 + 4mu_c + m_V^2)/J exactly",
      sp.simplify(w2_B4 - ((alp + bet) * kk**2 + 4 * mu_c + mV**2) / Jm) == 0, w2_B4)

# transverse 4x4 -> circular basis (u+, phi+) (+) (u-, phi-)
idxT = [0, 1, 3, 4]
KT = K6[idxT, idxT]
s2 = 1 / sp.sqrt(2)
C = sp.Matrix([
    [s2, 0, s2, 0],
    [sp.I * s2, 0, -sp.I * s2, 0],
    [0, s2, 0, s2],
    [0, sp.I * s2, 0, -sp.I * s2]])
Kc = sp.expand(C.H * KT * C)
offblock = sp.simplify(sp.Matrix([[Kc[0, 2], Kc[0, 3]], [Kc[1, 2], Kc[1, 3]]]))
check("circular basis block-diagonalizes transverse sector (chi3=0)",
      offblock == sp.zeros(2, 2))
K2 = Kc[0:2, 0:2]  # one circular polarization: basis (u+, phi+)
A_c = 2 * mu_c + mV**2 / 2
check("2x2 block structure: K_uu = mu k^2 + A_c k^2/2, K_pp = 2A_c + (beta+gamma)k^2/2, |K_up| = A_c k",
      sp.simplify(K2[0, 0] - (mu * kk**2 + A_c * kk**2 / 2)) == 0
      and sp.simplify(K2[1, 1] - (2 * A_c + (bet + gam) * kk**2 / 2)) == 0
      and sp.simplify(sp.Abs(K2[0, 1])**2 - A_c**2 * kk**2) == 0)

# dispersion roots of det(K2 - w2 diag(rho0, J)) = 0
w2s = sp.symbols("w2")
charpoly = sp.expand((K2[0, 0] - w2s * rho0) * (K2[1, 1] - w2s * Jm)
                     - K2[0, 1] * sp.conjugate(K2[0, 1]))
def clean_sqrt(e):
    # factor perfect squares under sqrt (sympy leaves sqrt((m_V^2+4mu_c)^2) unevaluated)
    return e.replace(lambda x: x.is_Pow and x.exp == sp.S.Half,
                     lambda x: sp.sqrt(sp.factor(x.base)))


roots = [sp.expand(clean_sqrt(r)) for r in sp.solve(charpoly, w2s)]
gap0 = [sp.simplify(clean_sqrt(sp.simplify(r.subs(kk, 0)))) for r in roots]
light, heavy = (roots[0], roots[1]) if gap0[0] == 0 else (roots[1], roots[0])
check("root selection: one gapless + one gapped root at k = 0", 0 in gap0, gap0)

light_ser = sp.simplify(clean_sqrt(sp.series(light, kk, 0, 3).removeO()))
heavy_ser = sp.expand(clean_sqrt(sp.series(heavy, kk, 0, 3).removeO()))
c2_B2 = sp.simplify(clean_sqrt(sp.simplify(light_ser / kk**2)))
gap_B3 = sp.simplify(clean_sqrt(sp.simplify(heavy_ser.subs(kk, 0))))
c2_B3 = sp.simplify(clean_sqrt(sp.simplify((heavy_ser - gap_B3) / kk**2)))

check("B2: gapless, c^2 = mu/rho0 exactly (Theorem II.2 IR content)",
      sp.simplify(c2_B2 - mu / rho0) == 0, f"c2_B2 = {c2_B2}")
check("B3: gap = (4 mu_c + m_V^2)/J (tier-1 G06 closed form)",
      sp.simplify(gap_B3 - (4 * mu_c + mV**2) / Jm) == 0, f"gap = {gap_B3}")
c2_B3_expect = (bet + gam) / (2 * Jm) + A_c / (2 * rho0)
check("B3 cone: c_psi^2 = (beta+gamma)/(2J) + (2mu_c + m_V^2/2)/(2rho0)",
      sp.simplify(c2_B3 - c2_B3_expect) == 0, f"c_psi^2 = {sp.simplify(c2_B3)}")

BENCH = {rho0: 1, Jm: 1, lam: 1, mu: 1, mu_c: 5, alp: sp.Rational(1, 2),
         bet: sp.Rational(1, 2), gam: sp.Rational(1, 2), mV: 1}
bench_vals = {
    "c_L^2": float(((lam + 2 * mu) / rho0).subs(BENCH)),
    "c_B2^2": float((mu / rho0).subs(BENCH)),
    "gap_B3=gap_B4": float(gap_B3.subs(BENCH)),
    "c_psi^2_B3": float(c2_B3.subs(BENCH)),
    "cone_B4=alpha+beta": float((alp + bet).subs(BENCH)),
}
print("benchmark values:", bench_vals)

# --- numeric cross-validation against the campaign's measured tier-1 numbers
fK = sp.lambdify((kk,), K6.subs(BENCH), "numpy")
Minv_sqrt = np.diag([1.0] * 3 + [1.0] * 3)  # rho0 = J = 1 at benchmark


def eig_at(kval):
    Km = np.array(fK(kval), dtype=complex)
    H = 0.5 * (Km + Km.conj().T)
    return np.sort(np.linalg.eigvalsh(H))


ev1 = eig_at(1.0)
# campaign table (fs-gum-cosserat RESULTS.md, eigensolve column, k = 1):
camp = {"B2": 0.947870595456, "B1": 1.732050807569, "B4": 4.690415759823,
        "B3": 5.181847289748}
w_all = np.sqrt(ev1)
ours = {
    "B2": float(np.sqrt(complex(light.subs(BENCH).subs(kk, 1)).real)),
    "B1": float(np.sqrt(bench_vals["c_L^2"])),
    "B4": float(np.sqrt(float(w2_B4.subs(BENCH).subs(kk, 1)))),
    "B3": float(np.sqrt(complex(heavy.subs(BENCH).subs(kk, 1)).real)),
}
dev_campaign = {b: abs(ours[b] - camp[b]) / camp[b] for b in camp}
check("closed-form dispersions reproduce campaign eigensolve at k=1 (rel < 1e-9)",
      max(dev_campaign.values()) < 1e-9, dev_campaign)
dev_eig = max(abs(np.sqrt(ev1[0]) - ours["B2"]), abs(np.sqrt(ev1[2]) - ours["B1"]),
              abs(np.sqrt(ev1[3]) - ours["B4"]), abs(np.sqrt(ev1[4]) - ours["B3"]))
check("closed forms match direct 6x6 eigensolve at k=1 (< 1e-10)", dev_eig < 1e-10, dev_eig)

RESULTS["part1_dispersions"] = {
    "conventions": "corpus W2 Case A (objective mass), chi3=0 baseline, k||z; "
                   "port of tier1-spectrum/spectrum.py to sympy",
    "B1": {"omega2": str(w2_B1), "in_E6_sum": False,
           "note": "longitudinal acoustic; corpus relegates it to the elliptic "
                   "constraint sector (VI B lists B2+-, B3, B4, knot band-edges)"},
    "B2": {"omega2_smallk": f"({c2_B2})*k^2 + O(k^4)", "N_s": 2, "gap": 0,
           "exact_root": str(sp.simplify(light))},
    "B3": {"omega2_smallk": f"{gap_B3} + ({sp.simplify(c2_B3)})*k^2 + O(k^4)", "N_s": 2},
    "B4": {"omega2": str(w2_B4), "N_s": 1, "exact_KG": True},
    "knot_band_edge": {"omega2": "omega_0^2 + c^2 k^2 + O(k^4 a^2)  (corpus (7.1))",
                       "statistics": "fermionic by corpus Theorem c''"},
    "benchmark_values": bench_vals,
    "campaign_cross_validation_rel_dev": dev_campaign,
}

# ======================================================================
# PART 2 — Gilkey reduction of the generic sector operator
# ======================================================================
print("\n=== PART 2: Laplace-type reduction and the a1 R-coefficient ===")

tt, xx, yy, zz = sp.symbols("tau x y z", real=True)
epsl = sp.symbols("varepsilon", real=True)
p_, q_, r_ = sp.symbols("p q r", real=True)
M0 = sp.symbols("M_0", nonnegative=True)


def christoffel(g, ginv, coords):
    n = len(coords)
    Gam = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nu in range(n):
                s = sum(ginv[l, a] * (sp.diff(g[a, nu], coords[m])
                                      + sp.diff(g[a, m], coords[nu])
                                      - sp.diff(g[m, nu], coords[a]))
                        for a in range(n))
                Gam[l][m][nu] = sp.together(s / 2)
    return Gam


def ricci_scalar(g, ginv, Gam, coords):
    n = len(coords)
    R = sp.Integer(0)
    for si in range(n):
        for nu in range(n):
            Rsn = sp.Integer(0)
            for m in range(n):
                Rsn += sp.diff(Gam[m][nu][si], coords[m]) - sp.diff(Gam[m][m][si], coords[nu])
                for l in range(n):
                    Rsn += Gam[m][m][l] * Gam[l][nu][si] - Gam[m][nu][l] * Gam[l][m][si]
            R += ginv[si, nu] * Rsn
    return R


# --- anchor 0: curvature convention on S^3 (radius rho): R must be +6/rho^2
chi_, th_, ph_, rho_ = sp.symbols("chi theta phi rho", positive=True)
g3 = sp.diag(rho_**2, rho_**2 * sp.sin(chi_)**2,
             rho_**2 * sp.sin(chi_)**2 * sp.sin(th_)**2)
g3inv = g3.inv()
Gam3 = christoffel(g3, g3inv, [chi_, th_, ph_])
R3 = sp.simplify(ricci_scalar(g3, g3inv, Gam3, [chi_, th_, ph_]))
check("curvature convention: R(S^3, rho) = +6/rho^2 (spheres positive)",
      sp.simplify(R3 - 6 / rho_**2) == 0, R3)

# --- the modulated sector background (linear order in eps; h static, arbitrary)
h = sp.Function("h")(xx, yy, zz)
coords = [tt, xx, yy, zz]
Wt = 1 + epsl * q_ * h          # sector inertia response  (w_t / w_t0)
Wx = 1 + epsl * p_ * h          # sector stiffness response (w_x / w_x0)
Mx = M0 * (1 + epsl * r_ * h)   # sector gap response
ginv = sp.diag(Wt, Wx, Wx, Wx)  # Gilkey metric is FORCED: g^{mu nu} = a^{mu nu}
g = sp.diag(1 / Wt, 1 / Wx, 1 / Wx, 1 / Wx)
Gam = christoffel(g, ginv, coords)
Rsc = ricci_scalar(g, ginv, Gam, coords)

Cvec = [sum(ginv[m, s] * Gam[nu][m][s] for m in range(4) for s in range(4))
        for nu in range(4)]


def gilkey_E(avec, bscal):
    omega = [sp.together(sum(g[d, nu] * (avec[nu] + Cvec[nu]) for nu in range(4)) / 2)
             for d in range(4)]
    E = bscal
    for m in range(4):
        for nu in range(4):
            term = (sp.diff(omega[nu], coords[m]) + omega[m] * omega[nu]
                    - sum(omega[s] * Gam[s][m][nu] for s in range(4)))
            E -= ginv[m, nu] * term
    return E


def lin(e):
    return sp.expand(sp.diff(e, epsl).subs(epsl, 0))


hxx, hyy, hzz = sp.diff(h, xx, 2), sp.diff(h, yy, 2), sp.diff(h, zz, 2)
hxy = sp.diff(h, xx, yy)

# --- anchor 1: minimal Laplacian => E = 0 identically (checked at linear order)
E_min = gilkey_E([-Cvec[nu] for nu in range(4)], sp.Integer(0))
check("anchor: minimal Laplacian -Delta_g has E = 0 (O(eps))", sp.simplify(lin(E_min)) == 0)
check("anchor: minimal Laplacian a1 = R/6, i.e. a1_ratio = +1/6 by construction", True)

# --- anchor 2: conformal operator -Delta_g + R/6 => a1 = 0
E_conf = gilkey_E([-Cvec[nu] for nu in range(4)], -Rsc / 6)
check("anchor: conformal operator (xi = 1/6) has a1-integrand = 0 (O(eps))",
      sp.simplify(lin(E_conf + Rsc / 6)) == 0)

# --- the physical (variational / divergence-form) sector operator
#     D = -Wt d_tau^2 - d_i Wx d_i + M  =  -(a^{mu nu} dd + a^s d_s + b)
avec_phys = [sp.Integer(0), sp.diff(Wx, xx), sp.diff(Wx, yy), sp.diff(Wx, zz)]
E_phys = gilkey_E(avec_phys, -Mx)
check("E_phys at eps = 0 equals -M0 (flat background)",
      sp.simplify(E_phys.subs(epsl, 0) + M0) == 0)

E_lin, R_lin = lin(E_phys), lin(Rsc)
a1_lin = sp.expand(E_lin + R_lin / 6)

cR = R_lin.coeff(hxx)
cA = sp.simplify(a1_lin.coeff(hxx))
check("isotropy: coefficients of h_xx = h_yy = h_zz in both R and a1; no mixed terms",
      sp.simplify(R_lin.coeff(hxx) - R_lin.coeff(hyy)) == 0
      and sp.simplify(a1_lin.coeff(hxx) - a1_lin.coeff(hzz)) == 0
      and R_lin.coeff(hxy) == 0 and a1_lin.coeff(hxy) == 0)
resid = sp.simplify(a1_lin - cA * (hxx + hyy + hzz) - (-M0 * r_ * h))
check("a1 linear structure: a1 = cA * Delta h - M0 r h exactly (mass part separates)",
      resid == 0, f"cA = {cA}")
check("R linear structure: R = (q + 2p) Delta h",
      sp.simplify(cR - (q_ + 2 * p_)) == 0, f"cR = {cR}")

a1_ratio = sp.simplify(cA / cR)
a1_expect = -(5 * p_ + q_) / (12 * (2 * p_ + q_))
check("MASTER RESULT: a1_ratio(p,q) = -(5p + q) / (12 (2p + q)) per scalar mode",
      sp.simplify(a1_ratio - a1_expect) == 0, f"a1_ratio = {a1_ratio}")

# --- anchor 3: unimodular ray q = -3p  <=>  sqrt(g) = const  (the exact-acoustic /
#     Unruh case, where the divergence-form operator IS the minimal Laplacian)
check("anchor: unimodular ray q = -3p gives a1_ratio = +1/6 (recovers minimal)",
      sp.simplify(a1_ratio.subs(q_, -3 * p_) - sp.Rational(1, 6)) == 0)

RAYS = {
    "cone-only (corpus's stated form g_s = gbar + 2 delta c^2 n(x)n; p=0)":
        (0, 1, sp.simplify(a1_ratio.subs({p_: 0, q_: 1}))),
    "stiffness-only (q=0)": (1, 0, sp.simplify(a1_ratio.subs({p_: 1, q_: 0}))),
    "conformal (p=q)": (1, 1, sp.simplify(a1_ratio.subs({p_: 1, q_: 1}))),
    "GR weak-field (q=-p)": (1, -1, sp.simplify(a1_ratio.subs({p_: 1, q_: -1}))),
    "unimodular / exact-acoustic (q=-3p)":
        (1, -3, sp.simplify(a1_ratio.subs({p_: 1, q_: -3}))),
}
print("named response rays, a1 per scalar mode:")
for name, (pv, qv, val) in RAYS.items():
    print(f"  {name}: a1 = {val}")

RESULTS["part2_reduction"] = {
    "operator": "D = -w_t(x) d_tau^2 - d_i w_x(x) d_i + M(x); Gilkey metric "
                "g^{mu nu} = diag(w_t, w_x, w_x, w_x) is forced by the leading symbol",
    "modulation": "w_x ~ 1 + eps p h, w_t ~ 1 + eps q h, M ~ M0(1 + eps r h), linear order",
    "R_linear": "R = (q + 2p) Delta h",
    "a1_linear": f"a1 = ({cA}) Delta h - M0 r h",
    "master_a1_ratio": str(a1_ratio),
    "positive_iff": "(5p + q)(2p + q) < 0  [opposite-sign response wedge only]",
    "pole": "q = -2p: R vanishes at linear order while the a1 term does not -- on this "
            "ray the induced term is not of the form a1*R at all (the (6.1) ansatz fails)",
    "rays": {n: {"p": v[0], "q": v[1], "a1_per_mode": str(v[2])} for n, v in RAYS.items()},
    "anchors": {"S3": str(R3), "minimal_E": "0", "conformal_a1": "0",
                "unimodular": "+1/6"},
    "mass_note": "literal reading (a1 includes -m^2): on flat/weakly-curved background "
                 "a1_literal = -m_s^2 < 0 for every gapped sector (B3, B4, knots) "
                 "whenever R < 6 m_s^2; benchmark: -21 (B3/B4).",
}

# ======================================================================
# PART 3 — the knot band-edge sector (fermionic by corpus Theorem c'')
# ======================================================================
print("\n=== PART 3: knot band-edge sector, statistics bookkeeping ===")
# Dirac-type reduction (Lichnerowicz, torsion-free import): Dslash^2 = -Delta + R/4,
# i.e. Gilkey E = -R/4 * 1_4;  tr(E + R/6) = 4(1/6 - 1/4) R = -R/3.
a1_dirac = 4 * (sp.Rational(1, 6) - sp.Rational(1, 4))
check("Dirac-type a1 = 4(1/6 - 1/4) = -1/3 per Dirac field", a1_dirac == sp.Rational(-1, 3))
# Physical loop bookkeeping: Gamma_scalar = +1/2 Tr ln D; Gamma_Dirac = -1/2 Tr ln Dslash^2
# => one Dirac field contributes like  (-1)*(-1/3)/(1/6) = +2 minimal scalar modes.
equiv_scalars = sp.simplify((-1) * a1_dirac / sp.Rational(1, 6))
check("statistics-signed Dirac contribution = +2 minimal-scalar equivalents",
      equiv_scalars == 2)

RESULTS["part3_knot_sector"] = {
    "K1_bosonic_KG_reading": "corpus (7.1) taken literally as a bosonic KG operator: "
                             "a1 follows the PART-2 master formula (minimal: +N/6; "
                             "cone-only ray: -N/12)",
    "K2_fermionic_corpus_literal": "Theorem c'' honored, w_s = N_s Lambda_s^2 a1^(s) read "
                                   "as printed with N_s > 0 a mode count: a1(Dirac) = -1/3 "
                                   "=> w_knot < 0: convexity fails as written",
    "K3_fermionic_statistics_signed": "loop sign included (Gamma_F = -Tr ln): contribution "
                                      "= +2 scalar-equivalents per Dirac field > 0 under "
                                      "the minimal reading; the printed formula then needs "
                                      "an erratum (N_s -> supertrace-signed count)",
    "import": "Lichnerowicz Dslash^2 = -Delta + R/4 (torsion-free)",
}

# ======================================================================
# PART 4 — deterministic numeric scans
# ======================================================================
print("\n=== PART 4: numeric scans ===")
rng = np.random.default_rng(SEED)
NS = 20000

# (i) Eringen-admissible moduli scan: paper's own positivity (Thm VII.4 'all moduli
# positive at the bond'; quadratic-form admissibility: mu>0, 3lam+2mu>0, mu_c>0,
# beta>0, gamma>0, 3alpha+beta>0, rho0>0, J>0, m_V^2>=0).
def lu(lo, hi, n):
    return np.exp(rng.uniform(np.log(lo), np.log(hi), n))


mu_n, muc_n, bet_n, gam_n, J_n, rho_n = (lu(1e-3, 1e3, NS) for _ in range(6))
lam_n = -2 * mu_n / 3 + lu(1e-3, 1e3, NS)          # 3 lam + 2 mu > 0, lam may be < 0
alp_n = -bet_n / 3 + lu(1e-3, 1e3, NS)             # 3 alp + bet > 0, alp may be < 0
mV2_n = np.where(rng.random(NS) < 0.1, 0.0, lu(1e-6, 1e3, NS))

cL2 = (lam_n + 2 * mu_n) / rho_n
cT2 = mu_n / rho_n
Ac_n = 2 * muc_n + mV2_n / 2
gapB3 = (4 * muc_n + mV2_n) / J_n
cpsi2 = (bet_n + gam_n) / (2 * J_n) + Ac_n / (2 * rho_n)
cB4 = (alp_n + bet_n) / J_n
adm = {
    "n_draws": NS,
    "cL2_min": float(cL2.min()), "cT2_min": float(cT2.min()),
    "cpsi2_min": float(cpsi2.min()), "cB4cone_min": float(cB4.min()),
    "gap_min": float(gapB3.min()),
}
check("admissible-moduli scan: every sector cone and gap positive (acoustic metrics exist)",
      min(cL2.min(), cT2.min(), cpsi2.min(), cB4.min()) > 0 and gapB3.min() >= 0, adm)
# and the a1_ratio is moduli-INDEPENDENT at linear order (only the response ray matters):
check("a1_ratio contains no moduli symbols (sign is moduli-independent, ray-dependent)",
      a1_ratio.free_symbols == {p_, q_})

# (ii) response-ray sign map
f_ratio = sp.lambdify((p_, q_), a1_ratio, "numpy")
thetas = np.linspace(0, 2 * np.pi, 1441)[:-1]
pv, qv = np.cos(thetas), np.sin(thetas)
den = 2 * pv + qv
vals = np.where(np.abs(den) > 1e-12, f_ratio(pv, qv), np.nan)
frac_pos = float(np.nanmean(vals > 0))
# analytic wedge measure: between q=-2p and q=-5p (both branches)
wedge = 2 * abs(np.arctan(-2.0) - np.arctan(-5.0)) / (2 * np.pi)
check("positive-sign wedge measure: numeric ray scan matches analytic wedge",
      abs(frac_pos - wedge) < 2e-3, f"numeric {frac_pos:.4f} vs analytic {wedge:.4f}")

# (iii) tier-1 dispersion cross-check on a k grid (closed forms vs eigensolve)
kgrid = np.logspace(-3, np.log10(3.0), 60)
f_light = sp.lambdify((kk,), light.subs(BENCH), "numpy")
f_heavy = sp.lambdify((kk,), heavy.subs(BENCH), "numpy")
f_B4 = sp.lambdify((kk,), w2_B4.subs(BENCH), "numpy")
maxdev = 0.0
for kv in kgrid:
    ev = eig_at(float(kv))
    pred = np.sort(np.array([
        3 * kv**2,                        # B1
        float(np.real(f_light(kv))), float(np.real(f_light(kv))),   # B2 doublet
        float(np.real(f_heavy(kv))), float(np.real(f_heavy(kv))),   # B3 doublet
        float(f_B4(kv))]))                # B4
    maxdev = max(maxdev, float(np.max(np.abs(pred - ev) / np.maximum(ev, 1e-30))))
check("closed-form spectrum == 6x6 eigensolve on full k grid (rel < 1e-8)",
      maxdev < 1e-8, f"max rel dev {maxdev:.2e}")
RESULTS["part4_scans"] = {
    "admissible_moduli_scan": adm,
    "ray_scan": {"n_rays": len(thetas), "fraction_positive": frac_pos,
                 "analytic_wedge_fraction": wedge,
                 "positive_wedge": "-5p < q < -2p (p>0) and mirror: (5p+q)(2p+q)<0"},
    "kgrid_crossvalidation_max_rel_dev": maxdev,
}

# ======================================================================
# SECTOR TABLE + VERDICT
# ======================================================================
ray_a1 = {n: v[2] for n, v in RAYS.items()}
cone_only = sp.Rational(-1, 12)


def sector_row(name, N, wt, wx, M, in_sum, extra=""):
    return {
        "sector": name, "N_s": N, "w_t (inertia)": wt, "w_x (stiffness)": wx,
        "gap M": M, "in_E6_sum": in_sum,
        "a1_minimal_reading": f"+{N}/6",
        "a1_constitutive_cone_only": str(N * cone_only),
        "a1_constitutive_general": f"{N} * (-(5p+q)/(12(2p+q)))",
        "a1_literal_flat": ("0 (massless)" if M == "0" else f"-{M} < 0"),
        "note": extra,
    }


SECTORS = [
    sector_row("B1 (long. acoustic)", 1, "rho0", "lambda + 2 mu", "0", False,
               "not in the E.6 sum (elliptic constraint sector); reference only"),
    sector_row("B2+- (light transverse doublet / photon)", 2, "rho0", "mu", "0", True,
               "c^2 = mu/rho0 exact; heavy-field elimination corrections are "
               "O(k^4/A_c), i.e. O(R/m~^2) here — subleading, truncated"),
    sector_row("B3 (heavy KG doublet)", 2, "J", "J c_psi^2 = (beta+gamma)/2 + A_c J/(2 rho0)",
               "(4 mu_c + m_V^2)/J", True,
               "gyroscopic u-phi mixing removed by branch projection, valid to O(k^2/gap)"),
    sector_row("B4 (longitudinal twist)", 1, "J", "alpha + beta", "(4 mu_c + m_V^2)/J", True,
               "exactly KG at quadratic order (no k^4 term)"),
    {
        "sector": "knot band-edges", "N_s": "per species", "w_t (inertia)": "1",
        "w_x (stiffness)": "c^2", "gap M": "omega_0^2", "in_E6_sum": True,
        "a1_bosonic_KG_reading": "master formula (minimal +N/6; cone-only -N/12)",
        "a1_fermionic_corpus_literal": "-N/3 < 0 (Lichnerowicz; w_s < 0 as printed)",
        "a1_fermionic_statistics_signed": "+2N/6 scalar-equivalent > 0 (minimal reading)",
        "note": "fermionic by corpus Theorem c''; the corpus's (7.1) writes a bosonic KG "
                "form — the tension is the audited point",
    },
]
RESULTS["sector_table"] = SECTORS

VERDICT = {
    "reduction_status": "REDUCIBLE per sector (Laplace-type after branch projection and "
                        "second-order truncation; declared error O(k^2/gap, R/gap)). The "
                        "FULL 6x6 operator is NOT Laplace-type (multi-cone leading symbol) "
                        "— the per-sector split of (6.1) is forced, and audited as such.",
    "headline": "E.6's blanket a1^(s) > 0 is NOT derivable from the corpus's operators. "
                "The sign is a function of an underdetermined constitutive response ray "
                "(p, q) that the corpus never fixes: a1 = -(5p+q)/(12(2p+q)) per mode. "
                "It is NEGATIVE on every same-sign response ray, including the corpus's "
                "own displayed perturbation form g_s = gbar + 2 delta_s c^2 n(x)n "
                "(cone-only, p = 0: a1 = -1/12 per mode), and positive only on the "
                "opposite-sign wedge -5p < q < -2p (8.5% of rays), which contains the "
                "unimodular/exact-acoustic ray q = -3p (a1 = +1/6) — a derivation the "
                "corpus would need and does not have.",
    "fermion_flank": "Under the corpus's own Theorem c'' the knot band-edge sector is "
                     "fermionic: as printed (w_s = N_s Lambda_s^2 a1^(s), N_s > 0), "
                     "a1(Dirac) = -1/3 gives w_knot < 0 in EVERY geometric reading; the "
                     "statistics-signed repair makes it +2 scalar-equivalents but "
                     "falsifies the printed formula.",
    "convexity_consequences": [
        "Reading (minimal coupling + statistics-signed fermions): all w_s > 0 — "
        "Theorem VI.1 convexity SOUND, but both ingredients are unforced repairs the "
        "corpus neither states nor derives (minimal coupling is FALSE for the medium's "
        "own variational operators except on the unimodular ray).",
        "Reading (corpus-literal formula + corpus's displayed cone-only modulation): all "
        "bosonic w_s < 0 AND w_knot < 0 — the ratio (6.2) survives formally (common sign "
        "divides out) but the induced 1/16 pi G = Sum w_s < 0: an inverted-sign Newton "
        "constant, contradicting VI.B's G ~ c^3/(hbar N Lambda^2) identification.",
        "Mixed readings (constitutive bosons + statistics-signed fermions, or minimal "
        "bosons + literal fermions): weights of BOTH signs — c_GW^2 exits the convex "
        "hull; the slaving bound |Delta c_GW/c| <= delta_UV fails as an inequality.",
    ],
    "verdict": "DEFECT-CANDIDATE (conditional) / GAP CONFIRMED AND SHARPENED: the one-line "
               "E.6 justification ('heat-kernel positivity') is false as a theorem about "
               "a1; per-sector computation shows the positivity holds only in a measure-"
               "0.085 wedge of constitutive response space (plus a statistics-sign "
               "erratum), and fails on the corpus's own displayed perturbation form. "
               "Region with negative weight: every response ray with (5p+q)(2p+q) > 0, "
               "including p=0 (cone-only), q=0 (stiffness-only), p=q (conformal), "
               "q=-p (GR weak field); plus the knot sector as printed in every reading.",
}
RESULTS["verdict"] = VERDICT
RESULTS["checks"] = CHECKS
RESULTS["all_checks_pass"] = all(c["ok"] for c in CHECKS)
RESULTS["runtime_s"] = round(time.time() - T0, 1)

out = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit/h25_results.json"
with open(out, "w") as f:
    json.dump(RESULTS, f, indent=1, default=str)
print(f"\nwrote {out}")
print(f"checks: {sum(c['ok'] for c in CHECKS)}/{len(CHECKS)} pass; "
      f"runtime {RESULTS['runtime_s']} s")
print("\nVERDICT:", VERDICT["verdict"])
