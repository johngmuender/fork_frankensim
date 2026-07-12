#!/usr/bin/env python3
"""
gstar_solve.py — Phase G3: is the halo-saturated constant G* equal to the
corpus's deep-BPS constant c0 = 128 sqrt(42)/(105 pi), or is there a pinned
gap?

ANSWER (derived here, then verified against the frozen campaign engine):

    G*  =  16 sqrt(2) / 9  =  2.514157444218836...   EXACTLY,

so G* != c0 = 2.514753626327608...; the gap is pinned in closed form:

    c0 - G* = 128 sqrt(42)/(105 pi) - 16 sqrt(2)/9 = 5.96182e-4
    c0 / G* = 24 sqrt(21) / (35 pi) = 1.000237134...

(c0 is transcendental — algebraic/pi — while G* is algebraic, so the
identity fails EXACTLY, not merely at some numerical precision.)

----------------------------------------------------------------------------
DERIVATION (all in the campaign's frozen conventions; axi_solve.py sectors)
----------------------------------------------------------------------------
The saturated closure objective at t -> 0 (axi3_solve.py, Section 4) is

    G = E6 + E0 - I/(16 pi),
    E6 = pi^3 INT b^2 dV,  E0 = (1/4pi) INT (1-q0) dV,
    I  = INT 2 (q1^2+q2^2) dV,

minimised over degree-1 fields.  Since E0 and I are pointwise functions of
the target point q and E6 is quadratic in the topological density b, G is
the sextic+potential (BPS-Skyrme-type) functional

    G = INT [ pi^3 b^2 + W(q) ] dV,
    W(q) = (1/4pi)(1-q0) - (1/8pi)(q1^2+q2^2)
         = (1/8pi) [ (1-cos psi)^2 + sin^2 psi cos^2 Theta ]  >= 0,

with target coordinates q0 = cos psi, (q1,q2,q3) = sin psi
(sin Theta cos phi, sin Theta sin phi, cos Theta).  NOTE: W depends on the
target DIRECTION through cos^2 Theta — the inertia discount -I/(16pi) is
exactly what turns the pure-potential (1/4pi)(1-q0) into this
Theta-dependent effective potential.

(1) SHARP LOWER BOUND (rigorous, all degree-1 fields, axisymmetric or not).
AM-GM pointwise:  pi^3 b^2 + W >= 2 sqrt(pi^3 W) |b|;  and since b dV is
1/(2 pi^2) times the pullback of the target volume form, for |deg| = 1

    G >= 2 sqrt(pi^3) (1/2pi^2) INT_{S^3} sqrt(W) dOmega
       = (1/sqrt2) INT INT sin^2 psi sin Theta
                    sqrt( (1-cos psi)^2 + sin^2 psi cos^2 Theta ) dpsi dTheta
       = 16 sqrt2/15  +  8 sqrt2 INT_0^{pi/2} sin^5 x cos x ln cot(x/2) dx
       = 16 sqrt2/15  +  8 sqrt2 * (4/45)          [by parts, exact]
       = 16 sqrt2 / 9.

(The bound quoted in axi3, pi/sqrt2 = 2.22144, is this bound with the
cos^2 Theta term dropped; keeping it is what sharpens 2.221 -> 2.514.)

(2) ATTAINMENT (explicit closed-form minimiser).  The bound is saturated
iff  pi^3 b^2 = W  pointwise with single-signed b.  Within the
direction-locked axisymmetric ansatz Theta = theta, F = F(r,theta) this is
the per-theta ODE  sin^2 F |dF/dr| = 2 sqrt(pi) r^2 sqrt(W(F,theta)),
which integrates in closed form.  With

    g(x) := arcsin x - x sqrt(1-x^2)      (monotone [0,1] -> [0, pi/2])

the solution is the implicit surface   g( cos(F/2) sin theta ) =
rho^3 / (6 sqrt2),  rho = r sin theta  (cylindrical radius), i.e.

    F(rho, z) = 2 arccos( min(1, r x / rho) ),   r = sqrt(rho^2+z^2),
    x = g^{-1}( min(rho^3/(6 sqrt2), pi/2) ).

Substituting back, the BPS identity reduces to the ALGEBRAIC identity
sqrt((1-cos F)^2 + sin^2 F cos^2 theta) = 2 sin(F/2) sqrt(1 - x^2)  with
x = cos(F/2) sin theta — exact.  Properties (all verified below):
  * on the symmetry axis it is EXACTLY the campaign's compacton
    F = 2 arccos(r / 2^{5/6});  support radius R(0) = 2^{5/6} = 1.78180;
  * the support is oblate: R(pi/2) = (3 sqrt2 pi)^{1/3} = 2.37098 at the
    equator (edge: g(sin theta) = rho^3/(6 sqrt2));
  * it is direction-locked (Theta = theta): the tilt freedom is NOT used
    by the true minimiser (the tilt modes in axi3's joint basis were
    compensating for the missing per-theta profile-shape freedom);
  * sector values in closed form (pushforward dV = dOmega/(2 sqrt(pi W))):
      E6 = 8 sqrt2/9,  E0 = 4 sqrt2/3,  I = 64 sqrt2 pi/9,
      E0 - I/(16 pi) = 8 sqrt2/9 = E6  (d-stationary at d = 1),
      G* = 16 sqrt2/9.
  * saturated-closure tuple in closed form:
      L = sqrt(8 pi) G* = 64 sqrt(pi)/9,   c_paper = (4/pi) G* = 64 sqrt2/(9 pi)
        = 3.20112...,   kappa_paper = 1/sqrt2,
      I_tot = 8 pi G* = 128 sqrt2 pi/9,  halo_D = 32 sqrt2 pi/9,
      E0_halo = 4 sqrt2/9,  E_stat_tot = 24 sqrt2/9 = (3/2) G*  (identity).
G is SDiff-invariant (E6 quadratic in a density, E0 and I derivative-free),
so the minimiser is a whole SDiff orbit; the field above is its locked
axisymmetric representative.  At t -> 0+ the degeneracy is lifted by
t(E2+E4) at O(t log t) in G — irrelevant to the t -> 0 value.

(3) VERIFICATION PROGRAM (this script):
  A. exact algebra gates: the bound integral by three independent routes;
     the attainer's E6/E0/I closed forms by target-space quadrature;
  B. field-level gates: Newton-inversion residual; the pointwise BPS
     identity on the quadrature nodes; analytic F_r vs finite differences;
  C. spectral spherical quadrature ladder of E6, E0, I on the explicit
     minimiser (Gauss-Legendre in mu and in the edge-regularising variable
     r = R(theta)(1-v^2)) -> closed forms to ~1e-9;
  D. FROZEN-ENGINE ladder: axi_solve.sector_pieces (REUSED, unmodified —
     the same engine that produced every reported G* in axi3/axi4) on an
     N-ladder, Richardson-extrapolated, against 16 sqrt2/9;  degree gate;
     d-stationarity gate;
  E. apples-to-apples: axi3's 18-mode extended optimum re-evaluated on the
     same engine at the same N as our minimiser;
  F. stationarity/minimality probes: contour-warp, tilt and dilation
     perturbations of the minimiser at engine level (dG = 0, d2G >= 0).

Deterministic, no RNG.  Outputs: gstar_results.json, gstar_fig.png
(gstar_run.log via tee; gstar_RESULTS.md written from the printed tables).
"""

import json
import math
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier2-closure"
sys.path.insert(0, OUTDIR)
import axi_solve as ax                    # frozen Step-2 engine (REUSED)
import axi3_solve as a3                   # Step-3 basis/objective (REUSED)

from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

SQ2 = math.sqrt(2.0)
SIXTEEN_PI = 16.0 * math.pi
GSTAR = 16.0 * SQ2 / 9.0                          # the claimed exact G*
C0 = 128.0 * math.sqrt(42.0) / (105.0 * math.pi)  # corpus deep-BPS constant
E6_EX = 8.0 * SQ2 / 9.0
E0_EX = 4.0 * SQ2 / 3.0
I_EX = 64.0 * SQ2 * math.pi / 9.0
R_AXIS = 2.0 ** (5.0 / 6.0)
R_EQ = (3.0 * SQ2 * math.pi) ** (1.0 / 3.0)
G_BOUND_OLD = math.pi / SQ2                       # axi3's (unsharpened) bound

QUICK = "--quick" in sys.argv


# ----------------------------------------------------------------------
# the closed-form minimiser
# ----------------------------------------------------------------------
_FACT = [math.factorial(2 * m + 1) for m in range(1, 14)]


def gfun(u):
    """g(sin u) = u - sin u cos u, evaluated stably (series for small u)."""
    u = np.asarray(u, float)
    out = u - np.sin(u) * np.cos(u)
    small = u < 0.5
    if np.any(small):
        us = u[small]
        u2 = us * us
        term = us.copy()
        s = np.zeros_like(us)
        for m in range(1, 14):
            term = term * u2                       # u^{2m+1}
            s = s + ((-1.0) ** (m + 1)) * (4.0 ** m) / _FACT[m - 1] * term
        out[small] = s
    return out


def u_of_y(y, iters=40):
    """Invert u - sin u cos u = y on [0, pi/2] (vectorised Newton on the
    cube-root-transformed equation, well-conditioned at both ends)."""
    y = np.asarray(y, float)
    y3 = np.cbrt(y)
    u = np.clip(np.cbrt(1.5 * y), 0.0, 0.5 * math.pi)
    for _ in range(iters):
        g = gfun(u)
        G = np.cbrt(g)
        dG = 2.0 * np.sin(u) ** 2 / (3.0 * np.maximum(g, 1e-300) ** (2.0 / 3.0))
        dG = np.where(u < 1e-12, (2.0 / 3.0) ** (1.0 / 3.0), dG)
        u = np.clip(u - (G - y3) / dG, 0.0, 0.5 * math.pi)
    return u


def x_of_rho(arho):
    """x = g^{-1}(rho^3/(6 sqrt2)), clipped at x = 1 (vacuum beyond)."""
    y = np.minimum(np.asarray(arho, float) ** 3 / (6.0 * SQ2), 0.5 * math.pi)
    return np.sin(u_of_y(y))


def F_from(r, arho):
    """The minimiser profile: F = 2 arccos(min(1, r x(rho)/rho))."""
    x = x_of_rho(arho)
    w = np.minimum(r * x / np.maximum(arho, 1e-300), 1.0)
    return 2.0 * np.arccos(w)


def qfun_min(R, Z):
    """Engine-format field builder for the exact minimiser (locked)."""
    r = np.hypot(R, Z)
    F = F_from(r, np.abs(R))
    th = np.arctan2(R, Z)
    sF = np.sin(F)
    q1 = sF * np.sin(th)
    return np.cos(F), q1, sF * np.cos(th), sF / r


def Redge(s):
    """Support radius R(theta), s = sin theta:  g(s) = (R s)^3/(6 sqrt2)."""
    s = np.asarray(s, float)
    return np.cbrt(6.0 * SQ2 * gfun(np.arcsin(s))) / s


def Fr_analytic(r, s, F, x):
    """dF/dr at fixed theta from implicit differentiation:
    F_r = - rho^2 / (sqrt2 g'(x) sin(F/2)),  g'(x) = 2 x^2/sqrt(1-x^2)."""
    rho = r * s
    gp = 2.0 * x * x / np.sqrt(np.maximum(1.0 - x * x, 1e-300))
    return -rho * rho / (SQ2 * gp * np.sin(0.5 * F))


# ----------------------------------------------------------------------
# A. exact-algebra gates
# ----------------------------------------------------------------------
def bound_quadratures():
    res = {}
    # (a) 2-D Gauss-Legendre of the target-space integral
    n = 800
    xp, wp = leggauss(n)
    psi = 0.5 * math.pi * (xp + 1.0)
    wpsi = 0.5 * math.pi * wp
    up, wu = leggauss(n)                     # u = cos Theta on [-1, 1]
    P, U = np.meshgrid(psi, up, indexing="ij")
    WW = np.outer(wpsi, wu)
    Q = (1.0 - np.cos(P)) ** 2 + np.sin(P) ** 2 * U * U
    res["G_2d"] = float(np.sum(np.sin(P) ** 2 * np.sqrt(Q) * WW)) / SQ2
    # (b) 1-D reduced form (inner Theta integral in closed form):
    # INT_{-1}^{1} sqrt(A + B u^2) du = sqrt(A+B) + (A/sqrt B) asinh(sqrt(B/A))
    c = np.cos(psi)
    inner = (2.0 * np.abs(np.sin(0.5 * psi))
             + (1.0 - c) ** 1.5 / np.sqrt(1.0 + c)
             * np.arcsinh(np.sqrt((1.0 + c) / (1.0 - c))))
    res["G_1d"] = float(np.sum(np.sin(psi) ** 2 * inner * wpsi)) / SQ2
    # (c) closed-form split 16 sqrt2/15 + 8 sqrt2 * (4/45)
    res["G_split"] = 16.0 * SQ2 / 15.0 + 8.0 * SQ2 * (4.0 / 45.0)
    # attainer sector closed forms by target-space pushforward quadrature
    sq = np.sqrt(Q)
    res["E0_tq"] = float(np.sum(np.sin(P) ** 2 * (1.0 - np.cos(P)) / sq * WW)) / SQ2
    coef = (1.0 / (2.0 * math.sqrt(math.pi))) * 2.0 * 2.0 * math.pi \
        * 2.0 * math.sqrt(2.0 * math.pi)
    res["I_tq"] = coef * float(np.sum(np.sin(P) ** 4 * (1.0 - U * U) / sq * WW))
    return res


# ----------------------------------------------------------------------
# B/C. spectral spherical quadrature of the minimiser
# ----------------------------------------------------------------------
def spectral_sectors(nmu, nv):
    """E6 (analytic F_r), E0, I of the explicit minimiser.
    Gauss-Legendre in mu = cos theta on [0,1] (x2 symmetry) and in v with
    r = R(theta)(1 - v^2) (regularises the sqrt edge).  Also returns the
    max pointwise BPS residual on the nodes."""
    xm, wm = leggauss(nmu)
    mu = 0.5 * (xm + 1.0)
    wmu = 0.5 * wm
    xv, wv = leggauss(nv)
    v = 0.5 * (xv + 1.0)
    wv = 0.5 * wv
    s = np.sqrt(1.0 - mu * mu)
    Rth = Redge(s)
    MU, V = np.meshgrid(mu, v, indexing="ij")
    S = np.sqrt(1.0 - MU * MU)
    RR = Rth[:, None] * np.ones_like(V)
    r = RR * (1.0 - V * V)
    jac = 2.0 * RR * V                          # |dr/dv|
    rho = r * S
    x = x_of_rho(rho)
    F = F_from(r, rho)
    Fr = Fr_analytic(r, S, F, x)
    W2 = np.outer(wmu, wv)
    dE6 = 0.5 * np.sin(F) ** 4 * Fr * Fr / (r * r)
    dPOT = 0.25 * r * r * ((1.0 - np.cos(F)) ** 2
                           + np.sin(F) ** 2 * MU * MU)
    E6 = 2.0 * float(np.sum(dE6 * jac * W2))
    E0 = 2.0 * float(np.sum(0.5 * (1.0 - np.cos(F)) * r * r * jac * W2))
    In = 2.0 * float(np.sum(4.0 * math.pi * np.sin(F) ** 2
                            * (1.0 - MU * MU) * r * r * jac * W2))
    # pointwise BPS residual on interior nodes (F > 0.05; at the support
    # edge both densities vanish like F^2 and the RATIO is pure roundoff)
    msk = F > 0.05
    bps = float(np.max(np.abs(dE6[msk] / dPOT[msk] - 1.0)))
    return dict(nmu=nmu, nv=nv, E6=E6, E0=E0, I=In,
                G=E6 + E0 - In / SIXTEEN_PI, bps_resid=bps)


def field_gates():
    g = {}
    # Newton residual of the inversion on a hard sample
    rho = np.concatenate([np.geomspace(1e-9, R_EQ, 4001),
                          R_EQ - np.geomspace(1e-12, 0.5, 1001)])
    y = np.minimum(rho ** 3 / (6.0 * SQ2), 0.5 * math.pi)
    u = u_of_y(y)
    g["newton_resid"] = float(np.max(np.abs(gfun(u) - y)))
    # axis limit: F(r, theta->0) -> compacton 2 arccos(r/2^{5/6})
    r = np.linspace(1e-3, R_AXIS - 1e-3, 400)
    th = 1e-7
    Fa = F_from(r, r * math.sin(th))
    Fc = 2.0 * np.arccos(np.clip(r / R_AXIS, 0.0, 1.0))
    g["axis_vs_compacton"] = float(np.max(np.abs(Fa - Fc)))
    # analytic F_r vs central finite differences (interior sample)
    thv = np.linspace(0.15, math.pi / 2 - 0.05, 9)
    errs = []
    for t in thv:
        s = math.sin(t)
        RT = float(Redge(np.array([s]))[0])
        rr = np.linspace(0.1 * RT, 0.9 * RT, 30)
        h = 1e-6
        Fp = F_from(rr + h, (rr + h) * s)
        Fm = F_from(rr - h, (rr - h) * s)
        fd = (Fp - Fm) / (2.0 * h)
        an = Fr_analytic(rr, s, F_from(rr, rr * s), x_of_rho(rr * s))
        errs.append(np.max(np.abs(fd / an - 1.0)))
    g["Fr_fd_vs_analytic"] = float(np.max(errs))
    return g


# ----------------------------------------------------------------------
# D. frozen-engine ladder + Richardson
# ----------------------------------------------------------------------
def engine_ladder(Ns, dom):
    rows = []
    for N in Ns:
        t0 = time.time()
        P = ax.sector_pieces(qfun_min, N, dom, dom)
        G = P["E6"] + P["E0"] - P["I"] / SIXTEEN_PI
        dstar = (P["E6"] / (P["E0"] - P["I"] / SIXTEEN_PI)) ** (1.0 / 6.0)
        rows.append(dict(N=N, G=G, E6=P["E6"], E0=P["E0"], I=P["I"],
                         BDEG=P["BDEG"], dstar=dstar,
                         G_err=G - GSTAR, secs=time.time() - t0))
        print("   N=%4d  G=%.9f  G-G*=%+.3e  E6-ex=%+.2e  E0-ex=%+.2e  "
              "I-ex=%+.2e  B=%.6f  d*-1=%+.1e  (%.1f s)"
              % (N, G, G - GSTAR, P["E6"] - E6_EX, P["E0"] - E0_EX,
                 P["I"] - I_EX, P["BDEG"], dstar - 1.0, rows[-1]["secs"]))
    return rows


def richardson(rows, nfit=5):
    """Fit G_N = Ginf + c N^-p on the last nfit rows; also fixed p = 2
    pairwise extrapolants.  Returns estimates + honest spread."""
    sub = rows[-nfit:]
    N = np.array([r["N"] for r in sub], float)
    G = np.array([r["G"] for r in sub], float)

    def resid(q):
        return G - (q[0] + q[1] * N ** (-q[2]))

    q0 = np.array([G[-1], (G[0] - G[-1]) * N[0] ** 2, 2.0])
    sol = least_squares(resid, q0, method="lm", xtol=1e-15, ftol=1e-15)
    Ginf, cc, p = sol.x
    fitres = float(np.max(np.abs(resid(sol.x))))
    pair = []
    for i in range(len(rows) - 1):
        N1, N2 = rows[i]["N"], rows[i + 1]["N"]
        G1, G2 = rows[i]["G"], rows[i + 1]["G"]
        pair.append(G2 + (G2 - G1) / ((N2 / N1) ** 2 - 1.0))
    est = [Ginf] + pair[-3:]
    spread = max(est) - min(est)
    budget = max(spread, fitres, abs(pair[-1] - Ginf))
    return dict(Ginf_fit=float(Ginf), p_fit=float(p), c_fit=float(cc),
                fit_resid=fitres, pairwise_p2=pair,
                estimates=[float(e) for e in est],
                spread=float(spread), budget=float(budget))


# ----------------------------------------------------------------------
# E. axi3 extended optimum on the same engine (apples-to-apples)
# ----------------------------------------------------------------------
def axi3_crosscheck(N):
    with open(OUTDIR + "/axi3_results.json") as fh:
        J = json.load(fh)
    ext = J["R1_saturated"]["extended"]
    rsup1 = ax.RSTAR * 1.02
    basis = a3.GenBasis(rsup1, K=3, a_l=(0, 2, 4), b_l=(1, 2, 3))
    qf = a3.make_qfun(ax.f_compacton, basis,
                      np.asarray(ext["a"]), np.asarray(ext["b"]))
    P = ax.sector_pieces(qf, N, rsup1 * 1.5, rsup1 * 1.5)
    d, G = a3.min_d_G(P, 0.0)
    return dict(N=N, G=float(G), d=float(d), G_json=ext["G"],
                gap_to_gstar=float(G - GSTAR))


# ----------------------------------------------------------------------
# F. stationarity / minimality probes at engine level
# ----------------------------------------------------------------------
def bump(r, c, w):
    x = np.abs((np.asarray(r, float) - c) / w)
    return np.where(x < 1.0, 1.0 - 1.5 * x * x * (1.0 - 0.5 * x),
                    0.25 * np.clip(2.0 - x, 0.0, None) ** 3)


def qfun_perturbed(mode, amp, c, w, l):
    """Perturbation families through the minimiser (amp = 0)."""
    def qf(R, Z):
        r = np.hypot(R, Z)
        aR = np.abs(R)
        s = aR / r
        cth = Z / r
        th = np.arctan2(R, Z)
        if mode == "contour":                    # F(r e^{A}, theta)
            A = amp * bump(r, c, w) * (a3.legC(l, cth))
            re = r * np.exp(A)
            F = F_from(re, re * s)
        else:
            F = F_from(r, aR)
        if mode == "tilt":
            Th = th + amp * (R / r) * cth * bump(r, c, w)
        else:
            Th = th
        sF = np.sin(F)
        q1 = sF * np.sin(Th)
        return np.cos(F), q1, sF * np.cos(Th), q1 / R
    return qf


def stationarity(N, dom, h=0.05):
    """Central-difference dG and d2G at the minimiser, at step h AND h/2:
    a true stationary point shows dG ~ h^2 (pure truncation, ratio ~ 4)."""
    modes = [("contour", 0.8, 0.8, 0), ("contour", 1.6, 0.8, 0),
             ("contour", 0.8, 0.8, 2), ("contour", 1.6, 0.8, 2),
             ("tilt", 0.8, 0.8, 1), ("tilt", 1.6, 0.8, 1)]
    P0 = ax.sector_pieces(qfun_min, N, dom, dom)
    G0 = P0["E6"] + P0["E0"] - P0["I"] / SIXTEEN_PI

    def Gat(mode, amp, c, w, l):
        P = ax.sector_pieces(qfun_perturbed(mode, amp, c, w, l), N, dom, dom)
        return P["E6"] + P["E0"] - P["I"] / SIXTEEN_PI

    rows = []
    for mode, c, w, l in modes:
        Gp, Gm = Gat(mode, h, c, w, l), Gat(mode, -h, c, w, l)
        Gp2, Gm2 = Gat(mode, 0.5 * h, c, w, l), Gat(mode, -0.5 * h, c, w, l)
        dG = (Gp - Gm) / (2.0 * h)
        dG2 = (Gp2 - Gm2) / h
        rows.append(dict(mode="%s[c=%.1f,l=%d]" % (mode, c, l), h=h,
                         dG=dG, dG_halfh=dG2,
                         trunc_ratio=dG / dG2 if dG2 else float("nan"),
                         d2G=(Gp - 2.0 * G0 + Gm) / (h * h),
                         Gp=Gp, Gm=Gm))
        print("   %-18s dG(h) = %+9.2e  dG(h/2) = %+9.2e (ratio %5.2f; "
              "4 = pure h^2 truncation)   d2G = %+9.5f   min(G+,G-)-G0 = "
              "%+9.2e"
              % (rows[-1]["mode"], dG, dG2, rows[-1]["trunc_ratio"],
                 rows[-1]["d2G"], min(Gp, Gm) - G0))
    # dilation direction (analytic in d from the sector scalings)
    dstar = (P0["E6"] / (P0["E0"] - P0["I"] / SIXTEEN_PI)) ** (1.0 / 6.0)
    print("   dilation: d* = %.7f (analytic; exact minimiser has d* = 1)"
          % dstar)
    return dict(G0=G0, N=N, rows=rows, dstar=dstar)


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
INK, GRID_C = "#333333", "#dddddd"
BLUE, AQUA, RED, VIOLET = "#2a78d6", "#1baf7a", "#e34948", "#4a3aa7"


def style_ax(a):
    a.set_facecolor("white")
    a.grid(color=GRID_C, lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.tick_params(colors=INK, labelsize=9)


def figure(rows, rich, path):
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), dpi=160)
    fig.patch.set_facecolor("white")

    a1 = axes[0]
    style_ax(a1)
    N = np.array([r["N"] for r in rows], float)
    err = np.array([abs(r["G_err"]) for r in rows])
    a1.loglog(N, err, "o-", color=BLUE, ms=7, mec="white", lw=1.6,
              label=r"$|G_N - 16\sqrt{2}/9|$ (frozen engine)")
    xg = np.array([N[0], N[-1]])
    a1.loglog(xg, err[0] * (xg / N[0]) ** -2.0, "--", color=INK, lw=1.0,
              label=r"$N^{-2}$ reference")
    a1.axhline(rich["budget"], color=RED, lw=1.0, ls=":")
    a1.text(N[0], rich["budget"] * 1.15,
            " extrapolation budget %.1e" % rich["budget"], color=RED,
            fontsize=8)
    a1.set_xlabel("engine grid N", color=INK)
    a1.set_ylabel("|G error|", color=INK)
    a1.legend(frameon=False, fontsize=8)
    a1.set_title("(a) engine convergence to $16\\sqrt{2}/9$", fontsize=10,
                 color=INK)

    a2 = axes[1]
    style_ax(a2)
    n = 460
    dom = 1.05 * R_EQ
    hx = dom / n
    xx = (np.arange(n) + 0.5) * hx
    Rg, Zg = np.meshgrid(xx, xx, indexing="ij")
    F = F_from(np.hypot(Rg, Zg), Rg)
    F0 = 2.0 * np.arccos(np.clip(np.hypot(Rg, Zg) / R_AXIS, 0.0, 1.0))
    zz = np.concatenate([-xx[::-1], xx])

    def mirror(A):
        return np.concatenate([A[:, ::-1], A], axis=1)

    lev = [np.pi / 6, np.pi / 3, np.pi / 2, 2 * np.pi / 3, 5 * np.pi / 6]
    a2.contour(xx, zz, mirror(F0).T, levels=lev, colors=RED,
               linewidths=1.0, linestyles="--")
    a2.contour(xx, zz, mirror(F).T, levels=lev, colors=BLUE, linewidths=1.6)
    th = np.linspace(1e-4, math.pi - 1e-4, 600)
    Rb = Redge(np.abs(np.sin(th)))
    a2.plot(Rb * np.sin(th), Rb * np.cos(th), color=AQUA, lw=1.8)
    a2.plot([], [], color=RED, ls="--", lw=1.0, label="spherical compacton")
    a2.plot([], [], color=BLUE, lw=1.6, label="exact $G$ minimiser")
    a2.plot([], [], color=AQUA, lw=1.8, label="support edge $R(\\theta)$")
    a2.set_aspect("equal")
    a2.legend(frameon=False, fontsize=8, loc="upper right")
    a2.set_xlabel(r"$\rho$", color=INK)
    a2.set_ylabel(r"$z$", color=INK)
    a2.set_title("(b) the oblate BPS compacton (closed form)", fontsize=10,
                 color=INK)

    a3x = axes[2]
    style_ax(a3x)
    marks = [(GSTAR, "$G^*=16\\sqrt{2}/9$ (exact)", BLUE),
             (2.5144238, "axi3 extended (18-mode family)", AQUA),
             (2.5144510, "axi3 joint (12-mode family)", AQUA),
             (2.515991, "axi4 locked ring-$G^*$ (radius-warp family)", VIOLET),
             (C0, "$\\mathfrak{c}_0 = 128\\sqrt{42}/(105\\pi)$", RED)]
    for i, (v, lab, col) in enumerate(marks):
        a3x.plot([v], [i], "o", ms=9, color=col, mec="white")
        a3x.text(v + 6e-5, i, "  %s = %.7f" % (lab, v), color=col,
                 fontsize=8, va="center")
    a3x.axvline(GSTAR, color=BLUE, lw=1.0, ls=":")
    a3x.axvline(C0, color=RED, lw=1.0, ls=":")
    a3x.set_ylim(-1.0, len(marks))
    a3x.set_xlim(2.5138, 2.5175)
    a3x.set_yticks([])
    a3x.set_xlabel("G (saturated objective, campaign convention)", color=INK)
    a3x.set_title("(c) $G^*$ vs $\\mathfrak{c}_0$: pinned gap "
                  "$5.962\\times 10^{-4}$", fontsize=10, color=INK)

    fig.suptitle("Phase G3: the saturated constant is $G^* = 16\\sqrt{2}/9$"
                 " exactly — not $\\mathfrak{c}_0$", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    out = dict(config=dict(quick=QUICK),
               exact=dict(GSTAR=GSTAR, C0=C0, gap=C0 - GSTAR,
                          gap_rel=(C0 - GSTAR) / C0,
                          ratio_c0_over_gstar=C0 / GSTAR,
                          ratio_closed_form=24.0 * math.sqrt(21.0)
                          / (35.0 * math.pi),
                          E6=E6_EX, E0=E0_EX, I=I_EX,
                          R_axis=R_AXIS, R_eq=R_EQ,
                          L=64.0 * math.sqrt(math.pi) / 9.0,
                          c_paper=64.0 * SQ2 / (9.0 * math.pi),
                          kappa_paper=1.0 / SQ2,
                          I_tot=8.0 * math.pi * GSTAR,
                          halo_D=32.0 * SQ2 * math.pi / 9.0,
                          E0_halo=4.0 * SQ2 / 9.0,
                          Estat_tot=24.0 * SQ2 / 9.0,
                          G_bound_unsharpened=G_BOUND_OLD))

    print("== EXACT CONSTANTS ==")
    print("   G*  = 16 sqrt2/9              = %.15f" % GSTAR)
    print("   c0  = 128 sqrt42/(105 pi)     = %.15f" % C0)
    print("   c0 - G*                       = %.9e   (rel %.4e)"
          % (C0 - GSTAR, (C0 - GSTAR) / C0))
    print("   c0 / G* = 24 sqrt21/(35 pi)   = %.12f  (check %.12f)"
          % (C0 / GSTAR, 24.0 * math.sqrt(21.0) / (35.0 * math.pi)))
    print("   c_paper at saturation (4/pi)G* = 64 sqrt2/(9 pi) = %.7f"
          % (64.0 * SQ2 / (9.0 * math.pi)))

    print("\n== A. EXACT-ALGEBRA GATES ==")
    bq = bound_quadratures()
    out["gates_algebra"] = {k: (v if not isinstance(v, float) else v)
                            for k, v in bq.items()}
    print("   sharp bound, 2-D target quadrature : %.15f  (err %+.1e)"
          % (bq["G_2d"], bq["G_2d"] - GSTAR))
    print("   sharp bound, 1-D reduced form      : %.15f  (err %+.1e)"
          % (bq["G_1d"], bq["G_1d"] - GSTAR))
    print("   sharp bound, closed-form split     : %.15f  (err %+.1e)"
          % (bq["G_split"], bq["G_split"] - GSTAR))
    print("   attainer E0 (pushforward quad)     : %.12f vs 4sqrt2/3 "
          "(err %+.1e)" % (bq["E0_tq"], bq["E0_tq"] - E0_EX))
    print("   attainer I  (pushforward quad)     : %.10f vs 64sqrt2pi/9 "
          "(err %+.1e)" % (bq["I_tq"], bq["I_tq"] - I_EX))

    print("\n== B. FIELD-LEVEL GATES (closed-form minimiser) ==")
    fg = field_gates()
    out["gates_field"] = fg
    print("   Newton inversion residual max      : %.2e" % fg["newton_resid"])
    print("   axis profile vs exact compacton    : %.2e"
          % fg["axis_vs_compacton"])
    print("   analytic F_r vs finite differences : %.2e (rel)"
          % fg["Fr_fd_vs_analytic"])

    print("\n== C. SPECTRAL QUADRATURE LADDER (E6 analytic-derivative) ==")
    lad = []
    for nn in ((48, 48), (96, 96), (192, 192), (320, 320)):
        s = spectral_sectors(*nn)
        lad.append(s)
        print("   n=%3dx%3d  G=%.12f (err %+.1e)  E6 %+.1e  E0 %+.1e  "
              "I %+.1e  BPS-resid %.1e"
              % (s["nmu"], s["nv"], s["G"], s["G"] - GSTAR,
                 s["E6"] - E6_EX, s["E0"] - E0_EX, s["I"] - I_EX,
                 s["bps_resid"]))
    out["spectral_ladder"] = lad

    print("\n== D. FROZEN-ENGINE LADDER (axi_solve.sector_pieces, reused) ==")
    dom = 1.06 * R_EQ
    Ns = [160, 224, 320, 448, 640, 720, 896] if not QUICK \
        else [112, 160, 224, 320]
    rows = engine_ladder(Ns, dom)
    rich = richardson(rows, nfit=5 if not QUICK else 4)
    out["engine_ladder"] = rows
    out["engine_richardson"] = rich
    print("   Richardson: Ginf(fit) = %.9f  p = %.2f  fit-resid %.1e"
          % (rich["Ginf_fit"], rich["p_fit"], rich["fit_resid"]))
    print("   pairwise p=2 extrapolants (last 3): %s"
          % ", ".join("%.9f" % e for e in rich["pairwise_p2"][-3:]))
    print("   Ginf estimates spread (budget)    : %.2e" % rich["budget"])
    print("   |Ginf - 16 sqrt2/9|               : %.2e   %s"
          % (abs(rich["Ginf_fit"] - GSTAR),
             "CONFIRMED within budget"
             if abs(rich["Ginf_fit"] - GSTAR) <= 3.0 * rich["budget"]
             else "TENSION"))

    print("\n== E. APPLES-TO-APPLES: axi3 extended optimum, same engine ==")
    xc = axi3_crosscheck(720 if not QUICK else 320)
    out["axi3_crosscheck"] = xc
    print("   axi3 extended-basis optimum: G = %.7f at N=%d "
          "(json value %.7f)" % (xc["G"], xc["N"], xc["G_json"]))
    print("   gap above exact G*: %+..3e".replace("..", ".")
          % xc["gap_to_gstar"])

    print("\n== F. STATIONARITY / MINIMALITY PROBES (engine N=320) ==")
    st = stationarity(320 if not QUICK else 160, dom)
    out["stationarity"] = st

    # adjudication summary numbers
    out["adjudication"] = dict(
        gstar_exact="16*sqrt(2)/9",
        gstar_value=GSTAR,
        c0_exact="128*sqrt(42)/(105*pi)",
        c0_value=C0,
        gap=C0 - GSTAR, gap_rel=(C0 - GSTAR) / C0,
        engine_budget=rich["budget"],
        gap_over_budget=(C0 - GSTAR) / rich["budget"],
        identity_holds=False,
        family_gaps=dict(axi3_joint=2.5144510380468388 - GSTAR,
                         axi3_extended=2.5144237706410673 - GSTAR,
                         axi4_locked_ring=2.515991 - GSTAR),
        note="G* is algebraic, c0 = algebraic/pi is transcendental: "
             "G* = c0 is impossible exactly; measured gap 5.9618e-4 "
             "pinned in closed form c0 - 16 sqrt2/9.")

    print("\n== VERDICT (analysis layer — this step does the math) ==")
    print("   G* = 16 sqrt2/9 = %.12f EXACTLY (rigorous sharp bound + "
          "explicit attainer," % GSTAR)
    print("   engine-verified to %.1e).  c0 = %.12f." % (rich["budget"], C0))
    print("   |G* - c0| = %.6e  =  %.0f x the engine budget."
          % (C0 - GSTAR, (C0 - GSTAR) / rich["budget"]))
    print("   IDENTITY G* = c0: FAILS (pinned gap, closed form "
          "c0 - 16 sqrt2/9; ratio 24 sqrt21/(35 pi) != 1).")

    out["runtime_s"] = time.time() - t_start
    with open(OUTDIR + "/gstar_results.json", "w") as fh:
        json.dump(a3.jsonable(out), fh, indent=2)
    figure(rows, rich, OUTDIR + "/gstar_fig.png")
    print("\nruntime: %.1f s" % out["runtime_s"])
    print("wrote gstar_results.json, gstar_fig.png")


if __name__ == "__main__":
    main()
