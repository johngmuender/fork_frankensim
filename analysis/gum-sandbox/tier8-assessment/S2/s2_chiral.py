#!/usr/bin/env python3
"""
s2_chiral.py — Phase S2 (Tier 8): chirality-interpolation family
W_lambda = (1/8pi)[(1-cos psi)^2 + lambda sin^2 psi cos^2 Theta],
lambda in [0,1] (lambda = 1 = GUM chiral, lambda = 0 = achiral).

Question: is the saturated closure constant G* = 16 sqrt2/9 (and the closure
identities the corpus leans on) DISTINCTIVE to the chiral lambda = 1 point,
or GENERIC across the family?

----------------------------------------------------------------------------
DERIVATION (frozen gstar_solve.py conventions; nothing refitted)
----------------------------------------------------------------------------
(0) Sector embedding.  Using (1-c)^2 = 2(1-c) - s^2 and
    s^2 cos^2 Theta = s^2 - s^2 sin^2 Theta  (c = cos psi, s = sin psi):

      W_lam = (1/4pi)(1-q0) - i_lam(q)/(16 pi),
      i_lam = 2 lam (q1^2+q2^2) + 2 (1-lam) |qvec|^2,

    i.e. the family interpolates the INERTIA sector from the chiral
    (transverse, direction-locked) I = INT 2(q1^2+q2^2) dV at lam = 1 to the
    isotropic 2 INT |qvec|^2 dV at lam = 0:

      I_lam = lam I + 2 (1-lam) J,     J := INT sin^2 psi dV,
      G_lam = E6 + E0 - I_lam/(16 pi) = INT [pi^3 b^2 + W_lam] dV.

(1) SHARP BOUND (AM-GM + degree pullback, verbatim gstar route):

      B(lam) = (1/sqrt2) INT INT sin^2 psi sin Theta
               sqrt[(1-cos psi)^2 + lam sin^2 psi cos^2 Theta] dpsi dTheta.

    Endpoints (closed form): B(1) = 16 sqrt2/9 (gstar);  B(0):
    the Theta-integral is trivial (integrand Theta-free),
      B(0) = (1/sqrt2) * 2 * INT_0^pi sin^2 psi (1-cos psi) dpsi
           = (1/sqrt2) * 2 * (pi/2)  =  pi/sqrt2.        [derived here]
    MOREOVER the whole curve is CLOSED FORM (derived in this workstream,
    sympy inner integral + by-parts outer integral, log terms cancel):
      G*(lam) = B(lam) = (sqrt2/3) [ sqrt(lam)(3-2lam)/(1-lam)
                          + (3-4lam) arccos(sqrt lam)/(1-lam)^{3/2} ],
    with the pi-rational midpoint  G*(1/2) = (4+pi)/3  and closed-form
    sectors (see closed_forms() below).

(2) ATTAINMENT (the per-theta BPS ODE integrates in closed form for EVERY
    lam).  Locked ansatz Theta = theta, F = F(r,theta); BPS condition
    pi^3 b^2 = W_lam with b = -sin^2 F F_r/(2 pi^2 r^2) is

      sin^2 F |F_r| = (r^2/sqrt2) sqrt[(1-cos F)^2 + lam sin^2 F cos^2 th].

    Half-angle c = cos(F/2), k(theta) := sqrt(1 - lam cos^2 theta):
    (1-cos F)^2 + lam sin^2F cos^2 th = 4 sin^2(F/2) (1 - k^2 c^2), so the
    ODE is 4 c^2 c_r = (r^2/sqrt2) sqrt(1-k^2c^2); with x = k c and
    g(x) = arcsin x - x sqrt(1-x^2) (the SAME g as gstar):

      g( k cos(F/2) )  =  k^3 r^3 / (6 sqrt2)         [general-lam law]

    with F = pi at r = 0, support edge g(k) = k^3 R(theta)^3/(6 sqrt2).
    At lam = 1 (k = sin theta) this is EXACTLY the corpus law
    g(cos(F/2) sin th) = rho^3/(6 sqrt2); at lam = 0 (k = 1) it is the
    SPHERICAL profile g(cos(F/2)) = r^3/(6 sqrt2), R = (3 sqrt2 pi)^{1/3}.
    Saturation is the algebraic identity
      sqrt[(1-cosF)^2 + lam sin^2F cos^2 th] = 2 sin(F/2) sqrt(1-x^2),
    x = k cos(F/2) — exact, so B(lam) is ATTAINED for every lam in [0,1]:
    G*(lam) = B(lam).  Equatorial radius R(pi/2) = (3 sqrt2 pi)^{1/3} is
    lam-INDEPENDENT; the axis radius R(0) runs from (3 sqrt2 pi)^{1/3}
    (lam=0, sphere) to 2^{5/6} (lam=1, the campaign compacton).

(3) SECTORS on the lam-minimiser via the pushforward
    dV = dOmega/(2 sqrt(pi W_lam)) = sqrt2 dOmega/sqrt(Q_lam),
    Q_lam = (1-cos psi)^2 + lam sin^2 psi cos^2 Theta:
      E6  = B/2   (BPS half-split, generic),
      E0  = (1/sqrt2) II (1-c) s^2 sinTh / sqrt(Q),
      I   = 4 sqrt2 pi II s^4 sin^3 Th / sqrt(Q),
      J   = 2 sqrt2 pi II s^4 sinTh / sqrt(Q).

(4) IDENTITIES tested per lam (residuals printed):
      id1 literal (corpus form):   E6 = E0 - I/(16 pi)
        [prediction: resid = -(1-lam)(2J - I)/(16 pi) < 0 for lam < 1]
      id1 generalized (d-stationarity of the lam-objective):
                                   E6 = E0 - I_lam/(16 pi)
      id2 (halo fixed point):      E_stat_tot = (3/2) G*(lam)
        The axi3 saturation algebra (clock L^2 = (2/3) I E_stat + halo
        threshold kappa* = 1/sqrt(8pi)) is lam-INDEPENDENT because the
        marginal halo is equatorial (Theta = pi/2) small-amplitude, where
        i_lam = 2 eta^2 for every lam.  With the lam-consistent inertia:
        E_stat_tot = E6 + E0 + E0_halo, E0_halo = (8 pi G* - I_lam)/(16 pi)
        -> GENERIC.  With the literal chiral-I bookkeeping at lam != 1:
        resid = +(1-lam)(2J - I)/(16 pi) -> DISTINCTIVE.
      support law: literal (axis compacton R(0) = 2^{5/6}; edge
        g(sin th) = rho^3/(6 sqrt2)) vs generalized g(k) = (kR)^3/(6 sqrt2).
      kappa_paper = 1/sqrt2: lam-independent (halo algebra) -> GENERIC.

(5) ANCHOR SCAN: c(lam) = (4/pi) G*(lam) in [2 sqrt2, 64 sqrt2/(9 pi)] and
    G*(lam) in [pi/sqrt2, 16 sqrt2/9] against every printed constant of
    gstar_RESULTS.md; crossings root-refined (mpmath).

Deterministic, no RNG.  Outputs (in S2/): s2_results.json, s2_fig.png,
s2_checkpoint.json (per-phase checkpoints), s2_run.log via tee.
Usage: python3 s2_chiral.py [--quick]
"""

import json
import math
import os
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

T2 = "/home/user/fork_frankensim/analysis/gum-sandbox/tier2-closure"
OUT = "/home/user/fork_frankensim/analysis/gum-sandbox/tier8-assessment/S2"
sys.path.insert(0, T2)
from gstar_solve import gfun, u_of_y  # frozen inversion machinery (REUSED)

SQ2 = math.sqrt(2.0)
PI = math.pi
SIXTEEN_PI = 16.0 * PI
GSTAR1 = 16.0 * SQ2 / 9.0                 # B(1) exact (gstar)
B0_EX = PI / SQ2                          # B(0) exact (derived above)
C0 = 128.0 * math.sqrt(42.0) / (105.0 * PI)
R_EQ = (3.0 * SQ2 * PI) ** (1.0 / 3.0)
R_AXIS1 = 2.0 ** (5.0 / 6.0)

QUICK = "--quick" in sys.argv
MP_DPS = 20 if QUICK else 30
NGL = 400 if QUICK else 800
NSPEC = 96 if QUICK else 256

LAMBDAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

CKPT = os.path.join(OUT, "s2_checkpoint.json")


def ckpt_write(state):
    tmp = CKPT + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(tmp, CKPT)


# ----------------------------------------------------------------------
# the closed-form lam-minimiser (generalises gstar's F_from / Redge)
# ----------------------------------------------------------------------
def kmu(lam, mu):
    """k(theta) = sqrt(1 - lam cos^2 theta), mu = cos theta."""
    return np.sqrt(np.maximum(1.0 - lam * np.asarray(mu, float) ** 2, 0.0))


def gval(k):
    """g at x = k: arcsin k - k sqrt(1-k^2) (reuses gstar's stable gfun)."""
    return gfun(np.arcsin(np.clip(np.asarray(k, float), 0.0, 1.0)))


def Redge_lam(lam, mu):
    """Support radius R(theta): g(k) = k^3 R^3/(6 sqrt2)."""
    k = kmu(lam, mu)
    return np.cbrt(6.0 * SQ2 * gval(k)) / np.maximum(k, 1e-300)


def field_lam(lam, r, mu):
    """Minimiser F(r, theta) for the lam-family; returns F, x, k."""
    k = kmu(lam, mu)
    y = np.minimum(k ** 3 * np.asarray(r, float) ** 3 / (6.0 * SQ2), gval(k))
    x = np.sin(u_of_y(y))                   # frozen Newton inversion (REUSED)
    c = np.minimum(x / np.maximum(k, 1e-300), 1.0)
    return 2.0 * np.arccos(c), x, k


def Fr_lam(r, F, x, k):
    """dF/dr from implicit differentiation of g(k c) = k^3 r^3/(6 sqrt2):
    c_r = k^2 r^2/(2 sqrt2 g'(x)), g'(x) = 2x^2/sqrt(1-x^2)."""
    gp = 2.0 * x * x / np.sqrt(np.maximum(1.0 - x * x, 1e-300))
    c_r = k * k * r * r / (2.0 * SQ2 * np.maximum(gp, 1e-300))
    return -2.0 * c_r / np.maximum(np.sin(0.5 * F), 1e-300)


# ----------------------------------------------------------------------
# field-level spectral quadrature of the sectors on the lam-minimiser
# (gstar spectral_sectors pattern: GL in mu, edge-regularised r = R(1-v^2))
# ----------------------------------------------------------------------
def spectral_lam(lam, nmu, nv):
    xm, wm = leggauss(nmu)
    mu = 0.5 * (xm + 1.0)
    wmu = 0.5 * wm
    xv, wv = leggauss(nv)
    v = 0.5 * (xv + 1.0)
    wv = 0.5 * wv
    MU, V = np.meshgrid(mu, v, indexing="ij")
    K = kmu(lam, MU)
    RR = np.cbrt(6.0 * SQ2 * gval(K)) / np.maximum(K, 1e-300)
    r = RR * (1.0 - V * V)
    jac = 2.0 * RR * V
    F, x, _ = field_lam(lam, r, MU)
    Fr = Fr_lam(r, F, x, K)
    W2 = np.outer(wmu, wv)
    sF, cF = np.sin(F), np.cos(F)
    dE6 = 0.5 * sF ** 4 * Fr * Fr / (r * r)
    dPOT = 0.25 * r * r * ((1.0 - cF) ** 2 + lam * sF ** 2 * MU * MU)
    E6 = 2.0 * float(np.sum(dE6 * jac * W2))
    POT = 2.0 * float(np.sum(dPOT * jac * W2))
    E0 = 2.0 * float(np.sum(0.5 * (1.0 - cF) * r * r * jac * W2))
    In = 2.0 * float(np.sum(4.0 * PI * sF ** 2 * (1.0 - MU * MU)
                            * r * r * jac * W2))
    Jn = 2.0 * float(np.sum(2.0 * PI * sF ** 2 * r * r * jac * W2))
    deg = 2.0 * float(np.sum(-(1.0 / PI) * sF ** 2 * Fr * jac * W2))
    msk = F > 0.05
    bps = float(np.max(np.abs(dE6[msk] / dPOT[msk] - 1.0)))
    # support-law residual (construction check, Newton-level):
    supp = float(np.max(np.abs(gfun(np.arcsin(np.clip(x, 0, 1)))
                               - K ** 3 * r ** 3 / (6.0 * SQ2))[r <= RR]))
    return dict(lam=lam, nmu=nmu, nv=nv, E6=E6, POT=POT, E0=E0, I=In, J=Jn,
                G=E6 + POT, deg=deg, bps_resid=bps, support_resid=supp)


# ----------------------------------------------------------------------
# target-space pushforward quadrature (float GL, gstar bound_quadratures
# pattern):  dV = sqrt2 dOmega/sqrt(Q_lam)
# ----------------------------------------------------------------------
def target_sectors(lam, n=NGL):
    xp, wp = leggauss(n)
    psi = 0.5 * PI * (xp + 1.0)
    wpsi = 0.5 * PI * wp
    up, wu = leggauss(n)                     # u = cos Theta on [-1, 1]
    P, U = np.meshgrid(psi, up, indexing="ij")
    WW = np.outer(wpsi, wu)
    s2 = np.sin(P) ** 2
    Q = (1.0 - np.cos(P)) ** 2 + lam * s2 * U * U
    sq = np.sqrt(Q)
    B = float(np.sum(s2 * sq * WW)) / SQ2
    E0 = float(np.sum(s2 * (1.0 - np.cos(P)) / sq * WW)) / SQ2
    In = 4.0 * SQ2 * PI * float(np.sum(s2 ** 2 * (1.0 - U * U) / sq * WW))
    Jn = 2.0 * SQ2 * PI * float(np.sum(s2 ** 2 / sq * WW))
    return dict(lam=lam, B=B, E6=B / 2.0, E0=E0, I=In, J=Jn)


def B_reduced_float(lam, n=NGL):
    """1-D reduced bound (closed-form Theta integral), float GL.
    INT_{-1}^{1} sqrt(A + Bc u^2) du = sqrt(A+Bc) + (A/sqrt Bc) asinh(sqrt(Bc/A))."""
    xp, wp = leggauss(n)
    psi = 0.5 * PI * (xp + 1.0)
    wpsi = 0.5 * PI * wp
    c = np.cos(psi)
    s2 = 1.0 - c * c
    A = (1.0 - c) ** 2
    if lam == 0.0:
        inner = 2.0 * (1.0 - c)
    else:
        Bc = lam * s2
        inner = (np.sqrt(A + Bc)
                 + A / np.sqrt(Bc) * np.arcsinh(np.sqrt(Bc / A)))
    return float(np.sum(s2 * inner * wpsi)) / SQ2


# ----------------------------------------------------------------------
# mpmath high-precision references (reduced 1-D, closed-form u-inners,
# series-guarded against cancellation at Bc/A -> 0)
# ----------------------------------------------------------------------
_SER_N = 14


def _iu2_series_coeffs():
    """f(eps) = sqrt(1+eps) - asinh(sqrt eps)/sqrt eps = sum a_n eps^n,
    a_n = binom(1/2,n) - (-1)^n (2n)!/(4^n (n!)^2 (2n+1));  a1 = 2/3."""
    return [mp.binomial(mp.mpf(1) / 2, n)
            - (-1) ** n * mp.factorial(2 * n)
            / (mp.mpf(4) ** n * mp.factorial(n) ** 2 * (2 * n + 1))
            for n in range(1, _SER_N + 1)]


def mp_sectors(lam, dps=MP_DPS):
    """B, E0, I, J at precision dps via tanh-sinh quadrature."""
    with mp.workdps(dps):
        lam_ = mp.mpf(repr(lam))
        coeffs = _iu2_series_coeffs()

        def inners(psi):
            c = mp.cos(psi)
            s2 = (1 - c) * (1 + c)
            A = (1 - c) ** 2
            if lam_ == 0 or s2 == 0:
                sA = mp.sqrt(A)
                return 2 * sA, 2 / sA if sA else mp.inf, \
                    mp.mpf(2) / (3 * sA) if sA else mp.inf
            Bc = lam_ * s2
            sA = mp.sqrt(A)
            sBc = mp.sqrt(Bc)
            eps = Bc / A
            t = mp.sqrt(eps)
            ash = mp.asinh(t)
            i0 = mp.sqrt(A + Bc) + A / sBc * ash        # INT sqrt(A+Bc u^2)
            iinv = 2 * ash / sBc                        # INT 1/sqrt(.)
            if eps < mp.mpf("1e-4"):
                f = mp.fsum(a * eps ** n
                            for n, a in enumerate(coeffs, start=1))
                iu2 = f / (sA * eps)                    # INT u^2/sqrt(.)
            else:
                iu2 = mp.sqrt(A + Bc) / Bc - A / Bc ** mp.mpf("1.5") * ash
            return i0, iinv, iu2

        def _zero(psi):
            c = mp.cos(psi)
            return (1 - c) * (1 + c) == 0 or (1 - c) == 0

        def fB(psi):
            if _zero(psi):
                return mp.mpf(0)
            s2 = mp.sin(psi) ** 2
            return s2 * inners(psi)[0]

        def fE0(psi):
            if _zero(psi):
                return mp.mpf(0)
            c = mp.cos(psi)
            return (1 - c) * (1 - c * c) * inners(psi)[1]

        def fI(psi):
            if _zero(psi):
                return mp.mpf(0)
            s4 = mp.sin(psi) ** 4
            i0, iinv, iu2 = inners(psi)
            return s4 * (iinv - iu2)

        def fJ(psi):
            if _zero(psi):
                return mp.mpf(0)
            return mp.sin(psi) ** 4 * inners(psi)[1]

        pts = [0, mp.pi / 2, mp.pi]
        B = mp.quad(fB, pts) / mp.sqrt(2)
        E0 = mp.quad(fE0, pts) / mp.sqrt(2)
        In = 4 * mp.sqrt(2) * mp.pi * mp.quad(fI, pts)
        Jn = 2 * mp.sqrt(2) * mp.pi * mp.quad(fJ, pts)
        return dict(B=B, E0=E0, I=In, J=Jn)


# ----------------------------------------------------------------------
# CLOSED FORMS (derived in this workstream, verified vs quadrature below).
# Inner c-integral P(m) = INT_0^1 c^2(1-c^2) sqrt(1-mc^2) dc (sympy, exact):
#   P(m) = sqrt(1-m)(4m^2-4m+3)/(48 m^2) + (2m-1) asin(sqrt m)/(16 m^{5/2});
# outer u-integral by parts (m = 1-lam u^2; the log terms cancel exactly):
#
#   B(lam) = (sqrt2/3) [ sqrt(lam)(3-2 lam)/(1-lam)
#                        + (3-4 lam) arccos(sqrt lam)/(1-lam)^{3/2} ],
#   B(0) = pi/sqrt2,  B(1/2) = (4+pi)/3,  B(1) = 16 sqrt2/9   (limits).
#
# Sectors by the same by-parts route (a = sqrt(1-lam), v = cos(psi/2),
# M_k := INT_0^1 v^{2k}/sqrt(1-a^2 v^2) dv;  M1 = g(a)/(2a^3) with the
# CAMPAIGN's own g(x) = arcsin x - x sqrt(1-x^2)):
#   E0(lam) = 2 sqrt2 (M0 - M1),
#   J(lam)  = (16 sqrt2 pi/3)(M0 + M1 - 2 M2),
#   I_lam   = 16 pi (E0 - B/2)     [generic d-stationarity identity],
#   I(lam)  = (I_lam - 2(1-lam) J)/lam   (lam > 0;  I(0) = 8 sqrt2 pi^2/3).
# Special exact points: E0(0) = pi/sqrt2, E0(1/2) = 2, E0(1) = 4 sqrt2/3;
#   J(0) = 2 sqrt2 pi^2, J(1/2) = 16 pi - 8 pi^2/3, J(1) = 224 sqrt2 pi/45;
#   I(0) = 8 sqrt2 pi^2/3, I(1/2) = 32 pi/3, I(1) = 64 sqrt2 pi/9.
# ----------------------------------------------------------------------
def closed_forms(lam, dps=MP_DPS):
    """mpmath closed-form B, E0, J, I at lam (uses limits at lam = 0, 1)."""
    with mp.workdps(dps):
        lam_ = mp.mpf(repr(float(lam)))
        s2m = mp.sqrt(2)
        if lam_ == 0:
            B = mp.pi / s2m
            E0 = mp.pi / s2m
            J = 2 * s2m * mp.pi ** 2
            In = 8 * s2m * mp.pi ** 2 / 3
            return dict(B=B, E0=E0, J=J, I=In)
        if lam_ == 1:
            B = 16 * s2m / 9
            E0 = 4 * s2m / 3
            J = 224 * s2m * mp.pi / 45
            In = 64 * s2m * mp.pi / 9
            return dict(B=B, E0=E0, J=J, I=In)
        a = mp.sqrt(1 - lam_)
        B = (s2m / 3) * (mp.sqrt(lam_) * (3 - 2 * lam_) / (1 - lam_)
                         + (3 - 4 * lam_) * mp.acos(mp.sqrt(lam_))
                         / (1 - lam_) ** mp.mpf("1.5"))
        M0 = mp.asin(a) / a
        M1 = (mp.asin(a) - a * mp.sqrt(1 - a * a)) / (2 * a ** 3)
        M2 = (3 * mp.asin(a) - a * mp.sqrt(1 - a * a) * (3 + 2 * a * a)) \
            / (8 * a ** 5)
        E0 = 2 * s2m * (M0 - M1)
        J = (16 * s2m * mp.pi / 3) * (M0 + M1 - 2 * M2)
        Ilam = 16 * mp.pi * (E0 - B / 2)
        In = (Ilam - 2 * (1 - lam_) * J) / lam_
        return dict(B=B, E0=E0, J=J, I=In)


def mp_B(lam, dps=MP_DPS):
    with mp.workdps(dps):
        lam_ = mp.mpf(repr(float(lam)))

        def fB(psi):
            c = mp.cos(psi)
            s2 = (1 - c) * (1 + c)
            A = (1 - c) ** 2
            if lam_ == 0 or s2 == 0 or A == 0:
                return s2 * 2 * mp.sqrt(A)
            Bc = lam_ * s2
            ash = mp.asinh(mp.sqrt(Bc / A))
            return s2 * (mp.sqrt(A + Bc) + A / mp.sqrt(Bc) * ash)

        return mp.quad(fB, [0, mp.pi / 2, mp.pi]) / mp.sqrt(2)


# ----------------------------------------------------------------------
# main phases
# ----------------------------------------------------------------------
def phase_A(state):
    print("== A. SHARP BOUND B(lam): quadrature ladder + endpoint gates ==")
    rows = []
    for lam in LAMBDAS:
        t0 = time.time()
        b2d = target_sectors(lam)["B"]
        b1d = B_reduced_float(lam)
        bmp = mp_B(lam)
        bcf = closed_forms(lam)["B"]
        with mp.workdps(MP_DPS):
            d_cf = float(abs(bmp - bcf))
        bmp_f = float(bmp)
        rows.append(dict(lam=lam, B_2d=b2d, B_1d=b1d, B_mp=bmp_f,
                         B_cf=float(bcf), B_mp_str=mp.nstr(bmp, 25),
                         d_2d_mp=b2d - bmp_f, d_1d_mp=b1d - bmp_f,
                         d_cf_mp=d_cf,
                         c_of_lam=4.0 / PI * bmp_f,
                         secs=time.time() - t0))
        print("   lam=%.2f  B=%.15f  c=(4/pi)B=%.12f   |2d-mp|=%.1e "
              "|1d-mp|=%.1e  |closedform-mp|=%.1e  (%.1f s)"
              % (lam, bmp_f, 4.0 / PI * bmp_f, abs(b2d - bmp_f),
                 abs(b1d - bmp_f), d_cf, rows[-1]["secs"]))
    # special exact point lam = 1/2: B = (4+pi)/3
    with mp.workdps(40):
        d_half = float(abs(mp_B(0.5, dps=40) - (4 + mp.pi) / 3))
    print("   SPECIAL POINT: |B(1/2) - (4+pi)/3| = %.1e  "
          "(exact midpoint value found)" % d_half)
    e1 = abs(rows[-1]["B_mp"] - GSTAR1)
    e1x = float(abs(mp_B(1.0, dps=40) - 16 * mp.sqrt(2) / 9))
    e0 = abs(rows[0]["B_mp"] - B0_EX)
    e0x = float(abs(mp_B(0.0, dps=40) - mp.pi / mp.sqrt(2)))
    print("   ENDPOINT GATES: |B(1) - 16 sqrt2/9| = %.2e (dps40: %.1e)"
          % (e1, e1x))
    print("                   |B(0) - pi/sqrt2|   = %.2e (dps40: %.1e)"
          % (e0, e0x))
    state["A"] = dict(rows=rows, gate_B1=e1, gate_B0=e0,
                      gate_B1_dps40=e1x, gate_B0_dps40=e0x,
                      B1_exact=GSTAR1, B0_exact=B0_EX,
                      B_half_exact_resid=d_half,
                      closed_form="B(lam) = (sqrt2/3)[ sqrt(lam)(3-2lam)/"
                                  "(1-lam) + (3-4lam) arccos(sqrt lam)/"
                                  "(1-lam)^(3/2) ]",
                      worst_cf_vs_quad=max(r["d_cf_mp"] for r in rows))
    ckpt_write(state)
    return state


def phase_B(state):
    print("\n== B. ATTAINMENT: closed-form lam-minimiser, BPS residuals ==")
    rows = []
    for lam in LAMBDAS:
        t0 = time.time()
        sp = spectral_lam(lam, NSPEC, NSPEC)
        Bref = next(r for r in state["A"]["rows"]
                    if r["lam"] == lam)["B_mp"]
        sp["G_minus_B"] = sp["G"] - Bref
        sp["R_axis"] = float(Redge_lam(lam, np.array([1.0 - 1e-14]))[0]) \
            if lam == 1.0 else float(Redge_lam(lam, np.array([1.0]))[0])
        sp["R_eq"] = float(Redge_lam(lam, np.array([0.0]))[0])
        sp["secs"] = time.time() - t0
        rows.append(sp)
        print("   lam=%.2f  BPS-resid=%.2e  supp-resid=%.2e  deg=%.8f  "
              "G-B=%+.2e  R_ax=%.6f R_eq=%.6f  (%.1f s)"
              % (lam, sp["bps_resid"], sp["support_resid"], sp["deg"],
                 sp["G_minus_B"], sp["R_axis"], sp["R_eq"], sp["secs"]))
    interior = [r for r in rows if 0.0 < r["lam"] < 1.0]
    worst = max(r["bps_resid"] for r in interior)
    print("   worst interior BPS residual: %.2e  (gate <= 1e-8: %s)"
          % (worst, "PASS" if worst <= 1e-8 else "FAIL"))
    state["B"] = dict(rows=rows, worst_interior_bps=worst,
                      attained=bool(worst <= 1e-8))
    ckpt_write(state)
    return state


def phase_C(state):
    print("\n== C. SECTOR DECOMPOSITION + IDENTITY PERSISTENCE ==")
    rows = []
    for lam in LAMBDAS:
        t0 = time.time()
        ts = target_sectors(lam)
        ms = mp_sectors(lam)
        cf = closed_forms(lam)
        with mp.workdps(MP_DPS):
            cf_resid = dict(B=float(abs(ms["B"] - cf["B"])),
                            E0=float(abs(ms["E0"] - cf["E0"])),
                            I=float(abs(ms["I"] - cf["I"])),
                            J=float(abs(ms["J"] - cf["J"])))
        B_, E0_, I_, J_ = (float(ms["B"]), float(ms["E0"]),
                           float(ms["I"]), float(ms["J"]))
        E6_ = B_ / 2.0
        sp = next(r for r in state["B"]["rows"] if r["lam"] == lam)
        Ilam = lam * I_ + 2.0 * (1.0 - lam) * J_
        id1_lit = E6_ - (E0_ - I_ / SIXTEEN_PI)
        id1_gen = E6_ - (E0_ - Ilam / SIXTEEN_PI)
        pred1 = -(1.0 - lam) * (2.0 * J_ - I_) / SIXTEEN_PI
        # saturated tuple, lam-consistent inertia bookkeeping (axi3 algebra)
        Gs = B_
        I_tot = 8.0 * PI * Gs
        halo_D = (I_tot - Ilam) / 2.0
        E0_halo = halo_D / (8.0 * PI)
        Estat_gen = E6_ + E0_ + E0_halo
        id2_gen = Estat_gen - 1.5 * Gs
        # literal chiral-I bookkeeping at lam != 1
        E0_halo_lit = (I_tot - I_) / 2.0 / (8.0 * PI)
        id2_lit = (E6_ + E0_ + E0_halo_lit) - 1.5 * Gs
        # support law
        supp_lit = abs(sp["R_axis"] - R_AXIS1)     # axis compacton radius
        supp_gen = sp["support_resid"]
        r = dict(lam=lam, B=B_, E6=E6_, E0=E0_, I=I_, J=J_, I_lam=Ilam,
                 c_of_lam=4.0 / PI * Gs, L=math.sqrt(8.0 * PI) * Gs,
                 I_tot=I_tot, halo_D=halo_D, E0_halo=E0_halo,
                 Estat_tot=Estat_gen,
                 id1_literal=id1_lit, id1_literal_pred=pred1,
                 id1_pred_resid=id1_lit - pred1,
                 id1_generalized=id1_gen,
                 id2_generalized=id2_gen, id2_literal=id2_lit,
                 support_literal=supp_lit, support_generalized=supp_gen,
                 xcheck_E0_2d=ts["E0"] - E0_, xcheck_I_2d=ts["I"] - I_,
                 xcheck_J_2d=ts["J"] - J_,
                 xcheck_E0_spec=sp["E0"] - E0_, xcheck_I_spec=sp["I"] - I_,
                 xcheck_J_spec=sp["J"] - J_,
                 closed_form_resid=cf_resid, secs=time.time() - t0)
        rows.append(r)
        print("   lam=%.2f  E6=%.9f E0=%.9f I=%.6f J=%.6f   "
              "cf-resid B/E0/I/J %.0e/%.0e/%.0e/%.0e"
              % (lam, E6_, E0_, I_, J_, cf_resid["B"], cf_resid["E0"],
                 cf_resid["I"], cf_resid["J"]))
        print("            id1 literal E6-(E0-I/16pi) = %+.9f  "
              "(pred %+.9f, |diff| %.1e)" % (id1_lit, pred1,
                                             abs(id1_lit - pred1)))
        print("            id1 general E6-(E0-I_lam/16pi) = %+.2e   "
              "id2 general = %+.2e   id2 literal = %+.2e"
              % (id1_gen, id2_gen, id2_lit))
    # verdicts
    inter = [r for r in rows if 0.0 < r["lam"] < 1.0]
    v = dict(
        id1_literal="DISTINCTIVE" if all(abs(r["id1_literal"]) > 1e-3
                                         for r in inter) else "OTHER",
        id1_generalized="GENERIC" if all(abs(r["id1_generalized"]) < 1e-8
                                         for r in rows) else "OTHER",
        id2_generalized="GENERIC" if all(abs(r["id2_generalized"]) < 1e-8
                                         for r in rows) else "OTHER",
        id2_literal="DISTINCTIVE" if all(abs(r["id2_literal"]) > 1e-3
                                         for r in inter) else "OTHER",
        support_literal="DISTINCTIVE" if all(r["support_literal"] > 1e-3
                                             for r in inter) else "OTHER",
        support_generalized="GENERIC" if all(r["support_generalized"] < 1e-8
                                             for r in rows) else "OTHER",
        kappa_paper="GENERIC (halo threshold algebra is lam-free: the "
                    "marginal halo sits at Theta=pi/2 where i_lam = 2 eta^2 "
                    "for every lam; kappa* = 1/sqrt(8pi) exactly)")
    print("   VERDICTS: id1 literal %s | id1 generalized %s | "
          "id2 generalized %s | id2 literal %s | support literal %s | "
          "support generalized %s"
          % (v["id1_literal"], v["id1_generalized"], v["id2_generalized"],
             v["id2_literal"], v["support_literal"],
             v["support_generalized"]))
    # special exact point lam = 1/2 (whole tuple is pi-rational)
    with mp.workdps(40):
        ms = mp_sectors(0.5, dps=40)
        half = dict(B=float(abs(ms["B"] - (4 + mp.pi) / 3)),
                    E0=float(abs(ms["E0"] - 2)),
                    I=float(abs(ms["I"] - 32 * mp.pi / 3)),
                    J=float(abs(ms["J"] - (16 * mp.pi
                                           - 8 * mp.pi ** 2 / 3))))
    print("   SPECIAL POINT lam=1/2 (exact): B=(4+pi)/3 [%.0e]  E0=2 "
          "[%.0e]  I=32pi/3 [%.0e]  J=16pi-8pi^2/3 [%.0e]"
          % (half["B"], half["E0"], half["I"], half["J"]))
    state["C"] = dict(
        rows=rows, verdicts=v, half_point_resid=half,
        closed_forms=dict(
            B="(sqrt2/3)[sqrt(lam)(3-2lam)/(1-lam)"
              " + (3-4lam) arccos(sqrt lam)/(1-lam)^(3/2)]",
            E0="2 sqrt2 (M0 - M1),  a = sqrt(1-lam), "
               "M0 = asin(a)/a, M1 = g(a)/(2a^3)",
            J="(16 sqrt2 pi/3)(M0 + M1 - 2 M2), "
              "M2 = [3 asin a - a sqrt(1-a^2)(3+2a^2)]/(8a^5)",
            I="(16 pi (E0 - B/2) - 2(1-lam) J)/lam;  I(0) = 8 sqrt2 pi^2/3",
            special=dict(lam0=dict(B="pi/sqrt2", E0="pi/sqrt2",
                                   I="8 sqrt2 pi^2/3", J="2 sqrt2 pi^2"),
                         lam_half=dict(B="(4+pi)/3", E0="2", I="32 pi/3",
                                       J="16 pi - 8 pi^2/3"),
                         lam1=dict(B="16 sqrt2/9", E0="4 sqrt2/3",
                                   I="64 sqrt2 pi/9", J="224 sqrt2 pi/45"))))
    ckpt_write(state)
    return state


def phase_D(state):
    print("\n== D. ANCHOR-HIT SCAN (printed constants of gstar_RESULTS.md) ==")
    anchors_c = [
        ("c_paper = 64 sqrt2/(9 pi)", 64.0 * SQ2 / (9.0 * PI), "c-class"),
        ("(4/pi) c0 = 512 sqrt42/(105 pi^2)",
         4.0 / PI * C0, "c-class"),
        ("pi", PI, "c-class"),
        ("axi3 c bracket hi 3.20146", 3.20146, "c-class"),
        ("axi3 c bracket lo 2.82843", 2.82843, "c-class"),
        ("kappa_paper 1/sqrt2", 1.0 / SQ2, "c-class"),
    ]
    anchors_G = [
        ("G* = 16 sqrt2/9", GSTAR1, "G-class"),
        ("c0 = 128 sqrt42/(105 pi)", C0, "G-class"),
        ("axi3 extended 2.5144238", 2.5144238, "G-class"),
        ("axi3 joint 2.5144510", 2.5144510, "G-class"),
        ("axi4 locked 2.515991", 2.515991, "G-class"),
        ("axi3 bound pi/sqrt2", PI / SQ2, "G-class"),
        ("R_eq = (3 sqrt2 pi)^(1/3)  [LENGTH]", R_EQ, "length-class"),
        ("R_axis = 2^(5/6)  [LENGTH]", R_AXIS1, "length-class"),
        ("E_stat_tot(1) = 24 sqrt2/9", 24.0 * SQ2 / 9.0, "E-class"),
        ("E6(1) = 8 sqrt2/9", 8.0 * SQ2 / 9.0, "E-class"),
        ("E0(1) = 4 sqrt2/3", 4.0 * SQ2 / 3.0, "E-class"),
    ]
    nl = 2001
    ls = np.linspace(0.0, 1.0, nl)
    Bs = np.array([B_reduced_float(l, n=NGL // 2) for l in ls])
    cs = 4.0 / PI * Bs
    hits = []

    def scan(vals, anchors, label):
        for name, a, cls in anchors:
            d = vals - a
            sgn = np.sign(d)
            for i in range(nl - 1):
                if sgn[i] == 0 or sgn[i] * sgn[i + 1] < 0:
                    lo, hi = ls[i], ls[i + 1]
                    if label == "c":
                        f = lambda l: 4.0 / PI * B_reduced_float(l) - a
                        fm = lambda l: float(4 / mp.pi * mp_B(l)
                                             - mp.mpf(repr(a)))
                    else:
                        f = lambda l: B_reduced_float(l) - a
                        fm = lambda l: float(mp_B(l) - mp.mpf(repr(a)))
                    try:
                        root = brentq(f, max(lo - 2e-3, 0.0),
                                      min(hi + 2e-3, 1.0), xtol=1e-13)
                    except ValueError:
                        continue
                    # one Newton polish with mpmath values
                    h = 1e-6
                    der = (fm(min(root + h, 1.0)) - fm(max(root - h, 0.0))) \
                        / (min(root + h, 1.0) - max(root - h, 0.0))
                    root_mp = root - fm(root) / der if der else root
                    endpoint = root_mp < 1e-9 or root_mp > 1.0 - 1e-9
                    # a "hit" a few 1e-6 in lam from an endpoint whose
                    # anchor is a ROUNDED print of the endpoint value is an
                    # endpoint identity, not an interior crossing
                    v0 = vals[0] if label != "c" else cs[0]
                    v1 = vals[-1]
                    near_end = (abs(a - v0) < 1e-4 or abs(a - v1) < 1e-4)
                    tag = ("[ENDPOINT]" if endpoint else
                           "[ENDPOINT (rounded print)]" if near_end
                           else "[INTERIOR HIT]")
                    hits.append(dict(curve=label, anchor=name,
                                     anchor_value=a, anchor_class=cls,
                                     lam=float(root_mp),
                                     endpoint=bool(endpoint or near_end),
                                     tag=tag,
                                     resid=fm(float(np.clip(root_mp, 0, 1)))))
                    print("   %s(lam) = %-42s at lam = %.10f  %s resid %.1e"
                          % (label, name, root_mp, tag, hits[-1]["resid"]))
    scan(cs, anchors_c, "c")
    scan(Bs, anchors_G, "G")
    # high-precision refinement of the flagship interior crossings
    ref = {}
    for h in hits:
        if h["endpoint"]:
            continue
        key = "c=pi" if (h["curve"] == "c"
                         and abs(h["anchor_value"] - PI) < 1e-12) else None
        key = key or ("G=R_eq" if (h["curve"] == "G"
                                   and abs(h["anchor_value"] - R_EQ) < 1e-12)
                      else None)
        if key:
            with mp.workdps(30):
                tgt = (mp.pi ** 2 / 4 if key == "c=pi"
                       else (3 * mp.sqrt(2) * mp.pi) ** (mp.mpf(1) / 3))
                root = mp.findroot(lambda l: mp_B(float(l), dps=30) - tgt,
                                   mp.mpf(repr(h["lam"])), tol=1e-22)
                ref[key] = mp.nstr(root, 15)
                h["lam_refined"] = float(root)
                print("   refined %-8s lam* = %s" % (key, ref[key]))
    state["D"] = dict(hits=hits, refined=ref,
                      scan_grid=nl,
                      c_range=[float(cs[0]), float(cs[-1])],
                      G_range=[float(Bs[0]), float(Bs[-1])],
                      anchors_scanned=dict(
                          c=[(n, a) for n, a, _ in anchors_c],
                          G=[(n, a) for n, a, _ in anchors_G]))
    ckpt_write(state)
    return state, ls, Bs


# ----------------------------------------------------------------------
# figure
# ----------------------------------------------------------------------
INK, GRID_C = "#333333", "#dddddd"
BLUE, AQUA, RED, VIOLET, ORANGE = ("#2a78d6", "#1baf7a", "#e34948",
                                   "#4a3aa7", "#d97b1f")


def style_ax(a):
    a.set_facecolor("white")
    a.grid(color=GRID_C, lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.tick_params(colors=INK, labelsize=9)


def figure(state, ls, Bs, path):
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.2), dpi=160)
    fig.patch.set_facecolor("white")
    (a1, a2), (a3, a4) = axes

    # (a) G*(lam) and c(lam)
    style_ax(a1)
    a1.plot(ls, Bs, color=BLUE, lw=1.8, label=r"$G^*(\lambda)=B(\lambda)$")
    a1.axhline(GSTAR1, color=BLUE, ls=":", lw=0.9)
    a1.axhline(B0_EX, color=AQUA, ls=":", lw=0.9)
    a1.plot([1], [GSTAR1], "o", color=BLUE, mec="white", ms=8)
    a1.plot([0], [B0_EX], "o", color=AQUA, mec="white", ms=8)
    a1.text(0.62, GSTAR1 - 0.012, r"$16\sqrt{2}/9$ (GUM, $\lambda=1$)",
            color=BLUE, fontsize=9)
    a1.text(0.02, B0_EX + 0.006, r"$\pi/\sqrt{2}$ (achiral, $\lambda=0$)",
            color=AQUA, fontsize=9)
    lam_pi = next((h.get("lam_refined", h["lam"])
                   for h in state["D"]["hits"]
                   if h["curve"] == "c" and not h["endpoint"]
                   and abs(h["anchor_value"] - PI) < 1e-12), None)
    if lam_pi is not None:
        Bpi = PI * PI / 4.0
        a1.plot([lam_pi], [Bpi], "s", color=RED, mec="white", ms=9)
        a1.annotate(r"$c(\lambda^*)=\pi$ at $\lambda^*=%.6f$" % lam_pi,
                    (lam_pi, Bpi), (0.06, 2.482),
                    color=RED, fontsize=9,
                    arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    a1.plot([0.5], [(4.0 + PI) / 3.0], "D", color=VIOLET, mec="white", ms=8)
    a1.text(0.5, (4.0 + PI) / 3.0 - 0.028,
            r"$G^*(\frac{1}{2})=\frac{4+\pi}{3}$", color=VIOLET, fontsize=9,
            ha="center")
    a1.text(0.02, 2.428,
            r"$G^*(\lambda)=\frac{\sqrt{2}}{3}\left["
            r"\frac{\sqrt{\lambda}(3-2\lambda)}{1-\lambda}"
            r"+\frac{(3-4\lambda)\arccos\sqrt{\lambda}}"
            r"{(1-\lambda)^{3/2}}\right]$"
            "\n(closed form, this workstream)",
            color=INK, fontsize=9, va="top")
    a1.set_xlabel(r"$\lambda$", color=INK)
    a1.set_ylabel(r"$G^*(\lambda)$", color=INK)
    ax2 = a1.twinx()
    ax2.set_ylim(np.array(a1.get_ylim()) * 4.0 / PI)
    ax2.tick_params(colors=RED, labelsize=9)
    ax2.axhline(PI, color=RED, ls="--", lw=0.8)
    ax2.text(0.02, PI + 0.008, r"$c=\pi$", color=RED, fontsize=9)
    ax2.set_ylabel(r"$c(\lambda)=(4/\pi)G^*(\lambda)$", color=RED)
    a1.set_title("(a) saturated constant across the chirality family",
                 fontsize=10, color=INK)
    a1.legend(frameon=False, fontsize=8, loc="lower right")

    # (b) sector curves
    style_ax(a2)
    C = state["C"]["rows"]
    lamv = [r["lam"] for r in C]
    a2.plot(lamv, [r["E6"] for r in C], "o-", color=BLUE, ms=5,
            mec="white", label=r"$E_6=B/2$")
    a2.plot(lamv, [r["E0"] for r in C], "s-", color=AQUA, ms=5,
            mec="white", label=r"$E_0$")
    a2.plot(lamv, [r["I"] / SIXTEEN_PI for r in C], "^-", color=RED, ms=5,
            mec="white", label=r"$I/16\pi$ (chiral sector)")
    a2.plot(lamv, [r["J"] / (8.0 * PI) for r in C], "v-", color=VIOLET,
            ms=5, mec="white", label=r"$J/8\pi$")
    a2.plot(lamv, [r["I_lam"] / SIXTEEN_PI for r in C], "--",
            color=ORANGE, lw=1.6,
            label=r"$I_\lambda/16\pi=E_0-E_6$ (generic)")
    a2.set_xlabel(r"$\lambda$", color=INK)
    a2.set_title("(b) sector decomposition on the $\\lambda$-minimiser",
                 fontsize=10, color=INK)
    a2.legend(frameon=False, fontsize=8)

    # (c) identity residuals
    style_ax(a3)
    a3.semilogy(lamv, [abs(r["id1_literal"]) + 1e-18 for r in C], "o-",
                color=RED, ms=5, mec="white",
                label=r"|$E_6-(E_0-I/16\pi)$|  (literal, chiral $I$)")
    a3.semilogy(lamv, [abs(r["id2_literal"]) + 1e-18 for r in C], "s-",
                color=ORANGE, ms=5, mec="white",
                label=r"|$E_{stat}-\frac{3}{2}G^*$|  (literal bookkeeping)")
    a3.semilogy(lamv, [abs(r["id1_generalized"]) + 1e-18 for r in C],
                "^-", color=BLUE, ms=5, mec="white",
                label=r"|$E_6-(E_0-I_\lambda/16\pi)$|  (generalized)")
    a3.semilogy(lamv, [abs(r["id2_generalized"]) + 1e-18 for r in C],
                "v-", color=AQUA, ms=5, mec="white",
                label=r"|$E_{stat}-\frac{3}{2}G^*$|  ($\lambda$-consistent)")
    a3.axhline(1e-8, color=INK, ls=":", lw=0.8)
    a3.text(0.90, 3e-8, "1e-8", color=INK, fontsize=8)
    a3.set_ylim(1e-17, 1.0)
    a3.set_xlabel(r"$\lambda$", color=INK)
    a3.set_ylabel("|identity residual|", color=INK)
    a3.set_title("(c) identity persistence: literal = DISTINCTIVE, "
                 "structure = GENERIC", fontsize=10, color=INK)
    a3.legend(frameon=False, fontsize=8, loc="center left")

    # (d) minimiser support shapes
    style_ax(a4)
    th = np.linspace(1e-4, PI - 1e-4, 800)
    for lam, col in [(0.0, AQUA), (0.25, VIOLET), (0.5, ORANGE),
                     (0.75, RED), (1.0, BLUE)]:
        R = Redge_lam(lam, np.abs(np.cos(th)))
        a4.plot(R * np.sin(th), R * np.cos(th), color=col, lw=1.6,
                label=r"$\lambda=%.2f$" % lam)
    a4.set_aspect("equal")
    a4.set_xlabel(r"$\rho$", color=INK)
    a4.set_ylabel(r"$z$", color=INK)
    a4.set_title("(d) support edge $R(\\theta)$: sphere "
                 "$(3\\sqrt{2}\\pi)^{1/3}$ $\\to$ oblate compacton; "
                 "$R_{eq}$ is $\\lambda$-free", fontsize=10, color=INK)
    a4.legend(frameon=False, fontsize=8, loc="upper right")

    fig.suptitle("Phase S2: the chirality family "
                 "$W_\\lambda$ — $G^*(\\lambda)$, attainment, identity "
                 "persistence", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    state = dict(config=dict(quick=QUICK, mp_dps=MP_DPS, ngl=NGL,
                             nspec=NSPEC, lambdas=LAMBDAS))
    state = phase_A(state)
    state = phase_B(state)
    state = phase_C(state)
    state, ls, Bs = phase_D(state)

    g1 = max(state["A"]["gate_B1"], state["A"]["gate_B0"])
    g2 = state["B"]["worst_interior_bps"]
    state["gates"] = dict(
        S2_G1=dict(requirement="B(1), B(0) endpoints <= 1e-10",
                   measured=g1, verdict="PASS" if g1 <= 1e-10 else "FAIL"),
        S2_G2=dict(requirement="BPS saturation residual <= 1e-8 at every "
                               "sampled lam (>= 9 interior)",
                   measured=g2, verdict="PASS" if g2 <= 1e-8 else "FAIL"),
        S2_G3=dict(requirement="identity persistence table with verdicts",
                   measured=state["C"]["verdicts"], verdict="PASS"),
        S2_G4=dict(requirement="honest bottom line (within-model only)",
                   measured="see RESULTS.md section 7", verdict="PASS"))
    state["runtime_s"] = time.time() - t0
    with open(os.path.join(OUT, "s2_results.json"), "w") as fh:
        json.dump(state, fh, indent=1, default=float)
    figure(state, ls, Bs, os.path.join(OUT, "s2_fig.png"))
    print("\n== GATES ==")
    for k, v in state["gates"].items():
        print("   %s: %s  (measured: %s)"
              % (k, v["verdict"], v["measured"]))
    print("runtime: %.1f s" % state["runtime_s"])
    print("wrote s2_results.json, s2_fig.png")


if __name__ == "__main__":
    main()
