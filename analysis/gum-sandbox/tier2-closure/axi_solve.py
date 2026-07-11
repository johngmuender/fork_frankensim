#!/usr/bin/env python3
"""
axi_solve.py — Tier-2b STEP 2: the spheroidal-ansatz isorotating closure at
eps = 0.05, at FIELD LEVEL, with ZERO calibrated constants.

Builds on step 1 (radial_solve.py / radial_results.json): hedgehog profile f(r)
solved at a6 = a0 = m = 1, a2 = a4 = t = 0.008276434949296802 (frozen), giving
the radial sector integrals E2, E4, E6, E0 and eps = t(E2+E4)/(E6+E0) = 0.049985.

----------------------------------------------------------------------------
FIELD, CURRENTS, DENSITIES (all constants DERIVED, then VERIFIED in gate V1)
----------------------------------------------------------------------------
Unit quaternion q = (q0, q) with P = q0*1 + i q.sigma; hedgehog q0 = cos f,
q = sin f * xhat. Axisymmetric winding-1 fields obey
q(rho,z,phi) = R_phi q(rho,z,0) R_phi^dag, R_phi = exp(i phi sigma3/2), so
everything is evaluated on the phi = 0 half-plane and
d_phi P = (i/2)[sigma3, P] is ANALYTIC (no phi differencing).

With L_i = P^dag d_i P = i a_i.sigma and the orthonormal-frame gradients
D = (L_rho, L_z, L_phi/rho), su(2) identities give
  -(1/2) Tr(D_i D_i)              = sum |a_i|^2
  -(1/8) Tr([D_i,D_j]^2)          = |a_i x a_j|^2
  b = -(1/(24 pi^2)) eps_ijk Tr(D_i D_j D_k)
    = -(1/(2 pi^2)) a_rho.(a_z x a_phi/rho)   [(rho,z,phi) ordering -> B=+1
       for f(0)=pi; sign verified by the degree in gate V1]
On the phi = 0 plane q2 = 0 and (from d_phi P): a_phi/rho = (q1/rho)(q3,-q0,-q1),
with q1/rho supplied analytically by each field builder (regular on the axis).

Sector normalisations are DERIVED by evaluating each density on the hedgehog
and matching the step-1 radial-integral convention (E_rad = E_3D / 4pi):
  E2 = (1/4pi) INT sum|a_i|^2 dV          -> INT (f'^2 r^2 + 2 sin^2 f) dr
  E4 = (1/4pi) INT sum_{i<j}|a_i x a_j|^2 -> INT sin^2 f (2f'^2 + s^2/r^2) dr
  E6 = pi^3    INT b^2 dV                 -> INT sin^4 f f'^2 / r^2 dr
  E0 = (1/4pi) INT (1 - q0) dV            -> INT (1-cos f) r^2 dr
  I  =         INT 2 (q1^2 + q2^2) dV     =  (16 pi/3) INT r^2 sin^2 f dr
(dV = 2 pi rho drho dz, doubled for z<0; z-parity P(rho,-z) = s1 P(rho,z) s1
makes every density even in z — used and stated, quarter-plane grid.)
These constants (1/4pi, 1/4pi, pi^3, 1/4pi, 1) are frozen BEFORE any closure;
gate V1/V2 verifies each against the step-1 numbers to < 0.5%. NOTHING is
fitted to make the closure work.

----------------------------------------------------------------------------
DEFORMATION FAMILIES
----------------------------------------------------------------------------
Family A — the task's literal deformation (an SDiff element + dilation):
  q_(lam,d)(x,y,z) = q_base(x/(e d), y/(e d), z/(p d)),
  e = lam^(-1/3), p = lam^(2/3)  (volume preserving, e^2 p = 1).
  Pulling the integrals back to base coordinates (X,Z) = (rho/(ed), z/(pd))
  gives EXACT closed-form (lam, d) dependence from ONE base evaluation:
    E2(lam,d) = (1/4pi) d [ lam^(2/3)(S11+S33) + lam^(-4/3) S22 ]
    E4(lam,d) = (1/4pi) (1/d) [ lam^(-2/3)(C12+C23) + lam^(4/3) C13 ]
    E6(lam,d) = pi^3 BQ / d^3            (lam-independent: SDiff-flat)
    E0(lam,d) = (1/4pi) POT d^3          (lam-independent: SDiff-flat)
    I (lam,d) = INR d^3                  (lam-independent — see FINDING)
  where S/C/BQ/POT/INR are the base (sphere) 2-D integrals of |a_rho|^2 etc.
  FINDING (exact theorem, verified numerically in V3/V4): every functional
  with no spatial gradients — the potential E0 AND the isorotation inertia
  I = INT 2(q1^2+q2^2) dV — is invariant under EVERY volume-preserving
  diffeomorphism (pure change of variables). Hence in the paper's own
  deformation class g(lam) == 1 identically: the claimed spheroidal inertia
  enhancement g_oblate(lam) is NOT SDiff-realisable at field level.

Family B — the paper's actual spheroidal-shell ansatz (NOT an SDiff image):
  q = ( cos f(m), sin f(m) xhat_lab ),  m^2 = (rho/e)^2 + (z/p)^2, then x -> x/d.
  Here the profile is deformed onto spheroidal shells but the iso-direction
  stays the LAB hedgehog xhat. This reproduces the paper's inertia curve
  exactly: I_B(lam)/I_B(1) = <sin^2 theta_lab>_spheroid / (2/3) = g_oblate(lam)
  (verified against closure_reduced.g_oblate in gate V4) — at the price that
  E6 is no longer flat (the ansatz leaves the SDiff orbit). Family B is the
  minimal field-level realisation of the paper's g-mechanism; both closures
  are run and adjudicated.
  Dilation scaling for BOTH families is exact: E2 ~ d, E4 ~ 1/d, E6 ~ 1/d^3,
  E0 ~ d^3, I ~ d^3 (q_(lam,d)(x) = q_(lam,1)(x/d)).

----------------------------------------------------------------------------
CLOSURE (per family)
----------------------------------------------------------------------------
E_static = t E2 + t E4 + E6 + E0 (radial-integral convention, as step 1),
I as the full 3-D integral (exactly as the task defines both).
Routhian R(lam,d;L) = E_static + L^2/(2I); minimise over (lam,d)
(lam-grid + golden refine; golden in ln d — deterministic, no RNG);
clock condition 2 L^2 / I* = R*(L) solved by bisection on L.
(The task's bracket L in [0.5,3] is paper-units; in our units the same
window maps to L ~ [2,25] via the unit map below — bracket [1,40] used.)
Outputs c = 2L, kappa = L/I*, and the convention-free invariants
lam*, g* = I(lam*,d*)/I(1,d*), V* = (d*/d0)^3 (d0 = static lam=1 optimum),
E_rot/E_tot.

UNIT MAP to the paper's (e0, i0) convention — ZERO calibration, fixed by the
exact eps=0 compacton f0 = 2 arccos(r/R*), R* = 2^(5/6), which both sides own:
  ours:  E6 = E0 = 16 sqrt2/15 (=(1/2) BPS bound), I0 = 128 pi R*^3/105
  paper: e0 = 64/(15 pi),  i0 = 256/(105 pi)
  =>  u_E = e0 / (2 E6)      = sqrt2/pi        (energy unit)
      u_I = i0 / I0          = 2^(-3/2)/pi^2   (inertia unit)
      c_paper     = c_ours * sqrt(u_E u_I) = c_ours / sqrt(2 pi^3)
      kappa_paper = kappa_ours * sqrt(u_E/u_I) = 2 sqrt(pi) kappa_ours
  Analytic identity check: the family-A t=0 closure gives, in closed form,
  s := L^2/(2I) = E6, V* = sqrt2, c_paper = 2 sqrt(e0 i0) = 2.0533,
  kappa_paper = 0.9354 — exactly the reduced model's g=1 endpoint. (The
  reduced model's quoted eps->0 endpoint c0 = 2.5147, kappa = sqrt(7/12)
  is its g = 3/2, lam -> 0 pancake endpoint.)

VALIDATION GATES (V1..V3 hard, abort on failure; V4 is the physics
comparison the task flags as informative):
  V1 sphere sectors vs step-1 radial integrals (< 0.5%), degree B = +1
  V2 sphere inertia vs (16 pi/3) INT r^2 sin^2 f dr (< 0.5%)
  V3 lam = 0.6 INDEPENDENT lab-grid run (direct central differences of the
     deformed field on a stretched (rho,z) grid): E6+E0 unchanged (< 0.2%),
     plus cross-check of the family-A closed forms for E2, E4, I
  V4 I(0.42)/I(1.0): family A -> 1.0000 exactly (see FINDING);
     family B -> g_oblate(0.42) = 1.2855 (paper's curve, few %)

EPS-SCAN: t re-dialled (step-1 machinery, warm-started fixed point) so the
achieved radial ratio hits eps in {0.02, 0.035, 0.05, 0.075, 0.1}; closure
re-run per point per family; deficit exponent of c(eps) fitted against the
t = 0 endpoint c0 of the same family (compacton profile; family-A c0 also
known in closed form — quadrature cross-check).

Deterministic. No RNG. Runtime ~6-10 min at N = 512 (quarter-plane 512x512
~ half-plane 512x1024, inside the task's suggested budget); one coarser run
(N = 336) repeats the gates and the eps = 0.05 closures for convergence.

Outputs: axi_results.json, axi_landscape.png, axi_epsscan.png (run log via
tee -> axi_run.log; axi_RESULTS.md written from the printed tables).
"""

import json
import math
import sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier2-closure"
sys.path.insert(0, OUTDIR)
import radial_solve as rs           # step-1 machinery (grid, solver, thomas)
import closure_reduced as cr        # only for g_oblate comparison in V4

# ----------------------------------------------------------------------
# frozen conventions / constants (NOTHING here is fitted)
# ----------------------------------------------------------------------
TFROZEN = 0.008276434949296802      # a2 = a4 from step 1 (frozen)
RSTAR = 2.0 ** (5.0 / 6.0)
FOURPI = 4.0 * np.pi
C6 = np.pi ** 3                     # E6 = pi^3 INT b^2 dV (derived, see header)
E6_COMPACTON = 16.0 * math.sqrt(2.0) / 15.0          # = E0_compacton
I_COMPACTON = 128.0 * math.pi * RSTAR ** 3 / 105.0   # exact
U_E = math.sqrt(2.0) / math.pi                        # paper energy unit
U_I = 2.0 ** (-1.5) / math.pi ** 2                    # paper inertia unit
CMAP_L = math.sqrt(U_E * U_I)       # L_paper = CMAP_L * L_ours = 1/sqrt(2 pi^3)
KMAP = math.sqrt(U_E / U_I)         # kappa_paper = KMAP * kappa_ours = 2 sqrt(pi)

EPS_LIST = [0.02, 0.035, 0.05, 0.075, 0.1]
QUICK = "--quick" in sys.argv


# ----------------------------------------------------------------------
# natural cubic spline (numpy-only; f = 0 beyond the last node)
# ----------------------------------------------------------------------
class CubicSpline1D:
    def __init__(self, x, y, right_zero=True):
        self.right_zero = right_zero
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        h = np.diff(x)
        diag = 2.0 * (h[:-1] + h[1:])
        off = h[1:-1]
        rhs = 6.0 * ((y[2:] - y[1:-1]) / h[1:] - (y[1:-1] - y[:-2]) / h[:-1])
        M = np.zeros(len(x))
        M[1:-1] = rs.thomas(diag, off, rhs)   # step-1 Thomas solver reused
        self.x, self.y, self.h, self.M = x, y, h, M

    def __call__(self, xq):
        xq = np.asarray(xq, float)
        i = np.clip(np.searchsorted(self.x, xq) - 1, 0, len(self.x) - 2)
        hh = self.h[i]
        tt = (xq - self.x[i]) / hh
        v = (self.y[i] * (1.0 - tt) + self.y[i + 1] * tt
             + ((1.0 - tt) ** 3 - (1.0 - tt)) * self.M[i] * hh * hh / 6.0
             + (tt ** 3 - tt) * self.M[i + 1] * hh * hh / 6.0)
        if self.right_zero:              # radial profile: f = 0 beyond support
            v = np.where(xq >= self.x[-1], 0.0, v)
        v = np.where(xq <= self.x[0], self.y[0], v)
        return v


# ----------------------------------------------------------------------
# step-1 radial machinery wrappers
# ----------------------------------------------------------------------
def radial_profile(t, f_init=None, r=None, N=4000, rmax=6.0):
    if r is None:
        r = rs.make_grid(N, rmax)
    f0 = rs.init_profile(r) if f_init is None else f_init
    f, E, sect, g, gmax, dE, it = rs.solve_profile(r, f0, t)
    return r, f, sect


def dial_eps(eps_target, t_guess, rtol=1.0e-3):
    """Fixed-point on t (ratio ~ linear in t), warm-started profiles."""
    t, r, f = t_guess, None, None
    ratio, sect = None, None
    for _ in range(12):
        r, f, sect = radial_profile(t, f_init=f, r=r)
        E2, E4, E6, E0 = sect
        ratio = t * (E2 + E4) / (E6 + E0)
        if abs(ratio / eps_target - 1.0) < rtol:
            break
        t *= eps_target / ratio
    return t, r, f, sect, ratio


def radial_inertia(r, f):
    rm = 0.5 * (r[1:] + r[:-1])
    fm = 0.5 * (f[1:] + f[:-1])
    return float((16.0 * np.pi / 3.0)
                 * np.sum(np.diff(r) * rm * rm * np.sin(fm) ** 2))


def f_compacton(rq):
    return 2.0 * np.arccos(np.clip(np.asarray(rq, float) / RSTAR, 0.0, 1.0))


# ----------------------------------------------------------------------
# 2-D axisymmetric engine (phi = 0 plane, quarter-plane, cell-centred)
# ----------------------------------------------------------------------
def make_famA(fint, lam, d):
    """Task's literal deformation: full pullback q_base(x/(ed), y/(ed), z/(pd))."""
    e = lam ** (-1.0 / 3.0)
    p = lam ** (2.0 / 3.0)

    def qfun(R, Z):
        X = R / (e * d)
        Y = Z / (p * d)
        rp = np.hypot(X, Y)
        F = fint(rp)
        s = np.sin(F)
        return np.cos(F), s * X / rp, s * Y / rp, s / (rp * e * d)

    return qfun


def make_famB(fint, lam, d):
    """Paper's spheroidal-shell ansatz: f on spheroids, iso-direction = lab xhat."""
    e = lam ** (-1.0 / 3.0)
    p = lam ** (2.0 / 3.0)

    def qfun(R, Z):
        m = np.hypot(R / (e * d), Z / (p * d))
        rl = np.hypot(R, Z)
        F = fint(m)
        s = np.sin(F)
        return np.cos(F), s * R / rl, s * Z / rl, s / rl

    return qfun


def sector_pieces(qfun, N, rho_max, z_max, Nz=None):
    """All 2-D sector integrals for an axisymmetric winding-1 field.

    Cell-centred quarter-plane grid (z >= 0 doubled; every density is even in
    z since P(rho,-z) = sigma1 P(rho,z) sigma1 for both families), one ghost
    ring on each side evaluated at its physical (possibly negative) coordinate
    — the builders are globally defined and carry the reflection parities
    exactly, so central differences are valid up to the boundary. L_rho, L_z
    from central differences of P; L_phi analytic (see module docstring)."""
    if Nz is None:
        Nz = N
    hr = rho_max / N
    hz = z_max / Nz
    rho = (np.arange(N + 2) - 0.5) * hr
    zz = (np.arange(Nz + 2) - 0.5) * hz
    R, Z = np.meshgrid(rho, zz, indexing="ij")
    q0, q1, q3, q1r = qfun(R, Z)
    S = np.s_[1:-1]
    q0i, q1i, q3i, q1ri = q0[S, S], q1[S, S], q3[S, S], q1r[S, S]

    def dr(A):
        return (A[2:, 1:-1] - A[:-2, 1:-1]) / (2.0 * hr)

    def dz(A):
        return (A[1:-1, 2:] - A[1:-1, :-2]) / (2.0 * hz)

    def current(p0, p1, p3):        # L = i a.sigma from P^dag dP, q2 = 0
        return (q0i * p1 - p0 * q1i,
                q3i * p1 - q1i * p3,
                q0i * p3 - p0 * q3i)

    ar = current(dr(q0), dr(q1), dr(q3))
    az = current(dz(q0), dz(q1), dz(q3))
    ap = (q1ri * q3i, -q1ri * q0i, -q1ri * q1i)      # a_phi / rho, analytic

    def dot(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

    def cr2(a, b):                   # |a x b|^2
        return dot(a, a) * dot(b, b) - dot(a, b) ** 2

    trip = (ar[0] * (az[1] * ap[2] - az[2] * ap[1])
            + ar[1] * (az[2] * ap[0] - az[0] * ap[2])
            + ar[2] * (az[0] * ap[1] - az[1] * ap[0]))
    b = -trip / (2.0 * np.pi ** 2)

    W = 2.0 * 2.0 * np.pi * R[S, S] * hr * hz        # x2 for the z<0 mirror

    def integ(D):
        return float(np.sum(D * W))

    P = dict(S11=integ(dot(ar, ar)), S22=integ(dot(az, az)),
             S33=integ(dot(ap, ap)),
             C12=integ(cr2(ar, az)), C13=integ(cr2(ar, ap)),
             C23=integ(cr2(az, ap)),
             BQ=integ(b * b), BDEG=integ(b),
             POT=integ(1.0 - q0i), INR=integ(2.0 * q1i * q1i))
    P["E2"] = (P["S11"] + P["S22"] + P["S33"]) / FOURPI
    P["E4"] = (P["C12"] + P["C13"] + P["C23"]) / FOURPI
    P["E6"] = C6 * P["BQ"]
    P["E0"] = P["POT"] / FOURPI
    P["I"] = P["INR"]
    return P


def support_radius(r, f, thresh=1.0e-10, margin=1.05):
    idx = np.nonzero(f > thresh)[0]
    return float(min(r[idx[-1]] * margin, r[-1]))


def run_lab(fint, family, lam, d, N, r_sup, iso=False):
    """One direct lab-grid evaluation.

    Extents follow the spheroid. iso=False: N x N cells stretched with the
    deformation (resolution-matched; used for the family-B lambda scans).
    iso=True: SQUARE cells h = rho_max/N (z-cell count adjusted) — for family A
    this sampling is INCOMMENSURATE with the pullback of the base grid, so the
    gate-V3/V4 checks are genuinely independent of the closed-form scalings
    (a commensurate grid would reproduce them to machine precision by
    construction, testing nothing)."""
    e = lam ** (-1.0 / 3.0)
    p = lam ** (2.0 / 3.0)
    mk = make_famA if family == "A" else make_famB
    rho_max, z_max = e * d * r_sup, p * d * r_sup
    if not iso:
        return sector_pieces(mk(fint, lam, d), N, rho_max, z_max)
    h = rho_max / N
    Nz = int(math.ceil(z_max / h))
    return sector_pieces(mk(fint, lam, d), N, rho_max, Nz * h, Nz=Nz)


# ----------------------------------------------------------------------
# closure: family functions at d = 1, exact d-scaling, minimise, clock
# ----------------------------------------------------------------------
def famA_funcs(P, lmin=0.15, lmax=1.60):
    T, S2 = P["S11"] + P["S33"], P["S22"]
    CA, CB = P["C12"] + P["C23"], P["C13"]
    e6c, e0c, iic = P["E6"], P["E0"], P["I"]
    return dict(
        e2=lambda lam: (T * lam ** (2.0 / 3.0) + S2 * lam ** (-4.0 / 3.0)) / FOURPI,
        e4=lambda lam: (CA * lam ** (-2.0 / 3.0) + CB * lam ** (4.0 / 3.0)) / FOURPI,
        e6=lambda lam: e6c, e0=lambda lam: e0c, iI=lambda lam: iic,
        lmin=lmin, lmax=lmax)


def famB_funcs(fint, r_sup, N, lam_grid):
    rows = {k: [] for k in ("E2", "E4", "E6", "E0", "I")}
    for lam in lam_grid:
        P = run_lab(fint, "B", lam, 1.0, N, r_sup)
        for k in rows:
            rows[k].append(P[k])
    spl = {k: CubicSpline1D(np.asarray(lam_grid), np.asarray(v),
                            right_zero=False)
           for k, v in rows.items()}
    F = dict(e2=lambda lam: float(spl["E2"](lam)),
             e4=lambda lam: float(spl["E4"](lam)),
             e6=lambda lam: float(spl["E6"](lam)),
             e0=lambda lam: float(spl["E0"](lam)),
             iI=lambda lam: float(spl["I"](lam)),
             lmin=float(lam_grid[0]), lmax=float(lam_grid[-1]))
    return F, rows


def golden(fun, a, b, iters=60):
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, fd = fun(c), fun(d)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = fun(d)
    x = 0.5 * (a + b)
    return x, fun(x)


def R_of(F, t, lam, d, L):
    e2 = F["e2"](lam) * d
    e4 = F["e4"](lam) / d
    e6 = F["e6"](lam) / d ** 3
    e0 = F["e0"](lam) * d ** 3
    iI = F["iI"](lam) * d ** 3
    Estat = t * (e2 + e4) + e6 + e0
    return Estat + (L * L / (2.0 * iI) if L else 0.0), (e2, e4, e6, e0, iI, Estat)


def min_over_d(F, t, lam, L, dlo=0.4, dhi=2.2):
    u, Ru = golden(lambda w: R_of(F, t, lam, math.exp(w), L)[0],
                   math.log(dlo), math.log(dhi), iters=48)
    return math.exp(u), Ru


def min_R(F, t, L, nlam=41):
    lg = np.linspace(F["lmin"], F["lmax"], nlam)
    Rs = [min_over_d(F, t, lam, L)[1] for lam in lg]
    i = int(np.argmin(Rs))
    a, b = lg[max(i - 1, 0)], lg[min(i + 1, nlam - 1)]
    lam, _ = golden(lambda l: min_over_d(F, t, l, L)[1], a, b, iters=44)
    d, R = min_over_d(F, t, lam, L)
    return lam, d, R


def clock_solve(F, t, Llo=1.0, Lhi=40.0, iters=48):
    def h(L):
        lam, d, R = min_R(F, t, L)
        iI = F["iI"](lam) * d ** 3
        return 2.0 * L * L / iI - R, (lam, d, R, iI)

    hlo, _ = h(Llo)
    hhi, _ = h(Lhi)
    assert hlo < 0.0 < hhi, "clock bracket failed: %.3g %.3g" % (hlo, hhi)
    for _ in range(iters):
        Lm = 0.5 * (Llo + Lhi)
        hm, _ = h(Lm)
        if hlo * hm <= 0.0:
            Lhi = Lm
        else:
            Llo, hlo = Lm, hm
    L = 0.5 * (Llo + Lhi)
    hm, (lam, d, R, iI) = h(L)
    return L, lam, d, R, iI, hm


def closure_report(F, t, name):
    L, lam, d, R, iI, resid = clock_solve(F, t)
    _, (e2, e4, e6, e0, iI2, Estat) = R_of(F, t, lam, d, L)
    d0, _ = min_over_d(F, t, 1.0, 0.0)          # static (L = 0) optimum, lam = 1
    Erot = L * L / (2.0 * iI)
    g = F["iI"](lam) / F["iI"](1.0)             # same d: pure shape enhancement
    return dict(name=name, t=t, L=L, c_ours=2.0 * L, kappa_ours=L / iI,
                c_paper=2.0 * L * CMAP_L, kappa_paper=(L / iI) * KMAP,
                lam=lam, d=d, d0=d0, V=(d / d0) ** 3, g=g, R=R, Estat=Estat,
                Erot=Erot, erot_frac=Erot / R,
                eps_at_solution=t * (e2 + e4) / (e6 + e0),
                clock_resid=resid, E2=e2, E4=e4, E6=e6, E0=e0, I=iI)


# ----------------------------------------------------------------------
# pipeline at one resolution (gates + eps = 0.05 closures)
# ----------------------------------------------------------------------
def pipeline(N, fspl, r, f, sectN, targets, lam_grid, verbose=True, hard=True):
    out = {}
    r_sup = support_radius(r, f)

    # ---- gate V1: sphere sectors vs step-1 radial integrals ------------
    P1 = sector_pieces(make_famA(fspl, 1.0, 1.0), N, r_sup, r_sup)
    v1 = {k: P1[k] / targets[k] - 1.0 for k in ("E2", "E4", "E6", "E0")}
    v1_pass = max(abs(x) for x in v1.values()) < 5.0e-3
    # ---- gate V2: inertia vs radial formula ----------------------------
    Irad = radial_inertia(r, f)
    v2 = P1["I"] / Irad - 1.0
    v2_pass = abs(v2) < 5.0e-3

    # ---- gate V3: independent lab run at lam = 0.6 (isotropic cells,
    # incommensurate with the pullback -> a real differentiation test) ----
    P06 = run_lab(fspl, "A", 0.6, 1.0, N, r_sup, iso=True)
    bps = (P06["E6"] + P06["E0"]) / (P1["E6"] + P1["E0"]) - 1.0
    v3_pass = abs(bps) < 2.0e-3
    FA = famA_funcs(P1)
    xchk = dict(E2=P06["E2"] / FA["e2"](0.6) - 1.0,
                E4=P06["E4"] / FA["e4"](0.6) - 1.0,
                E6=P06["E6"] / P1["E6"] - 1.0,
                E0=P06["E0"] / P1["E0"] - 1.0,
                I=P06["I"] / P1["I"] - 1.0)

    # ---- gate V4: inertia curve at lam = 0.42 --------------------------
    P42A = run_lab(fspl, "A", 0.42, 1.0, N, r_sup, iso=True)
    P42B = run_lab(fspl, "B", 0.42, 1.0, N, r_sup)
    gtarget = cr.g_oblate(0.42) / cr.g_oblate(1.0)
    v4 = dict(gA=P42A["I"] / P1["I"], gB=P42B["I"] / P1["I"], g_paper=gtarget)

    out["gates"] = dict(V1=v1, V1_pass=v1_pass, B_degree=P1["BDEG"],
                        V2=v2, V2_pass=v2_pass, I2D=P1["I"], Irad=Irad,
                        V3_bps=bps, V3_pass=v3_pass, V3_cross=xchk,
                        V4=v4, r_sup=r_sup)
    out["P1"] = P1
    if verbose:
        print("  V1 sphere sectors (2D vs step-1 radial), rel. err:")
        for k in ("E2", "E4", "E6", "E0"):
            print("     %s: %+.3e  (2D %.6f | radial %.6f)"
                  % (k, v1[k], P1[k], targets[k]))
        print("     degree B (2D) = %+.6f ;  V1 %s"
              % (P1["BDEG"], "PASS" if v1_pass else "FAIL"))
        print("  V2 inertia: I_2D = %.5f vs radial %.5f, rel %+.3e  %s"
              % (P1["I"], Irad, v2, "PASS" if v2_pass else "FAIL"))
        print("  V3 lab run lam=0.6: (E6+E0) rel shift %+.3e  %s"
              % (bps, "PASS" if v3_pass else "FAIL"))
        print("     cross-check lab vs closed-form A: "
              + "  ".join("%s %+..1e".replace("..", ".") % (k, v)
                          for k, v in xchk.items()))
        print("  V4 I(0.42)/I(1): family A = %.5f (SDiff pullback; exact 1)"
              % v4["gA"])
        print("     family B = %.5f vs paper g_oblate = %.5f (rel %+.3e)"
              % (v4["gB"], gtarget, v4["gB"] / gtarget - 1.0))
    if not (v1_pass and v2_pass and v3_pass):
        if hard:
            print("!! HARD GATE FAILURE — aborting before closure.")
            with open(OUTDIR + "/axi_results.json", "w") as fh:
                json.dump(_jsonable(out), fh, indent=2)
            sys.exit(1)
        print("   (coarse-grid gate tolerance exceeded — expected at reduced "
              "N; convergence deltas below quantify it)")

    # ---- closures at eps = 0.05 (t frozen) ------------------------------
    FB, rowsB = famB_funcs(fspl, r_sup, N, lam_grid)
    out["solA"] = closure_report(FA, TFROZEN, "A(SDiff)")
    out["solB"] = closure_report(FB, TFROZEN, "B(shell)")
    out["FA"], out["FB"], out["rowsB"] = FA, FB, rowsB
    return out


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()
                if not callable(v) and k not in ("FA", "FB", "P1")}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    return obj


def print_sol(s, extra=""):
    print("   %-9s L=%8.4f  c_our=%8.4f c_pap=%6.4f | kap_our=%7.5f "
          "kap_pap=%6.4f | lam*=%6.4f g*=%6.4f V*=%6.4f Erot/E=%8.6f | "
          "d*=%6.4f d0=%6.4f eps@sol=%7.5f%s"
          % (s["name"], s["L"], s["c_ours"], s["c_paper"], s["kappa_ours"],
             s["kappa_paper"], s["lam"], s["g"], s["V"], s["erot_frac"],
             s["d"], s["d0"], s["eps_at_solution"], extra))


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------
INK, GRID_C = "#333333", "#dddddd"
BLUE, AQUA, RED, VIOLET = "#2a78d6", "#1baf7a", "#e34948", "#4a3aa7"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
       "#1c5cab", "#104281", "#0d366b"]
CMAP_BLUES = LinearSegmentedColormap.from_list("seqblue", SEQ[::-1])  # dark=low


def style_ax(ax):
    ax.set_facecolor("white")
    ax.grid(color=GRID_C, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(colors=INK, labelsize=9)


def fig_landscape(res):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), dpi=160)
    fig.patch.set_facecolor("white")
    for ax, (F, sol, ttl) in zip(axes, [
            (res["FA"], res["solA"],
             "family A — task's SDiff pullback (I flat in $\\lambda$)"),
            (res["FB"], res["solB"],
             "family B — paper's spheroidal-shell ansatz")]):
        lam_v = np.linspace(F["lmin"], min(F["lmax"], 1.45), 90)
        d_v = np.linspace(0.85, 1.55, 80)
        L = sol["L"]
        RR = np.array([[R_of(F, TFROZEN, lam, dd, L)[0] for lam in lam_v]
                       for dd in d_v])
        lo = RR.min()
        lev = lo + (RR.max() - lo) * np.linspace(0.0, 1.0, 18) ** 1.6
        cf = ax.contourf(lam_v, d_v, RR, levels=lev, cmap=CMAP_BLUES)
        ax.contour(lam_v, d_v, RR, levels=lev[1:12], colors="white",
                   linewidths=0.5, alpha=0.55)
        ax.plot(sol["lam"], sol["d"], marker="*", ms=17, color=RED,
                mec="white", mew=1.2, ls="none", zorder=5)
        ax.annotate("$(\\lambda^*, d^*) = (%.3f, %.3f)$\n$R^*=%.4f$"
                    % (sol["lam"], sol["d"], sol["R"]),
                    (sol["lam"], sol["d"]), xytext=(-10, 12), ha="right",
                    textcoords="offset points", fontsize=9, color="white")
        ax.axhline(sol["d0"], color=RED, lw=1.0, ls=":", alpha=0.8)
        ax.text(lam_v[0] + 0.02, sol["d0"] + 0.008, "$d_0$ (static)",
                color=RED, fontsize=8, va="bottom")
        cb = fig.colorbar(cf, ax=ax, shrink=0.92)
        cb.set_label("Routhian $R$ (darker = lower)", fontsize=8, color=INK)
        cb.ax.tick_params(labelsize=7, colors=INK)
        ax.set_xlabel("oblateness $\\lambda$ (polar/equatorial)", color=INK)
        ax.set_ylabel("dilation $d$", color=INK)
        ax.set_title(ttl + "  [$L^*=%.3f$]" % L, fontsize=10, color=INK)
        ax.tick_params(colors=INK, labelsize=9)
    fig.suptitle("Tier-2b step 2: Routhian landscapes at the solved spin, "
                 "$\\epsilon = 0.05$ ($t$ frozen from step 1)",
                 fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUTDIR + "/axi_landscape.png", facecolor="white")
    plt.close(fig)


def fig_epsscan(scan, endA, endB, fits):
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6), dpi=160)
    fig.patch.set_facecolor("white")
    eps = np.array([row["eps"] for row in scan])
    cA = np.array([row["A"]["c_paper"] for row in scan])
    cB = np.array([row["B"]["c_paper"] for row in scan])
    c0A, c0B = endA["c_paper"], endB["c_paper"]

    ax = axes[0, 0]
    style_ax(ax)
    ax.plot(eps, cA, "-o", color=BLUE, lw=2, ms=6, mec="white")
    ax.plot(eps, cB, "-s", color=AQUA, lw=2, ms=6, mec="white")
    ax.axhline(c0A, color=BLUE, lw=1.2, ls=":")
    ax.axhline(c0B, color=AQUA, lw=1.2, ls=":")
    ax.axhline(2.5147, color=VIOLET, lw=1.0, ls="--")
    ax.errorbar([0.05], [2.37], yerr=[0.09], fmt="D", color=RED, ms=6,
                capsize=4, mec="white")
    ax.text(0.052, 2.37, " paper $\\langle r_1\\rangle$\n $2.37\\pm0.09$",
            color=RED, fontsize=8, va="center")
    ax.text(eps[-1], c0A, "  $c_0$ A", color=BLUE, fontsize=8, va="bottom")
    ax.text(eps[-1], c0B, "  $c_0$ B", color=AQUA, fontsize=8, va="bottom")
    ax.text(eps[0], 2.5147, " paper $\\epsilon\\to0$ endpoint (g=3/2): 2.5147",
            color=VIOLET, fontsize=8, va="bottom")
    ax.annotate("A (SDiff pullback)", (eps[3], cA[3]), xytext=(6, -16),
                textcoords="offset points", color=BLUE, fontsize=9)
    ax.annotate("B (spheroidal shell)", (eps[0], cB[0]), xytext=(-2, 10),
                textcoords="offset points", color=AQUA, fontsize=9)
    ax.set_xlabel("$\\epsilon$ (achieved sector ratio)", color=INK)
    ax.set_ylabel("$c = 2L$  (paper units via compacton map)", color=INK)
    ax.set_title("(a) clock charge vs $\\epsilon$", fontsize=10, color=INK)

    ax = axes[0, 1]
    style_ax(ax)
    for c, c0, col, lab in ((cA, c0A, BLUE, "A"), (cB, c0B, AQUA, "B")):
        dfc = np.abs(1.0 - c / c0)
        nu, amp, sgn = fits[lab]
        xe = np.array([eps[0] * 0.9, eps[-1] * 1.1])
        if lab == "A":
            ax.loglog(eps, dfc, "o", color=col, ms=7, mec="white")
            ax.loglog(xe, amp * xe ** nu, "-", color=col, lw=1.6, alpha=0.8)
        else:   # A and B deficits nearly coincide: open squares, dashed fit
            ax.loglog(eps, dfc, "s", color=col, ms=8, mfc="none", mew=1.6)
            ax.loglog(xe, amp * xe ** nu, "--", color=col, lw=1.6, alpha=0.9)
        ax.text(0.03, 0.93 if lab == "A" else 0.85,
                "%s: $|1-c/c_0| \\propto \\epsilon^{%.3f}$ (%s)"
                % (lab, nu, "$c>c_0$" if sgn < 0 else "$c<c_0$"),
                transform=ax.transAxes, color=col, fontsize=9)
    ax.text(0.03, 0.78, "(A and B nearly coincide)", transform=ax.transAxes,
            color=INK, fontsize=8)
    xe = np.array([0.015, 0.11])
    ax.loglog(xe, 0.42 * xe ** (2.0 / 3.0), "--", color=VIOLET, lw=1.2)
    ax.text(xe[-1], 0.42 * xe[-1] ** (2. / 3.), "corpus $0.42\\,\\epsilon^{2/3}$ ",
            color=VIOLET, fontsize=8, va="bottom", ha="right")
    ax.set_xlabel("$\\epsilon$", color=INK)
    ax.set_ylabel("$|1 - c/c_0|$", color=INK)
    ax.set_title("(b) deficit law and fitted exponents", fontsize=10, color=INK)

    ax = axes[1, 0]
    style_ax(ax)
    lA = [row["A"]["lam"] for row in scan]
    lB = [row["B"]["lam"] for row in scan]
    ax.plot(eps, lA, "-o", color=BLUE, lw=2, ms=6, mec="white")
    ax.plot(eps, lB, "-s", color=AQUA, lw=2, ms=6, mec="white")
    ax.axhspan(0.37, 0.47, color=RED, alpha=0.10, lw=0)
    ax.axhline(0.42, color=RED, lw=1.0, ls="--")
    ax.text(eps[0], 0.425, " paper $\\lambda^* = 0.42\\pm0.05$", color=RED,
            fontsize=8)
    ax.annotate("A", (eps[-2], lA[-2]), xytext=(0, 8),
                textcoords="offset points", color=BLUE, fontsize=10)
    ax.annotate("B", (eps[-2], lB[-2]), xytext=(0, 8),
                textcoords="offset points", color=AQUA, fontsize=10)
    ax.set_xlabel("$\\epsilon$", color=INK)
    ax.set_ylabel("$\\lambda^*$", color=INK)
    ax.set_ylim(0.0, 1.15)
    ax.set_title("(c) optimal oblateness", fontsize=10, color=INK)

    ax = axes[1, 1]
    style_ax(ax)
    gA = [row["A"]["g"] for row in scan]
    gB = [row["B"]["g"] for row in scan]
    ax.plot(eps, gA, "-o", color=BLUE, lw=2, ms=6, mec="white")
    ax.plot(eps, gB, "-s", color=AQUA, lw=2, ms=6, mec="white")
    ax.axhspan(1.27, 1.35, color=RED, alpha=0.10, lw=0)
    ax.axhline(1.31, color=RED, lw=1.0, ls="--")
    ax.text(eps[0], 1.315, " paper $g^* = 1.31\\pm0.04$", color=RED, fontsize=8)
    ax.annotate("A (exactly 1: SDiff-invariant inertia)", (eps[1], gA[1]),
                xytext=(0, 7), textcoords="offset points", color=BLUE,
                fontsize=9)
    ax.annotate("B", (eps[-2], gB[-2]), xytext=(0, 8),
                textcoords="offset points", color=AQUA, fontsize=10)
    ax.set_ylim(0.96, 1.38)
    ax.set_xlabel("$\\epsilon$", color=INK)
    ax.set_ylabel("$g^* = I(\\lambda^*,d^*)/I(1,d^*)$", color=INK)
    ax.set_title("(d) shape enhancement of the inertia", fontsize=10, color=INK)

    fig.suptitle("Tier-2b step 2: $\\epsilon$-scan of the field-level closure "
                 "(zero calibrated constants)", fontsize=11, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTDIR + "/axi_epsscan.png", facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    N = 160 if QUICK else 512
    NC = 112 if QUICK else 336
    lam_grid = np.arange(0.15, 1.1001, 0.05 if QUICK else 0.025)

    with open(OUTDIR + "/radial_results.json") as fh:
        RJ = json.load(fh)
    targets = {k: RJ["N_run"][k] for k in ("E2", "E4", "E6", "E0")}

    print("== step-1 radial re-solve at frozen t = %.15g ==" % TFROZEN)
    r, f, sect = radial_profile(TFROZEN)
    print("   sectors vs radial_results.json (should be ~1e-9): "
          + "  ".join("%s %+.1e" % (k, sect[i] / targets[k] - 1.0)
                      for i, k in enumerate(("E2", "E4", "E6", "E0"))))
    fspl = CubicSpline1D(r, f)

    print("\n== VALIDATION GATES + eps=0.05 closures  (N = %d) ==" % N)
    res = pipeline(N, fspl, r, f, sect, targets, lam_grid)
    print("\n   closures at eps = 0.05 (t frozen):")
    print_sol(res["solA"])
    print_sol(res["solB"])

    # ---- t = 0 endpoints (compacton profile) ---------------------------
    print("\n== t = 0 endpoints (exact eps=0 compacton profile) ==")
    r_supc = RSTAR * 1.02
    P1c = sector_pieces(make_famA(f_compacton, 1.0, 1.0), N, r_supc, r_supc)
    print("   compacton engine check: E6 %+.2e  E0 %+.2e  I %+.2e  (vs exact)"
          % (P1c["E6"] / E6_COMPACTON - 1.0, P1c["E0"] / E6_COMPACTON - 1.0,
             P1c["I"] / I_COMPACTON - 1.0))
    FA0 = famA_funcs(P1c)
    FB0, _ = famB_funcs(f_compacton, r_supc, N, lam_grid)
    endA = closure_report(FA0, 0.0, "A0")
    endB = closure_report(FB0, 0.0, "B0")
    # closed-form family-A t=0 solution (uses exact compacton values)
    L0 = math.sqrt(2.0 * I_COMPACTON * E6_COMPACTON)
    exact = dict(c_ours=2.0 * L0, c_paper=2.0 * L0 * CMAP_L,
                 kappa_paper=KMAP * L0 / (I_COMPACTON * math.sqrt(2.0)),
                 V=math.sqrt(2.0))
    print_sol(endA, "   [analytic: c_pap %.4f kap_pap %.4f V %.6f]"
              % (exact["c_paper"], exact["kappa_paper"], exact["V"]))
    print("   NOTE: family-A t=0 lambda direction is EXACTLY FLAT (E6, E0, I "
          "all SDiff-invariant); the reported lam* is an arbitrary point of "
          "the degenerate valley, NOT a physical optimum. The task expected "
          "a runaway to extreme oblate here; at field level there is no "
          "gradient in lambda at all.")
    print_sol(endB)
    if endB["lam"] - FB0["lmin"] < 0.02:
        print("   NOTE: family-B t=0 optimum at the lambda-grid edge %.2f "
              "(family cut-off)" % FB0["lmin"])

    # ---- eps-scan --------------------------------------------------------
    print("\n== eps-scan (t re-dialled per point) ==")
    scan = []
    eps_list = [0.02, 0.05, 0.1] if QUICK else EPS_LIST
    for eps in eps_list:
        if abs(eps - 0.05) < 1e-12:
            te, re_, fe, secte, ratio = (TFROZEN, r, f, sect,
                                         TFROZEN * (sect[0] + sect[1])
                                         / (sect[2] + sect[3]))
            solA, solB = res["solA"], res["solB"]
        else:
            te, re_, fe, secte, ratio = dial_eps(eps, TFROZEN * eps / 0.05)
            fspl_e = CubicSpline1D(re_, fe)
            r_sup_e = support_radius(re_, fe)
            P1e = sector_pieces(make_famA(fspl_e, 1.0, 1.0), N,
                                r_sup_e, r_sup_e)
            FAe = famA_funcs(P1e)
            FBe, _ = famB_funcs(fspl_e, r_sup_e, N, lam_grid)
            solA = closure_report(FAe, te, "A")
            solB = closure_report(FBe, te, "B")
        print("   eps %.4f (achieved %.5f, t = %.7g):" % (eps, ratio, te))
        print_sol(solA)
        print_sol(solB)
        scan.append(dict(eps=float(ratio), t=float(te), A=solA, B=solB))

    # deficit exponent fits vs same-family t=0 endpoint
    fits = {}
    for lab, end in (("A", endA), ("B", endB)):
        e_arr = np.array([row["eps"] for row in scan])
        c_arr = np.array([row[lab]["c_ours"] for row in scan])
        dfc = 1.0 - c_arr / end["c_ours"]
        sgn = float(np.sign(dfc.mean()))
        nu, lnA = np.polyfit(np.log(e_arr), np.log(np.abs(dfc)), 1)
        fits[lab] = (float(nu), float(math.exp(lnA)), sgn)
        print("   family %s deficit fit: 1 - c/c0 = %+.4f eps^%.4f  "
              "(deficit sign %s; c0_ours = %.4f)"
              % (lab, sgn * math.exp(lnA), nu,
                 "negative (c > c0)" if sgn < 0 else "positive", end["c_ours"]))

    # ---- convergence: coarser rerun of gates + eps = 0.05 closures ------
    print("\n== convergence: coarse rerun (N = %d vs %d) ==" % (NC, N))
    resC = pipeline(NC, fspl, r, f, sect, targets, lam_grid, verbose=False,
                    hard=False)
    conv = {}
    for fam in ("solA", "solB"):
        for k in ("c_ours", "kappa_ours", "lam", "g", "V", "R"):
            key = "%s.%s" % (fam, k)
            a, b = res[fam][k], resC[fam][k]
            conv[key] = dict(fine=a, coarse=b,
                             rel=abs(b - a) / (abs(a) if a else 1.0))
            print("   %-16s fine %12.6f  coarse %12.6f  rel delta %.2e"
                  % (key, a, b, conv[key]["rel"]))
    gate_conv = dict(
        V1_maxerr_fine=max(abs(x) for x in res["gates"]["V1"].values()),
        V1_maxerr_coarse=max(abs(x) for x in resC["gates"]["V1"].values()),
        V3_fine=res["gates"]["V3_bps"], V3_coarse=resC["gates"]["V3_bps"],
        gB42_fine=res["gates"]["V4"]["gB"], gB42_coarse=resC["gates"]["V4"]["gB"])
    print("   gates: V1 max err %.2e -> %.2e | V3 %.2e -> %.2e | "
          "gB(0.42) %.5f -> %.5f"
          % (gate_conv["V1_maxerr_coarse"], gate_conv["V1_maxerr_fine"],
             gate_conv["V3_coarse"], gate_conv["V3_fine"],
             gate_conv["gB42_coarse"], gate_conv["gB42_fine"]))

    # ---- figures + json --------------------------------------------------
    fig_landscape(res)
    fig_epsscan(scan, endA, endB, fits)
    out = dict(N=N, N_coarse=NC, t_frozen=TFROZEN,
               unit_map=dict(u_E=U_E, u_I=U_I, c_map=CMAP_L, kappa_map=KMAP),
               gates=res["gates"], solA=res["solA"], solB=res["solB"],
               endA=endA, endB=endB, endA_analytic=exact,
               scan=scan, fits=fits, convergence=conv, gate_conv=gate_conv,
               lam_gridB=[float(x) for x in lam_grid],
               rowsB_eps005=res["rowsB"],
               runtime_s=time.time() - t_start)
    with open(OUTDIR + "/axi_results.json", "w") as fh:
        json.dump(_jsonable(out), fh, indent=2)

    # ---- verdict table ----------------------------------------------------
    print("\n== VERDICTS vs benchmark targets (family B = paper's ansatz; "
          "family A = task's literal SDiff family) ==")
    sB, sA = res["solB"], res["solA"]
    rows = [
        ("lambda*  (0.42 +/- 0.05)", sA["lam"], sB["lam"]),
        ("g*       (1.31 +/- 0.04)", sA["g"], sB["g"]),
        ("V*       (1.409 +/- 0.010)", sA["V"], sB["V"]),
        ("Erot/E   (0.2500)", sA["erot_frac"], sB["erot_frac"]),
        ("c_paper  (2.37 +/- 0.09)", sA["c_paper"], sB["c_paper"]),
        ("kappa_pap(0.802 +/- 0.018)", sA["kappa_paper"], sB["kappa_paper"]),
        ("exponent (2/3-class; 2a: 0.75)", fits["A"][0], fits["B"][0]),
    ]
    print("   %-32s %12s %12s" % ("target", "family A", "family B"))
    for name, a, b in rows:
        print("   %-32s %12.4f %12.4f" % (name, a, b))
    print("\nruntime: %.1f s" % (time.time() - t_start))
    print("wrote axi_results.json, axi_landscape.png, axi_epsscan.png")


if __name__ == "__main__":
    main()
