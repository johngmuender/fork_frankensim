#!/usr/bin/env python3
"""
axi3_solve.py — Tier-2b STEP 3: enlarged-basis isorotating closure that
adjudicates FINDING F-R4.

Question: over a configuration family that can TILT the field's internal
direction relative to the spatial angle (the freedom neither Step-2 family
had), does the isorotating closure at eps -> 0 approach the paper's
shape-exact endpoint (g -> 3/2, c_paper -> 2.515), or stay near the scaling
rung (g ~ 1, c_paper ~ 2.05)?  And at eps = 0.05, does it approach the <r1>
benchmark (c = 2.37 +/- 0.09, kappa = 0.802 +/- 0.018, g* = 1.31 +/- 0.04,
V = 1.409 +/- 0.010)?

----------------------------------------------------------------------------
ENLARGED ANSATZ (general axisymmetric unit-winding field)
----------------------------------------------------------------------------
On the phi = 0 half-plane, in spherical coordinates (r, theta), mu = cos(theta):

  q = ( cos F,  sin F sin Theta,  0,  sin F cos Theta ),        with
  F(r,theta)     = f_base( r_eff ),   r_eff = r * exp( sum_kl a_kl B_k(r) C_l(mu) )
  Theta(r,theta) = theta + sum_kl b_kl B_k(r) S_l(theta)

  B_k(r) : K normalised cubic B-spline bumps spanning the soliton support
           (centres r_sup (k+1/2)/K, width r_sup/K; peak value 1)
  C_l(mu): even Legendre polynomials  {1, P2, [P4]}   (l = 0, 2, [4])
  S_l    : odd-in-mu tilt functions vanishing at theta = 0, pi/2, pi:
           {sin t cos t, sin t cos^3 t, [sin t cos^5 t]}   (l = 1, 2, [3])

Boundary conditions EXACT by construction: Theta(r,0) = 0, Theta(r,pi) = pi,
Theta(r,pi/2) = pi/2, so the field points along +/-z on the axis and keeps
the hedgehog-class z -> -z symmetry (C_l even, S_l odd in mu).  The
phi-dependence is the standard unit winding (analytic L_phi), so the whole
Step-2 quadrature engine (axi_solve.sector_pieces) is REUSED unchanged; the
builder uses signed sin(theta) = rho/r and theta = arctan2(rho, z), which
extends it to the ghost rings with the exact reflection parities.  These are
continuous deformations of the hedgehog => degree stays 1; VERIFIED at every
accepted optimum (gate: |B - 1| < 0.5%).  Overall dilation d is analytic
(E2 ~ d, E4 ~ 1/d, E6 ~ 1/d^3, E0 ~ d^3, I ~ d^3): quadratures at d = 1,
d minimised per evaluation (golden in ln d).

----------------------------------------------------------------------------
CLOSURE (derived; verified against the Step-2 solutions in gate G0)
----------------------------------------------------------------------------
From E_tot = E_static + L^2/(2I) and the clock c kappa = E_tot (c = 2L,
kappa = L/I):  2L^2/I = E_static + L^2/(2I)  =>  L^2 = (2/3) I E_static at
the solution.  NOT substituted into R (that would ignore inertia): the
physical solution minimises R(p; L) = E_static(p) + L^2/(2 I(p)) over p at
FIXED L; the clock then selects L.  Iteration: (i) L from the hedgehog
closure; (ii) Nelder-Mead minimisation of R(.; L) (deterministic seeds:
hedgehog, spheroidal-shell projections, tilt seeds); (iii) L updated by the
closed form L^2 = (2/3) I(p*) E_static(p*); repeat to |dL/L| < 1e-7.  After
convergence the (d, L) fixed point is re-solved exactly at fixed p*, so the
reported clock residual is machine-level.

----------------------------------------------------------------------------
THE STEP-3 FINDING THIS SCRIPT ESTABLISHES AND VERIFIES (halo saturation)
----------------------------------------------------------------------------
Once tilt modes exist, the closure has a runaway channel none of the Step-2
families contained.  Consider adding to any configuration a low-amplitude
"halo": F = eta (small) over a large region far from the core, with the
internal direction tilted to the equator (Theta ~ pi/2, allowed off-axis).
Exact small-eta energetics (m = a0 = 1, our unit conventions):
    dE_static = dE0 = (1/4pi) INT (1 - cos F) dV ~ (1/8pi) INT eta^2 dV
    dI        = INT 2 sin^2 F sin^2 Theta dV     ~  2      INT eta^2 dV
    (dE2, dE4 ~ t k^2 eta^2 -> 0 for soft halos; dE6 ~ eta^6 negligible)
so  dR = dE_static - (L^2/2I^2) dI = [INT eta^2 dV] (1/(8pi) - kappa_ours^2),
kappa_ours = L/I.  A spinning configuration is stable against halo growth
ONLY if kappa_ours <= 1/sqrt(8pi) = 0.19947, i.e. in paper units (KMAP =
2 sqrt(pi)):  kappa_paper <= 1/sqrt(2) = 0.70711 — the isorotation channel's
mass threshold (the classic omega <= m_meson criterion; here m_iso = 1/sqrt2
in paper units from (1-cos f) ~ f^2/2).  EVERY Step-2 closure solution and
EVERY corpus tuple sits ABOVE it (solA 0.9414, solB 0.9069, endpoints
0.9354/0.9009, paper deep-BPS 0.7638, <r1> 0.802) => all are saddle points
of the fixed-L Routhian, unstable to tilted-halo growth.  The halo grows
until L^2/I^2 = 1/(8pi):

  SATURATED CLOSURE (exact fixed point, halo is G-neutral):
    kappa_ours* = 1/sqrt(8pi)   (kappa_paper* = 1/sqrt2, EXACT)
    core minimises  G := E_static - I/(16 pi)   (Legendre-reduced problem;
      dG = 0 for the ideal halo, dG > 0 for everything else, so G is
      bounded below and the reduced problem is WELL-POSED)
    L = sqrt(8pi) G*,  E_static_tot = (3/2) G*,  I_tot = 8pi G*,
    R = 2 G*,  E_rot/E_tot = 1/4 (automatic),
    c_ours = 2 sqrt(8pi) G*,   c_paper = c_ours/sqrt(2 pi^3) = (4/pi) G*.
  RIGOROUS LOWER BOUND (any t >= 0, any degree-1 field): pointwise
  sin^2 Theta <= 1 gives G >= E6 + (1/8pi) INT (1-cos F)^2 dV, and the
  sextic+potential Bogomolny bound (normalisation verified against the
  step-1 compacton bound 32 sqrt2/15) gives
    G >= 2 INT_0^pi sin^2 f sqrt(U_eff) df,  U_eff = (1-cos f)^2 / 2
      = pi/sqrt2   (closed form; = 16 sqrt2 INT_0^{pi/2} sin^4 u cos^2 u du)
    =>  c_paper >= (4/pi)(pi/sqrt2) = 2 sqrt2 = 2.82843.
  The paper's shape-exact endpoint (c0 = 2.5147, kappa = sqrt(7/12) =
  0.7638) and the <r1> benchmark (2.37, 0.802) are therefore EXCLUDED as
  unrestricted-closure limits: c0 < 2 sqrt2 and both kappas != 1/sqrt2.

What the script does about it:
  * runs the task-literal DIRECT closure (R1 2x2 ablation, extended basis,
    R2) and reports the boundary-hitting (the family-too-small signal);
  * scans the parametrisation box size (boundary scan) showing c rising and
    kappa falling monotonically toward the saturated values;
  * solves the WELL-POSED saturated closure (G-minimisation over the same
    enlarged basis, same 2x2 ablation) whose optimum IS interior;
  * verifies the halo energetics with an explicit engine-level probe (core
    + tilted shell, measured kappa_crit vs the derived 1/sqrt(8pi));
  * computes the marginal-mode spectrum (d2E/dp^2, dI/dp) at the compacton
    and at the R1 optimum (the F-R4-relevant spectrum).

Deterministic, no RNG.  Outputs: axi3_results.json, axi3_modes.png
(axi3_run.log via tee; axi3_RESULTS.md written from the printed tables).
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
import axi_solve as ax                    # Step-2 engine (REUSED, not modified)

# ----------------------------------------------------------------------
# frozen conventions (identical to Step 2; NOTHING fitted)
# ----------------------------------------------------------------------
TFROZEN = ax.TFROZEN
RSTAR = ax.RSTAR
CMAP_L = ax.CMAP_L                        # c_paper   = 2 L * CMAP_L
KMAP = ax.KMAP                            # kappa_pap = (L/I) * KMAP
E6_COMPACTON = ax.E6_COMPACTON
I_COMPACTON = ax.I_COMPACTON
SIXTEEN_PI = 16.0 * math.pi
KAPPA_CRIT_OURS = 1.0 / math.sqrt(8.0 * math.pi)
KAPPA_CRIT_PAPER = KAPPA_CRIT_OURS * KMAP            # = 1/sqrt(2)
G_BOUND = math.pi / math.sqrt(2.0)                   # closed form (verified)
C_PAPER_LOWER = 4.0 / math.pi * G_BOUND              # = 2 sqrt(2)

QUICK = "--quick" in sys.argv
ABOX, BBOX = 0.60, 1.50                   # soft parameter boxes (report if hit)
AMAXCAP = 0.45                            # soft cap on the shrink side of A

if QUICK:
    N_OPT, N_POL, N_FIN, N_COARSE = 112, 176, 256, 176
    B_TRI, B_COLD, B_WARM, B_POL = 120, 500, 180, 100
    B_SCAN = 350
else:
    N_OPT, N_POL, N_FIN, N_COARSE = 240, 480, 720, 360
    B_TRI, B_COLD, B_WARM, B_POL = 500, 3200, 800, 300
    B_SCAN = 1500

TARGETS_R2 = dict(c=(2.37, 0.09), kappa=(0.802, 0.018), g=(1.31, 0.04),
                  V=(1.409, 0.010))


# ----------------------------------------------------------------------
# basis
# ----------------------------------------------------------------------
def legC(l, mu):
    if l == 0:
        return np.ones_like(mu)
    if l == 2:
        return 0.5 * (3.0 * mu * mu - 1.0)
    if l == 4:
        m2 = mu * mu
        return (35.0 * m2 * m2 - 30.0 * m2 + 3.0) / 8.0
    raise ValueError(l)


def tiltS(l, sth, cth):
    if l == 1:
        return sth * cth
    if l == 2:
        return sth * cth ** 3
    if l == 3:
        return sth * cth ** 5
    raise ValueError(l)


class GenBasis:
    def __init__(self, r_sup, K=3, a_l=(0, 2), b_l=(1, 2)):
        self.r_sup = float(r_sup)
        self.K = K
        s = self.r_sup / K
        self.cent = s * (np.arange(K) + 0.5)
        self.w = s
        self.a_modes = [(k, l) for k in range(K) for l in a_l]
        self.b_modes = [(k, l) for k in range(K) for l in b_l]
        self.na, self.nb = len(self.a_modes), len(self.b_modes)
        self.n = self.na + self.nb
        self.labels = (["a[%d,%d]" % m for m in self.a_modes]
                       + ["b[%d,%d]" % m for m in self.b_modes])
        self.rd = np.linspace(0.0, self.r_sup * 1.05, 601)
        self.Bd = np.array([self.bump(k, self.rd) for k in range(K)])
        mud = np.linspace(-1.0, 1.0, 81)
        self.Cd = {l: legC(l, mud)
                   for l in sorted(set(l for _, l in self.a_modes))}

    def bump(self, k, r):
        x = np.abs((np.asarray(r, float) - self.cent[k]) / self.w)
        return np.where(x < 1.0, 1.0 - 1.5 * x * x * (1.0 - 0.5 * x),
                        0.25 * np.clip(2.0 - x, 0.0, None) ** 3)

    def afield(self, r, cth, a):
        A, Bk, Cl = None, {}, {}
        for (k, l), coef in zip(self.a_modes, a):
            if coef == 0.0:
                continue
            if k not in Bk:
                Bk[k] = self.bump(k, r)
            if l not in Cl:
                Cl[l] = legC(l, cth)
            term = coef * Bk[k] * Cl[l]
            A = term if A is None else A + term
        return A

    def bfield(self, r, sth, cth, b):
        T, Bk, Sl = None, {}, {}
        for (k, l), coef in zip(self.b_modes, b):
            if coef == 0.0:
                continue
            if k not in Bk:
                Bk[k] = self.bump(k, r)
            if l not in Sl:
                Sl[l] = tiltS(l, sth, cth)
            term = coef * Bk[k] * Sl[l]
            T = term if T is None else T + term
        return T

    def a_extremes(self, a):
        A = np.zeros((self.rd.size, 81))
        for (k, l), coef in zip(self.a_modes, a):
            if coef != 0.0:
                A += coef * self.Bd[k][:, None] * self.Cd[l][None, :]
        return float(A.min()), float(A.max())


def make_qfun(fint, basis, a, b):
    use_a = any(v != 0.0 for v in a)
    use_b = any(v != 0.0 for v in b)

    def qfun(R, Z):
        r = np.hypot(R, Z)
        sth = R / r
        cth = Z / r
        if use_a:
            F = fint(r * np.exp(basis.afield(r, cth, a)))
        else:
            F = fint(r)
        th = np.arctan2(R, Z)
        Th = th + basis.bfield(r, sth, cth, b) if use_b else th
        sF = np.sin(F)
        q1 = sF * np.sin(Th)
        return np.cos(F), q1, sF * np.cos(Th), q1 / R

    return qfun


# ----------------------------------------------------------------------
# closure pieces (all d-dependence exact/closed form)
# ----------------------------------------------------------------------
def R_from(P, t, L, d):
    e2, e4 = P["E2"] * d, P["E4"] / d
    e6, e0 = P["E6"] / d ** 3, P["E0"] * d ** 3
    I3 = P["I"] * d ** 3
    Estat = t * (e2 + e4) + e6 + e0
    return Estat + (L * L / (2.0 * I3) if L else 0.0), Estat, I3


def min_d(P, t, L, dlo=0.35, dhi=2.8):
    w, R = ax.golden(lambda u: R_from(P, t, L, math.exp(u))[0],
                     math.log(dlo), math.log(dhi), iters=56)
    return math.exp(w), R


def G_from(P, t, d):
    return (t * (P["E2"] * d + P["E4"] / d) + P["E6"] / d ** 3
            + (P["E0"] - P["I"] / SIXTEEN_PI) * d ** 3)


def min_d_G(P, t, dlo=0.35, dhi=2.8):
    w, G = ax.golden(lambda u: G_from(P, t, math.exp(u)),
                     math.log(dlo), math.log(dhi), iters=56)
    return math.exp(w), G


def dl_fixed_point(P, t, L0):
    """Exact (d, L) clock fixed point at frozen sectors P (closed form)."""
    L = max(L0, 1e-6)
    for _ in range(200):
        d, R = min_d(P, t, L)
        _, Estat, I3 = R_from(P, t, L, d)
        Ln = math.sqrt((2.0 / 3.0) * I3 * Estat)
        if abs(Ln - L) < 1e-14 * Ln:
            L = Ln
            break
        L = Ln
    d, R = min_d(P, t, L)
    _, Estat, I3 = R_from(P, t, L, d)
    return L, d, R, Estat, I3


# ----------------------------------------------------------------------
# objective (fixed domain => smooth deterministic function of p)
# ----------------------------------------------------------------------
class Objective:
    def __init__(self, fint, basis, t, N, dom, active, stretch_lim,
                 mode="direct", amaxcap=AMAXCAP):
        self.fint, self.basis, self.t = fint, basis, t
        self.N, self.dom = N, dom
        self.active = list(active)
        self.stretch_lim = stretch_lim
        self.amaxcap = amaxcap
        self.mode = mode
        self.L = None
        self.nfev = 0

    def split(self, x):
        p = np.zeros(self.basis.n)
        p[self.active] = x
        return p[:self.basis.na], p[self.basis.na:]

    def sectors(self, x, N=None):
        a, b = self.split(x)
        NN = self.N if N is None else N
        return ax.sector_pieces(make_qfun(self.fint, self.basis, a, b),
                                NN, self.dom, self.dom)

    def penalty(self, a, b):
        amin, amax = self.basis.a_extremes(a)
        if -amin > self.stretch_lim:                 # hard: support > domain
            return None, -amin, amax
        pen = 400.0 * max(0.0, -amin - (self.stretch_lim - 0.02)) ** 2
        pen += 400.0 * max(0.0, amax - self.amaxcap) ** 2
        for v in a:
            pen += 200.0 * max(0.0, abs(v) - ABOX) ** 2
        for v in b:
            pen += 200.0 * max(0.0, abs(v) - BBOX) ** 2
        return pen, -amin, amax

    def __call__(self, x):
        self.nfev += 1
        a, b = self.split(x)
        pen, stretch, amax = self.penalty(a, b)
        if pen is None:
            return 1.0e6 * (1.0 + stretch)
        P = self.sectors(x)
        if self.mode == "direct":
            _, val = min_d(P, self.t, self.L)
        else:
            _, val = min_d_G(P, self.t)
        return val + pen


# ----------------------------------------------------------------------
# deterministic Nelder-Mead (scipy absent)
# ----------------------------------------------------------------------
def nelder_mead(f, x0, steps, maxfev=1000, ftol=1e-10, xtol=1e-9):
    n = len(x0)
    if n == 0:
        return np.asarray(x0, float), f(np.asarray(x0, float)), 1
    sim = [np.asarray(x0, float)]
    fx = [f(sim[0])]
    for i in range(n):
        xi = sim[0].copy()
        xi[i] += steps[i]
        sim.append(xi)
        fx.append(f(xi))
    nfev = n + 1
    while nfev < maxfev:
        o = np.argsort(fx)
        sim = [sim[j] for j in o]
        fx = [fx[j] for j in o]
        if (fx[-1] - fx[0] < ftol * (abs(fx[0]) + 1e-12)
                and max(np.max(np.abs(s - sim[0])) for s in sim[1:]) < xtol):
            break
        cen = np.mean(sim[:-1], axis=0)
        xr = cen + (cen - sim[-1])
        fr = f(xr)
        nfev += 1
        if fr < fx[0]:
            xe = cen + 2.0 * (cen - sim[-1])
            fe = f(xe)
            nfev += 1
            if fe < fr:
                sim[-1], fx[-1] = xe, fe
            else:
                sim[-1], fx[-1] = xr, fr
        elif fr < fx[-2]:
            sim[-1], fx[-1] = xr, fr
        else:
            xc = cen + 0.5 * (sim[-1] - cen)
            fc = f(xc)
            nfev += 1
            if fc < fx[-1]:
                sim[-1], fx[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    sim[i] = sim[0] + 0.5 * (sim[i] - sim[0])
                    fx[i] = f(sim[i])
                nfev += n
    i = int(np.argmin(fx))
    return sim[i], fx[i], nfev


def optimize(obj, seeds, steps, b_tri, b_main):
    best = None
    if len(seeds) > 1:
        for s in seeds:
            x, fv, _ = nelder_mead(obj, s, steps, maxfev=b_tri)
            if best is None or fv < best[1]:
                best = (x, fv)
    else:
        x0 = np.asarray(seeds[0], float)
        best = (x0, obj(x0))
    x, fv, _ = nelder_mead(obj, best[0], steps, maxfev=b_main)
    x, fv, _ = nelder_mead(obj, x, steps * 0.25, maxfev=max(b_main // 3, 60))
    return x, fv


# ----------------------------------------------------------------------
# diagnostics
# ----------------------------------------------------------------------
def moments(qfun, N, dom):
    h = dom / N
    x = (np.arange(N) + 0.5) * h
    R, Z = np.meshgrid(x, x, indexing="ij")
    q0 = qfun(R, Z)[0]
    W = 4.0 * np.pi * R * h * h
    dens = (1.0 - q0) * W
    rr = float(np.sum(dens * R * R))
    zz = float(np.sum(dens * Z * Z))
    return float(np.sum(dens)), rr, zz, math.sqrt(2.0 * zz / rr)


def sph_check(fint, basis, a, b, dom, n=1400):
    """Independent SPHERICAL-coordinate quadrature for E0 and I."""
    hr, hth = dom / n, math.pi / n
    r = (np.arange(n) + 0.5) * hr
    th = (np.arange(n) + 0.5) * hth
    Rg, Tg = np.meshgrid(r, th, indexing="ij")
    sth, cth = np.sin(Tg), np.cos(Tg)
    A = basis.afield(Rg, cth, a)
    F = fint(Rg * np.exp(A)) if A is not None else fint(Rg)
    Tt = basis.bfield(Rg, sth, cth, b)
    Th = Tg + (Tt if Tt is not None else 0.0)
    w = 2.0 * np.pi * Rg * Rg * sth * hr * hth
    E0 = float(np.sum((1.0 - np.cos(F)) * w)) / (4.0 * np.pi)
    I = float(np.sum(2.0 * np.sin(F) ** 2 * np.sin(Th) ** 2 * w))
    return E0, I


def paperize(L, I3):
    return dict(c_ours=2.0 * L, c_paper=2.0 * L * CMAP_L,
                kappa_ours=L / I3, kappa_paper=(L / I3) * KMAP)


def boundary_report(obj, a, b):
    pen, stretch, amax = obj.penalty(a, b)
    return dict(penalty=pen, stretch_used=stretch,
                stretch_frac=stretch / obj.stretch_lim, amax=amax,
                amax_frac=amax / obj.amaxcap,
                flags=dict(
                    stretch=bool(stretch > 0.90 * obj.stretch_lim),
                    amax=bool(amax > 0.90 * obj.amaxcap),
                    abox=bool(np.max(np.abs(a), initial=0.0) > 0.90 * ABOX),
                    bbox=bool(np.max(np.abs(b), initial=0.0) > 0.90 * BBOX)))


# ----------------------------------------------------------------------
# finalisation (direct closure and saturated closure)
# ----------------------------------------------------------------------
def finalize_direct(obj, x, L_guess, N, P0N):
    P = obj.sectors(x, N=N)
    L, d, R, Estat, I3 = dl_fixed_point(P, obj.t, L_guess)
    d0, _ = min_d(P0N, obj.t, 0.0)
    a, b = obj.split(x)
    _, _, _, lam_eff = moments(make_qfun(obj.fint, obj.basis, a, b), N,
                               obj.dom)
    rep = dict(N=N, L=L, d=d, R=R, Estat=Estat, I=I3, I1=P["I"],
               Erot=L * L / (2.0 * I3),
               erot_frac=(L * L / (2.0 * I3)) / R,
               clock_resid_rel=(2.0 * L * L / I3 - R) / R,
               g=P["I"] / P0N["I"], d0=d0, V=(d / d0) ** 3,
               V_E0=(P["E0"] * d ** 3) / (P0N["E0"] * d0 ** 3),
               lam_eff=lam_eff, degree=P["BDEG"],
               degree_ok=bool(abs(P["BDEG"] - 1.0) < 5e-3),
               eps_at_solution=obj.t * (P["E2"] * d + P["E4"] / d)
               / (P["E6"] / d ** 3 + P["E0"] * d ** 3),
               a=list(map(float, a)), b=list(map(float, b)))
    rep.update(paperize(L, I3))
    rep.update(boundary_report(obj, a, b))
    return rep


def finalize_G(obj, x, N, P0N):
    """Saturated closure: core minimises G; halo closes the fixed point."""
    P = obj.sectors(x, N=N)
    d, Gv = min_d_G(P, obj.t)
    L = math.sqrt(8.0 * math.pi) * Gv
    I_tot = 8.0 * math.pi * Gv
    Estat_tot = 1.5 * Gv
    I_core = P["I"] * d ** 3
    halo_D = (I_tot - I_core) / 2.0
    E0_halo = halo_D / (8.0 * math.pi)
    _, Estat_core, _ = R_from(P, obj.t, 0.0, d)
    d0, _ = min_d(P0N, obj.t, 0.0)
    a, b = obj.split(x)
    _, _, _, lam_eff = moments(make_qfun(obj.fint, obj.basis, a, b), N,
                               obj.dom)
    rep = dict(N=N, G=Gv, d=d, d0=d0, L=L,
               c_ours=2.0 * L, c_paper=4.0 / math.pi * Gv,
               kappa_ours=KAPPA_CRIT_OURS, kappa_paper=KAPPA_CRIT_PAPER,
               Estat_tot=Estat_tot, I_tot=I_tot, R=2.0 * Gv,
               erot_frac=0.25,
               Estat_core=Estat_core, I_core=I_core, halo_D=halo_D,
               E0_halo=E0_halo, halo_ok=bool(halo_D >= 0.0),
               consistency=(Estat_core + E0_halo) / Estat_tot - 1.0,
               g_total=I_tot / (P0N["I"] * d0 ** 3),
               g_core=P["I"] / P0N["I"],
               g_core_scaled=I_core / (P0N["I"] * d0 ** 3),
               V_core=(d / d0) ** 3,
               V_E0_tot=(P["E0"] * d ** 3 + E0_halo)
               / (P0N["E0"] * d0 ** 3),
               lam_eff=lam_eff, degree=P["BDEG"],
               degree_ok=bool(abs(P["BDEG"] - 1.0) < 5e-3),
               a=list(map(float, a)), b=list(map(float, b)))
    rep.update(boundary_report(obj, a, b))
    return rep


# ----------------------------------------------------------------------
# drivers
# ----------------------------------------------------------------------
def run_closure(tag, fint, basis, t, active, seeds, dom, stretch_lim,
                P0_cache, b_cold=None, max_outer=12, log=True,
                amaxcap=AMAXCAP):
    """Task-literal direct closure: min_p R(p; L), clock update on L."""
    t0 = time.time()
    obj = Objective(fint, basis, t, N_OPT, dom, active, stretch_lim,
                    amaxcap=amaxcap)
    na = basis.na
    steps = np.array([0.06 if i < na else 0.15 for i in obj.active])
    b_cold = B_COLD if b_cold is None else b_cold

    def P0(N):
        if N not in P0_cache:
            P0_cache[N] = ax.sector_pieces(
                make_qfun(fint, basis, np.zeros(na), np.zeros(basis.nb)),
                N, dom, dom)
        return P0_cache[N]

    L, _, _, _, _ = dl_fixed_point(P0(N_OPT), t, 5.0)
    x = np.zeros(len(obj.active))
    hist = []
    outer_conv = True
    if len(obj.active) == 0:
        rep = finalize_direct(obj, x, L, N_FIN, P0(N_FIN))
        rep["coarse"] = finalize_direct(obj, x, L, N_COARSE, P0(N_COARSE))
        rep["hist"], rep["nfev"], rep["outer_converged"] = [], 0, True
        rep["runtime_s"] = time.time() - t0
        if log:
            print("   %-10s (0 params)  c_pap=%.4f kap_pap=%.4f g=%.4f "
                  "V=%.4f" % (tag, rep["c_paper"], rep["kappa_paper"],
                              rep["g"], rep["V"]))
        return rep, x, obj

    for it in range(max_outer):
        obj.L = L
        if it == 0:
            x, fv = optimize(obj, seeds, steps, B_TRI, b_cold)
        else:
            x, fv, _ = nelder_mead(obj, x, steps * 0.2, maxfev=B_WARM)
        P = obj.sectors(x)
        Lx, d, R, Estat, I3 = dl_fixed_point(P, t, L)
        rel = abs(Lx - L) / Lx
        hist.append(dict(it=it, L=L, L_new=Lx, R=fv,
                         c_paper=2.0 * Lx * CMAP_L,
                         kappa_paper=Lx / I3 * KMAP, nfev=obj.nfev))
        if log:
            print("   %-10s outer %2d: L %.6f -> %.6f (dL/L %.2e)  R=%.6f  "
                  "nfev=%d" % (tag, it, L, Lx, rel, fv, obj.nfev))
        L = Lx
        if rel < 1.0e-7 and it >= 2:
            break
    else:
        outer_conv = False

    objP = Objective(fint, basis, t, N_POL, dom, active, stretch_lim,
                     amaxcap=amaxcap)
    for _ in range(2):
        objP.L = L
        x, fv, _ = nelder_mead(objP, x, steps * 0.08, maxfev=B_POL)
        P = objP.sectors(x)
        L, d, R, Estat, I3 = dl_fixed_point(P, t, L)
    if log:
        print("   %-10s polish (N=%d): L=%.6f  R=%.6f  nfev=%d"
              % (tag, N_POL, L, R, objP.nfev))

    rep = finalize_direct(obj, x, L, N_FIN, P0(N_FIN))
    rep["coarse"] = finalize_direct(obj, x, L, N_COARSE, P0(N_COARSE))
    rep["conv_c_rel"] = abs(rep["coarse"]["c_paper"] / rep["c_paper"] - 1.0)
    rep["hist"], rep["nfev"] = hist, obj.nfev + objP.nfev
    rep["outer_converged"] = outer_conv
    rep["runtime_s"] = time.time() - t0
    if log:
        print("   %-10s FINAL (N=%d): c_pap=%.4f kap_pap=%.4f g*=%.4f "
              "lam_eff=%.4f V=%.4f V_E0=%.4f Erot/E=%.6f resid=%.1e "
              "B=%.5f [coarse: c=%.4f]  (%.0f s)"
              % (tag, N_FIN, rep["c_paper"], rep["kappa_paper"], rep["g"],
                 rep["lam_eff"], rep["V"], rep["V_E0"], rep["erot_frac"],
                 rep["clock_resid_rel"], rep["degree"],
                 rep["coarse"]["c_paper"], rep["runtime_s"]))
        fl = rep["flags"]
        if any(fl.values()):
            print("   %-10s !! BOUNDARY: %s (stretch %.0f%%, amax %.0f%% of "
                  "caps) — family-too-small signal"
                  % (tag, [k for k, v in fl.items() if v],
                     100 * rep["stretch_frac"], 100 * rep["amax_frac"]))
    return rep, x, obj


def run_gmin(tag, fint, basis, t, active, seeds, dom, stretch_lim, P0_cache,
             log=True):
    """Well-posed saturated closure: minimise G = E_static - I/(16 pi)."""
    t0 = time.time()
    obj = Objective(fint, basis, t, N_OPT, dom, active, stretch_lim,
                    mode="G")
    na = basis.na
    steps = np.array([0.06 if i < na else 0.15 for i in obj.active])

    def P0(N):
        if N not in P0_cache:
            P0_cache[N] = ax.sector_pieces(
                make_qfun(fint, basis, np.zeros(na), np.zeros(basis.nb)),
                N, dom, dom)
        return P0_cache[N]

    if len(obj.active) == 0:
        x = np.zeros(0)
    else:
        x, fv = optimize(obj, seeds, steps, B_TRI, B_COLD)
        objP = Objective(fint, basis, t, N_POL, dom, active, stretch_lim,
                         mode="G")
        x, fv, _ = nelder_mead(objP, x, steps * 0.08, maxfev=B_POL)
        obj.nfev += objP.nfev
    rep = finalize_G(obj, x, N_FIN, P0(N_FIN))
    rep["coarse"] = finalize_G(obj, x, N_COARSE, P0(N_COARSE))
    rep["conv_c_rel"] = abs(rep["coarse"]["c_paper"] / rep["c_paper"] - 1.0)
    rep["nfev"] = obj.nfev
    rep["runtime_s"] = time.time() - t0
    if log:
        print("   %-10s G*=%.6f  ->  c_pap=%.4f (kappa_pap=1/sqrt2)  "
              "g_tot=%.3f g_core=%.4f V_E0=%.4f halo_D=%.3f B=%.5f "
              "[coarse: c=%.4f]  (%.0f s, nfev=%d)"
              % (tag, rep["G"], rep["c_paper"], rep["g_total"],
                 rep["g_core_scaled"], rep["V_E0_tot"], rep["halo_D"],
                 rep["degree"], rep["coarse"]["c_paper"], rep["runtime_s"],
                 rep["nfev"]))
        fl = rep["flags"]
        if any(fl.values()):
            print("   %-10s !! BOUNDARY: %s" % (tag,
                                                [k for k, v in fl.items() if v]))
    return rep, x, obj


# ----------------------------------------------------------------------
# seeds
# ----------------------------------------------------------------------
def bseed_alpha2(lam):
    mu = np.linspace(-1.0, 1.0, 201)
    h = 0.5 * np.log(lam ** (2.0 / 3.0) * (1.0 - mu * mu)
                     + lam ** (-4.0 / 3.0) * mu * mu)
    X = np.stack([np.ones_like(mu), legC(2, mu)], axis=1)
    coef, *_ = np.linalg.lstsq(X, h, rcond=None)
    return float(coef[1])


def seed_vec(basis, active, a_map=None, b_map=None):
    p = np.zeros(basis.n)
    if a_map:
        for i, (k, l) in enumerate(basis.a_modes):
            p[i] = a_map.get(l, 0.0)
    if b_map:
        for i, (k, l) in enumerate(basis.b_modes):
            p[basis.na + i] = b_map.get(l, 0.0)
    return p[list(active)]


# ----------------------------------------------------------------------
# marginal-mode spectrum (R3)
# ----------------------------------------------------------------------
def spectrum(fint, basis, p_full, t, dstar, N, dom, stretch_lim):
    obj = Objective(fint, basis, t, N, dom, list(range(basis.n)),
                    stretch_lim)
    x0 = np.asarray(p_full, float)
    P = obj.sectors(x0, N=N)
    _, E0stat, I0 = R_from(P, t, 0.0, dstar)
    rows = []
    for i in range(basis.n):
        h = 0.02 if i < basis.na else 0.04
        vals = []
        for sgn in (1.0, -1.0):
            xp = x0.copy()
            xp[i] += sgn * h
            Pp = obj.sectors(xp, N=N)
            _, Es, I3 = R_from(Pp, t, 0.0, dstar)
            vals.append((Es, I3))
        (Ep, Ip), (Em, Im) = vals
        d2E = (Ep - 2.0 * E0stat + Em) / (h * h)
        dI = (Ip - Im) / (2.0 * h)
        rows.append(dict(mode=basis.labels[i], h=h,
                         dE=(Ep - Em) / (2.0 * h), d2E=d2E, dI=dI,
                         gain=dI / d2E if d2E != 0 else float("inf")))
    return dict(Estat=E0stat, I=I0, dstar=dstar, rows=rows)


def print_spectrum(sp, title):
    print("   %s  (E_static=%.6f, I=%.4f, d*=%.4f)"
          % (title, sp["Estat"], sp["I"], sp["dstar"]))
    print("   %-9s %12s %12s %12s %12s" % ("mode", "dE/dp", "d2E/dp2",
                                           "dI/dp", "dI/d2E"))
    for r in sp["rows"]:
        print("   %-9s %+12.5f %12.5f %+12.4f %+12.4f"
              % (r["mode"], r["dE"], r["d2E"], r["dI"], r["gain"]))


# ----------------------------------------------------------------------
# explicit halo probe (engine-level verification of the threshold)
# ----------------------------------------------------------------------
def halo_probe(fint, rsup, t, kappa_sol, label, N=None):
    """Core hedgehog + explicit low-amplitude tilted shell far outside the
    core.  Measures kappa_crit = sqrt(2 dE_static / dI) with the FULL engine
    (all sectors) and compares with the derived ideal 1/sqrt(8 pi).  The
    shell tilt uses the parity-exact S1 = sin(t)cos(t) (partial tilt), so
    the measured kappa_crit sits slightly above the ideal."""
    N = (N_POL if not QUICK else 160) if N is None else N
    rc, w = 2.1 * rsup, 0.45 * rsup
    dom = 3.2 * rsup

    def bumpf(r):
        x = np.abs((r - rc) / w)
        return np.where(x < 1.0, 1.0 - 1.5 * x * x * (1.0 - 0.5 * x),
                        0.25 * np.clip(2.0 - x, 0.0, None) ** 3)

    def mk(eta, beta):
        def qf(R, Z):
            r = np.hypot(R, Z)
            chi = bumpf(r)
            F = fint(r) + eta * chi
            th = np.arctan2(R, Z)
            sth, cth = R / r, Z / r
            Th = th + beta * sth * cth * chi
            sF = np.sin(F)
            q1 = sF * np.sin(Th)
            return np.cos(F), q1, sF * np.cos(Th), q1 / R
        return qf

    P00 = ax.sector_pieces(mk(0.0, 0.0), N, dom, dom)
    _, E00, _ = R_from(P00, t, 0.0, 1.0)
    res = {}
    for eta in (0.02, 0.04):
        P = ax.sector_pieces(mk(eta, 1.0), N, dom, dom)
        _, Es, _ = R_from(P, t, 0.0, 1.0)
        dE, dI = Es - E00, P["I"] - P00["I"]
        res["eta_%.2f" % eta] = dict(dE=dE, dI=dI,
                                     kcrit=math.sqrt(2.0 * dE / dI))
    k1 = res["eta_0.02"]["kcrit"]
    k2 = res["eta_0.04"]["kcrit"]
    quad = res["eta_0.04"]["dI"] / res["eta_0.02"]["dI"]
    print("   %s: kappa_crit(shell) = %.5f / %.5f (eta 0.02/0.04; ideal "
          "1/sqrt(8pi) = %.5f)" % (label, k1, k2, KAPPA_CRIT_OURS))
    print("       dI ratio eta->2eta: %.3f (expect ~4 for quadratic)  |  "
          "solution kappa_ours = %.5f -> UNSTABLE: %s"
          % (quad, kappa_sol, kappa_sol > k1))
    res.update(kappa_solution=kappa_sol, unstable=bool(kappa_sol > k1),
               kappa_crit_ideal=KAPPA_CRIT_OURS, N=N, rc=rc, w=w, beta=1.0)
    return res


# ----------------------------------------------------------------------
# figure (house palette, validated)
# ----------------------------------------------------------------------
INK, GRID_C = "#333333", "#dddddd"
BLUE, AQUA, RED, VIOLET = "#2a78d6", "#1baf7a", "#e34948", "#4a3aa7"


def style_ax(a):
    a.set_facecolor("white")
    a.grid(color=GRID_C, lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.tick_params(colors=INK, labelsize=9)


def fig_modes(fint1, basis1, repG1, repG2, sp_comp, scan, gfam, path):
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.2), dpi=160)
    fig.patch.set_facecolor("white")

    def fields(fint, basis, rep, n=340, fac=1.25):
        dom = basis.r_sup * fac
        h = dom / n
        x = (np.arange(n) + 0.5) * h
        R, Z = np.meshgrid(x, x, indexing="ij")
        r = np.hypot(R, Z)
        sth, cth = R / r, Z / r
        A = basis.afield(r, cth, np.asarray(rep["a"]))
        F = fint(r * np.exp(A)) if A is not None else fint(r)
        T = basis.bfield(r, sth, cth, np.asarray(rep["b"]))
        return x, F, fint(r), (T if T is not None else np.zeros_like(r))

    def mirror(A):
        return np.concatenate([A[:, ::-1], A], axis=1)

    # (a) F contours at the R1 saturated-core optimum
    a1 = axes[0, 0]
    style_ax(a1)
    x, F, F0, dTh = fields(fint1, basis1, repG1)
    zz = np.concatenate([-x[::-1], x])
    lev = [np.pi / 6, np.pi / 3, np.pi / 2, 2 * np.pi / 3, 5 * np.pi / 6]
    a1.contour(x, zz, mirror(F0).T, levels=lev, colors=RED, linewidths=1.0,
               linestyles="--")
    a1.contour(x, zz, mirror(F).T, levels=lev, colors=BLUE, linewidths=1.6)
    a1.plot([], [], color=RED, ls="--", lw=1.0, label="hedgehog compacton")
    a1.plot([], [], color=BLUE, lw=1.6, label="R1 core optimum")
    a1.set_aspect("equal")
    a1.legend(frameon=False, fontsize=8, loc="upper right")
    a1.set_xlabel(r"$\rho$", color=INK)
    a1.set_ylabel(r"$z$", color=INK)
    a1.set_title("(a) $F$ contours, R1 ($t=0$) core optimum ($d=1$ frame)",
                 fontsize=10, color=INK)

    # (b) tilt field
    a2 = axes[0, 1]
    style_ax(a2)
    div = LinearSegmentedColormap.from_list(
        "div", ["#104281", "#3987e5", "#9ec5f4", "#f5f5f2", "#f2a09f",
                "#e34948", "#8f1d1c"])
    m = max(float(np.max(np.abs(dTh))), 1e-6)
    im = a2.pcolormesh(x, zz, mirror(dTh).T, cmap=div, vmin=-m, vmax=m,
                       shading="auto")
    a2.contour(x, zz, mirror(F).T, levels=[np.pi / 2], colors=INK,
               linewidths=0.8)
    cb = fig.colorbar(im, ax=a2, shrink=0.9)
    cb.set_label(r"$\Delta\Theta=\Theta-\theta$ (rad); $>0$ = equatorward",
                 fontsize=8, color=INK)
    cb.ax.tick_params(labelsize=7, colors=INK)
    a2.set_aspect("equal")
    a2.set_xlabel(r"$\rho$", color=INK)
    a2.set_ylabel(r"$z$", color=INK)
    a2.set_title("(b) internal-direction tilt at the R1 core optimum\n"
                 "(black: $F=\\pi/2$ shell)", fontsize=10, color=INK)

    # (c) optimal coefficients
    a3 = axes[0, 2]
    style_ax(a3)
    lab = basis1.labels
    v1 = np.asarray(repG1["a"] + repG1["b"])
    v2 = np.asarray(repG2["a"] + repG2["b"])
    ix = np.arange(len(lab))
    cols = [BLUE] * basis1.na + [AQUA] * basis1.nb
    a3.bar(ix - 0.2, v1, width=0.38, color=cols)
    a3.bar(ix + 0.2, v2, width=0.38, color=cols, alpha=0.45)
    a3.axhline(0.0, color=INK, lw=0.8)
    a3.set_xticks(ix, lab, rotation=60, fontsize=7)
    from matplotlib.patches import Patch
    a3.legend(handles=[Patch(color=BLUE, label="contour modes $a_{kl}$"),
                       Patch(color=AQUA, label="tilt modes $b_{kl}$"),
                       Patch(color="#999999", alpha=0.45,
                             label="lighter: R2 ($\\epsilon=0.05$)")],
              frameon=False, fontsize=8)
    a3.set_ylabel("coefficient", color=INK)
    a3.set_title("(c) saturated-core mode amplitudes (solid: R1)",
                 fontsize=10, color=INK)

    # (d) marginal spectrum at the compacton
    a4 = axes[1, 0]
    style_ax(a4)
    for r in sp_comp["rows"]:
        is_b = r["mode"].startswith("b")
        a4.plot(r["d2E"], abs(r["dI"]), marker="o" if is_b else "s", ms=8,
                mec="white", mew=0.8, ls="none",
                color=AQUA if is_b else BLUE)
        a4.annotate(r["mode"], (r["d2E"], abs(r["dI"])), xytext=(5, 4),
                    textcoords="offset points", fontsize=7, color=INK)
    xg = np.linspace(0.0, max(r["d2E"] for r in sp_comp["rows"]) * 1.05, 50)
    a4.plot([], [], "s", color=BLUE, label="contour modes (a)")
    a4.plot([], [], "o", color=AQUA, label="tilt modes (b)")
    a4.legend(frameon=False, fontsize=8, loc="upper left")
    a4.set_xlabel(r"quadratic cost $\partial^2 E_{\rm stat}/\partial p^2$",
                  color=INK)
    a4.set_ylabel(r"inertia gain $|\partial I/\partial p|$", color=INK)
    a4.set_title("(d) marginal mode spectrum at the compacton ($t=0$):\n"
                 "every mode buys inertia at finite quadratic cost",
                 fontsize=10, color=INK)

    # (e) boundary scan: c
    a5 = axes[1, 1]
    style_ax(a5)
    sl = [row["amaxcap"] for row in scan]
    cv = [row["c_paper"] for row in scan]
    a5.plot(sl, cv, "-o", color=BLUE, lw=2, ms=7, mec="white")
    a5.axhline(4.0 / math.pi * gfam, color=AQUA, lw=1.4, ls="-")
    a5.text(sl[0], 4.0 / math.pi * gfam,
            " saturated closure $(4/\\pi)G^*_{\\rm fam}$", color=AQUA,
            fontsize=8, va="bottom")
    a5.axhline(C_PAPER_LOWER, color=INK, lw=1.2, ls="-.")
    a5.text(sl[0], C_PAPER_LOWER, " rigorous bound $2\\sqrt{2}$", color=INK,
            fontsize=8, va="bottom")
    a5.axhline(2.5147, color=VIOLET, lw=1.0, ls="--")
    a5.text(sl[-1], 2.5147, " paper 2.5147", color=VIOLET, fontsize=8,
            va="bottom", ha="right")
    a5.axhline(2.0533, color=RED, lw=1.0, ls=":")
    a5.text(sl[-1], 2.0533, " scaling rung 2.0533", color=RED, fontsize=8,
            va="bottom", ha="right")
    a5.set_xlabel("family size (cap on $\\ln r_{\\rm eff}/r$; "
                  "stretch wall dialled jointly)", color=INK)
    a5.set_ylabel("$c_{\\rm paper}$ at the direct closure", color=INK)
    a5.set_title("(e) direct closure vs family size (R1, $t=0$):\n"
                 "no interior optimum — runs toward halo saturation",
                 fontsize=10, color=INK)

    # (f) boundary scan: kappa
    a6 = axes[1, 2]
    style_ax(a6)
    kv = [row["kappa_paper"] for row in scan]
    a6.plot(sl, kv, "-o", color=BLUE, lw=2, ms=7, mec="white")
    a6.axhline(KAPPA_CRIT_PAPER, color=AQUA, lw=1.4)
    a6.text(sl[0], KAPPA_CRIT_PAPER, " saturation $1/\\sqrt{2} = 0.70711$",
            color=AQUA, fontsize=8, va="bottom")
    a6.axhline(0.7638, color=VIOLET, lw=1.0, ls="--")
    a6.text(sl[-1], 0.7638, " paper $\\sqrt{7/12}$", color=VIOLET,
            fontsize=8, va="bottom", ha="right")
    a6.axhline(0.802, color=RED, lw=1.0, ls=":")
    a6.text(sl[-1], 0.802, " $\\langle r_1\\rangle$ 0.802", color=RED,
            fontsize=8, va="bottom", ha="right")
    a6.set_xlabel("family size (cap on $\\ln r_{\\rm eff}/r$)",
                  color=INK)
    a6.set_ylabel("$\\kappa_{\\rm paper}$", color=INK)
    a6.set_title("(f) $\\kappa$ falls toward the isorotation mass\n"
                 "threshold $1/\\sqrt{2}$ as the family grows", fontsize=10,
                 color=INK)

    fig.suptitle("Tier-2b step 3: enlarged-basis (contour + tilt) closure — "
                 "F-R4 adjudicated by halo saturation", fontsize=12,
                 color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
# json helper
# ----------------------------------------------------------------------
def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return o


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    out = dict(config=dict(N_OPT=N_OPT, N_POL=N_POL, N_FIN=N_FIN,
                           N_COARSE=N_COARSE, ABOX=ABOX, BBOX=BBOX,
                           AMAXCAP=AMAXCAP,
                           budgets=[B_TRI, B_COLD, B_WARM, B_POL, B_SCAN],
                           quick=QUICK, t_frozen=TFROZEN),
               theory=dict(kappa_crit_ours=KAPPA_CRIT_OURS,
                           kappa_crit_paper=KAPPA_CRIT_PAPER,
                           G_bound=G_BOUND, c_paper_lower=C_PAPER_LOWER,
                           statement="halo saturation: kappa_paper=1/sqrt2,"
                           " c_paper=(4/pi) min(E_static - I/(16 pi))"))

    # ------------------------------------------------------------------
    print("== GATES ==")
    gates = {}
    with open(OUTDIR + "/axi_results.json") as fh:
        A2 = json.load(fh)
    g0 = {}
    for k in ("solA", "solB", "endA", "endB"):
        s = A2[k]
        g0[k] = math.sqrt((2.0 / 3.0) * s["I"] * s["Estat"]) / s["L"] - 1.0
    gates["G0_clock_identity_vs_step2"] = g0
    g0max = max(abs(v) for v in g0.values())
    print("  G0 clock closed form L=sqrt(2/3 I E) vs step-2 solutions: "
          "max rel %.2e  %s" % (g0max, "PASS" if g0max < 1e-4 else "FAIL"))

    # G0b: Bogomolny normalisation for the G bound (verified vs step-1)
    fq = np.linspace(0.0, math.pi, 200001)
    bps_num = 2.0 * float(np.trapezoid(np.sin(fq) ** 2
                                       * np.sqrt(1.0 - np.cos(fq)), fq))
    gb_num = 2.0 * float(np.trapezoid(np.sin(fq) ** 2 * math.sqrt(0.5)
                                      * (1.0 - np.cos(fq)), fq))
    g0b = dict(bps_vs_known=bps_num / (32.0 * math.sqrt(2.0) / 15.0) - 1.0,
               gbound_vs_closed=gb_num / G_BOUND - 1.0)
    gates["G0b_bound_normalisation"] = g0b
    print("  G0b Bogomolny normalisation: BPS %+.1e vs 32sqrt2/15; "
          "G_bound quadrature %+.1e vs pi/sqrt2 = %.6f  %s"
          % (g0b["bps_vs_known"], g0b["gbound_vs_closed"], G_BOUND,
             "PASS" if max(abs(v) for v in g0b.values()) < 1e-6 else "FAIL"))

    with open(OUTDIR + "/radial_results.json") as fh:
        RJ = json.load(fh)
    targets = {k: RJ["N_run"][k] for k in ("E2", "E4", "E6", "E0")}
    r1d, f1d, sect1d = ax.radial_profile(TFROZEN)
    fspl = ax.CubicSpline1D(r1d, f1d)
    rsup2 = ax.support_radius(r1d, f1d)
    rsup1 = RSTAR * 1.02

    DOMFAC0 = 1.5
    SLIM0 = math.log(DOMFAC0) - 0.03
    basis1 = GenBasis(rsup1)
    basis2 = GenBasis(rsup2)
    dom1, dom2 = rsup1 * DOMFAC0, rsup2 * DOMFAC0
    z6, z6b = np.zeros(6), np.zeros(6)

    P0h = ax.sector_pieces(make_qfun(fspl, basis2, z6, z6b), N_FIN, dom2,
                           dom2)
    g1 = {k: P0h[k] / targets[k] - 1.0 for k in ("E2", "E4", "E6", "E0")}
    Irad = ax.radial_inertia(r1d, f1d)
    g2 = P0h["I"] / Irad - 1.0
    gates["G1_sectors_vs_radial"] = g1
    gates["G1_degree"] = P0h["BDEG"]
    gates["G2_inertia_vs_radial"] = g2
    ok1 = (max(abs(v) for v in g1.values()) < 5e-3
           and abs(P0h["BDEG"] - 1) < 5e-3)
    print("  G1 p=0 sectors vs step-1 radial (N=%d): %s  degree %.5f  %s"
          % (N_FIN, "  ".join("%s %+.1e" % kv for kv in g1.items()),
             P0h["BDEG"], "PASS" if ok1 else "FAIL"))
    print("  G2 p=0 inertia vs radial: rel %+.2e  %s"
          % (g2, "PASS" if abs(g2) < 5e-3 else "FAIL"))

    P0c = ax.sector_pieces(make_qfun(ax.f_compacton, basis1, z6, z6b),
                           N_FIN, dom1, dom1)
    g3 = dict(E6=P0c["E6"] / E6_COMPACTON - 1.0,
              E0=P0c["E0"] / E6_COMPACTON - 1.0,
              I=P0c["I"] / I_COMPACTON - 1.0, degree=P0c["BDEG"])
    gates["G3_compacton_vs_exact"] = g3
    ok3 = max(abs(g3[k]) for k in ("E6", "E0", "I")) < 2e-3
    print("  G3 compacton p=0 vs exact: E6 %+.1e  E0 %+.1e  I %+.1e  "
          "degree %.5f  %s" % (g3["E6"], g3["E0"], g3["I"], g3["degree"],
                               "PASS" if ok3 else "FAIL"))

    at = np.array([0.08, -0.05, 0.06, -0.04, 0.05, 0.03])
    bt = np.array([0.30, -0.15, 0.22, -0.10, 0.15, 0.10])
    Pt = ax.sector_pieces(make_qfun(ax.f_compacton, basis1, at, bt), N_FIN,
                          dom1, dom1)
    E0s, Is = sph_check(ax.f_compacton, basis1, at, bt, dom1)
    g4 = dict(E0=Pt["E0"] / E0s - 1.0, I=Pt["I"] / Is - 1.0,
              degree=Pt["BDEG"])
    gates["G4_spherical_crosscheck"] = g4
    ok4 = (abs(g4["E0"]) < 2e-3 and abs(g4["I"]) < 2e-3
           and abs(Pt["BDEG"] - 1.0) < 5e-3)
    print("  G4 deformed test config vs independent spherical quadrature: "
          "E0 %+.1e  I %+.1e  degree %.5f  %s"
          % (g4["E0"], g4["I"], g4["degree"], "PASS" if ok4 else "FAIL"))

    Lh0, _, _, _, _ = dl_fixed_point(P0h, TFROZEN, 5.0)
    g5 = dict(c_paper=2.0 * Lh0 * CMAP_L,
              vs_solA=2.0 * Lh0 * CMAP_L / A2["solA"]["c_paper"] - 1.0)
    gates["G5_hedgehog_closure_vs_solA"] = g5
    print("  G5 p=0 closure at eps=0.05: c_paper %.5f vs step-2 solA %.5f "
          "(rel %+.1e)  %s" % (g5["c_paper"], A2["solA"]["c_paper"],
                               g5["vs_solA"],
                               "PASS" if abs(g5["vs_solA"]) < 3e-3
                               else "FAIL"))
    out["gates"] = gates
    if not (g0max < 1e-4 and ok1 and abs(g2) < 5e-3 and ok3 and ok4
            and abs(g5["vs_solA"]) < 3e-3):
        print("!! GATE FAILURE — aborting.")
        with open(OUTDIR + "/axi3_results.json", "w") as fh:
            json.dump(jsonable(out), fh, indent=2)
        sys.exit(1)

    # ------------------------------------------------------------------
    # THEORY banner (the derivation the numerics will verify)
    # ------------------------------------------------------------------
    print("\n== HALO THRESHOLD (derived; verified numerically below) ==")
    print("  dR(tilted halo) = [INT eta^2 dV] (1/(8pi) - kappa_ours^2)")
    print("  kappa_crit_ours = 1/sqrt(8pi) = %.5f;  kappa_crit_paper = "
          "1/sqrt2 = %.5f" % (KAPPA_CRIT_OURS, KAPPA_CRIT_PAPER))
    for k in ("solA", "solB", "endA", "endB"):
        print("    step-2 %s: kappa_ours = %.5f  -> halo-UNSTABLE (%.2fx "
              "critical)" % (k, A2[k]["kappa_ours"],
                             A2[k]["kappa_ours"] / KAPPA_CRIT_OURS))
    print("  saturated closure: kappa_paper = 1/sqrt2 EXACT, c_paper = "
          "(4/pi) G*, G* = min(E_static - I/16pi)")
    print("  rigorous bound: G >= pi/sqrt2 -> c_paper >= 2 sqrt2 = %.5f"
          % C_PAPER_LOWER)

    # ------------------------------------------------------------------
    # engine-level halo verification
    # ------------------------------------------------------------------
    print("\n== HALO PROBE (explicit shell, full engine) ==")
    hp1 = halo_probe(ax.f_compacton, rsup1, 0.0,
                     A2["endA"]["kappa_ours"], "t=0    vs endA ")
    hp2 = halo_probe(fspl, rsup2, TFROZEN,
                     A2["solA"]["kappa_ours"], "t=eps  vs solA ")
    out["halo_probe"] = dict(t0=hp1, teps=hp2)

    # ------------------------------------------------------------------
    # R1: t = 0 direct closure (task-literal), 2x2 ablation
    # ------------------------------------------------------------------
    print("\n== R1: t = 0 DIRECT closure (compacton base), 2x2 ablation ==")
    alpha07, alpha045 = bseed_alpha2(0.7), bseed_alpha2(0.45)
    P0c_cache = {N_FIN: P0c}
    n1a, n1b = basis1.na, basis1.nb
    act_a = list(range(n1a))
    act_b = list(range(n1a, n1a + n1b))
    act_j = act_a + act_b

    rep_hh, _, _ = run_closure("R1-none", ax.f_compacton, basis1, 0.0, [],
                               [np.zeros(0)], dom1, SLIM0, P0c_cache)
    seeds_a = [np.zeros(n1a),
               seed_vec(basis1, act_a, a_map={2: alpha07}),
               seed_vec(basis1, act_a, a_map={2: alpha045})]
    rep_a, x_a, _ = run_closure("R1-a", ax.f_compacton, basis1, 0.0, act_a,
                                seeds_a, dom1, SLIM0, P0c_cache)
    seeds_b = [np.zeros(n1b),
               seed_vec(basis1, act_b, b_map={1: 0.35}),
               seed_vec(basis1, act_b, b_map={1: 0.70})]
    rep_b, x_b, _ = run_closure("R1-b", ax.f_compacton, basis1, 0.0, act_b,
                                seeds_b, dom1, SLIM0, P0c_cache)
    seeds_j = [np.concatenate([x_a, x_b]),
               np.concatenate([x_a, np.zeros(n1b)]),
               np.concatenate([np.zeros(n1a), x_b]),
               seed_vec(basis1, act_j, a_map={2: alpha07},
                        b_map={1: 0.35}),
               np.zeros(n1a + n1b)]
    rep_j, x_j, obj_j = run_closure("R1-joint", ax.f_compacton, basis1, 0.0,
                                    act_j, seeds_j, dom1, SLIM0, P0c_cache)

    # fixed-L 2x2 Routhian table at L = L*_joint
    print("\n   fixed-L 2x2 ablation at L* = %.6f (R1 joint):" % rep_j["L"])
    Lj = rep_j["L"]
    _, Rn = min_d(P0c_cache[N_OPT], 0.0, Lj)
    tab = {"none": Rn}
    for name, act, xw in (("a-only", act_a, x_a), ("b-only", act_b, x_b),
                          ("joint", act_j, x_j)):
        o = Objective(ax.f_compacton, basis1, 0.0, N_OPT, dom1, act, SLIM0)
        o.L = Lj
        st = np.array([0.02 if i < n1a else 0.05 for i in act])
        _, fv, _ = nelder_mead(o, xw, st, maxfev=1200)
        tab[name] = fv
    for k in ("none", "a-only", "b-only", "joint"):
        print("      R(%-7s; L*) = %.6f   (dR vs none: %+.6f)"
              % (k, tab[k], tab[k] - Rn))
    out["R1_ablation_fixedL"] = dict(L=Lj, R=tab)

    # enlarged basis
    print("\n   enlarged-basis rerun (a_l={0,2,4}, b_l={1,2,3}; 18 params):")
    basis1x = GenBasis(rsup1, K=3, a_l=(0, 2, 4), b_l=(1, 2, 3))
    embed = np.zeros(basis1x.n)
    for i, m in enumerate(basis1.a_modes):
        embed[basis1x.a_modes.index(m)] = rep_j["a"][i]
    for i, m in enumerate(basis1.b_modes):
        embed[basis1x.na + basis1x.b_modes.index(m)] = rep_j["b"][i]
    P0x_cache = {}
    rep_x, x_x, _ = run_closure("R1-ext", ax.f_compacton, basis1x, 0.0,
                                list(range(basis1x.n)),
                                [embed, np.zeros(basis1x.n)], dom1, SLIM0,
                                P0x_cache)
    out["R1_direct"] = dict(none=rep_hh, a_only=rep_a, b_only=rep_b,
                            joint=rep_j, extended=rep_x)

    # ------------------------------------------------------------------
    # boundary scan (family size dial): the no-interior-optimum signal
    # ------------------------------------------------------------------
    print("\n== R1 boundary scan: direct closure vs family size "
          "(both walls dialled together) ==")
    scan = [dict(amaxcap=AMAXCAP, stretch_lim=SLIM0,
                 c_paper=rep_j["c_paper"],
                 kappa_paper=rep_j["kappa_paper"], V=rep_j["V"],
                 g=rep_j["g"], flags=rep_j["flags"])]
    for domfac, acap in ((1.30, 0.30), (1.75, 0.65)):
        slim = math.log(domfac) - 0.03
        cache = {}
        rp, _, _ = run_closure("R1-c%.2f" % acap, ax.f_compacton, basis1,
                               0.0, act_j, [x_j], rsup1 * domfac, slim,
                               cache, b_cold=B_SCAN, max_outer=8,
                               amaxcap=acap)
        scan.append(dict(amaxcap=acap, stretch_lim=slim,
                         c_paper=rp["c_paper"],
                         kappa_paper=rp["kappa_paper"], V=rp["V"],
                         g=rp["g"], flags=rp["flags"]))
    scan.sort(key=lambda s: s["amaxcap"])
    print("   %-8s %-12s %10s %12s %8s %8s" % ("amaxcap", "stretch_lim",
                                               "c_paper", "kappa_paper",
                                               "V", "g*"))
    for s in scan:
        print("   %-8.2f %-12.4f %10.4f %12.4f %8.3f %8.3f" % (
            s["amaxcap"], s["stretch_lim"], s["c_paper"],
            s["kappa_paper"], s["V"], s["g"]))
    out["R1_boundary_scan"] = scan

    # ------------------------------------------------------------------
    # R1 saturated closure (well-posed): minimise G, 2x2 ablation
    # ------------------------------------------------------------------
    print("\n== R1 SATURATED closure: G = E_static - I/(16 pi) minimised ==")
    gcacheC = {N_FIN: P0c}
    repG_hh, _, _ = run_gmin("G1-none", ax.f_compacton, basis1, 0.0, [],
                             [np.zeros(0)], dom1, SLIM0, gcacheC)
    repG_a, gx_a, _ = run_gmin("G1-a", ax.f_compacton, basis1, 0.0, act_a,
                               seeds_a, dom1, SLIM0, gcacheC)
    repG_b, gx_b, _ = run_gmin("G1-b", ax.f_compacton, basis1, 0.0, act_b,
                               seeds_b, dom1, SLIM0, gcacheC)
    seeds_gj = [np.concatenate([gx_a, gx_b]),
                np.concatenate([gx_a, np.zeros(n1b)]),
                np.concatenate([np.zeros(n1a), gx_b]),
                np.zeros(n1a + n1b)]
    repG_j, gx_j, _ = run_gmin("G1-joint", ax.f_compacton, basis1, 0.0,
                               act_j, seeds_gj, dom1, SLIM0, gcacheC)
    embedG = np.zeros(basis1x.n)
    for i, m in enumerate(basis1.a_modes):
        embedG[basis1x.a_modes.index(m)] = repG_j["a"][i]
    for i, m in enumerate(basis1.b_modes):
        embedG[basis1x.na + basis1x.b_modes.index(m)] = repG_j["b"][i]
    repG_x, _, _ = run_gmin("G1-ext", ax.f_compacton, basis1x, 0.0,
                            list(range(basis1x.n)),
                            [embedG, np.zeros(basis1x.n)], dom1, SLIM0,
                            P0x_cache)
    print("   G ablation: none %.6f | a %.6f | b %.6f | joint %.6f | "
          "ext %.6f | rigorous bound %.6f"
          % (repG_hh["G"], repG_a["G"], repG_b["G"], repG_j["G"],
             repG_x["G"], G_BOUND))
    out["R1_saturated"] = dict(none=repG_hh, a_only=repG_a, b_only=repG_b,
                               joint=repG_j, extended=repG_x,
                               G_bound=G_BOUND,
                               c_paper_bracket=[C_PAPER_LOWER,
                                                4.0 / math.pi * repG_x["G"]])

    # ------------------------------------------------------------------
    # R3: marginal-mode spectrum
    # ------------------------------------------------------------------
    print("\n== R3: marginal-mode spectrum (d2E_static/dp^2, dI/dp) ==")
    sp_comp = spectrum(ax.f_compacton, basis1, np.zeros(basis1.n), 0.0,
                       1.0, N_POL, dom1, SLIM0)
    print_spectrum(sp_comp, "at the COMPACTON (p = 0, d = 1):")
    p_j = np.asarray(rep_j["a"] + rep_j["b"])
    sp_opt = spectrum(ax.f_compacton, basis1, p_j, 0.0, rep_j["d"], N_POL,
                      dom1, SLIM0)
    print_spectrum(sp_opt, "at the R1 DIRECT-JOINT optimum (d* = %.4f):"
                   % rep_j["d"])
    print("   NOTE: at the compacton every dE/dp = 0 (BPS minimum) while "
          "dI/dp != 0 for tilt")
    print("   modes -> R has a strictly downhill direction whenever "
          "kappa > kappa_crit (F-R4's")
    print("   'free' SDiff dial is replaced by a PAID dial that the closure "
          "can always afford).")
    out["R3_spectrum_compacton"] = sp_comp
    out["R3_spectrum_R1optimum"] = sp_opt

    # ------------------------------------------------------------------
    # R2: eps = 0.05
    # ------------------------------------------------------------------
    print("\n== R2: eps = 0.05 (t frozen), DIRECT closure ==")
    P0h_cache = {N_FIN: P0h}
    seeds2 = [np.asarray(rep_j["a"] + rep_j["b"]),
              seed_vec(basis2, act_j, a_map={2: alpha07}, b_map={1: 0.35}),
              np.concatenate([np.zeros(n1a), np.asarray(rep_j["b"])]),
              np.zeros(basis2.n)]
    rep2, x2, _ = run_closure("R2-joint", fspl, basis2, TFROZEN, act_j,
                              seeds2, dom2, SLIM0, P0h_cache)

    print("\n== R2 SATURATED closure ==")
    repG2_hh, _, _ = run_gmin("G2-none", fspl, basis2, TFROZEN, [],
                              [np.zeros(0)], dom2, SLIM0, P0h_cache)
    seeds_g2 = [np.asarray(repG_j["a"] + repG_j["b"]),
                np.zeros(basis2.n),
                seed_vec(basis2, act_j, a_map={2: alpha07},
                         b_map={1: 0.35})]
    repG2, gx2, _ = run_gmin("G2-joint", fspl, basis2, TFROZEN, act_j,
                             seeds_g2, dom2, SLIM0, P0h_cache)
    out["R2"] = dict(direct=rep2, saturated=repG2,
                     saturated_hedgehog=repG2_hh, targets=TARGETS_R2)

    # ------------------------------------------------------------------
    # verdicts
    # ------------------------------------------------------------------
    print("\n== VERDICTS ==")
    print("  R1 (t -> 0), DIRECT closure over the enlarged family "
          "(boundary-limited):")
    print("   %-26s %9s %9s %8s %8s %8s  %s" % ("config", "c_paper",
          "kap_pap", "g*", "lam_eff", "V", "boundary"))
    for nm, rp in (("hedgehog (scaling rung)", rep_hh),
                   ("a-only (contour)", rep_a), ("b-only (tilt)", rep_b),
                   ("joint", rep_j), ("extended (18 modes)", rep_x)):
        print("   %-26s %9.4f %9.4f %8.4f %8.4f %8.3f  %s"
              % (nm, rp["c_paper"], rp["kappa_paper"], rp["g"],
                 rp["lam_eff"], rp["V"],
                 [k for k, v in rp["flags"].items() if v] or "-"))
    print("  R1 SATURATED closure (well-posed limit of the same problem):")
    print("   %-26s %9s %9s %8s %8s" % ("config", "G*", "c_paper",
                                        "g_tot", "V_E0"))
    for nm, rp in (("hedgehog core", repG_hh), ("a-only core", repG_a),
                   ("b-only core", repG_b), ("joint core", repG_j),
                   ("extended core", repG_x)):
        print("   %-26s %9.5f %9.4f %8.3f %8.3f"
              % (nm, rp["G"], rp["c_paper"], rp["g_total"], rp["V_E0_tot"]))
    print("   kappa_paper = 1/sqrt2 = %.5f EXACT at saturation"
          % KAPPA_CRIT_PAPER)
    print("   c_paper(t->0) bracket: [%.5f (rigorous), %.5f (family)]"
          % (C_PAPER_LOWER, 4.0 / math.pi * repG_x["G"]))
    print("   vs paper shape-exact endpoint (2.5147, 0.7638, g=3/2): "
          "c EXCLUDED from below by 2sqrt2 = %.5f; kappa != 1/sqrt2"
          % C_PAPER_LOWER)
    print("   vs scaling rung (2.0533, 0.9354, g=1): left behind as soon "
          "as tilt modes exist")

    print("\n  R2 (eps = 0.05) vs <r1> benchmark "
          "(c 2.37+/-0.09, kap 0.802+/-0.018, g 1.31+/-0.04, "
          "V 1.409+/-0.010):")
    tt = TARGETS_R2
    print("   direct (boundary-limited): c=%.4f (%+.1f sig)  kap=%.4f "
          "(%+.1f sig)  g=%.3f  V=%.3f  V_E0=%.3f"
          % (rep2["c_paper"], (rep2["c_paper"] - tt["c"][0]) / tt["c"][1],
             rep2["kappa_paper"],
             (rep2["kappa_paper"] - tt["kappa"][0]) / tt["kappa"][1],
             rep2["g"], rep2["V"], rep2["V_E0"]))
    print("   saturated (unrestricted limit): c=%.4f (%+.1f sig)  "
          "kap=%.5f (%+.1f sig)  g_tot=%.3f (%+.1f sig)  V_E0=%.3f "
          "(%+.1f sig)"
          % (repG2["c_paper"],
             (repG2["c_paper"] - tt["c"][0]) / tt["c"][1],
             KAPPA_CRIT_PAPER,
             (KAPPA_CRIT_PAPER - tt["kappa"][0]) / tt["kappa"][1],
             repG2["g_total"], (repG2["g_total"] - tt["g"][0]) / tt["g"][1],
             repG2["V_E0_tot"],
             (repG2["V_E0_tot"] - tt["V"][0]) / tt["V"][1]))
    print("   -> the <r1> tuple is NOT the unrestricted-closure solution: "
          "it sits at a halo-unstable saddle.")

    # ------------------------------------------------------------------
    out["runtime_s"] = time.time() - t_start
    with open(OUTDIR + "/axi3_results.json", "w") as fh:
        json.dump(jsonable(out), fh, indent=2)
    fig_modes(ax.f_compacton, basis1, repG_j, repG2, sp_comp, scan,
              repG_x["G"], OUTDIR + "/axi3_modes.png")
    print("\nruntime: %.1f s" % out["runtime_s"])
    print("wrote axi3_results.json, axi3_modes.png")


if __name__ == "__main__":
    main()
