#!/usr/bin/env python3
"""
axi4_locked_solve.py — Tier-4C: the DIRECTION-LOCKED (hedgehog-locked) closure.

Hypothesis under test (axi3_ADJUDICATION.md §3): the corpus's frozen <r1>
benchmark code plausibly direction-locks the field — Theta(r,theta) == theta
exactly (hedgehog direction n-hat fixed), only the profile F(r,theta) free.
Under that lock the tilt-halo channel of Step 3 is removed and a weaker
F-halo channel remains; the benchmark may then be a genuine constrained
minimum.  This script builds the locked problem honestly and measures it.

----------------------------------------------------------------------------
LOCKED FIELD AND ITS HALO CHANNEL(S) — derivation done independently here
----------------------------------------------------------------------------
q = (cos F, sin F sin(theta), 0, sin F cos(theta)) on the phi = 0 half-plane
(standard unit winding in phi).  Add a far-field profile halo
F = eta * chi(r) * s(theta) with the core untouched (compact support => the
regions are disjoint).  Exact small-eta energetics in our conventions
(E0 = (1/4pi) INT (1-cos F) dV,  I = INT 2 sin^2 F sin^2 theta dV):

    dE0 = (1/8pi) INT eta^2 chi^2 s^2 dV               [1 - cos x ~ x^2/2]
    dI  = 2 <sin^2 theta>_w INT eta^2 chi^2 s^2 dV     [w = s^2 weight]
    dE2, dE4: the radial-gradient and sin^2F/r^2 pieces vanish per unit halo
      norm as the shell radius grows (t-weighted anyway); dE6 = O(eta^4).
    dR = dE0 - (kappa_ours^2/2) dI
       = [INT eta^2 chi^2 s^2 dV] * ( 1/(8pi) - kappa_ours^2 <sin^2 th>_w )

  =>  kappa_crit_ours(s) = 1/sqrt(8 pi <sin^2 theta>_w)
      kappa_crit_paper(s) = KMAP * that = 1/sqrt(2 <sin^2 theta>_w).

Angle-uniform halo s = 1:  <sin^2 theta> = 2/3
      kappa_crit_ours = sqrt(3/(16 pi)) = 0.2443014
      kappa_crit_paper = sqrt(3)/2      = 0.8660254   (the task's value  OK)

FLAGGED DISCREPANCY vs the single-threshold expectation: the lock freezes
Theta but NOT the angular shape of F, so the locked problem retains a
one-parameter FAMILY of halo channels labelled by the angular profile s:
      s = 1        -> kappa_paper = sqrt(3)/2      = 0.86603
      s = sin th   -> 1/sqrt(8/5)                  = 0.79057
      s = sin^2 th -> 1/sqrt(12/7) = sqrt(7/12)    = 0.76376  (= deep-BPS k0!)
      s = sin^4 th -> 1/sqrt(20/11)                = 0.74162
      s^2 = |cos|  -> 1/sqrt(1)                    = 1.00000  (exact)
      s = |cos th| -> 1/sqrt(4/5)                  = 1.11803
      equatorial-ring limit s^2 -> delta(th-pi/2)  -> 1/sqrt2 = 0.70711
        (the Step-3 tilt value: the lock does not remove it, it only makes
         it harder to reach — the descent needs equator-concentrated far
         halos instead of arbitrary tilted ones).
At t = 0 the angular concentration is exactly free at O(eta^2) (E0 and I are
pointwise densities); at t > 0 the angular-gradient cost per unit halo norm
scales as 1/r_c^2 and vanishes in the far-field limit.  All of this is
MEASURED below with explicit engine-level shell probes (two amplitudes,
quadratic-scaling checks, radius scans).

----------------------------------------------------------------------------
LOCKED CLOSURE
----------------------------------------------------------------------------
Basis (task 4C): F(r,theta) = f_base(r * exp(SUM_kl a_kl B_k(r) C_l(mu))),
K = 4 (default) or 5 (ladder) cubic-B-spline bumps x even Legendre {1,P2,P4}
=> 12 / 15 parameters + analytic dilation d.  Theta == theta ALWAYS (the
lock).  Fixed-L Routhian minimisation (deterministic Nelder-Mead) + clock
iteration L^2 = (2/3) I E_static to |dL/L| <= 1e-8, warm-started; caps are
soft quartic walls and the ACCEPTANCE criterion is that no optimum sits
within 80% of any cap — if one does, all caps are widened (ladder below) and
the solve repeats.  Everything else (quadrature engine, sector constants,
unit map, clock) is REUSED unchanged from axi_solve.py / axi3_solve.py.

Deterministic, no RNG.  Python3 + numpy only.
Outputs: axi4_locked_results.json, axi4_locked_landscape.png,
axi4_locked_epsscan.png (log via tee -> axi4_locked_run.log).
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
import axi_solve as ax                    # Step-2 engine (REUSED, unmodified)
import axi3_solve as a3                   # Step-3 machinery (REUSED, unmodified)

# ----------------------------------------------------------------------
# frozen conventions (identical to Steps 1-3; NOTHING fitted)
# ----------------------------------------------------------------------
TFROZEN = ax.TFROZEN
RSTAR = ax.RSTAR
CMAP_L = ax.CMAP_L                        # c_paper     = 2 L * CMAP_L
KMAP = ax.KMAP                            # kappa_paper = (L/I) * KMAP = 2 sqrt(pi) k_ours
E6_COMPACTON = ax.E6_COMPACTON
I_COMPACTON = ax.I_COMPACTON

KLOCK_OURS = math.sqrt(3.0 / (16.0 * math.pi))       # uniform locked F-halo
KLOCK_PAPER = KLOCK_OURS * KMAP                      # = sqrt(3)/2
KRING_OURS = 1.0 / math.sqrt(8.0 * math.pi)          # ring limit (= tilt value)
KRING_PAPER = KRING_OURS * KMAP                      # = 1/sqrt2

TARGETS = dict(c=(2.37, 0.09), kappa=(0.802, 0.018), g=(1.31, 0.04),
               lam=(0.42, 0.05), V=(1.409, 0.010), erot=0.2500)
CORPUS_ONSET = (1.000, 0.004)                        # App. I.4 emission onset
CORPUS_TUPLES = [("scaling rung", 0.9354), ("<r1> benchmark", 0.802),
                 ("deep-BPS endpoint", 0.7638)]

QUICK = "--quick" in sys.argv

if QUICK:
    N_POL, N_FIN, N_COARSE, N_PROBE = 144, 192, 128, 144
    CAP_LEVELS = [dict(name="L0", domfac=1.5, amax=0.45, abox=0.80, N_opt=96),
                  dict(name="L1", domfac=2.0, amax=0.80, abox=1.40, N_opt=112)]
    BUD = dict(tri=60, cold=260, warm=90, pol=60)
    BUD_SCAN = dict(tri=60, cold=140, warm=70, pol=50)
    MAX_OUTER, MAX_OUTER_SCAN = 5, 4
    EPS_LIST = [0.02, 0.05, 0.10]
else:
    N_POL, N_FIN, N_COARSE, N_PROBE = 400, 720, 360, 480
    CAP_LEVELS = [dict(name="L0", domfac=1.5, amax=0.45, abox=0.80, N_opt=200),
                  dict(name="L1", domfac=2.0, amax=0.80, abox=1.40, N_opt=240),
                  dict(name="L2", domfac=2.7, amax=1.40, abox=2.40, N_opt=288)]
    BUD = dict(tri=260, cold=2200, warm=420, pol=170)
    BUD_SCAN = dict(tri=180, cold=1100, warm=340, pol=140)
    MAX_OUTER, MAX_OUTER_SCAN = 16, 12
    EPS_LIST = [0.02, 0.035, 0.05, 0.075, 0.1]

for lv in CAP_LEVELS:
    lv["stretch_lim"] = math.log(lv["domfac"]) - 0.03


# ----------------------------------------------------------------------
# basis: axi3 contour machinery, richer (K = 4/5, {1, P2, P4}), no tilt
# ----------------------------------------------------------------------
class LockedBasis(a3.GenBasis):
    """a-modes only (Theta == theta is the lock).  The dense grid used for
    the stretch/shrink extremes is extended to cover the FULL bump support
    (a3's grid stopped at 1.05 r_sup while bump tails reach (1+1.5/K) r_sup)."""

    def __init__(self, r_sup, K=4, a_l=(0, 2, 4)):
        a3.GenBasis.__init__(self, r_sup, K=K, a_l=a_l, b_l=())
        span = r_sup * (1.0 + 1.5 / K) * 1.02
        self.rd = np.linspace(0.0, span, 801)
        self.Bd = np.array([self.bump(k, self.rd) for k in range(K)])
        self._mud = np.linspace(-1.0, 1.0, 161)
        self.Cd = {l: a3.legC(l, self._mud)
                   for l in sorted(set(l for _, l in self.a_modes))}

    def a_extremes(self, a):
        A = np.zeros((self.rd.size, self._mud.size))
        for (k, l), coef in zip(self.a_modes, a):
            if coef != 0.0:
                A += coef * self.Bd[k][:, None] * self.Cd[l][None, :]
        return float(A.min()), float(A.max())


ZB = np.zeros(0)                          # the empty tilt vector (the lock)


def qfun_locked(fint, basis, a):
    return a3.make_qfun(fint, basis, np.asarray(a, float), ZB)


# ----------------------------------------------------------------------
# closure helpers (wider dilation bracket than a3's defaults)
# ----------------------------------------------------------------------
DLO, DHI = 0.30, 3.6


def MIN_D(P, t, L):
    w, R = ax.golden(lambda u: a3.R_from(P, t, L, math.exp(u))[0],
                     math.log(DLO), math.log(DHI), iters=60)
    return math.exp(w), R


def MIN_DG(P, t):
    w, G = ax.golden(lambda u: a3.G_from(P, t, math.exp(u)),
                     math.log(DLO), math.log(DHI), iters=60)
    return math.exp(w), G


def DLFP(P, t, L0):
    """Exact (d, L) clock fixed point at frozen sectors P."""
    L = max(L0, 1e-6)
    for _ in range(200):
        d, _ = MIN_D(P, t, L)
        _, Estat, I3 = a3.R_from(P, t, L, d)
        Ln = math.sqrt((2.0 / 3.0) * I3 * Estat)
        if abs(Ln - L) < 1e-14 * Ln:
            L = Ln
            break
        L = Ln
    d, _ = MIN_D(P, t, L)
    R, Estat, I3 = a3.R_from(P, t, L, d)
    return L, d, R, Estat, I3


# ----------------------------------------------------------------------
# objective with explicit, dial-able caps (soft quartic walls)
# ----------------------------------------------------------------------
class LockedObj:
    def __init__(self, fint, basis, t, N, dom, caps, mode="direct"):
        self.fint, self.basis, self.t = fint, basis, t
        self.N, self.dom = N, dom
        self.abox = caps["abox"]
        self.amaxcap = caps["amax"]
        self.stretch_lim = caps["stretch_lim"]
        self.mode = mode
        self.L = None
        self.nfev = 0

    def sectors(self, x, N=None):
        return ax.sector_pieces(qfun_locked(self.fint, self.basis, x),
                                self.N if N is None else N, self.dom, self.dom)

    def penalty(self, a):
        amin, amax = self.basis.a_extremes(a)
        if -amin > self.stretch_lim:          # hard: support beyond domain
            return None, -amin, amax
        pen = 400.0 * max(0.0, -amin - (self.stretch_lim - 0.02)) ** 2
        pen += 400.0 * max(0.0, amax - self.amaxcap) ** 2
        for v in a:
            pen += 200.0 * max(0.0, abs(v) - self.abox) ** 2
        return pen, -amin, amax

    def fracs(self, a):
        amin, amax = self.basis.a_extremes(a)
        return dict(abox=float(np.max(np.abs(a)) / self.abox),
                    amax=max(amax, 0.0) / self.amaxcap,
                    stretch=max(-amin, 0.0) / self.stretch_lim)

    def __call__(self, x):
        self.nfev += 1
        a = np.asarray(x, float)
        pen, stretch, _ = self.penalty(a)
        if pen is None:
            return 1.0e6 * (1.0 + stretch)
        P = self.sectors(a)
        if self.mode == "direct":
            _, val = MIN_D(P, self.t, self.L)
        else:
            _, val = MIN_DG(P, self.t)
        return val + pen


# ----------------------------------------------------------------------
# finalisation
# ----------------------------------------------------------------------
def finalize_locked(obj, x, L_guess, N, P0N):
    P = obj.sectors(x, N=N)
    L, d, R, Estat, I3 = DLFP(P, obj.t, L_guess)
    d0, _ = MIN_D(P0N, obj.t, 0.0)
    a = np.asarray(x, float)
    _, _, _, lam_eff = a3.moments(qfun_locked(obj.fint, obj.basis, a), N,
                                  obj.dom)
    fr = obj.fracs(a)
    maxfrac = max(fr.values())
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
               c_ours=2.0 * L, c_paper=2.0 * L * CMAP_L,
               kappa_ours=L / I3, kappa_paper=(L / I3) * KMAP,
               cap_fracs=fr, max_cap_frac=maxfrac,
               interior=bool(maxfrac <= 0.80),
               a=[float(v) for v in a])
    return rep


def finalize_locked_G(obj, x, N, P0N):
    """Ring-saturated closure over the LOCKED family (the ideal limit of the
    equator-concentrated F-halo channel; same algebra as axi3 Section 4)."""
    P = obj.sectors(x, N=N)
    d, Gv = MIN_DG(P, obj.t)
    L = math.sqrt(8.0 * math.pi) * Gv
    I_tot = 8.0 * math.pi * Gv
    I_core = P["I"] * d ** 3
    _, Estat_core, _ = a3.R_from(P, obj.t, 0.0, d)
    d0, _ = MIN_D(P0N, obj.t, 0.0)
    a = np.asarray(x, float)
    fr = obj.fracs(a)
    rep = dict(N=N, G=Gv, d=d, d0=d0, L=L, c_ours=2.0 * L,
               c_paper=4.0 / math.pi * Gv,
               kappa_ours=KRING_OURS, kappa_paper=KRING_PAPER,
               Estat_tot=1.5 * Gv, I_tot=I_tot, R=2.0 * Gv, erot_frac=0.25,
               Estat_core=Estat_core, I_core=I_core,
               halo_D=(I_tot - I_core) / 2.0,
               consistency=(Estat_core + (I_tot - I_core) / (16.0 * math.pi))
               / (1.5 * Gv) - 1.0,
               g_total=I_tot / (P0N["I"] * d0 ** 3),
               g_core=P["I"] / P0N["I"],
               degree=P["BDEG"], cap_fracs=fr,
               max_cap_frac=max(fr.values()),
               interior=bool(max(fr.values()) <= 0.80),
               a=[float(v) for v in a])
    return rep


def uniform_sat_reference(P0, t):
    """Analytic reference: hedgehog core + ANGLE-UNIFORM saturated F-halo.
    At saturation of the uniform channel kappa_ours = KLOCK_OURS and the halo
    is neutral for G_L := E_static - (3/(32 pi)) I; core relaxes in d only.
    (The G_L minimisation over the full locked family is NOT well-posed —
    finite-amplitude equatorial rings drive it to -inf — so only this
    dilation-relaxed hedgehog-core fixed point is quoted.)"""
    def GL(d):
        _, Estat, I3 = a3.R_from(P0, t, 0.0, d)
        return Estat - 3.0 / (32.0 * math.pi) * I3
    u, _ = ax.golden(lambda w: GL(math.exp(w)), math.log(DLO), math.log(DHI),
                     iters=60)
    d = math.exp(u)
    G = GL(d)
    c_ours = 2.0 * G / KLOCK_OURS
    return dict(G=G, d=d, c_ours=c_ours, c_paper=c_ours * CMAP_L,
                kappa_ours=KLOCK_OURS, kappa_paper=KLOCK_PAPER)


# ----------------------------------------------------------------------
# drivers
# ----------------------------------------------------------------------
def run_locked(tag, fint, basis, t, dom, caps, N_opt, seeds, bud, P0_cache,
               max_outer=MAX_OUTER, log=True):
    """Task-literal direct locked closure: min_a R(a; L) + clock on L, to
    |dL/L| <= 1e-8 (outer), then N_POL polish repeated to the same tol."""
    t0 = time.time()
    obj = LockedObj(fint, basis, t, N_opt, dom, caps)
    steps0 = np.full(basis.n, 0.06)

    def P0(N):
        key = (N, round(dom, 9))
        if key not in P0_cache:
            P0_cache[key] = ax.sector_pieces(
                qfun_locked(fint, basis, np.zeros(basis.n)), N, dom, dom)
        return P0_cache[key]

    L, _, _, _, _ = DLFP(P0(N_opt), t, 5.0)
    x, hist, rel, outer_conv = None, [], 1.0, False
    for it in range(max_outer):
        obj.L = L
        if it == 0:
            best = None
            if len(seeds) > 1:
                for s in seeds:
                    xs, fv, _ = a3.nelder_mead(obj, np.asarray(s, float),
                                               steps0, maxfev=bud["tri"])
                    if best is None or fv < best[1]:
                        best = (xs, fv)
            else:
                s0 = np.asarray(seeds[0], float)
                best = (s0, obj(s0))
            x, fv, _ = a3.nelder_mead(obj, best[0], steps0,
                                      maxfev=bud["cold"])
            x, fv, _ = a3.nelder_mead(obj, x, steps0 * 0.25,
                                      maxfev=max(bud["cold"] // 3, 40))
        else:
            sc = max(0.02, min(0.2, 50.0 * rel))
            x, fv, _ = a3.nelder_mead(obj, x, steps0 * sc,
                                      maxfev=bud["warm"])
        P = obj.sectors(x)
        Lx, d, R, Estat, I3 = DLFP(P, t, L)
        rel = abs(Lx - L) / Lx
        hist.append(dict(it=it, L=L, L_new=Lx, rel=rel, R=fv,
                         c_paper=2.0 * Lx * CMAP_L,
                         kappa_paper=Lx / I3 * KMAP, nfev=obj.nfev))
        if log:
            print("   %-12s outer %2d: L %.8f -> %.8f (dL/L %.2e)  R=%.7f"
                  "  nfev=%d" % (tag, it, L, Lx, rel, fv, obj.nfev))
        L = Lx
        if rel < 1.0e-8 and it >= 3:
            outer_conv = True
            break

    objP = LockedObj(fint, basis, t, N_POL, dom, caps)
    pol_rel = None
    for k in range(4):
        objP.L = L
        x, fv, _ = a3.nelder_mead(objP, x, steps0 * 0.06, maxfev=bud["pol"])
        P = objP.sectors(x)
        Lx, d, R, Estat, I3 = DLFP(P, t, L)
        pol_rel = abs(Lx - L) / Lx
        L = Lx
        if pol_rel < 1.0e-8 and k >= 1:
            break
    if log:
        print("   %-12s polish (N=%d): L=%.8f  dL/L=%.2e  nfev=%d"
              % (tag, N_POL, L, pol_rel, objP.nfev))

    rep = finalize_locked(objP, x, L, N_FIN, P0(N_FIN))
    rep["coarse"] = finalize_locked(objP, x, L, N_COARSE, P0(N_COARSE))
    rep["conv_c_rel"] = abs(rep["coarse"]["c_paper"] / rep["c_paper"] - 1.0)
    rep["hist"] = hist
    rep["nfev"] = obj.nfev + objP.nfev
    rep["outer_converged"] = bool(outer_conv)
    rep["outer_last_rel"] = rel
    rep["polish_last_rel"] = pol_rel
    rep["caps"] = dict(caps)
    rep["N_opt"] = N_opt
    rep["runtime_s"] = time.time() - t0
    if log:
        fr = rep["cap_fracs"]
        print("   %-12s FINAL (N=%d): c_pap=%.5f kap_pap=%.5f g*=%.4f "
              "lam_eff=%.4f V=%.4f V_E0=%.4f Erot/E=%.6f B=%.5f "
              "[coarse c=%.5f]  (%.0f s)"
              % (tag, N_FIN, rep["c_paper"], rep["kappa_paper"], rep["g"],
                 rep["lam_eff"], rep["V"], rep["V_E0"], rep["erot_frac"],
                 rep["degree"], rep["coarse"]["c_paper"], rep["runtime_s"]))
        print("   %-12s caps: abox %.0f%%  amax %.0f%%  stretch %.0f%%  -> %s"
              % (tag, 100 * fr["abox"], 100 * fr["amax"],
                 100 * fr["stretch"],
                 "INTERIOR (all <= 80%)" if rep["interior"]
                 else "!! CAP-LIMITED (>= 80%)"))
    return rep, x


def run_locked_G(tag, fint, basis, t, dom, caps, N_opt, seeds, bud, P0_cache,
                 log=True):
    t0 = time.time()
    obj = LockedObj(fint, basis, t, N_opt, dom, caps, mode="G")
    steps0 = np.full(basis.n, 0.06)

    def P0(N):
        key = (N, round(dom, 9))
        if key not in P0_cache:
            P0_cache[key] = ax.sector_pieces(
                qfun_locked(fint, basis, np.zeros(basis.n)), N, dom, dom)
        return P0_cache[key]

    best = None
    if len(seeds) > 1:
        for s in seeds:
            xs, fv, _ = a3.nelder_mead(obj, np.asarray(s, float), steps0,
                                       maxfev=bud["tri"])
            if best is None or fv < best[1]:
                best = (xs, fv)
    else:
        s0 = np.asarray(seeds[0], float)
        best = (s0, obj(s0))
    x, fv, _ = a3.nelder_mead(obj, best[0], steps0, maxfev=bud["cold"])
    objP = LockedObj(fint, basis, t, N_POL, dom, caps, mode="G")
    x, fv, _ = a3.nelder_mead(objP, x, steps0 * 0.08, maxfev=bud["pol"])
    rep = finalize_locked_G(objP, x, N_FIN, P0(N_FIN))
    rep["coarse"] = finalize_locked_G(objP, x, N_COARSE, P0(N_COARSE))
    rep["conv_c_rel"] = abs(rep["coarse"]["c_paper"] / rep["c_paper"] - 1.0)
    rep["nfev"] = obj.nfev + objP.nfev
    rep["runtime_s"] = time.time() - t0
    if log:
        print("   %-12s G*=%.6f -> c_pap=%.5f (kap_pap=1/sqrt2)  g_core=%.4f"
              "  halo_D=%.3f  B=%.5f  caps max %.0f%%  (%.0f s)"
              % (tag, rep["G"], rep["c_paper"], rep["g_core"], rep["halo_D"],
                 rep["degree"], 100 * rep["max_cap_frac"], rep["runtime_s"]))
    return rep, x


# ----------------------------------------------------------------------
# seeds
# ----------------------------------------------------------------------
def spheroid_seed(basis, lam, scale=1.0 / 1.45):
    """Project A = ln(m/r) of the spheroidal-shell map onto the C_l set; the
    1/1.45 undoes the interior bump-partition sum (seed quality only)."""
    mu = np.linspace(-1.0, 1.0, 241)
    h = 0.5 * np.log(lam ** (2.0 / 3.0) * (1.0 - mu * mu)
                     + lam ** (-4.0 / 3.0) * mu * mu)
    ls = sorted(set(l for _, l in basis.a_modes))
    X = np.stack([a3.legC(l, mu) for l in ls], axis=1)
    coef, *_ = np.linalg.lstsq(X, h, rcond=None)
    amap = dict(zip(ls, coef))
    p = np.zeros(basis.n)
    for i, (k, l) in enumerate(basis.a_modes):
        p[i] = amap.get(l, 0.0) * scale
    return p


def skirt_seed(basis, amp=-0.15):
    p = np.zeros(basis.n)
    for i, (k, l) in enumerate(basis.a_modes):
        if k == basis.K - 1 and l == 0:
            p[i] = amp
    return p


def project_field(b_from, x_from, b_to):
    """Least-squares embed of the warp field A(r, mu) into another basis."""
    r = np.linspace(0.0, b_to.r_sup * (1.0 + 1.5 / b_to.K), 260)
    mu = np.linspace(-1.0, 1.0, 91)
    Rg, Mg = np.meshgrid(r, mu, indexing="ij")
    A = b_from.afield(Rg, Mg, np.asarray(x_from, float))
    if A is None:
        return np.zeros(b_to.n)
    cols = [(b_to.bump(k, Rg) * a3.legC(l, Mg)).ravel()
            for (k, l) in b_to.a_modes]
    coef, *_ = np.linalg.lstsq(np.stack(cols, axis=1), A.ravel(), rcond=None)
    return coef


# ----------------------------------------------------------------------
# halo probes (Task A / D): explicit shells on the full engine
# ----------------------------------------------------------------------
SPROFILES = [
    ("s=1 (uniform)", lambda s, c: np.ones_like(s)),
    ("s=sin", lambda s, c: np.abs(s)),
    ("s=sin^2", lambda s, c: s * s),
    ("s=sin^4", lambda s, c: s ** 4),
    ("s^2~|cos| (sm.)", lambda s, c: (c * c + 0.01) ** 0.25),
    ("s=|cos|", lambda s, c: np.abs(c)),
]


def sprofile_moment(sfun):
    u = np.linspace(-1.0, 1.0, 20001)
    s = np.sqrt(np.clip(1.0 - u * u, 0.0, None))
    w = sfun(s, u) ** 2
    return float(np.trapezoid(w * (1.0 - u * u), u) / np.trapezoid(w, u))


def halo_probe(fint, rsup, t, sfun, rc_fac, w_fac, N, etas=(0.02, 0.04),
               base_cache=None):
    rc, w = rc_fac * rsup, w_fac * rsup
    dom = (rc + 2.0 * w) * 1.02

    def mk(eta):
        def qf(R, Z):
            r = np.hypot(R, Z)
            xx = np.abs((r - rc) / w)
            chi = np.where(xx < 1.0, 1.0 - 1.5 * xx * xx * (1.0 - 0.5 * xx),
                           0.25 * np.clip(2.0 - xx, 0.0, None) ** 3)
            sth, cth = R / r, Z / r
            F = fint(r) + eta * chi * sfun(sth, cth)
            sF = np.sin(F)
            return np.cos(F), sF * sth, sF * cth, sF / r
        return qf

    key = (round(dom, 9), N)
    if base_cache is not None and key in base_cache:
        P00 = base_cache[key]
    else:
        P00 = ax.sector_pieces(mk(0.0), N, dom, dom)
        if base_cache is not None:
            base_cache[key] = P00
    _, E00, _ = a3.R_from(P00, t, 0.0, 1.0)
    dEs, dIs, kco = [], [], []
    for eta in etas:
        P = ax.sector_pieces(mk(eta), N, dom, dom)
        _, Es, _ = a3.R_from(P, t, 0.0, 1.0)
        dE, dI = Es - E00, P["I"] - P00["I"]
        dEs.append(dE)
        dIs.append(dI)
        kco.append(math.sqrt(2.0 * dE / dI) if dI > 0 else float("nan"))
    return dict(rc_fac=rc_fac, w_fac=w_fac, N=N, dom=dom,
                etas=list(etas), dE=dEs, dI=dIs,
                kcrit_ours=kco, kcrit_paper=[k * KMAP for k in kco],
                dI_ratio=dIs[1] / dIs[0] if dIs[0] != 0 else float("nan"),
                dE_ratio=dEs[1] / dEs[0] if dEs[0] != 0 else float("nan"))


def placement(kpap, name):
    ko = kpap / KMAP
    return dict(name=name, kappa_paper=kpap, kappa_ours=ko,
                x_lock=ko / KLOCK_OURS, x_ring=ko / KRING_OURS,
                req_conc=1.0 / (8.0 * math.pi * ko * ko),
                stable_uniform=bool(ko < KLOCK_OURS),
                stable_ring=bool(ko < KRING_OURS))


# ----------------------------------------------------------------------
# figures (house palette; validated for this campaign)
# ----------------------------------------------------------------------
INK, GRID_C = "#333333", "#dddddd"
BLUE, AQUA, RED, VIOLET = "#2a78d6", "#1baf7a", "#e34948", "#4a3aa7"


def style_ax(a):
    a.set_facecolor("white")
    a.grid(color=GRID_C, lw=0.6)
    for sp in ("top", "right"):
        a.spines[sp].set_visible(False)
    a.tick_params(colors=INK, labelsize=9)


def fig_landscape(fint, basis, rep, ladder, k5rep, probes_t0, probes_eps,
                  derived, path):
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.2), dpi=160)
    fig.patch.set_facecolor("white")

    # fields of the accepted optimum
    n, fac = 340, 1.30
    dom = basis.r_sup * fac
    h = dom / n
    xg = (np.arange(n) + 0.5) * h
    Rg, Zg = np.meshgrid(xg, xg, indexing="ij")
    r = np.hypot(Rg, Zg)
    A = basis.afield(r, Zg / r, np.asarray(rep["a"]))
    F = fint(r * np.exp(A)) if A is not None else fint(r)
    F0 = fint(r)
    Afield = A if A is not None else np.zeros_like(r)
    zz = np.concatenate([-xg[::-1], xg])

    def mirror(M):
        return np.concatenate([M[:, ::-1], M], axis=1)

    a1 = axes[0, 0]
    style_ax(a1)
    lev = [np.pi / 6, np.pi / 3, np.pi / 2, 2 * np.pi / 3, 5 * np.pi / 6]
    a1.contour(xg, zz, mirror(F0).T, levels=lev, colors=RED,
               linewidths=1.0, linestyles="--")
    a1.contour(xg, zz, mirror(F).T, levels=lev, colors=BLUE, linewidths=1.6)
    a1.plot([], [], color=RED, ls="--", lw=1.0, label="hedgehog profile")
    a1.plot([], [], color=BLUE, lw=1.6, label="locked optimum")
    a1.set_aspect("equal")
    a1.legend(frameon=False, fontsize=8, loc="upper right")
    a1.set_xlabel(r"$\rho$", color=INK)
    a1.set_ylabel(r"$z$", color=INK)
    a1.set_title("(a) $F$ contours, locked closure at $\\epsilon=0.05$ "
                 "($d=1$ frame)", fontsize=10, color=INK)

    a2 = axes[0, 1]
    style_ax(a2)
    div = LinearSegmentedColormap.from_list(
        "div", ["#104281", "#3987e5", "#9ec5f4", "#f5f5f2", "#f2a09f",
                "#e34948", "#8f1d1c"])
    m = max(float(np.max(np.abs(Afield))), 1e-6)
    im = a2.pcolormesh(xg, zz, mirror(Afield).T, cmap=div, vmin=-m, vmax=m,
                       shading="auto")
    a2.contour(xg, zz, mirror(F).T, levels=[np.pi / 2], colors=INK,
               linewidths=0.8)
    cb = fig.colorbar(im, ax=a2, shrink=0.9)
    cb.set_label(r"warp $A=\ln(r_{\rm eff}/r)$; $<0$ = local stretch",
                 fontsize=8, color=INK)
    cb.ax.tick_params(labelsize=7, colors=INK)
    a2.set_aspect("equal")
    a2.set_xlabel(r"$\rho$", color=INK)
    a2.set_ylabel(r"$z$", color=INK)
    a2.set_title("(b) contour-warp field of the locked optimum\n"
                 "(black: $F=\\pi/2$ shell)", fontsize=10, color=INK)

    # (c) onset ladder
    a3x = axes[0, 2]
    style_ax(a3x)
    xs = np.linspace(0.30, 1.0, 200)
    a3x.plot(xs, 1.0 / np.sqrt(2.0 * xs), color=INK, lw=1.4,
             label=r"derived $1/\sqrt{2\langle\sin^2\theta\rangle_w}$")
    for probes, col, mk, lab in ((probes_t0, BLUE, "s", "measured, $t=0$"),
                                 (probes_eps, AQUA, "o",
                                  "measured, $\\epsilon=0.05$")):
        px = [d["conc"] for d in probes]
        py = [d["kcrit_paper"] for d in probes]
        a3x.plot(px, py, mk, color=col, ms=8, mec="white", mew=0.8,
                 ls="none", label=lab)
    for y, lab, col in ((1.000, "corpus onset 1.000", VIOLET),
                        (0.9354, "rung 0.9354", INK),
                        (0.86603, r"$\sqrt3/2$ (uniform)", INK),
                        (0.802, r"$\langle r_1\rangle$ 0.802", RED),
                        (0.76376, r"$\sqrt{7/12}$ (deep-BPS)", INK),
                        (0.70711, r"$1/\sqrt2$ (ring limit)", INK)):
        a3x.axhline(y, color=col, lw=0.9,
                    ls="--" if col != INK else ":", alpha=0.85)
        a3x.text(1.0, y, " " + lab, color=col, fontsize=7, va="bottom",
                 ha="right")
    a3x.set_xlabel(r"halo angular concentration "
                   r"$\langle\sin^2\theta\rangle_{s^2}$", color=INK)
    a3x.set_ylabel(r"$\kappa_{\rm paper}$ at onset", color=INK)
    a3x.legend(frameon=False, fontsize=8, loc="upper right")
    a3x.set_title("(c) locked F-halo onset family: derivation vs engine",
                  fontsize=10, color=INK)

    # (d, e) cap/basis ladder
    labs = [row["label"] for row in ladder]
    cv = [row["c_paper"] for row in ladder]
    kv = [row["kappa_paper"] for row in ladder]
    ix = np.arange(len(labs))
    a4 = axes[1, 0]
    style_ax(a4)
    a4.plot(ix, cv, "-o", color=BLUE, lw=2, ms=8, mec="white")
    a4.axhspan(2.37 - 0.09, 2.37 + 0.09, color=RED, alpha=0.10, lw=0)
    a4.axhline(2.37, color=RED, lw=1.0, ls="--")
    a4.text(ix[0], 2.372, r" $\langle r_1\rangle$: $2.37\pm0.09$",
            color=RED, fontsize=8, va="bottom")
    a4.axhline(2.15824, color=INK, lw=0.9, ls=":")
    a4.text(ix[-1], 2.159, " hedgehog (solA) 2.158", color=INK, fontsize=8,
            va="bottom", ha="right")
    a4.set_xticks(ix, labs, fontsize=8)
    a4.set_ylabel(r"$c_{\rm paper}$", color=INK)
    a4.set_title("(d) locked closure vs cap level and basis size "
                 "($\\epsilon=0.05$)", fontsize=10, color=INK)

    a5 = axes[1, 1]
    style_ax(a5)
    a5.plot(ix, kv, "-o", color=BLUE, lw=2, ms=8, mec="white")
    a5.axhspan(0.802 - 0.018, 0.802 + 0.018, color=RED, alpha=0.10, lw=0)
    a5.axhline(0.802, color=RED, lw=1.0, ls="--")
    a5.text(ix[0], 0.803, r" $\langle r_1\rangle$: $0.802\pm0.018$",
            color=RED, fontsize=8, va="bottom")
    for y, lab in ((KLOCK_PAPER, r" $\sqrt3/2$ uniform-halo onset"),
                   (KRING_PAPER, r" $1/\sqrt2$ ring limit")):
        a5.axhline(y, color=AQUA, lw=1.2)
        a5.text(ix[0], y, lab, color=AQUA, fontsize=8, va="bottom")
    a5.set_xticks(ix, labs, fontsize=8)
    a5.set_ylabel(r"$\kappa_{\rm paper}$", color=INK)
    a5.set_title("(e) $\\kappa$ placement vs the locked thresholds",
                 fontsize=10, color=INK)

    # (f) mode amplitudes K=4 vs K=5
    a6 = axes[1, 2]
    style_ax(a6)
    lab4 = basis.labels
    v4 = np.asarray(rep["a"])
    a6.bar(np.arange(len(lab4)) - 0.18, v4, width=0.34, color=BLUE,
           label="K=4 (12 params)")
    if k5rep is not None:
        v5 = np.asarray(k5rep["a"])
        a6.bar(np.arange(len(v5)) * (len(lab4) - 1.0)
               / max(len(v5) - 1.0, 1.0) + 0.18, v5, width=0.34,
               color=AQUA, alpha=0.75, label="K=5 (15 params, rescaled axis)")
    a6.axhline(0.0, color=INK, lw=0.8)
    a6.set_xticks(np.arange(len(lab4)), lab4, rotation=60, fontsize=7)
    a6.legend(frameon=False, fontsize=8)
    a6.set_ylabel("coefficient", color=INK)
    a6.set_title("(f) locked-optimum mode amplitudes", fontsize=10,
                 color=INK)

    fig.suptitle("Tier-4C: direction-locked closure — halo-onset family, "
                 "cap/basis ladder, optimum", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def fig_epsscan(scan, end, fit, path):
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.6), dpi=160)
    fig.patch.set_facecolor("white")
    eps = np.array([row["eps"] for row in scan])
    cp = np.array([row["rep"]["c_paper"] for row in scan])
    c0 = end["c_paper"]

    a1 = axes[0, 0]
    style_ax(a1)
    a1.plot(eps, cp, "-o", color=BLUE, lw=2, ms=6, mec="white",
            label="locked closure")
    a1.axhline(c0, color=BLUE, lw=1.2, ls=":")
    a1.text(eps[-1], c0, "  locked $c_0$ ($t=0$)", color=BLUE, fontsize=8,
            va="bottom", ha="right")
    xe = np.linspace(eps[0] * 0.7, eps[-1] * 1.05, 100)
    a1.plot(xe, c0 * (1.0 - 0.42 * xe ** (2.0 / 3.0)), "--", color=VIOLET,
            lw=1.4, label=r"corpus law $c_0[1-0.42\,\epsilon^{2/3}]$")
    a1.errorbar([0.05], [2.37], yerr=[0.09], fmt="D", color=RED, ms=6,
                capsize=4, mec="white")
    a1.text(0.052, 2.37, r" $\langle r_1\rangle$ $2.37\pm0.09$", color=RED,
            fontsize=8, va="center")
    a1.legend(frameon=False, fontsize=8, loc="lower right")
    a1.set_xlabel(r"$\epsilon$ (achieved sector ratio)", color=INK)
    a1.set_ylabel(r"$c_{\rm paper}$", color=INK)
    a1.set_title("(a) clock charge vs $\\epsilon$, locked closure",
                 fontsize=10, color=INK)

    a2 = axes[0, 1]
    style_ax(a2)
    dfc = 1.0 - cp / c0
    a2.loglog(eps, np.abs(dfc), "o", color=BLUE, ms=7, mec="white")
    xe = np.array([eps[0] * 0.9, eps[-1] * 1.1])
    a2.loglog(xe, fit["amp"] * xe ** fit["p"], "-", color=BLUE, lw=1.6,
              alpha=0.85)
    a2.text(0.04, 0.93, r"locked: $|1-c/c_0| \propto \epsilon^{%.3f}$ (%s)"
            % (fit["p"], "$c<c_0$" if fit["sign"] > 0 else "$c>c_0$"),
            transform=a2.transAxes, color=BLUE, fontsize=9)
    a2.loglog(xe, 0.42 * xe ** (2.0 / 3.0), "--", color=VIOLET, lw=1.2)
    a2.text(xe[-1], 0.42 * xe[-1] ** (2.0 / 3.0),
            r"corpus $0.42\,\epsilon^{2/3}$ ", color=VIOLET, fontsize=8,
            va="bottom", ha="right")
    a2.set_xlabel(r"$\epsilon$", color=INK)
    a2.set_ylabel(r"$|1 - c/c_0|$", color=INK)
    a2.set_title("(b) deficit law: sign and exponent", fontsize=10,
                 color=INK)

    a3x = axes[1, 0]
    style_ax(a3x)
    kp = np.array([row["rep"]["kappa_paper"] for row in scan])
    a3x.plot(eps, kp, "-o", color=BLUE, lw=2, ms=6, mec="white")
    a3x.axhspan(0.802 - 0.018, 0.802 + 0.018, color=RED, alpha=0.10, lw=0)
    a3x.axhline(0.802, color=RED, lw=1.0, ls="--")
    a3x.text(eps[0], 0.803, r" $\langle r_1\rangle$ $0.802\pm0.018$",
             color=RED, fontsize=8)
    for y, lab in ((KLOCK_PAPER, r" $\sqrt3/2$ uniform-halo onset"),
                   (KRING_PAPER, r" $1/\sqrt2$ ring limit")):
        a3x.axhline(y, color=AQUA, lw=1.2)
        a3x.text(eps[0], y, lab, color=AQUA, fontsize=8, va="bottom")
    a3x.set_xlabel(r"$\epsilon$", color=INK)
    a3x.set_ylabel(r"$\kappa_{\rm paper}$", color=INK)
    a3x.set_title("(c) $\\kappa(\\epsilon)$ vs locked thresholds",
                  fontsize=10, color=INK)

    a4 = axes[1, 1]
    style_ax(a4)
    gv = [row["rep"]["g"] for row in scan]
    a4.plot(eps, gv, "-o", color=BLUE, lw=2, ms=6, mec="white")
    a4.axhspan(1.27, 1.35, color=RED, alpha=0.10, lw=0)
    a4.axhline(1.31, color=RED, lw=1.0, ls="--")
    a4.text(eps[0], 1.315, r" paper $g^* = 1.31\pm0.04$", color=RED,
            fontsize=8)
    a4.set_xlabel(r"$\epsilon$", color=INK)
    a4.set_ylabel(r"$g^* = I(p^*)/I(0)$ at $d=1$", color=INK)
    a4.set_title("(d) inertia shape enhancement, locked", fontsize=10,
                 color=INK)

    fig.suptitle("Tier-4C: locked-closure $\\epsilon$-scan (t re-dialled "
                 "per point; zero calibrated constants)", fontsize=11,
                 color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, facecolor="white")
    plt.close(fig)


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t_start = time.time()
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    out = dict(config=dict(N_POL=N_POL, N_FIN=N_FIN, N_COARSE=N_COARSE,
                           N_PROBE=N_PROBE, cap_levels=CAP_LEVELS,
                           budgets=BUD, budgets_scan=BUD_SCAN,
                           quick=QUICK, t_frozen=TFROZEN),
               theory=dict(kappa_lock_ours=KLOCK_OURS,
                           kappa_lock_paper=KLOCK_PAPER,
                           kappa_ring_ours=KRING_OURS,
                           kappa_ring_paper=KRING_PAPER,
                           statement="locked F-halo family: kappa_crit_ours"
                           "(s) = 1/sqrt(8 pi <sin^2 th>_{s^2}); uniform s=1"
                           " gives sqrt(3/16pi) (paper sqrt3/2); equatorial"
                           " concentration degrades it toward 1/sqrt(8pi)"
                           " (paper 1/sqrt2)"))

    # ------------------------------------------------------------------
    print("== GATES (engine revalidation before any new physics) ==")
    gates = {}
    with open(OUTDIR + "/axi_results.json") as fh:
        A2 = json.load(fh)
    with open(OUTDIR + "/radial_results.json") as fh:
        RJ = json.load(fh)
    targets1 = {k: RJ["N_run"][k] for k in ("E2", "E4", "E6", "E0")}

    r1d, f1d, sect1d = ax.radial_profile(TFROZEN)
    fspl = ax.CubicSpline1D(r1d, f1d)
    rsup_e = ax.support_radius(r1d, f1d)
    rsup_c = RSTAR * 1.02

    basis_e = LockedBasis(rsup_e, K=4)
    basis_c = LockedBasis(rsup_c, K=4)
    dom_e0 = rsup_e * 1.5
    dom_c0 = rsup_c * 1.5

    P0h = ax.sector_pieces(qfun_locked(fspl, basis_e, np.zeros(basis_e.n)),
                           N_FIN, dom_e0, dom_e0)
    g1 = {k: P0h[k] / targets1[k] - 1.0 for k in ("E2", "E4", "E6", "E0")}
    Irad = ax.radial_inertia(r1d, f1d)
    g2 = P0h["I"] / Irad - 1.0
    ok1 = (max(abs(v) for v in g1.values()) < 5e-3
           and abs(P0h["BDEG"] - 1.0) < 5e-3)
    print("  G1 p=0 sectors vs step-1 radial (N=%d): %s  degree %.5f  %s"
          % (N_FIN, "  ".join("%s %+.1e" % kv for kv in g1.items()),
             P0h["BDEG"], "PASS" if ok1 else "FAIL"))
    print("  G2 p=0 inertia vs radial: rel %+.2e  %s"
          % (g2, "PASS" if abs(g2) < 5e-3 else "FAIL"))

    P0c = ax.sector_pieces(qfun_locked(ax.f_compacton, basis_c,
                                       np.zeros(basis_c.n)),
                           N_FIN, dom_c0, dom_c0)
    g3 = dict(E6=P0c["E6"] / E6_COMPACTON - 1.0,
              E0=P0c["E0"] / E6_COMPACTON - 1.0,
              I=P0c["I"] / I_COMPACTON - 1.0, degree=P0c["BDEG"])
    ok3 = max(abs(g3[k]) for k in ("E6", "E0", "I")) < 2e-3
    print("  G3 compacton p=0 vs exact: E6 %+.1e  E0 %+.1e  I %+.1e  "
          "degree %.5f  %s" % (g3["E6"], g3["E0"], g3["I"], g3["degree"],
                               "PASS" if ok3 else "FAIL"))

    at = np.array([0.08, -0.05, 0.04, 0.06, -0.04, 0.03,
                   0.05, 0.03, -0.02, -0.04, 0.02, 0.02])
    Pt = ax.sector_pieces(qfun_locked(ax.f_compacton, basis_c, at), N_FIN,
                          dom_c0, dom_c0)
    E0s, Is = a3.sph_check(ax.f_compacton, basis_c, at, ZB, dom_c0)
    g4 = dict(E0=Pt["E0"] / E0s - 1.0, I=Pt["I"] / Is - 1.0,
              degree=Pt["BDEG"])
    ok4 = (abs(g4["E0"]) < 2e-3 and abs(g4["I"]) < 2e-3
           and abs(Pt["BDEG"] - 1.0) < 5e-3)
    print("  G4L locked deformed config (K=4, P4 active) vs independent "
          "spherical quadrature: E0 %+.1e  I %+.1e  degree %.5f  %s"
          % (g4["E0"], g4["I"], g4["degree"], "PASS" if ok4 else "FAIL"))

    Lh0, _, _, _, _ = DLFP(P0h, TFROZEN, 5.0)
    g5 = dict(c_paper=2.0 * Lh0 * CMAP_L,
              vs_solA=2.0 * Lh0 * CMAP_L / A2["solA"]["c_paper"] - 1.0)
    print("  G5 p=0 closure at eps=0.05: c_paper %.5f vs step-2 solA %.5f "
          "(rel %+.1e)  %s" % (g5["c_paper"], A2["solA"]["c_paper"],
                               g5["vs_solA"],
                               "PASS" if abs(g5["vs_solA"]) < 3e-3
                               else "FAIL"))
    gates.update(G1=g1, G1_degree=P0h["BDEG"], G2=g2, G3=g3, G4L=g4, G5=g5)
    out["gates"] = gates
    if not (ok1 and abs(g2) < 5e-3 and ok3 and ok4
            and abs(g5["vs_solA"]) < 3e-3):
        print("!! GATE FAILURE — aborting.")
        with open(OUTDIR + "/axi4_locked_results.json", "w") as fh:
            json.dump(a3.jsonable(out), fh, indent=2)
        sys.exit(1)

    # ------------------------------------------------------------------
    # TASK A: locked halo threshold — derivation table + engine probes
    # ------------------------------------------------------------------
    print("\n== TASK A: LOCKED F-HALO THRESHOLD (derived, then measured) ==")
    print("  derivation: dR = [INT eta^2 s^2 dV] (1/(8pi) "
          "- kappa_ours^2 <sin^2 th>_{s^2})")
    print("  kappa_crit_ours(s) = 1/sqrt(8 pi <sin^2 th>_w);  "
          "kappa_crit_paper(s) = 1/sqrt(2 <sin^2 th>_w)")
    print("  uniform s=1: kappa_ours = sqrt(3/16pi) = %.5f, kappa_paper = "
          "sqrt(3)/2 = %.5f   [task's expected locked threshold: CONFIRMED "
          "for angle-uniform halos]" % (KLOCK_OURS, KLOCK_PAPER))
    print("  FLAG: the lock does NOT remove the angular freedom of F; "
          "equator-concentrated far halos degrade the onset toward the ring "
          "limit 1/sqrt(8pi) = %.5f (paper 1/sqrt2 = %.5f). Measured below."
          % (KRING_OURS, KRING_PAPER))
    derived = []
    for nm, sf in SPROFILES:
        conc = sprofile_moment(sf)
        derived.append(dict(name=nm, conc=conc,
                            kcrit_ours=1.0 / math.sqrt(8.0 * math.pi * conc),
                            kcrit_paper=1.0 / math.sqrt(2.0 * conc)))
    print("   %-18s %12s %14s %14s" % ("angular profile", "<sin^2 th>_w",
                                       "k_crit ours", "k_crit paper"))
    for d in derived:
        print("   %-18s %12.6f %14.5f %14.5f"
              % (d["name"], d["conc"], d["kcrit_ours"], d["kcrit_paper"]))
    out["theory"]["derived_profiles"] = derived

    print("\n  -- engine probes: far shell r_c = 2.1 r_sup, w = 0.45 r_sup, "
          "eta = 0.02 / 0.04 --")
    probes = {"t0": [], "eps": []}
    for tag, fint, rsup, tt in (("t0", ax.f_compacton, rsup_c, 0.0),
                                ("eps", fspl, rsup_e, TFROZEN)):
        cache = {}
        print("   [%s profile]" % ("t = 0 compacton" if tag == "t0"
                                   else "eps = 0.05"))
        for (nm, sf), d in zip(SPROFILES, derived):
            pr = halo_probe(fint, rsup, tt, sf, 2.1, 0.45, N_PROBE,
                            base_cache=cache)
            pr.update(name=nm, conc=d["conc"],
                      kcrit_paper_ideal=d["kcrit_paper"],
                      kcrit_paper=pr["kcrit_paper"][0])
            probes[tag].append(pr)
            print("    %-18s k_crit ours %.5f / %.5f (ideal %.5f)  paper "
                  "%.5f (ideal %.5f)  dI ratio %.3f"
                  % (nm, pr["kcrit_ours"][0], pr["kcrit_ours"][1],
                     d["kcrit_ours"], pr["kcrit_paper"],
                     d["kcrit_paper"], pr["dI_ratio"]))
        # radius scan (uniform + sin^2): approach to the far-field ideal
        rad = []
        for rcf, prof_i in ((1.3, 0), (3.2, 0), (3.2, 2)):
            nm, sf = SPROFILES[prof_i]
            pr = halo_probe(fint, rsup, tt, sf, rcf, 0.45, N_PROBE,
                            base_cache=cache)
            pr.update(name=nm, rc_fac=rcf, conc=derived[prof_i]["conc"],
                      kcrit_paper_ideal=derived[prof_i]["kcrit_paper"])
            rad.append(pr)
            print("    radius scan %-12s r_c=%.1f r_sup: k_ours %.5f "
                  "(ideal %.5f)" % (nm, rcf, pr["kcrit_ours"][0],
                                    derived[prof_i]["kcrit_ours"]))
        probes[tag + "_radius"] = rad
        # mid-field / overlap probe (Task D): shell straddling the edge
        rc_mid = 1.0 / 1.02 if tag == "t0" else 0.95
        prm = halo_probe(fint, rsup, tt, SPROFILES[0][1], rc_mid, 0.25,
                         N_PROBE, base_cache=cache)
        prm.update(name="edge shell (overlap)")
        probes[tag + "_edge"] = prm
        print("    edge shell r_c=%.2f r_sup: dI(0.02)=%.4f dI(0.04)=%.4f "
              "ratio %.3f (LINEAR -> profile-relaxation direction, not an "
              "onset channel); dE ratio %.3f (quadratic)"
              % (rc_mid, prm["dI"][0], prm["dI"][1], prm["dI_ratio"],
                 prm["dE_ratio"]))
    out["probes"] = probes

    print("\n  -- placement of the corpus tuples (measure, don't assume) --")
    plc = [placement(k, nm) for nm, k in CORPUS_TUPLES]
    print("   %-20s %8s %9s %11s %11s %10s %10s %8s"
          % ("tuple", "k_paper", "k_ours", "x uniform", "x ring",
             "st.unif", "st.ring", "req<s2>"))
    for p in plc:
        print("   %-20s %8.4f %9.5f %11.3f %11.3f %10s %10s %8.4f"
              % (p["name"], p["kappa_paper"], p["kappa_ours"], p["x_lock"],
                 p["x_ring"], "STABLE" if p["stable_uniform"] else "unstab",
                 "STABLE" if p["stable_ring"] else "unstab", p["req_conc"]))
    print("   (req<s2> = minimal halo concentration <sin^2 th>_w that opens "
          "descent; uniform = 2/3, ring = 1; deep-BPS needs %.6f = 6/7 "
          "exactly — it IS the s=sin^2 threshold sqrt(7/12))"
          % plc[2]["req_conc"])
    out["placement"] = plc

    # ------------------------------------------------------------------
    # TASK B: locked closure at eps = 0.05, cap ladder, K = 4 then K = 5
    # ------------------------------------------------------------------
    print("\n== TASK B: LOCKED CLOSURE at eps = 0.05 (t frozen), cap "
          "ladder ==")
    cache_e = {}
    seeds0 = [np.zeros(basis_e.n),
              spheroid_seed(basis_e, 0.7),
              spheroid_seed(basis_e, 0.45),
              spheroid_seed(basis_e, 0.45, scale=1.0),
              skirt_seed(basis_e)]
    ladder = []
    x_e, rep_e, lev_used = None, None, None
    for il, lev in enumerate(CAP_LEVELS):
        caps = dict(abox=lev["abox"], amax=lev["amax"],
                    stretch_lim=lev["stretch_lim"])
        dom = rsup_e * lev["domfac"]
        seeds = seeds0 if x_e is None else [x_e] + seeds0[:2]
        rep, x = run_locked("R-eps %s" % lev["name"], fspl, basis_e, TFROZEN,
                            dom, caps, lev["N_opt"], seeds,
                            BUD if il == 0 else BUD_SCAN, cache_e,
                            max_outer=MAX_OUTER)
        rep["level"] = lev["name"]
        ladder.append(dict(label="K4/" + lev["name"],
                           c_paper=rep["c_paper"],
                           kappa_paper=rep["kappa_paper"], g=rep["g"],
                           V=rep["V"], interior=rep["interior"],
                           max_cap_frac=rep["max_cap_frac"]))
        rep_e, x_e, lev_used = rep, x, il
        if rep["interior"]:
            break
    out["ladder_eps005"] = ladder
    lev = CAP_LEVELS[lev_used]
    caps_f = dict(abox=lev["abox"], amax=lev["amax"],
                  stretch_lim=lev["stretch_lim"])
    domfac_f, Nopt_f = lev["domfac"], lev["N_opt"]
    print("   accepted cap level: %s (domfac %.2f, amax %.2f, abox %.2f)  "
          "interior=%s" % (lev["name"], lev["domfac"], lev["amax"],
                           lev["abox"], rep_e["interior"]))

    # basis-size ladder: K = 5 (15 params) warm-started at the same caps
    basis_e5 = LockedBasis(rsup_e, K=5)
    x5_seed = project_field(basis_e, x_e, basis_e5)
    cache_e5 = {}
    rep_e5, x_e5 = run_locked("R-eps K5", fspl, basis_e5, TFROZEN,
                              rsup_e * domfac_f, caps_f, Nopt_f,
                              [x5_seed, np.zeros(basis_e5.n)], BUD_SCAN,
                              cache_e5, max_outer=MAX_OUTER)
    ladder.append(dict(label="K5/" + lev["name"], c_paper=rep_e5["c_paper"],
                       kappa_paper=rep_e5["kappa_paper"], g=rep_e5["g"],
                       V=rep_e5["V"], interior=rep_e5["interior"],
                       max_cap_frac=rep_e5["max_cap_frac"]))
    basis_shift = dict(
        c=abs(rep_e5["c_paper"] / rep_e["c_paper"] - 1.0),
        kappa=abs(rep_e5["kappa_paper"] / rep_e["kappa_paper"] - 1.0),
        g=abs(rep_e5["g"] / rep_e["g"] - 1.0),
        V=abs(rep_e5["V"] / rep_e["V"] - 1.0))
    print("   basis ladder K=4 -> K=5 relative shifts: c %.2e  kappa %.2e  "
          "g %.2e  V %.2e" % (basis_shift["c"], basis_shift["kappa"],
                              basis_shift["g"], basis_shift["V"]))
    out["locked_eps005"] = dict(K4=rep_e, K5=rep_e5,
                                basis_shift=basis_shift)

    # stability placement of the locked solution itself
    for nm, rp in (("K4", rep_e), ("K5", rep_e5)):
        ko = rp["kappa_ours"]
        print("   locked solution %s: kappa_ours = %.5f -> vs uniform-halo "
              "threshold %.5f: %s | vs ring/tilt threshold %.5f: %s | "
              "req. destabilising concentration <sin^2 th>_w = %.4f"
              % (nm, ko, KLOCK_OURS,
                 "BELOW (stable)" if ko < KLOCK_OURS else "ABOVE (unstable)",
                 KRING_OURS,
                 "below" if ko < KRING_OURS else "ABOVE",
                 1.0 / (8.0 * math.pi * ko * ko)))

    # saturated references (where the residual ring channel would run)
    print("\n   -- saturated references --")
    repG_e, _ = run_locked_G("G-eps", fspl, basis_e, TFROZEN,
                             rsup_e * domfac_f, caps_f, Nopt_f,
                             [x_e, np.zeros(basis_e.n),
                              spheroid_seed(basis_e, 0.7)], BUD_SCAN,
                             cache_e)
    uref_e = uniform_sat_reference(P0h, TFROZEN)
    print("   uniform-channel hedgehog-core reference (eps=0.05): "
          "G_L*=%.6f  c_paper=%.5f at kappa_paper=%.5f"
          % (uref_e["G"], uref_e["c_paper"], uref_e["kappa_paper"]))
    out["saturated"] = dict(ring_eps=repG_e, uniform_ref_eps=uref_e)

    # ------------------------------------------------------------------
    # TASK C: eps-scan + t = 0 endpoint (same caps, same discipline)
    # ------------------------------------------------------------------
    print("\n== TASK C: LOCKED eps-SCAN (t re-dialled per point) ==")
    scan = [dict(eps=float(TFROZEN * (sect1d[0] + sect1d[1])
                           / (sect1d[2] + sect1d[3])),
                 t=TFROZEN, rep=rep_e)]
    x_by_eps = {0.05: x_e}
    for eps in [e for e in EPS_LIST if abs(e - 0.05) > 1e-12]:
        te, re_, fe, secte, ratio = ax.dial_eps(eps, TFROZEN * eps / 0.05)
        fspl_x = ax.CubicSpline1D(re_, fe)
        rsup_x = ax.support_radius(re_, fe)
        basis_x = LockedBasis(rsup_x, K=4)
        src = min(x_by_eps.keys(), key=lambda k: abs(k - eps))
        cache_x = {}
        rep, xx = run_locked("R-e%.3f" % eps, fspl_x, basis_x, te,
                             rsup_x * domfac_f, caps_f, Nopt_f,
                             [x_by_eps[src]], BUD_SCAN, cache_x,
                             max_outer=MAX_OUTER_SCAN)
        x_by_eps[eps] = xx
        scan.append(dict(eps=float(ratio), t=float(te), rep=rep))
    scan.sort(key=lambda s: s["eps"])

    print("\n   t = 0 endpoint (compacton base), same caps:")
    cache_c = {}
    src = min(x_by_eps.keys())
    end_rep, x_c = run_locked("R-t0", ax.f_compacton, basis_c, 0.0,
                              rsup_c * domfac_f, caps_f, Nopt_f,
                              [x_by_eps[src], np.zeros(basis_c.n),
                               spheroid_seed(basis_c, 0.7)],
                              BUD_SCAN, cache_c, max_outer=MAX_OUTER)
    endpoint_escalated = False
    while (not end_rep["interior"]) and lev_used + 1 < len(CAP_LEVELS):
        endpoint_escalated = True
        lev_used += 1
        lev = CAP_LEVELS[lev_used]
        caps_f = dict(abox=lev["abox"], amax=lev["amax"],
                      stretch_lim=lev["stretch_lim"])
        domfac_f, Nopt_f = lev["domfac"], lev["N_opt"]
        print("   t=0 endpoint cap-limited -> escalating ALL closures to "
              "%s and re-running scan chain" % lev["name"])
        end_rep, x_c = run_locked("R-t0 %s" % lev["name"], ax.f_compacton,
                                  basis_c, 0.0, rsup_c * domfac_f, caps_f,
                                  Nopt_f, [x_c], BUD_SCAN, cache_c,
                                  max_outer=MAX_OUTER)
    if endpoint_escalated:
        # re-run every scan point at the final caps (consistent family)
        scan2 = []
        for row in scan:
            te = row["t"]
            if abs(te - TFROZEN) < 1e-15:
                fspl_x, rsup_x = fspl, rsup_e
            else:
                _, re_, fe, _, _ = ax.dial_eps(row["eps"], te)
                fspl_x = ax.CubicSpline1D(re_, fe)
                rsup_x = ax.support_radius(re_, fe)
            basis_x = LockedBasis(rsup_x, K=4)
            rep, xx = run_locked("R2-e%.3f" % row["eps"], fspl_x, basis_x,
                                 te, rsup_x * domfac_f, caps_f, Nopt_f,
                                 [np.asarray(row["rep"]["a"])], BUD_SCAN,
                                 {}, max_outer=MAX_OUTER_SCAN)
            scan2.append(dict(eps=row["eps"], t=te, rep=rep))
        scan = scan2
        rep_e = [s["rep"] for s in scan
                 if abs(s["t"] - TFROZEN) < 1e-15][0]
        ladder.append(dict(label="K4/" + lev["name"],
                           c_paper=rep_e["c_paper"],
                           kappa_paper=rep_e["kappa_paper"], g=rep_e["g"],
                           V=rep_e["V"], interior=rep_e["interior"],
                           max_cap_frac=rep_e["max_cap_frac"]))
        out["locked_eps005"]["K4"] = rep_e

    # K = 5 at t = 0 (endpoint basis systematic)
    basis_c5 = LockedBasis(rsup_c, K=5)
    end_rep5, _ = run_locked("R-t0 K5", ax.f_compacton, basis_c5, 0.0,
                             rsup_c * domfac_f, caps_f, Nopt_f,
                             [project_field(basis_c, x_c, basis_c5)],
                             BUD_SCAN, {}, max_outer=MAX_OUTER_SCAN)
    repG_c, _ = run_locked_G("G-t0", ax.f_compacton, basis_c, 0.0,
                             rsup_c * domfac_f, caps_f, Nopt_f,
                             [x_c, np.zeros(basis_c.n)], BUD_SCAN, cache_c)
    uref_c = uniform_sat_reference(P0c, 0.0)
    print("   uniform-channel hedgehog-core reference (t=0): G_L*=%.6f  "
          "c_paper=%.5f at kappa_paper=%.5f"
          % (uref_c["G"], uref_c["c_paper"], uref_c["kappa_paper"]))
    out["saturated"]["ring_t0"] = repG_c
    out["saturated"]["uniform_ref_t0"] = uref_c
    out["endpoint_t0"] = dict(K4=end_rep, K5=end_rep5,
                              escalated=endpoint_escalated)
    out["scan"] = scan

    print("\n   scan summary (locked, K=4, caps %s):" %
          CAP_LEVELS[lev_used]["name"])
    print("   %-8s %10s %12s %8s %8s %8s %10s %9s"
          % ("eps", "c_paper", "kappa_paper", "g*", "lam_eff", "V*",
             "maxcap%", "interior"))
    for row in scan:
        rp = row["rep"]
        print("   %-8.5f %10.5f %12.5f %8.4f %8.4f %8.4f %9.0f%% %9s"
              % (row["eps"], rp["c_paper"], rp["kappa_paper"], rp["g"],
                 rp["lam_eff"], rp["V"], 100 * rp["max_cap_frac"],
                 "yes" if rp["interior"] else "NO"))
    print("   t=0 endpoint: c_paper %.5f (K5: %.5f)  kappa %.5f  g %.4f  "
          "V %.4f  maxcap %.0f%%  interior=%s"
          % (end_rep["c_paper"], end_rep5["c_paper"],
             end_rep["kappa_paper"], end_rep["g"], end_rep["V"],
             100 * end_rep["max_cap_frac"], end_rep["interior"]))

    # deficit fit vs the locked t = 0 endpoint
    e_arr = np.array([row["eps"] for row in scan])
    c_arr = np.array([row["rep"]["c_ours"] for row in scan])
    dfc = 1.0 - c_arr / end_rep["c_ours"]
    sgn = float(np.sign(dfc.mean()))
    pfit, lnA = np.polyfit(np.log(e_arr), np.log(np.abs(dfc)), 1)
    fit = dict(p=float(pfit), amp=float(math.exp(lnA)), sign=sgn,
               c0_ours=end_rep["c_ours"], c0_paper=end_rep["c_paper"],
               deficits=[float(v) for v in dfc])
    out["fit"] = fit
    print("\n   deficit fit vs locked t=0 endpoint:  1 - c/c0 = %+.4f "
          "eps^%.4f   (sign: %s;  corpus law: +0.42 eps^0.667, c "
          "DECREASING)" % (sgn * math.exp(lnA), pfit,
                           "positive, c < c0 (corpus sign)" if sgn > 0
                           else "NEGATIVE, c > c0 (anti-corpus)"))

    # ------------------------------------------------------------------
    # TASK D: over-spin onset table in paper units
    # ------------------------------------------------------------------
    print("\n== TASK D: LOCKED OVER-SPIN ONSETS (paper units) vs corpus "
          "1.000 +/- 0.004 and 0.866 ==")
    onset = []
    for d, pt0, pe in zip(derived, probes["t0"], probes["eps"]):
        onset.append(dict(channel="F-halo " + d["name"], derived=
                          d["kcrit_paper"], measured_t0=pt0["kcrit_paper"],
                          measured_eps=pe["kcrit_paper"]))
    onset.append(dict(channel="F-halo ring limit (s^2 -> delta)",
                      derived=KRING_PAPER, measured_t0=None,
                      measured_eps=None))
    onset.append(dict(channel="tilt halo (EXCLUDED by the lock)",
                      derived=KRING_PAPER, measured_t0=0.21625 * KMAP,
                      measured_eps=0.21688 * KMAP,
                      note="axi3 measured values, for reference"))
    print("   %-38s %10s %13s %13s" % ("channel", "derived",
                                       "measured t=0", "measured eps"))
    for o in onset:
        print("   %-38s %10.5f %13s %13s"
              % (o["channel"], o["derived"],
                 "%.5f" % o["measured_t0"] if o["measured_t0"] else "—",
                 "%.5f" % o["measured_eps"] if o["measured_eps"] else "—"))
    print("   corpus over-spin onset quote:  1.000 +/- 0.004")
    print("   locked uniform-halo onset:     %.5f  (= sqrt3/2; the task's "
          "0.866)" % KLOCK_PAPER)
    print("   locked s^2~|cos| onset:        1.00000 EXACT in the ideal "
          "far-field limit (measured smoothed variant above)")
    print("   edge shells: dI linear in eta (ratios %.3f / %.3f at t=0 / "
          "eps) -> first-order profile relaxation, no onset threshold "
          "at the compacton edge" % (probes["t0_edge"]["dI_ratio"],
                                     probes["eps_edge"]["dI_ratio"]))
    out["onset_table"] = onset

    # ------------------------------------------------------------------
    # final tuple vs benchmark
    # ------------------------------------------------------------------
    print("\n== LOCKED TUPLE vs <r1> BENCHMARK (sigma pulls) ==")
    tt = TARGETS
    pulls = {}
    for nm, rp in (("K4", rep_e), ("K5", rep_e5)):
        pl = dict(c=(rp["c_paper"] - tt["c"][0]) / tt["c"][1],
                  kappa=(rp["kappa_paper"] - tt["kappa"][0]) / tt["kappa"][1],
                  g=(rp["g"] - tt["g"][0]) / tt["g"][1],
                  lam=(rp["lam_eff"] - tt["lam"][0]) / tt["lam"][1],
                  V=(rp["V"] - tt["V"][0]) / tt["V"][1])
        pulls[nm] = pl
        print("   %s: c=%.4f (%+.1f sig)  kappa=%.4f (%+.1f sig)  g*=%.3f "
              "(%+.1f sig)  lam_eff=%.3f (%+.1f sig)  V=%.3f (%+.1f sig)  "
              "Erot/E=%.6f (target 0.2500)"
              % (nm, rp["c_paper"], pl["c"], rp["kappa_paper"], pl["kappa"],
                 rp["g"], pl["g"], rp["lam_eff"], pl["lam"], rp["V"],
                 pl["V"], rp["erot_frac"]))
    print("   (lam_eff = sqrt(2<z^2>/<rho^2>), (1-q0)-weighted; = 1 for a "
          "sphere — NOT the paper's ansatz lambda, reported for shape "
          "reference only)")
    out["pulls"] = pulls
    out["targets"] = TARGETS

    # ------------------------------------------------------------------
    out["runtime_s"] = time.time() - t_start
    with open(OUTDIR + "/axi4_locked_results.json", "w") as fh:
        json.dump(a3.jsonable(out), fh, indent=2)
    fig_landscape(fspl, basis_e, rep_e, ladder, rep_e5, probes["t0"],
                  probes["eps"], derived,
                  OUTDIR + "/axi4_locked_landscape.png")
    fig_epsscan(scan, end_rep, fit, OUTDIR + "/axi4_locked_epsscan.png")
    print("\nruntime: %.1f s" % out["runtime_s"])
    print("wrote axi4_locked_results.json, axi4_locked_landscape.png, "
          "axi4_locked_epsscan.png")


if __name__ == "__main__":
    main()
