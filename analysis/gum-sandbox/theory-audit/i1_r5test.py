#!/usr/bin/env python3
"""
i1_r5test.py — Phase I1 (ROADMAP_v6_SYNTHESIS): the <r5> corroboration test.

QUESTION. The corpus's frozen quartic-regime closure <r5> reports
c_q = 3.1 +/- 0.2 (used at 6 sigma to kill the shallow-knot neutrino,
paper Sec. VII.G; I.8 manifest: "r5 (quartic-regime closure
c_q = 3.1 +/- 0.2; scaling audit)").  The campaign's repaired (saturated)
closure gives exactly 64 sqrt2/(9 pi) = 3.2011247 at eps -> 0, and the
F-R5 halo threshold kappa_ours = 1/sqrt(8 pi) is t-INDEPENDENT, so the
saturated branch exists at every eps.  Was <r5> measuring the saturated
branch?

WHAT THIS SCRIPT COMPUTES (all in the campaign's frozen conventions,
a6 = a0 = m = 1, a2 = a4 = t, compacton-anchored unit map, zero calibrated
constants; every engine REUSED unmodified from tier2-closure):

  1. The eps-dial ladder: t(eps) for eps in {0.02, 0.05, 0.1, 0.2, 0.5,
     1.0, 1.5, 2.0}, eps = t(E2+E4)/(E6+E0) on the step-1 hedgehog radial
     solution (the campaign's dial, radial_solve.py).  Note at eps >= 1 the
     kinetic FRACTION of the static energy is eps/(1+eps) -> 0.5+.
  2. Reading (a), eps-dial continuation:
     (i)  the saturated curve c_sat(eps) = (4/pi) G_t*, where
          G_t = t(E2+E4) + E6 + E0 - I/(16 pi) is minimised over the core
          (hedgehog base profile at that t, + the axi3 12-mode
          contour/tilt basis + analytic dilation).  Derivation of the
          t > 0 saturated objective: axi3_solve.py Sec. "saturated
          closure" — the halo ledger is t-blind because the soft halo's
          E2/E4 cost vanishes in the delocalised limit while its E0 and I
          contributions keep the exact ratio 1/(16 pi); the factorisation
          (core min G; kappa_paper = 1/sqrt2; c_paper = (4/pi) G*) holds
          verbatim with G carrying the t-sectors of the CORE.
     (ii) the direct (restricted, saddle-branch) closures on the corpus's
          method-class families: family A (hedgehog + dilation, the
          scaling rung's continuation) and family B (spheroidal shell,
          the paper's G.5 ansatz) — interior fixed points of the
          fixed-L Routhian + clock, i.e. what a restricted Newton-Krylov
          would find.
  3. Reading (b), pure quartic limit a6 = 0, E = t(E2+E4) + E0 with
     t = 1 (with a6 = 0, t is pure scale: E and G scale as t^{3/2}, so
     the model has a canonical normalisation only at a2 = a4 = 1):
     hedgehog radial solve, direct closure, saturated closure over the
     same 12-mode basis, plus the analytic Bogomolny-type floors:
       - Faddeev-type:      E2 + E4 >= 3 pi |B|      (these conventions)
       - quartic+potential (classic BPS-submodel bound):
            E4 + E0 >= 16*2^{1/4} G(7/4)G(3/2)/G(13/4) = 6.07912
       - saturated quartic floor:
            G_q >= (2^{-1/4}/pi) INT_{S^3} [(1-cos psi)^2
                    + sin^2 psi cos^2 Theta]^{1/4} dOmega   (quadrature)
  4. Analytic floors for the eps-dial saturated curve (rigorous):
       G_t >= 16 sqrt2/9  (drop t(E2+E4) >= 0)  =>  c_sat(eps) >= 3.20112
       G_t >= 16 sqrt2/9 + 3 pi t  =>  c_sat(eps) >= 3.20112 + 12 t(eps).
  5. Convergence ladders (radial N and rmax; engine N; mode-count) for
     every quoted number; deterministic throughout (no RNG).

Outputs: i1_results.json (+ printed tables -> i1_RESULTS.md written from
them).  Deterministic; runtime ~10-20 min.
"""

import json
import math
import sys
import time

import numpy as np

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit"
T2DIR = "/home/user/fork_frankensim/analysis/gum-sandbox/tier2-closure"
sys.path.insert(0, T2DIR)

import radial_solve as rs          # step-1 machinery (REUSED)
import axi_solve as ax             # step-2 engine + families (REUSED)
import axi3_solve as a3            # step-3 basis/objective/NM (REUSED)
import gstar_solve as gs           # G3 exact minimiser (REUSED)

from numpy.polynomial.legendre import leggauss

SQ2 = math.sqrt(2.0)
PI = math.pi
SIXTEEN_PI = 16.0 * PI
GSTAR = 16.0 * SQ2 / 9.0                    # exact t->0 saturated constant
CSAT0 = 64.0 * SQ2 / (9.0 * PI)             # = (4/pi) GSTAR = 3.2011247
R5_C, R5_SIG = 3.1, 0.2                     # the corpus's <r5> release
C0_CORPUS = 128.0 * math.sqrt(42.0) / (105.0 * PI)   # corpus deep-BPS 2.5148
RUNG = 256.0 / (15.0 * math.sqrt(7.0) * PI)          # scaling rung 2.0533
CMAP_L, KMAP = ax.CMAP_L, ax.KMAP
EPS_LIST = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
QUICK = "--quick" in sys.argv

# NM budgets (reduced vs axi3 full: convergence quantified by ladders below)
if QUICK:
    N_OPT, N_POL, N_FIN, N_COARSE = 112, 176, 256, 176
    B_TRI, B_COLD, B_POL = 120, 500, 100
else:
    N_OPT, N_POL, N_FIN, N_COARSE = 192, 320, 480, 320
    B_TRI, B_COLD, B_POL = 220, 1500, 260
N_FAM, N_FAM_C = 320, 224                   # direct-family engine grids


def pull(x):
    return (x - R5_C) / R5_SIG


# ----------------------------------------------------------------------
# 1. the eps-dial (campaign convention: ratio on the step-1 radial solve)
# ----------------------------------------------------------------------
def rmax_for(eps):
    if eps <= 0.11:
        return 6.0
    if eps <= 0.21:
        return 8.0
    if eps <= 0.51:
        return 10.0
    return 12.0


def dial(eps_target, t_guess, N=4000, rmax=6.0, rtol=3.0e-4, iters=14):
    """Fixed point on t so that t(E2+E4)/(E6+E0) = eps_target on the solved
    hedgehog profile (rs machinery reused; same convention as axi_solve)."""
    r = rs.make_grid(N, rmax)
    f = rs.init_profile(r)
    t = t_guess
    ratio, sect, gmax = None, None, None
    for _ in range(iters):
        f, E, sect, g, gmax, dE, it = rs.solve_profile(r, f, t)
        E2, E4, E6, E0 = sect
        ratio = t * (E2 + E4) / (E6 + E0)
        if abs(ratio / eps_target - 1.0) < rtol:
            break
        t *= eps_target / ratio
    return dict(t=t, ratio=ratio, sect=list(sect), r=r, f=f, gmax=gmax,
                N=N, rmax=rmax)


def dial_ladder(eps, sol):
    """Convergence of the dial: re-solve at same t with N x2 and rmax +2."""
    out = {}
    for tag, N, rmax in (("N2x", 2 * sol["N"], sol["rmax"]),
                         ("rmax+", sol["N"], sol["rmax"] + 2.0)):
        r = rs.make_grid(N, rmax)
        f, E, sect, g, gmax, dE, it = rs.solve_profile(
            r, rs.init_profile(r), sol["t"])
        E2, E4, E6, E0 = sect
        ratio = sol["t"] * (E2 + E4) / (E6 + E0)
        out[tag] = dict(ratio=ratio, d_ratio_rel=ratio / sol["ratio"] - 1.0,
                        dE0_rel=E0 / sol["sect"][3] - 1.0)
    return out


# ----------------------------------------------------------------------
# 2a. direct (restricted) closures — the saddle branch, method-class of
#     the corpus (interior optimum of a restricted family)
# ----------------------------------------------------------------------
def closure_direct(F, t, Lhi):
    """ax.closure_report with a wide clock bracket (large-eps c is big)."""
    L, lam, d, R, iI, resid = ax.clock_solve(F, t, Llo=0.5, Lhi=Lhi)
    _, (e2, e4, e6, e0, iI2, Estat) = ax.R_of(F, t, lam, d, L)
    d0, _ = ax.min_over_d(F, t, 1.0, 0.0)
    Erot = L * L / (2.0 * iI)
    return dict(t=t, L=L, c_ours=2.0 * L, kappa_ours=L / iI,
                c_paper=2.0 * L * CMAP_L, kappa_paper=(L / iI) * KMAP,
                lam=lam, d=d, d0=d0, V=(d / d0) ** 3,
                g=F["iI"](lam) / F["iI"](1.0), R=R, Estat=Estat,
                Erot=Erot, erot_frac=Erot / R,
                eps_at_solution=t * (e2 + e4) / (e6 + e0),
                clock_resid=resid, E2=e2, E4=e4, E6=e6, E0=e0, I=iI,
                above_halo_threshold=bool((L / iI) * KMAP > 1.0 / SQ2))


def direct_pair(fspl, t, dom_tail, N, lam_grid, Lhi):
    """Family A (hedgehog+dilation) and family B (spheroidal shell)."""
    P1 = ax.sector_pieces(ax.make_famA(fspl, 1.0, 1.0), N, dom_tail,
                          dom_tail)
    FA = ax.famA_funcs(P1)
    solA = closure_direct(FA, t, Lhi)
    solA["degree"] = P1["BDEG"]
    FB, _ = ax.famB_funcs(fspl, dom_tail, N, lam_grid)
    solB = closure_direct(FB, t, Lhi)
    return solA, solB, P1


# ----------------------------------------------------------------------
# 2b. saturated closure at t > 0 — G_t minimised over hedgehog base +
#     axi3 12-mode contour/tilt basis + analytic dilation (axi3 REUSED)
# ----------------------------------------------------------------------
def saturated_at(fspl, t, r_sup_b, dom, tag, log=True):
    basis = a3.GenBasis(r_sup_b, K=3, a_l=(0, 2), b_l=(1, 2))
    active = list(range(basis.n))
    stretch_lim = 0.50
    P0_cache = {}

    def P0(N):
        if N not in P0_cache:
            P0_cache[N] = ax.sector_pieces(
                a3.make_qfun(fspl, basis, np.zeros(basis.na),
                             np.zeros(basis.nb)), N, dom, dom)
        return P0_cache[N]

    # hedgehog-core reference (0 modes)
    obj0 = a3.Objective(fspl, basis, t, N_OPT, dom, [], stretch_lim,
                        mode="G")
    rep0 = a3.finalize_G(obj0, np.zeros(0), N_FIN, P0(N_FIN))
    rep0["coarse"] = a3.finalize_G(obj0, np.zeros(0), N_COARSE,
                                   P0(N_COARSE))

    # 12-mode optimisation (deterministic NM, fixed seeds)
    obj = a3.Objective(fspl, basis, t, N_OPT, dom, active, stretch_lim,
                       mode="G")
    steps = np.array([0.06 if i < basis.na else 0.15 for i in active])
    seeds = [np.zeros(basis.n),
             a3.seed_vec(basis, active, a_map={0: 0.0, 2: 0.12}),
             a3.seed_vec(basis, active, a_map={0: 0.0, 2: -0.12},
                         b_map={1: 0.35, 2: 0.15}),
             a3.seed_vec(basis, active, b_map={1: 0.45, 2: 0.20})]
    x, fv = a3.optimize(obj, seeds, steps, B_TRI, B_COLD)
    objP = a3.Objective(fspl, basis, t, N_POL, dom, active, stretch_lim,
                        mode="G")
    x, fv, _ = a3.nelder_mead(objP, x, steps * 0.08, maxfev=B_POL)
    rep = a3.finalize_G(obj, x, N_FIN, P0(N_FIN))
    rep["coarse"] = a3.finalize_G(obj, x, N_COARSE, P0(N_COARSE))
    rep["conv_c_rel"] = abs(rep["coarse"]["c_paper"] / rep["c_paper"] - 1.0)
    rep["nfev"] = obj.nfev + objP.nfev
    # the core's own sector ratio and kinetic fraction of E_stat_tot
    P = a3.make_qfun(fspl, basis, np.asarray(rep["a"]), np.asarray(rep["b"]))
    PS = ax.sector_pieces(P, N_FIN, dom, dom)
    d = rep["d"]
    e24 = t * (PS["E2"] * d + PS["E4"] / d)
    e60 = PS["E6"] / d ** 3 + PS["E0"] * d ** 3
    rep["eps_core"] = e24 / e60
    rep["kin_frac_tot"] = e24 / rep["Estat_tot"]
    if log:
        print("   %-8s G*(12-mode)=%.6f -> c_sat=%.5f  [hedgehog-core "
              "%.5f]  d=%.3f  B=%.5f  conv %.1e  eps_core=%.4f  nfev=%d"
              % (tag, rep["G"], rep["c_paper"], rep0["c_paper"], rep["d"],
                 rep["degree"], rep["conv_c_rel"], rep["eps_core"],
                 rep["nfev"]))
        fl = rep["flags"]
        if any(fl.values()):
            print("   %-8s boundary flags: %s"
                  % (tag, [k for k, v in fl.items() if v]))
    return rep, rep0


# ----------------------------------------------------------------------
# 3. reading (b): pure quartic model a6 = 0, t = 1
# ----------------------------------------------------------------------
def assemble_gen(r, f, t2, t4, c6, c0):
    """rs.assemble with general sector coefficients (E = t2 E2 + t4 E4
    + c6 E6 + c0 E0); midpoint rule, exact gradient + tridiag Hessian.
    Verbatim generalisation of radial_solve.assemble (gate GQ1 checks it
    against the original at t2=t4=t, c6=c0=1)."""
    h = np.diff(r)
    rm = 0.5 * (r[1:] + r[:-1])
    rm2 = rm * rm
    fm = 0.5 * (f[1:] + f[:-1])
    d = np.diff(f) / h
    d2 = d * d
    s = np.sin(fm)
    c = np.cos(fm)
    s2, s3, s4 = s * s, s * s * s, s * s * s * s
    sin2f = 2.0 * s * c
    cos2f = c * c - s2
    L2 = d2 * rm2 + 2.0 * s2
    L4 = s2 * (2.0 * d2 + s2 / rm2)
    L6 = s4 * d2 / rm2
    L0 = (1.0 - np.cos(fm)) * rm2
    E2 = float(np.sum(h * L2))
    E4 = float(np.sum(h * L4))
    E6 = float(np.sum(h * L6))
    E0 = float(np.sum(h * L0))
    E = t2 * E2 + t4 * E4 + c6 * E6 + c0 * E0
    L_d = 2.0 * t2 * d * rm2 + 4.0 * t4 * s2 * d + 2.0 * c6 * s4 * d / rm2
    L_f = (t2 * 2.0 * sin2f + t4 * (2.0 * sin2f * d2 + 4.0 * s3 * c / rm2)
           + c6 * 4.0 * d2 * s3 * c / rm2 + c0 * s * rm2)
    gL = 0.5 * h * L_f - L_d
    gR = 0.5 * h * L_f + L_d
    g = gR[:-1] + gL[1:]
    L_dd = 2.0 * t2 * rm2 + 4.0 * t4 * s2 + 2.0 * c6 * s4 / rm2
    L_fd = 2.0 * t4 * d * sin2f * 2.0 + 8.0 * c6 * s3 * c * d / rm2
    L_ff = (t2 * 4.0 * cos2f + t4 * (4.0 * cos2f * d2
                                     + 4.0 * (3.0 * s2 * c * c - s4) / rm2)
            + c6 * 4.0 * d2 * (3.0 * s2 * c * c - s4) / rm2
            + c0 * c * rm2)
    HLL = 0.25 * h * L_ff - L_fd + L_dd / h
    HRR = 0.25 * h * L_ff + L_fd + L_dd / h
    HLR = 0.25 * h * L_ff - L_dd / h
    diag = HRR[:-1] + HLL[1:]
    off = HLR[1:-1]
    return E, (E2, E4, E6, E0), g, diag, off


def solve_profile_gen(r, f0, t2, t4, c6, c0, gtol=1e-12):
    """flow + damped Newton (rs pattern) on assemble_gen."""
    f = f0.copy()
    E, _, g, diagd, _ = assemble_gen(r, f, t2, t4, c6, c0)
    dt = 0.5
    for _ in range(8000):                        # preconditioned flow
        if np.max(np.abs(g)) < 1e-2:
            break
        prec = np.maximum(diagd, 1e-3 * np.max(diagd))
        fn = f.copy()
        fn[1:-1] = f[1:-1] - dt * g / prec
        np.clip(fn, 0.0, np.pi, out=fn)
        En, _, gn, diagn, _ = assemble_gen(r, fn, t2, t4, c6, c0)
        if En < E:
            f, E, g, diagd = fn, En, gn, diagn
            dt = min(dt * 1.1, 5.0)
        else:
            dt *= 0.5
            if dt < 1e-9:
                break
    E, sect, g, diagd, off = assemble_gen(r, f, t2, t4, c6, c0)
    lam = 1e-3
    for _ in range(300):                         # Newton
        gmax = float(np.max(np.abs(g)))
        if gmax < gtol:
            break
        scale = np.abs(diagd) + 1e-30
        accepted = False
        while lam < 1e14:
            step = rs.thomas(diagd + lam * scale, off, -g)
            if np.all(np.isfinite(step)):
                fn = f.copy()
                fn[1:-1] += step
                En, sectn, gn, diagn, offn = assemble_gen(r, fn, t2, t4,
                                                          c6, c0)
                if np.isfinite(En) and float(np.max(np.abs(gn))) < gmax:
                    f, E, sect, g, diagd, off = fn, En, sectn, gn, diagn, offn
                    lam = max(lam * 0.25, 1e-14)
                    accepted = True
                    break
            lam *= 10.0
        if not accepted:
            break
    return f, E, sect, float(np.max(np.abs(g)))


class QuarticObjective(a3.Objective):
    """axi3 Objective with the sextic sector removed (a6 = 0)."""

    def sectors(self, x, N=None):
        P = dict(a3.Objective.sectors(self, x, N=N))
        P["E6_raw"] = P["E6"]
        P["E6"] = 0.0
        return P


def quartic_bounds():
    """Analytic floors for the a6=0, t=1 model (derivations in header)."""
    out = {}
    # Faddeev-type: E2+E4 >= 3 pi |B|   (exact constant in our conventions)
    out["faddeev_E2E4"] = 3.0 * PI
    # classic quartic+potential BPS-submodel bound, closed form
    from math import gamma
    out["quartic_pot_closed"] = (16.0 * 2.0 ** 0.25 * gamma(7.0 / 4.0)
                                 * gamma(3.0 / 2.0) / gamma(13.0 / 4.0))
    # same by quadrature: 4 INT sin^2 psi (1-cos psi)^{1/4} dpsi
    n = 4000
    xp, wp = leggauss(n)
    psi = 0.5 * PI * (xp + 1.0)
    w = 0.5 * PI * wp
    out["quartic_pot_quad"] = 4.0 * float(
        np.sum(np.sin(psi) ** 2 * (1.0 - np.cos(psi)) ** 0.25 * w))
    # saturated quartic floor: (2^{-1/4}/pi) INT_{S3} W~^{1/4} dOmega,
    # W~ = (1-cos psi)^2 + sin^2 psi cos^2 Theta
    up, wu = leggauss(n // 4)
    TH = 0.5 * PI * (up + 1.0)
    wth = 0.5 * PI * wu
    Ps, Th = np.meshgrid(psi, TH, indexing="ij")
    WW = np.outer(w, wth)
    Wt = (1.0 - np.cos(Ps)) ** 2 + np.sin(Ps) ** 2 * np.cos(Th) ** 2
    integ = float(np.sum(np.sin(Ps) ** 2 * np.sin(Th) * Wt ** 0.25 * WW))
    out["saturated_floor"] = 2.0 ** (-0.25) / PI * (2.0 * PI * integ)
    return out


def quartic_model():
    """The a6 = 0, t = 1 model end-to-end (reading b)."""
    print("\n== READING (b): pure quartic limit a6 = 0, E = E2+E4+E0 "
          "(t = 1 canonical: with a6 = 0, t scales out as t^{3/2}) ==")
    res = {}
    # gate GQ1: assemble_gen == rs.assemble at a6 = a0 = 1, t = TFROZEN
    r = rs.make_grid(3000, 6.0)
    f = rs.init_profile(r)
    Ea, sa, ga, da, oa = assemble_gen(r, f, ax.TFROZEN, ax.TFROZEN, 1.0, 1.0)
    Eb, sb, gb, db, ob = rs.assemble(r, f, ax.TFROZEN)
    res["gate_GQ1"] = dict(dE=abs(Ea - Eb),
                           dg=float(np.max(np.abs(ga - gb))),
                           dH=float(np.max(np.abs(da - db))))
    print("   gate GQ1 (assemble_gen vs radial_solve.assemble): dE=%.1e "
          "dgrad=%.1e dHess=%.1e" % (res["gate_GQ1"]["dE"],
                                     res["gate_GQ1"]["dg"],
                                     res["gate_GQ1"]["dH"]))
    # radial solve of the quartic model (N-ladder)
    prof = {}
    for N, rmax in ((3000, 12.0), (6000, 12.0), (6000, 16.0)):
        r = rs.make_grid(N, rmax)
        f, E, sect, gmax = solve_profile_gen(r, rs.init_profile(r),
                                             1.0, 1.0, 0.0, 1.0)
        prof[(N, rmax)] = (r, f, E, sect, gmax)
        E2, E4, E6, E0 = sect
        print("   radial solve N=%d rmax=%.0f: E=%.8f (E2 %.5f E4 %.5f "
              "E0 %.5f; E6[coeff 0] %.5f)  gmax=%.1e  virial "
              "tE2-tE4+3E0=%.2e"
              % (N, rmax, E, E2, E4, E0, E6, gmax, E2 - E4 + 3.0 * E0))
    r, f, E, sect, gmax = prof[(6000, 16.0)]
    res["radial"] = dict(E=E, sect=list(sect), gmax=gmax,
                         ladder={"%d_%g" % k: v[2] for k, v in prof.items()})
    fspl = ax.CubicSpline1D(r, f)
    dom_tail = float(r[np.nonzero(f > 2e-5)[0][-1]]) * 1.05
    r_sup_b = float(r[np.nonzero(f > 0.02)[0][-1]])
    res["dom_tail"], res["r_sup_b"] = dom_tail, r_sup_b

    # direct closure, family A (hedgehog + dilation), E6 zeroed
    P1 = ax.sector_pieces(ax.make_famA(fspl, 1.0, 1.0), N_FAM, dom_tail,
                          dom_tail)
    res["gate_engine_vs_radial"] = dict(
        E2=P1["E2"] / sect[0] - 1.0, E4=P1["E4"] / sect[1] - 1.0,
        E0=P1["E0"] / sect[3] - 1.0, B=P1["BDEG"])
    print("   engine vs radial sectors: E2 %+.1e  E4 %+.1e  E0 %+.1e  "
          "B=%.5f" % (res["gate_engine_vs_radial"]["E2"],
                      res["gate_engine_vs_radial"]["E4"],
                      res["gate_engine_vs_radial"]["E0"], P1["BDEG"]))
    P1q = dict(P1)
    P1q["E6"] = 0.0
    FA = ax.famA_funcs(P1q)
    solA = closure_direct(FA, 1.0, 600.0)
    res["directA"] = solA
    print("   direct closure (famA): c_paper=%.4f  kappa_paper=%.4f  "
          "E_stat=%.4f  (E2+E4 measured %.4f vs Faddeev floor 3pi=%.4f)"
          % (solA["c_paper"], solA["kappa_paper"], solA["Estat"],
             solA["E2"] + solA["E4"], 3.0 * PI))

    # saturated closure over the 12-mode basis (E6 zeroed)
    basis = a3.GenBasis(r_sup_b, K=3, a_l=(0, 2), b_l=(1, 2))
    active = list(range(basis.n))
    obj = QuarticObjective(fspl, basis, 1.0, N_OPT, dom_tail, active, 0.50,
                           mode="G")
    steps = np.array([0.06 if i < basis.na else 0.15 for i in active])
    seeds = [np.zeros(basis.n),
             a3.seed_vec(basis, active, a_map={2: 0.12}),
             a3.seed_vec(basis, active, b_map={1: 0.45, 2: 0.20})]
    x, fv = a3.optimize(obj, seeds, steps, B_TRI, B_COLD)
    objP = QuarticObjective(fspl, basis, 1.0, N_POL, dom_tail, active, 0.50,
                            mode="G")
    x, fv, _ = a3.nelder_mead(objP, x, steps * 0.08, maxfev=B_POL)
    P0N = dict(ax.sector_pieces(a3.make_qfun(fspl, basis,
                                             np.zeros(basis.na),
                                             np.zeros(basis.nb)),
                                N_FIN, dom_tail, dom_tail))
    P0N["E6"] = 0.0
    repq = a3.finalize_G(obj, x, N_FIN, P0N)
    P0C = dict(ax.sector_pieces(a3.make_qfun(fspl, basis,
                                             np.zeros(basis.na),
                                             np.zeros(basis.nb)),
                                N_COARSE, dom_tail, dom_tail))
    P0C["E6"] = 0.0
    repq["coarse"] = a3.finalize_G(obj, x, N_COARSE, P0C)
    repq["conv_c_rel"] = abs(repq["coarse"]["c_paper"] / repq["c_paper"]
                             - 1.0)
    obj0 = QuarticObjective(fspl, basis, 1.0, N_OPT, dom_tail, [], 0.50,
                            mode="G")
    rep0q = a3.finalize_G(obj0, np.zeros(0), N_FIN, P0N)
    res["saturated"] = repq
    res["saturated_hedgehog_core"] = rep0q
    print("   saturated closure: G_q*(12-mode)=%.5f -> c=%.4f  "
          "[hedgehog-core %.4f]  conv %.1e  B=%.5f"
          % (repq["G"], repq["c_paper"], rep0q["c_paper"],
             repq["conv_c_rel"], repq["degree"]))
    bnd = quartic_bounds()
    res["bounds"] = bnd
    print("   analytic floors: E2+E4 >= 3pi = %.5f;  E4+E0 >= %.5f "
          "(closed form; quad check %.5f);  G_q >= %.5f (saturated floor)"
          % (bnd["faddeev_E2E4"], bnd["quartic_pot_closed"],
             bnd["quartic_pot_quad"], bnd["saturated_floor"]))
    print("   => c_q(reading b) floors: direct >= (unit-map) ...; "
          "saturated c >= (4/pi)*3pi = 12 exactly.")
    res["c_floor_saturated"] = 12.0
    return res


# ----------------------------------------------------------------------
# gates (frozen references)
# ----------------------------------------------------------------------
def gates():
    out = {}
    print("== GATES (frozen references) ==")
    # G-A: engine on the G3 exact minimiser vs 16 sqrt2/9
    N = 448 if not QUICK else 224
    P = ax.sector_pieces(gs.qfun_min, N, 1.06 * gs.R_EQ, 1.06 * gs.R_EQ)
    G = P["E6"] + P["E0"] - P["I"] / SIXTEEN_PI
    out["GA_gstar"] = dict(N=N, G=G, err=G - GSTAR, BDEG=P["BDEG"])
    print("   G-A engine@G3-minimiser N=%d: G=%.9f (G-G* = %+.2e; gstar "
          "table at 448: -1.13e-5)  B=%.6f" % (N, G, G - GSTAR, P["BDEG"]))
    # G-B: dial at eps=0.05 must reproduce the frozen t
    d = dial(0.05, ax.TFROZEN)
    out["GB_dial"] = dict(t=d["t"], t_frozen=ax.TFROZEN,
                          rel=d["t"] / ax.TFROZEN - 1.0, ratio=d["ratio"])
    print("   G-B dial(0.05): t=%.9g vs frozen %.9g (rel %+.1e), "
          "ratio=%.6f" % (d["t"], ax.TFROZEN, out["GB_dial"]["rel"],
                          d["ratio"]))
    # G-C: frozen axi/axi3 reference numbers for eps=0.05 (read from JSON)
    with open(T2DIR + "/axi_results.json") as fh:
        AJ = json.load(fh)
    with open(T2DIR + "/axi3_results.json") as fh:
        A3J = json.load(fh)
    out["refs"] = dict(solA_c=AJ["solA"]["c_paper"],
                       solB_c=AJ["solB"]["c_paper"],
                       axi3_R2_sat_c=A3J["R2"]["saturated"]["c_paper"],
                       axi3_R2_sat_hh_c=A3J["R2"]["saturated_hedgehog"]
                       ["c_paper"])
    print("   G-C frozen refs: step-2 solA c=%.5f, solB c=%.5f, axi3 R2 "
          "saturated c=%.5f (hedgehog core %.5f)"
          % (out["refs"]["solA_c"], out["refs"]["solB_c"],
             out["refs"]["axi3_R2_sat_c"], out["refs"]["axi3_R2_sat_hh_c"]))
    return out


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    out = dict(config=dict(quick=QUICK, eps_list=EPS_LIST,
                           N_OPT=N_OPT, N_POL=N_POL, N_FIN=N_FIN,
                           N_COARSE=N_COARSE, N_FAM=N_FAM,
                           B_TRI=B_TRI, B_COLD=B_COLD, B_POL=B_POL),
               constants=dict(GSTAR=GSTAR, CSAT0=CSAT0, R5=[R5_C, R5_SIG],
                              C0_corpus=C0_CORPUS, rung=RUNG,
                              pull_CSAT0=pull(CSAT0)))
    out["gates"] = gates()

    # ---------------- eps ladder ----------------
    print("\n== EPS-DIAL LADDER (campaign convention: "
          "eps = t(E2+E4)/(E6+E0) on the step-1 hedgehog solve) ==")
    rows = []
    t_guess = ax.TFROZEN * EPS_LIST[0] / 0.05
    for eps in EPS_LIST:
        sol = dial(eps, t_guess, rmax=rmax_for(eps))
        t_guess = sol["t"] * (EPS_LIST[min(EPS_LIST.index(eps) + 1,
                                           len(EPS_LIST) - 1)] / eps)
        lad = dial_ladder(eps, sol)
        E2, E4, E6, E0 = sol["sect"]
        Estat = sol["t"] * (E2 + E4) + E6 + E0
        rows.append(dict(eps=eps, t=sol["t"], ratio=sol["ratio"],
                         sect=sol["sect"], Estat_rad=Estat,
                         kin_frac=sol["ratio"] / (1.0 + sol["ratio"]),
                         mu_true=1.0 / math.sqrt(2.0 * sol["t"]),
                         gmax=sol["gmax"], rmax=sol["rmax"],
                         ladder=lad, _r=sol["r"], _f=sol["f"]))
        print("   eps=%.3f: t=%.7g  achieved=%.5f  kin_frac=%.4f  "
              "E_stat_rad=%.5f  mu=%.3f  [conv: dratio N2x %+.1e, "
              "rmax+ %+.1e]"
              % (eps, sol["t"], sol["ratio"], rows[-1]["kin_frac"],
                 Estat, rows[-1]["mu_true"],
                 lad["N2x"]["d_ratio_rel"], lad["rmax+"]["d_ratio_rel"]))
    out["dial"] = [{k: v for k, v in r.items()
                    if not k.startswith("_")} for r in rows]

    # ---------------- closures per eps (reading a) ----------------
    print("\n== READING (a): eps-dial continuation — direct (saddle) and "
          "saturated branches ==")
    lam_grid = np.arange(0.15, 1.4501, 0.05)
    curve = []
    for row in rows:
        eps, t = row["eps"], row["t"]
        r, f = row["_r"], row["_f"]
        fspl = ax.CubicSpline1D(r, f)
        idx = np.nonzero(f > 2e-5)[0]
        dom_tail = float(min(r[idx[-1]] * 1.05, r[-1]))
        r_sup_b = float(r[np.nonzero(f > 0.02)[0][-1]])
        dom_sat = float(min(max(1.5 * r_sup_b, dom_tail), r[-1]))
        Lhi = 60.0 * (1.0 + eps) ** 1.5
        tA = time.time()
        solA, solB, P1 = direct_pair(fspl, t, dom_tail, N_FAM, lam_grid,
                                     Lhi)
        # coarse-grid rerun of the direct pair (convergence)
        solAc, solBc, _ = direct_pair(fspl, t, dom_tail, N_FAM_C, lam_grid,
                                      Lhi)
        convA = abs(solAc["c_paper"] / solA["c_paper"] - 1.0)
        convB = abs(solBc["c_paper"] / solB["c_paper"] - 1.0)
        print("   eps=%.3f DIRECT: famA c=%.5f (kap=%.4f, conv %.1e) | "
              "famB c=%.5f (kap=%.4f, lam*=%.3f, conv %.1e)  (%.0f s)"
              % (eps, solA["c_paper"], solA["kappa_paper"], convA,
                 solB["c_paper"], solB["kappa_paper"], solB["lam"], convB,
                 time.time() - tA))
        tS = time.time()
        repS, rep0 = saturated_at(fspl, t, r_sup_b, dom_sat,
                                  "eps=%.3f" % eps)
        floor_c = CSAT0 + 12.0 * t                # rigorous analytic floor
        curve.append(dict(
            eps=eps, t=t,
            c_directA=solA["c_paper"], kappa_directA=solA["kappa_paper"],
            c_directA_conv=convA,
            c_directB=solB["c_paper"], kappa_directB=solB["kappa_paper"],
            lamB=solB["lam"], c_directB_conv=convB,
            directA=solA, directB=solB,
            c_sat=repS["c_paper"], G_sat=repS["G"],
            c_sat_hedgehog_core=rep0["c_paper"],
            c_sat_conv=repS["conv_c_rel"],
            sat=repS, sat_core0=rep0,
            c_sat_floor_analytic=floor_c,
            floor_ok=bool(repS["c_paper"] >= floor_c - 5e-3),
            pull_sat=pull(repS["c_paper"]),
            pull_directA=pull(solA["c_paper"]),
            pull_directB=pull(solB["c_paper"]),
            secs_direct=time.time() - tA - (time.time() - tS),
            secs_sat=time.time() - tS))
        if not curve[-1]["floor_ok"]:
            print("   !! eps=%.3f saturated value BELOW analytic floor "
                  "%.5f — investigate" % (eps, floor_c))
    out["curve"] = [{k: v for k, v in c.items()
                     if k not in ("directA", "directB", "sat", "sat_core0")}
                    for c in curve]
    out["curve_full"] = [dict(eps=c["eps"],
                              directA=c["directA"], directB=c["directB"],
                              sat=c["sat"], sat_core0=c["sat_core0"])
                         for c in curve]

    # ---------------- reading (b) ----------------
    out["quartic_model"] = quartic_model()

    # ---------------- verdict tables ----------------
    print("\n== VERDICT TABLE (reading a; frozen compacton unit map; "
          "<r5>: c_q = 3.1 +/- 0.2; saturated anchor 64sqrt2/9pi = "
          "%.6f) ==" % CSAT0)
    print("   %-6s %-9s | %-9s %-7s | %-9s %-7s | %-9s %-7s | %s"
          % ("eps", "t", "c_dirA", "pull", "c_dirB", "pull", "c_sat",
             "pull", "c_sat floor"))
    for c in curve:
        print("   %-6.3f %-9.3g | %-9.5f %+-7.2f | %-9.5f %+-7.2f | "
              "%-9.5f %+-7.2f | %.4f"
              % (c["eps"], c["t"], c["c_directA"], c["pull_directA"],
                 c["c_directB"], c["pull_directB"], c["c_sat"],
                 c["pull_sat"], c["c_sat_floor_analytic"]))

    # normalisation-reading enumeration (the 'scaling audit' question)
    print("\n== NORMALISATION READINGS (candidate c_q conventions at "
          "eps >= 1; measured, no fitting) ==")
    print("   N0 raw c(eps); N1 c/(1+eps); N2 c/E_growth "
          "(E_stat_rad(eps)/E_stat_rad(eps->0), rad anchor 32sqrt2/15*"
          "(1+eps) exact at 0); N3 eps->0 intercept of branch")
    E0_anchor = 32.0 * SQ2 / 15.0     # E6+E0 at the compacton (BPS bound)
    readings = []
    for c in curve:
        row = next(rr for rr in rows if rr["eps"] == c["eps"])
        Egrow = row["Estat_rad"] / E0_anchor
        rd = dict(eps=c["eps"],
                  sat_N0=c["c_sat"], sat_N1=c["c_sat"] / (1.0 + c["eps"]),
                  sat_N2=c["c_sat"] / Egrow,
                  dirB_N0=c["c_directB"],
                  dirB_N1=c["c_directB"] / (1.0 + c["eps"]),
                  dirB_N2=c["c_directB"] / Egrow,
                  E_growth=Egrow)
        for k in ("sat_N0", "sat_N1", "sat_N2", "dirB_N0", "dirB_N1",
                  "dirB_N2"):
            rd["pull_" + k] = pull(rd[k])
        readings.append(rd)
        print("   eps=%.3f E_growth=%.4f | sat: N0 %.4f (%+.1f s) N1 %.4f "
              "(%+.1f s) N2 %.4f (%+.1f s) | dirB: N0 %.4f (%+.1f s) N1 "
              "%.4f (%+.1f s) N2 %.4f (%+.1f s)"
              % (rd["eps"], Egrow, rd["sat_N0"], rd["pull_sat_N0"],
                 rd["sat_N1"], rd["pull_sat_N1"], rd["sat_N2"],
                 rd["pull_sat_N2"], rd["dirB_N0"], rd["pull_dirB_N0"],
                 rd["dirB_N1"], rd["pull_dirB_N1"], rd["dirB_N2"],
                 rd["pull_dirB_N2"]))
    out["readings"] = readings

    # branch crossings of the r5 band
    eps_arr = np.array([c["eps"] for c in curve])
    for name, key in (("directA", "c_directA"), ("directB", "c_directB"),
                      ("saturated", "c_sat")):
        cv = np.array([c[key] for c in curve])
        cross = None
        for i in range(len(cv) - 1):
            if (cv[i] - R5_C) * (cv[i + 1] - R5_C) < 0.0:
                w = (R5_C - cv[i]) / (cv[i + 1] - cv[i])
                cross = float(eps_arr[i] + w * (eps_arr[i + 1] - eps_arr[i]))
        out.setdefault("crossings", {})[name] = cross
        print("   branch %-9s crosses c = 3.1 at eps ~ %s (linear interp "
              "on the grid)" % (name, "%.3f" % cross if cross else "never"))

    out["runtime_s"] = time.time() - t0
    with open(OUTDIR + "/i1_results.json", "w") as fh:
        json.dump(a3.jsonable(out), fh, indent=2)
    print("\nruntime %.1f s — wrote i1_results.json" % out["runtime_s"])


if __name__ == "__main__":
    main()
