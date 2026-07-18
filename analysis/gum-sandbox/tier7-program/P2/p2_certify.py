#!/usr/bin/env python3
"""
WORKSTREAM P2 -- T4-W5 certification pilot: exact-field, high-precision
backward paths (F-T7-P2).  Upgrades O1 (F-T7-O1) one rung toward certified.

KEY IDEA (pre-registered, ROADMAP_v9 Phase P / P2).  The L5'/L7b/O1 engine's
field is a FINITE trigonometric sum: the state is exactly factorized,
    phi(x, z, t) = f(x, t) sin(pi z) e^{-i pi^2 t / 2},
and f(., t) is the free evolution of the Nx = 2048 grid datum on the
periodic box [-15, 45].  The FFT coefficients of that datum ARE the exact
spectral coefficients of the band-limited field the engine evolves:
    f(x, t) = sum_k  C_k  e^{i kappa_k (x - xmin)}  e^{-i kappa_k^2 t / 2},
    C_k = F0[k] / Nx,   kappa_k = (pi/30) m_k,  m_k in [-1024, 1023],
which the grid engine reproduces EXACTLY at grid nodes (ifft of the phased
spectrum).  Hence f, f' -- and rho, j, v -- can be evaluated exactly (to
floating-point rounding) at ARBITRARY (x, z, t) by direct mode summation:
the x-part is the 2048-mode sum above, the z-part is analytic
(rho = 2|f|^2 sin^2(pi z);  v_x = Im(f'/f) - n_y pi cot(pi z),
v_z = n_y Re(f'/f)).  This removes O1's dominant error source -- linear
interpolation of grid velocity fields -- entirely.  The engine's velocity
regularization (a, b clipped to +-500; denominator floored at
1e-14 x max_x rho1d(., t); z clipped to [1e-9, 1-1e-9] inside the cot) is
mirrored EXACTLY and monitored: on every certified path we record whether
any clip/floor ever binds (none does on class-(a) paths).

WHAT IS RE-RUN (stratified subset of O1's 64,000 kill-bin backward paths,
per family n_y in {+1, +0.75}; the -n_y mirrors are numerically identical
to the +n_y families under z -> 1-z, verified in O1 and re-verified here
from o1_raw.npz at the 3.5e-8 level):
  - ALL 16 box-exit/wall-adjacent paths per family: every O1 BOX_EXIT path
    (5 for |n_y| = 1, 3 for 0.75) + the extreme-wall-row (z = 0.0025 or
    0.9975) class-(a) paths with the EARLIEST certified crossing epoch
    s_pc (longest certified segment = hardest wall cases), filled to 16;
  - the 100 smallest-margin class-(a) paths per family;
  - 100 random class-(a) paths per family (seed 20260718);
  = 216 per family, 432 total.
Integration: scipy DOP853, rtol 1e-12, atol 1e-14 (adaptive; the adaptive
step ~7e-4 resolves the fastest spectral beat kappa_max^2/2 ~ 5.7e3 rad/t),
dense output; margin = max_s X_x - d refined by local optimization on the
order-7 dense interpolant; s_pc by root-finding X_x = d + 1e-3.

PRECISION LADDER (>= 20 paths: the worst margins + a spread; 24 here):
re-integration in 200-bit (~60 decimal digit >= 50-digit) GMP/MPFR
arithmetic via gmpy2 -- the same backend mpmath uses; a cross-check of one
field evaluation against mpmath at dps = 60 is recorded.  Method: classical
fixed-step RK4, h = t0/ceil(t0/5e-4) (~5e-4; float64 calibration vs the
adaptive reference shows the h = 5e-4 endpoint error is 1e-10 class), with
the exact per-mode phase advanced by precomputed half-step recurrences
(no per-step trigonometry).  Step-halving (h -> 2.5e-4) on the 3
worst-margin paths bounds the mp-side truncation independently of float64.
Margins refined by local h/50 re-integration + 3-point parabola; s_pc by
h/50 bracket + linear interpolation.

GATES (pre-registered, ROADMAP_v9 Phase P / P2):
  G1 exact-vs-grid field agreement at nodes <= 1e-12 rel, >= 1000 random
     nodes x several times (both f and the velocity components a, b;
     "rel" = sup-normalized over all sampled nodes AND per-node relative
     on well-conditioned nodes).
  G2 every re-run path keeps its O1 classification; margins agree with O1
     within O1's own error estimate (2.8e-5 class at p95): pre-registered
     numeric criteria  p95|dm| <= 5.6e-5 (2x O1's Nx-doubling p95) and
     max|dm| <= 1.2e-2 (2x O1's max refinement diff), over the paths O1
     recorded uncapped; capped paths must re-certify margin >= cap = 2.0.
  G3 precision ladder: float64-DOP853 vs 50-digit endpoints <= 1e-8;
     margins stable to >= 6 digits (rel diff <= 1e-6) on the ladder subset.
  G4 the worst-margin path certified at ladder grade (full record printed).
  G5 honest scope: precision-certification, NOT formal interval arithmetic;
     what formal certification still requires; the analytic proof remains
     the terminal open.

Engine constants inherited verbatim from o1_run.py / l7b_run.py /
l5prime_run.py.  hbar = m = 1.  Usage: p2_certify.py [g1|select|dop|ladder|
gates|fig|all]  (default: all remaining phases, checkpointed).
"""
import json
import os
import sys
import time

import numpy as np
from scipy import fft as sfft
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar

# ----------------------------------------------------------------------
# Parameters (inherited verbatim from O1 / L7b / L5prime)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
VCLAMP = 500.0
RHO_EPS_REL = 1e-14
XMIN, XMAX = -15.0, 45.0
NX = 2048                       # O1 reference resolution: its grid DEFINES
                                # the band-limited field we evaluate exactly
KILL_BINS = [(6.4, 6.8), (6.8, 7.2)]
NY_FAMS = [1.0, 0.75]           # mirrors -1, -0.75 numerically identical
RHO_FLOOR_REL = 1e-10
MARGIN_EPS = 1e-9
MARGIN_CLASS = 1e-3
MARGIN_CAP = 2.0
ZCLIP = 1e-9

N_WALL, N_SMALL, N_RAND = 16, 100, 100    # per family -> 216 x 2 = 432
SEED_SELECT = 20260718

RTOL, ATOL = 1e-12, 1e-14                 # DOP853 (pre-registered)
H_MP = 5e-4                               # mp ladder nominal step
MP_BITS = 200                             # ~60 decimal digits (>= 50)
N_LADDER = 24
N_HALVE = 3                               # step-halving checks (worst margins)

OUTDIR = os.path.dirname(os.path.abspath(__file__))
O1DIR = os.path.join(os.path.dirname(OUTDIR), "O1")

CLS_NAME = {0: "A_PRECROSSED", 1: "B_RHOFLOOR", 2: "C_FRESH_IN_SUPP",
            3: "FRESH_OUT_SUPP", 4: "BOX_EXIT", 5: "GRAZE_UNCERTAIN"}


# ----------------------------------------------------------------------
# Initial profile (verbatim) and exact spectral representation
# ----------------------------------------------------------------------
def smoothstep(t):
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def taper(x):
    u = np.abs(np.asarray(x, float) - X0) / SIGMA_X
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    return W


def amp(x):
    d2 = (np.asarray(x, float) - X0) ** 2 / SIGMA_X ** 2
    return np.exp(-d2 / 2.0) * taper(x)


def build_spectrum():
    dx = (XMAX - XMIN) / NX
    xg = XMIN + dx * np.arange(NX)
    f0 = amp(xg).astype(complex) * np.exp(1j * K0 * xg)
    f0 /= np.sqrt(np.sum(np.abs(f0) ** 2) * dx)
    kx = 2.0 * np.pi * sfft.fftfreq(NX, dx)
    F0 = sfft.fft(f0)
    return xg, dx, f0, kx, F0


XG, DX, F0_GRID, KX, F0 = build_spectrum()
CFT = F0 / NX                    # exact spectral coefficients C_k
IKC = 1j * KX * CFT              # coefficients of f'

# rho1d max table (for the engine's relative denominator floor; the floor
# never binds on certified paths -- monitored -- but is mirrored exactly)
_TS = np.arange(0.0, 7.2 + 2e-3, 1e-3)
_RMAX = np.array([(np.abs(sfft.ifft(F0 * np.exp(-0.5j * KX ** 2 * t))) ** 2
                   ).max() for t in _TS])


def fields_exact(x, t):
    """Exact band-limited f, f' at arbitrary (x, t) by direct mode sum."""
    ph = np.exp(1j * (KX * (x - XMIN) - 0.5 * KX ** 2 * t))
    return (CFT * ph).sum(), (IKC * ph).sum()


def rhs_factory(nyv):
    """Backward-ODE right-hand side on the exact field, with the engine's
    regularization mirrored exactly (clips/floor; they never bind on
    class-(a) paths -- monitored post hoc)."""
    def rhs(s, y):
        x, z = y
        ze = min(max(z, ZCLIP), L - ZCLIP)
        f, fp = fields_exact(x, s)
        r = (f.real * f.real + f.imag * f.imag)
        den = max(r, RHO_EPS_REL * np.interp(s, _TS, _RMAX))
        q = np.conj(f) * fp
        a = min(max(q.imag / den, -VCLAMP), VCLAMP)
        b = min(max(q.real / den, -VCLAMP), VCLAMP)
        vx = min(max(a - nyv * np.pi / np.tan(np.pi * ze / L), -VCLAMP),
                 VCLAMP)
        return np.array([vx, nyv * b])
    return rhs


# ----------------------------------------------------------------------
# Phase G1: exact evaluator vs grid engine at grid nodes
# ----------------------------------------------------------------------
def phase_g1():
    """Exact-vs-grid agreement at grid nodes.  Both engines compute the
    SAME finite trigonometric sum (FFT synthesis vs direct summation), so
    the difference is pure floating-point rounding.  Metric (fixed here,
    before the production run): (i) sup-normalized field agreement
    |df| / max|f| <= 1e-12 over ALL sampled nodes; (ii) per-node RELATIVE
    agreement <= 1e-12 for f on well-conditioned nodes (|f| >= 3e-2 max)
    and for the velocity components a, b on well-conditioned nodes
    (rho >= 3e-2 max, velocity scale max(1, |a|)); (iii) on ALL nodes the
    a, b deviation must be explained by rounding amplification: the
    conditioning-normalized deviation |da| * (rho / rho_max) / max(1,|a|)
    <= 1e-12 (near-floor nodes amplify the ~1e-15 rounding by
    rho_max / rho -- the same regularized regime O1 flux-bounds rather
    than trusts).  Per-time uniform sample 1200 nodes + 300 targeted
    well-conditioned nodes, 6 times."""
    rng = np.random.default_rng(11)
    times = [0.5, 1.7, 3.3, 5.1, 6.6, 7.195]
    per_time = []
    for t in times:
        Ft = F0 * np.exp(-0.5j * KX ** 2 * t)
        fg = sfft.ifft(Ft)
        fpg = sfft.ifft(1j * KX * Ft)
        rho = np.abs(fg) ** 2
        rmax = rho.max()
        den = np.maximum(rho, RHO_EPS_REL * rmax)
        w = np.conj(fg) * fpg
        ag = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
        bg = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
        fmax = np.abs(fg).max()
        idx_u = rng.integers(0, NX, 1200)
        well_pool = np.where(rho >= 3e-2 * rmax)[0]
        idx_w = rng.choice(well_pool, 300, replace=True)
        idx = np.concatenate([idx_u, idx_w])
        sup_f = rel_f = rel_ab = cond_ab = 0.0
        n_well_f = n_well_ab = 0
        for i in idx:
            fe, fpe = fields_exact(XG[i], t)
            re_ = (fe.real ** 2 + fe.imag ** 2)
            dene = max(re_, RHO_EPS_REL * rmax)
            qe = np.conj(fe) * fpe
            ae = min(max(qe.imag / dene, -VCLAMP), VCLAMP)
            be = min(max(qe.real / dene, -VCLAMP), VCLAMP)
            da = max(abs(ae - ag[i]), abs(be - bg[i])) / max(1.0, abs(ag[i]))
            sup_f = max(sup_f, abs(fe - fg[i]) / fmax)
            cond_ab = max(cond_ab, da * min(1.0, rho[i] / rmax))
            if abs(fg[i]) >= 3e-2 * fmax:
                rel_f = max(rel_f, abs(fe - fg[i]) / abs(fg[i]))
                n_well_f += 1
            if rho[i] >= 3e-2 * rmax:
                rel_ab = max(rel_ab, da)
                n_well_ab += 1
        per_time.append(dict(t=t, n_nodes=int(idx.size),
                             n_well_f=n_well_f, n_well_ab=n_well_ab,
                             sup_rel_f=float(sup_f),
                             pernode_rel_f_well=float(rel_f),
                             pernode_rel_ab_well=float(rel_ab),
                             cond_normalized_ab_all=float(cond_ab)))
        print("  t=%.3f  sup|df|/max|f|=%.2e  rel f(well,n=%d)=%.2e  "
              "rel ab(well,n=%d)=%.2e  cond-norm ab(all)=%.2e"
              % (t, sup_f, n_well_f, rel_f, n_well_ab, rel_ab, cond_ab),
              flush=True)
    worst = dict(
        sup_rel_f=max(p["sup_rel_f"] for p in per_time),
        pernode_rel_f_well=max(p["pernode_rel_f_well"] for p in per_time),
        pernode_rel_ab_well=max(p["pernode_rel_ab_well"] for p in per_time),
        cond_normalized_ab_all=max(p["cond_normalized_ab_all"]
                                   for p in per_time),
        n_nodes_total=sum(p["n_nodes"] for p in per_time),
        n_well_ab_total=sum(p["n_well_ab"] for p in per_time))
    ok = (worst["sup_rel_f"] <= 1e-12
          and worst["pernode_rel_f_well"] <= 1e-12
          and worst["pernode_rel_ab_well"] <= 1e-12
          and worst["cond_normalized_ab_all"] <= 1e-12)
    out = dict(times=times, per_time=per_time, worst=worst,
               metric=("(i) sup-normalized f over all nodes; (ii) per-node "
                       "relative f and a,b on well-conditioned nodes "
                       "(3e-2 cut); (iii) conditioning-normalized a,b "
                       "deviation on all nodes; each <= 1e-12"),
               verdict="PASS" if ok else "FAIL")
    json.dump(out, open(os.path.join(OUTDIR, "_g1.json"), "w"), indent=2)
    print("G1:", out["verdict"], json.dumps(worst), flush=True)
    return out


# ----------------------------------------------------------------------
# Phase SELECT: stratified subset from O1's raw classification
# ----------------------------------------------------------------------
def phase_select():
    raw = np.load(os.path.join(O1DIR, "o1_raw.npz"))
    ny, tf, zf = raw["ny_flat"], raw["t_flat"], raw["z_flat"]
    cls, mar, spc = raw["cls_ref"], raw["margin_ref"], raw["s_pc_ref"]

    # mirror re-verification (+ny vs -ny under z -> 1-z)
    mirror = {}
    for fam in NY_FAMS:
        def keyed(nyv):
            m = ny == nyv
            zz = zf[m] if nyv > 0 else 1.0 - zf[m]
            order = np.lexsort((np.round(zz, 6), np.round(tf[m], 6)))
            return mar[m][order], cls[m][order]
        m1, c1 = keyed(fam)
        m2, c2 = keyed(-fam)
        mirror["%+.2f" % fam] = dict(
            max_abs_dmargin=float(np.nanmax(np.abs(m1 - m2))),
            cls_identical=bool((c1 == c2).all()))

    rng = np.random.default_rng(SEED_SELECT)
    sel_idx, sel_stratum = [], []
    for fam in NY_FAMS:
        famsel = ny == fam
        # -- wall-adjacent stratum: all BOX_EXIT + earliest-s_pc wall rows
        box = np.where(famsel & (cls == 4))[0]
        wallrow = np.where(famsel & (cls == 0)
                           & ((zf == zf.min()) | (zf == zf.max())))[0]
        wall_fill = wallrow[np.argsort(spc[wallrow])][:N_WALL - box.size]
        wall = np.concatenate([box, wall_fill])
        # -- smallest-margin class-(a) stratum
        a_all = np.where(famsel & (cls == 0))[0]
        a_rest = np.setdiff1d(a_all, wall, assume_unique=False)
        small = a_rest[np.argsort(mar[a_rest])][:N_SMALL]
        # -- random class-(a) stratum
        a_rest2 = np.setdiff1d(a_rest, small, assume_unique=False)
        rand = rng.choice(a_rest2, N_RAND, replace=False)
        for arr, name in ((wall, "wall"), (small, "small"), (rand, "rand")):
            sel_idx.extend(arr.tolist())
            sel_stratum.extend([name] * arr.size)
        print("  family %+.2f: wall=%d (box-exit %d) small=%d rand=%d"
              % (fam, wall.size, box.size, small.size, rand.size), flush=True)
    sel_idx = np.array(sel_idx)
    np.savez_compressed(
        os.path.join(OUTDIR, "_selection.npz"),
        idx=sel_idx, stratum=np.array(sel_stratum),
        ny=ny[sel_idx], t0=tf[sel_idx], z0=zf[sel_idx],
        cls_o1=cls[sel_idx], margin_o1=mar[sel_idx], s_pc_o1=spc[sel_idx])
    json.dump(dict(mirror=mirror, n_selected=int(sel_idx.size)),
              open(os.path.join(OUTDIR, "_select.json"), "w"), indent=2)
    print("selected %d paths; mirror check: %s"
          % (sel_idx.size, json.dumps(mirror)), flush=True)


# ----------------------------------------------------------------------
# Phase DOP: float64 DOP853 backward integration of the subset
# ----------------------------------------------------------------------
def _ev_left(s, y):
    return y[0] - (XMIN + 0.5)


def _ev_right(s, y):
    return y[0] - (XMAX - 0.5)


_ev_left.terminal = True
_ev_right.terminal = True


def integrate_path(args):
    nyv, t0, z0 = args
    rhs = rhs_factory(nyv)
    tw = time.time()
    sol = solve_ivp(rhs, (t0, 0.0), [D_NEAR, z0], method="DOP853",
                    rtol=RTOL, atol=ATOL, dense_output=True,
                    events=[_ev_left, _ev_right])
    s_end = sol.t[-1]
    n_s = max(2, int(np.ceil((t0 - s_end) / 2.5e-4)) + 1)
    ss = np.linspace(t0, s_end, n_s)
    ys = sol.sol(ss)
    xs = ys[0]
    # refined margin: local optimization on the dense interpolant
    im = int(np.argmax(xs))
    lo, hi = ss[min(im + 1, n_s - 1)], ss[max(im - 1, 0)]
    if hi > lo:
        r = minimize_scalar(lambda s: -sol.sol(s)[0], bounds=(lo, hi),
                            method="bounded",
                            options=dict(xatol=1e-12))
        margin = max(xs[im], -r.fun) - D_NEAR
    else:
        margin = xs[im] - D_NEAR
    # s_pc: first backward crossing of d + 1e-3
    s_pc = np.nan
    above = xs >= D_NEAR + MARGIN_CLASS
    if above.any():
        j = int(np.argmax(above))          # first sample (largest s) above
        if j == 0:
            s_pc = ss[0]
        else:
            s_pc = brentq(lambda s: sol.sol(s)[0] - (D_NEAR + MARGIN_CLASS),
                          ss[j], ss[j - 1], xtol=1e-13)
    # monitors on a subsample: rho floor, clip headroom (z-clipped as engine)
    sub = ss[::4]
    min_rho_rel = np.inf
    s_floor = np.nan
    max_a = max_b = 0.0
    for s, x, z in zip(sub, sol.sol(sub)[0], sol.sol(sub)[1]):
        f, fp = fields_exact(x, s)
        r = f.real ** 2 + f.imag ** 2
        rm = np.interp(s, _TS, _RMAX)
        ze = min(max(z, ZCLIP), L - ZCLIP)
        rho_rel = r * np.sin(np.pi * ze) ** 2 / rm
        if rho_rel < min_rho_rel:
            min_rho_rel = rho_rel
        if np.isnan(s_floor) and rho_rel < RHO_FLOOR_REL:
            s_floor = s
        q = np.conj(f) * fp
        max_a = max(max_a, abs(q.imag / max(r, RHO_EPS_REL * rm)))
        max_b = max(max_b, abs(q.real / max(r, RHO_EPS_REL * rm)))
    # classification (O1 semantics)
    boxexit = sol.status == 1
    floored_first = (not np.isnan(s_floor)) and \
        (np.isnan(s_pc) or s_floor > s_pc)
    if not np.isnan(s_pc) and not floored_first:
        c = 0
    elif floored_first:
        c = 1
    elif boxexit:
        c = 4
    elif MARGIN_EPS < margin < MARGIN_CLASS:
        c = 5
    elif abs(sol.y[0, -1] - X0) <= W_OUT * SIGMA_X:
        c = 2
    else:
        c = 3
    return dict(margin=float(margin), s_pc=float(s_pc), cls=int(c),
                x_end=float(sol.y[0, -1]), z_end=float(sol.y[1, -1]),
                s_end=float(s_end), boxexit=bool(boxexit),
                min_rho_rel=float(min_rho_rel), s_floor=float(s_floor),
                max_abs_a=float(max_a), max_abs_b=float(max_b),
                nfev=int(sol.nfev), nsteps=int(sol.t.size - 1),
                elapsed=float(time.time() - tw))


def phase_dop():
    from multiprocessing import Pool
    sel = np.load(os.path.join(OUTDIR, "_selection.npz"))
    args = list(zip(sel["ny"], sel["t0"], sel["z0"]))
    tw = time.time()
    with Pool(4) as pool:
        recs = []
        for i, r in enumerate(pool.imap(integrate_path, args, chunksize=4)):
            recs.append(r)
            if (i + 1) % 40 == 0:
                print("  dop %d/%d  elapsed %.0f s" % (i + 1, len(args),
                      time.time() - tw), flush=True)
    out = {k: np.array([r[k] for r in recs]) for k in recs[0]}
    np.savez_compressed(os.path.join(OUTDIR, "_dop.npz"), **out)
    print("dop done: %.0f s total, mean nfev %.0f"
          % (time.time() - tw, out["nfev"].mean()), flush=True)


# ----------------------------------------------------------------------
# Phase LADDER: multiprecision (200-bit ~ 60 >= 50 digit) re-integration
# ----------------------------------------------------------------------
def ladder_selection():
    sel = np.load(os.path.join(OUTDIR, "_selection.npz"))
    dop = np.load(os.path.join(OUTDIR, "_dop.npz"))
    cls_ok = dop["cls"] == 0
    idx_all = np.arange(sel["idx"].size)
    # 12 globally smallest margins among class-(a)
    order = np.argsort(np.where(cls_ok, dop["margin"], np.inf))
    worst12 = order[:12]
    chosen = list(worst12)
    # 6 spread over family +1 class-(a) by s_pc quantiles
    fam1 = idx_all[(sel["ny"] == 1.0) & cls_ok & ~np.isin(idx_all, chosen)]
    q = np.quantile(dop["s_pc"][fam1], np.linspace(0, 1, 6))
    for v in q:
        j = fam1[np.argmin(np.abs(dop["s_pc"][fam1] - v))]
        if j not in chosen:
            chosen.append(j)
    # 6 spread over family +0.75 class-(a) by margin quantiles
    fam2 = idx_all[(sel["ny"] == 0.75) & cls_ok & ~np.isin(idx_all, chosen)]
    q = np.quantile(dop["margin"][fam2], np.linspace(0, 1, 6))
    for v in q:
        j = fam2[np.argmin(np.abs(dop["margin"][fam2] - v))]
        if j not in chosen:
            chosen.append(j)
    # top up (quantile collisions) with next-worst margins
    k = 12
    while len(chosen) < N_LADDER:
        j = order[k]
        if j not in chosen:
            chosen.append(j)
        k += 1
    return np.array(chosen[:N_LADDER])


def mp_run_path(args):
    """Fixed-step RK4 backward integration in 200-bit gmpy2 arithmetic on
    the exact mode sum (phase recurrences, no per-step trigonometry).
    Returns endpoint, refined margin, s_pc, all as strings + floats."""
    nyv, t0f, z0f, h_nom = args
    import gmpy2
    from gmpy2 import mpfr, mpc
    gmpy2.get_context().precision = MP_BITS
    tw = time.time()

    one = mpfr(1)
    pi = gmpy2.const_pi()
    base = pi / 30                       # kappa spacing (pi/30 exactly)
    xminm = mpfr(XMIN)
    nym = mpfr(nyv)
    t0 = mpfr(t0f)
    d_class = mpfr(D_NEAR) + mpfr(MARGIN_CLASS)

    mvals = list(range(0, NX // 2)) + list(range(-NX // 2, 0))
    Cm = [mpc(mpfr(c.real), mpfr(c.imag)) for c in CFT]
    IKm = [mpc(mpfr(0), base * m) for m in mvals]    # i kappa_k

    def phase_array(theta):
        """[e^{-i theta m^2} for m in mvals] via quadratic recurrence
        (u^{(m+1)^2} = u^{m^2} * u^{2m+1}; one exp total)."""
        u = gmpy2.exp(mpc(mpfr(0), -theta))
        u2 = u * u
        sq = [None] * (NX // 2 + 1)      # u^{m^2}, m = 0 .. NX/2
        q = mpc(1, 0)
        r = u
        for m in range(NX // 2 + 1):
            sq[m] = q
            q = q * r
            r = r * u2
        return [sq[abs(k if k < NX // 2 else k - NX)] for k in range(NX)]

    n = int(np.ceil(t0f / h_nom))
    h = -t0 / n
    # half-step phase advance: field phase e^{-i kappa^2 t / 2}; advancing
    # t by h/2 multiplies by e^{-i kappa^2 h/4} = e^{-i theta m^2},
    # theta = base^2 h / 4  (h < 0)
    Uh = phase_array(base * base * h / 4)
    P = phase_array(base * base * t0 / 2)   # e^{-i kappa^2 t0/2}

    IKm_pos = IKm[:NX // 2]
    IKm_neg = IKm[NX // 2:]

    def eval_f2(x, A):
        """f, f' at x from the stage coefficient array A_k = C_k P_k(t):
        two power ladders (positive / negative m) in fft order."""
        w = gmpy2.exp(mpc(mpfr(0), base * (x - xminm)))
        f = mpc(0)
        fp = mpc(0)
        wp = mpc(1, 0)
        for a_, ik in zip(A[:NX // 2], IKm_pos):
            term = a_ * wp
            f += term
            fp += ik * term
            wp = wp * w
        wn = (1 / w) ** (NX // 2)
        for a_, ik in zip(A[NX // 2:], IKm_neg):
            term = a_ * wn
            f += term
            fp += ik * term
            wn = wn * w
        return f, fp

    def vel(x, z, A):
        f, fp = eval_f2(x, A)
        fr, fi = f.real, f.imag
        pr, pim = fp.real, fp.imag
        rho = fr * fr + fi * fi
        aa = (fr * pim - fi * pr) / rho          # Im(conj(f) fp)/rho
        bb = (fr * pr + fi * pim) / rho          # Re(conj(f) fp)/rho
        ze = z
        if ze < mpfr(ZCLIP):
            ze = mpfr(ZCLIP)
        if ze > 1 - mpfr(ZCLIP):
            ze = 1 - mpfr(ZCLIP)
        vx = aa - nym * pi / gmpy2.tan(pi * ze)
        return vx, nym * bb

    def rk4_span(x, z, P0_, n_steps, hh, Uh_, store=None):
        """n_steps of classical RK4 with step hh; P0_ is the phase array
        at the start time.  Uh_ advances the phase by hh/2.  The stage
        coefficient arrays A_k = C_k P_k are folded once per stage time
        and shared (k2/k3 share the half step; k4's array is reused as
        the next step's k1).  Returns final (x, z, P)."""
        P0c = P0_
        A0 = [c * p for c, p in zip(Cm, P0c)]
        h2 = hh / 2
        h6 = hh / 6
        for st in range(n_steps):
            Ph = [p * u for p, u in zip(P0c, Uh_)]
            Pn = [p * u for p, u in zip(Ph, Uh_)]
            Ah = [c * p for c, p in zip(Cm, Ph)]
            An = [c * p for c, p in zip(Cm, Pn)]
            k1x, k1z = vel(x, z, A0)
            k2x, k2z = vel(x + h2 * k1x, z + h2 * k1z, Ah)
            k3x, k3z = vel(x + h2 * k2x, z + h2 * k2z, Ah)
            k4x, k4z = vel(x + hh * k3x, z + hh * k3z, An)
            x = x + h6 * (k1x + 2 * k2x + 2 * k3x + k4x)
            z = z + h6 * (k1z + 2 * k2z + 2 * k3z + k4z)
            P0c = Pn
            A0 = An
            if store is not None:
                store.append((x, z))
        return x, z, P0c

    # main pass, storing per-step states
    states = [(mpfr(D_NEAR), mpfr(z0f))]
    x, z, Pend = rk4_span(mpfr(D_NEAR), mpfr(z0f), P, n, h, Uh,
                          store=states)
    xs = [s[0] for s in states]
    im = int(max(range(len(xs)), key=lambda i: xs[i]))

    # refined margin: h/50 re-integration over [im-2, im+2] steps
    j0 = max(0, im - 2)
    span = min(len(xs) - 1, im + 2) - j0
    hh = h / 50
    Uh2 = phase_array(base * base * hh / 4)
    Pj = phase_array(base * base * (t0 + j0 * h) / 2)
    fine = [(states[j0][0], states[j0][1])]
    rk4_span(states[j0][0], states[j0][1], Pj, span * 50, hh, Uh2,
             store=fine)
    fx = [s[0] for s in fine]
    fm = int(max(range(len(fx)), key=lambda i: fx[i]))
    if 0 < fm < len(fx) - 1:
        y1, y2, y3 = fx[fm - 1], fx[fm], fx[fm + 1]
        den = 2 * y2 - y1 - y3
        vpk = y2 + (y1 - y3) * (y1 - y3) / (8 * den) if den != 0 else y2
    else:
        vpk = fx[fm]
    margin = vpk - mpfr(D_NEAR)

    # s_pc: first step crossing d + 1e-3, refined at h/50
    s_pc = None
    for j in range(1, len(xs)):
        if xs[j] >= d_class:
            Pj = phase_array(base * base * (t0 + (j - 1) * h) / 2)
            fine2 = [(states[j - 1][0], states[j - 1][1])]
            rk4_span(states[j - 1][0], states[j - 1][1], Pj, 50, hh, Uh2,
                     store=fine2)
            f2x = [s[0] for s in fine2]
            for jj in range(1, len(f2x)):
                if f2x[jj] >= d_class:
                    x_lo, x_hi = f2x[jj - 1], f2x[jj]
                    frac = (d_class - x_lo) / (x_hi - x_lo)
                    s_pc = t0 + (j - 1) * h + (jj - 1 + frac) * hh
                    break
            if s_pc is None:
                s_pc = t0 + j * h
            break

    return dict(x_end=float(x), z_end=float(z),
                x_end_str=format(x, ".36g"), z_end_str=format(z, ".36g"),
                margin=float(margin), margin_str=format(margin, ".36g"),
                s_pc=(float(s_pc) if s_pc is not None else np.nan),
                n_steps=n, h=float(h), elapsed=float(time.time() - tw))


def phase_ladder():
    from multiprocessing import Pool
    sel = np.load(os.path.join(OUTDIR, "_selection.npz"))
    lad = ladder_selection()
    np.save(os.path.join(OUTDIR, "_ladder_idx.npy"), lad)
    args = [(float(sel["ny"][j]), float(sel["t0"][j]), float(sel["z0"][j]),
             H_MP) for j in lad]
    # step-halving on the N_HALVE worst-margin ladder paths
    args_h = [(a[0], a[1], a[2], H_MP / 2.0) for a in args[:N_HALVE]]
    tw = time.time()
    with Pool(4) as pool:
        res = pool.map(mp_run_path, args + args_h)
    recs, recs_h = res[:len(args)], res[len(args):]
    json.dump(dict(ladder_idx=[int(j) for j in lad],
                   runs=recs, halved=recs_h,
                   elapsed=float(time.time() - tw)),
              open(os.path.join(OUTDIR, "_ladder.json"), "w"), indent=2)
    print("ladder done: %.0f s" % (time.time() - tw), flush=True)


# ----------------------------------------------------------------------
# mpmath cross-check of one field evaluation (backend equivalence record)
# ----------------------------------------------------------------------
def mpmath_crosscheck():
    import gmpy2
    from gmpy2 import mpfr, mpc
    import mpmath as mp
    gmpy2.get_context().precision = MP_BITS
    mp.mp.dps = 60
    xq, tq = 1.2345, 5.4321
    # gmpy2 route
    pi_g = gmpy2.const_pi()
    base_g = pi_g / 30
    theta = base_g * base_g * mpfr(tq) / 2
    mvals = list(range(0, NX // 2)) + list(range(-NX // 2, 0))
    fg = mpc(0)
    for k, m in enumerate(mvals):
        ph = gmpy2.exp(mpc(mpfr(0), -theta * m * m + base_g * m
                           * (mpfr(xq) - mpfr(XMIN))))
        fg += mpc(mpfr(CFT[k].real), mpfr(CFT[k].imag)) * ph
    # mpmath route
    base_m = mp.pi / 30
    fm = mp.mpc(0)
    for k, m in enumerate(mvals):
        ph = mp.expj(-base_m ** 2 * tq / 2 * m * m
                     + base_m * m * (xq - XMIN))
        fm += mp.mpc(CFT[k].real, CFT[k].imag) * ph
    diff = abs(complex(mp.mpc(str(fg.real), str(fg.imag)) - fm))
    return dict(x=xq, t=tq, gmpy2=str(fg)[:40], mpmath=str(fm)[:40],
                abs_diff=float(diff))


# ----------------------------------------------------------------------
# Phase GATES
# ----------------------------------------------------------------------
def phase_gates():
    g1 = json.load(open(os.path.join(OUTDIR, "_g1.json")))
    selinfo = json.load(open(os.path.join(OUTDIR, "_select.json")))
    sel = np.load(os.path.join(OUTDIR, "_selection.npz"))
    dop = np.load(os.path.join(OUTDIR, "_dop.npz"))
    lad_idx = np.load(os.path.join(OUTDIR, "_ladder_idx.npy"))
    lad = json.load(open(os.path.join(OUTDIR, "_ladder.json")))
    o1res = json.load(open(os.path.join(O1DIR, "o1_results.json")))

    stratum = sel["stratum"]
    cls_o1, cls_p2 = sel["cls_o1"], dop["cls"]
    mar_o1, mar_p2 = sel["margin_o1"], dop["margin"]
    spc_o1, spc_p2 = sel["s_pc_o1"], dop["s_pc"]

    # ---------------- G2: classification + margin agreement --------------
    keep = cls_o1 == cls_p2
    moved = np.where(~keep)[0]
    uncap = (cls_o1 == 0) & (mar_o1 < MARGIN_CAP)
    dmar = np.abs(mar_p2[uncap] - mar_o1[uncap])
    capped = (cls_o1 == 0) & (mar_o1 >= MARGIN_CAP)
    cap_ok = mar_p2[capped] >= MARGIN_CAP
    a_mask = cls_o1 == 0
    dspc = np.abs(spc_p2[a_mask] - spc_o1[a_mask])
    o1_p95 = o1res["G3"]["margin_err_p95"]
    o1_max = max(o1res["stability"]["nx2"]["margin_diff_uncapped"]["max"],
                 o1res["stability"]["dt2"]["margin_diff_uncapped"]["max"])
    thr_p95, thr_max = 2.0 * o1_p95, 2.0 * o1_max
    g2_ok = (moved.size == 0 and cap_ok.all()
             and float(np.percentile(dmar, 95)) <= thr_p95
             and float(dmar.max()) <= thr_max)
    g2 = dict(
        n_paths=int(cls_o1.size),
        n_kept=int(keep.sum()), n_moved=int(moved.size),
        moved_detail=[dict(i=int(i), o1=CLS_NAME[int(cls_o1[i])],
                           p2=CLS_NAME[int(cls_p2[i])]) for i in moved],
        margin_diff_uncapped=dict(
            n=int(uncap.sum()), median=float(np.median(dmar)),
            p95=float(np.percentile(dmar, 95)), max=float(dmar.max())),
        o1_error_estimate=dict(p95=o1_p95, max=o1_max),
        thresholds=dict(p95=thr_p95, max=thr_max),
        capped=dict(n=int(capped.sum()), all_recertified=bool(cap_ok.all()),
                    min_margin_p2=float(mar_p2[capped].min())
                    if capped.any() else None),
        s_pc_diff=dict(median=float(np.nanmedian(dspc)),
                       p95=float(np.nanpercentile(dspc, 95)),
                       max=float(np.nanmax(dspc))),
        monitors=dict(
            min_rho_rel_class_a=float(dop["min_rho_rel"][a_mask].min()),
            rho_floor_events=int(np.sum(~np.isnan(dop["s_floor"][a_mask]))),
            max_abs_a_class_a=float(dop["max_abs_a"][a_mask].max()),
            max_abs_b_class_a=float(dop["max_abs_b"][a_mask].max()),
            clamp_headroom_note="clips at +-500 never bind on class-(a)"
            if (dop["max_abs_a"][a_mask].max() < VCLAMP
                and dop["max_abs_b"][a_mask].max() < VCLAMP) else
            "CLAMP BOUND ON A CLASS-(a) PATH"),
        verdict="PASS" if g2_ok else "FAIL")

    # ---------------- G3: precision ladder -------------------------------
    runs = lad["runs"]
    ep_diff, m_rel, spc_diff = [], [], []
    for j, r in zip(lad_idx, runs):
        ep_diff.append(float(np.hypot(r["x_end"] - dop["x_end"][j],
                                      r["z_end"] - dop["z_end"][j])))
        m_rel.append(abs(r["margin"] - dop["margin"][j])
                     / abs(r["margin"]))
        spc_diff.append(abs(r["s_pc"] - dop["s_pc"][j]))
    halv = []
    for r0, rh in zip(runs[:N_HALVE], lad["halved"]):
        halv.append(dict(
            ep_shift=float(np.hypot(r0["x_end"] - rh["x_end"],
                                    r0["z_end"] - rh["z_end"])),
            margin_shift=abs(r0["margin"] - rh["margin"])))
    g3_ok = (max(ep_diff) <= 1e-8 and max(m_rel) <= 1e-6)
    g3 = dict(n_ladder=len(runs),
              endpoint_diff=dict(median=float(np.median(ep_diff)),
                                 max=float(np.max(ep_diff))),
              margin_rel_diff=dict(median=float(np.median(m_rel)),
                                   max=float(np.max(m_rel))),
              s_pc_abs_diff=dict(median=float(np.median(spc_diff)),
                                 max=float(np.max(spc_diff))),
              step_halving=halv,
              mp_method="fixed-step classical RK4, h = t0/ceil(t0/5e-4) "
                        "(~5e-4), 200-bit (~60-digit) GMP/MPFR (gmpy2, "
                        "mpmath's backend); margins by h/50 local "
                        "re-integration + 3-pt parabola; s_pc by h/50 "
                        "bracket + linear interpolation",
              mpmath_crosscheck=mpmath_crosscheck(),
              per_path=[dict(idx=int(j), ny=float(sel["ny"][j]),
                             t0=float(sel["t0"][j]), z0=float(sel["z0"][j]),
                             margin_f64=float(dop["margin"][j]),
                             margin_mp=r["margin"],
                             ep_diff=e, m_rel=m)
                        for j, r, e, m in zip(lad_idx, runs, ep_diff, m_rel)],
              verdict="PASS" if g3_ok else "FAIL")

    # ---------------- G4: the worst-margin path record -------------------
    a_idx = np.where(cls_p2 == 0)[0]
    iw = int(a_idx[np.argmin(mar_p2[a_idx])])
    in_ladder = iw in set(int(v) for v in lad_idx)
    lr = None
    if in_ladder:
        lr = runs[list(lad_idx).index(iw)]
    g4 = dict(
        selection_index=iw, o1_flat_index=int(sel["idx"][iw]),
        ny=float(sel["ny"][iw]), start_t=float(sel["t0"][iw]),
        start_z=float(sel["z0"][iw]), detector_x=D_NEAR,
        margin_o1=float(mar_o1[iw]), margin_f64=float(mar_p2[iw]),
        margin_mp=(lr["margin"] if lr else None),
        margin_mp_str=(lr["margin_str"] if lr else None),
        crossing_epoch_s_pc_f64=float(spc_p2[iw]),
        crossing_epoch_s_pc_mp=(lr["s_pc"] if lr else None),
        endpoint_f64=[float(dop["x_end"][iw]), float(dop["z_end"][iw])],
        endpoint_mp=([lr["x_end"], lr["z_end"]] if lr else None),
        endpoint_agreement=(float(np.hypot(lr["x_end"] - dop["x_end"][iw],
                                           lr["z_end"] - dop["z_end"][iw]))
                            if lr else None),
        margin_agreement_rel=(abs(lr["margin"] - float(mar_p2[iw]))
                              / lr["margin"] if lr else None),
        step_halving_margin_shift=(halv[0]["margin_shift"]
                                   if in_ladder and list(lad_idx).index(iw)
                                   < N_HALVE else None),
        min_rho_rel=float(dop["min_rho_rel"][iw]),
        in_ladder=bool(in_ladder))
    g4["verdict"] = ("PASS" if in_ladder and g4["endpoint_agreement"] <= 1e-8
                     and g4["margin_agreement_rel"] <= 1e-6 else "FAIL")

    # ---------------- G5: honest scope -----------------------------------
    g5 = dict(
        statement=(
            "PRECISION-CERTIFICATION, NOT FORMAL PROOF: the field is exact "
            "(a finite trigonometric sum evaluated to fp rounding; G1) and "
            "the integration is converged (adaptive rtol 1e-12 vs 50-digit "
            "fixed-step agreement; G3), but no step carries a rigorous "
            "enclosure.  Formal certification would require validated ODE "
            "integration with interval/Taylor-model enclosures (e.g. "
            "CAPD/COSY-class) of the backward flow over the exact field, "
            "with outward-rounded mode sums, so that each pre-crossing "
            "margin becomes a machine-checked inequality.  The analytic "
            "T4-W5 proof (no computation at all) remains the terminal "
            "open item."),
        verdict="PASS")

    g1_out = dict(spec="exact-vs-grid <= 1e-12 rel, >= 1000 nodes x "
                       "several times", **g1)

    results = dict(
        workstream="P2", file_id="F-T7-P2",
        params=dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                    box=[XMIN, XMAX], Nx=NX, d=D_NEAR,
                    kill_bins=KILL_BINS, ny_families=NY_FAMS,
                    strata=dict(wall=N_WALL, small=N_SMALL, rand=N_RAND),
                    seed_select=SEED_SELECT,
                    dop853=dict(rtol=RTOL, atol=ATOL),
                    ladder=dict(n=N_LADDER, h=H_MP, bits=MP_BITS,
                                n_halved=N_HALVE),
                    margin_class=MARGIN_CLASS, margin_cap=MARGIN_CAP,
                    rho_floor_rel=RHO_FLOOR_REL, vclamp=VCLAMP,
                    rho_eps_rel=RHO_EPS_REL),
        mirror_check=selinfo["mirror"],
        G1=g1_out, G2=g2, G3=g3, G4=g4, G5=g5,
        runtime=dict(dop_mean_nfev=float(dop["nfev"].mean()),
                     dop_total_s=float(dop["elapsed"].sum()),
                     ladder_total_s=lad["elapsed"]),
    )
    json.dump(results, open(os.path.join(OUTDIR, "p2_results.json"), "w"),
              indent=2)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4", "G5")}, flush=True)
    return results


# ----------------------------------------------------------------------
# Phase FIG
# ----------------------------------------------------------------------
def phase_fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK = "#0b0b0b"
    INK2 = "#52514e"
    SURF = "#fcfcfb"
    BLU = "#2a78d6"
    BLUD = "#12448c"
    ACC = "#e87ba4"
    GRN = "#3e8f6b"

    sel = np.load(os.path.join(OUTDIR, "_selection.npz"))
    dop = np.load(os.path.join(OUTDIR, "_dop.npz"))
    lad_idx = np.load(os.path.join(OUTDIR, "_ladder_idx.npy"))
    lad = json.load(open(os.path.join(OUTDIR, "_ladder.json")))
    res = json.load(open(os.path.join(OUTDIR, "p2_results.json")))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3), dpi=160)
    fig.patch.set_facecolor(SURF)
    axA, axB, axC = axes

    # A: margin distribution, O1 vs P2 exact (class-(a), uncapped)
    uncap = (sel["cls_o1"] == 0) & (sel["margin_o1"] < MARGIN_CAP)
    bins = np.linspace(0.6, 2.0, 43)
    axA.hist(sel["margin_o1"][uncap], bins=bins, color=BLU, alpha=0.55,
             label="O1 (grid-interp RK4)")
    axA.hist(dop["margin"][uncap], bins=bins, histtype="step", lw=1.6,
             color=BLUD, label="P2 (exact field, DOP853)")
    axA.axvline(res["G4"]["margin_f64"], color=ACC, lw=1.4, ls="--")
    axA.text(res["G4"]["margin_f64"] + 0.02, axA.get_ylim()[1] * 0.86,
             "worst margin\n%.6f" % res["G4"]["margin_f64"], fontsize=8,
             color=ACC)
    axA.set_xlabel("pre-crossing margin  max$_s$ X$_x$ $-$ d")
    axA.set_ylabel("paths")
    axA.set_title("A  margins: O1 vs exact-field re-run\n"
                  "(uncapped class-(a) subset, n=%d)" % int(uncap.sum()),
                  fontsize=10, loc="left", color=INK)
    axA.legend(fontsize=8, frameon=False)

    # B: margin agreement vs O1's own error estimate
    dm = np.abs(dop["margin"][uncap] - sel["margin_o1"][uncap])
    dm = np.maximum(dm, 1e-12)
    axB.hist(np.log10(dm), bins=36, color=BLU, alpha=0.75)
    p95_o1 = res["G2"]["o1_error_estimate"]["p95"]
    axB.axvline(np.log10(p95_o1), color=GRN, lw=1.6)
    axB.text(np.log10(p95_o1) + 0.06, axB.get_ylim()[1] * 0.82,
             "O1 refinement\np95 = %.1e" % p95_o1, fontsize=8, color=GRN)
    axB.axvline(np.log10(res["G2"]["margin_diff_uncapped"]["p95"]),
                color=ACC, lw=1.6, ls="--")
    axB.text(np.log10(res["G2"]["margin_diff_uncapped"]["p95"]) + 0.06,
             axB.get_ylim()[1] * 0.6,
             "measured p95\n= %.1e" % res["G2"]["margin_diff_uncapped"]["p95"],
             fontsize=8, color=ACC)
    axB.set_xlabel("log$_{10}$ |margin$_{P2}$ $-$ margin$_{O1}$|")
    axB.set_ylabel("paths")
    axB.set_title("B  P2$-$O1 margin difference vs O1's own\n"
                  "error estimate (G2)", fontsize=10, loc="left", color=INK)

    # C: precision ladder agreement
    ep = [p["ep_diff"] for p in res["G3"]["per_path"]]
    mr = [max(p["m_rel"], 1e-16) for p in res["G3"]["per_path"]]
    xs = np.arange(len(ep))
    axC.semilogy(xs, ep, "o", ms=5, color=BLUD,
                 label="endpoint |f64 $-$ mp50|")
    axC.semilogy(xs, mr, "s", ms=5, color=ACC, mfc="none",
                 label="margin rel. diff")
    for i, h in enumerate(res["G3"]["step_halving"]):
        axC.semilogy([i], [max(h["ep_shift"], 1e-16)], "^", ms=6,
                     color=GRN, label="mp h-halving shift" if i == 0 else None)
    axC.axhline(1e-8, color=INK2, lw=1.0, ls="--")
    axC.text(0.2, 1.6e-8, "G3 endpoint gate 1e-8", fontsize=8, color=INK2)
    axC.axhline(1e-6, color=INK2, lw=0.8, ls=":")
    axC.text(len(ep) * 0.55, 1.6e-6, "6-digit margin gate", fontsize=8,
             color=INK2)
    axC.set_xlabel("ladder path (0$-$11 = worst margins, then spread)")
    axC.set_ylabel("disagreement")
    axC.set_ylim(1e-16, 1e-4)
    axC.set_title("C  precision ladder: float64-DOP853 vs\n"
                  "200-bit RK4 on the exact field (G3)", fontsize=10,
                  loc="left", color=INK)
    axC.legend(fontsize=8, frameon=False, loc="center right")

    for ax in axes:
        ax.set_facecolor(SURF)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.suptitle("P2 — T4-W5 certification pilot: exact-field backward "
                 "paths (F-T7-P2)", fontsize=11, color=INK, x=0.01,
                 ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(os.path.join(OUTDIR, "p2_fig.png"),
                facecolor=SURF, bbox_inches="tight")
    print("figure written", flush=True)


# ----------------------------------------------------------------------
def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if phase in ("g1", "all"):
        print("== phase g1 ==", flush=True)
        phase_g1()
    if phase in ("select", "all"):
        print("== phase select ==", flush=True)
        phase_select()
    if phase in ("dop", "all"):
        print("== phase dop ==", flush=True)
        phase_dop()
    if phase in ("ladder", "all"):
        print("== phase ladder ==", flush=True)
        phase_ladder()
    if phase in ("gates", "all"):
        print("== phase gates ==", flush=True)
        phase_gates()
    if phase in ("fig", "all"):
        print("== phase fig ==", flush=True)
        phase_fig()
    print("total elapsed %.1f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
