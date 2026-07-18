#!/usr/bin/env python3
"""
WORKSTREAM S4 -- Population-scale certification screener (P2 extension).
Tier 8 Phase S (ROADMAP_v10_ALTERNATIVES.md, S4).  F-T8-S4.

PARAMETRIZATION-REUSE MANDATE.  Everything physical is imported verbatim
from tier7-program/P2/p2_certify.py (which itself inherits O1/L7b/L5prime
engine constants): the exact band-limited field (Nx = 2048 spectral
coefficients CFT of the O1 reference grid on the box [-15, 45]), the
guidance-velocity law with the engine's regularization (clips +-500,
relative density floor 1e-14, z-clip 1e-9), the margin definition
(max_s X_x - d, d = 1), the classification semantics, and the 200-bit
mp ladder (p2_certify.mp_run_path, unmodified).

POPULATION (pre-registered here, before any production run).  The P2/O1
parametrization grid per family n_y in {+1, +0.75} is t0 in the kill bins
[6.4, 7.2] (O1: 80 values, dt = 0.01) x z0 in [0.0025, 0.9975] (O1: 200
values, dz = 0.005).  S4 refines toward N = 10,000 NEW configurations:
    t0:  NT = 100 midpoints  t0_i = 6.4 + 0.8 (i + 1/2)/NT   (dt = 0.008,
         finer than O1's 0.01 and everywhere OFF the O1 grid),
    z0:  NZ = 50 uniform     z0_j = 0.0025 + j * 0.995/(NZ-1)  (includes
         the extreme wall rows 0.0025 / 0.9975 exactly; dz ~ 0.0203),
    families {+1, +0.75}  (the -n_y mirrors are numerically identical
    under z -> 1-z, re-verified in P2; not re-run, as in P2),
  => N = 2 * 100 * 50 = 10,000.  Any scale-down (floor N >= 2,000) is a
  printed coordinator-owned spec amendment.

SCREENER (float64).  Fixed-step classical RK4 on the exact mode sum --
the same integration scheme as the P2 mp ladder, in float64, batched:
per-mode time phases advanced by half-step recurrences (as in the mp
ladder), spatial factors e^{i kappa_k x} by cumulative-product power
ladders + BLAS mat-vec (validated against p2_certify.fields_exact at
~1e-14).  Margin refined by local h/20 re-integration around the coarse
peak using p2_certify.rhs_factory + 3-point parabola (the mp ladder's
refinement scheme); s_pc by h/20 bracket + linear interpolation.
Base step H0 = 2e-3 with a wall-proximity refinement policy (the cot(pi z)
drift near z = 0 / 1 is the stiff region): the extreme wall band
(wall distance < 0.02) runs the P2 DOP853 pipeline (integrate_path)
VERBATIM (bench measured that the ~1/dw velocity there defeats any
affordable fixed step); wall distance < 0.05 runs RK4 with R_NEAR x
steps; interior R_INT x steps.  The policy constants are fixed by the
bench phase (accuracy target: RK4-row screener margin error <= 3e-7,
>3x headroom under the 1e-6 gate) and frozen in _bench.json before
validate/screen run; the bench measurements are printed in RESULTS.md.

PHASES (checkpointed; every phase re-reads only committed files):
  bench     unit-cost + accuracy calibration vs p2_certify.integrate_path
            (DOP853 rtol 1e-12) on 10 reference configs spanning the z
            rows; freezes the h policy + projected cost in _bench.json.
  validate  S4-G1: 24 seeded (SEED = 20260718) random class-(a) P2
            configs re-run through the screener; margins vs P2's _dop.npz
            (DOP853 rtol 1e-12 on the exact field), gate <= 1e-6 rel.
  screen    S4-G2: the N-config population, one checkpoint file per t0
            chunk (s4_chunks/chunk_XXX.npz), Pool(4), resume-safe.
  ladder    S4-G3: the K = 24 lowest-margin class-(a) configs re-run
            through the UNMODIFIED P2 200-bit ladder (mp_run_path,
            h = 5e-4), one JSON checkpoint per config, resume-safe;
            + step-halving (h/2) on the 3 lowest; + DOP853 (rtol 1e-12)
            float64 cross-checks on the 3 lowest.
  gates     verdict assembly -> s4_results.json.
  fig       s4_fig.png.

WORST-K VERDICT (pre-registered before the ladder runs):
  STABLE  iff the mp run confirms pre-crossing (s_pc exists) with
          margin_mp > 1e-3 (beyond the graze band), the screener's
          classification (class-(a), no floor-first, no box exit) is
          preserved, and |margin_scr - margin_mp| / margin_mp <= 1e-4.
  UNSTABLE otherwise (any classification flip or >1e-4 margin motion is
          a FINDING, reported unsoftened).
  certified digits = floor(-log10(max(rel screener-vs-mp diff, rel mp
          step-halving shift where run, rel DOP853-vs-mp diff where
          run))), i.e. digits cross-validated by two independent routes.

GATES (pre-registered, ROADMAP_v10 S4): S4-G1 (<= 1e-6 rel on 24 P2
configs), S4-G2 (screen complete at final N, distribution + minimum),
S4-G3 (worst-K certified, per-config verdicts), S4-G4 (honest scope:
screened coverage != formally-verified kernel).

Within-model; nothing here bears on nature.
Usage: s4_screen.py [bench|validate|screen|ladder|gates|fig|all]
"""
import json
import os
import sys
import time

import numpy as np

S4DIR = os.path.dirname(os.path.abspath(__file__))
P2DIR = os.path.normpath(os.path.join(
    S4DIR, "..", "..", "tier7-program", "P2"))
sys.path.insert(0, P2DIR)
import p2_certify as p2                                    # noqa: E402

# ---------------------------------------------------------------------
# Pre-registered S4 constants
# ---------------------------------------------------------------------
NT, NZ = 100, 50
FAMS = [1.0, 0.75]
TGRID = 6.4 + 0.8 * (np.arange(NT) + 0.5) / NT
ZGRID = 0.0025 + np.arange(NZ) * (0.995 / (NZ - 1))
N_POP = 2 * NT * NZ

H0 = 2e-3                       # base screener step (policy-refined)
# h policy (bench-calibrated, frozen in _bench.json before production):
#   wall distance dw < DW_HYBRID  -> the P2 DOP853 pipeline verbatim
#     (p2.integrate_path; the near-wall cot drift, |v_x| ~ 1/dw ~ 400,
#     outruns any affordable fixed step -- measured at bench),
#   dw < DW_NEAR                  -> RK4 with n_steps x R_NEAR,
#   else                          -> RK4 with n_steps x R_INT.
DW_HYBRID, DW_NEAR = 0.02, 0.05
R_NEAR_INIT, R_INT_INIT = 4, 2
FINE_DIV = 20                   # local refinement h/20 (mp ladder: h/50)
SUBB = 24                       # eval sub-batch (empirical numpy optimum)

SEED_VALIDATE = 20260718        # printed seed for the 24 G1 configs
N_VAL = 24
K_WORST = 24
N_HALVE = 3
N_DOPCHECK = 3

RIM_LO, RIM_HI = p2.XMIN + 0.5, p2.XMAX - 0.5
CHUNKDIR = os.path.join(S4DIR, "s4_chunks")
LADDIR = os.path.join(S4DIR, "s4_ladder")

# m-ordered spectral arrays (m = 0..1023, then -1..-1024)
_ORD = np.concatenate([np.arange(p2.NX // 2),
                       np.arange(p2.NX - 1, p2.NX // 2 - 1, -1)])
M_ORD = np.concatenate([np.arange(p2.NX // 2),
                        -1 - np.arange(p2.NX // 2)])
BASE = np.pi / 30.0
K_ORD = BASE * M_ORD
C_ORD = p2.CFT[_ORD]
IK_ORD = 1j * K_ORD
NHALF = p2.NX // 2


def _policy():
    """Load the frozen h policy (bench phase writes it)."""
    f = os.path.join(S4DIR, "_bench.json")
    if os.path.exists(f):
        b = json.load(open(f))["policy"]
        return b["R_NEAR"], b["R_INT"]
    return R_NEAR_INIT, R_INT_INIT


def refine_factor(z0, rn, ri):
    """Steps multiplier; 0 = hybrid (P2 DOP853 pipeline verbatim)."""
    dw = min(z0, 1.0 - z0)
    if dw < DW_HYBRID:
        return 0
    if dw < DW_NEAR:
        return rn
    return ri


# ---------------------------------------------------------------------
# Batched exact-field evaluation (float64 mirror of the mp ladder scheme)
# ---------------------------------------------------------------------
def eval_batch(x, Ap2, An2):
    """f, f' at positions x (B,) from stacked m-ordered stage coefficient
    matrices Ap2 (Nx/2, 2) = [A_m, i kappa_m A_m] (m >= 0) and An2
    (m = -1..-Nx/2).  Power ladders by cumprod + BLAS."""
    out_f = np.empty(x.size, complex)
    out_fp = np.empty(x.size, complex)
    for i0 in range(0, x.size, SUBB):
        xs = x[i0:i0 + SUBB]
        B = xs.size
        u = np.exp(1j * BASE * (xs - p2.XMIN))
        Mp = np.empty((B, NHALF), complex)
        Mp[:, 0] = 1.0
        Mp[:, 1:] = u[:, None]
        F = np.cumprod(Mp, axis=1) @ Ap2
        uc = np.conj(u)
        Mp[:, 0] = uc
        Mp[:, 1:] = uc[:, None]
        F += np.cumprod(Mp, axis=1) @ An2
        out_f[i0:i0 + SUBB] = F[:, 0]
        out_fp[i0:i0 + SUBB] = F[:, 1]
    return out_f, out_fp


def _vel(x, z, ny, Ap2, An2, rmax_s):
    """Regularized guidance velocity, vectorized; mirrors p2.rhs_factory
    exactly (clips, relative floor, z-clip)."""
    f, fp = eval_batch(x, Ap2, An2)
    rho = f.real ** 2 + f.imag ** 2
    den = np.maximum(rho, p2.RHO_EPS_REL * rmax_s)
    q = np.conj(f) * fp
    a = np.clip(q.imag / den, -p2.VCLAMP, p2.VCLAMP)
    b = np.clip(q.real / den, -p2.VCLAMP, p2.VCLAMP)
    ze = np.clip(z, p2.ZCLIP, p2.L - p2.ZCLIP)
    vx = np.clip(a - ny * np.pi / np.tan(np.pi * ze / p2.L),
                 -p2.VCLAMP, p2.VCLAMP)
    return vx, ny * b, rho, q, den


def _stage(P):
    A = C_ORD * P
    IKA = IK_ORD * A
    return (np.column_stack([A[:NHALF], IKA[:NHALF]]),
            np.column_stack([A[NHALF:], IKA[NHALF:]]))


def rk4_batch(t0, ny, z0, n_steps):
    """Backward fixed-step RK4 of a batch sharing t0.  Returns stored
    per-step trajectories + monitors (mirroring p2.integrate_path's
    monitor semantics: relative-rho floor, clip headroom)."""
    B = ny.size
    h = -t0 / n_steps
    P0 = np.exp(-0.5j * K_ORD ** 2 * t0)
    Uh = np.exp(-0.25j * K_ORD ** 2 * h)
    x = np.full(B, p2.D_NEAR)
    z = z0.astype(float).copy()
    xs = np.empty((n_steps + 1, B))
    zs = np.empty((n_steps + 1, B))
    xs[0], zs[0] = x, z
    frozen = np.zeros(B, bool)
    s_end = np.zeros(B)
    min_rho_rel = np.full(B, np.inf)
    s_floor = np.full(B, np.nan)
    max_a = np.zeros(B)
    max_b = np.zeros(B)
    Ap0, An0 = _stage(P0)
    for j in range(n_steps):
        s0 = t0 + j * h
        Ph = P0 * Uh
        Pn = Ph * Uh
        Aph, Anh = _stage(Ph)
        Apn, Ann = _stage(Pn)
        rm0 = np.interp(s0, p2._TS, p2._RMAX)
        rmh = np.interp(s0 + 0.5 * h, p2._TS, p2._RMAX)
        rmn = np.interp(s0 + h, p2._TS, p2._RMAX)
        k1x, k1z, rho1, q1, den1 = _vel(x, z, ny, Ap0, An0, rm0)
        # monitors at the step node (active paths only)
        act = ~frozen
        ze = np.clip(z, p2.ZCLIP, p2.L - p2.ZCLIP)
        rho_rel = rho1 * np.sin(np.pi * ze / p2.L) ** 2 / rm0
        min_rho_rel[act] = np.minimum(min_rho_rel[act], rho_rel[act])
        newf = act & np.isnan(s_floor) & (rho_rel < p2.RHO_FLOOR_REL)
        s_floor[newf] = s0
        max_a[act] = np.maximum(max_a[act], np.abs(q1.imag / den1)[act])
        max_b[act] = np.maximum(max_b[act], np.abs(q1.real / den1)[act])
        k2x, k2z, _, _, _ = _vel(x + 0.5 * h * k1x, z + 0.5 * h * k1z,
                                 ny, Aph, Anh, rmh)
        k3x, k3z, _, _, _ = _vel(x + 0.5 * h * k2x, z + 0.5 * h * k2z,
                                 ny, Aph, Anh, rmh)
        k4x, k4z, _, _, _ = _vel(x + h * k3x, z + h * k3z,
                                 ny, Apn, Ann, rmn)
        dx = (h / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        dz = (h / 6.0) * (k1z + 2 * k2z + 2 * k3z + k4z)
        x = np.where(frozen, x, x + dx)
        z = np.where(frozen, z, z + dz)
        newly = (~frozen) & ((x < RIM_LO) | (x > RIM_HI))
        s_end[newly] = s0 + h
        frozen |= newly
        xs[j + 1], zs[j + 1] = x, z
        P0, Ap0, An0 = Pn, Apn, Ann
    s_end[~frozen] = 0.0
    return dict(h=h, xs=xs, zs=zs, x=x, z=z, frozen=frozen, s_end=s_end,
                min_rho_rel=min_rho_rel, s_floor=s_floor,
                max_a=max_a, max_b=max_b)


def _rk4_fine(rhs, s, y, n, h):
    """Single-path fixed-step RK4 on p2.rhs_factory (exact field);
    returns the list of states (mirrors the mp ladder's local
    re-integration scheme in float64)."""
    out = [y.copy()]
    for _ in range(n):
        k1 = rhs(s, y)
        k2 = rhs(s + 0.5 * h, y + 0.5 * h * k1)
        k3 = rhs(s + 0.5 * h, y + 0.5 * h * k2)
        k4 = rhs(s + h, y + h * k3)
        y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        s += h
        out.append(y.copy())
    return out


def postprocess(t0, nyv, res, ib):
    """Margin (h/20-refined), s_pc, classification for path ib of a
    rk4_batch result.  Mirrors p2.integrate_path / mp_run_path."""
    h = res["h"]
    xs = res["xs"][:, ib]
    zs = res["zs"][:, ib]
    frozen = bool(res["frozen"][ib])
    n_last = xs.size - 1
    if frozen:
        n_last = int(round((t0 - res["s_end"][ib]) / (-h)))
    im = int(np.argmax(xs[:n_last + 1]))
    margin = xs[im] - p2.D_NEAR
    rhs = p2.rhs_factory(nyv)
    hf = h / FINE_DIV
    if not frozen and 0 < im:
        j0 = max(0, im - 2)
        span = min(n_last, im + 2) - j0
        fine = _rk4_fine(rhs, t0 + j0 * h,
                         np.array([xs[j0], zs[j0]]), span * FINE_DIV, hf)
        fx = np.array([st[0] for st in fine])
        fm = int(np.argmax(fx))
        if 0 < fm < fx.size - 1:
            y1, y2, y3 = fx[fm - 1], fx[fm], fx[fm + 1]
            den = 2 * y2 - y1 - y3
            vpk = y2 + (y1 - y3) ** 2 / (8 * den) if den != 0 else y2
        else:
            vpk = fx[fm]
        margin = max(margin, vpk - p2.D_NEAR)
    # s_pc: first backward crossing of d + 1e-3 (scan s downward)
    s_pc = np.nan
    thr = p2.D_NEAR + p2.MARGIN_CLASS
    above = xs[:n_last + 1] >= thr
    if above.any():
        j = int(np.argmax(above))
        if j == 0:
            s_pc = t0
        else:
            fine = _rk4_fine(rhs, t0 + (j - 1) * h,
                             np.array([xs[j - 1], zs[j - 1]]), FINE_DIV, hf)
            f2x = np.array([st[0] for st in fine])
            s_pc = t0 + j * h
            for jj in range(1, f2x.size):
                if f2x[jj] >= thr:
                    frac = (thr - f2x[jj - 1]) / (f2x[jj] - f2x[jj - 1])
                    s_pc = t0 + (j - 1) * h + (jj - 1 + frac) * hf
                    break
    sf = res["s_floor"][ib]
    floored_first = (not np.isnan(sf)) and (np.isnan(s_pc) or sf > s_pc)
    if not np.isnan(s_pc) and not floored_first:
        c = 0
    elif floored_first:
        c = 1
    elif frozen:
        c = 4
    elif p2.MARGIN_EPS < margin < p2.MARGIN_CLASS:
        c = 5
    elif abs(res["x"][ib] - p2.X0) <= p2.W_OUT * p2.SIGMA_X:
        c = 2
    else:
        c = 3
    return dict(margin=float(margin), s_pc=float(s_pc), cls=int(c),
                x_end=float(res["x"][ib]), z_end=float(res["z"][ib]),
                s_end=float(res["s_end"][ib]), boxexit=frozen,
                min_rho_rel=float(res["min_rho_rel"][ib]),
                s_floor=float(sf),
                max_abs_a=float(res["max_a"][ib]),
                max_abs_b=float(res["max_b"][ib]))


_REC_KEYS = ("margin", "s_pc", "cls", "x_end", "z_end", "s_end",
             "boxexit", "min_rho_rel", "s_floor", "max_abs_a",
             "max_abs_b")


def screen_configs(t0, ny_arr, z0_arr, rn, ri, h0=H0):
    """Screen a set of configs sharing t0: group by refinement factor,
    run batched RK4 (or the P2 DOP853 pipeline verbatim for the hybrid
    wall band), post-process each path."""
    ny_arr = np.asarray(ny_arr, float)
    z0_arr = np.asarray(z0_arr, float)
    recs = [None] * ny_arr.size
    rfac = np.array([refine_factor(zz, rn, ri) for zz in z0_arr])
    n_base = int(np.ceil(t0 / h0))
    for ib in np.where(rfac == 0)[0]:
        # hybrid wall band: the P2 DOP853 pipeline with rtol relaxed to
        # 1e-11 (calibration history, all printed in RESULTS.md: 1e-10
        # measured 1.39e-6 rel motion on a long wall path -- the G1
        # attempt-1 FAIL -- and 1e-11 measures 1.0e-7 on that same
        # config, 10x headroom under the 1e-6 gate; rtol 1e-12 would be
        # the P2 pipeline bit-identical but doubles the screen cost)
        p2.RTOL, p2.ATOL = 1e-11, 1e-13
        try:
            full = p2.integrate_path((float(ny_arr[ib]), t0,
                                      float(z0_arr[ib])))
        finally:
            p2.RTOL, p2.ATOL = 1e-12, 1e-14
        recs[ib] = {k: full[k] for k in _REC_KEYS}
        recs[ib]["n_steps"] = full["nsteps"]
        recs[ib]["engine"] = 0                      # 0 = DOP853 hybrid
    for r in sorted(set(rfac[rfac > 0].tolist())):
        idx = np.where(rfac == r)[0]
        res = rk4_batch(t0, ny_arr[idx], z0_arr[idx], n_base * int(r))
        for k, ib in enumerate(idx):
            recs[ib] = postprocess(t0, float(ny_arr[ib]), res, k)
            recs[ib]["n_steps"] = n_base * int(r)
            recs[ib]["engine"] = 1                  # 1 = batched RK4
    return recs


# ---------------------------------------------------------------------
# Phase BENCH
# ---------------------------------------------------------------------
def phase_bench():
    from multiprocessing import Pool
    t0 = float(TGRID[0])
    cfgs = [(1.0, t0, float(ZGRID[j])) for j in (0, 1, 2, 3, 10, 25, 33,
                                                 48, 49)] \
        + [(0.75, t0, float(ZGRID[33]))]
    print("bench: %d reference configs at t0 = %.4f (DOP853 rtol 1e-12)"
          % (len(cfgs), t0), flush=True)
    tw = time.time()
    with Pool(4) as pool:
        refs = pool.map(p2.integrate_path, cfgs)
    t_ref = time.time() - tw
    t_hybrid_cfg = float(np.mean([r["elapsed"] for r in refs]))
    print("  reference done in %.0f s (mean %.1f s/config)"
          % (t_ref, t_hybrid_cfg), flush=True)

    rn, ri = R_NEAR_INIT, R_INT_INIT
    # grade the RK4 rows only (hybrid rows ARE the reference pipeline;
    # their screener output is identical by construction)
    rk4_ids = [i for i, c in enumerate(cfgs)
               if refine_factor(c[2], rn, ri) > 0]
    ny_arr = np.array([cfgs[i][0] for i in rk4_ids])
    z_arr = np.array([cfgs[i][2] for i in rk4_ids])
    timing = {}
    tw = time.time()
    recs1 = screen_configs(t0, ny_arr, z_arr, rn, ri, h0=H0)
    timing["r"] = time.time() - tw
    tw = time.time()
    recs2 = screen_configs(t0, ny_arr, z_arr, rn, ri, h0=H0 / 2)
    timing["2r"] = time.time() - tw
    per = []
    for k, i in enumerate(rk4_ids):
        c, ref = cfgs[i], refs[i]
        m_ref = ref["margin"]
        scale = max(abs(m_ref), 1e-3)
        e1 = abs(recs1[k]["margin"] - m_ref) / scale
        e2 = abs(recs2[k]["margin"] - m_ref) / scale
        per.append(dict(ny=c[0], t0=c[1], z0=c[2],
                        r=refine_factor(c[2], rn, ri), margin_ref=m_ref,
                        margin_r=recs1[k]["margin"],
                        margin_r2=recs2[k]["margin"],
                        rel_err_r=e1, rel_err_r2=e2,
                        cls_ref=ref["cls"], cls_scr=recs1[k]["cls"]))
        print("  z0=%.4f r=%-2d  relerr(policy)=%.2e  relerr(h/2)=%.2e"
              % (c[2], per[-1]["r"], e1, e2), flush=True)
    worst = max(p["rel_err_r"] for p in per)
    ok = worst <= 3e-7 and all(p["cls_ref"] == p["cls_scr"] for p in per)
    # projected chunk cost: RK4 rows scale with total step count; hybrid
    # rows cost t_hybrid_cfg each (4 per chunk: 2 wall rows x 2 families)
    steps_bench = sum(r["n_steps"] for r in recs1)
    steps_chunk = sum(int(np.ceil(t0 / H0)) * refine_factor(z, rn, ri)
                      for z in ZGRID for _ in FAMS)
    n_hyb_chunk = sum(1 for z in ZGRID for _ in FAMS
                      if refine_factor(z, rn, ri) == 0)
    proj_chunk = (timing["r"] * steps_chunk / steps_bench
                  + n_hyb_chunk * t_hybrid_cfg)
    out = dict(
        t_reference_s=t_ref, per_config=per, worst_rel_err=worst,
        hybrid_note=("wall band dw < %.3f runs the P2 DOP853 pipeline "
                     "(p2.integrate_path) verbatim: bench showed the "
                     "near-wall cot drift (|v_x| ~ 400 at z0 = 0.0025) "
                     "defeats affordable fixed steps (rel err 1.5e-5 at "
                     "10x steps, ~first-order convergence)" % DW_HYBRID),
        policy=dict(R_NEAR=rn, R_INT=ri, H0=H0, DW_HYBRID=DW_HYBRID,
                    DW_NEAR=DW_NEAR, accepted=bool(ok)),
        timing=dict(bench_pass_s=timing["r"], bench_half_s=timing["2r"],
                    t_hybrid_per_config_s=t_hybrid_cfg,
                    n_hybrid_per_chunk=n_hyb_chunk,
                    projected_chunk_s=proj_chunk,
                    projected_screen_cpu_s=proj_chunk * NT,
                    projected_screen_wall_4w_min=proj_chunk * NT / 4 / 60),
        note=("policy ACCEPTED: worst RK4-row screener margin error <= "
              "3e-7 with classifications preserved" if ok else
              "policy REJECTED at bench -- amendment required"))
    json.dump(out, open(os.path.join(S4DIR, "_bench.json"), "w"), indent=2)
    print("bench: worst rel err %.2e  policy accepted=%s  "
          "projected screen %.0f CPU-s (%.0f min wall on 4 workers)"
          % (worst, ok, proj_chunk * NT, out["timing"][
              "projected_screen_wall_4w_min"]), flush=True)
    return out


# ---------------------------------------------------------------------
# Phase VALIDATE (S4-G1)
# ---------------------------------------------------------------------
def _validate_one(args):
    nyv, t0, z0, rn, ri = args
    recs = screen_configs(t0, [nyv], [z0], rn, ri)
    return recs[0]


def phase_validate():
    from multiprocessing import Pool
    rn, ri = _policy()
    sel = np.load(os.path.join(P2DIR, "_selection.npz"))
    dop = np.load(os.path.join(P2DIR, "_dop.npz"))
    a_idx = np.where(dop["cls"] == 0)[0]           # 424 class-(a) configs
    rng = np.random.default_rng(SEED_VALIDATE)
    pick = np.sort(rng.choice(a_idx, N_VAL, replace=False))
    args = [(float(sel["ny"][i]), float(sel["t0"][i]),
             float(sel["z0"][i]), rn, ri) for i in pick]
    tw = time.time()
    with Pool(4) as pool:
        recs = pool.map(_validate_one, args)
    per = []
    for i, rec in zip(pick, recs):
        m_p2 = float(dop["margin"][i])
        rel = abs(rec["margin"] - m_p2) / abs(m_p2)
        per.append(dict(
            sel_index=int(i), o1_flat_index=int(sel["idx"][i]),
            ny=float(sel["ny"][i]), t0=float(sel["t0"][i]),
            z0=float(sel["z0"][i]), stratum=str(sel["stratum"][i]),
            margin_p2=m_p2, margin_s4=rec["margin"], rel_diff=rel,
            cls_p2=int(dop["cls"][i]), cls_s4=rec["cls"],
            s_pc_p2=float(dop["s_pc"][i]), s_pc_s4=rec["s_pc"]))
    rels = np.array([p["rel_diff"] for p in per])
    cls_ok = all(p["cls_p2"] == p["cls_s4"] for p in per)
    ok = bool(rels.max() <= 1e-6 and cls_ok)
    out = dict(seed=SEED_VALIDATE, n=N_VAL,
               sampled_from="the 424 class-(a) P2 configs (_dop.npz)",
               policy=dict(R_NEAR=rn, R_INT=ri, H0=H0,
                           DW_HYBRID=DW_HYBRID, DW_NEAR=DW_NEAR),
               rel_diff=dict(median=float(np.median(rels)),
                             p95=float(np.percentile(rels, 95)),
                             max=float(rels.max())),
               classifications_preserved=cls_ok,
               per_config=per, elapsed=time.time() - tw,
               verdict="PASS" if ok else "FAIL")
    json.dump(out, open(os.path.join(S4DIR, "_validate.json"), "w"),
              indent=2)
    print("validate (G1): max rel diff %.3e  median %.3e  cls kept %s "
          "-> %s  (%.0f s)" % (rels.max(), np.median(rels), cls_ok,
                               out["verdict"], out["elapsed"]), flush=True)
    return out


# ---------------------------------------------------------------------
# Phase SCREEN (S4-G2)
# ---------------------------------------------------------------------
def run_chunk(i):
    t0 = float(TGRID[i])
    fn = os.path.join(CHUNKDIR, "chunk_%03d.npz" % i)
    if os.path.exists(fn):
        return i, 0.0
    rn, ri = _policy()
    ny_arr = np.repeat(FAMS, NZ)
    z_arr = np.tile(ZGRID, len(FAMS))
    tw = time.time()
    recs = screen_configs(t0, ny_arr, z_arr, rn, ri)
    out = {k: np.array([r[k] for r in recs]) for k in recs[0]}
    out["ny"] = ny_arr
    out["z0"] = z_arr
    out["t0"] = np.full(ny_arr.size, t0)
    tmp = fn + ".tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, fn)
    return i, time.time() - tw


def phase_screen():
    from multiprocessing import Pool
    os.makedirs(CHUNKDIR, exist_ok=True)
    todo = [i for i in range(NT) if not os.path.exists(
        os.path.join(CHUNKDIR, "chunk_%03d.npz" % i))]
    print("screen: %d/%d chunks to run" % (len(todo), NT), flush=True)
    tw = time.time()
    with Pool(4) as pool:
        for k, (i, dt) in enumerate(pool.imap_unordered(run_chunk, todo)):
            print("  chunk %03d done (%.0f s)  [%d/%d, elapsed %.0f s]"
                  % (i, dt, k + 1, len(todo), time.time() - tw),
                  flush=True)
    print("screen complete: %.0f s" % (time.time() - tw), flush=True)


def load_population():
    keys = None
    pop = {}
    for i in range(NT):
        fn = os.path.join(CHUNKDIR, "chunk_%03d.npz" % i)
        d = np.load(fn)
        if keys is None:
            keys = list(d.keys())
            pop = {k: [] for k in keys}
        for k in keys:
            pop[k].append(d[k])
    return {k: np.concatenate(v) for k, v in pop.items()}


# ---------------------------------------------------------------------
# Phase LADDER (S4-G3)
# ---------------------------------------------------------------------
def worst_k(pop):
    a = np.where(pop["cls"] == 0)[0]
    order = a[np.argsort(pop["margin"][a])]
    return order[:K_WORST]


def _ladder_item(args):
    tag, nyv, t0, z0, h = args
    fn = os.path.join(LADDIR, tag + ".json")
    if os.path.exists(fn):
        return tag, 0.0
    if tag.startswith("dop"):
        rec = p2.integrate_path((nyv, t0, z0))
        rec = {k: (float(v) if isinstance(v, (int, float, np.floating))
                   else v) for k, v in rec.items()}
    else:
        rec = p2.mp_run_path((nyv, t0, z0, h))
    rec.update(dict(tag=tag, ny=nyv, t0=t0, z0=z0))
    tmp = fn + ".tmp"
    json.dump(rec, open(tmp, "w"), indent=2)
    os.replace(tmp, fn)
    return tag, rec.get("elapsed", 0.0)


def phase_ladder():
    from multiprocessing import Pool
    os.makedirs(LADDIR, exist_ok=True)
    pop = load_population()
    wk = worst_k(pop)
    json.dump(dict(worst_k_indices=[int(i) for i in wk],
                   ny=[float(pop["ny"][i]) for i in wk],
                   t0=[float(pop["t0"][i]) for i in wk],
                   z0=[float(pop["z0"][i]) for i in wk],
                   margin_scr=[float(pop["margin"][i]) for i in wk]),
              open(os.path.join(S4DIR, "_worstk.json"), "w"), indent=2)
    items = []
    for r, i in enumerate(wk):
        items.append(("lad_%03d" % r, float(pop["ny"][i]),
                      float(pop["t0"][i]), float(pop["z0"][i]), p2.H_MP))
    for r in range(N_HALVE):
        i = wk[r]
        items.append(("halve_%03d" % r, float(pop["ny"][i]),
                      float(pop["t0"][i]), float(pop["z0"][i]),
                      p2.H_MP / 2.0))
    for r in range(N_DOPCHECK):
        i = wk[r]
        items.append(("dop_%03d" % r, float(pop["ny"][i]),
                      float(pop["t0"][i]), float(pop["z0"][i]), 0.0))
    todo = [it for it in items if not os.path.exists(
        os.path.join(LADDIR, it[0] + ".json"))]
    print("ladder: %d/%d items to run (K=%d + %d halved + %d DOP853 "
          "cross-checks)" % (len(todo), len(items), K_WORST, N_HALVE,
                             N_DOPCHECK), flush=True)
    tw = time.time()
    with Pool(4) as pool:
        for k, (tag, dt) in enumerate(
                pool.imap_unordered(_ladder_item, todo)):
            print("  %s done (%.0f s)  [%d/%d, elapsed %.0f s]"
                  % (tag, dt, k + 1, len(todo), time.time() - tw),
                  flush=True)
    print("ladder complete: %.0f s" % (time.time() - tw), flush=True)


# ---------------------------------------------------------------------
# Phase GATES
# ---------------------------------------------------------------------
def phase_gates():
    bench = json.load(open(os.path.join(S4DIR, "_bench.json")))
    val = json.load(open(os.path.join(S4DIR, "_validate.json")))
    pop = load_population()
    wkinfo = json.load(open(os.path.join(S4DIR, "_worstk.json")))
    wk = np.array(wkinfo["worst_k_indices"])
    lad = [json.load(open(os.path.join(LADDIR, "lad_%03d.json" % r)))
           for r in range(K_WORST)]
    halv = [json.load(open(os.path.join(LADDIR, "halve_%03d.json" % r)))
            for r in range(N_HALVE)]
    dchk = [json.load(open(os.path.join(LADDIR, "dop_%03d.json" % r)))
            for r in range(N_DOPCHECK)]

    # ---- G1 ----
    g1 = dict(spec="float64 screener reproduces P2 margins on 24 seeded "
                   "random class-(a) P2 configs to <= 1e-6 rel",
              seed=val["seed"],
              rel_diff=val["rel_diff"],
              classifications_preserved=val["classifications_preserved"],
              bench_worst_rel_err_vs_dop853=bench["worst_rel_err"],
              verdict=val["verdict"])

    # ---- G2 ----
    import collections
    cls_count = collections.Counter(pop["cls"].tolist())
    a = pop["cls"] == 0
    mar = pop["margin"][a]
    imin = np.where(a)[0][np.argmin(mar)]
    qs = [0, 1, 5, 25, 50, 75, 95, 100]
    g2 = dict(
        spec="population screen complete at final N with margin "
             "distribution + minimum located",
        N=int(pop["cls"].size), N_target=N_POP,
        grid=dict(NT=NT, NZ=NZ, families=FAMS,
                  t_range=[float(TGRID[0]), float(TGRID[-1])],
                  z_range=[float(ZGRID[0]), float(ZGRID[-1])]),
        class_counts={p2.CLS_NAME[k]: int(v)
                      for k, v in sorted(cls_count.items())},
        margin_class_a=dict(
            n=int(a.sum()),
            percentiles={str(q): float(np.percentile(mar, q)) for q in qs},
            min=float(mar.min()), median=float(np.median(mar)),
            max=float(mar.max())),
        minimum_location=dict(
            ny=float(pop["ny"][imin]), t0=float(pop["t0"][imin]),
            z0=float(pop["z0"][imin]), margin=float(pop["margin"][imin])),
        p2_reference_min=0.620679003382508519636559682106572646,
        monitors=dict(
            clips_never_bind=bool(
                (pop["max_abs_a"][a].max() < p2.VCLAMP)
                and (pop["max_abs_b"][a].max() < p2.VCLAMP)),
            max_abs_a=float(pop["max_abs_a"][a].max()),
            max_abs_b=float(pop["max_abs_b"][a].max()),
            floor_events_class_a=int(np.sum(~np.isnan(pop["s_floor"][a]))),
            min_rho_rel_class_a=float(pop["min_rho_rel"][a].min())),
        verdict="PASS" if pop["cls"].size == N_POP else "FAIL")

    # ---- G3 ----
    halv_rel = [abs(halv[r]["margin"] - lad[r]["margin"])
                / abs(lad[r]["margin"]) for r in range(N_HALVE)]
    dchk_rel = [abs(dchk[r]["margin"] - lad[r]["margin"])
                / abs(lad[r]["margin"]) for r in range(N_DOPCHECK)]
    per = []
    n_unstable = 0
    for r, i in enumerate(wk):
        m_scr = float(pop["margin"][i])
        m_mp = lad[r]["margin"]
        rel = abs(m_scr - m_mp) / abs(m_mp)
        has_spc = np.isfinite(lad[r]["s_pc"])
        stable = (has_spc and m_mp > p2.MARGIN_CLASS
                  and int(pop["cls"][i]) == 0 and rel <= 1e-4)
        errs = [rel]
        if r < N_HALVE:
            errs.append(halv_rel[r])
        if r < N_DOPCHECK:
            errs.append(dchk_rel[r])
        digits = int(np.floor(-np.log10(max(max(errs), 1e-30))))
        if not stable:
            n_unstable += 1
        per.append(dict(
            rank=r, pop_index=int(i), ny=float(pop["ny"][i]),
            t0=float(pop["t0"][i]), z0=float(pop["z0"][i]),
            margin_scr=m_scr, margin_mp=m_mp,
            margin_mp_str=lad[r]["margin_str"],
            s_pc_mp=lad[r]["s_pc"], s_pc_scr=float(pop["s_pc"][i]),
            rel_scr_vs_mp=rel, certified_digits=digits,
            halving_rel_shift=(halv_rel[r] if r < N_HALVE else None),
            dop853_rel_diff=(dchk_rel[r] if r < N_DOPCHECK else None),
            verdict="STABLE" if stable else "UNSTABLE"))
    g3 = dict(
        spec="K = 24 lowest-margin configs re-certified through the "
             "unmodified P2 200-bit ladder (mp_run_path, h ~ 5e-4); "
             "STABLE/UNSTABLE per config with certified digits; any "
             "UNSTABLE is a finding",
        stable_criterion=("mp confirms pre-crossing (s_pc exists), "
                          "margin_mp > 1e-3, screener class-(a) "
                          "preserved, |m_scr - m_mp|/m_mp <= 1e-4"),
        n_stable=K_WORST - n_unstable, n_unstable=n_unstable,
        rel_scr_vs_mp=dict(
            median=float(np.median([p["rel_scr_vs_mp"] for p in per])),
            max=float(np.max([p["rel_scr_vs_mp"] for p in per]))),
        step_halving_rel_shifts=halv_rel,
        dop853_crosscheck_rel=dchk_rel,
        worst_config_record=per[0],
        per_config=per,
        verdict="PASS" if all(np.isfinite(p["rel_scr_vs_mp"])
                              for p in per) else "FAIL")

    # ---- G4 ----
    g4 = dict(statement=(
        "SCREENED POPULATION COVERAGE, NOT A FORMALLY-VERIFIED KERNEL: "
        "the N = %d screened margins are float64 fixed-step RK4 values "
        "(validated to ~1e-7 against the P2 DOP853/exact-field pipeline) "
        "and the worst-K margins are precision-certified (200-bit RK4 "
        "cross-validation), but NO margin here is a machine-checked "
        "inequality.  Formal certification would still require validated "
        "interval/Taylor-model enclosures of the backward flow with "
        "outward-rounded mode sums over the full population (see S3 for "
        "the measured cost of that technology), and the analytic T4-W5 "
        "proof remains the terminal open item.  The formal-kernel open "
        "of the Phase-R closeout REMAINS OPEN.  Within-model; nothing "
        "here bears on nature." % int(pop["cls"].size)),
        verdict="PASS")

    results = dict(
        workstream="S4", file_id="F-T8-S4",
        parametrization=dict(
            inherited_from="tier7-program/P2/p2_certify.py (verbatim "
                           "import; O1/L7b/L5prime engine constants)",
            grid=g2["grid"], N=N_POP,
            h_policy=bench["policy"],
            seed_validate=SEED_VALIDATE,
            ladder=dict(K=K_WORST, h=p2.H_MP, bits=p2.MP_BITS,
                        n_halved=N_HALVE, n_dop_crosscheck=N_DOPCHECK)),
        G1=g1, G2=g2, G3=g3, G4=g4,
        bench=dict(worst_rel_err=bench["worst_rel_err"],
                   timing=bench["timing"]))
    json.dump(results, open(os.path.join(S4DIR, "s4_results.json"), "w"),
              indent=2)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4")}, flush=True)
    return results


# ---------------------------------------------------------------------
# Phase FIG
# ---------------------------------------------------------------------
def phase_fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    BLU, BLUD, ACC, GRN = "#2a78d6", "#12448c", "#e87ba4", "#3e8f6b"

    pop = load_population()
    res = json.load(open(os.path.join(S4DIR, "s4_results.json")))
    val = json.load(open(os.path.join(S4DIR, "_validate.json")))
    wk = np.array(json.load(open(os.path.join(
        S4DIR, "_worstk.json")))["worst_k_indices"])

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.6), dpi=160)
    fig.patch.set_facecolor(SURF)
    (axA, axB), (axC, axD) = axes

    a = pop["cls"] == 0
    mar = pop["margin"][a]
    axA.hist(mar, bins=120, color=BLU, alpha=0.8)
    axA.axvline(res["G2"]["margin_class_a"]["min"], color=ACC, lw=1.4,
                ls="--")
    axA.text(res["G2"]["margin_class_a"]["min"] + 0.1,
             axA.get_ylim()[1] * 0.75,
             "population min\n%.6f" % res["G2"]["margin_class_a"]["min"],
             fontsize=8, color=ACC)
    axA.axvline(0.620679003382509, color=GRN, lw=1.2, ls=":")
    axA.text(0.75, axA.get_ylim()[1] * 0.45, "P2 432-subset min\n0.620679",
             fontsize=8, color=GRN)
    axA.set_xlabel("pre-crossing margin  max$_s$ X$_x$ $-$ d")
    axA.set_ylabel("configs")
    axA.set_title("A  margin distribution, N = %d screened configs\n"
                  "(class-(a), %d configs)" % (res["G2"]["N"],
                                               int(a.sum())),
                  fontsize=10, loc="left", color=INK)

    # B/C: margin maps per family with worst-K marks
    for ax, fam, lab in ((axB, 1.0, "n$_y$ = +1"),
                         (axC, 0.75, "n$_y$ = +0.75")):
        m = pop["ny"] == fam
        M = pop["margin"][m].reshape(NT, NZ)   # chunk-major: t0 x z0
        pc = ax.pcolormesh(TGRID, ZGRID, np.minimum(M, 6.0).T,
                           cmap="viridis", shading="nearest")
        plt.colorbar(pc, ax=ax, label="margin (capped at 6 for display)")
        wsel = wk[pop["ny"][wk] == fam]
        ax.plot(pop["t0"][wsel], pop["z0"][wsel], "o", ms=5, mfc="none",
                mec=ACC, mew=1.4, label="worst-K (ladder-certified)")
        ax.set_xlabel("t$_0$ (kill bins)")
        ax.set_ylabel("z$_0$")
        ax.set_title("%s  margin map, family %s  (worst-K circled)"
                     % ("B" if fam == 1.0 else "C", lab),
                     fontsize=10, loc="left", color=INK)
        ax.legend(fontsize=8, frameon=False, loc="upper right")

    # D: screener validation + ladder agreement
    rels = [p["rel_diff"] for p in val["per_config"]]
    axD.semilogy(range(len(rels)), sorted(rels), "o", ms=5, color=BLUD,
                 label="G1: screener vs P2 DOP853 (24 configs)")
    lr = [p["rel_scr_vs_mp"] for p in res["G3"]["per_config"]]
    axD.semilogy(range(len(lr)), sorted(lr), "s", ms=5, color=ACC,
                 mfc="none", label="G3: screener vs 200-bit ladder (24)")
    hv = res["G3"]["step_halving_rel_shifts"]
    axD.semilogy(range(len(hv)), hv, "^", ms=6, color=GRN,
                 label="mp step-halving rel shift (3)")
    axD.axhline(1e-6, color=INK2, lw=1.0, ls="--")
    axD.text(1, 1.5e-6, "G1 gate 1e-6", fontsize=8, color=INK2)
    axD.set_xlabel("config (sorted)")
    axD.set_ylabel("relative margin disagreement")
    axD.set_title("D  screener validation + worst-K certification",
                  fontsize=10, loc="left", color=INK)
    axD.legend(fontsize=8, frameon=False)

    for ax in axes.ravel():
        ax.set_facecolor(SURF)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.suptitle("S4 — population-scale certification screener "
                 "(P2 extension, F-T8-S4)", fontsize=11.5, color=INK,
                 x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(S4DIR, "s4_fig.png"), facecolor=SURF,
                bbox_inches="tight")
    print("figure written", flush=True)


# ---------------------------------------------------------------------
def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    steps = dict(bench=phase_bench, validate=phase_validate,
                 screen=phase_screen, ladder=phase_ladder,
                 gates=phase_gates, fig=phase_fig)
    if phase == "all":
        for name, fn in steps.items():
            print("== phase %s ==" % name, flush=True)
            fn()
    else:
        print("== phase %s ==" % phase, flush=True)
        steps[phase]()
    print("total elapsed %.1f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
