#!/usr/bin/env python3
"""
WORKSTREAM S3 -- O-1 pilot: validated FORWARD flow-map covering
(Q1 technology extension).  File as F-T8-S3.

Q2's ANALYTIC_RECON.md SS5 names obligation O-1 (early-time flow control on
[0, T0] x supp rho_0) as THE dominant obligation of Conjecture C, with "the
only identified technology" being the Q1 validated-enclosure covering.  Q1
certified 7 single BACKWARD paths.  This pilot extends the same engine --
mpmath.iv/disc-interval band-limited field, verified Picard a-priori boxes,
Lohner QR mean-value propagation -- from thin backward tubes to FAT FORWARD
boxes, and runs a pre-registered 24-box covering of the Q2 T_ball band.

A pilot is NOT an O-1 discharge.  The deliverable is the FEASIBILITY DATUM:
fraction of the subdomain certified crossed/localized, per-box validated-step
cost, wrapping/width growth factors, and the honest extrapolated cost of a
full O-1 covering.

THE CERTIFIED OBJECT (identical to O1/P2/Q1): the exact Nx = 2048
band-limited field f(x,t) = sum_k C_k e^{i kappa_k (x - xmin)}
e^{-i kappa_k^2 t/2} with the analytic z-factor; guidance velocity
    v_x = Im(f'/f) - n_y pi cot(pi z),   v_z = n_y Re(f'/f).
All interval semantics (scalar outward-rounded float64 pairs, complex disc
arithmetic, mpmath.iv prec-120 transcendentals, Lohner QR with rigorous 2x2
adjugate inverse) are IMPORTED UNCHANGED from tier7-program/Q1/q1_enclose.py.
The only changes are (i) the time direction: forward step s0 -> s1 = s0 + h
via the exact identity
    X(s1) = X0 + h v(X0, s0) + int_{s0}^{s1} (s1 - u) D(X(u), u) du,
    (s1 - u) in [0, h]  =>  X(s1) in X0 + h v(X0,s0) + (h^2/2) hull(D(B,I_s))
with S = I + h J(B, I_s) in the Lohner mean-value form (Q1 used the mirrored
backward identity with S = I - h J); (ii) FAT initial data: the Lohner frame
starts at A = I with q = box - center (an interval RECTANGLE, not a point);
(iii) the Picard a-priori pad is h-scaled (0.35 w(h v) + 4e-3 h + 1e-14,
inflation x2 on failure) instead of Q1's tube-width-scaled pad -- for fat
boxes Q1's 0.55*width pad would inflate B by ~55% of the box size and
artificially accelerate wrapping; the pad is heuristic, the Picard
CONTAINMENT CHECK (the rigorous part) is unchanged.

SUBDOMAIN (pre-registered, ROADMAP_v10 S3): n_y = 1; the Q2 T_ball band
w in [-1.2, -0.6], z in [0.05, 0.95] at t0 = 1.0, mapped to (x, z) via the
Q2 lag/lead label w = (x - x0)/t - k0 (x0 = -5, k0 = 2), i.e. at t0 = 1.0:
    x = x0 + (k0 + w) t0 = w - 3   =>   x in [-4.2, -3.6].
COORDINATE CLARIFICATION (printed per spec): the (w, z) rectangle maps to an
(x, z) rectangle exactly because t0 is a single time slice; the pilot box
grid is 6 x 4 over x in [-4.2, -3.6] (width 0.1 each) x z in [0.05, 0.95]
(width 0.225 each), M = 24.  Forward horizon T0 = 5.14 (~ the measured
T_max(1) = 5.1310 of Conjecture C; L5' tau_max = 5.130950).

CLASSIFICATION per box (pre-registered):
  (i)  CROSSED    -- enclosure wholly past x = d = 1: lower(X_x) > d at some
                     s <= T0 (all member trajectories have crossed by s);
  (ii) LOCALIZED  -- wholly x < d at T0: upper(X_x) < d at s = T0,
                     enclosure printed;
  (iii) BLOWUP    -- enclosure width > 1 before T0; also grouped here the
                     engine-terminal failures FAIL_WALL (z interval reaches
                     the wall lane, cot/csc singular), FAIL_FDISC (f-disc
                     contains 0), FAIL_PICARD (no verified a-priori box).
  A box reaching T0 with width <= 1 but straddling x = d is UNRESOLVED
  (counted with FAIL in the fraction accounting -- spec clarification).
  Failure-class parents are subdivided ONCE 2x2, children re-run, then stop.

FEASIBILITY DIAGNOSTIC (in addition to the 24-box covering; labeled as
diagnostic, not part of the pre-registered covering): a width LADDER at the
band center -- single boxes of decreasing width run to terminal status --
locating the maximal initial width w* the engine certifies to T0, which
anchors the full-O-1 extrapolation.

GATES (pre-registered): S3-G1 forward replay of Q1 path 134 contains the Q1
anchor (containment, not distance); S3-G2 covering executed, zero
unaccounted subdomain; S3-G3 feasibility datum incl. full-O-1 extrapolation;
S3-G4 scope honesty (pilot != discharge; kill-bin theorem remains open).

Usage: s3_cover.py [g1|bench|cover|ladder|gates|fig|produce]
Per-box results _s3_box_<tag>.json; running checkpoints _s3_ck_<tag>.json
(atomic, resumable); production under nohup.  hbar = m = 1.
Within-model; nothing here bears on nature.
"""
import json
import math
import os
import sys
import time

import numpy as np

S3DIR = os.path.dirname(os.path.abspath(__file__))
Q1DIR = os.path.normpath(os.path.join(
    S3DIR, "..", "..", "tier7-program", "Q1"))
sys.path.insert(0, Q1DIR)

import q1_enclose as q1                                        # noqa: E402
from q1_enclose import (                                       # noqa: E402
    ipt, iadd, isub, imul, iscale, ineg, ihull, iwid, imid, isubset,
    i_isect, _lo, _hi, INF,
    vel_pack, mat_ivmul_float, mat_ivmul_iv, matvec_iv, matvec_float_iv,
    inv2_iv, qr2, D_NEAR)

# ----------------------------------------------------------------------
# S3 parameters (pre-registered)
# ----------------------------------------------------------------------
NY = 1.0
T_START = 1.0                       # exact float
T0_END = 5.14                       # fl(5.14); classification epoch
D = D_NEAR                          # detector x = 1.0
X_LO_BAND, X_HI_BAND = -4.2, -3.6   # w in [-1.2,-0.6] at t0=1  ->  x = w-3
Z_LO_BAND, Z_HI_BAND = 0.05, 0.95
NBX, NBZ = 6, 4                     # M = 24
H_COVER = 2.0 ** -12                # fixed by bench (printed there)
WIDTH_CAP_S3 = 1.0                  # pre-registered blow-up threshold
CK_EVERY = 4000
FAIL_CLASS = ("BLOWUP", "FAIL_WALL", "FAIL_FDISC", "FAIL_PICARD",
              "FAIL_TIMEOUT", "FAIL_EXC")
LADDER_WIDTHS = [0.05, 0.02, 0.01, 0.005, 0.0025, 0.00125,
                 0.000625, 0.0003125]

Q1_PATH_FOR_G1 = 134                # shallowest Q1 certificate (cheap replay)


def parent_grid():
    dx = (X_HI_BAND - X_LO_BAND) / NBX
    dz = (Z_HI_BAND - Z_LO_BAND) / NBZ
    out = []
    for i in range(NBX):
        for j in range(NBZ):
            out.append(dict(
                tag="p%d%d" % (i, j), level=0, i=i, j=j,
                x_lo=X_LO_BAND + i * dx, x_hi=X_LO_BAND + (i + 1) * dx,
                z_lo=Z_LO_BAND + j * dz, z_hi=Z_LO_BAND + (j + 1) * dz,
                h=H_COVER))
    return out


def children_of(cfg):
    xm = 0.5 * (cfg["x_lo"] + cfg["x_hi"])
    zm = 0.5 * (cfg["z_lo"] + cfg["z_hi"])
    out = []
    for ci, (xl, xh) in enumerate([(cfg["x_lo"], xm), (xm, cfg["x_hi"])]):
        for cj, (zl, zh) in enumerate([(cfg["z_lo"], zm), (zm, cfg["z_hi"])]):
            out.append(dict(
                tag="%s_c%d%d" % (cfg["tag"], ci, cj), level=1,
                i=cfg["i"], j=cfg["j"],
                x_lo=xl, x_hi=xh, z_lo=zl, z_hi=zh, h=cfg["h"]))
    return out


# ----------------------------------------------------------------------
# forward Lohner box integration
# ----------------------------------------------------------------------
def forward_box(cfg, quiet=True, max_steps=None, max_hours=1.5,
                classify=True, t_ref=None, k_off=0, n_total=None,
                t_final=None, resfile=True):
    """Integrate the initial rectangle [x_lo,x_hi]x[z_lo,z_hi] FORWARD.

    Time grid: s(k) = t_ref + (k - k_off) h  (one rounding; exact when
    t_ref, h binary and exponents align -- true for the covering grid
    t_ref = 1.0).  After the last full step a single partial step lands
    exactly on t_final (outward-rounded time intervals keep this rigorous).
    Default (classify=True): t_ref = T_START, t_final = T0_END.
    """
    tag, h = cfg["tag"], cfg["h"]
    if t_ref is None:
        t_ref = T_START
    if t_final is None and classify:
        t_final = T0_END
    if n_total is None:
        n_total = int(math.floor((t_final - t_ref) / h)) if t_final else 0
    ckfile = os.path.join(S3DIR, "_s3_ck_%s.json" % tag)
    REC_EVERY = max(1, int(round(2.0 ** -6 / h)))    # ~64 samples / unit s

    x_int = (cfg["x_lo"], cfg["x_hi"])
    z_int = (cfg["z_lo"], cfg["z_hi"])
    w_init = max(iwid(x_int), iwid(z_int))
    c = [imid(x_int), imid(z_int)]
    A = [[1.0, 0.0], [0.0, 1.0]]
    q = [isub(x_int, ipt(c[0])), isub(z_int, ipt(c[1]))]
    box = [x_int, z_int]
    k = 0
    rec = []
    mon = dict(min_rho_x_lo=INF, max_g_up=0.0, z_lo=z_int[0], z_hi=z_int[1],
               max_wid=w_init, picard_retries=0)
    wall0 = time.time()
    if os.path.exists(ckfile):
        ck = json.load(open(ckfile))
        k = ck["k"]
        c = [float.fromhex(v) for v in ck["c"]]
        A = [[float.fromhex(v) for v in row] for row in ck["A"]]
        q = [tuple(float.fromhex(v) for v in pair) for pair in ck["q"]]
        box = [tuple(float.fromhex(v) for v in pair) for pair in ck["box"]]
        rec = ck["rec"]
        mon = ck["mon"]
        if not quiet:
            print("  [%s] resumed at step %d" % (tag, k), flush=True)

    status, s_term = None, None
    vg = None                        # warm-started Picard velocity guess
    while True:
        last = (k >= n_total)
        s0 = t_ref + (k - k_off) * h
        s1 = (t_final if last else t_ref + (k + 1 - k_off) * h)
        if last and (t_final is None or s1 <= s0):
            status = "DONE_GRID"     # replay mode: grid exhausted
            s_term = s0
            break
        s0_i = (_lo(s0), _hi(s0))
        s1_i = (_lo(s1), _hi(s1))
        is_i = ihull(s0_i, s1_i)
        hh = s1 - s0                 # = h except on the final partial step
        if max_steps is not None and k >= max_steps:
            status = "BENCH_STOP"
            s_term = s0
            break
        if (time.time() - wall0) > max_hours * 3600:
            status = "FAIL_TIMEOUT"
            s_term = s0
            break

        try:
            # thin velocity at the box center
            pc = vel_pack(ipt(c[0]), ipt(c[1]), s0_i, NY, 1)
            if vg is None:
                pb = vel_pack(box[0], box[1], s0_i, NY, 1)
                vg = [pb["vx"], pb["vz"]]
            # verified Picard a-priori box (forward: [0,h] * +v)
            ok = False
            infl = 1.0
            for attempt in range(8):
                B = []
                for d_ in range(2):
                    mv = ihull(ipt(0.0), iscale(vg[d_], hh))
                    Bd = iadd(box[d_], mv)
                    pad = infl * (0.35 * iwid(mv) + 4.0 * hh * 1e-3 + 1e-14)
                    B.append((_lo(Bd[0] - pad), _hi(Bd[1] + pad)))
                pf = vel_pack(B[0], B[1], is_i, NY, 3)
                V = [pf["vx"], pf["vz"]]
                cont = all(isubset(iadd(box[d_],
                                        ihull(ipt(0.0), iscale(V[d_], hh))),
                                   B[d_]) for d_ in range(2))
                if cont:
                    vg = V
                    ok = True
                    break
                vg = V
                infl *= 2.0
                mon["picard_retries"] += 1
            if not ok:
                status = "FAIL_PICARD"
                s_term = s0
                break
        except ZeroDivisionError:
            status = "FAIL_FDISC"
            s_term = s0
            break
        except RuntimeError as e:
            status = ("FAIL_WALL" if "z interval" in str(e)
                      or "sin(pi z)" in str(e) else "FAIL_EXC")
            s_term = s0
            break

        # forward step: u_new = c + h v(c,s0) + (h^2/2) D(B, I_s)
        # (h^2/2 applied as two outward-rounded scalings: 0.5*hh is an
        #  exact power-of-two scaling, so no pre-rounded h^2 constant)
        u_new = [
            iadd(iadd(ipt(c[0]), iscale(pc["vx"], hh)),
                 iscale(iscale(pf["Dx"], 0.5 * hh), hh)),
            iadd(iadd(ipt(c[1]), iscale(pc["vz"], hh)),
                 iscale(iscale(pf["Dz"], 0.5 * hh), hh)),
        ]
        # S = I + h J(B, I_s)      (Jzz = 0 exactly: v_z = n_y b(x,t))
        S = [[iadd(ipt(1.0), iscale(pf["Jxx"], hh)),
              iscale(pf["Jxz"], hh)],
             [iscale(pf["Jzx"], hh), ipt(1.0)]]
        M = mat_ivmul_float(S, A)
        Aq = matvec_float_iv(A, q)
        direct = [iadd(u_new[d_], matvec_iv(S, Aq)[d_]) for d_ in range(2)]
        Mm = [[imid(M[i][j]) for j in range(2)] for i in range(2)]
        A2 = qr2(Mm)
        Ainv = inv2_iv(A2)
        c2 = [imid(u_new[0]), imid(u_new[1])]
        e = [isub(u_new[d_], ipt(c2[d_])) for d_ in range(2)]
        BQ = mat_ivmul_iv(Ainv, M)
        q2v = [iadd(matvec_iv(BQ, q)[d_], matvec_iv(Ainv, e)[d_])
               for d_ in range(2)]
        lohner = [iadd(ipt(c2[d_]), matvec_float_iv(A2, q2v)[d_])
                  for d_ in range(2)]
        box2 = [i_isect(lohner[d_], direct[d_]) for d_ in range(2)]
        c, A, q, box = c2, A2, q2v, box2
        k += 1

        w0, w1 = iwid(box[0]), iwid(box[1])
        mon["max_wid"] = max(mon["max_wid"], w0, w1)
        mon["min_rho_x_lo"] = min(mon["min_rho_x_lo"], pf["rho_x_lo"])
        mon["max_g_up"] = max(mon["max_g_up"], pf["g_up"])
        mon["z_lo"] = min(mon["z_lo"], box[1][0])
        mon["z_hi"] = max(mon["z_hi"], box[1][1])
        if k % REC_EVERY == 0 or last:
            rec.append([s1, box[0][0], box[0][1], box[1][0], box[1][1],
                        w0, w1])

        if classify:
            if box[0][0] > D:
                status = "CROSSED"
                s_term = s1
                break
            if max(w0, w1) > WIDTH_CAP_S3:
                status = "BLOWUP"
                s_term = s1
                break
            if last:
                status = "LOCALIZED" if box[0][1] < D else "UNRESOLVED"
                s_term = s1
                break
        elif last:
            status = "DONE_GRID"
            s_term = s1
            break

        if k % CK_EVERY == 0:
            _ck = dict(k=k, c=[v.hex() for v in c],
                       A=[[v.hex() for v in row] for row in A],
                       q=[[v[0].hex(), v[1].hex()] for v in q],
                       box=[[v[0].hex(), v[1].hex()] for v in box],
                       rec=rec, mon=mon)
            json.dump(_ck, open(ckfile + ".tmp", "w"))
            os.replace(ckfile + ".tmp", ckfile)
            if not quiet:
                print("  [%s] step %d s=%.4f x=[%.4f,%.4f] w=%.2e %.0fs"
                      % (tag, k, s1, box[0][0], box[0][1],
                         max(w0, w1), time.time() - wall0), flush=True)

    out = dict(
        tag=tag, level=cfg.get("level", 0), h=h, ny=NY,
        init_box=[[cfg["x_lo"], cfg["x_hi"]], [cfg["z_lo"], cfg["z_hi"]]],
        w_init=w_init, status=status, n_steps=k, s_term=s_term,
        final_box=[[box[0][0], box[0][1]], [box[1][0], box[1][1]]],
        final_wid=[iwid(box[0]), iwid(box[1])],
        growth_factor=mon["max_wid"] / w_init if w_init > 0 else None,
        mon=mon, rec=rec, elapsed=time.time() - wall0)
    if resfile:
        fn = os.path.join(S3DIR, "_s3_box_%s.json" % tag)
        json.dump(out, open(fn + ".tmp", "w"))
        os.replace(fn + ".tmp", fn)
    if os.path.exists(ckfile):
        os.remove(ckfile)
    if not quiet:
        print("  [%s] %s steps=%d s_term=%s max_wid=%.3g (%.0f s)"
              % (tag, status, k, s_term, mon["max_wid"],
                 time.time() - wall0), flush=True)
    return out


# ----------------------------------------------------------------------
# phase G1: forward replay of Q1 path 134 -> must contain the anchor
# ----------------------------------------------------------------------
def phase_g1():
    fn = os.path.join(Q1DIR, "_q1_path_%d.json" % Q1_PATH_FOR_G1)
    p = json.load(open(fn))
    assert p["status"] == "CERTIFIED"
    K = p["n_steps"]
    hq = p["h"]
    t0 = p["t0"]
    fb = p["final_box"]
    anchor = (D_NEAR, p["z0"])
    print("G1: forward replay of Q1 path %d: t0=%.6g z0=%.6g h=2^%d "
          "K=%d s_stop=%.10f" % (Q1_PATH_FOR_G1, t0, p["z0"],
                                 round(math.log2(hq)), K, p["s_stop"]),
          flush=True)
    print("    start box (Q1 endpoint enclosure): x=[%.12f,%.12f] "
          "z=[%.12f,%.12f]" % (fb[0][0], fb[0][1], fb[1][0], fb[1][1]),
          flush=True)
    cfg = dict(tag="g1replay", level=0, x_lo=fb[0][0], x_hi=fb[0][1],
               z_lo=fb[1][0], z_hi=fb[1][1], h=hq)
    # time grid s(k) = t0 - (K - k) h : bitwise-identical floats to Q1's
    r = forward_box(cfg, quiet=False, classify=False, t_ref=t0, k_off=K,
                    n_total=K, t_final=None, resfile=False)
    bx, bz = r["final_box"]
    cont = (bx[0] <= anchor[0] <= bx[1]) and (bz[0] <= anchor[1] <= bz[1])
    out = dict(
        q1_path=Q1_PATH_FOR_G1, t0=t0, z0=p["z0"], h=hq, n_steps=r["n_steps"],
        start_box=fb, anchor=list(anchor),
        end_box=r["final_box"], end_wid=r["final_wid"],
        anchor_contained=bool(cont),
        growth_factor=r["growth_factor"], elapsed=r["elapsed"],
        note=("forward integration of the time-reverse of the certified "
              "backward path: the true trajectory through the anchor lies "
              "in the Q1 endpoint enclosure at s_stop, so a valid forward "
              "enclosure of that box at t0 MUST contain the anchor"),
        verdict="PASS" if cont else "FAIL")
    json.dump(out, open(os.path.join(S3DIR, "_s3_g1.json"), "w"), indent=2)
    print("G1 %s: end box x=[%.8f,%.8f] z=[%.8f,%.8f] contains anchor "
          "(%.4f, %.4f): %s  (wid %.2e/%.2e, growth %.1fx, %.0f s)"
          % (out["verdict"], bx[0], bx[1], bz[0], bz[1], anchor[0],
             anchor[1], cont, r["final_wid"][0], r["final_wid"][1],
             r["growth_factor"], r["elapsed"]), flush=True)
    return out


# ----------------------------------------------------------------------
# phase BENCH: one box, limited steps -> unit cost + h decision
# ----------------------------------------------------------------------
def phase_bench():
    grid = parent_grid()
    cfg = dict(grid[13])             # p31: x [-3.9,-3.8], z [0.275,0.5]
    cfg["tag"] = "bench"
    nb = 1000
    t0w = time.time()
    r = forward_box(cfg, quiet=False, max_steps=nb, resfile=False)
    dt = (time.time() - t0w) / max(r["n_steps"], 1)
    n_full = int(math.floor((T0_END - T_START) / H_COVER)) + 1
    est_full = dt * n_full
    est_cover = est_full * (NBX * NBZ) / 4          # 4 workers, worst case
    out = dict(bench_box=cfg["tag"], base=grid[13]["tag"], h=H_COVER,
               n_steps=r["n_steps"], status=r["status"],
               per_step_ms=dt * 1e3,
               s_reached=r["s_term"], max_wid=r["mon"]["max_wid"],
               growth_factor=r["growth_factor"],
               n_steps_full_horizon=n_full,
               est_seconds_full_box=est_full,
               est_seconds_cover_4workers=est_cover,
               decision=("h = 2^-12 retained: per-step cost %.1f ms, full "
                         "4.14-unit horizon %d steps = %.0f s/box worst "
                         "case; blow-up-terminated boxes cost far less. "
                         "24 parents + <= 96 children fit the wall-clock "
                         "budget on 4 workers." % (dt * 1e3, n_full,
                                                   est_full)))
    json.dump(out, open(os.path.join(S3DIR, "_s3_bench.json"), "w"),
              indent=2)
    print("BENCH:", json.dumps({k: v for k, v in out.items()
                                if k != "decision"}), flush=True)
    print("BENCH decision:", out["decision"], flush=True)
    return out


# ----------------------------------------------------------------------
# phase COVER: 24 parents (Pool 4), subdivide failures once, children
# ----------------------------------------------------------------------
def _worker(cfg):
    try:
        return forward_box(cfg, quiet=True)
    except Exception as exc:                        # noqa: BLE001
        import traceback
        traceback.print_exc()
        out = dict(tag=cfg["tag"], level=cfg.get("level", 0),
                   status="FAIL_EXC_" + str(exc)[:60], n_steps=0,
                   init_box=[[cfg["x_lo"], cfg["x_hi"]],
                             [cfg["z_lo"], cfg["z_hi"]]],
                   w_init=max(cfg["x_hi"] - cfg["x_lo"],
                              cfg["z_hi"] - cfg["z_lo"]),
                   s_term=None, final_box=None, final_wid=None,
                   growth_factor=None, mon={}, rec=[], elapsed=0.0)
        fn = os.path.join(S3DIR, "_s3_box_%s.json" % cfg["tag"])
        json.dump(out, open(fn, "w"))
        return out


def _run_pool(cfgs, label):
    from multiprocessing import Pool
    todo = [c for c in cfgs if not os.path.exists(
        os.path.join(S3DIR, "_s3_box_%s.json" % c["tag"]))]
    print("%s: %d boxes, %d to run" % (label, len(cfgs), len(todo)),
          flush=True)
    if todo:
        with Pool(min(4, len(todo))) as pool:
            for r in pool.imap_unordered(_worker, todo):
                print("done %s -> %s (steps %d, s_term %s, max_wid %s)"
                      % (r["tag"], r["status"], r["n_steps"], r["s_term"],
                         (None if not r.get("mon") else
                          "%.3g" % r["mon"]["max_wid"])), flush=True)
    return {c["tag"]: json.load(open(os.path.join(
        S3DIR, "_s3_box_%s.json" % c["tag"]))) for c in cfgs}


def phase_cover():
    grid = parent_grid()
    parents = _run_pool(grid, "COVER parents")
    subdiv = [c for c in grid
              if any(parents[c["tag"]]["status"].startswith(f)
                     for f in FAIL_CLASS)]
    kids = []
    for c in subdiv:
        kids.extend(children_of(c))
    print("COVER: %d parents in failure class -> %d children"
          % (len(subdiv), len(kids)), flush=True)
    if kids:
        _run_pool(kids, "COVER children")
    print("COVER complete", flush=True)


# ----------------------------------------------------------------------
# phase LADDER: max feasible initial width at band center (diagnostic)
# ----------------------------------------------------------------------
def phase_ladder():
    xc, zc = -3.9, 0.5
    for wx in LADDER_WIDTHS:
        wz = 2.25 * wx               # keep the 6x4 grid aspect ratio
        tag = "lad%s" % ("%g" % wx).replace(".", "p")
        fn = os.path.join(S3DIR, "_s3_box_%s.json" % tag)
        if os.path.exists(fn):
            r = json.load(open(fn))
        else:
            cfg = dict(tag=tag, level=9, x_lo=xc - wx / 2, x_hi=xc + wx / 2,
                       z_lo=zc - wz / 2, z_hi=zc + wz / 2, h=H_COVER)
            r = forward_box(cfg, quiet=False)
        print("LADDER wx=%g wz=%g -> %s (steps %d, s_term %s, growth %s)"
              % (wx, wz, r["status"], r["n_steps"], r["s_term"],
                 ("%.1f" % r["growth_factor"]) if r["growth_factor"]
                 else "-"), flush=True)
        if r["status"] in ("CROSSED", "LOCALIZED"):
            print("LADDER: first certifying width w* = %g (x) found; stop"
                  % wx, flush=True)
            break
    print("LADDER complete", flush=True)


# ----------------------------------------------------------------------
# phase TRUEGROWTH: float-layer diagnostic (NOT rigorous, labeled as such)
# true flow-map amplification ||DPhi_{1->t}|| vs enclosure growth
# = the wrapping-effect factor the roadmap's feasibility datum asks for
# ----------------------------------------------------------------------
def phase_truegrowth():
    from scipy.integrate import solve_ivp
    KV = (np.pi / 30) * q1.MVALS

    def f012(x, t):
        ph = np.exp(1j * (KV * (x - q1.XMIN) - 0.5 * KV ** 2 * t))
        f = np.sum(q1.CFT * ph)
        f1 = np.sum(q1.CFT * 1j * KV * ph)
        f2 = np.sum(q1.CFT * (1j * KV) ** 2 * ph)
        return f, f1, f2

    def rhs(t, y):
        x, z = y[0], y[1]
        V = y[2:].reshape(2, 2)
        f, f1, f2 = f012(x, t)
        g = f1 / f
        G2 = f2 / f - g * g
        zc = min(max(z, 1e-6), 1 - 1e-6)
        vx = g.imag - np.pi / np.tan(np.pi * zc)
        vz = g.real
        J = np.array([[G2.imag, np.pi ** 2 / np.sin(np.pi * zc) ** 2],
                      [G2.real, 0.0]])
        return np.concatenate([[vx, vz], (J @ V).ravel()])

    # ladder death times + enclosure growth at death (from the rec logs)
    lad = []
    for wx in LADDER_WIDTHS:
        tag = "lad%s" % ("%g" % wx).replace(".", "p")
        r = _load(tag)
        if r is None:
            continue
        lad.append(dict(wx=wx, t_death=r["s_term"],
                        G_enc=r["growth_factor"], status=r["status"]))

    pts = [(-3.9, 0.5)] + [(x, z) for x in (-4.2, -3.9, -3.6)
                           for z in (0.1, 0.5, 0.9) if (x, z) != (-3.9, 0.5)]
    rows = []
    for (x0, z0) in pts:
        y0 = np.array([x0, z0, 1.0, 0.0, 0.0, 1.0])
        tev = sorted(set([l["t_death"] for l in lad] + [T0_END]))
        try:
            sol = solve_ivp(rhs, (T_START, T0_END), y0, method="DOP853",
                            rtol=1e-10, atol=1e-12, t_eval=tev,
                            dense_output=False, max_step=0.05)
            ok = sol.success and sol.t[-1] >= T0_END - 1e-9
            amps = {}
            wall = False
            for ti, yi in zip(sol.t, sol.y.T):
                if yi[1] > 0.999:
                    wall = True
                V = yi[2:].reshape(2, 2)
                amps["%.6f" % ti] = float(np.linalg.norm(V, 2))
            rows.append(dict(x0=x0, z0=z0, ok=bool(ok), wall_lane=wall,
                             xz_final=[float(sol.y[0, -1]),
                                       float(sol.y[1, -1])],
                             A_true=amps))
        except Exception as exc:                    # noqa: BLE001
            rows.append(dict(x0=x0, z0=z0, ok=False, error=str(exc)[:80]))

    # wrapping factor at the band center vs each ladder death time
    ctr = rows[0]
    wrap = []
    for l in lad:
        key = "%.6f" % l["t_death"]
        if ctr.get("A_true") and key in ctr["A_true"]:
            at = ctr["A_true"][key]
            wrap.append(dict(wx=l["wx"], t_death=l["t_death"],
                             G_enc=l["G_enc"], A_true_center=at,
                             wrapping_factor=l["G_enc"] / at))
    a_finals = [r["A_true"]["%.6f" % T0_END] for r in rows
                if r.get("ok") and r.get("A_true")]
    out = dict(
        note=("FLOAT-LAYER DIAGNOSTIC, not rigorous: DOP853 rtol 1e-10 "
              "trajectories + variational matrix V' = J V; A_true = "
              "||DPhi_{1->t}||_2. Enclosure growth divided by A_true = "
              "wrapping/interval overhead of the pilot engine."),
        points=rows, ladder_wrap=wrap,
        A_true_T0=dict(n=len(a_finals),
                       min=min(a_finals) if a_finals else None,
                       median=(sorted(a_finals)[len(a_finals) // 2]
                               if a_finals else None),
                       max=max(a_finals) if a_finals else None))
    json.dump(out, open(os.path.join(S3DIR, "_s3_truegrowth.json"), "w"),
              indent=1)
    print("TRUEGROWTH: A_true(1 -> 5.14) over %d pts: min %.3g med %.3g "
          "max %.3g" % (len(a_finals), out["A_true_T0"]["min"],
                        out["A_true_T0"]["median"], out["A_true_T0"]["max"]),
          flush=True)
    for w in wrap:
        print("  wx=%g t_death=%.3f G_enc=%.1f A_true=%.2f wrap=%.1fx"
              % (w["wx"], w["t_death"], w["G_enc"], w["A_true_center"],
                 w["wrapping_factor"]), flush=True)
    return out


# ----------------------------------------------------------------------
# phase GATES
# ----------------------------------------------------------------------
def _load(tag):
    fn = os.path.join(S3DIR, "_s3_box_%s.json" % tag)
    return json.load(open(fn)) if os.path.exists(fn) else None


def phase_gates():
    g1 = json.load(open(os.path.join(S3DIR, "_s3_g1.json")))
    bench = json.load(open(os.path.join(S3DIR, "_s3_bench.json")))
    grid = parent_grid()

    # leaf accounting: every parent either terminal or replaced by its 4
    # children; leaf areas must sum to the exact band area
    band_area = (X_HI_BAND - X_LO_BAND) * (Z_HI_BAND - Z_LO_BAND)
    leaves = []
    missing = []
    for c in grid:
        p = _load(c["tag"])
        if p is None:
            missing.append(c["tag"])
            continue
        if any(p["status"].startswith(f) for f in FAIL_CLASS):
            for ch in children_of(c):
                r = _load(ch["tag"])
                if r is None:
                    missing.append(ch["tag"])
                else:
                    r["parent"] = c["tag"]
                    leaves.append(r)
            p["subdivided"] = True
            p["is_leaf"] = False
        else:
            p["is_leaf"] = True
            leaves.append(p)
    area = sum((r["init_box"][0][1] - r["init_box"][0][0])
               * (r["init_box"][1][1] - r["init_box"][1][0])
               for r in leaves)

    def _cls(r):
        s = r["status"]
        if s == "CROSSED":
            return "CROSSED"
        if s == "LOCALIZED":
            return "LOCALIZED"
        if s == "UNRESOLVED":
            return "UNRESOLVED"
        return "FAIL"

    frac = {}
    for r in leaves:
        a = ((r["init_box"][0][1] - r["init_box"][0][0])
             * (r["init_box"][1][1] - r["init_box"][1][0]))
        frac[_cls(r)] = frac.get(_cls(r), 0.0) + a / band_area
    counts = {}
    for r in leaves:
        counts[_cls(r)] = counts.get(_cls(r), 0) + 1
    statuses = {}
    for r in leaves:
        statuses[r["status"]] = statuses.get(r["status"], 0) + 1

    steps = sorted(r["n_steps"] for r in leaves)
    walls = sorted(r["elapsed"] for r in leaves)
    growths = sorted(r["growth_factor"] for r in leaves
                     if r.get("growth_factor"))
    med = lambda v: v[len(v) // 2] if v else None    # noqa: E731

    # ladder summary
    ladder = []
    wstar = None
    for wx in LADDER_WIDTHS:
        tag = "lad%s" % ("%g" % wx).replace(".", "p")
        r = _load(tag)
        if r is None:
            continue
        ladder.append(dict(wx=wx, wz=2.25 * wx, status=r["status"],
                           n_steps=r["n_steps"], s_term=r["s_term"],
                           growth_factor=r.get("growth_factor"),
                           final_wid=r.get("final_wid"),
                           elapsed=r["elapsed"]))
        if r["status"] in ("CROSSED", "LOCALIZED") and wstar is None:
            wstar = wx

    # full O-1 extrapolation (honest; assumptions printed)
    supp_x, z_span = 7.0, 0.9        # supp rho_0 x in [-8.5,-1.5]; z band
    extrap = dict(assumptions=(
        "full O-1 needs [0, T0] x supp rho_0: x-span 7.0 (supp A0 width, "
        "[-8.5,-1.5]) x z in [0.05, 0.95] (the same near-wall-strip "
        "exception class as O-2/O1 flux bounds; the strips are NOT covered "
        "here), horizon [0, 5.14] vs the pilot's [1.0, 5.14] (~1.24x "
        "steps), per certified-|n_y| family (>= 2: 0.75, 1.0; mirrors by "
        "symmetry; a continuum-n_y treatment would add an interval "
        "dimension on top)"))
    if wstar is not None:
        lad_ok = next(l for l in ladder if l["wx"] == wstar)
        n_boxes = math.ceil(supp_x / wstar) * math.ceil(z_span
                                                        / (2.25 * wstar))
        cost_box = lad_ok["elapsed"] * (5.14 / 4.14)
        extrap.update(
            w_star_x=wstar, w_star_z=2.25 * wstar,
            n_boxes_one_family=n_boxes,
            certifying_box_wall_s=lad_ok["elapsed"],
            est_core_hours_one_family=n_boxes * cost_box / 3600,
            est_core_hours_two_families=2 * n_boxes * cost_box / 3600,
            note=("assumes the band-center certifying width w* holds "
                  "across the domain -- optimistic near field nodes and "
                  "the wall lane where widths must shrink further; an "
                  "adaptive covering would redistribute but not remove "
                  "the wrapping-driven scale"))
    else:
        # nothing certified: build the honest bracketed extrapolation from
        # the ladder death-time sequence t_d(w0) and the float-layer
        # true-growth diagnostic (if computed)
        deaths = [(l["wx"], l["s_term"]) for l in ladder
                  if l["s_term"] is not None]
        rates = []
        for (w1, t1), (w2, t2) in zip(deaths, deaths[1:]):
            if t2 > t1 and w2 < w1:
                rates.append(dict(t_mid=0.5 * (t1 + t2),
                                  lam=math.log(w1 / w2) / (t2 - t1)))
        lam_last = rates[-1]["lam"] if rates else None
        tg = None
        fn = os.path.join(S3DIR, "_s3_truegrowth.json")
        if os.path.exists(fn):
            tg = json.load(open(fn))
        proj = {}
        if lam_last and deaths:
            w_end, t_end_d = deaths[-1][0], deaths[-1][1]
            dt_rem = T0_END - t_end_d
            # model A: enclosure growth rate stays at the last measured
            # value (pessimistic: rate was still decreasing)
            wA = w_end * math.exp(-lam_last * dt_rem)
            # model B: rate decays ~ 1/t from the last measured point
            # (optimistic: matches the dispersive smoothing of the field)
            wB = w_end * math.exp(-lam_last * t_end_d
                                  * math.log(T0_END / t_end_d))
            for nm, wv in (("model_A_const_rate", wA),
                           ("model_B_rate_1_over_t", wB)):
                nb = (math.ceil(supp_x / wv)
                      * math.ceil(z_span / (2.25 * wv)))
                proj[nm] = dict(w_star_x_projected=wv, n_boxes=nb,
                                core_hours_at_107s_per_box=nb * 107 / 3600)
        extrap.update(
            w_star_x=None,
            ladder_death_rates=rates,
            projection=proj,
            A_true=None if tg is None else tg["A_true_T0"],
            wrapping_factors=None if tg is None else tg["ladder_wrap"],
            ideal_adaptive_bound=(None if tg is None or
                                  not tg["A_true_T0"]["max"] else dict(
                note=("wrapping-free lower bound: ANY uniform interval "
                      "covering must start below ~0.25/A_true to keep "
                      "final widths under the f-disc scale; adaptive "
                      "splitting pays this locally instead of globally "
                      "but the pilot engine's wrapping overhead (the "
                      "measured G_enc/A_true factors) comes on top"),
                w0=0.25 / tg["A_true_T0"]["max"],
                n_boxes=(math.ceil(supp_x / (0.25
                                             / tg["A_true_T0"]["max"]))
                         * math.ceil(z_span / (2.25 * 0.25
                                               / tg["A_true_T0"]["max"]))))),
            note=(
                "NO ladder width certified to T0: every rung down to "
                "wx = %g died on the f-disc-zero blow-up mode before "
                "t = 2.0; the projected certifying width and box count "
                "for ONE family are bracketed by models A/B above; "
                "full O-1 by uniform covering with THIS engine is "
                "infeasible -- the pilot's conclusion is that O-1 needs "
                "higher-order (Taylor-model) enclosures and adaptive "
                "subdivision, not more compute" % LADDER_WIDTHS[-1]))

    g2_ok = (not missing) and abs(area - band_area) < 1e-12
    g3_ok = True                     # datum printed below regardless
    gates = dict(
        S3_G1=dict(spec=("forward replay of one Q1-certified backward "
                         "path contains the Q1 anchor (containment)"),
                   measured=g1, verdict=g1["verdict"]),
        S3_G2=dict(spec=("covering executed: all 24 parents (+ <= 1 round "
                         "2x2 subdivisions) to T0 or terminal blow-up, "
                         "checkpointed, zero unaccounted subdomain"),
                   n_parents=len(grid), n_leaves=len(leaves),
                   missing=missing, leaf_area=area, band_area=band_area,
                   verdict="PASS" if g2_ok else "FAIL"),
        S3_G3=dict(spec=("feasibility datum: fractions, median step count "
                         "and wall cost, width growth, full-O-1 "
                         "extrapolation"),
                   verdict="PASS" if g3_ok else "FAIL"),
        S3_G4=dict(spec="scope honesty",
                   statement=("PILOT, NOT AN O-1 DISCHARGE: the kill-bin "
                              "theorem (T4-W5 analytic proof, Conjecture C "
                              "obligations O-1..O-5) remains OPEN; no "
                              "register or grade language changes; "
                              "within-model throughout"),
                   verdict="PASS"),
    )
    results = dict(
        workstream="S3", file_id="F-T8-S3",
        params=dict(ny=NY, t_start=T_START, T0=T0_END, d=D,
                    band_x=[X_LO_BAND, X_HI_BAND],
                    band_z=[Z_LO_BAND, Z_HI_BAND],
                    band_w=[-1.2, -0.6], grid=[NBX, NBZ], h=H_COVER,
                    width_cap=WIDTH_CAP_S3,
                    coordinate_note=("w = (x-x0)/t - k0 at t0=1.0 -> "
                                     "x = w - 3; rectangle maps to "
                                     "rectangle exactly on a time slice")),
        bench=bench,
        gates=gates,
        feasibility=dict(
            fractions_area=frac, counts=counts, statuses=statuses,
            n_leaves=len(leaves),
            median_steps_per_box=med(steps),
            median_wall_s_per_box=med(walls),
            total_wall_s=sum(walls),
            growth_factors=dict(
                min=growths[0] if growths else None,
                median=med(growths), max=growths[-1] if growths else None),
            ladder=ladder, w_star_x=wstar,
            full_O1_extrapolation=extrap),
        leaves=[{k: v for k, v in r.items() if k != "rec"}
                for r in leaves],
    )
    json.dump(results, open(os.path.join(S3DIR, "s3_results.json"), "w"),
              indent=1)
    print("GATES:", json.dumps({g: gates[g]["verdict"] for g in gates}),
          flush=True)
    print("fractions (area):", json.dumps(frac), flush=True)
    print("counts:", json.dumps(counts), " statuses:", json.dumps(statuses),
          flush=True)
    print("median steps/box:", med(steps), " median wall s/box:",
          None if not walls else round(med(walls), 1), flush=True)
    print("growth factors: min %s med %s max %s"
          % (growths[0] if growths else None, med(growths),
             growths[-1] if growths else None), flush=True)
    print("ladder w*:", wstar, flush=True)
    print("extrapolation:", json.dumps(extrap), flush=True)
    return results


# ----------------------------------------------------------------------
# phase FIG
# ----------------------------------------------------------------------
def phase_fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    res = json.load(open(os.path.join(S3DIR, "s3_results.json")))
    leaves = res["leaves"]
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    COL = dict(CROSSED="#2a78d6", LOCALIZED="#3e8f6b", BLOWUP="#c8503c",
               FAIL_WALL="#e0913f", FAIL_FDISC="#8d5bb8",
               FAIL_PICARD="#8d5bb8", UNRESOLVED="#9a9890")

    def colof(s):
        for k, v in COL.items():
            if s.startswith(k):
                return v
        return "#9a9890"

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.5), dpi=160)
    fig.patch.set_facecolor(SURF)
    axA, axB, axC = axes

    # A: subdomain map
    for r in leaves:
        (xl, xh), (zl, zh) = r["init_box"]
        axA.add_patch(Rectangle((xl, zl), xh - xl, zh - zl,
                                facecolor=colof(r["status"]),
                                edgecolor=SURF, lw=0.8, alpha=0.9))
        if r["level"] == 0:
            axA.add_patch(Rectangle((xl, zl), xh - xl, zh - zl,
                                    facecolor="none", edgecolor=INK,
                                    lw=0.5))
    axA.set_xlim(X_LO_BAND - 0.02, X_HI_BAND + 0.02)
    axA.set_ylim(Z_LO_BAND - 0.03, Z_HI_BAND + 0.03)
    axA.set_xlabel("x at t$_0$ = 1  (lag band w $\\in$ [-1.2, -0.6])")
    axA.set_ylabel("z")
    axA.set_title("A  covering classification at T$_0$ = 5.14\n"
                  "(children shown where parents were subdivided)",
                  fontsize=10, loc="left", color=INK)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="s", ls="", markersize=8,
                      color=v, label=k) for k, v in COL.items()
               if any(r["status"].startswith(k) for r in leaves)]
    axA.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper "
               "left", bbox_to_anchor=(0.0, -0.14), ncol=3)

    # B: width histories (need rec from per-box files)
    for r in leaves:
        fn = os.path.join(S3DIR, "_s3_box_%s.json" % r["tag"])
        rr = json.load(open(fn))
        if not rr.get("rec"):
            continue
        a = np.array(rr["rec"])
        axB.semilogy(a[:, 0], np.maximum(a[:, 5], a[:, 6]),
                     color=colof(r["status"]), lw=0.8, alpha=0.7)
    for tag, lab in [("lad0p01", "ladder 0.01"),
                     ("lad0p005", "ladder 0.005"),
                     ("lad0p0025", "ladder 0.0025")]:
        fn = os.path.join(S3DIR, "_s3_box_%s.json" % tag)
        if os.path.exists(fn):
            rr = json.load(open(fn))
            if rr.get("rec"):
                a = np.array(rr["rec"])
                axB.semilogy(a[:, 0], np.maximum(a[:, 5], a[:, 6]),
                             color=INK, lw=1.3, ls="--", label=lab)
    axB.axhline(WIDTH_CAP_S3, color=INK2, lw=1.0, ls=":")
    axB.text(1.1, 1.15, "blow-up cap (width 1)", fontsize=7, color=INK2)
    axB.set_xlabel("t")
    axB.set_ylabel("max enclosure width")
    axB.set_title("B  width growth along the forward flow\n"
                  "(dashed: band-center width ladder)", fontsize=10,
                  loc="left", color=INK)
    if axB.get_legend_handles_labels()[0]:
        axB.legend(fontsize=6.5, frameon=False)

    # C: ladder bar: terminal status vs initial width
    lad = res["feasibility"]["ladder"]
    if lad:
        xs = np.arange(len(lad))
        gs = [l["growth_factor"] or 0 for l in lad]
        cols = [colof(l["status"]) for l in lad]
        axC.bar(xs, gs, color=cols, width=0.62)
        axC.set_yscale("log")
        axC.set_xticks(xs, ["%g" % l["wx"] for l in lad])
        for x, l in zip(xs, lad):
            axC.text(x, max(gs) * 1.1, l["status"][:9], rotation=90,
                     fontsize=6, ha="center", color=INK)
        axC.set_xlabel("initial x-width (z-width = 2.25x)")
        axC.set_ylabel("width growth factor at termination")
        axC.set_title("C  width ladder at band center\n"
                      "(bar color = terminal classification)", fontsize=10,
                      loc="left", color=INK)

    for ax in axes:
        ax.set_facecolor(SURF)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.suptitle("S3 — O-1 pilot: validated forward flow-map covering "
                 "(F-T8-S3; pilot, NOT an O-1 discharge)", fontsize=11,
                 color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(os.path.join(S3DIR, "s3_fig.png"), facecolor=SURF,
                bbox_inches="tight")
    print("figure written", flush=True)


# ----------------------------------------------------------------------
def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "produce"
    t0 = time.time()
    if phase in ("g1", "produce") and not os.path.exists(
            os.path.join(S3DIR, "_s3_g1.json")):
        print("== phase g1 ==", flush=True)
        phase_g1()
    if phase in ("bench", "produce") and not os.path.exists(
            os.path.join(S3DIR, "_s3_bench.json")):
        print("== phase bench ==", flush=True)
        phase_bench()
    if phase in ("cover", "produce"):
        print("== phase cover ==", flush=True)
        phase_cover()
    if phase in ("ladder", "produce"):
        print("== phase ladder ==", flush=True)
        phase_ladder()
    if phase in ("truegrowth", "produce") and not os.path.exists(
            os.path.join(S3DIR, "_s3_truegrowth.json")):
        print("== phase truegrowth ==", flush=True)
        phase_truegrowth()
    if phase in ("gates", "produce"):
        print("== phase gates ==", flush=True)
        phase_gates()
    if phase in ("fig", "produce"):
        print("== phase fig ==", flush=True)
        phase_fig()
    print("total elapsed %.1f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
