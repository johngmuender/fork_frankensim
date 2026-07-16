#!/usr/bin/env python3
"""
J1 — App I.1's "p-grid over the (B.4) window", executed.

WITHIN-MODEL ONLY: this asks whether the corpus's printed p = 0.84 +- 0.03
solves the corpus's printed equations under a stated normalization
convention.  Nothing here validates physics/nature.

Corpus spec (corpus2/01-GUM-Omega-Paper-v3.0-ext.md):
  App B.4:  V = m~^2 (1-sigma) + c2 (1-sigma)^2, window c2 in [-m~^2/2, m~^2/6)
  App B.5:  tail mu^2 = m~^2/(2 a_psi)
  App I.1:  "Profile shooting/relaxation: virial <= 3e-4; INT b = 1.000;
             p-grid over the (B.4) window."
  Sec VII:  "Yukawa tail f -> p (1+mu r) e^{-mu r}/(mu r)^2 * mu,
             p = 0.84 +- 0.03 <r1> [CAL]"  (literal reading: p = A/mu).

The Tier-2b pilot (tier2-closure/radial_solve.py, certified: all I.1 gates
PASS) solved c2 = 0 only.  This script adds the c2 term to the SAME
machinery (midpoint energy, graded grid, flow + exact tridiagonal Newton,
eps-dial to 0.05) and scans c2 over the printed window, extracting the tail
amplitude under BOTH candidate conventions identified by I3:
  raw/naive        p = A e^{-mu R*} / (2 pi)   (pilot's nearest, 0.805 at c2=0)
  edge-referenced  p = A e^{-mu R*} / 6        (I3's joint-solve pick, 0.843)
plus /(15pi/8), the literal A/mu, and f_fit(R*) for the record.
R* = R*(c2) is the eps = 0 compacton radius from the Bogomolny saturation
quadrature R*^3 = 3 INT_0^pi sin^2 f / sqrt(V(f)) df (2^{5/6} at c2 = 0).

Deterministic, no RNG.  Runtime ~1 min.  Writes j1_results.json
STAGE-FLUSHED (partial JSON after every grid point) and prints a PASS/FAIL
gate table.
"""

import json
import time

import numpy as np

# ----------------------------------------------------------------------
# constants / conventions (identical to the certified radial_solve.py)
# ----------------------------------------------------------------------
M = 1.0                      # potential mass m~
A6 = 1.0
A0 = 1.0
EPS_TARGET = 0.05
EPS_TOL = 5.0e-5
RMAX = 6.0
N_BASE = 4000                # base resolution (pilot's N)
C2_LO = -0.5 * M * M         # B.4 window, closed lower edge
C2_HI = M * M / 6.0          # B.4 window, OPEN upper edge

# the B.4 window grid: 9 points, upper point just inside the open edge
C2_GRID = [-0.500, -0.400, -0.300, -0.200, -0.100, 0.000, 0.060, 0.120, 0.160]
SENTINELS_4N = [-0.500, 0.000, 0.160]   # extra 4N resolution at these

P_TARGET = (0.84, 0.03)     # corpus p band
B_TARGET = (42.0, 6.0)      # corpus b_eff band
# best available converged C_d: G5b relaxed cap-1200 (i3_results.json inputs)
CD_G5B = (426075768.6944963, 82599726.60227346)     # value, sigma (19%)
# superseded cap-400 value, printed only for the caveat
CD_G4RELAX = (7865282.617028015, 3256457.565093463)

OUTDIR = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit"

RESULTS = {"spec": {
    "window": [C2_LO, "m^2/6 (open)"], "c2_grid": C2_GRID,
    "eps_target": EPS_TARGET, "Rmax": RMAX, "N_base": N_BASE,
    "p_target": P_TARGET, "b_target": B_TARGET,
    "Cd_G5b_relaxed_1200": CD_G5B,
}, "points": {}}


def flush():
    with open(OUTDIR + "/j1_results.json", "w") as fh:
        json.dump(RESULTS, fh, indent=1)


# ----------------------------------------------------------------------
# B.4 potential and the eps = 0 compacton geometry
# ----------------------------------------------------------------------
def Vpot(f, c2):
    """V(f) = m^2 (1-cos f) + c2 (1-cos f)^2  (sigma = cos f)."""
    one_minus = 1.0 - np.cos(f)
    return M * M * one_minus + c2 * one_minus ** 2


def compacton(c2, nquad=2000):
    """R*(c2) and BPS bound from Gauss-Legendre quadrature.

    Saturation of E6 + E0 (a6 = 1): sin^2 f f'/r = -sqrt(V) r integrates to
    R*^3 = 3 INT_0^pi sin^2 f / sqrt(V(f)) df; bound = 2 INT sin^2 f sqrt(V) df.
    Positivity in-window: V = 2 s^2 (m^2 + 2 c2 s^2) >= 0 for c2 >= -m^2/2.
    """
    x, w = np.polynomial.legendre.leggauss(nquad)
    f = 0.5 * np.pi * (x + 1.0)
    w = 0.5 * np.pi * w
    V = Vpot(f, c2)
    s2f = np.sin(f) ** 2
    I1 = float(np.sum(w * s2f / np.sqrt(np.maximum(V, 1e-300))))
    bound = 2.0 * float(np.sum(w * s2f * np.sqrt(np.maximum(V, 0.0))))
    return (3.0 * I1) ** (1.0 / 3.0), bound


# ----------------------------------------------------------------------
# grid / profile init (as certified, with c2-dependent edge location)
# ----------------------------------------------------------------------
def make_grid(N, rmax, rstar):
    s = np.linspace(0.0, 1.0, 200001)
    rr = s * rmax
    w = (0.55
         + 2.0 * np.exp(-(rr / 0.5) ** 2)
         + 3.0 * np.exp(-((rr - rstar) / 0.45) ** 2)
         + 0.8 * np.exp(-((rr - (rstar + 1.2)) / 0.9) ** 2))
    cw = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(rr))])
    cw /= cw[-1]
    r = np.interp(np.linspace(0.0, 1.0, N + 1), cw, rr)
    r[0] = 0.0
    r[-1] = rmax
    return r


def init_profile(r, rstar):
    u = r / rstar
    x = u / (1.0 + u ** 8) ** 0.125
    f = 2.0 * np.arccos(np.clip(x, 0.0, 1.0))
    f[0] = np.pi
    f[-1] = 0.0
    return f


# ----------------------------------------------------------------------
# discretized energy with the c2 term; exact gradient + tridiag Hessian
# ----------------------------------------------------------------------
def assemble(r, f, t, c2, need_hess=True):
    h = np.diff(r)
    rm = 0.5 * (r[1:] + r[:-1])
    rm2 = rm * rm
    fm = 0.5 * (f[1:] + f[:-1])
    d = np.diff(f) / h
    d2 = d * d
    s = np.sin(fm)
    c = np.cos(fm)
    s2 = s * s
    s3 = s2 * s
    s4 = s2 * s2
    sin2f = 2.0 * s * c
    cos2f = c * c - s2
    om = 1.0 - c                       # 1 - cos f

    L2 = d2 * rm2 + 2.0 * s2
    L4 = s2 * (2.0 * d2 + s2 / rm2)
    L6 = s4 * d2 / rm2
    L0 = (M * M * om + c2 * om * om) * rm2
    E2 = float(np.sum(h * L2))
    E4 = float(np.sum(h * L4))
    E6 = float(np.sum(h * L6))
    E0 = float(np.sum(h * L0))
    E = t * (E2 + E4) + A6 * E6 + A0 * E0

    L_d = 2.0 * t * d * rm2 + 4.0 * t * s2 * d + 2.0 * s4 * d / rm2
    L_f = (2.0 * t * sin2f * (1.0 + d2)
           + 4.0 * (t + d2) * s3 * c / rm2
           + (M * M * s + 2.0 * c2 * om * s) * rm2)

    gL = 0.5 * h * L_f - L_d
    gR = 0.5 * h * L_f + L_d
    g = gR[:-1] + gL[1:]

    if not need_hess:
        return E, (E2, E4, E6, E0), g, None, None

    L_dd = 2.0 * t * rm2 + 4.0 * t * s2 + 2.0 * s4 / rm2
    L_fd = 4.0 * t * d * sin2f + 8.0 * s3 * c * d / rm2
    L_ff = (4.0 * t * cos2f * (1.0 + d2)
            + 4.0 * (t + d2) * (3.0 * s2 * c * c - s4) / rm2
            + (M * M * c + 2.0 * c2 * (s2 + om * c)) * rm2)

    HLL = 0.25 * h * L_ff - L_fd + L_dd / h
    HRR = 0.25 * h * L_ff + L_fd + L_dd / h
    HLR = 0.25 * h * L_ff - L_dd / h
    diag = HRR[:-1] + HLL[1:]
    off = HLR[1:-1]
    return E, (E2, E4, E6, E0), g, diag, off


def thomas(diag, off, b):
    n = diag.size
    dd = diag.tolist()
    ee = off.tolist()
    bb = b.tolist()
    for i in range(1, n):
        mfac = ee[i - 1] / dd[i - 1]
        dd[i] -= mfac * ee[i - 1]
        bb[i] -= mfac * bb[i - 1]
    x = [0.0] * n
    x[-1] = bb[-1] / dd[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (bb[i] - ee[i] * x[i + 1]) / dd[i]
    return np.asarray(x)


def flow(r, f, t, c2, gmax_target=1e-2, maxsteps=6000):
    E, _, g, diag, _ = assemble(r, f, t, c2)
    dt = 0.5
    for _ in range(maxsteps):
        if np.max(np.abs(g)) < gmax_target:
            break
        prec = np.maximum(diag, 1e-3 * np.max(diag))
        fn = f.copy()
        fn[1:-1] = f[1:-1] - dt * g / prec
        np.clip(fn, 0.0, np.pi, out=fn)
        En, _, gn, diagn, _ = assemble(r, fn, t, c2)
        if En < E:
            f, E, g, diag = fn, En, gn, diagn
            dt = min(dt * 1.1, 5.0)
        else:
            dt *= 0.5
            if dt < 1e-9:
                break
    return f


def newton(r, f, t, c2, gtol=1e-13, maxit=300):
    E, sect, g, diag, off = assemble(r, f, t, c2)
    lam = 1e-3
    it = 0
    while it < maxit:
        gmax = float(np.max(np.abs(g)))
        if gmax < gtol:
            break
        scale = np.abs(diag) + 1e-30
        accepted = False
        while lam < 1e14:
            step = thomas(diag + lam * scale, off, -g)
            if np.all(np.isfinite(step)):
                fn = f.copy()
                fn[1:-1] += step
                En, sectn, gn, diagn, offn = assemble(r, fn, t, c2)
                gmaxn = float(np.max(np.abs(gn)))
                if np.isfinite(En) and gmaxn < gmax:
                    f, E, sect, g, diag, off = fn, En, sectn, gn, diagn, offn
                    lam = max(lam * 0.25, 1e-14)
                    accepted = True
                    break
            lam *= 10.0
        if not accepted:
            break
        it += 1
    return f, E, sect, float(np.max(np.abs(g))), it


def solve_profile(r, f0, t, c2):
    f = flow(r, f0, t, c2, gmax_target=1e-2)
    return newton(r, f, t, c2)


def degree(r, f):
    fm = 0.5 * (f[1:] + f[:-1])
    return float(-(2.0 / np.pi) * np.sum(np.sin(fm) ** 2 * np.diff(f)))


def tail_fit(r, f, mu_lo, mu_hi, rstar):
    mask = (f > 1e-6) & (f < 1e-2) & (r > rstar)
    rw = r[mask]
    y = np.log(f[mask])

    def misfit(mu):
        x = mu * rw
        model = -x + np.log(1.0 / x + 1.0 / (x * x))
        lgA = float(np.mean(y - model))
        res = y - model - lgA
        return float(np.sqrt(np.mean(res ** 2))), lgA

    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    a, b = mu_lo, mu_hi
    cc = b - invphi * (b - a)
    dd = a + invphi * (b - a)
    fc, _ = misfit(cc)
    fd, _ = misfit(dd)
    for _ in range(120):
        if fc < fd:
            b, dd, fd = dd, cc, fc
            cc = b - invphi * (b - a)
            fc, _ = misfit(cc)
        else:
            a, cc, fc = cc, dd, fd
            dd = a + invphi * (b - a)
            fd, _ = misfit(dd)
    mu = 0.5 * (a + b)
    rms, lgA = misfit(mu)
    return mu, float(np.exp(lgA)), rms, int(mask.sum())


# ----------------------------------------------------------------------
# per-point pipeline
# ----------------------------------------------------------------------
def conventions(A, mu, rstar):
    env = np.exp(-mu * rstar)
    x = mu * rstar
    return {
        "p_naive_2pi": float(A * env / (2.0 * np.pi)),
        "p_edge_6": float(A * env / 6.0),
        "p_edge_15pi8": float(A * env / (15.0 * np.pi / 8.0)),
        "p_literal_A_over_mu": float(A / mu),
        "f_fit_at_Rstar": float(A * env * (1.0 / x + 1.0 / x ** 2)),
    }


def run_point(c2, N, t0, f_init_from=None):
    """Solve the eps-dial at resolution N for this c2; return dict + (t, r, f)."""
    rstar, bps_bound = compacton(c2, 2000)
    rstar_hi, _ = compacton(c2, 4000)
    r = make_grid(N, RMAX, rstar)
    if f_init_from is None:
        f = init_profile(r, rstar)
    else:
        r_prev, f_prev = f_init_from
        f = np.interp(r, r_prev, f_prev)
        f[0] = np.pi
        f[-1] = 0.0
    t = t0
    for _ in range(20):
        f, E, sect, gmax, nit = solve_profile(r, f, t, c2)
        E2, E4, E6, E0 = sect
        ratio = t * (E2 + E4) / (A6 * E6 + A0 * E0)
        if abs(ratio - EPS_TARGET) < EPS_TOL:
            break
        t *= EPS_TARGET / ratio
    K = degree(r, f)
    vir_rel = abs(t * E2 - t * E4 - 3.0 * A6 * E6 + 3.0 * A0 * E0) / E
    mu_true = M / np.sqrt(2.0 * t)      # B.5: mu^2 = m~^2/(2 a_psi), c2-blind
    mu_fit, A_fit, rms, npts = tail_fit(r, f, 0.3 * mu_true, 2.5 * mu_true, rstar)
    i_chk = np.searchsorted(r, 0.97 * RMAX)
    f_far = float(np.max(np.abs(f[i_chk:])))
    D = {
        "c2": c2, "N": N, "t": float(t), "Rstar": float(rstar),
        "Rstar_quad_delta": float(abs(rstar_hi - rstar)),
        "BPS_bound": float(bps_bound),
        "E2": E2, "E4": E4, "E6": E6, "E0": E0, "E_total": float(E),
        "eps_ratio": float(ratio), "degree_K": float(K),
        "virial_rel": float(vir_rel), "grad_max": float(gmax),
        "mu_fit": float(mu_fit), "mu_derived": float(mu_true),
        "mu_rel_err": float(abs(mu_fit / mu_true - 1.0)),
        "A_fit": float(A_fit), "tail_rms": float(rms), "tail_npts": npts,
        "BPS_deficit": float((E6 + E0) / bps_bound - 1.0),
        "f_far": f_far,
        "monotone": bool(np.all(np.diff(f) <= 1e-12)),
    }
    D.update(conventions(A_fit, mu_fit, rstar))
    return D, t, r, f


def gates(D):
    return {
        "virial<=3e-4": D["virial_rel"] <= 3e-4,
        "degree=1(1e-3)": abs(D["degree_K"] - 1.0) <= 1e-3,
        "eps=0.05(1e-3)": abs(D["eps_ratio"] - EPS_TARGET) <= 1e-3,
        "muB.5(1e-4)": D["mu_rel_err"] <= 1e-4,
        "tail_rms<=1e-4": D["tail_rms"] <= 1e-4,
        "contained": D["f_far"] <= 1e-8,
        "monotone": D["monotone"],
    }


def main():
    t_wall = time.time()
    # walk outward from c2 = 0 in both directions for warm starts
    order = sorted(C2_GRID, key=lambda c: (0 if c >= 0 else 1, abs(c)))
    warm = {}    # c2 -> (t, r, f) at N_BASE
    for c2 in order:
        # nearest already-solved point as warm start
        t0, init = 0.008, None
        if warm:
            cnear = min(warm, key=lambda c: abs(c - c2))
            tprev, rprev, fprev = warm[cnear]
            t0, init = tprev, (rprev, fprev)
        D1, t1, r1, f1 = run_point(c2, N_BASE, t0, init)
        warm[c2] = (t1, r1, f1)
        D2, t2, r2, f2 = run_point(c2, 2 * N_BASE, t1, (r1, f1))
        entry = {"N": D1, "2N": D2, "gates_N": gates(D1), "gates_2N": gates(D2)}
        if c2 in SENTINELS_4N:
            D4, _, _, _ = run_point(c2, 4 * N_BASE, t2, (r2, f2))
            entry["4N"] = D4
            entry["gates_4N"] = gates(D4)
        RESULTS["points"][f"{c2:+.3f}"] = entry
        flush()
        print(f"[c2 = {c2:+.3f}]  t = {D2['t']:.8f}  mu = {D2['mu_fit']:.5f}  "
              f"R* = {D2['Rstar']:.6f}  A = {D2['A_fit']:.5e}  "
              f"p/2pi = {D2['p_naive_2pi']:.4f}  p/6 = {D2['p_edge_6']:.4f}  "
              f"({time.time()-t_wall:.1f}s)")

    # ------------------------------------------------------------------
    # analysis: where does 0.84 +- 0.03 land under each convention
    # ------------------------------------------------------------------
    cs = sorted(C2_GRID)
    keyorder = [f"{c:+.3f}" for c in cs]
    pt2 = {k: RESULTS["points"][k]["2N"] for k in keyorder}
    pmid, psig = P_TARGET
    summary = {}
    for conv in ["p_naive_2pi", "p_edge_6", "p_edge_15pi8",
                 "p_literal_A_over_mu", "f_fit_at_Rstar"]:
        vals = np.array([pt2[k][conv] for k in keyorder])
        valsN = np.array([RESULTS["points"][k]["N"][conv] for k in keyorder])
        inband = (vals >= pmid - psig) & (vals <= pmid + psig)
        # linear-interpolated window locations where the curve crosses 0.84
        crossings = []
        for i in range(len(cs) - 1):
            a, b = vals[i] - pmid, vals[i + 1] - pmid
            if a == 0.0 or a * b < 0.0:
                crossings.append(cs[i] + (cs[i + 1] - cs[i]) * a / (a - b))
        summary[conv] = {
            "values_2N": [float(v) for v in vals],
            "min": float(vals.min()), "max": float(vals.max()),
            "mean_over_window": float(vals.mean()),
            "value_at_c2_0": float(pt2["+0.000"][conv]),
            "n_inband": int(inband.sum()), "n_grid": len(cs),
            "band_hit_anywhere": bool(inband.any()),
            "crossings_of_0.84": [float(c) for c in crossings],
            "max_abs_N_vs_2N": float(np.max(np.abs(vals - valsN))),
            "spread_over_window_pct": float(
                100.0 * (vals.max() - vals.min()) / vals.mean())
            if conv != "p_literal_A_over_mu" else None,
        }
    RESULTS["p_grid_summary"] = summary
    RESULTS["c2_grid_sorted"] = cs
    flush()

    # ------------------------------------------------------------------
    # joint b_eff check at the c2 = 0 benchmark (two-knot ran there)
    # ------------------------------------------------------------------
    D0 = pt2["+0.000"]
    mu0, R0 = D0["mu_fit"], D0["Rstar"]
    env0 = np.exp(-mu0 * R0)
    bmid, bsig = B_TARGET
    joint = {}
    for cd_name, (cd, cd_sig) in [("G5b_relaxed_1200", CD_G5B),
                                  ("G4_relaxed_400_SUPERSEDED", CD_G4RELAX)]:
        rows = {}
        for p_name in ["p_naive_2pi", "p_edge_6", "p_edge_15pi8"]:
            p = D0[p_name]
            for fb_name, fb in [("2pi", 2.0 * np.pi), ("2", 2.0),
                                ("4pi^2", 4.0 * np.pi ** 2),
                                ("8pi^2", 8.0 * np.pi ** 2)]:
                b = cd * mu0 * env0 / (fb * p * p)
                sig_eff = np.sqrt(bsig ** 2 + (b * cd_sig / cd) ** 2)
                rows[f"{p_name}|Fb={fb_name}"] = {
                    "b_eff": float(b), "pull": float((b - bmid) / sig_eff)}
        joint[cd_name] = rows
    RESULTS["joint_b_eff_c2_0"] = {
        "formula": "b = C_d mu e^{-mu R*} / (F_b p^2)  [I3 gb=1 edge envelope]",
        "mu": mu0, "Rstar": R0, "rows": joint}
    flush()

    # ------------------------------------------------------------------
    # gate table
    # ------------------------------------------------------------------
    print("\n===== PASS/FAIL gate table =====")
    gate_names = list(gates(D0).keys())
    hdr = "c2      res  " + "  ".join(f"{g:>15s}" for g in gate_names)
    print(hdr)
    all_ok = True
    for k in keyorder:
        for res in ["N", "2N"] + (["4N"] if "4N" in RESULTS["points"][k] else []):
            gg = RESULTS["points"][k][f"gates_{res}"]
            all_ok &= all(gg.values())
            print(f"{k:7s} {res:>3s}  " + "  ".join(
                f"{'PASS' if gg[g] else 'FAIL':>15s}" for g in gate_names))
    print(f"\nALL SOLVER GATES: {'PASS' if all_ok else 'FAIL'}")
    RESULTS["all_solver_gates_pass"] = bool(all_ok)

    print("\n===== p over the B.4 window (2N values) =====")
    print("c2      " + "  ".join(f"{c:>12s}" for c in
                                 ["p/2pi", "p/6", "p/(15pi/8)", "mu", "R*"]))
    for k in keyorder:
        d = pt2[k]
        print(f"{k:7s} {d['p_naive_2pi']:12.5f} {d['p_edge_6']:12.5f} "
              f"{d['p_edge_15pi8']:12.5f} {d['mu_fit']:12.5f} "
              f"{d['Rstar']:12.6f}")

    print("\n===== verdicts vs p = 0.84 +- 0.03 =====")
    verd = {}
    for conv, tag in [("p_naive_2pi", "raw/naive  A e^-muR*/2pi"),
                      ("p_edge_6", "edge-ref.  A e^-muR*/6")]:
        s = summary[conv]
        hit = s["band_hit_anywhere"]
        verd[conv] = hit
        print(f"{tag:30s}: range [{s['min']:.4f}, {s['max']:.4f}]  "
              f"in-band {s['n_inband']}/{s['n_grid']} grid pts  "
              f"0.84 crossed in-window: {s['crossings_of_0.84']}  "
              f"-> {'CONFIRM (0.84 reachable in-window)' if hit else 'NOT RECOVERED in-window'}")
    RESULTS["verdicts_p"] = verd

    print("\n===== joint b_eff (c2 = 0, G5b relaxed cap-1200) =====")
    for key, row in joint["G5b_relaxed_1200"].items():
        inb = abs(row["pull"]) <= 1.5
        print(f"{key:32s} b = {row['b_eff']:10.2f}  pull = {row['pull']:+6.2f}"
              f"  {'in 1.5 sigma' if inb else 'OUT'}")
    flush()
    print(f"\ntotal runtime {time.time()-t_wall:.1f} s; wrote j1_results.json")


if __name__ == "__main__":
    main()
