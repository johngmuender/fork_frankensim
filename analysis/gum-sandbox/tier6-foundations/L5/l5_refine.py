#!/usr/bin/env python3
"""
L5 refinement ladder (G2 retry).

The 2-D run measured z-mode leakage = 0.0 exactly: the state remains
phi(x,z,t) = f(x,t) sin(pi z / L) e^{-i pi^2 t/2}
for all t (initial data is the pure ground sine mode; the free
Hamiltonian is diagonal in the sine basis).  The guidance law then
reduces EXACTLY to

  a(x,t) = Im(f'/f)                (convective x-velocity)
  b(x,t) = Re(f'/f) = dx ln|f|     (transverse spin z-velocity)
  transverse:  dx/dt = a - pi cot(pi z),  dz/dt = b
  axial:       dx/dt = a,                 dz/dt = 0

so the z-dependence is ANALYTIC (no z-grid, no near-wall interpolation
error).  f(x,t) is propagated exactly in 1-D Fourier space.  This engine
is validated against the full 2-D engine at matched (Nx, dt), then used
for a refinement ladder: double Nx, successively halve dt.
"""
import json
import os

import numpy as np
from scipy import fft as sfft
from scipy import stats

import l5_arrival as m

OUTDIR = os.path.dirname(os.path.abspath(__file__))
VCLAMP = m.VCLAMP


def build_f(Nx):
    dx = (m.XMAX - m.XMIN) / Nx
    xg = m.XMIN + dx * np.arange(Nx)
    f = np.exp(-(xg - m.X0) ** 2 / (4.0 * m.SIGMA_X ** 2)).astype(complex)
    f[np.abs(xg - m.X0) > m.TRUNC] = 0.0
    f *= np.exp(1j * m.K0 * xg)
    f /= np.sqrt(np.sum(np.abs(f) ** 2) * dx)   # 1-D norm (z-factor separately normalized)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    return xg, dx, f, kx


def lin_interp(V, px, xg0, dx, Nx):
    g = np.clip((px - xg0) / dx, 0.0, Nx - 1.000001)
    i = g.astype(np.int64)
    fr = g - i
    return V[i] * (1 - fr) + V[i + 1] * fr


def run_fast(Nx, dt, P0, t_final=m.T_FINAL, rho_eps_rel=m.RHO_EPS_REL):
    xg, dx, f0, kx = build_f(Nx)
    F = sfft.fft(f0)
    Eh = np.exp(-0.5j * kx ** 2 * (dt / 2.0))
    nsteps = int(round(t_final / dt))
    N = P0.shape[0]

    def fields(F):
        f = sfft.ifft(F)
        fp = sfft.ifft(1j * kx * F)
        rho = np.abs(f) ** 2
        den = np.maximum(rho, rho_eps_rel * rho.max())
        w = np.conj(f) * fp
        a = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
        b = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
        return a, b

    ens = {"T": P0.copy(), "A": P0.copy()}
    active = {"T": np.ones(N, bool), "A": np.ones(N, bool)}
    tau = {s: {d: np.full(N, np.nan) for d in ("near", "far")} for s in "TA"}
    dets = {"near": m.D_NEAR, "far": m.D_FAR}

    def vel(ab, s, px, pz):
        a, b = ab
        va = lin_interp(a, px, xg[0], dx, Nx)
        if s == "A":
            return va, np.zeros_like(va)
        vb = lin_interp(b, px, xg[0], dx, Nx)
        vx = np.clip(va - np.pi / np.tan(np.pi * pz / m.L), -VCLAMP, VCLAMP)
        return vx, vb

    ab = fields(F)
    for n in range(nsteps):
        t = n * dt
        F = F * Eh
        ab_h = fields(F)
        F = F * Eh
        ab_n = fields(F)
        for s in "TA":
            P = ens[s]
            act = active[s]
            if not act.any():
                continue
            idx = np.where(act)[0]
            px, pz = P[idx, 0], P[idx, 1]
            k1x, k1z = vel(ab, s, px, pz)
            k2x, k2z = vel(ab_h, s, px + 0.5 * dt * k1x,
                           np.clip(pz + 0.5 * dt * k1z, 1e-9, m.L - 1e-9))
            k3x, k3z = vel(ab_h, s, px + 0.5 * dt * k2x,
                           np.clip(pz + 0.5 * dt * k2z, 1e-9, m.L - 1e-9))
            k4x, k4z = vel(ab_n, s, px + dt * k3x,
                           np.clip(pz + dt * k3z, 1e-9, m.L - 1e-9))
            nx = px + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz = np.clip(pz + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
                         1e-9, m.L - 1e-9)
            for dname, dpos in dets.items():
                ta = tau[s][dname]
                cross = np.isnan(ta[idx]) & (px < dpos) & (nx >= dpos)
                if cross.any():
                    ci = idx[cross]
                    frac = (dpos - px[cross]) / (nx[cross] - px[cross])
                    ta[ci] = t + dt * frac
            P[idx, 0] = nx
            P[idx, 1] = nz
            esc = (nx < m.XMIN + 0.5) | (nx > m.XMAX - 0.5)
            if esc.any():
                active[s][idx[esc]] = False
        ab = ab_n
    return dict(tau=tau, active=active)


def summarize(run, label):
    out = {"label": label}
    for s in "TA":
        for d in ("near", "far"):
            a = run["tau"][s][d]
            fin = a[np.isfinite(a)]
            out["n_%s_%s" % (s, d)] = int(fin.size)
            out["max_%s_%s" % (s, d)] = float(fin.max()) if fin.size else None
    tT = run["tau"]["T"]["near"]
    fin = np.isfinite(tT)
    srt = np.sort(tT[fin])
    out["T_near_top5"] = [float(v) for v in srt[-5:]]
    out["tau_1249"] = float(tT[1249]) if np.isfinite(tT[1249]) else None
    out["tau_653"] = float(tT[653]) if np.isfinite(tT[653]) else None
    return out


def main():
    rng = np.random.default_rng(m.SEED)
    xg, dx, _, _ = build_f(m.NX)
    P0 = m.sample_ensemble(rng, xg, dx)     # identical ensemble (rng-only)

    ladder = [
        ("Nx1024_dt2e-3", 1024, 2.0e-3),    # validation vs 2-D engine
        ("Nx1024_dt1e-3", 1024, 1.0e-3),    # validation vs 2-D refined run
        ("Nx2048_dt1e-3", 2048, 1.0e-3),
        ("Nx2048_dt5e-4", 2048, 5.0e-4),
        ("Nx2048_dt2.5e-4", 2048, 2.5e-4),
    ]
    runs = {}
    for label, Nx, dt in ladder:
        print("run", label, flush=True)
        runs[label] = run_fast(Nx, dt, P0)
        print(json.dumps(summarize(runs[label], label)), flush=True)

    # rho-floor sensitivity for the edge outlier (finest run, looser floor)
    print("run eps-sensitivity (Nx2048_dt2.5e-4, eps=1e-10)", flush=True)
    eps_run = run_fast(2048, 2.5e-4, P0, rho_eps_rel=1e-10)
    print(json.dumps(summarize(eps_run, "eps1e-10")), flush=True)

    # cross-engine validation vs stored 2-D results
    raw = np.load(os.path.join(OUTDIR, "L5_raw.npz"))
    val = {}
    for label, key in (("Nx1024_dt2e-3", "tau_T_near"),
                       ("Nx1024_dt1e-3", "tau_T_near_refined")):
        t2d = raw[key]
        t1d = runs[label]["tau"]["T"]["near"]
        both = np.isfinite(t2d) & np.isfinite(t1d)
        diff = np.abs(t2d[both] - t1d[both])
        val[label] = dict(n_common=int(both.sum()),
                          median_abs_diff=float(np.median(diff)),
                          p99_abs_diff=float(np.percentile(diff, 99)),
                          max_abs_diff=float(diff.max()),
                          ks_p=float(stats.ks_2samp(t2d[np.isfinite(t2d)],
                                                    t1d[np.isfinite(t1d)]).pvalue))
    print("cross-engine validation:", json.dumps(val), flush=True)

    # ---- final gate numbers from the two finest runs ----
    fine = runs["Nx2048_dt2.5e-4"]
    half = runs["Nx2048_dt5e-4"]

    def get(run, s, d):
        a = run["tau"][s][d]
        return a[np.isfinite(a)]

    tT = get(fine, "T", "near")
    tA = get(fine, "A", "near")
    fT = get(fine, "T", "far")
    fA = get(fine, "A", "far")

    tau_max_fine = float(tT.max())
    tau_max_half = float(get(half, "T", "near").max())
    rel = abs(tau_max_fine - tau_max_half) / tau_max_fine
    iqr = float(np.subtract(*np.percentile(tT, [75, 25])))

    # bulk (excluding the flagged edge particle 1249)
    mask = np.ones(m.N_PART, bool)
    mask[1249] = False
    tT_all = fine["tau"]["T"]["near"]
    tT_bulk = tT_all[mask & np.isfinite(tT_all)]
    tTh_all = half["tau"]["T"]["near"]
    tTh_bulk = tTh_all[mask & np.isfinite(tTh_all)]
    bulk_max_fine = float(tT_bulk.max())
    bulk_max_half = float(tTh_bulk.max())
    bulk_rel = abs(bulk_max_fine - bulk_max_half) / bulk_max_fine

    axial_max = float(tA.max())
    axial_max_half = float(get(half, "A", "near").max())
    frac_axial_late_bulk = float(np.mean(tA > bulk_max_fine))
    frac_axial_late_full = float(np.mean(tA > tau_max_fine))
    ks = stats.ks_2samp(fT, fA)
    rngb = np.random.default_rng(777)
    ratio = float(tT.mean() / tA.mean())
    boots = np.array([rngb.choice(tT, tT.size).mean()
                      / rngb.choice(tA, tA.size).mean() for _ in range(2000)])

    out = dict(
        ladder={lab: summarize(runs[lab], lab) for lab, _, _ in ladder},
        eps_sensitivity=summarize(eps_run, "eps1e-10"),
        validation=val,
        fine=dict(tau_max=tau_max_fine, tau_max_half=tau_max_half,
                  tau_max_rel_change=rel, iqr=iqr,
                  bulk_max=bulk_max_fine, bulk_max_half=bulk_max_half,
                  bulk_rel_change=bulk_rel,
                  axial_max=axial_max, axial_max_half=axial_max_half,
                  frac_axial_beyond_bulk_cutoff=frac_axial_late_bulk,
                  frac_axial_beyond_full_max=frac_axial_late_full,
                  n_T_near=int(tT.size), n_A_near=int(tA.size),
                  n_T_far=int(fT.size), n_A_far=int(fA.size),
                  ks_far_stat=float(ks.statistic), ks_far_p=float(ks.pvalue),
                  mean_T=float(tT.mean()), mean_A=float(tA.mean()),
                  ratio=ratio, ratio_err=float(boots.std(ddof=1))),
    )
    with open(os.path.join(OUTDIR, "L5_refine.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    np.savez_compressed(
        os.path.join(OUTDIR, "L5_refine_raw.npz"),
        tau_T_near=fine["tau"]["T"]["near"], tau_A_near=fine["tau"]["A"]["near"],
        tau_T_far=fine["tau"]["T"]["far"], tau_A_far=fine["tau"]["A"]["far"],
        tau_T_near_half=half["tau"]["T"]["near"],
        tau_A_near_half=half["tau"]["A"]["near"])
    print("FINE:", json.dumps(out["fine"], indent=2))


if __name__ == "__main__":
    main()
