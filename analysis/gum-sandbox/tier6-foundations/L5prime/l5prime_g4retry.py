#!/usr/bin/env python3
"""
L5' G4 refinement retry (coverage, not physics tuning).

At T = 24 the far-field KS compares differently-truncated ensembles
(non-arrivals at d_far = 25: transverse 30/2000 vs axial 290/2000 -- the
transverse spin current advances the slow half of the guide, so the T = 24
horizon cuts the axial slow tail but not the transverse one).  Per campaign
rules this is a coverage artifact, so one refinement retry: extend the
far-field evolution to T_far = 40 for the IDENTICAL ensemble and re-time the
far crossings.  The box is extended to [-15, 225] at the SAME reference dx
(Nx = 8192) so the packet front (~ x0 + (k0 + 3 sigma_v) T ~ 159) stays
contained (periodic FFT box).  Extra detector planes d = 10..35 give the
KS-distance-vs-d trend for the diagnosis if the gate still fails.

Companion run at dt = 5e-4 checks stability of the retry KS numbers.
Outputs: L5prime_g4retry.json, L5prime_g4retry_raw.npz, and an updated G4
block (with final verdict) merged into L5prime_results.json.
"""
import json
import os

import numpy as np
from scipy import stats

import l5prime_run as m

OUTDIR = os.path.dirname(os.path.abspath(__file__))

T_FAR = 40.0
XMIN_R, XMAX_R = -15.0, 225.0          # length 240; Nx=8192 -> dx = ref dx
NX_R = 8192
DETS = {"near": m.D_NEAR, "d10": 10.0, "d15": 15.0, "d20": 20.0,
        "far": m.D_FAR, "d30": 30.0, "d35": 35.0}


def fin(a):
    return a[np.isfinite(a)]


def ks_table(run, t_horizon):
    out = {}
    for name, dpos in sorted(DETS.items(), key=lambda kv: kv[1]):
        tT = run["tau"]["T"][name]
        tA = run["tau"]["A"][name]
        aT = fin(tT); aT = aT[aT <= t_horizon]
        aA = fin(tA); aA = aA[aA <= t_horizon]
        ks = stats.ks_2samp(aT, aA)
        out[name] = dict(d=dpos, n_T=int(aT.size), n_A=int(aA.size),
                         nonarr_T=int(m.N_PART - aT.size),
                         nonarr_A=int(m.N_PART - aA.size),
                         ks_stat=float(ks.statistic), ks_p=float(ks.pvalue),
                         median_T=float(np.median(aT)),
                         median_A=float(np.median(aA)))
    return out


def main():
    P0 = m.sample_ensemble(m.SEED, "literal")

    print("retry run: Nx=%d box=[%g,%g] dt=2.5e-4 T=%g" %
          (NX_R, XMIN_R, XMAX_R, T_FAR), flush=True)
    run = m.run_1d(NX_R, 2.5e-4, P0, XMIN_R, XMAX_R, "literal",
                   t_final=T_FAR, dets=DETS)
    print("  elapsed %.1f s" % run["elapsed"], flush=True)

    print("companion run: dt=5e-4", flush=True)
    run_h = m.run_1d(NX_R, 5.0e-4, P0, XMIN_R, XMAX_R, "literal",
                     t_final=T_FAR, dets=DETS)
    print("  elapsed %.1f s" % run_h["elapsed"], flush=True)

    # consistency with the original T=24 window (same detectors, bigger box)
    tab24 = ks_table(run, 24.0)
    tab40 = ks_table(run, T_FAR)
    tab40_h = ks_table(run_h, T_FAR)

    # gate quantity: far detector, full retry horizon
    g4 = tab40["far"]
    g4_pass = g4["ks_p"] > 0.05

    # near-field invariants re-checked on the retry run (box/horizon sanity)
    tTn = fin(run["tau"]["T"]["near"])
    tau_max_retry = float(tTn[tTn <= m.T_NEAR].max())
    n_T_gap = int(np.sum(tTn > tau_max_retry))

    out = dict(
        config=dict(Nx=NX_R, dx=(XMAX_R - XMIN_R) / NX_R, dt=2.5e-4,
                    dt_companion=5e-4, box=[XMIN_R, XMAX_R], T_far=T_FAR,
                    detectors={k: v for k, v in DETS.items()},
                    N=m.N_PART, seed=m.SEED),
        ks_vs_d_T24=tab24,
        ks_vs_d_T40=tab40,
        ks_vs_d_T40_dt5e4=tab40_h,
        gate=dict(ks_stat=g4["ks_stat"], ks_p=g4["ks_p"],
                  n_far_T=g4["n_T"], n_far_A=g4["n_A"],
                  nonarrival_far_T=g4["nonarr_T"],
                  nonarrival_far_A=g4["nonarr_A"],
                  verdict="PASS" if g4_pass else "FAIL"),
        near_field_sanity=dict(tau_max=tau_max_retry,
                               transverse_arrivals_beyond=n_T_gap),
    )
    with open(os.path.join(OUTDIR, "L5prime_g4retry.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    np.savez_compressed(
        os.path.join(OUTDIR, "L5prime_g4retry_raw.npz"),
        **{"tau_%s_%s" % (s, d): run["tau"][s][d]
           for s in "TA" for d in DETS})

    # ---- merge final G4 into L5prime_results.json ----
    res = json.load(open(os.path.join(OUTDIR, "L5prime_results.json")))
    res["G4"] = dict(
        ks_stat_T24=res["G4"]["ks_stat"], ks_p_T24=res["G4"]["ks_p"],
        n_far_T_T24=res["G4"]["n_far_T"], n_far_A_T24=res["G4"]["n_far_A"],
        ks_p_2d_engine_T24=res["G4"]["ks_p_2d_engine"],
        coverage_T24=dict(nonarrival_far_T=30, nonarrival_far_A=290),
        retry=out["gate"],
        retry_ks_p_dt5e4=tab40_h["far"]["ks_p"],
        ks_stat_vs_d_T40={k: dict(d=v["d"], stat=v["ks_stat"], p=v["ks_p"],
                                  nonarr_T=v["nonarr_T"],
                                  nonarr_A=v["nonarr_A"])
                          for k, v in tab40.items()},
        verdict="PASS" if g4_pass else "FAIL",
    )
    res["L5_comparison"]["L5prime"]["ks_far_p"] = g4["ks_p"]
    with open(os.path.join(OUTDIR, "L5prime_results.json"), "w") as fh:
        json.dump(res, fh, indent=2)

    print(json.dumps(out["gate"], indent=2))
    print("KS vs d (T=40):")
    for k, v in sorted(tab40.items(), key=lambda kv: kv[1]["d"]):
        print("  d=%5.1f stat=%.4f p=%.4f nonarr T/A = %d/%d"
              % (v["d"], v["ks_stat"], v["ks_p"],
                 v["nonarr_T"], v["nonarr_A"]))


if __name__ == "__main__":
    main()
