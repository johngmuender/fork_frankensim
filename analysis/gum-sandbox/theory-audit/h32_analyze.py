#!/usr/bin/env python3
"""Phase H3.2 — analysis layer over h32_results.json (produced by
h32_solve/src/bin/h32_selection.rs).

Computes, deterministically, from the raw per-cycle trajectories:
  * the kappa_paper = 1 crossing (the corpus's claimed stopping regime):
    interpolated cycle/iteration, and whether ANYTHING in the dynamics
    distinguishes it — local dkappa/dcycle, dR/diter, halo growth rate,
    I-growth rate in windows before/at/after the crossing, plus a plateau
    test (is there any local slowdown of |dkappa| near the crossing?);
  * the same battery at the corpus benchmark kappa = 0.802;
  * the endpoint: linear-in-kappa rate fit dkappa/dcyc = -lam*(kappa-kinf)
    (exponential relaxation), giving the flow's attractor kinf, fitted on
    the full trajectory and on the tail half (stability check);
  * the control runs: endpoint, drift, whether the projection binds
    (dI_cut per cycle), and whether any constrained flow stabilises near
    the corpus benchmark.

Appends an "analysis" block into h32_results.json (idempotent: recomputed
from the raw trajectories each invocation) and prints a summary.

Usage: python3 h32_analyze.py [path/to/h32_results.json]
"""

import json
import sys

SQRT2INV = 2.0 ** -0.5


def traj_arrays(run):
    cyc = run["cycles"]
    return {
        "cyc": [c["cyc"] for c in cyc],
        # kappa of the clock-consistent state at the end of each cycle
        "kappa": [c["kappa_paper_next"] for c in cyc],
        # kappa of the block (L held during descent)
        "kappa_blk": [c["kappa_paper"] for c in cyc],
        "I": [c["I"] for c in cyc],
        "I_pre": [c["I_preproj"] for c in cyc],
        "E": [c["Estat"] for c in cyc],
        "L": [c["L"] for c in cyc],
        "Lnext": [c["L_next"] for c in cyc],
        "halo": [c["halo_pre"] for c in cyc],
        "if150": [c["ifrac_r150"] for c in cyc],
        "if125": [c["ifrac_r125"] for c in cyc],
        "if100": [c["ifrac_r100"] for c in cyc],
        "rrms": [c["r_rms_over_rstar"] for c in cyc],
        "dRdit": [(c["R_end"] - c["R0"]) / max(c["iters"], 1) for c in cyc],
        "erot": [c["Erot_over_Etot"] for c in cyc],
        "arrests": [c["arrests"] for c in cyc],
        "cpaper": [c["c_paper_next"] for c in cyc],
    }


def crossing(t, level):
    """First downward crossing of kappa through `level`: interpolated cycle."""
    k = t["kappa"]
    for n in range(1, len(k)):
        if k[n - 1] > level >= k[n]:
            frac = (k[n - 1] - level) / (k[n - 1] - k[n])
            return (n - 1) + frac + 1.0  # cycles are 1-based
    return None


def window_stats(t, c0, half):
    """Mean dkappa/dcyc, dR/diter, dhalo/dcyc, dlnI/dcyc in cycles [c0-half, c0+half]."""
    k, h, i_arr, dr = t["kappa"], t["halo"], t["I"], t["dRdit"]
    lo = max(1, int(round(c0)) - half)
    hi = min(len(k) - 1, int(round(c0)) + half)
    if hi <= lo:
        return None
    dk = (k[hi] - k[lo]) / (hi - lo)
    dh = (h[hi] - h[lo]) / (hi - lo)
    import math

    dlni = (math.log(i_arr[hi]) - math.log(i_arr[lo])) / (hi - lo)
    drm = sum(dr[lo : hi + 1]) / (hi - lo + 1)
    return {"dkappa_dcyc": dk, "dR_diter": drm, "dhalo_dcyc": dh, "dlnI_dcyc": dlni,
            "window_cycles": [lo + 1, hi + 1]}


def plateau_test(t, c0, half):
    """Any local slowdown near c0?  min |dkappa/dcyc| in the window vs the
    median |dkappa/dcyc| in the surrounding +-3*half band."""
    k = t["kappa"]
    dk = [k[n] - k[n - 1] for n in range(1, len(k))]
    n0 = int(round(c0)) - 1
    lo, hi = max(0, n0 - half), min(len(dk), n0 + half)
    inner = [abs(v) for v in dk[lo:hi]]
    blo, bhi = max(0, n0 - 3 * half), min(len(dk), n0 + 3 * half)
    band = sorted(abs(v) for v in dk[blo:bhi])
    med = band[len(band) // 2] if band else float("nan")
    return {
        "min_abs_dkappa_in_window": min(inner) if inner else float("nan"),
        "median_abs_dkappa_band": med,
        "ratio_min_over_median": (min(inner) / med) if inner and med else float("nan"),
        "plateau": bool(inner) and med > 0 and min(inner) < 0.1 * med,
    }


def rate_fit(t, lo_frac=0.0):
    """Least-squares fit dkappa/dcyc = -lam*(kappa - kinf) on cycles beyond
    lo_frac of the trajectory.  Returns (lam, kinf, rms)."""
    k = t["kappa"]
    n0 = int(len(k) * lo_frac)
    xs, ys = [], []
    for n in range(max(1, n0), len(k)):
        xs.append(0.5 * (k[n] + k[n - 1]))
        ys.append(k[n] - k[n - 1])
    m = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    den = m * sxx - sx * sx
    a = (m * sxy - sx * sy) / den  # slope   = -lam
    b = (sy * sxx - sx * sxy) / den  # intercept = lam*kinf
    lam = -a
    kinf = b / lam if lam != 0 else float("nan")
    rms = (sum((a * x + b - y) ** 2 for x, y in zip(xs, ys)) / m) ** 0.5
    return {"lambda_per_cycle": lam, "kappa_infinity": kinf, "rms_residual": rms,
            "fit_cycles": [max(1, n0) + 1, len(k)]}


def endpoint(t, tail=20):
    k = t["kappa"]
    n = len(k)
    tailk = k[n - tail :]
    drift = (k[-1] - k[n - tail]) / tail
    return {
        "kappa_last": k[-1],
        "kappa_tail_mean": sum(tailk) / len(tailk),
        "dkappa_dcyc_tail": drift,
        "halo_last": t["halo"][-1],
        "ifrac_r150_last": t["if150"][-1],
        "ifrac_r125_last": t["if125"][-1],
        "ifrac_r100_last": t["if100"][-1],
        "r_rms_last": t["rrms"][-1],
        "I_last": t["I"][-1],
        "dI_cut_last": t["I_pre"][-1] - t["I"][-1],
        "dR_diter_last": t["dRdit"][-1],
        "c_paper_last": t["cpaper"][-1],
        "erot_frac_last": t["erot"][-1],
    }


def analyze_run(run, name):
    t = traj_arrays(run)
    out = {"name": name}
    # crossings of the two narrative levels
    for label, level in [("kappa_eq_1", 1.0), ("benchmark_0p802", 0.802),
                         ("deep_bps_0p7638", (7.0 / 12.0) ** 0.5),
                         ("saturation_0p7071", SQRT2INV)]:
        c = crossing(t, level)
        blk = {"level": level, "crossing_cycle": c}
        if c is not None:
            blk["before"] = window_stats(t, c - 12, 6)
            blk["at"] = window_stats(t, c, 6)
            blk["after"] = window_stats(t, c + 12, 6)
            blk["plateau_test"] = plateau_test(t, c, 8)
        out[label] = blk
    out["rate_fit_full"] = rate_fit(t, 0.0)
    out["rate_fit_tail_half"] = rate_fit(t, 0.5)
    out["endpoint"] = endpoint(t)
    out["max_erot_frac_dev_from_quarter"] = max(abs(v - 0.25) for v in t["erot"])
    out["total_arrests"] = sum(t["arrests"])
    return out


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "h32_results.json"
    with open(path) as fh:
        data = json.load(fh)

    runs = [
        ("free", data["free"]),
        ("control_halo_1p5rstar", data["control_halo_1p5rstar"]),
        ("control_support_1p25rstar", data["control_support_1p25rstar"]),
    ]
    analysis = {"script": "h32_analyze.py",
                "kappa_convention": "kappa_paper_next = 2 sqrt(pi) L_next / I (the clock-consistent state at each cycle end)"}
    for name, run in runs:
        analysis[name] = analyze_run(run, name)

    # cross-run verdict facts
    free = analysis["free"]
    ch = analysis["control_halo_1p5rstar"]
    cs = analysis["control_support_1p25rstar"]
    analysis["verdict_facts"] = {
        "rigid_rung_start_kappa": data["rigid_rung"]["kappa_paper"],
        "free_crosses_1": free["kappa_eq_1"]["crossing_cycle"] is not None,
        "free_kappa1_plateau": free["kappa_eq_1"].get("plateau_test", {}).get("plateau"),
        "free_crosses_benchmark": free["benchmark_0p802"]["crossing_cycle"] is not None,
        "free_benchmark_plateau": free["benchmark_0p802"].get("plateau_test", {}).get("plateau"),
        "free_attractor_kinf_full": free["rate_fit_full"]["kappa_infinity"],
        "free_attractor_kinf_tail": free["rate_fit_tail_half"]["kappa_infinity"],
        "saturation_locus": SQRT2INV,
        "free_kinf_minus_saturation_tail": free["rate_fit_tail_half"]["kappa_infinity"] - SQRT2INV,
        "free_kinf_minus_benchmark_tail": free["rate_fit_tail_half"]["kappa_infinity"] - 0.802,
        "control_halo_attractor_kinf_tail": ch["rate_fit_tail_half"]["kappa_infinity"],
        "control_support_attractor_kinf_tail": cs["rate_fit_tail_half"]["kappa_infinity"],
        "control_support_endpoint_kappa": cs["endpoint"]["kappa_tail_mean"],
        "control_support_tail_drift": cs["endpoint"]["dkappa_dcyc_tail"],
        "control_support_dI_cut_last": cs["endpoint"]["dI_cut_last"],
        "control_halo_endpoint_kappa": ch["endpoint"]["kappa_tail_mean"],
        "corpus_benchmark": [0.802, 0.018],
    }

    data["analysis"] = analysis
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1)

    vf = analysis["verdict_facts"]
    print("=== H3.2 analysis summary ===")
    for name in ("free", "control_halo_1p5rstar", "control_support_1p25rstar"):
        a = analysis[name]
        print(f"\n[{name}]")
        for lab in ("kappa_eq_1", "benchmark_0p802", "deep_bps_0p7638", "saturation_0p7071"):
            b = a[lab]
            cc = b["crossing_cycle"]
            msg = f"  cross {b['level']:.4f}: " + (f"cycle {cc:.1f}" if cc else "NOT crossed")
            if cc:
                pt = b["plateau_test"]
                msg += (f"; dk/dcyc before/at/after = {b['before']['dkappa_dcyc']:+.2e}/"
                        f"{b['at']['dkappa_dcyc']:+.2e}/{b['after']['dkappa_dcyc']:+.2e}"
                        f"; plateau={pt['plateau']} (min/med={pt['ratio_min_over_median']:.2f})")
            print(msg)
        rf, rt = a["rate_fit_full"], a["rate_fit_tail_half"]
        print(f"  rate fit dk/dcyc=-lam(k-kinf): full kinf={rf['kappa_infinity']:.4f} lam={rf['lambda_per_cycle']:.2e} rms={rf['rms_residual']:.1e}")
        print(f"                                 tail kinf={rt['kappa_infinity']:.4f} lam={rt['lambda_per_cycle']:.2e} rms={rt['rms_residual']:.1e}")
        e = a["endpoint"]
        print(f"  endpoint: kappa={e['kappa_last']:.5f} (tail mean {e['kappa_tail_mean']:.5f}, drift {e['dkappa_dcyc_tail']:+.2e}/cyc); halo={e['halo_last']:.4f} if(1.25R*)={e['ifrac_r125_last']:.4f} r_rms={e['r_rms_last']:.3f}R*; dI_cut={e['dI_cut_last']:+.3e}; dR/dit={e['dR_diter_last']:+.2e}")
        print(f"  max |Erot/E - 1/4| = {a['max_erot_frac_dev_from_quarter']:.2e}; arrests total = {a['total_arrests']}")
    print("\n=== verdict facts ===")
    for k, v in vf.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
