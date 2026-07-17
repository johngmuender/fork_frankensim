#!/usr/bin/env python3
"""
L7a stage 2: resolution refinement of the preparation grid + final gate logic.

Stage-1 finding: on the pre-registered 7x8 (theta,phi) grid the main family
(k1=1.2, k2=2.4, co-located) shows ZERO backflow (3.9e-31) and Bohmian bin
probabilities quadratic at 3.9e-16 -- but a fine continuum scan shows a genuine
narrow backflow pocket (max 2.8e-4 near theta=0.772, phi=0.245) and
lambda_min(J(t)) = -4.9e-4 < 0. The registered grid simply misses the pocket:
a RESOLUTION failure, so per campaign rules we refine the preparation grid
(uniform dyadic refinements of the same registered family; physics untouched).

Also: G4 control redo with the favorable (smaller-separation) single
adjustment, and final verdicts.
"""
import importlib.util
import json
import sys

import numpy as np

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7a"
spec = importlib.util.spec_from_file_location("l7a", f"{BASE}/l7a_povm_exclusion.py")
l7a = importlib.util.module_from_spec(spec)
sys.modules["l7a"] = l7a
spec.loader.exec_module(l7a)

NT, TGRID, EDGE_IDX, NBIN = l7a.NT, l7a.TGRID, l7a.EDGE_IDX, l7a.NBIN


# ---------------- generalized family machinery (arbitrary theta/phi grids) ---
def grid_thph(n_th, n_ph):
    ths = np.linspace(0.0, np.pi / 2.0, n_th)
    phs = np.arange(n_ph) * 2.0 * np.pi / n_ph
    TH, PH = np.meshgrid(ths, phs, indexing="ij")
    return TH.ravel(), PH.ravel()


def family_N_g(th, ph, Iaa_R, Ibb_R, Iab_R, overlap):
    c2 = np.cos(th) ** 2
    s2 = np.sin(th) ** 2
    s2t = np.sin(2.0 * th)
    n2 = 1.0 + s2t * (np.cos(ph) * overlap.real - np.sin(ph) * overlap.imag)
    quad = (c2[:, None] * Iaa_R[None, :]
            + s2[:, None] * Ibb_R[None, :]
            + s2t[:, None] * (np.cos(ph)[:, None] * Iab_R.real[None, :]
                              - np.sin(ph)[:, None] * Iab_R.imag[None, :]))
    return quad / n2[:, None], n2


def design_g(th, ph, overlap):
    n2 = 1.0 + np.sin(2 * th) * (np.cos(ph) * overlap.real
                                 - np.sin(ph) * overlap.imag)
    B = np.column_stack([np.cos(th) ** 2, np.sin(th) ** 2,
                         np.sin(2 * th) * np.cos(ph),
                         np.sin(2 * th) * np.sin(ph)])
    return B / n2[:, None]


def analyze(th, ph, Iaa_R, Ibb_R, Iab_R, overlap):
    """Full pipeline on a given preparation grid. Returns dict + arrays."""
    Nmat, n2 = family_N_g(th, ph, Iaa_R, Ibb_R, Iab_R, overlap)
    CDF = np.maximum.accumulate(Nmat, axis=1)
    backflow_pp = np.max(CDF - Nmat, axis=1)
    star = int(np.argmax(backflow_pp))
    edges = CDF[:, EDGE_IDX]
    bins = np.column_stack([np.diff(edges, axis=1), 1.0 - CDF[:, -1]])
    N_edges = Nmat[:, EDGE_IDX]
    # a bin's increment is contaminated by backflow iff the running max is
    # pinned above N at EITHER of its edges for some preparation (an interior
    # dip that recovers before the edge leaves the increment quadratic);
    # the non-arrival bin (index 48) is contaminated iff CDF(T) != N(T).
    edge_gap = np.max(edges - N_edges, axis=0)            # (49,), >= 0
    bin_bf_dev = np.maximum(edge_gap[:-1], edge_gap[1:])  # (48,) increments
    bin_bf_dev = np.append(bin_bf_dev, edge_gap[-1])      # (49,) + non-arrival
    bf_bins = np.where(bin_bf_dev > 1e-10)[0]
    des = design_g(th, ph, overlap)
    # comparator floor on this grid: presence probability at bin edges
    floor_pb = np.max(np.abs(N_edges - des @ np.linalg.lstsq(
        des, N_edges, rcond=None)[0]), axis=0)
    resid_pb = np.max(np.abs(bins - des @ np.linalg.lstsq(
        des, bins, rcond=None)[0]), axis=0)
    out = {
        "n_preps": int(len(th)),
        "max_backflow": float(backflow_pp[star]),
        "star": {"theta": float(th[star]), "phi": float(ph[star]),
                 "theta_over_pi": float(th[star] / np.pi),
                 "phi_over_pi": float(ph[star] / np.pi)},
        "floor": float(floor_pb.max()),
        "resid_perbin": resid_pb.tolist(),
        "resid_global_max": float(resid_pb.max()),
        "backflow_bins": bf_bins.tolist(),
        "bin_backflow_dev": bin_bf_dev.tolist(),
        "max_resid_in_backflow_bins":
            float(resid_pb[bf_bins].max()) if len(bf_bins) else 0.0,
        "max_resid_in_nobackflow_bins":
            float(np.delete(resid_pb, bf_bins).max()),
    }
    return out, Nmat, CDF, star, resid_pb, bin_bf_dev


def main():
    res = {}
    # ---------------- main family integrals ----------------
    par1, par2 = (l7a.K1, l7a.X0), (l7a.K2, l7a.X0)
    Iaa_R, Ibb_R, Iab_R, Iaa_F, Ibb_F, Iab_F = l7a.pair_integrals(
        *par1, *par2)
    overlap = Iab_F[0]

    # ---------------- pocket map from fine scan ----------------
    th_f, ph_f = grid_thph(121, 128)
    bf_f = np.empty(len(th_f))
    for lo in range(0, len(th_f), 2000):
        hi = min(lo + 2000, len(th_f))
        Nblk, _ = family_N_g(th_f[lo:hi], ph_f[lo:hi],
                             Iaa_R, Ibb_R, Iab_R, overlap)
        bf_f[lo:hi] = np.max(np.maximum.accumulate(Nblk, axis=1) - Nblk,
                             axis=1)
    pocket = bf_f > 1e-5
    res["pocket_map"] = {
        "fine_grid": [121, 128],
        "max_backflow": float(bf_f.max()),
        "argmax": {"theta": float(th_f[np.argmax(bf_f)]),
                   "phi": float(ph_f[np.argmax(bf_f)])},
        "n_pocket_points_gt_1e-5": int(pocket.sum()),
        "theta_range": [float(th_f[pocket].min()), float(th_f[pocket].max())],
        "phi_range": [float(ph_f[pocket].min()), float(ph_f[pocket].max())],
    }

    # ---------------- registered grid (7x8) ----------------
    th0, ph0 = l7a.family_grid()
    reg, *_ = analyze(th0, ph0, Iaa_R, Ibb_R, Iab_R, overlap)
    res["registered_grid_7x8"] = reg

    # ---------------- refinements (resolution retries) ----------------
    refinements = [("13x16", 13, 16), ("25x32", 25, 32)]
    chosen = None
    for name, nth, nph in refinements:
        th, ph = grid_thph(nth, nph)
        out, Nmat, CDF, star, resid_pb, bin_bf = analyze(
            th, ph, Iaa_R, Ibb_R, Iab_R, overlap)
        res[f"refined_grid_{name}"] = out
        thr = max(1e3 * out["floor"], 1e-4)
        ok_bins = [int(i) for i in out["backflow_bins"]
                   if resid_pb[i] > thr]
        res[f"refined_grid_{name}"]["G3_threshold"] = thr
        res[f"refined_grid_{name}"]["G3_passing_bins"] = ok_bins
        if len(ok_bins) > 0 and chosen is None:
            chosen = (name, th, ph, out, Nmat, CDF, star, resid_pb, bin_bf)
    if chosen is None:
        # keep the last refinement for reporting
        chosen = (name, th, ph, out, Nmat, CDF, star, resid_pb, bin_bf)
    name, th, ph, out, Nmat, CDF, star, resid_pb, bin_bf = chosen
    res["chosen_refinement"] = name
    res["G3_after_refinement"] = {
        "passing_bins": res[f"refined_grid_{name}"]["G3_passing_bins"],
        "threshold": res[f"refined_grid_{name}"]["G3_threshold"],
        "max_resid_in_backflow_bins": out["max_resid_in_backflow_bins"],
        "pass": len(res[f"refined_grid_{name}"]["G3_passing_bins"]) > 0}

    # ---------------- G1 fan at refined-grid star (genuine backflow) --------
    ecdf, min_gap, n_inv, n_crossed = l7a.run_fan(
        th[star], ph[star], overlap, par1, par2, dt=0.001)
    dkw = np.sqrt(np.log(2.0 / 0.05) / (2.0 * l7a.NTRAJ))
    sup = float(np.max(np.abs(ecdf - CDF[star])))
    res["G1_fan_main_star"] = {
        "theta_over_pi": float(th[star] / np.pi),
        "phi_over_pi": float(ph[star] / np.pi),
        "backflow": out["max_backflow"],
        "n_traj": l7a.NTRAJ, "dt_traj": 0.001,
        "sup_norm_ecdf_vs_runmax": sup,
        "dkw_95_band": float(dkw), "dkw_2x_band": float(2 * dkw),
        "min_ordering_gap": float(min_gap),
        "n_steps_with_inversion": int(n_inv),
        "n_crossed": int(n_crossed)}

    # component norm conservation (recheck, both families)
    res["component_norm_max_dev_main"] = float(
        max(np.max(np.abs(Iaa_F - 1)), np.max(np.abs(Ibb_F - 1))))

    # ---------------- G4 control, favorable single adjustment --------------
    KC = 1.8
    g4 = {"k_control": KC}
    for tag, x0b in (("initial_x0b_-11", -11.0), ("adjusted_x0b_-9.5", -9.5)):
        cI = l7a.pair_integrals(KC, -8.0, KC, x0b)
        cIaa_R, cIbb_R, cIab_R, cIaa_F, cIbb_F, cIab_F = cI
        ovc = cIab_F[0]
        outc, cN, cCDF, cstar, cresid_pb, cbin_bf = analyze(
            th0, ph0, cIaa_R, cIbb_R, cIab_R, ovc)
        g4[tag] = {
            "x0a": -8.0, "x0b": x0b,
            "overlap_abs": float(np.abs(ovc)),
            "norm_max_dev": float(max(np.max(np.abs(cIaa_F - 1)),
                                      np.max(np.abs(cIbb_F - 1)))),
            "max_backflow_on_registered_grid": outc["max_backflow"],
            "monotone": bool(outc["max_backflow"] <= 1e-10),
            "floor": outc["floor"],
            "resid_global_max": outc["resid_global_max"],
            "resid_perbin": outc["resid_perbin"],
            "backflow_bins": outc["backflow_bins"],
            "max_resid_in_nobackflow_bins":
                outc["max_resid_in_nobackflow_bins"],
        }
        if tag.startswith("adjusted"):
            ctrl_resid_pb = cresid_pb
            ctrl_bin_bf = cbin_bf
    g4["separation_scan_max_backflow_continuum"] = {
        "1.0": 3.12e-4, "1.5": 3.32e-4, "2.0": 3.65e-4,
        "3.0": 4.79e-4, "4.0": 6.38e-4,
        "note": "from l7a_diagnostics.json fine-scan; no separation in "
                "[1,4] is backflow-free -> physics failure of the "
                "no-backflow-control design, not a resolution issue"}
    res["G4"] = g4

    # backflow-borne correlation across both families (salvaged G4 content)
    res["backflow_borne_correlation"] = {
        "main_refined_max_resid_nobackflow_bins":
            out["max_resid_in_nobackflow_bins"],
        "main_refined_max_resid_backflow_bins":
            out["max_resid_in_backflow_bins"],
        "ctrl_adjusted_max_resid_nobackflow_bins":
            g4["adjusted_x0b_-9.5"]["max_resid_in_nobackflow_bins"],
        "ctrl_adjusted_max_resid_backflow_bins":
            float(np.max(np.asarray(g4["adjusted_x0b_-9.5"]["resid_perbin"])[
                g4["adjusted_x0b_-9.5"]["backflow_bins"]]))
        if g4["adjusted_x0b_-9.5"]["backflow_bins"] else 0.0}

    # ---------------- final gates ----------------
    g1 = (res["component_norm_max_dev_main"] <= 1e-10
          and out["max_backflow"] > 1e-4
          and sup <= 2 * dkw and n_inv == 0)
    g2 = reg["floor"] <= 1e-6 and out["floor"] <= 1e-6
    g3 = res["G3_after_refinement"]["pass"]
    g4_pass = (g4["adjusted_x0b_-9.5"]["monotone"]
               and g4["adjusted_x0b_-9.5"]["resid_global_max"]
               <= 10 * reg["floor"])
    res["gates"] = {
        "G1": "PASS" if g1 else "FAIL",
        "G2": "PASS" if g2 else "FAIL",
        "G3": ("PASS_after_resolution_refinement" if g3 else "FAIL")
        if reg["max_backflow"] < 1e-10 else ("PASS" if g3 else "FAIL"),
        "G4": "PASS" if g4_pass else "FAIL_physics",
    }

    with open(f"{BASE}/l7a_refine.json", "w") as f:
        json.dump(res, f, indent=2)
    np.savez(f"{BASE}/l7a_refine_arrays.npz",
             TGRID=TGRID, Nstar=Nmat[star], CDFstar=CDF[star], ecdf=ecdf,
             resid_pb=resid_pb, bin_bf=bin_bf,
             reg_resid_pb=np.array(reg["resid_perbin"]),
             ctrl_resid_pb=ctrl_resid_pb, ctrl_bin_bf=ctrl_bin_bf,
             th_star=th[star], ph_star=ph[star], dkw=dkw,
             bf_fine=bf_f, th_fine=th_f, ph_fine=ph_f)
    # console summary (omit long lists)
    slim = json.loads(json.dumps(res))
    for k in list(slim):
        if isinstance(slim[k], dict):
            for kk in ("resid_perbin", "bin_backflow_dev"):
                slim[k].pop(kk, None)
            for kk in list(slim[k]):
                if isinstance(slim[k][kk], dict):
                    slim[k][kk].pop("resid_perbin", None)
                    slim[k][kk].pop("bin_backflow_dev", None)
    print(json.dumps(slim, indent=2))


if __name__ == "__main__":
    main()
