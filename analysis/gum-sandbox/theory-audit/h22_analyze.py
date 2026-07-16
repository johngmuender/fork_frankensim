#!/usr/bin/env python3
"""Phase H2.2 — NR-D2 P4 two-texture reconciliation test: analysis.

Reads the deterministic h22_solve run JSONs (h22_runs/) and assembles the
verdict-relevant numbers into h22_results.json + a printed table:

  * naive mismatch wall dE_naive (the engine analog of the corpus's f^2 p
    per-crossing cost): E_seed(repulse) - E_seed(align), decomposed into
    the SDiff-frozen value part D(E6+E0) and the eps-suppressed kinetic
    part t*D(E2+E4);
  * cost_full  = E_fin(repulse, guarded ANF) - E_fin(align, guarded ANF);
  * cost_sdiff = E_fin(repulse, SDiff-projected descent) - E_fin(align);
  * cost_twist = min_alpha E(seed o T_alpha) - E_fin(align) over the exact
    volume-preserving interpolating-twist family;
  * the corpus screening scale S = eps * dE_naive at the engine dial
    eps = 0.05 (unit map stated in h22_RESULTS.md), and the ratios to it;
  * the F-R4 floor checks: (i) value-sector freeze under SDiff moves
    (measured drift of E6+E0 and degree along the projected descent and
    the twist scan), (ii) whether every SDiff endpoint respects the
    rigorous value floor E >= (E6+E0)_seed, (iii) whether the eps-lifted
    cost is a lower bound for both relaxations.

Epistemic notice: within-model numbers on a speculative theory's
functional; measured facts only, adjudication is the coordinator's.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "h22_runs")
OUT = os.path.join(HERE, "h22_results.json")

EPS_ENGINE = 0.05  # the frozen eps dial: t(E2+E4)/(E6+E0) at the hedgehog
T_FROZEN = 0.008276434949296802
CORPUS = {
    "fp_in_f_units": 7.6e7,          # f*p (NR-D2 section 2)
    "naive_wall": "f^2 p per crossing (K = f^2, mismatch K/p^2 over p^3)",
    "screened": "eps * f^2 p (P4's claimed SDiff-screened cost)",
    "eps_band": [1.0e-5, 6.0e-4],    # corpus frustration-sector band
    "eps_mapped_here": EPS_ENGINE,   # per task: map eps = 0.05, engine units
}


def load(name):
    with open(os.path.join(RUNS, name)) as fh:
        return json.load(fh)


def sect(o):
    """corner-scheme bare sector summary of an Out dict."""
    return {
        "estat": o["estat"],
        "e2": o["e2"],
        "e4": o["e4"],
        "e6": o["e6"],
        "e0": o["e0"],
        "value": o["e6"] + o["e0"],
        "kinetic_t": T_FROZEN * (o["e2"] + o["e4"]),
        "deg": o["deg"],
    }


def tail_flatness(series, frac=0.2):
    """relative estat change over the last `frac` of the recorded series."""
    if not series:
        return None
    k = max(1, int(len(series) * frac))
    tail = series[-k:]
    e0, e1 = tail[0]["estat"], tail[-1]["estat"]
    return abs(e1 - e0) / max(abs(e1), 1e-300)


def main():
    single = load("h22_single.json")
    al = load("h22_pair_align_m09.json")
    rp = load("h22_pair_repulse_m09.json")
    sd = load("h22_sdiff_repulse_m09.json")
    tw = load("h22_twist_repulse_m09.json")
    sc = load("h22_seedscan_m06_m20.json")

    s_al_seed = sect(al["seed"]["corner"])
    s_rp_seed = sect(rp["seed"]["corner"])
    s_al_fin = sect(al["final"]["corner"])
    s_rp_fin = sect(rp["final"]["corner"])
    s_sd_fin = sect(sd["final"]["corner"])
    s_single = sect(single["final"]["corner"])

    e_ref = s_al_fin["estat"]          # reconciled reference
    e_single2 = 2.0 * s_single["estat"]

    # --- the naive wall (engine analog of f^2 p) and its decomposition -----
    de_naive = s_rp_seed["estat"] - s_al_seed["estat"]
    dv_naive = s_rp_seed["value"] - s_al_seed["value"]
    dk_naive = s_rp_seed["kinetic_t"] - s_al_seed["kinetic_t"]

    # --- costs --------------------------------------------------------------
    cost_full = s_rp_fin["estat"] - e_ref
    cost_sdiff = s_sd_fin["estat"] - e_ref

    scan = tw["scan"]
    tw_estat = [r["corner"]["estat"] for r in scan]
    tw_value = [r["corner"]["e6"] + r["corner"]["e0"] for r in scan]
    tw_deg = [r["corner"]["deg"] for r in scan]
    i_min = min(range(len(scan)), key=lambda i: tw_estat[i])
    cost_twist = tw_estat[i_min] - e_ref
    tw_value_spread = max(tw_value) - min(tw_value)
    tw_deg_spread = max(tw_deg) - min(tw_deg)

    # --- screening scale (the corpus's eps * f^2 p, engine units) ----------
    # the engine's naive wall is NEGATIVE (see results: the aligned pair is
    # the costlier seed); the corpus scale is a magnitude, so |dE_naive|
    s_scale = EPS_ENGINE * abs(de_naive)
    # measured dial at the relaxed single knot (cross-check of eps = 0.05)
    eps_measured = s_single["kinetic_t"] / s_single["value"]

    # --- SDiff verdict quantities -------------------------------------------
    sdiff_total_reduction = s_rp_seed["estat"] - s_sd_fin["estat"]
    full_total_reduction = s_rp_seed["estat"] - s_rp_fin["estat"]
    sdiff_value_drift = sd["descent"]["value_drift_final"]
    sdiff_deg_drift = sd["descent"]["deg_drift_final"]
    sdiff_kinetic_reduction = s_rp_seed["kinetic_t"] - s_sd_fin["kinetic_t"]
    # rigorous SDiff value floor: any SDiff image of the mismatch seed has
    # E >= (E6+E0)(seed) since t(E2+E4) >= 0
    value_floor_seed = s_rp_seed["value"]
    sdiff_above_value_floor = s_sd_fin["estat"] - value_floor_seed
    full_final_value_change = s_rp_fin["value"] - s_rp_seed["value"]
    # reduction beyond the measured constraint leak (the continuum-legal
    # part of the SDiff reduction; <= 0 means SDiff moves achieved nothing
    # beyond the numerical volume-preservation error)
    sdiff_reduction_beyond_leak = sdiff_total_reduction + sdiff_value_drift

    # --- seed scan: sign structure of the naive wall across separations ----
    seed_scan = [
        {"m": r["m"], "x_sep": r["x_sep"],
         "dE_rep_minus_align": r["repulse"]["estat"] - r["align"]["estat"],
         "dValue": (r["repulse"]["e6"] + r["repulse"]["e0"])
                   - (r["align"]["e6"] + r["align"]["e0"])}
        for r in sc["scan"]
    ]
    no_positive_wall = all(r["dE_rep_minus_align"] < 0 for r in seed_scan)

    tw_monotone = all(b >= a for a, b in zip(tw_estat, tw_estat[1:]))
    tw_e0_spread = (max(r["corner"]["e0"] for r in scan)
                    - min(r["corner"]["e0"] for r in scan))

    results = {
        "phase": "H2.2 — NR-D2 P4 two-texture reconciliation test",
        "date": "2026-07-16",
        "corpus_scales": CORPUS,
        "protocol": {
            "grid": "Field3 n=96 cubic, half=4.5, h=0.09375 (frozen)",
            "separation": {"m": 9, "d": rp["protocol"]["d"],
                           "x_sep_d_over_rstar": rp["protocol"]["x_sep"]},
            "mismatch_channel": "repulse (theta=pi about the separation axis)",
            "reconciled_reference": "align (identity relative rotation), full guarded ANF",
            "guards": rp["guards"],
            "anf": {"align": al["anf"], "repulse": rp["anf"],
                    "single": single["anf"]},
            "sdiff_descent": sd["descent"],
            "eps_engine_dial": EPS_ENGINE,
            "eps_measured_single": eps_measured,
        },
        "states": {
            "single_final": s_single,
            "align_seed": s_al_seed, "align_final": s_al_fin,
            "repulse_seed": s_rp_seed, "repulse_final": s_rp_fin,
            "sdiff_final": s_sd_fin,
            "twist_min": {"alpha": scan[i_min]["alpha"],
                          **sect(scan[i_min]["corner"])},
            "twist_alpha_pi": sect(scan[-1]["corner"]),
        },
        "costs": {
            "reference_E": e_ref,
            "reference_2x_single": e_single2,
            "dE_naive_seed_mismatch": de_naive,
            "dE_naive_value_part": dv_naive,
            "dE_naive_kinetic_part": dk_naive,
            "cost_full": cost_full,
            "cost_sdiff_projected": cost_sdiff,
            "cost_twist_best": cost_twist,
            "E_int_align_vs_2single": s_al_fin["estat"] - e_single2,
            "E_int_repulse_vs_2single": s_rp_fin["estat"] - e_single2,
        },
        "screening": {
            "screening_scale_S_eq_eps_dE_naive": s_scale,
            "cost_full_over_S": cost_full / s_scale,
            "cost_sdiff_over_S": cost_sdiff / s_scale,
            "eps_eff_full_eq_cost_full_over_dE_naive": cost_full / de_naive,
            "corpus_P4_prediction_eps_eff": EPS_ENGINE,
        },
        "sdiff_projection": {
            "total_reduction_sdiff": sdiff_total_reduction,
            "total_reduction_full": full_total_reduction,
            "sdiff_capture_fraction_of_full":
                (sdiff_total_reduction / full_total_reduction)
                if full_total_reduction != 0 else None,
            "kinetic_reduction_sdiff": sdiff_kinetic_reduction,
            "value_drift_sdiff_final": sdiff_value_drift,
            "deg_drift_sdiff_final": sdiff_deg_drift,
            "value_change_full_final": full_final_value_change,
            "reduction_beyond_leak": sdiff_reduction_beyond_leak,
            "status": sd["descent"]["status"],
            "budget_hits": sd["descent"].get("budget_hits"),
            "detj_max": sd["final"].get("detj_max"),
            "detj_mean": sd["final"].get("detj_mean"),
            "twist_value_spread_over_scan": tw_value_spread,
            "twist_e0_spread_over_scan": tw_e0_spread,
            "twist_deg_spread_over_scan": tw_deg_spread,
            "twist_estat_monotone_increasing": tw_monotone,
            "sdiff_map_dmax": sd["final"]["map_dmax"],
            "sdiff_map_dl2": sd["final"]["map_dl2"],
        },
        "seed_scan_sign_structure": {
            "rows": seed_scan,
            "no_positive_wall_regime_found": no_positive_wall,
            "note": "E_seed(repulse) < E_seed(align) at every overlapping "
                    "separation; the corpus's positive naive wall does not "
                    "reproduce in sign in this sector, and the difference "
                    "is 94-98% value-sector (SDiff-frozen).",
        },
        "floor_tests": {
            "value_floor_seed_E6E0": value_floor_seed,
            "sdiff_final_minus_value_floor": sdiff_above_value_floor,
            "sdiff_respects_value_floor": sdiff_above_value_floor >= 0.0,
            "cost_full_geq_S": cost_full >= s_scale,
            "cost_sdiff_geq_S": cost_sdiff >= s_scale,
            "cost_twist_geq_S": cost_twist >= s_scale,
        },
        "caveats": {
            "align_endpoint_penalties": {"epen": al["final"]["corner"]["epen"],
                                         "efpen": al["final"]["corner"]["efpen"]},
            "repulse_endpoint_penalties": {"epen": rp["final"]["corner"]["epen"],
                                           "efpen": rp["final"]["corner"]["efpen"]},
            "align_anf_tail_flatness": tail_flatness(al["series"]),
            "repulse_anf_tail_flatness": tail_flatness(rp["series"]),
            "d_eff": {"align_seed": al["seed"]["d_eff"],
                      "align_final": al["final"]["d_eff"],
                      "repulse_seed": rp["seed"]["d_eff"],
                      "repulse_final": rp["final"]["d_eff"]},
            "boundary_tails": {"align_final": al["final"]["boundary_tail"],
                               "repulse_final": rp["final"]["boundary_tail"],
                               "sdiff_final": sd["final"]["boundary_tail"]},
            "dist_to_align_seed": {
                "repulse_seed": rp["seed"]["dist_to_align_seed"],
                "repulse_final": rp["final"]["dist_to_align_seed"],
                "sdiff_final": sd["final"]["dist_to_align_seed"],
                "twist_min": scan[i_min]["dist_to_align_seed"],
            },
        },
        "twist_scan": [
            {"alpha": r["alpha"], "estat": r["corner"]["estat"],
             "e2": r["corner"]["e2"], "e4": r["corner"]["e4"],
             "e6": r["corner"]["e6"], "e0": r["corner"]["e0"],
             "deg": r["corner"]["deg"],
             "dist_to_align_seed": r["dist_to_align_seed"]}
            for r in scan
        ],
    }

    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=1)
    print("wrote", OUT)

    # ---- printed table -------------------------------------------------------
    def row(k, v, note=""):
        print("  %-42s %+14.6f   %s" % (k, v, note))

    print("\n== H2.2 cost table (engine units, corner scheme, bare E_static) ==")
    row("E_ref (align pair, relaxed)", e_ref)
    row("2 x E(single, relaxed)", e_single2)
    row("dE_naive (seed mismatch wall ~ f^2 p)", de_naive)
    row("  value part D(E6+E0)", dv_naive)
    row("  kinetic part t D(E2+E4)", dk_naive)
    row("cost_full (ANF endpoint - E_ref)", cost_full)
    row("cost_sdiff_projected", cost_sdiff)
    row("cost_twist_best (alpha=%.3f)" % scan[i_min]["alpha"], cost_twist)
    row("screening scale S = 0.05 * dE_naive", s_scale)
    print("\n== ratios ==")
    row("cost_full / S", cost_full / s_scale)
    row("cost_sdiff / S", cost_sdiff / s_scale)
    row("eps_eff = cost_full / dE_naive", cost_full / de_naive,
        "(corpus P4 predicts ~ %.3f)" % EPS_ENGINE)
    print("\n== SDiff projection ==")
    row("reduction achieved (sdiff)", sdiff_total_reduction,
        "status=%s budget_hits=%s" % (sd["descent"]["status"],
                                      sd["descent"].get("budget_hits")))
    row("reduction achieved (full ANF)", full_total_reduction)
    row("value drift along sdiff (E6+E0)", sdiff_value_drift,
        "(frozen-value check; 0 in continuum)")
    row("reduction beyond measured leak", sdiff_reduction_beyond_leak,
        "(<= 0: nothing beyond numerical leak)")
    row("value change along full ANF", full_final_value_change)
    row("twist-scan E6+E0 spread", tw_value_spread,
        "(exact SDiff family; E0 spread %.2e)" % tw_e0_spread)
    print("\n== floor tests ==")
    row("sdiff endpoint - value floor", sdiff_above_value_floor,
        ">= 0 required by F-R4 lemma")
    print("  cost_full >= S: %s | cost_sdiff >= S: %s | cost_twist >= S: %s"
          % (cost_full >= s_scale, cost_sdiff >= s_scale, cost_twist >= s_scale))
    return 0


if __name__ == "__main__":
    sys.exit(main())
