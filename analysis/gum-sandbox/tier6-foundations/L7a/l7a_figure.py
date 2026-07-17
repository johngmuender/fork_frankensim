#!/usr/bin/env python3
"""L7a: final figure + merged L7a_results.json."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7a"

ref = np.load(f"{BASE}/l7a_refine_arrays.npz")
diag = np.load(f"{BASE}/l7a_diag_arrays.npz")
with open(f"{BASE}/l7a_refine.json") as f:
    R = json.load(f)
with open(f"{BASE}/l7a_diagnostics.json") as f:
    DG = json.load(f)
with open(f"{BASE}/L7a_results.json") as f:
    S1 = json.load(f)
if "normalization" in S1:          # already-merged file from a previous run
    S1 = dict(S1["normalization"], params=S1["params"],
              n2_family_min_max=S1["normalization"]["n2_family_min_max"])

T = ref["TGRID"]
Nstar, CDFstar, ecdf = ref["Nstar"], ref["CDFstar"], ref["ecdf"]
sup_ecdf_cdf = float(np.max(np.abs(ecdf - CDFstar)))
sup_ecdf_N = float(np.max(np.abs(ecdf - Nstar)))

# ---------------- palette (fixed categorical order, CVD-safe) ----------------
BLUE, ORANGE, RED, TEAL, GRAY = ("#4269d0", "#efb118", "#ff725c",
                                 "#6cc5b0", "#9498a0")
INK, MUTED = "#1f2430", "#6b7280"
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": "#d7dbe0", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": "#e8eaee", "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 9.5,
    "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "legend.frameon": False, "legend.fontsize": 8.5,
})

fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.2))
fig.subplots_adjust(hspace=0.42, wspace=0.30, left=0.075, right=0.975,
                    top=0.90, bottom=0.08)
fig.suptitle("L7a - POVM exclusion via quadraticity of Bohmian first-arrival "
             "statistics (within-model testbed, F-T6-L7a)",
             fontsize=12, fontweight="bold", x=0.075, ha="left", color=INK)

# ---- (a) presence probability vs running-max CDF at the star preparation ----
ax = axs[0, 0]
ax.plot(T, Nstar, color=BLUE, lw=1.6, label=None)
ax.plot(T, CDFstar, color=ORANGE, lw=1.6, ls="--", label=None)
ax.set_title("(a) Backflow at the star preparation")
ax.set_xlabel("t"); ax.set_ylabel("probability beyond detector")
ax.text(8.2, 0.42, "N(t) = P(x>0)", color=BLUE, fontsize=9)
ax.text(8.2, 0.30, "CDF(t) = running max", color="#c98d09", fontsize=9)
ax.text(0.3, 0.60, r"$\theta=0.229\pi,\ \phi=0.0625\pi$"
        "\nmax backflow = 2.39e-4", fontsize=8.5, color=MUTED)
# inset zoom on the dip
axi = ax.inset_axes([0.55, 0.12, 0.42, 0.38])
m = (T >= 3.8) & (T <= 5.0)
axi.plot(T[m], Nstar[m], color=BLUE, lw=1.4)
axi.plot(T[m], CDFstar[m], color=ORANGE, lw=1.4, ls="--")
axi.set_xlim(3.8, 5.0)
lo = Nstar[m].min(); hi = CDFstar[m].max()
axi.set_ylim(lo - 1e-4, hi + 1e-4)
axi.tick_params(labelsize=7)
axi.set_title("dip, zoom", fontsize=7.5, color=MUTED)
ax.indicate_inset_zoom(axi, edgecolor=MUTED, alpha=0.5)

# ---- (b) per-bin residual profile, main family ----
ax = axs[0, 1]
bins = np.arange(49)
resid_ref = ref["resid_pb"]
resid_reg = ref["reg_resid_pb"]
bf_bins = np.where(ref["bin_bf"] > 1e-10)[0]
for b in bf_bins:
    ax.axvspan(b - 0.5, b + 0.5, color="#f6e7c8", zorder=0, lw=0)
ax.semilogy(bins, np.maximum(resid_ref, 1e-18), color=RED, lw=1.6,
            marker="o", ms=2.6, mew=0)
ax.semilogy(bins, np.maximum(resid_reg, 1e-18), color=BLUE, lw=1.4,
            marker="o", ms=2.4, mew=0)
ax.axhline(2.9e-15, color=GRAY, lw=1.2, ls=":")
ax.axhline(1e-4, color=INK, lw=0.9, ls="--", alpha=0.5)
ax.set_title("(b) POVM-form residual per time bin (main family)")
ax.set_xlabel("bin (0.25 wide; 48 = non-arrival)")
ax.set_ylabel("max |data - quadratic fit|")
ax.set_ylim(1e-18, 3e-3)
ax.text(1, 4e-4, "refined 25x32 grid", color=RED, fontsize=9)
ax.text(1, 6e-15, "registered 7x8 grid", color=BLUE, fontsize=9)
ax.text(35.5, 4.5e-15, "G2 floor", color=MUTED, fontsize=8)
ax.text(35.5, 1.6e-4, "G3 abs. threshold", color=MUTED, fontsize=8)
ax.text(14, 1e-17, "shaded = backflow bins", color="#b98a1e", fontsize=8)

# ---- (c) trajectory-fan validation of the running-max CDF ----
ax = axs[1, 0]
ax.plot(T, (CDFstar - Nstar) / 1e-4, color=RED, lw=1.8)
ax.plot(T, (ecdf - CDFstar) / 1e-4, color=BLUE, lw=0.9, alpha=0.85)
ax.set_ylim(-2.2, 3.2)
ax.set_title("(c) 4000-trajectory fan vs exact CDF (star prep.)")
ax.set_xlabel("t"); ax.set_ylabel("deviation  [1e-4]")
ax.text(5.6, 2.55, "backflow signal CDF(t) - N(t)  (max 2.39e-4)",
        color=RED, fontsize=9)
ax.text(0.3, -1.95,
        "fan error ECDF - runmax CDF: sup = 1.25e-4 = 1/2n\n"
        "(stratified-quantile resolution; 2x DKW band 4.3e-2 off scale)",
        color=BLUE, fontsize=9)

# ---- (d) control family: violation is backflow-borne ----
ax = axs[1, 1]
c_init = np.array(R["G4"]["initial_x0b_-11"]["resid_perbin"])
c_adj = np.array(R["G4"]["adjusted_x0b_-9.5"]["resid_perbin"])
for b in R["G4"]["initial_x0b_-11"]["backflow_bins"]:
    ax.axvspan(b - 0.5, b + 0.5, color="#f6e7c8", zorder=0, lw=0)
ax.semilogy(bins, np.maximum(c_init, 1e-18), color=TEAL, lw=1.6,
            marker="o", ms=2.6, mew=0)
ax.semilogy(bins, np.maximum(c_adj, 1e-18), color=GRAY, lw=1.4,
            marker="o", ms=2.4, mew=0)
ax.axhline(1.33e-15, color=GRAY, lw=1.2, ls=":")
ax.set_ylim(1e-18, 3e-3)
ax.set_title("(d) Control family (same k = 1.8, displaced centers)")
ax.set_xlabel("bin (0.25 wide; 48 = non-arrival)")
ax.set_ylabel("max |data - quadratic fit|")
ax.text(1, 4e-4, "sep 3.0: backflow 4.3e-4", color="#3c9c86", fontsize=9)
ax.text(1, 6e-15, "sep 1.5 (adjusted): monotone, at floor",
        color=MUTED, fontsize=9)
ax.text(14, 1e-17, "shaded = backflow bins (sep 3.0)",
        color="#b98a1e", fontsize=8)

fig.savefig(f"{BASE}/L7a_povm_exclusion.png", dpi=170)
print("figure written")

# ---------------- merged final results JSON ----------------
final = {
    "workstream": "F-T6-L7a",
    "schema_note": ("A fixed physical detector induces one preparation-"
                    "independent POVM {E(dt)} via Naimark; over a family of "
                    "preparations, POVM-realizable arrival statistics must be "
                    "a sesquilinear functional <psi|E(B)|psi>. A single "
                    "preparation excludes nothing; all content is in the "
                    "family."),
    "params": S1["params"],
    "normalization": {
        "overlap_analytic": S1["overlap_analytic"],
        "overlap_numeric_t0": S1["overlap_numeric_t0"],
        "overlap_abs": S1["overlap_abs"],
        "overlap_max_dev_from_analytic": S1["overlap_max_dev_from_analytic"],
        "overlap_max_dev_from_const_in_t": S1["overlap_max_dev_from_const_in_t"],
        "component_norm_max_dev": S1["component_norm_max_dev"],
        "n2_family_min_max": S1["n2_family_min_max"],
        "note": ("<psi1|psi2> = 0.698, NOT zero; family members are "
                 "normalized explicitly, and the quadratic-form fit uses "
                 "basis functions divided by n^2(theta,phi) so that a POVM "
                 "functional is exactly linear in the fitted parameters.")},
    "registered_grid_7x8": R["registered_grid_7x8"],
    "backflow_pocket_map": R["pocket_map"],
    "current_matrix_diagnostics": {
        "main_lambda_min": DG["main_lambda_min_over_t"],
        "main_lambda_min_argmin_t": DG["main_lambda_min_argmin_t"],
        "local_dk_at_detector": DG["local_dk_at_detector"],
        "control_separation_scan": DG["control_separation_scan"],
        "ctrl_continuity_max_dev": DG["ctrl_continuity_max_dev"],
        "ctrl_min_current": DG["ctrl_min_current"]},
    "refined_grid_13x16": R["refined_grid_13x16"],
    "refined_grid_25x32": R["refined_grid_25x32"],
    "chosen_refinement": R["chosen_refinement"],
    "G1_fan_main_star": dict(R["G1_fan_main_star"],
                             sup_norm_ecdf_vs_N=sup_ecdf_N),
    "G1_fan_ctrl_star_sep3": DG["ctrl_fan"],
    "G2": {"floor_registered": R["registered_grid_7x8"]["floor"],
           "floor_refined_25x32": R["refined_grid_25x32"]["floor"],
           "spec": "<= 1e-6"},
    "G3": dict(R["G3_after_refinement"],
               registered_grid_result={
                   "max_backflow": R["registered_grid_7x8"]["max_backflow"],
                   "resid_global_max":
                       R["registered_grid_7x8"]["resid_global_max"],
                   "verdict_at_registered_resolution":
                       "no exclusion (grid misses the backflow pocket)"}),
    "G4": R["G4"],
    "backflow_borne_correlation": R["backflow_borne_correlation"],
    "gates": R["gates"],
}
with open(f"{BASE}/L7a_results.json", "w") as f:
    json.dump(final, f, indent=2)
print("L7a_results.json merged")
print("sup|ECDF-CDF| =", sup_ecdf_cdf, " sup|ECDF-N| =", sup_ecdf_N)
