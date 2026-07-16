#!/usr/bin/env python3
"""
K1 — the /6 vs /(15pi/8) discriminator (J1's residual split), executed.

WITHIN-MODEL ONLY.  J1 recovered the edge-referenced amplitude convention
AS A CLASS: p = A e^{-mu R*}/F_p with F_p in {6, 15pi/8}.  The two
candidates differ by the exact factor 48/(15pi) = 1.018592 in p (1.86%).
This script (a) enumerates every archive quantity that could be sensitive
to the choice, with its sensitivity exponent and the archive's printed
precision; (b) executes every channel that has any discrimination power;
(c) prints a PASS/FAIL gate table.  A verified null (nothing in the
archive separates the candidates at >= 2 sigma) is a legitimate terminal
state and is what the gates below test for.

Inputs (frozen records only; nothing re-run, nothing modified):
  tier2-closure/radial_results.json                     A, mu, R*, f(R*)
  theory-audit/j1_results.json                          the B.4 p-grid
  theory-audit/i3_results.json                          targets + factor family
  gum-core/fs-gum-twoknot/twoknot_results.json          C_d fits, running b_eff
  gum-core/fs-gum-twoknot/twoknot_refine_results.json   G5b C_d, extrap arrays
  gum-core/fs-gum-twoknot/runs/run_single.json          record-content scan

Deterministic, no RNG, runtime ~2 s.  Writes k1_results.json STAGE-FLUSHED.
"""

import json
import math
import os
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SBOX = os.path.dirname(HERE)

RADIAL = os.path.join(SBOX, "tier2-closure", "radial_results.json")
J1 = os.path.join(HERE, "j1_results.json")
I3 = os.path.join(HERE, "i3_results.json")
TWOKNOT = os.path.join(SBOX, "gum-core", "fs-gum-twoknot", "twoknot_results.json")
REFINE = os.path.join(SBOX, "gum-core", "fs-gum-twoknot", "twoknot_refine_results.json")
RUNSINGLE = os.path.join(SBOX, "gum-core", "fs-gum-twoknot", "runs", "run_single.json")

OUT = os.path.join(HERE, "k1_results.json")

# ---------------------------------------------------------------- targets
P_TGT, P_SIG = 0.84, 0.03            # corpus VII.C / I.1 / session map
B_TGT, B_SIG = 42.0, 6.0             # corpus VII.D / I.2
X0_LOOP, X0_LOOP_SIG = 1.90, 0.05    # corpus bond equation at b = 42
X0_MEAS, X0_MEAS_SIG = 1.92, 0.08    # corpus measured bond location
BINF = 4.5                           # corpus "b_inf ~ 4.5" (no error printed)

F6 = 6.0
F15 = 15.0 * math.pi / 8.0
SPLIT_RATIO = F6 / F15               # p_15pi8 / p_6 = 48/(15 pi) = 1.018592
SPLIT_LN = math.log(SPLIT_RATIO)     # 0.018421 (1.86% in p)
SPLIT_B = SPLIT_RATIO ** 2           # 1.037530 (3.75% in b, b ~ p^-2)

RESULTS = {"spec": {
    "task": "K1 discriminate p = A e^{-muR*}/6 vs /(15pi/8)",
    "F6": F6, "F15pi8": F15, "split_ratio_p": SPLIT_RATIO,
    "split_pct_p": 100.0 * (SPLIT_RATIO - 1.0),
    "split_pct_b": 100.0 * (SPLIT_B - 1.0),
    "discrimination_bar_sigma": 2.0,
}}


def flush():
    with open(OUT, "w") as fh:
        json.dump(RESULTS, fh, indent=1, default=float)


GATES = []


def gate(name, ok, measured, threshold, note=""):
    GATES.append(dict(gate=name, ok=bool(ok), measured=measured,
                      threshold=threshold, note=note))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {measured} vs {threshold}"
          f"{'  -- ' + note if note else ''}")


# ---------------------------------------------------------------- inputs
rad_all = json.load(open(RADIAL))
rad = rad_all["2N_run"]
j1 = json.load(open(J1))
i3 = json.load(open(I3))
two = json.load(open(TWOKNOT))
ref = json.load(open(REFINE))

A = rad["A_fit"]
MU = rad["mu_analytic_derived"]
RSTAR = rad_all["RSTAR"]
MURS = MU * RSTAR
ENV = math.exp(-MURS)

D0 = j1["points"]["+0.000"]["2N"]        # J1 benchmark point, 2N
P6 = D0["p_edge_6"]
P15 = D0["p_edge_15pi8"]

CD = {
    "G4-relaxed-corner(cap400,SUPERSEDED)": (
        two["fits"]["corner"]["relaxed_anchored"]["C_d"],
        two["fits"]["corner"]["relaxed_anchored"]["C_d_err"]),
    "G4-seed-corner": (
        two["fits"]["corner"]["seed_ansatz"]["C_d"],
        two["fits"]["corner"]["seed_ansatz"]["C_d_err"]),
    "G4-seed-central4": (
        two["fits"]["central4"]["seed_ansatz"]["C_d"],
        two["fits"]["central4"]["seed_ansatz"]["C_d_err"]),
    "G5b-relaxed-1200": (
        ref["fit_hi"]["C_d"], math.sqrt(ref["fit_hi"]["cov_CB"][0][0])),
}

FACTORS = [(n, v) for n, v in i3["factor_family"]]

X_GRID = [1.4732, 1.6837, 1.8942, 2.1046, 2.4203, 2.7360, 3.0517]
BEFF_SEED = two["bond"]["corner"]["beff_running_seed"]
BEFF_RELAX = two["bond"]["corner"]["beff_running_relaxed"]
X_GRID_EXTRAP = [1.6837, 1.8942, 2.1046, 2.4203]
BEFF_EXTRAP = ref["bond"]["beff_running_extrap"]
P_OLD, F_OLD = 0.84, 2.0 * math.pi   # normalization the arrays were built with

print("=" * 88)
print("K1 — /6 vs /(15pi/8): candidate values at the c2 = 0 benchmark (J1, 2N)")
print("=" * 88)
print(f"  p_6     = {P6:.6f}   (pull {(P6 - P_TGT) / P_SIG:+.3f} vs 0.84 +- 0.03)")
print(f"  p_15pi8 = {P15:.6f}   (pull {(P15 - P_TGT) / P_SIG:+.3f})")
print(f"  exact split 48/(15pi) = {SPLIT_RATIO:.6f}  ->  "
      f"{100 * (SPLIT_RATIO - 1):.3f}% in p, {100 * (SPLIT_B - 1):.3f}% in b\n")

# ======================================================================
# STAGE A — enumeration: every archive quantity vs amplitude sensitivity
# ======================================================================
print("STAGE A — enumeration (sensitivity exponent n: Q ~ p^n; archive precision)")

# split (in sigma of the archive band) = |n| * (SPLIT_RATIO-1) * Q / sigma_Q
def power(n, frac_sigma):
    """Discrimination power in sigma for a quantity Q ~ p^n with archive
    fractional precision frac_sigma (instrument error folded by caller)."""
    return abs(n) * (SPLIT_RATIO - 1.0) / frac_sigma

CATALOG = [
    dict(q="p = 0.84 +- 0.03", src="VII.C(509)/I.1(777)/map47/Course15",
         n=1, frac_sigma=P_SIG / P_TGT, digits="2 dp + band",
         note="the calibrand itself; band is the corpus's own precision"),
    dict(q="b_eff = 42 +- 6", src="VII.D(523)/I.2/A.4",
         n=-2, frac_sigma=B_SIG / B_TGT, digits="2 sf + band",
         note="b ~ C_d/(F_b p^2); instrument C_d error and F_b identity fold in"),
    dict(q="x0(bond eq at b=42) = 1.90 +- 0.05", src="VII.D(523)",
         n=None, frac_sigma=X0_LOOP_SIG / X0_LOOP, digits="2 dp + band",
         note="sensitivity dx0/dln b measured in stage B3 from archived arrays"),
    dict(q="x0(measured) = 1.92 +- 0.08", src="VII.D(523)",
         n=0, frac_sigma=X0_MEAS_SIG / X0_MEAS, digits="2 dp + band",
         note="well location from E_int(x) shape: no amplitude enters; n = 0 exact"),
    dict(q="b_inf ~ 4.5", src="VII.D(523)/Course16",
         n=-2, frac_sigma=None, digits="2 sf, '~', NO band",
         note="no printed error; F_b unpinned; far zone unmeasured by the instrument"),
    dict(q="B.4 window edges [-m^2/2, m^2/6)", src="B.4(741)",
         n=0, frac_sigma=None, digits="exact rationals",
         note="potential-space statement; amplitude readout cannot enter (stage B4)"),
    dict(q="C6 = 64 Lm/15pi", src="(7.3)(503)/map47",
         n=0, frac_sigma=None, digits="exact",
         note="energy bound constant, amplitude-free (stage B5)"),
    dict(q="mu^2 = m^2/(2 a_psi)", src="B.5(741)",
         n=0, frac_sigma=2e-6, digits="formula",
         note="replicated to 2e-6; convention-free (J1 gate <=1e-4 all window)"),
    dict(q="I.1 gate constants (virial 3e-4, INT b=1.000, eps=0.05)", src="I.1(777)",
         n=0, frac_sigma=None, digits="exact thresholds",
         note="dimensionless shape/energy gates; amplitude-free"),
    dict(q="p_dipole coefficient: p ~ p/mu^2", src="(7.4)(511)/(15.6)",
         n=1, frac_sigma=None, digits="NO coefficient printed",
         note="proportionality only: zero archive digits to test against"),
    dict(q="x0_lock = 2.42 +- 0.12, shift 0.50 +- 0.14", src="VII.D(525)",
         n=0, frac_sigma=0.14 / 0.50, digits="2 dp + band",
         note="ln2/slope shift within one law: p^2 prefactor cancels; also <r6> superseded (B.6)"),
    dict(q="C_d (pair-law amplitude)", src="(7.5)(521)",
         n=None, frac_sigma=None, digits="NEVER printed numerically",
         note="corpus prints the form only; our measurement enters via b_eff channel"),
]
for c in CATALOG:
    c["power_sigma"] = (power(c["n"], c["frac_sigma"])
                        if (c["n"] not in (None, 0) and c["frac_sigma"]) else
                        (0.0 if c["n"] == 0 else None))
RESULTS["stageA_catalog"] = CATALOG
flush()
for c in CATALOG:
    p = c["power_sigma"]
    ptxt = "n/a (stage B)" if p is None else f"{p:.3f} sigma"
    print(f"  {c['q']:55s} n={str(c['n']):>4s}  power={ptxt}")
print()

# ======================================================================
# STAGE B1 — the p channel (band + printed-digit granularity)
# ======================================================================
print("STAGE B1 — p channel")
split_abs = P15 - P6
pow_p = split_abs / P_SIG
# rebuild the grid summary from the archived per-point records.
# DEFECT (recorded): the on-disk j1_results.json is a PARTIAL flush -- it
# holds only the four c2 >= 0 points and no summary block, although
# j1_run.log and j1_RESULTS.md carry the full 9-point run.  K1 therefore
# (i) restricts grid arithmetic to the archived points (both 0.84 crossings
# lie in the archived c2 >= 0 range) and (ii) integrity-checks the benchmark
# values against j1_run.log.
cs = sorted(c for c in j1["spec"]["c2_grid"]
            if f"{c:+.3f}" in j1["points"])
pt2 = {c: j1["points"][f"{c:+.3f}"]["2N"] for c in cs}
RESULTS["defect_j1_json_partial"] = dict(
    archived_points=sorted(j1["points"].keys()),
    spec_grid=j1["spec"]["c2_grid"],
    note="j1_results.json on disk is a partial stage-flush (4/9 points, no "
         "summary) inconsistent with j1_run.log / j1_RESULTS.md; K1 uses the "
         "archived points and log cross-check only.")
logline = next(l for l in open(os.path.join(HERE, "j1_run.log"))
               if l.startswith("[c2 = +0.000]"))
p6_log = float(logline.split("p/6 =")[1].split()[0])

def crossings_of(conv, thr=P_TGT):
    vals = [pt2[c][conv] for c in cs]
    out = []
    for i in range(len(cs) - 1):
        a, b = vals[i] - thr, vals[i + 1] - thr
        if a == 0.0 or a * b < 0.0:
            out.append(cs[i] + (cs[i + 1] - cs[i]) * a / (a - b))
    return out

b1 = dict(
    p6=P6, p15=P15, split_abs=split_abs, split_sigma=pow_p,
    ratio_measured=P15 / P6, ratio_exact=SPLIT_RATIO,
    required_band_for_2sigma=split_abs / 2.0,
    crossings_c2={"p_edge_6": crossings_of("p_edge_6"),
                  "p_edge_15pi8": crossings_of("p_edge_15pi8")},
)
# slope dp/dc2 near benchmark (J1 grid, 2N, /6 convention) -> c2-offset equivalent
i0 = cs.index(0.0)
slope = (pt2[cs[i0 + 1]]["p_edge_6"] - pt2[cs[i0]]["p_edge_6"]) / (cs[i0 + 1] - cs[i0])
b1["dp_dc2_at_benchmark"] = slope
b1["c2_offset_equivalent"] = split_abs / abs(slope)
# printed-digit (rounding) channel: "0.84" = half-ulp window [0.835, 0.845)
lo, hi = 0.835, 0.845
b1["rounding_window"] = [lo, hi]
b1["p6_in_window"] = bool(lo <= P6 < hi)
b1["p15_in_window"] = bool(lo <= P15 < hi)
sig_round = 0.005 / math.sqrt(3.0)          # uniform half-ulp, sigma-equivalent
b1["conditional_sigma_equiv_p15_excluded"] = abs(P15 - P_TGT) / sig_round
b1["conditional_sigma_equiv_p6"] = abs(P6 - P_TGT) / sig_round
b1["conditions"] = [
    "benchmark potential is c2 = 0 exactly (corpus never prints its c2)",
    "printed 0.84 is the pipeline central value rounded to 2 dp, not the center "
    "of a c2-scan band (J1: +-0.03 corresponds to a c2 half-width ~0.07, so a "
    "scan-summary reading is live)",
    "cross-code amplitude agreement A_ours = A_frozen to << 0.5% (unverifiable: "
    "the frozen pipeline is not available -- this is the original OPEN itself)",
]
RESULTS["stageB1_p_channel"] = b1
flush()
gate("B1 input integrity: benchmark p/6 matches j1_run.log", abs(P6 - p6_log) < 5e-5,
     f"json {P6:.6f} vs log {p6_log:.4f}", "< 5e-5",
     "guards against the partial-flush defect recorded above")
gate("B1 p-band cannot split candidates (power < 2 sigma)", pow_p < 2.0,
     f"{pow_p:.3f} sigma", "< 2", f"needs +-{split_abs/2:.4f}, archive prints +-0.03")
gate("B1 rounding channel: p6 rounds to printed 0.84", b1["p6_in_window"],
     f"{P6:.4f} in [0.835,0.845)", "True", "CONDITIONAL channel, see conditions")
gate("B1 rounding channel: p15pi8 does NOT round to 0.84", not b1["p15_in_window"],
     f"{P15:.4f} in window = {b1['p15_in_window']}", "False",
     f"would print 0.86; {b1['conditional_sigma_equiv_p15_excluded']:.1f} sigma-equiv IF conditions held")
print()

# ======================================================================
# STAGE B2 — the b_eff channel (both candidates x 4 C_d x factor family)
# ======================================================================
print("STAGE B2 — b_eff channel")
mu0, R0 = D0["mu_fit"], D0["Rstar"]
env0 = math.exp(-mu0 * R0)
rows = []
max_dpull = 0.0
max_dpull_row = None
for est, (cd, cderr) in CD.items():
    for fn, fv in FACTORS:
        out = {}
        for tag, p in [("6", P6), ("15pi8", P15)]:
            b = cd * mu0 * env0 / (fv * p * p)
            sig = math.hypot(B_SIG, b * cderr / cd)
            out[tag] = (b, (b - B_TGT) / sig)
        # only factor choices that put at least one candidate near band matter
        if min(abs(out["6"][1]), abs(out["15pi8"][1])) > 3.0:
            continue
        dp = abs(out["6"][1] - out["15pi8"][1])
        row = dict(estimator=est, Fb=fn, Fb_value=fv,
                   b_6=out["6"][0], pull_6=out["6"][1],
                   b_15pi8=out["15pi8"][0], pull_15pi8=out["15pi8"][1],
                   delta_pull=dp)
        rows.append(row)
        if dp > max_dpull:
            max_dpull, max_dpull_row = dp, row
rows.sort(key=lambda r: -r["delta_pull"])
# F_b relabel compensation: can a named-factor swap absorb the candidate switch?
best_comp = None
for (n1, v1), (n2, v2) in combinations(FACTORS, 2):
    r = abs(math.log(v2 / v1))
    resid = abs(r - 2.0 * SPLIT_LN)      # b shifts by 2 ln(ratio) under the switch
    if best_comp is None or resid < best_comp[2]:
        best_comp = (n1, n2, resid)
comp_resid_sigma = best_comp[2] / (B_SIG / B_TGT)
b2 = dict(
    mu=mu0, Rstar=R0,
    n_rows_kept=len(rows), rows_top=rows[:12],
    max_delta_pull=max_dpull, max_delta_pull_row=max_dpull_row,
    required_total_frac_sigma_for_2sigma=(SPLIT_B - 1.0) / 2.0,
    archive_band_frac=B_SIG / B_TGT,
    cd_frac_errors={k: v[1] / v[0] for k, v in CD.items()},
    cd_cap_systematic_x54="cap-400 -> cap-1200 moved C_d x54 (unconverged; i3 caveat 2)",
    fb_relabel_best=dict(pair=[best_comp[0], best_comp[1]],
                         residual_ln=best_comp[2],
                         residual_sigma_of_band=comp_resid_sigma),
)
RESULTS["stageB2_b_channel"] = b2
flush()
print(f"  b-split exact: x{SPLIT_B:.5f}; archive band +-{100*B_SIG/B_TGT:.1f}%; "
      f"needs +-{100*(SPLIT_B-1)/2:.2f}% total for 2 sigma")
for r in rows[:6]:
    print(f"  {r['estimator']:38s} Fb={r['Fb']:>6s}  b6={r['b_6']:7.2f} ({r['pull_6']:+.2f})  "
          f"b15={r['b_15pi8']:7.2f} ({r['pull_15pi8']:+.2f})  d_pull={r['delta_pull']:.3f}")
gate("B2 b_eff cannot split candidates (max delta-pull < 2 sigma)",
     max_dpull < 2.0, f"{max_dpull:.3f} sigma max over C_d x F_b", "< 2",
     "on top of the unpinned F_b identity and the x54 C_d systematic")
gate("B2 F_b relabel absorbs the switch below band resolution",
     comp_resid_sigma < 1.0,
     f"pair ({best_comp[0]},{best_comp[1]}) residual {best_comp[2]:.4f} ln "
     f"= {comp_resid_sigma:.2f} sigma_band", "< 1",
     "structural degeneracy: factor identity itself unmeasured (i3 caveat 1)")
# b_inf ~ 4.5: three independent blockers, printed as measured facts
binf_blockers = dict(
    no_error_printed=True,
    far_zone_instrument=("archived running b_eff at largest x (extrap array): "
                         + ", ".join(f"{v:.2e}" for v in BEFF_EXTRAP)
                         + "  (sign-indefinite: far zone unmeasured)"),
    fb_unpinned="split 3.75% < factor-family granularity 6.7% median (i3)",
)
RESULTS["stageB2_b_inf"] = binf_blockers
flush()
gate("B2 b_inf~4.5 unusable (no band + far zone unmeasured + F_b unpinned)",
     True, "3 independent blockers recorded", "informational",
     "even a +-0.05 ulp-reading would be voided by the F_b degeneracy")
print()

# ======================================================================
# STAGE B3 — the bond-equation x0 loop (reading-B crossing sensitivity)
# ======================================================================
print("STAGE B3 — x0 bond-equation loop channel")

def crossing(xs, beff, thr):
    for i in range(len(xs) - 1):
        if beff[i] < thr <= beff[i + 1]:
            f = (thr - beff[i]) / (beff[i + 1] - beff[i])
            return xs[i] + f * (xs[i + 1] - xs[i])
    return None

b3_rows = []
max_dx0 = 0.0
for fbn, fbv in [("2", 2.0), ("4pi^2", 4 * math.pi ** 2), ("8pi^2", 8 * math.pi ** 2)]:
    for tag, xs, arr in [("seed", X_GRID, BEFF_SEED),
                         ("relaxed", X_GRID, BEFF_RELAX),
                         ("extrap", X_GRID_EXTRAP, BEFF_EXTRAP)]:
        res = {}
        for pn, p in [("6", P6), ("15pi8", P15)]:
            scale = ENV * F_OLD * P_OLD ** 2 / (fbv * p * p)
            res[pn] = crossing(xs, [v * scale for v in arr], B_TGT)
        if res["6"] is not None and res["15pi8"] is not None:
            dx0 = abs(res["6"] - res["15pi8"])
            max_dx0 = max(max_dx0, dx0)
        else:
            dx0 = None
        b3_rows.append(dict(Fb=fbn, array=tag, x0_6=res["6"],
                            x0_15pi8=res["15pi8"], delta_x0=dx0))
pow_x0 = max_dx0 / X0_LOOP_SIG
RESULTS["stageB3_x0_loop"] = dict(
    rows=b3_rows, max_delta_x0=max_dx0, power_sigma=pow_x0,
    note="x0(measured)=1.92+-0.08 is amplitude-free by construction "
         "(well location of E_int(x); no p enters the energy differences)")
flush()
for r in b3_rows:
    t = (f"{r['delta_x0']:.4f}" if r["delta_x0"] is not None else "n/a")
    x6 = f"{r['x0_6']:.4f}" if r["x0_6"] else "none"
    x15 = f"{r['x0_15pi8']:.4f}" if r["x0_15pi8"] else "none"
    print(f"  Fb={r['Fb']:>6s} {r['array']:8s} x0(/6)={x6:>7s} x0(/15pi8)={x15:>7s}  dx0={t}")
gate("B3 bond-equation x0 loop cannot split (power < 2 sigma)",
     pow_x0 < 2.0, f"max dx0 = {max_dx0:.4f} = {pow_x0:.3f} sigma of +-0.05", "< 2",
     "the crossing sits in the steep rise: 3.75% in b moves x0 by ~0.01")
print()

# ======================================================================
# STAGE B4 — the B.4 window edges re-derived; amplitude-independence
# ======================================================================
print("STAGE B4 — B.4 window edges [-m^2/2, m^2/6)")

def V(f, c2):
    x = 1.0 - np.cos(f)
    return x + c2 * x * x          # m = 1

# lower edge: positivity of V on f in [0, pi]; V(pi) = 2 + 4 c2
fgrid = np.linspace(0.0, np.pi, 200001)
vmin_at_edge = float(np.min(V(fgrid, -0.5)))
vmin_below = float(np.min(V(fgrid, -0.5 - 1e-4)))
gate("B4 lower edge = positivity: min V(c2=-1/2) = 0 (at f=pi)",
     abs(vmin_at_edge) < 1e-12 and abs(V(np.array([np.pi]), -0.5)[0]) < 1e-12,
     f"min V = {vmin_at_edge:.2e}", "= 0", "V(pi) = 2 + 4 c2 -> c2 >= -1/2 exact")
gate("B4 lower edge sharp: V < 0 just below", vmin_below < 0,
     f"min V(c2=-1/2-1e-4) = {vmin_below:.2e}", "< 0")

# upper edge: quartic Taylor coefficient of V about the vacuum f = 0
# analytic: V = f^2/2 + (c2/4 - 1/24) f^4 + O(f^6) -> zero at c2 = 1/6
def a4_numeric(c2):
    ff = np.linspace(1e-3, 0.15, 400)
    y = V(ff, c2) - 0.5 * ff ** 2
    M = np.vstack([ff ** 4, ff ** 6, ff ** 8, ff ** 10]).T
    coef, *_ = np.linalg.lstsq(M, y, rcond=None)
    return float(coef[0])

lo_c, hi_c = 0.0, 0.4
for _ in range(80):                      # bisection on the numeric a4
    mid = 0.5 * (lo_c + hi_c)
    if a4_numeric(mid) < 0.0:
        lo_c = mid
    else:
        hi_c = mid
c2_root_numeric = 0.5 * (lo_c + hi_c)
a4_analytic_err = max(abs(a4_numeric(c2) - (c2 / 4.0 - 1.0 / 24.0))
                      for c2 in [-0.5, -0.2, 0.0, 1.0 / 6.0, 0.3])
gate("B4 upper edge = quartic-coefficient sign change at c2 = 1/6 exactly",
     abs(c2_root_numeric - 1.0 / 6.0) < 1e-9,
     f"numeric root {c2_root_numeric:.12f}", "1/6 = 0.166666666667",
     "a4(c2) = c2/4 - 1/24; matches B.4's 'amplitude-mode coupling' language")
gate("B4 numeric a4 matches analytic c2/4 - 1/24",
     a4_analytic_err < 1e-9, f"max |err| = {a4_analytic_err:.2e}", "< 1e-9",
     "fit-window truncation limited; 7+ orders below any physical scale here")
RESULTS["stageB4_window_edges"] = dict(
    lower_edge=dict(statement="V(pi) = 2 m^2 + 4 c2 >= 0", exact_edge=-0.5,
                    vmin_at_edge=vmin_at_edge, vmin_below=vmin_below),
    upper_edge=dict(statement="quartic Taylor coeff a4 = c2/4 - m^2/24 crosses 0",
                    exact_edge=1.0 / 6.0, numeric_root=c2_root_numeric,
                    analytic_check_max_err=a4_analytic_err),
    amplitude_independence="both edge conditions are functionals of V(f; m, c2) "
        "alone; the dimensionless amplitude p enters ONLY the exterior tail "
        "readout f -> p(1+mu r)e^{-mu r}/(mu r)^2 mu (VII.C). Sensitivity to "
        "F_p is identically zero: the edges re-derive bit-identically under "
        "either candidate because the computation takes no amplitude input.",
)
flush()
gate("B4 edges amplitude-free (sensitivity exactly 0)", True,
     "edge computation consumes (m, c2) only", "structural",
     "re-derives identically under BOTH candidates: zero discrimination power")
print()

# ======================================================================
# STAGE B5 — C6 Haar constant and mu formula (amplitude-free anchors)
# ======================================================================
print("STAGE B5 — C6 and mu (amplitude-free anchors, re-verified)")
x, w = np.polynomial.legendre.leggauss(2000)
chi = 0.5 * np.pi * (x + 1.0)
wq = 0.5 * np.pi * w
c6 = float((2.0 / np.pi) * np.sum(wq * 2.0 * np.sin(chi / 2.0) * np.sin(chi) ** 2))
c6_exact = 64.0 / (15.0 * np.pi)
gate("B5 C6 Haar quadrature = 64/(15pi)", abs(c6 - c6_exact) < 1e-12,
     f"{c6:.15f}", f"{c6_exact:.15f}",
     "energy-bound constant: no amplitude enters; 15pi here is NOT the tail F_p")
mu_errs = [j1["points"][k]["2N"]["mu_rel_err"] for k in j1["points"]]
gate("B5 mu = m/sqrt(2 a_psi) convention-free across window",
     max(mu_errs) <= 1e-4,
     f"max rel err {max(mu_errs):.2e} ({len(mu_errs)} archived pts; "
     "9/9 <= 1e-4 in j1_run.log)", "<= 1e-4")
RESULTS["stageB5_anchors"] = dict(C6_quadrature=c6, C6_exact=c6_exact,
                                  mu_rel_err_max=max(mu_errs))
flush()
print()

# ======================================================================
# STAGE B6 — frozen-record content scan (can A_3D be re-fit? NO)
# ======================================================================
print("STAGE B6 — two-knot frozen records: is a 3-D tail-amplitude fit possible?")
rs = json.load(open(RUNSINGLE))

def max_array_len(o):
    if isinstance(o, dict):
        return max([0] + [max_array_len(v) for v in o.values()])
    if isinstance(o, list):
        inner = max([0] + [max_array_len(v) for v in o
                           if isinstance(v, (dict, list))])
        return max(len(o), inner)
    return 0

mal = max_array_len(rs)
ncells = rs["protocol"]["nx"] * rs["protocol"]["nt"] ** 2
gate("B6 no field snapshot in frozen run records",
     mal < 1000 and mal < ncells,
     f"largest array in run_single.json = {mal} entries", f"<< {ncells} cells",
     "scalar diagnostics + iteration series only: the 3-D tail amplitude "
     "cannot be re-fit from the archive; its amplitude info enters via C_d (B2)")
RESULTS["stageB6_record_scan"] = dict(
    max_array_len=mal, n_cells=ncells,
    boundary_tail_scalar=rs["final"]["boundary_tail"],
    conclusion="A_3D not re-fittable from frozen records; re-running the 3-D "
               "solver would add OUR precision, not ARCHIVE precision -- the "
               "bottleneck is the printed +-0.03/+-6 bands, which no re-run moves.")
flush()
print()

# ======================================================================
# STAGE C — verdict
# ======================================================================
print("STAGE C — verdict")
powers = dict(p_band=pow_p, b_eff=max_dpull, x0_loop=pow_x0,
              x0_measured=0.0, window_edges=0.0, C6=0.0, mu=0.0,
              gate_constants=0.0, dipole_coefficient=None, b_inf=None)
max_uncond = max(v for v in powers.values() if v is not None)
verified_null = max_uncond < 2.0
gate("C VERIFIED NULL: no archive quantity separates candidates at >= 2 sigma",
     verified_null, f"max unconditional power = {max_uncond:.3f} sigma (p band)",
     "< 2", "conditional rounding pointer favors /6 but rests on unprinted assumptions")
RESULTS["stageC_verdict"] = dict(
    powers_sigma=powers,
    max_unconditional_power=max_uncond,
    verified_null=verified_null,
    verdict=("INDISTINGUISHABLE WITHIN THE ARCHIVE (verified null). Best "
             "unconditional channel is the p band itself at "
             f"{pow_p:.2f} sigma; every other amplitude-bearing quantity is "
             "either amplitude-free (n = 0), unbanded, degenerate with the "
             "unpinned F_b identity, or orders of magnitude short in "
             "precision. The CONDITIONAL printed-digit channel (0.84 rounds "
             "from 0.8434, not 0.8591) favors /6 at 6.6 sigma-equivalent IF "
             "three unprinted assumptions hold; it is not archive-backed and "
             "does not meet the >= 2 sigma unconditional bar."),
    prior_figure_note=("J1/ROADMAP carried the split as '0.7 sigma'; the "
                       f"measured split is {pow_p:.2f} sigma "
                       f"({split_abs:.4f}/0.03) - the prior figure mildly "
                       "OVERSTATES the split; the null is stronger than "
                       "advertised."),
)
RESULTS["gates"] = GATES
RESULTS["all_gates_pass"] = all(g["ok"] for g in GATES)
flush()
print(f"\nALL GATES: {'PASS' if RESULTS['all_gates_pass'] else 'FAIL'}"
      f"  ({sum(g['ok'] for g in GATES)}/{len(GATES)})")
print("wrote", OUT)
