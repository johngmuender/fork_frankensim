#!/usr/bin/env python3
"""
Phase I3 — the joint convention solve (p x b_eff).

Two spec-underdeterminations survived the Tier-2b/G4-G5b campaigns:

  p     = 0.84 +- 0.03   (App B.4 / I.1 radial tail amplitude; nearest
                          single-constraint candidate A e^{-mu R*}/2pi = 0.805,
                          ~1.2 sigma low; five natural candidates 0.365-5.06)
  b_eff = 42 +- 6        (App I.2 bond constant; shape closed in band,
                          ABSOLUTE normalization ~7 orders unrecovered under
                          the App-A reading b = C_d mu / (2 pi p^2), p = 0.84)

They are COUPLED: b's formula contains p^2.  This script enumerates the
normalization-convention space programmatically —

  p_cand = base * exp(-gp * mu R*) * mu^ap * R*^bp * Fp^sp
           base in {A (fitted k1 amplitude), f_fit(R*) (edge value)}
  b_cand = C_d * mu^cb * R*^eb * exp(-gb * mu R*) / (Fb^sb * p_cand^2)
           C_d in {G4 relaxed-anchored corner, G4 seed corner,
                   G4 seed central4, G5b cap-1200 relaxed}

with F drawn from the named 2^i 3^j pi^k factor family (1 .. ~4pi^2), and
scores every assignment jointly against (0.84 +- 0.03, 42 +- 6).  For the
top assignments it re-checks the bond-equation reading-B closed loop
(the running b_eff(x) crossing +42, corpus x0 = 1.90 +- 0.05) under the
re-scaled convention, using the archived running-estimator arrays.

All inputs are read from the frozen campaign records:
  tier2-closure/radial_results.json                (A, mu, R*, f(R*))
  gum-core/fs-gum-twoknot/twoknot_results.json     (C_d fits, running b_eff)
  gum-core/fs-gum-twoknot/twoknot_refine_results.json (G5b C_d, extrap b_eff)

Deterministic; ~1 s.  Writes i3_results.json.  No repo files are modified.
"""

import json
import math
import os
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
SBOX = os.path.dirname(HERE)  # analysis/gum-sandbox

RADIAL = os.path.join(SBOX, "tier2-closure", "radial_results.json")
TWOKNOT = os.path.join(SBOX, "gum-core", "fs-gum-twoknot", "twoknot_results.json")
REFINE = os.path.join(SBOX, "gum-core", "fs-gum-twoknot", "twoknot_refine_results.json")

# ---------------------------------------------------------------- targets
P_TGT, P_SIG = 0.84, 0.03          # corpus <r1> tail amplitude (App B.4/VII.C)
B_TGT, B_SIG = 42.0, 6.0           # corpus bond constant (App I.2 / VII.D)
X0_BOND, X0_BOND_SIG = 1.90, 0.05  # corpus bond equation at b = 42

# ---------------------------------------------------------------- inputs
rad = json.load(open(RADIAL))["2N_run"]
two = json.load(open(TWOKNOT))
ref = json.load(open(REFINE))

A = rad["A_fit"]                          # 5.2333e6 (k1-normalized tail amp)
MU = rad["mu_analytic_derived"]           # 7.77254686 (frozen referee mu)
RSTAR = json.load(open(RADIAL))["RSTAR"]  # 2^(5/6)
FRS = rad["f_fit_at_Rstar"]               # 0.39177 (fitted f at the edge)
MURS = MU * RSTAR                         # 13.8491 — THE large number
ENV = math.exp(-MURS)                     # 9.670e-7 — the edge envelope

CD_ESTIMATORS = {
    # name: (C_d, sigma_C_d, status note)
    "G4-relaxed-corner": (
        two["fits"]["corner"]["relaxed_anchored"]["C_d"],
        two["fits"]["corner"]["relaxed_anchored"]["C_d_err"],
        "cap-400; G5b showed under-converged (moved 54x up at cap 1200)",
    ),
    "G4-seed-corner": (
        two["fits"]["corner"]["seed_ansatz"]["C_d"],
        two["fits"]["corner"]["seed_ansatz"]["C_d_err"],
        "iteration-free product ansatz; O(overlap^2) bias in window",
    ),
    "G4-seed-central4": (
        two["fits"]["central4"]["seed_ansatz"]["C_d"],
        two["fits"]["central4"]["seed_ansatz"]["C_d_err"],
        "4th-order instrument on seed states (cross-check)",
    ),
    "G5b-relaxed-1200": (
        ref["fit_hi"]["C_d"],
        math.sqrt(ref["fit_hi"]["cov_CB"][0][0]),
        "deepest relaxation available; dof = 1",
    ),
}

# running b_eff(x) arrays as archived (computed with mu/(2 pi 0.84^2) norm)
X_GRID = [1.4732, 1.6837, 1.8942, 2.1046, 2.4203, 2.7360, 3.0517]
BEFF_SEED = two["bond"]["corner"]["beff_running_seed"]
BEFF_RELAX = two["bond"]["corner"]["beff_running_relaxed"]
X_GRID_EXTRAP = [1.6837, 1.8942, 2.1046, 2.4203]
BEFF_EXTRAP = ref["bond"]["beff_running_extrap"]
P_OLD, F_OLD = 0.84, 2.0 * math.pi  # normalization the arrays were built with

# ---------------------------------------------------------------- factor family
def factor_family():
    """Named factors m * pi^k, m in a small integer/rational/surd set."""
    ms = [
        ("1", 1.0), ("sqrt2", math.sqrt(2)), ("3/2", 1.5), ("2", 2.0),
        ("2sqrt2", 2 * math.sqrt(2)), ("3", 3.0), ("4", 4.0), ("6", 6.0),
        ("8", 8.0), ("16/15", 16 / 15), ("4/3", 4 / 3), ("15/8", 15 / 8),
        ("16", 16.0),
    ]
    out = []
    for mn, mv in ms:
        for k, kn in [(0, ""), (1, "pi"), (2, "pi^2")]:
            v = mv * math.pi ** k
            name = kn if mn == "1" and k > 0 else (mn if k == 0 else f"{mn}{kn}")
            out.append((name or "1", v))
    # dedup by value (keep first/simplest name)
    seen, ded = set(), []
    for n, v in sorted(out, key=lambda t: (round(math.log(t[1]), 9), len(t[0]))):
        key = round(math.log(v), 9)
        if key not in seen:
            seen.add(key)
            ded.append((n, v))
    return ded

FACTORS = factor_family()

# ---------------------------------------------------------------- p candidates
def enumerate_p(pull_keep=3.0):
    cands = []
    bases = [("A", A, False), ("f(R*)", FRS, True)]  # f(R*) already edge-valued
    for (bn, bv, edge_built_in), gp, ap, bp, (fn, fv), sp in product(
        bases, (0, 1), (-2, -1, 0, 1), (-1, 0, 1), FACTORS, (1, -1)
    ):
        if edge_built_in and gp != 0:
            continue  # f(R*) already contains the envelope physically
        val = bv * math.exp(-gp * MURS) * MU ** ap * RSTAR ** bp / (fv ** sp)
        if not (0.1 < val < 8.0):
            continue
        pull = (val - P_TGT) / P_SIG
        if abs(pull) > pull_keep:
            continue
        op = "/" if sp == 1 else "*"
        desc = bn
        if gp:
            desc += " e^{-muR*}"
        if ap:
            desc += f" mu^{ap}"
        if bp:
            desc += f" R*^{bp}"
        if fv != 1.0:
            desc += f" {op} {fn}"
        # complexity: departures from a bare amplitude statement
        comp = abs(ap) + abs(bp) + (0 if fv == 1.0 else 1) + (0 if sp == 1 else 1)
        cands.append(dict(desc=desc, value=val, pull=pull, comp=comp,
                          base=bn, gp=(1 if (gp or edge_built_in) else 0),
                          ap=ap, bp=bp, F=fn, s=sp))
    # dedup identical values, keep lowest complexity
    ded = {}
    for c in cands:
        key = round(math.log(c["value"]), 9)
        if key not in ded or (c["comp"], len(c["desc"])) < (ded[key]["comp"], len(ded[key]["desc"])):
            ded[key] = c
    return sorted(ded.values(), key=lambda c: (abs(c["pull"]), c["comp"]))

# ---------------------------------------------------------------- b assignments
def enumerate_joint(p_cands, p_gate=2.0, keep_maxpull=3.0):
    rows = []
    for pc in p_cands:
        if abs(pc["pull"]) > p_gate:
            continue
        p2 = pc["value"] ** 2
        for (est, (cd, cderr, note)), gb, cb, eb, (fn, fv), sb in product(
            CD_ESTIMATORS.items(), (0, 1, 2), (-1, 0, 1), (-1, 0, 1),
            FACTORS, (1, -1)
        ):
            b = cd * MU ** cb * RSTAR ** eb * math.exp(-gb * MURS) / (fv ** sb * p2)
            if not (1.0 < b < 2000.0):
                continue
            sig_eff = math.hypot(B_SIG, b * cderr / cd)
            pull_b = (b - B_TGT) / sig_eff
            if abs(pull_b) > keep_maxpull:
                continue
            op = "/" if sb == 1 else "*"
            desc = f"C_d[{est}]"
            if cb:
                desc += f" mu^{cb}"
            if eb:
                desc += f" R*^{eb}"
            if gb:
                desc += f" e^{{-{gb}muR*}}"
            desc += f" {op} ({fn} p^2)" if fv != 1.0 else " / p^2"
            comp_b = abs(cb - 1) + abs(eb) + (0 if abs(fv - 2 * math.pi) < 1e-12 else 1) \
                + (0 if sb == 1 else 1)
            maxpull = max(abs(pc["pull"]), abs(pull_b))
            rows.append(dict(
                p_desc=pc["desc"], p_value=pc["value"], p_pull=pc["pull"],
                b_desc=desc, b_value=b, b_pull=pull_b, estimator=est,
                gb=gb, cb=cb, eb=eb, Fb=fn, Fb_value=fv, sb=sb,
                comp=pc["comp"] + comp_b, maxpull=maxpull,
                joint_rss=math.hypot(pc["pull"], pull_b),
            ))
    # rank economy-first: with a factor family this dense, thousands of
    # assignments land inside 1.5 sigma — raw pull-ranking rewards numerology.
    rows.sort(key=lambda r: (r["comp"], r["maxpull"], r["b_desc"]))
    return rows

# ---------------------------------------------------------------- gb scan (the discrete pin)
def envelope_requirement():
    """For each estimator and each envelope count gb, the residual factor
    F_req = C_d mu e^{-gb muR*} / (42 p^2) that the convention family must
    supply (baseline cb=1, eb=0, p = 0.84).  Shows gb = 1 is forced."""
    out = {}
    p2 = P_TGT ** 2
    for est, (cd, _, _) in CD_ESTIMATORS.items():
        out[est] = {
            f"gb={gb}": cd * MU * math.exp(-gb * MURS) / (B_TGT * p2)
            for gb in (0, 1, 2)
        }
    return out

# ---------------------------------------------------------------- canonical readings
def canonical_readings():
    """The STRICT App-A formula b = B/(2 pi p^2) with B = C_d mu e^{-gb muR*}
    (gb = 1 forced by the envelope-requirement scan), scored for the two
    lowest-complexity in-band p candidates.  Shows whether the literal 2pi
    normalization closes without any extra factor."""
    out = []
    for pname, pval in [("A e^{-muR*} / 2pi", A * ENV / (2 * math.pi)),
                        ("A e^{-muR*} / 6", A * ENV / 6.0)]:
        for est, (cd, cderr, _) in CD_ESTIMATORS.items():
            b = cd * MU * ENV / (2 * math.pi * pval ** 2)
            sig = math.hypot(B_SIG, b * cderr / cd)
            out.append(dict(p_desc=pname, p_value=pval,
                            p_pull=(pval - P_TGT) / P_SIG,
                            estimator=est, b_value=b,
                            b_pull=(b - B_TGT) / sig))
    return out

# ---------------------------------------------------------------- reading-B validation
def crossing(xs, beff, thr):
    """First upward crossing of thr (linear interp), scanning in x."""
    for i in range(len(xs) - 1):
        if beff[i] < thr <= beff[i + 1]:
            f = (thr - beff[i]) / (beff[i + 1] - beff[i])
            return xs[i] + f * (xs[i + 1] - xs[i])
    return None

def reading_b_check(row):
    """Rescale the archived running b_eff(x) arrays to the candidate
    convention and locate the +42 crossing (bond-equation reading B)."""
    scale = (MU ** (row["cb"] - 1) * RSTAR ** row["eb"]
             * math.exp(-row["gb"] * MURS)
             * F_OLD * P_OLD ** 2 / (row["Fb_value"] ** row["sb"] * row["p_value"] ** 2))
    res = {}
    for tag, xs, arr in [("seed", X_GRID, BEFF_SEED),
                         ("relaxed", X_GRID, BEFF_RELAX),
                         ("extrap", X_GRID_EXTRAP, BEFF_EXTRAP)]:
        x0 = crossing(xs, [v * scale for v in arr], B_TGT)
        res[tag] = dict(
            x0=x0,
            pull=(None if x0 is None else (x0 - X0_BOND) / X0_BOND_SIG),
        )
    res["scale_vs_archived"] = scale
    return res

# ---------------------------------------------------------------- run
def main():
    p_cands = enumerate_p()
    joint = enumerate_joint(p_cands)

    n15 = sum(1 for r in joint if r["maxpull"] <= 1.5)
    n10 = sum(1 for r in joint if r["maxpull"] <= 1.0)
    envreq = envelope_requirement()

    # look-elsewhere density: factor-family multiplicative granularity
    lv = sorted(math.log(v) for _, v in FACTORS)
    gaps = [b - a for a, b in zip(lv, lv[1:])]
    med_gap = sorted(gaps)[len(gaps) // 2]

    top = joint[:24]
    # best (lowest-comp, then lowest-pull) assignment per C_d estimator
    best_per_est = {}
    for r in joint:
        if r["maxpull"] <= 1.5 and r["estimator"] not in best_per_est:
            best_per_est[r["estimator"]] = r
    for r in top + [r for r in best_per_est.values() if r not in top]:
        r["reading_B"] = reading_b_check(r)
    canon = canonical_readings()

    # ---- console report
    print("=" * 100)
    print("I3 — joint convention solve: p = 0.84 +- 0.03  x  b_eff = 42 +- 6")
    print("=" * 100)
    print(f"inputs: A = {A:.6g}, mu = {MU:.8g}, R* = {RSTAR:.8g}, "
          f"mu R* = {MURS:.4f}, e^-muR* = {ENV:.4e}, f(R*) = {FRS:.5f}")
    print(f"C_d estimators: " + ", ".join(
        f"{k} = {v[0]:.3e} ({100*v[1]/v[0]:.0f}%)" for k, v in CD_ESTIMATORS.items()))
    print()
    print("-- envelope requirement (residual factor needed at baseline b = C_d mu e^{-gb muR*}/(42 p^2), p=0.84):")
    for est, d in envreq.items():
        print(f"   {est:20s} " + "  ".join(f"{k}: {v:11.4e}" for k, v in d.items()))
    print()
    print(f"-- p candidates within 3 sigma of 0.84 (of the full enumeration): {len(p_cands)}")
    for c in p_cands[:12]:
        print(f"   p = {c['value']:.4f}  pull {c['pull']:+.2f}  comp {c['comp']}   {c['desc']}")
    print()
    print(f"-- joint assignments: {len(joint)} kept (|pull_b| <= 3, |pull_p| <= 2); "
          f"{n15} land BOTH within 1.5 sigma, {n10} within 1.0 sigma")
    print(f"-- factor-family granularity: median multiplicative gap "
          f"{math.exp(med_gap):.3f}x -> hits of this multiplicity are EXPECTED (see caveats)")
    print()
    print("-- STRICT App-A canonical reading b = C_d mu e^{-muR*}/(2pi p^2) "
          "(gb = 1 forced; no extra factor allowed):")
    for c in canon:
        print(f"   p = {c['p_desc']:18s} ({c['p_value']:.4f}, {c['p_pull']:+.2f}s)  "
              f"{c['estimator']:20s} b = {c['b_value']:8.2f}  pull {c['b_pull']:+6.2f}")
    print()

    def row_txt(r):
        rb = r["reading_B"]
        xs_, xe_ = rb["seed"]["x0"], rb["extrap"]["x0"]
        xtxt = (f"{xs_:.2f}" if xs_ else "none") + "/" + (f"{xe_:.2f}" if xe_ else "none")
        return (f"{r['p_desc']:26s} {r['p_value']:6.3f} {r['p_pull']:+5.2f}  "
                f"{r['b_desc']:56s} {r['b_value']:7.2f} {r['b_pull']:+5.2f} "
                f"{r['comp']:2d} {xtxt:>10s}")

    hdr = (f"{'p convention':26s} {'p':>6s} {'pull':>5s}  "
           f"{'b convention':56s} {'b':>7s} {'pull':>5s} {'cx':>2s} {'x0(B) s/e':>10s}")
    print("-- joint table, economy-ranked (cx = convention complexity; "
          "x0(B) = reading-B 42-crossing, seed/extrap arrays, corpus 1.90 +- 0.05):")
    print(hdr)
    print("-" * len(hdr))
    for r in top:
        print(row_txt(r))
    print()
    print("-- most economical in-band (<= 1.5 sigma both) assignment per C_d estimator:")
    for est, r in best_per_est.items():
        print(row_txt(r))

    out = dict(
        inputs=dict(A=A, mu=MU, Rstar=RSTAR, muRstar=MURS, env=ENV, f_at_Rstar=FRS,
                    C_d_estimators={k: dict(C_d=v[0], sigma=v[1], note=v[2])
                                    for k, v in CD_ESTIMATORS.items()}),
        targets=dict(p=[P_TGT, P_SIG], b_eff=[B_TGT, B_SIG],
                     x0_bond_eq=[X0_BOND, X0_BOND_SIG]),
        envelope_requirement=envreq,
        factor_family=[[n, v] for n, v in FACTORS],
        factor_family_median_gap=math.exp(med_gap),
        p_candidates_top=[{k: c[k] for k in ("desc", "value", "pull", "comp")}
                          for c in p_cands[:15]],
        n_joint_kept=len(joint),
        n_joint_within_1p5=n15,
        n_joint_within_1p0=n10,
        canonical_appA_reading=canon,
        joint_top=top,
        best_per_estimator=best_per_est,
    )
    with open(os.path.join(HERE, "i3_results.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(HERE, 'i3_results.json')}")

if __name__ == "__main__":
    main()
