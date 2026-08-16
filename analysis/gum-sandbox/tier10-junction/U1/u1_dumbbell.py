#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
U1 -- Junction-antijunction (J-Jbar) dumbbell spectrum from printed
corpus parameters.  [F-T10-U1]
ROADMAP_v12_JUNCTION.md workstream U1 (frozen spec; gates pre-registered).

Deterministic, no RNG.  Extends the Tier-9 T1 machinery (variable-mass
radial Schroedinger solver; frozen cutoff f* = 1.9724 from T1's ONE
validation calibration; Zhu-Kroemer ordering central, BenDaniel-Duke
systematic) to the minimal knot-free junction state: two order-3
junctions joined by three tubes (STAR's "baryonium glueball" class,
Csorgo-Gyulassy-Kharzeev 2004; Frenklakh-Kharzeev-Rossi-Veneziano
JHEP 07 (2024) 262 intercept ~1/2 [IM]).

  Dumbbell configuration coordinate: junction separation L.
  Dimensionless units x = sqrt(sigma) L, eps = E/sqrt(sigma),
  mhat = m_J/sqrt(sigma).

    V_N(x)  = 3x + 2*mhat + (pi*N - pi/4) * (1 - exp(-f x)) / x
    mu(x)   = 3x + 2*mhat          (configuration inertia = tube mass
                                    3*sigma*L + junction masses 2*m_J,
                                    imitating T1's mu(rho) = 2*pi*sigma*rho)

  3x           : three parallel tubes, energy 3*sigma*L at leading order
                 (Y-law end geometry collapses to parallel tubes on the
                 J-Jbar axis at leading order -- AMENDMENT U1-A1).
  pi*N/x       : transverse phonons of an open tube with both ends fixed
                 on junctions (Dirichlet), omega_n = n*pi/L; N = sum of
                 mode numbers over 3 tubes x 2 transverse polarizations.
  -pi/4 /x     : zeta-regularized zero point, 6 polarization channels x
                 (pi/L)*zeta(-1)/... : per channel (1/2)sum n*pi/L ->
                 (pi/2L)*zeta(-1) = -pi/(24 L); total -pi/(4L) = the
                 three-tube open-string Luscher term -3*(D-2)*pi/(24L).
  (1-e^{-fx})  : T1's short-distance smoothing transferred as
                 junction-end smoothing at the FROZEN f* (U1-A3).

Junction mass m_J: explicit parameter, NOT fitted; evaluated at the three
pre-registered values mhat in {0, 0.1355 [IM 2+1D lattice anchor,
arXiv:2508.00608 -- dimensional-transfer caveat], 0.39 = f_q central
0.17 GeV / sqrt(sigma) 0.436 GeV}.

VALIDATION FIRST (U1-G1): (a) two-body linear potential vs exact Airy
zeros to <= 1e-6 relative; (b) the same generic solver fed T1's loop
potential/inertia must reproduce T1's published eps(0++) = 3.5827 (and
0-+ 8.2749) at the frozen f* to <= 1e-3 relative.

Epistemic register (binding): within-model; nothing here bears on nature.
Seal q-theta: every mass number confronting a hadron mass is reportable
ONLY at V.F grade (dimensionally secure / structurally plausible /
quantitatively unclaimed).  STAR observables (1.84, 0.64, 1.04) are NOT
used anywhere in this workstream -- no fitting to STAR, no junction
observable enters any calibration.  Zero fitted parameters: f* frozen by
T1 BEFORE Tier 10 existed.
"""

import json
import math
import os
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import ai_zeros

HERE = os.path.dirname(os.path.abspath(__file__))
T1_JSON = os.path.normpath(os.path.join(
    HERE, "..", "..", "tier9-glueball", "T1", "t1_results.json"))

SIGMA_GUM = 0.19                      # GeV^2 [IM]; sqrt = 0.43589 GeV
SQS = math.sqrt(SIGMA_GUM)
MHAT_SET = [0.0, 0.1355, 0.39]        # pre-registered; NOT fitted
ZPD = -math.pi / 4.0                  # 6 channels x zeta-reg (-pi/24)

RESULTS = {"workstream": "U1", "finding": "F-T10-U1",
           "register": ("within-model; V.F-graded (dimensionally secure / "
                        "structurally plausible / quantitatively unclaimed); "
                        "nothing bears on nature; STAR observables never "
                        "used as targets"),
           "amendments": [], "self_checks": {}, "validation": {},
           "spectrum": {}, "confrontation": {}}


def amend(tag, text):
    line = f"[AMENDMENT {tag}] {text}"
    print(line)
    RESULTS["amendments"].append(line)


# =====================================================================
# Generic variable-mass radial Schroedinger solver (T1 algorithm,
# generalized to arbitrary V(x), mu(x); Dirichlet at 0 and xmax).
#   zk : H = mu^{-1/2} (-1/2 d^2/dx^2) mu^{-1/2} + V   (Zhu-Kroemer)
#   bdd: H = -1/2 d/dx [(1/mu) d/dx] + V               (BenDaniel-Duke)
# =====================================================================
def eigen_vm(Vfun, mufun, n=0, ordering="zk", npts=16000, xmax=18.0,
             nmax=None):
    h = xmax / (npts + 1)
    x = h * np.arange(1, npts + 1)
    V = Vfun(x)
    top = max(n, 0) if nmax is None else max(nmax, n)
    if ordering == "zk":
        d = 1.0 / np.sqrt(mufun(x))
        diag = d * d / (h * h) + V
        off = -(d[:-1] * d[1:]) / (2.0 * h * h)
    elif ordering == "bdd":
        a_l = 1.0 / mufun(x - h / 2.0)
        a_r = 1.0 / mufun(x + h / 2.0)
        diag = (a_l + a_r) / (2.0 * h * h) + V
        off = -a_r[:-1] / (2.0 * h * h)
    else:
        raise ValueError(ordering)
    vals = eigh_tridiagonal(diag, off, select="i", select_range=(0, top))[0]
    return float(vals[n])


# =====================================================================
# U1-G1a -- Airy validation: H = -(1/(2 mu)) d^2/dr^2 + k r on (0, inf),
# Dirichlet at 0  =>  E_n = (k^2/(2 mu))^(1/3) |a_n|  (a_n = Ai zeros).
# Textbook parameters mu = 1/2, k = 1  =>  E_n = |a_n| exactly.
# Same code path as production (eigen_vm, constant mass); Richardson
# h^2-extrapolation over two grids.
# =====================================================================
def validate_airy():
    mu_c, k = 0.5, 1.0
    zeros = np.abs(ai_zeros(6)[0])
    Vf = lambda x: k * x
    Mf = lambda x: mu_c * np.ones_like(np.asarray(x, dtype=float))
    xmax = 16.0
    errs = []
    for n in range(6):
        e1 = eigen_vm(Vf, Mf, n=n, ordering="zk", npts=25000, xmax=xmax)
        e2 = eigen_vm(Vf, Mf, n=n, ordering="zk", npts=50000, xmax=xmax)
        # grids h1 = xmax/25001, h2 = xmax/50001: near-exact 4x ratio
        r = (e1 - e2) / ((xmax / 25001.0) ** 2 - (xmax / 50001.0) ** 2)
        e_rich = e2 - r * (xmax / 50001.0) ** 2
        errs.append(abs(e_rich - zeros[n]) / zeros[n])
    err = float(max(errs))
    # constant-mass consistency: bdd == zk for constant mu
    e_zk = eigen_vm(Vf, Mf, n=0, ordering="zk", npts=25000, xmax=xmax)
    e_bd = eigen_vm(Vf, Mf, n=0, ordering="bdd", npts=25000, xmax=xmax)
    RESULTS["validation"]["airy"] = {
        "case": "mu = 1/2, k = 1 (textbook): E_n = |a_n| exactly",
        "airy_zeros": [float(z) for z in zeros],
        "max_rel_err_richardson": err,
        "zk_bdd_constant_mass_agreement": abs(e_zk - e_bd) / e_zk,
        "gate_U1_G1a": "PASS" if err <= 1e-6 else "FAIL"}
    print(f"[U1-G1a] Airy validation: max rel err = {err:.2e} over first 6 "
          f"levels (gate <= 1e-6): "
          f"{RESULTS['validation']['airy']['gate_U1_G1a']}")
    print(f"          zk/bdd constant-mass agreement: "
          f"{abs(e_zk - e_bd)/e_zk:.2e}")
    assert err <= 1e-6


# =====================================================================
# U1-G1b -- T1-limit reproduction: feed the generic solver T1's loop
# potential V = 2 pi x + (Mph - 13/12)(1 - e^{-f* x})/x and inertia
# mu = 2 pi x at T1's grid (npts = 12000, xmax = 14), ZK ordering.
# =====================================================================
def validate_t1_limit(fstar, t1_eps0, t1_eps0mp):
    ZP_T1 = -13.0 / 12.0

    def Vt1(Mph):
        return lambda x: (2.0 * math.pi * x
                          + (Mph + ZP_T1) * (1.0 - np.exp(-fstar * x)) / x)
    Mu = lambda x: 2.0 * math.pi * x
    e0 = eigen_vm(Vt1(0), Mu, n=0, ordering="zk", npts=12000, xmax=14.0)
    e4 = eigen_vm(Vt1(4), Mu, n=0, ordering="zk", npts=12000, xmax=14.0)
    r0 = abs(e0 - t1_eps0) / t1_eps0
    r4 = abs(e4 - t1_eps0mp) / t1_eps0mp
    RESULTS["validation"]["t1_limit"] = {
        "fstar_frozen": fstar,
        "t1_0pp_target": t1_eps0, "reproduced_0pp": e0, "rel_err_0pp": r0,
        "t1_0mp_target": t1_eps0mp, "reproduced_0mp": e4, "rel_err_0mp": r4,
        "gate_U1_G1b": "PASS" if max(r0, r4) <= 1e-3 else "FAIL"}
    print(f"[U1-G1b] T1-limit: 0++ eps = {e0:.10f} vs {t1_eps0:.10f} "
          f"(rel {r0:.2e}); 0-+ eps = {e4:.10f} vs {t1_eps0mp:.10f} "
          f"(rel {r4:.2e})  gate <= 1e-3: "
          f"{RESULTS['validation']['t1_limit']['gate_U1_G1b']}")
    assert max(r0, r4) <= 1e-3


# =====================================================================
# Dumbbell machinery
# =====================================================================
def V_dumbbell(N, mhat, f):
    return lambda x: (3.0 * x + 2.0 * mhat
                      + (math.pi * N + ZPD) * (1.0 - np.exp(-f * x)) / x)


def mu_dumbbell(mhat):
    return lambda x: 3.0 * x + 2.0 * mhat


def mu_stretch(mhat):
    # printed SENSITIVITY variant only (not in the frozen envelope):
    # physically derived symmetric-stretch inertia of the rigid-end
    # dumbbell (ends at +-L/2, linear velocity profile on the tubes):
    # T = 1/2 [m_J/2 + sigma L/4] Ldot^2
    return lambda x: 0.25 * x + 0.5 * mhat


def grid_convergence(fstar):
    Vf, Mf = V_dumbbell(0, 0.0, fstar), mu_dumbbell(0.0)
    Vh, _ = V_dumbbell(2, 0.39, fstar), None
    Mh = mu_dumbbell(0.39)
    e1 = eigen_vm(Vf, Mf, 0, "zk", 16000, 18.0)
    e2 = eigen_vm(Vf, Mf, 0, "zk", 32000, 18.0)
    e3 = eigen_vm(Vf, Mf, 0, "zk", 32000, 24.0)
    g1 = eigen_vm(Vh, Mh, 1, "zk", 16000, 18.0)
    g2 = eigen_vm(Vh, Mh, 1, "zk", 32000, 24.0)
    rel = max(abs(e2 - e1) / e1, abs(e3 - e2) / e2, abs(g2 - g1) / g1)
    RESULTS["self_checks"]["dumbbell_grid_convergence_relerr"] = rel
    print(f"[self-check] dumbbell grid convergence rel err = {rel:.2e}")
    assert rel < 5e-6, "grid not converged"


# =====================================================================
# MAIN
# =====================================================================
def main():
    print("=" * 72)
    print("U1 -- J-Jbar dumbbell spectrum: validation first (U1-G1)")
    print("=" * 72)

    with open(T1_JSON) as fh:
        t1 = json.load(fh)
    fstar = t1["route_A_IP"]["validation"]["fstar"]
    t1_levels = {L["label"]: L for L in t1["route_A_IP"]["levels_at_sigma_0.19"]}
    t1_eps0 = t1_levels["0++ (ground loop)"]["eps"]
    t1_eps0mp = [v for k, v in
                 ((L["label"], L["eps"]) for L in
                  t1["route_A_IP"]["levels_at_sigma_0.19"])
                 if k.startswith("0-+")][0]
    ng_n1 = [e["m_over_sqrtsigma"] for e in
             t1["route_B_NG"]["levels_at_sigma_0.19"] if e["N"] == 1][0]
    print(f"loaded T1 frozen machinery: f* = {fstar:.10f}, "
          f"T1 eps(0++) = {t1_eps0:.6f}, eps(0-+) = {t1_eps0mp:.6f}, "
          f"NG N=1 = {ng_n1:.6f}")

    validate_airy()
    validate_t1_limit(fstar, t1_eps0, t1_eps0mp)
    grid_convergence(fstar)

    amend("U1-A1", "Y-law end geometry idealized: at leading order in the "
          "adiabatic separation coordinate the three tubes run PARALLEL "
          "along the J-Jbar axis (energy exactly 3*sigma*L); the Y-law "
          "120-degree opening (a banked corpus structural consistency, "
          "favored by lattice [IM]) lives in the transverse junction "
          "structure that is collapsed to a point here.  Transverse "
          "junction shape, tube-tube interaction, and junction-core "
          "excitation are all beyond this leading-order dumbbell.")
    amend("U1-A2", "Configuration inertia mu(L) = 3*sigma*L + 2*m_J (the "
          "full configuration mass), imitating T1's mu(rho) = 2*pi*sigma*"
          "rho pattern per the frozen spec.  The physically derived "
          "symmetric-stretch inertia (linear velocity profile, ends at "
          "+-L/2) is sigma*L/4 + m_J/2 and is carried as a printed "
          "SENSITIVITY variant outside the frozen envelope; for the ring "
          "breathing mode the two coincide, for the dumbbell they do not "
          "-- a genuine transfer ambiguity, priced below.")
    amend("U1-A3", "Transverse-phonon content per T1 census conventions "
          "transferred to open tubes: both tube ends fixed on junctions "
          "(Dirichlet), omega_n = n*pi/L, 3 tubes x 2 transverse "
          "polarizations; zeta-regularized zero point = 6 x (-pi/24L) = "
          "-pi/(4L), i.e. three times the single-open-string D=4 Luscher "
          "term.  T1's short-distance factor (1 - e^{-f x}) transferred "
          "as junction-end smoothing at the FROZEN f* = 1.9724; scheme "
          "envelope = T1's frozen set {ZK x f in {f*/1.5, f*, 1.5 f*}} "
          "u {BDD x f in {8, 32}}.  No new parameter introduced and "
          "nothing recalibrated.")
    amend("U1-A4", "m_J enters both the statics (+2 m_J) and the inertia; "
          "evaluated ONLY at the three pre-registered values.  The 0.1355 "
          "anchor is a 2+1D SU(3) lattice determination (JHEP 12 (2025) "
          "019, arXiv:2508.00608 [IM]) -- the 3+1D value is unmeasured; "
          "dimensional transfer is an untested assumption, printed.  "
          "0.39 = f_q(central 0.17 GeV)/sqrt(sigma) is a corpus-internal "
          "scale analogy, not a junction mass.  Junction energy remains "
          "[CJ-new]: the corpus prints NO junction energy/mass/inertia.")

    # ---------------- U1-G2: spectrum at the three mhat values ----------
    print("\n" + "=" * 72)
    print("U1-G2 -- dumbbell spectrum, zero fitted parameters "
          f"(sigma = {SIGMA_GUM} GeV^2 [IM], sqrt(sigma) = {SQS:.5f} GeV)")
    print("=" * 72)
    schemes = [("zk", fstar / 1.5), ("zk", fstar), ("zk", fstar * 1.5),
               ("bdd", 8.0), ("bdd", 32.0)]
    levels = [
        ("ground (scalar-class head, N=0)",                    0, 0),
        ("first phonon level (N=1: 6 channels, 3 tubes x 2 pol)", 1, 0),
        ("separation excitation (N=0, n_r=1)",                 0, 1),
        ("second phonon level (N=2: one n=2 or two n=1)",      2, 0),
        ("phonon + separation (N=1, n_r=1)",                   1, 1),
    ]
    spec = {}
    for mhat in MHAT_SET:
        Mf = mu_dumbbell(mhat)
        out = []
        for label, N, nr in levels:
            e_c = eigen_vm(V_dumbbell(N, mhat, fstar), Mf, n=nr,
                           ordering="zk", npts=16000, xmax=18.0)
            band = [eigen_vm(V_dumbbell(N, mhat, fv), Mf, n=nr,
                             ordering=o, npts=16000, xmax=18.0)
                    for o, fv in schemes]
            lo, hi = min(band), max(band)
            out.append({"label": label, "N": N, "n_r": nr,
                        "eps": e_c, "eps_band": [lo, hi],
                        "mass_GeV_VF": e_c * SQS,
                        "mass_band_GeV_VF": [lo * SQS, hi * SQS]})
        # sensitivity variant (ground only; printed, not in envelope)
        e_s = eigen_vm(V_dumbbell(0, mhat, fstar), mu_stretch(mhat), n=0,
                       ordering="zk", npts=16000, xmax=18.0)
        spec[str(mhat)] = {
            "mhat": mhat, "levels": out,
            "ground_stretch_inertia_sensitivity": {
                "eps": e_s,
                "note": ("mu = sigma*L/4 + m_J/2 variant (AMENDMENT U1-A2); "
                         "lighter inertia raises eps; outside the frozen "
                         "envelope, printed as transfer-ambiguity price")}}
        print(f"\n  m_J/sqrt(sigma) = {mhat}:")
        for L in out:
            print(f"    {L['label']:55s} eps = {L['eps']:8.4f}  "
                  f"band [{L['eps_band'][0]:.4f}, {L['eps_band'][1]:.4f}]  "
                  f"m = {L['mass_GeV_VF']:.3f} GeV [V.F]")
        print(f"    {'(sensitivity: stretch-inertia ground)':55s} "
              f"eps = {e_s:8.4f}")
    RESULTS["spectrum"] = spec
    RESULTS["spectrum_conventions"] = {
        "potential": "V_N(x) = 3x + 2*mhat + (pi*N - pi/4)(1 - e^{-f x})/x",
        "inertia": "mu(x) = 3x + 2*mhat  (spec-frozen, T1 pattern)",
        "band_definition": ("envelope over T1's frozen scheme set "
                            "{ZK x f in [f*/1.5, f*, 1.5 f*]} + "
                            "{BDD x f in [8, 32]}; model systematics only, "
                            "nothing fitted"),
        "jpc_note": ("adiabatic model assignment only: ground = scalar-"
                     "class; the N=1 sextet decomposes under the dumbbell "
                     "symmetry (S3 tube permutation x polarization axial "
                     "doublet J_z = +-1 x end swap) into 1 + 2 per "
                     "polarization -- exact degeneracy at this level is an "
                     "adiabatic artifact (same caveat class as T1 M=2); "
                     "FKRV [IM] give J-Jbar Regge trajectories with "
                     "intercept ~1/2, not used here")}

    # ---------------- U1-G3: confrontation (V.F) -----------------------
    print("\n" + "=" * 72)
    print("U1-G3 -- dimensionless confrontation (V.F locked)")
    print("=" * 72)
    anchors = {"X2370_m_over_sqrtsigma": [5.41, 5.45],
               "lattice_0mp_over_sqrtsigma": [5.5, 5.9],
               "lattice_0pp_over_sqrtsigma": 3.5,
               "T1_loop_IP_0pp": t1_eps0, "T1_loop_IP_0mp": t1_eps0mp,
               "T1_loop_NG_N1": ng_n1}
    conf = {"anchors_IM_and_T1": anchors, "per_mhat": {}}
    t1_ip_all = sorted(L["eps"] for L in
                       t1["route_A_IP"]["levels_at_sigma_0.19"])
    for mhat in MHAT_SET:
        out = spec[str(mhat)]["levels"]
        g = out[0]
        # sector-relation classification vs the T1 loop ladder
        below_loop_ground = g["eps"] < t1_eps0
        dup = abs(g["eps"] - t1_eps0) / t1_eps0 < 0.05
        # which dumbbell levels' envelopes touch the X / lattice-0-+ window
        win = [5.41, 5.9]
        in_win = [L["label"] for L in out
                  if L["eps_band"][0] <= win[1] and L["eps_band"][1] >= win[0]]
        # interleaving: count dumbbell centrals falling between consecutive
        # T1 loop levels
        interleaved = []
        for L in out:
            for a, b in zip(t1_ip_all[:-1], t1_ip_all[1:]):
                if a < L["eps"] < b:
                    interleaved.append((L["label"], round(a, 3), round(b, 3)))
        conf["per_mhat"][str(mhat)] = {
            "ground_eps": g["eps"], "ground_band": g["eps_band"],
            "ground_below_T1_loop_0pp": bool(below_loop_ground),
            "ground_duplicates_T1_loop_0pp_5pct": bool(dup),
            "levels_touching_X_lattice0mp_window_5.41_5.9": in_win,
            "dumbbell_levels_interleaving_T1_loop_gaps": interleaved}
        print(f"  mhat = {mhat}: ground eps = {g['eps']:.4f} "
              f"band [{g['eps_band'][0]:.4f}, {g['eps_band'][1]:.4f}]; "
              f"below T1 loop 0++ ({t1_eps0:.3f})? {below_loop_ground}; "
              f"levels touching 5.41-5.9 window: {len(in_win)}")
    # sector-relation verdict (computed)
    all_below = all(conf["per_mhat"][str(m)]["ground_below_T1_loop_0pp"]
                    for m in MHAT_SET)
    any_dup = any(conf["per_mhat"][str(m)]["ground_duplicates_T1_loop_0pp_5pct"]
                  for m in MHAT_SET)
    n_inter = sum(len(conf["per_mhat"][str(m)]
                      ["dumbbell_levels_interleaving_T1_loop_gaps"])
                  for m in MHAT_SET)
    verdict = []
    if all_below:
        verdict.append("EXTENDS DOWNWARD: at every pre-registered junction "
                       "mass the J-Jbar ground state lies BELOW the loop "
                       "sector's 0++ head -- the junction sector opens a "
                       "new, lighter floor under the closed-loop spectrum")
    if any_dup:
        verdict.append("near-duplication of the loop 0++ head at the "
                       "heaviest mhat")
    verdict.append(f"INTERLEAVES: {n_inter} dumbbell level placements fall "
                   "inside gaps of the T1 loop ladder across the mhat set "
                   "-- the two sectors would overlay as one interleaved "
                   "spectrum, not two separated bands")
    # inertia-convention caveat, computed (AMENDMENT U1-A2 price)
    s_eps = {m: spec[str(m)]["ground_stretch_inertia_sensitivity"]["eps"]
             for m in MHAT_SET}
    if any(v > t1_eps0 for v in s_eps.values()):
        verdict.append(
            "CONVENTION-FRAGILE: under the physically derived stretch "
            "inertia (U1-A2 sensitivity variant) the ground rises to eps "
            "= " + ", ".join(f"{s_eps[m]:.2f} (mhat={m})" for m in MHAT_SET)
            + " -- ABOVE the loop 0++ head; the downward-extension half of "
            "this verdict is owned by the spec-frozen T1-pattern inertia, "
            "not by the dumbbell physics itself, and the interleaving half "
            "survives both conventions")
    conf["sector_relation_verdict"] = "; ".join(verdict)
    conf["tetraquark_kill_note"] = (
        "n = 3 ONLY: every network here is built from order-3 junctions "
        "(the corpus's Q-1 inventory class with the tetraquark kill armed, "
        "junction order owned by WS1-H3 -- an ansatz, [CAL]-at-best, seam "
        "F-K0-1 carried).  No 4-fold junction and no baryonium-TETRAQUARK "
        "(quark-terminated, i.e. knot-terminated) state is constructed or "
        "priced; the J-Jbar dumbbell is the knot-FREE minimal network.  "
        "Nothing here bears on the tetraquark kill's standing.")
    conf["vf_statement"] = (
        "All GeV values are sqrt(sigma)-rescalings of dimensionless model "
        "output with sigma = 0.19 GeV^2 [IM]: dimensionally secure / "
        "structurally plausible / quantitatively unclaimed.  No junction-"
        "sector state is claimed as any observed hadron.")
    print("\n  sector-relation verdict:", conf["sector_relation_verdict"])
    RESULTS["confrontation"] = conf

    # ---------------- U1-G4 register ------------------------------------
    RESULTS["register_G4"] = {
        "seal_q_theta": ("V.F on every hadron-mass-confronting number; "
                         "pipeline frozen by ROADMAP_v12 before execution"),
        "numerology_preemption": (
            "The corpus's dimensionless closure benchmark c-frak = 2.37 "
            "+- 0.09 shares digits with GeV-denominated numbers (e.g. "
            "m_X ~ 2.37 GeV) by man-made-unit coincidence only; no digit "
            "identity is cited anywhere in U1 as structure, and may never "
            "be."),
        "junction_energy_status": ("[CJ-new]: the corpus prints no junction "
                                   "energy/mass/inertia; m_J here is an "
                                   "explicit scanned parameter, never "
                                   "fitted, anchored only by the [IM] "
                                   "2+1D lattice value and a corpus scale "
                                   "analogy"),
        "no_star_fitting": ("STAR observables (1.84, 0.64, 1.04) appear "
                            "nowhere in this workstream's computations; "
                            "mechanical self-scan below"),
        "no_grade_motion": True, "no_stake": True}
    # mechanical self-scan: STAR anchor values must not appear as numeric
    # tokens anywhere outside the disclaimer sentences.  Needles built by
    # concatenation so the scan does not trip on itself; disclaimer triple
    # "(v1, v2, v3)" removed before scanning; token-boundary regex so
    # computed values like 2.064 cannot false-positive.
    import re
    needles = ["1." + "84", "0." + "64", "1." + "04"]
    triple = ", ".join(needles)
    src = open(os.path.abspath(__file__)).read().replace(triple, "")
    jtxt = json.dumps(RESULTS).replace(triple, "")
    hits = []
    for s in needles:
        pat = re.compile(r"(?<![\d.])" + re.escape(s) + r"(?![\d])")
        for name, txt in (("source", src), ("json", jtxt)):
            if pat.search(txt):
                hits.append(f"{s} in {name}")
    RESULTS["register_G4"]["star_value_scan_hits_outside_disclaimer"] = hits
    print(f"\n[U1-G4] STAR-value self-scan (outside disclaimer text): "
          f"{hits if hits else 'clean'}")

    with open(os.path.join(HERE, "u1_results.json"), "w") as fh:
        json.dump(RESULTS, fh, indent=1)
    print("wrote u1_results.json")
    make_figure(spec, anchors, t1, t1_eps0, t1_eps0mp, ng_n1)
    return RESULTS


# ---------------------------------------------------------------- figure
def make_figure(spec, anchors, t1, t1_eps0, t1_eps0mp, ng_n1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
    C_DB, C_T1 = "#2a78d6", "#eb6834"   # categorical slots 1, 2 (validated)
    fig, ax = plt.subplots(figsize=(9.4, 6.6), dpi=170)
    fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)

    # [IM] anchor bands, neutral grays, full width
    ax.axhspan(*anchors["lattice_0mp_over_sqrtsigma"], color="#e2e0da",
               zorder=0)
    ax.axhspan(*anchors["X2370_m_over_sqrtsigma"], color="#c9c6bd", zorder=1)
    ax.axhline(anchors["lattice_0pp_over_sqrtsigma"], color="#a5a29a",
               lw=1.4, ls=(0, (5, 3)), zorder=1)
    ax.text(4.68, sum(anchors["lattice_0mp_over_sqrtsigma"]) / 2 + 0.30,
            "lattice 0$^{-+}$ 5.5–5.9 [IM]", ha="right", va="bottom",
            fontsize=8.2, color=INK2)
    ax.text(4.68, anchors["X2370_m_over_sqrtsigma"][0] - 0.08,
            "X(2370) 5.41–5.45 [IM]", ha="right", va="top",
            fontsize=8.2, color=INK2)
    ax.text(4.68, anchors["lattice_0pp_over_sqrtsigma"] + 0.05,
            "lattice 0$^{++}$ ≈ 3.5 [IM]", ha="right", va="bottom",
            fontsize=8.2, color=INK2)

    def draw(xc, entries, color, w=0.30):
        for k, (lab, eps, band, dy) in enumerate(entries):
            ax.plot([xc - w, xc + w], [eps, eps], color=color, lw=2.2,
                    zorder=3, solid_capstyle="butt")
            if band[1] > band[0]:
                xw = xc - w - 0.05 - 0.05 * (k % 3)
                ax.plot([xw, xw], band, color=color, lw=1.3, zorder=2)
                for b in band:
                    ax.plot([xw - 0.02, xw + 0.02], [b, b], color=color,
                            lw=1.3, zorder=2)
            if lab:
                ax.annotate(lab, (xc + w + 0.02, eps + dy), fontsize=7.9,
                            color=INK, va="center", ha="left", zorder=4)

    short = {0: "N=0 ground", 1: "N=1 phonon ×6", 2: "n$_r$=1",
             3: "N=2", 4: "N=1, n$_r$=1"}
    xs = [0.62, 1.78, 2.94]
    for xc, mhat in zip(xs, MHAT_SET):
        entries = []
        lv = spec[str(mhat)]["levels"]
        order = sorted(range(len(lv)), key=lambda i: lv[i]["eps"])
        for rank, i in enumerate(order):
            L = lv[i]
            lab = short[i] if xc == xs[0] else ""
            dy = 0.0
            # de-collide close labels in first column
            if xc == xs[0] and rank > 0:
                prev = lv[order[rank - 1]]["eps"]
                if L["eps"] - prev < 0.30:
                    dy = 0.18
            entries.append((lab, L["eps"], L["eps_band"], dy))
        draw(xc, entries, C_DB)

    t1entries = [("loop 0$^{++}$ 3.58", t1_eps0, [0, 0], -0.34),
                 ("NG N=1 4.80", ng_n1, [0, 0], 0),
                 ("loop 0$^{-+}$ 8.27", t1_eps0mp, [0, 0], 0)]
    draw(4.10, t1entries, C_T1)

    ax.set_xlim(0.1, 4.72); ax.set_ylim(0, 11.0)
    ax.set_xticks(xs + [4.10])
    ax.set_xticklabels(
        ["J–J̄ dumbbell\n$m_J/\\sqrt{\\sigma}$ = 0",
         "J–J̄ dumbbell\n$m_J/\\sqrt{\\sigma}$ = 0.1355\n(2+1D lattice [IM])",
         "J–J̄ dumbbell\n$m_J/\\sqrt{\\sigma}$ = 0.39\n($f_q$/$\\sqrt{\\sigma}$)",
         "Tier-9 T1 loop\n(reference)"],
        fontsize=8.6, color=INK)
    ax.set_ylabel(r"$m/\sqrt{\sigma}$   (dimensionless; at $\sigma$ = 0.19 "
                  r"GeV$^2$ [IM]: $m$ = value × 0.436 GeV)",
                  fontsize=9, color=INK)
    ax.tick_params(colors=INK2, labelsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#a5a29a")
    ax.grid(axis="y", color="#eceae4", lw=0.8, zorder=-2)
    ax.set_axisbelow(True)
    ax.set_title("U1 — junction–antijunction dumbbell vs the Tier-9 loop "
                 "sector and [IM] anchors\nwithin-model, V.F grade: "
                 "dimensionally secure / structurally plausible / "
                 "quantitatively unclaimed", fontsize=10, color=INK, pad=12)
    fig.text(0.01, 0.012,
             "whiskers: frozen T1 scheme envelope (ZK×f ∪ BDD×f) — model "
             "systematics only; zero fitted parameters; m$_J$ scanned at "
             "three pre-registered values, never fitted.\nJunction energy "
             "is [CJ-new] (unprinted); inertia convention dominates beyond "
             "the drawn envelope (stretch variant ground ≈ 5.5–5.8; "
             "RESULTS §5).",
             fontsize=7.6, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = os.path.join(HERE, "u1_fig.png")
    fig.savefig(out, facecolor=SURF)
    print("wrote", out)


if __name__ == "__main__":
    main()
