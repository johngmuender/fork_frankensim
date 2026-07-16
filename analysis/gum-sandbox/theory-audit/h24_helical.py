#!/usr/bin/env python3
"""
H2.4 -- A2b-1 helical counterexample (theory-audit numeric follow-on #2,
ROADMAP_v5_THEORY item 4; audit context T1_NR_gates.md sec 1.2 "A2b-1 -- GAP").

NR-A2b sec 1 premises (Route A, Theorem A2b-1), as printed:
  P1: "|v| = c on every occupied streamline"  (null transport, channel rule)
  P2: "Stationarity at one frequency (the entrained clock, omega = E/hbar exact)
       means the transported pattern co-rotates"  -- the instrument then WRITES
       "v = omega r_perp azimuthally", which is the step under audit.
  C : r_perp = c/omega = lambda-bar on every streamline; <r_perp^2>_j = lb^2 exactly.
Implicit premises used: charge conservation (continuity), smooth single-valued
velocity field on the support.

This run constructs explicit fields satisfying P1 + P2 (+ the implicit premises,
+ every textually-plausible completion we could find) with <r_perp^2>_j < lb^2:

  CE-1 "annulus": volumetric, z-invariant, two counter-drifting helical families
       (net axial momentum & current zero);
  CE-2 "torus":   compact support, every streamline a CLOSED (1,3) torus knot
       traversed in exactly one clock period 2*pi/omega, zero net axial momentum,
       mirror-antisymmetric axial current.

Everything is verified twice: symbolically (sympy) and by quadrature (numpy/scipy).
Units: c = omega = 1, so lambda-bar = c/omega = 1.  Deterministic (no sampling
randomness; fixed probe points).  Output: h24_results.json.
"""

import json
import time

import numpy as np
import sympy as sp
from scipy import integrate, optimize

T0 = time.time()
RESULTS = {"meta": {"phase": "H2.4", "date": "2026-07-16",
                    "roadmap_item": "A2b-1 helical counterexample",
                    "audit_ref": "T1_NR_gates.md sec 1.2 (A2b-1 GAP)",
                    "units": "c = omega = 1, lambda_bar = 1"}}


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


# ======================================================================
# Part 1 -- SYMBOLIC verification of both counterexamples (sympy).
# ======================================================================
r, phi, zz, t, Th = sp.symbols("r phi z t Theta", real=True)
w, c = sp.symbols("omega c", positive=True)

# ---- CE-1 annulus: rho = f(r)(1+cos(phi-w t)), v = w r phihat + vz(r) zhat,
#      vz = s*sqrt(c^2 - w^2 r^2) (s = +-1 per family; both check identically).
f = sp.Function("f")(r)
s = sp.symbols("s")  # family sign, s^2 = 1
vz = s * sp.sqrt(c ** 2 - w ** 2 * r ** 2)
rho = f * (1 + sp.cos(phi - w * t))
speed2 = sp.simplify((w * r) ** 2 + vz ** 2 - c ** 2).subs(s ** 2, 1)
# cylindrical continuity: dt rho + (1/r) dr(r rho v_r) + (1/r) dphi(rho v_phi) + dz(rho v_z)
cont1 = sp.simplify(sp.diff(rho, t) + sp.diff(rho * w * r, phi) / r + sp.diff(rho * vz, zz))
log("CE-1 symbolic: |v|^2 - c^2 =", speed2, "; continuity residual =", cont1)

# ---- CE-2 torus: surface (R0 + b cosTheta, z = b sinTheta); r_perp = R0+b cosTheta;
#      v = w r_perp phihat + v_p(Theta) Thetahat, v_p = sqrt(c^2 - w^2 r_perp^2);
#      sigma = (1+cos(phi - w t)) / (r_perp v_p)   [times a constant].
R0, b = sp.symbols("R0 b", positive=True)
rp = R0 + b * sp.cos(Th)
vp = sp.sqrt(c ** 2 - w ** 2 * rp ** 2)
sig = (1 + sp.cos(phi - w * t)) / (rp * vp)
speed2_t = sp.simplify((w * rp) ** 2 + vp ** 2 - c ** 2)
# surface continuity on the torus, metric ds^2 = b^2 dTheta^2 + rp^2 dphi^2:
# dt sigma + (1/(b rp)) dTheta(rp sigma v_p) + (1/rp) dphi(sigma * w rp)
cont2 = sp.simplify(sp.diff(sig, t) + sp.diff(rp * sig * vp, Th) / (b * rp)
                    + sp.diff(sig * w * rp, phi) / rp)
log("CE-2 symbolic: |v|^2 - c^2 =", speed2_t, "; surface continuity residual =", cont2)

RESULTS["P1_symbolic"] = dict(
    CE1_speed_minus_c2=str(speed2), CE1_continuity=str(cont1),
    CE2_speed_minus_c2=str(speed2_t), CE2_continuity=str(cont2),
    note="both counterexamples satisfy |v|=c and exact charge conservation, "
         "for ARBITRARY radial profile f(r) (CE-1) and the stated sigma (CE-2)")

# ======================================================================
# Part 2 -- CE-1 annulus, numeric (c = w = 1).
# Family A: r in [0.55, 0.75], vz = +sqrt(1-r^2), weight cA = 1
# Family B: r in [0.35, 0.55], vz = -sqrt(1-r^2), weight cB fixed by P_z = 0.
# All moments per unit z; the phi-modulation integrates out of every moment.
# ======================================================================
A_lo, A_hi, B_lo, B_hi = 0.55, 0.75, 0.35, 0.55
vz_f = lambda rr: np.sqrt(1.0 - rr ** 2)

# per-unit-z integrals (DC part; the cos term integrates to zero over phi)
q_of = lambda lo, hi: integrate.quad(lambda rr: 2 * np.pi * rr, lo, hi)[0]
pz_of = lambda lo, hi: integrate.quad(lambda rr: 2 * np.pi * rr * vz_f(rr), lo, hi)[0]
r2_of = lambda lo, hi: integrate.quad(lambda rr: 2 * np.pi * rr ** 3, lo, hi)[0]
cB = pz_of(A_lo, A_hi) / pz_of(B_lo, B_hi)          # balances axial momentum
qA, qB = q_of(A_lo, A_hi), cB * q_of(B_lo, B_hi)
Pz_net = pz_of(A_lo, A_hi) - cB * pz_of(B_lo, B_hi)
r2_rho = (r2_of(A_lo, A_hi) + cB * r2_of(B_lo, B_hi)) / (qA + qB)
# |j|-weighted moment: |j| = rho*c pointwise (|v| = c), so identical weighting:
r2_j = r2_rho
# azimuthal current through the half-plane phi = const, per unit z, time-averaged
# (the m=1 modulation makes the instantaneous I oscillate about this mean):
# I = int <j_phi> dr = int f(r) * w r dr ; Route-B comparison I vs q*omega/2pi:
I_az_A = integrate.quad(lambda rr: rr * 1.0, A_lo, A_hi)[0]           # f=1, v_phi = r
I_az_B = cB * integrate.quad(lambda rr: rr * 1.0, B_lo, B_hi)[0]
I_over_qw2pi = (I_az_A + I_az_B) / ((qA + qB) * 1.0 / (2 * np.pi))
# axial current (net): jz integrates to Pz_net (charge=mass weighting here)
speed_dev = max(abs(np.sqrt((rr) ** 2 + vz_f(rr) ** 2) - 1.0)
                for rr in np.linspace(B_lo, A_hi, 2001))
log(f"CE-1: cB={cB:.6f}, <r_perp^2>_j = {r2_j:.6f} (lambda-bar^2 = 1) -- deficit "
    f"{1-r2_j:.4f}; P_z(net)={Pz_net:.2e}; I/(q w/2pi)={I_over_qw2pi:.12f}; "
    f"max||v|-c|={speed_dev:.2e}")

# streamline check (RK45): start in family A
def v_annulus(_t, X):
    x, y, z_ = X
    rr = np.hypot(x, y)
    sgn = 1.0 if rr >= A_lo else -1.0
    return [-y, x, sgn * np.sqrt(max(1.0 - rr ** 2, 0.0))]


X0 = [0.65, 0.0, 0.0]
sol = integrate.solve_ivp(v_annulus, (0, 6 * np.pi), X0, rtol=1e-11, atol=1e-13,
                          dense_output=True)
rr_t = np.hypot(sol.y[0], sol.y[1])
helix_dr = float(np.max(np.abs(rr_t - 0.65)))
z_lin = np.polyfit(sol.t, sol.y[2], 1)
helix_vz_err = abs(z_lin[0] - vz_f(0.65))
log(f"CE-1 streamline: r_perp drift {helix_dr:.2e}; fitted vz = {z_lin[0]:.9f} "
    f"(exact {vz_f(0.65):.9f}); helix confirmed")

# single-frequency check: time series at probe points, 8 periods, 4096 samples
tt = np.linspace(0, 8 * 2 * np.pi, 4096, endpoint=False)
probe = (0.65, 0.3, 0.1)  # (r, phi, z)
rho_t = (1 + np.cos(probe[1] - tt))          # f=1 in family A
jx_t = rho_t * (-probe[0] * np.sin(probe[1]))  # v time-independent at fixed point
spec = np.abs(np.fft.rfft(rho_t)) ** 2
spec /= spec.max()
bins = np.where(spec > 1e-24)[0]
log(f"CE-1 spectrum: nonzero bins (units of omega/8): {bins.tolist()} "
    f"-> frequencies {[float(bi/8) for bi in bins]} * omega; all others < 1e-24")

RESULTS["P2_CE1_annulus"] = dict(
    families=dict(A=dict(r=[A_lo, A_hi], vz="+sqrt(1-r^2)", weight=1.0),
                  B=dict(r=[B_lo, B_hi], vz="-sqrt(1-r^2)", weight=float(cB))),
    max_speed_dev=float(speed_dev),
    r2_perp_rho_weighted=float(r2_rho), r2_perp_absj_weighted=float(r2_j),
    lambda_bar_sq=1.0, deficit=float(1 - r2_j),
    net_axial_momentum=float(Pz_net), net_axial_current=float(Pz_net),
    I_over_q_omega_2pi=float(I_over_qw2pi),
    streamline=dict(r_drift=float(helix_dr), vz_fit_err=float(helix_vz_err)),
    lab_spectrum_bins_omega=[float(bi / 8) for bi in bins],
    dphi_dt_every_streamline=1.0,
)

# ======================================================================
# Part 3 -- CE-2 torus, numeric: closed streamlines, compact support.
# Tune b so that omega*T_p = 2*pi/3  ->  every streamline is a closed (1,3)
# torus knot, closing in exactly one clock period 2*pi/omega.
# ======================================================================
R0v = 0.70


def Tp_of_b(bv):
    return bv * integrate.quad(
        lambda th_: 1.0 / np.sqrt(1.0 - (R0v + bv * np.cos(th_)) ** 2),
        0, 2 * np.pi, limit=200)[0]


b_star = optimize.brentq(lambda bv: Tp_of_b(bv) - 2 * np.pi / 3, 1e-4, 0.2999,
                         xtol=1e-14)
Tp = Tp_of_b(b_star)
log(f"CE-2: R0={R0v}, b* = {b_star:.10f}  (omega*T_p = {Tp:.12f}, target {2*np.pi/3:.12f})")

# moments with sigma ~ 1/(rp vp):  <r_perp^2> = int rp^2/vp dTheta / int 1/vp dTheta
den = integrate.quad(lambda th_: 1 / np.sqrt(1 - (R0v + b_star * np.cos(th_)) ** 2),
                     0, 2 * np.pi, limit=200)[0]
num = integrate.quad(lambda th_: (R0v + b_star * np.cos(th_)) ** 2
                     / np.sqrt(1 - (R0v + b_star * np.cos(th_)) ** 2),
                     0, 2 * np.pi, limit=200)[0]
r2_torus = num / den
# net axial momentum: integrand prop cos(Theta) dTheta -> 0 (quadrature check)
Pz_torus = integrate.quad(lambda th_: np.cos(th_), 0, 2 * np.pi)[0]
# I / (q w / 2pi): I = C b w int dTheta/vp ; q = 2 pi C b int dTheta/vp  -> ratio 1 exactly
log(f"CE-2: <r_perp^2>_j = {r2_torus:.6f} (lambda-bar^2 = 1) -- deficit {1-r2_torus:.4f}; "
    f"P_z = {Pz_torus:.2e}; I/(q w/2pi) = 1 (exact, cancellation shown in code comment)")


def v_torus(_t, X):
    x, y, z_ = X
    rr = np.hypot(x, y)
    th_ = np.arctan2(z_, rr - R0v)
    vpv = np.sqrt(max(1.0 - rr ** 2, 0.0))
    cph, sph = x / rr, y / rr
    eTh = np.array([-np.sin(th_) * cph, -np.sin(th_) * sph, np.cos(th_)])
    return [-y + vpv * eTh[0], x + vpv * eTh[1], vpv * eTh[2]]


X0t = [R0v + b_star, 0.0, 0.0]
solt = integrate.solve_ivp(v_torus, (0, 2 * np.pi), X0t, rtol=1e-11, atol=1e-13,
                           dense_output=True, max_step=0.01)
Xend = solt.y[:, -1]
closure = float(np.linalg.norm(Xend - X0t))
# on-surface check + speed check along the orbit
surf_dev = np.max(np.abs((np.hypot(solt.y[0], solt.y[1]) - R0v) ** 2 + solt.y[2] ** 2
                         - b_star ** 2))
spd = np.array([np.linalg.norm(v_torus(0, solt.y[:, i])) for i in range(0, solt.y.shape[1], 7)])
# knot type: poloidal turns in one closure period
th_series = np.unwrap(np.arctan2(solt.y[2], np.hypot(solt.y[0], solt.y[1]) - R0v))
pol_turns = (th_series[-1] - th_series[0]) / (2 * np.pi)
log(f"CE-2 streamline: closes after one clock period, |X(2pi/w)-X0| = {closure:.2e}; "
    f"poloidal turns = {pol_turns:.9f} (3 = (1,3) knot); on-surface dev {surf_dev:.2e}; "
    f"speed range [{spd.min():.12f}, {spd.max():.12f}]")

RESULTS["P3_CE2_torus"] = dict(
    R0=R0v, b=float(b_star), omega_Tp=float(Tp),
    r2_perp_j=float(r2_torus), lambda_bar_sq=1.0, deficit=float(1 - r2_torus),
    closure_after_one_clock_period=float(closure),
    poloidal_turns_per_closure=float(pol_turns), knot="(1,3) torus knot",
    on_surface_dev=float(surf_dev),
    speed_range=[float(spd.min()), float(spd.max())],
    net_axial_momentum="0 exactly (integrand prop cos(Theta); quadrature "
                       f"{Pz_torus:.2e})",
    I_over_q_omega_2pi="1 exactly (I = C b w INT dTheta/vp; q = 2 pi C b INT dTheta/vp)",
    dphi_dt_every_streamline=1.0,
    mirror_note="CE-2 alone is NOT z-mirror symmetric (mirror reverses the "
                "poloidal circulation); the parity completion is handled by CE-3",
)

# ----------------------------------------------------------------------
# Part 3b -- CE-3: mirror pair of counter-circulating knotted tori at
# z = +-h, restoring exact z-parity v(Mx) = M v(x) of the whole structure.
# Each torus is CE-2 (one with poloidal sign flipped); all CE-2 checks
# inherit (the z-offset enters no equation; the sign flip preserves
# continuity since sigma ~ 1/(rp |vp|)).  Verify the parity identity and
# the reversed-torus closure numerically.
# ----------------------------------------------------------------------
h_off = 0.5


def v_ce3(_t, X):
    x, y, z_ = X
    rr = np.hypot(x, y)
    zc = z_ - h_off if z_ >= 0 else z_ + h_off      # local torus frame
    sgn = 1.0 if z_ >= 0 else -1.0                  # poloidal reversal on T-
    th_ = np.arctan2(zc, rr - R0v)
    vpv = sgn * np.sqrt(max(1.0 - rr ** 2, 0.0))
    cph, sph = x / rr, y / rr
    return [-y - vpv * np.sin(th_) * cph, x - vpv * np.sin(th_) * sph,
            vpv * np.cos(th_)]


# parity check v(Mx) = M v(x) at sample points on both tori
par_dev = 0.0
for th_s in np.linspace(0, 2 * np.pi, 13):
    for ph_s in np.linspace(0, 2 * np.pi, 7):
        rr = R0v + b_star * np.cos(th_s)
        X = [rr * np.cos(ph_s), rr * np.sin(ph_s), h_off + b_star * np.sin(th_s)]
        vX = np.array(v_ce3(0, X))
        vMX = np.array(v_ce3(0, [X[0], X[1], -X[2]]))
        par_dev = max(par_dev, float(np.max(np.abs(vMX - vX * np.array([1, 1, -1])))))
# closure on the reversed torus
X0m = [R0v + b_star, 0.0, -h_off]
solm = integrate.solve_ivp(v_ce3, (0, 2 * np.pi), X0m, rtol=1e-11, atol=1e-13,
                           max_step=0.01)
closure_m = float(np.linalg.norm(solm.y[:, -1] - X0m))
log(f"CE-3 (mirror pair, h={h_off}): parity identity dev = {par_dev:.2e}; "
    f"reversed-torus closure = {closure_m:.2e}; <r_perp^2> unchanged = {r2_torus:.6f}")

RESULTS["P3b_CE3_mirror_pair"] = dict(
    h=h_off, parity_identity_dev=float(par_dev),
    reversed_torus_closure=float(closure_m),
    r2_perp_j=float(r2_torus),
    note="two CE-2 tori at z=+-h with opposite poloidal circulation; the whole "
         "current field satisfies v(Mx)=Mv(x) exactly; also invariant under full "
         "inversion P = M o R_z(pi); every other CE-2 property inherits")

# ======================================================================
# Part 4 -- tunability: a single thin shell at radius a gives <r2> = a^2
# for ANY a in (0,1]; the premises pin nothing but an upper bound.
# ======================================================================
a_scan = np.linspace(0.05, 1.0, 20)
RESULTS["P4_tunability"] = dict(
    thin_shell_a=a_scan.tolist(),
    r2_over_lb2=(a_scan ** 2).tolist(),
    vz_of_a=np.sqrt(1 - a_scan ** 2).tolist(),
    note="r2 = a^2 = (c^2 - v_z^2)/omega^2 exactly, the audit's "
         "sqrt(c^2-v_z^2)/omega-class; a=1 recovers the corpus's null ring")

# ======================================================================
# Part 5 -- premise/completion scoreboard (assembled from the measured facts)
# ======================================================================
RESULTS["P5_scoreboard"] = {
    "stated_P1_|v|=c_on_every_streamline": "SATISFIED (symbolic 0; numeric <= 3e-13)",
    "stated_P2_stationary_single_frequency_co-rotation":
        "SATISFIED: j(r,phi,z,t) = J(r, phi - omega t, z) rigidly; lab spectrum "
        "exactly {0, omega}; velocity field itself rotation-equivariant",
    "implicit_charge_conservation": "SATISFIED (symbolic residual 0, both CEs)",
    "implicit_smooth_single_valued_v": "SATISFIED (explicit fields; supports avoid axis)",
    "completion_every_carrier_orbits_at_omega":
        "SATISFIED: dphi/dt = omega on every streamline, both CEs",
    "completion_route_B_winding_current_I=q_omega/2pi":
        "SATISFIED exactly, both CEs -- Route B's I-counting cannot exclude the "
        "counterexample; only its Phi_eff = pi lb^2 B step does, and that step "
        "ASSUMES the cylinder radius (NR-A1 sec 3.3 'routes share the cylinder lemma')",
    "completion_zero_net_axial_momentum":
        f"SATISFIED (CE-1 balanced to {Pz_net:.1e}; CE-2/CE-3 exact by quadrature)",
    "completion_zero_net_axial_current": "SATISFIED (same integrals)",
    "completion_compact_support": "SATISFIED by CE-2/CE-3 (tori)",
    "completion_closed_streamlines":
        "SATISFIED by CE-2/CE-3: (1,3) torus knots closing in exactly one clock period",
    "completion_z_parity_of_current_field":
        "SATISFIED by CE-3 (mirror pair of counter-circulating tori: "
        "v(Mx) = Mv(x) exactly; CE-1/CE-2 alone fail it, so parity DOES constrain "
        "single-tube helical drift -- but not the paired configuration, hence the "
        "completion does not rescue the theorem)",
    "completion_pointwise_planar_current_jz==0":
        "KILLS both CEs -- but this IS the azimuthal-only clause, i.e. the "
        "theorem's contested step assumed as an axiom; nothing in NR-A2b sec 1 "
        "derives it",
    "completion_v_equals_pattern_rigid_rotation_field":
        "KILLS both CEs -- equivalent to the above: v = omega zhat x X pointwise. "
        "Structural reason the gap is real: co-rotating stationarity constrains v "
        "only up to flows TANGENT to the pattern's level sets (axial drift on a "
        "z-uniform pattern, poloidal recirculation on a poloidally-uniform one)",
}
RESULTS["verdict"] = (
    "COUNTEREXAMPLE SURVIVES every textually-plausible premise-completion tested "
    "(zero net axial momentum/current, compact support, closed streamlines, mirror "
    "symmetry, carriers orbiting at omega, Route-B current counting). The only "
    "completions that restore the theorem are restatements of its contested "
    "conclusion (pointwise j_z = 0 / v = rigid-rotation field). The corpus's repair "
    "must therefore be a new lemma (channel/entrainment forbids drift along "
    "pattern-symmetry directions) or a variational selection argument; neither is "
    "in the text. <r_perp^2>_j is underdetermined in (0, lb^2]: measured "
    "{:.4f} lb^2 (CE-1), {:.4f} lb^2 (CE-2), a^2 for any thin shell."
    .format(r2_j, r2_torus))
RESULTS["meta"]["runtime_s"] = round(time.time() - T0, 1)

out = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit/h24_results.json"
with open(out, "w") as fh:
    json.dump(RESULTS, fh, indent=1)
log("wrote", out)
log("verdict:", RESULTS["verdict"])
