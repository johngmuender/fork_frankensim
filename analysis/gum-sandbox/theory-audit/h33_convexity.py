#!/usr/bin/env python3
"""
H3.3 Part A (T-H8): Lemma II.1 / App. F.1 pairs-only convexity claim.

Corpus statements audited (01-GUM-Omega-Paper-v2.0.1.md):

  Lemma II.1 (main text):  "For all admissible u, deg R~[u] = 0 (Diff_c(R^3)
  contractible; degree homotopy-invariant), hence deg P~ = deg Q~ == K: ...
  can change only where det F -> 0"

  App. F.1 (the printed proof): "Admissible deformations (det F > 0, identity
  at infinity) form a convex, hence contractible, set; R~[u] is a continuous
  map into Maps(R^3, SU(2)) from a contractible domain => deg R~ == deg R~[0]
  = 0; deg P~ = deg Q~. Contrapositive: Delta-Sigma-K != 0 => det F -> 0."

  Theorem VIII'.1 uses exactly the contrapositive ("Lemma II.1 contrapositive:
  single-knot creation is forbidden in tear-free matter").

This script machine-checks:
  C1  det(grad phi_twist) == 1 symbolically for a GENERIC twist profile
      theta(r)  (so the twist map is admissible: det F = 1 > 0, id at inf).
  C2  the straight segment (1-t)*id + t*phi_twist has det F = (1-2t)^2 at the
      core when theta = pi there: det -> 0 at t = 1/2  => the admissible set
      is NOT convex (and not star-shaped about u = 0).
  C3  strict counterexample (open-condition robust): phi = twist o stretch is
      admissible (det > 0 verified symbolically by chain rule + numerically
      on a grid), F(0) = diag(-2,-1,1), and the segment to the identity has
      det F_t(0) = (1-3t)(1-2t) < 0 on the OPEN interval t in (1/3, 1/2).
  C4  the degree integrator is validated on a hedgehog control map
      (|deg| = 1), and the twist texture's SU(2) lift has degree 0 along the
      whole ADMISSIBLE path s |-> theta_s = s*theta (integrand identically 0)
      -- the history-path repair operating exactly where convexity fails.

Deterministic; python3 + sympy + numpy. Writes its section of h33_summary.json.
"""

import json
import pathlib

import numpy as np
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
SUMMARY = HERE / "h33_summary.json"

checks = []            # (name, passed: bool, detail)


def check(name, passed, detail=""):
    checks.append({"name": name, "pass": bool(passed), "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}  {detail}")


# ----------------------------------------------------------------------
# C1. The twist map phi(x) = R_z(theta(r)) x is admissible: det F == 1.
#     (Each sphere r = const rotates rigidly; generic profile theta(r).)
# ----------------------------------------------------------------------
x, y, z, t = sp.symbols("x y z t", real=True)
r = sp.sqrt(x**2 + y**2 + z**2)
theta = sp.Function("theta")

c, s = sp.cos(theta(r)), sp.sin(theta(r))
Rz = sp.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
phi_twist = Rz * sp.Matrix([x, y, z])
F_twist = phi_twist.jacobian([x, y, z])
det_twist = sp.simplify(F_twist.det())
check("C1 det(grad phi_twist) == 1 for generic theta(r) [sympy]",
      sp.simplify(det_twist - 1) == 0, f"det = {det_twist}")

# identity at infinity: theta has compact support by construction (theta = pi
# for r <= 1, smoothstep down to 0 for r >= 2) -- profile used numerically:


def theta_prof(rr, amp=np.pi):
    """C^2 profile: amp for r<=1, quintic smoothstep to 0 on [1,2], 0 beyond."""
    u = np.clip(2.0 - rr, 0.0, 1.0)
    sstep = 6 * u**5 - 15 * u**4 + 10 * u**3
    return amp * sstep


check("C1b theta profile: theta(0.5)=pi, theta(2.5)=0",
      abs(theta_prof(0.5) - np.pi) < 1e-15 and theta_prof(2.5) == 0.0)

# ----------------------------------------------------------------------
# C2. Segment identity -> twist: at a core point (theta == pi exactly),
#     F_t = (1-t) I + t R_z(pi);  det = (1-2t)^2 -> 0 at t = 1/2.
#     The admissible set {det F > 0} is therefore NOT convex, and NOT even
#     star-shaped about the reference u = 0.
# ----------------------------------------------------------------------
th = sp.symbols("th", real=True)
Rz_th = sp.Matrix([[sp.cos(th), -sp.sin(th), 0],
                   [sp.sin(th), sp.cos(th), 0],
                   [0, 0, 1]])
seg = (1 - t) * sp.eye(3) + t * Rz_th
det_seg = sp.simplify(seg.det())
det_seg_pi = sp.simplify(det_seg.subs(th, sp.pi))
check("C2 det((1-t)I + t R_z(pi)) == (1-2t)^2 [sympy]",
      sp.simplify(det_seg_pi - (1 - 2 * t) ** 2) == 0, f"= {det_seg_pi}")
check("C2b segment det at t=1/2 is 0 (admissibility det>0 fails)",
      det_seg_pi.subs(t, sp.Rational(1, 2)) == 0)
vals = [float(det_seg_pi.subs(t, tv)) for tv in (0, 0.25, 0.5, 0.75, 1)]
check("C2c T4 numeric row reproduced (+1, .25, 0, .25, +1)",
      np.allclose(vals, [1, .25, 0, .25, 1]), f"{vals}")

# ----------------------------------------------------------------------
# C3. Strict (open-interval, det < 0) counterexample.
#     stretch: psi(x) = x + eta(r) x1 e1, eta(r) = exp(-r^2/R^2), R = 3
#       -> grad psi = I + e1 (eta' x1 x/r + eta e1)^T,
#          det grad psi = 1 + eta + eta' x1^2 / r   (rank-one update)
#     twist as above with theta == pi near 0.  phi = twist o stretch:
#       det grad phi = det(grad twist) * det(grad stretch) = det grad psi > 0,
#       F(0) = R_z(pi) diag(2,1,1) = diag(-2,-1,1).
#     Segment to identity at the origin: det F_t(0) = (1-3t)(1-2t) < 0 for
#     t in (1/3, 1/2): the segment leaves the admissible set through an OPEN
#     region of negative determinant (robust to perturbation).
# ----------------------------------------------------------------------
Rbig = sp.Symbol("R", positive=True)
eta = sp.exp(-r**2 / Rbig**2)
psi = sp.Matrix([x + eta * x, y, z])
det_psi = sp.simplify(psi.jacobian([x, y, z]).det())
# rank-one-update formula: det(I + e1 v^T) = 1 + v1,
# v1 = d(eta x)/dx = eta + eta'(r) x^2/r with eta'(r) = -2 r eta / R^2
det_psi_expected = sp.simplify(1 + eta + (-2 * r / Rbig**2) * eta * x**2 / r)
check("C3a det(grad stretch) == 1 + eta + eta' x^2/r [sympy]",
      sp.simplify(det_psi - det_psi_expected) == 0, f"det = {det_psi}")

# numeric global positivity of det(grad stretch) with R = 3 on a grid + rays
Rv = 3.0
N = 61
L = 12.0
g1 = np.linspace(-L, L, N)
X, Y, Z = np.meshgrid(g1, g1, g1, indexing="ij")
RR = np.sqrt(X**2 + Y**2 + Z**2) + 1e-300
ETA = np.exp(-RR**2 / Rv**2)
DETA = -2 * RR / Rv**2 * ETA
det_stretch = 1 + ETA + DETA * X**2 / RR
# analytic minimum: on the x-axis det = g(u) = 1 + e^-u (1-2u), u = x^2/R^2;
# g' = 0 at u = 3/2, g_min = 1 - 2 e^{-3/2} = 0.55374 > 0.
gmin_analytic = 1 - 2 * np.exp(-1.5)
check("C3b det(grad stretch) > 0 on [-12,12]^3 grid (min = 1-2e^-1.5)",
      det_stretch.min() > 0.55 and
      abs(det_stretch.min() - gmin_analytic) < 2e-3,
      f"min = {det_stretch.min():.6f}, analytic = {gmin_analytic:.6f}")
# worst radial line x-axis, dense:
xs = np.linspace(1e-9, 40, 400001)
det_line = 1 + np.exp(-xs**2 / Rv**2) * (1 - 2 * xs**2 / Rv**2)
check("C3c det(grad stretch) > 0 on x-axis dense scan (min = 1-2e^-1.5)",
      abs(det_line.min() - gmin_analytic) < 1e-6,
      f"min = {det_line.min():.6f}")

# F(0) of the composition and the segment determinant
F0 = np.array([[-2.0, 0, 0], [0, -1.0, 0], [0, 0, 1.0]])
# verify by finite differences of the composed map at the origin
eps = 1e-6


def stretch_map(p):
    rr = np.sqrt((p**2).sum())
    return p + np.exp(-rr**2 / Rv**2) * np.array([p[0], 0, 0])


def twist_map(p):
    rr = np.sqrt((p**2).sum())
    a = theta_prof(rr)
    ca, sa = np.cos(a), np.sin(a)
    return np.array([ca * p[0] - sa * p[1], sa * p[0] + ca * p[1], p[2]])


def composed(p):
    return twist_map(stretch_map(p))


Ffd = np.zeros((3, 3))
for j in range(3):
    e = np.zeros(3)
    e[j] = eps
    Ffd[:, j] = (composed(e) - composed(-e)) / (2 * eps)
check("C3d F(0) of twist o stretch = diag(-2,-1,1) [finite diff]",
      np.allclose(Ffd, F0, atol=1e-8), f"F(0) =\n{Ffd}")

tt = sp.Symbol("tt", real=True)
Fseg = (1 - tt) * sp.eye(3) + tt * sp.Matrix(F0)
det_Fseg = sp.expand(Fseg.det())
check("C3e segment det at origin == (1-3t)(1-2t) [sympy]",
      sp.simplify(det_Fseg - (1 - 3 * tt) * (1 - 2 * tt)) == 0,
      f"det = {sp.factor(det_Fseg)}")
d04 = float(det_Fseg.subs(tt, 0.4))
check("C3f det < 0 at t = 0.4 (open negative interval (1/3,1/2))",
      d04 < -0.03, f"det(0.4) = {d04:.4f}")

# composed map admissibility on the whole grid (chain rule: product of dets;
# twist det == 1 symbolically, so det grad phi = det grad stretch > 0):
check("C3g det(grad(twist o stretch)) = 1 * det(grad stretch) > 0 everywhere",
      det_stretch.min() > 0.55, "chain rule + C1 + C3b")

# ----------------------------------------------------------------------
# C4. Degree machinery: validate on hedgehog; twist lift degree = 0 along
#     the admissible path (integrand identically zero: image is a circle arc).
#     deg = (1/2 pi^2) int det[q, dq/dx, dq/dy, dq/dz] d^3x   (q in S^3)
# ----------------------------------------------------------------------
Nd = 120
Ld = 4.0
gd = (np.arange(Nd) + 0.5) / Nd * 2 * Ld - Ld   # cell-centred: avoids origin
dx = gd[1] - gd[0]
Xd, Yd, Zd = np.meshgrid(gd, gd, gd, indexing="ij")
Rd = np.sqrt(Xd**2 + Yd**2 + Zd**2)

# hedgehog control: q = (cos f, sin f * xhat), f = pi exp(-r^2): degree +-1
f = np.pi * np.exp(-Rd**2)
q0 = np.cos(f)
sf = np.sin(f) / Rd
q1, q2, q3 = sf * Xd, sf * Yd, sf * Zd
Q = np.stack([q0, q1, q2, q3], axis=-1)          # (N,N,N,4)
dQx = np.gradient(Q, dx, axis=0)
dQy = np.gradient(Q, dx, axis=1)
dQz = np.gradient(Q, dx, axis=2)
Mdet = np.stack([Q, dQx, dQy, dQz], axis=-2)     # (N,N,N,4,4)
deg_hedgehog = np.linalg.det(Mdet).sum() * dx**3 / (2 * np.pi**2)
check("C4a degree integrator control: hedgehog |deg| = 1 (tol 2%)",
      abs(abs(deg_hedgehog) - 1) < 0.02, f"deg = {deg_hedgehog:.4f}")

# twist texture lift along the admissible path theta_s = s*theta(r):
# q(x) = (cos(s*theta/2), 0, 0, sin(s*theta/2)) -- image is 1-dimensional,
# so the 4x4 determinant is identically zero: degree 0 for every s.
max_integrand = 0.0
for s_path in (0.25, 0.5, 0.75, 1.0):
    a2 = 0.5 * s_path * theta_prof(Rd)
    Qt = np.stack([np.cos(a2), 0 * Rd, 0 * Rd, np.sin(a2)], axis=-1)
    dQtx = np.gradient(Qt, dx, axis=0)
    dQty = np.gradient(Qt, dx, axis=1)
    dQtz = np.gradient(Qt, dx, axis=2)
    Mt = np.stack([Qt, dQtx, dQty, dQtz], axis=-2)
    max_integrand = max(max_integrand, np.abs(np.linalg.det(Mt)).max())
check("C4b twist-lift degree integrand == 0 along admissible path s in (0,1]",
      max_integrand < 1e-12, f"max |integrand| = {max_integrand:.2e}")

# ----------------------------------------------------------------------
# Verdict assembly
# ----------------------------------------------------------------------
n_pass = sum(1 for c in checks if c["pass"])
print(f"\ninternal checks: {n_pass}/{len(checks)} PASS")

readings = {
    "L1_printed_convexity": {
        "statement": "A = {u : det(1+grad u) > 0 pointwise, u -> 0 at inf} "
                     "is convex (App. F.1 as printed)",
        "verdict": "FALSE",
        "witness": "twist map (C1-C2): segment det -> 0 at t=1/2; strict "
                   "version twist o stretch (C3): det = (1-3t)(1-2t) < 0 on "
                   "(1/3,1/2). Not even star-shaped about u = 0 (C2).",
    },
    "L1p_contractibility_local_diffeos": {
        "statement": "same A (pointwise det F > 0, interpenetration allowed) "
                     "is contractible",
        "verdict": "NOT-PROVEN as printed; plausibly FALSE",
        "witness": "Gromov h-principle heuristic (T4 W7 repair-2, sketch "
                   "grade): pi_0 of the local-diffeo class ~ pi_3(SO(3)) = Z;"
                   " on nontrivial components deg R~ != 0.",
    },
    "L2_main_text_Diff_c": {
        "statement": "admissible = global orientation-preserving diffeos "
                     "decaying to identity; Diff_c(R^3) contractible",
        "verdict": "CONCLUSION TRUE (import: Cerf pi_0 = 0 / Hatcher "
                   "contractible), PRINTED PROOF STILL FALSE",
        "witness": "the twist map IS a global diffeo and the convex segment "
                   "still exits the class (C2/C3): convexity is not the "
                   "reason; a deep theorem is.",
    },
    "L4_history_path": {
        "statement": "along any continuous physical history u(t) with "
                     "det F(t,x) > 0, deg R~[u(t)] is constant",
        "verdict": "TRUE (elementary); sufficient for Theorem VIII'.1",
        "witness": "local constancy of degree under uniform convergence + "
                   "continuity of polar decomposition on GL+; demonstrated "
                   "on the admissible twist path (C4b) where the convex "
                   "path fails (C2).",
    },
}

section = {
    "task": "T-H8: Lemma II.1 / App. F.1 convexity",
    "verdict_one_line": (
        "printed convexity claim FALSE (machine-verified strict "
        "counterexample, det < 0 on an open t-interval); statement "
        "ill-posed as to the admissible class; conclusion true for the "
        "global-diffeo class (deep import) and, as actually used by "
        "Theorem VIII'.1, true under the elementary history-path repair"),
    "counterexample": {
        "map": "phi = twist(theta: pi inside r<=1, smoothstep to 0 on [1,2])"
               " o stretch(x -> x + exp(-r^2/9) x1 e1)",
        "F_at_origin": [[-2, 0, 0], [0, -1, 0], [0, 0, 1]],
        "segment_det_at_origin": "(1-3t)(1-2t)",
        "negative_interval": [1 / 3, 1 / 2],
        "det_at_t_0.4": d04,
    },
    "readings": readings,
    "repair": (
        "Lemma II.1' (corrected): let t -> u(t) be a continuous history "
        "with u(t) -> 0 at infinity and det F(t,x) >= delta > 0 on compacts. "
        "Then deg R~[u(t)] is constant in t; with u(0) = 0, deg R~ == 0 and "
        "deg P~ = deg Q~ along the history. Contrapositive: Delta Sigma K "
        "!= 0 along a history => inf_x det F -> 0. Proof: continuity of "
        "polar decomposition on GL+ + local constancy of Brouwer degree "
        "under uniform convergence. No convexity, contractibility, or "
        "connectedness of the admissible set is used; Theorem VIII'.1 "
        "(pairs-only) follows verbatim."),
    "checks": checks,
    "n_pass": n_pass,
    "n_total": len(checks),
}

summary = json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {}
summary.setdefault("phase", "H3.3")
summary["part_A_convexity"] = section
SUMMARY.write_text(json.dumps(summary, indent=1, sort_keys=False))
print(f"wrote {SUMMARY}")
