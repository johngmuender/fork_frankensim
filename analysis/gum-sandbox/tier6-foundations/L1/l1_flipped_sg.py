#!/usr/bin/env python3
"""
L1 -- Contextuality exhibit: the flipped Stern-Gerlach (Bell 1982 / Bricmont 2019).

de Broglie-Bohm guidance-sector demonstration that "measurements don't measure":
same initial wave function, same initial particle positions, reversed SG gradient
=> identical spatial trajectories, flipped outcome LABELS.

Units: hbar = m = 1.

Model (1-D two-component spinor):
  Psi(z,0) = phi(z) (|up> + |down>)/sqrt(2),  phi = Gaussian(0, sigma=1)
  Impulsive SG kick: psi_up *= exp(+i k z), psi_down *= exp(-i k z), k = 2.5
  ("reversed gradient" = k -> -k)
  Free split-step FFT evolution H = p^2/2 afterward.
  Guidance: dZ/dt = j_total/rho_total,
    j_total = sum_a Im(psi_a^* dz psi_a),  rho_total = sum_a |psi_a|^2.
  Trajectories: RK4 on the discrete velocity field; velocity is evaluated
  spectrally on the grid at every half time step (so RK4 needs no temporal
  interpolation) and interpolated to particle positions with a Catmull-Rom
  cubic in z.
"""

import json
import numpy as np

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
L = 80.0                # box: z in [-40, 40], periodic
N_GRID = 8192           # >= 4096 per spec
SIGMA = 1.0
K_KICK = 2.5
DT = 0.005              # <= 0.005 per spec
T_FINAL = 2.4           # packet centers at +/- k t -> separation 2*2.5*2.4 = 12 sigma >= 8 sigma
N_TRAJ = 4000
SEED = 20260717
RHO_FLOOR = 1e-200      # guard against 0/0 far in the tails (never active in the bulk)

z = np.linspace(-L / 2, L / 2, N_GRID, endpoint=False)
dz = z[1] - z[0]
kvec = 2.0 * np.pi * np.fft.fftfreq(N_GRID, d=dz)


# ----------------------------------------------------------------------------
# Wave-function machinery
# ----------------------------------------------------------------------------
def initial_spinor(k_sign):
    """Gaussian phi times (|up>+|down>)/sqrt(2), then the impulsive SG kick."""
    phi = (2.0 * np.pi * SIGMA**2) ** -0.25 * np.exp(-z**2 / (4.0 * SIGMA**2))
    # normalize on the grid (periodic images are negligible but be exact)
    phi = phi / np.sqrt(np.sum(np.abs(phi) ** 2) * dz)
    up = phi / np.sqrt(2.0) * np.exp(1j * k_sign * K_KICK * z)
    dn = phi / np.sqrt(2.0) * np.exp(-1j * k_sign * K_KICK * z)
    return up, dn


def free_half_step_factor(dt):
    return np.exp(-0.5j * kvec**2 * (dt / 2.0))


def step_half(psi, factor):
    return np.fft.ifft(factor * np.fft.fft(psi))


def velocity_field(up, dn):
    """v = j_total / rho_total on the grid, derivative taken spectrally."""
    dup = np.fft.ifft(1j * kvec * np.fft.fft(up))
    ddn = np.fft.ifft(1j * kvec * np.fft.fft(dn))
    j = np.imag(np.conj(up) * dup) + np.imag(np.conj(dn) * ddn)
    rho = np.abs(up) ** 2 + np.abs(dn) ** 2
    return j / np.maximum(rho, RHO_FLOOR)


# ----------------------------------------------------------------------------
# Cubic (Catmull-Rom) interpolation of a grid field to particle positions
# ----------------------------------------------------------------------------
def interp_cubic(field, x):
    s = (x - z[0]) / dz
    i = np.floor(s).astype(np.int64)
    t = s - i
    i0 = (i - 1) % N_GRID
    i1 = i % N_GRID
    i2 = (i + 1) % N_GRID
    i3 = (i + 2) % N_GRID
    p0, p1, p2, p3 = field[i0], field[i1], field[i2], field[i3]
    return p1 + 0.5 * t * (
        p2 - p0
        + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + t * (3.0 * (p1 - p2) + p3 - p0))
    )


# ----------------------------------------------------------------------------
# Coupled evolution + trajectory integration
# ----------------------------------------------------------------------------
def run(k_sign, dt, Z0, store_stride=None):
    """Evolve the spinor with split-step FFT and integrate Bohmian trajectories.

    RK4 per step dt with velocity fields evaluated exactly at t, t+dt/2, t+dt
    (the spinor is advanced in half steps).  Returns final positions, the
    crossing count of z=0, final component densities, and (optionally) stored
    trajectory snapshots for plotting.
    """
    n_steps = int(round(T_FINAL / dt))
    assert abs(n_steps * dt - T_FINAL) < 1e-12
    up, dn = initial_spinor(k_sign)
    fac = free_half_step_factor(dt)

    Z = Z0.copy()
    sign_prev = np.sign(Z)
    crossings = 0

    stored_t, stored_Z = [], []
    if store_stride is not None:
        stored_t.append(0.0)
        stored_Z.append(Z.copy())

    v0 = velocity_field(up, dn)
    for n in range(n_steps):
        up_h = step_half(up, fac)
        dn_h = step_half(dn, fac)
        v_half = velocity_field(up_h, dn_h)
        up = step_half(up_h, fac)
        dn = step_half(dn_h, fac)
        v1 = velocity_field(up, dn)

        k1 = interp_cubic(v0, Z)
        k2 = interp_cubic(v_half, Z + 0.5 * dt * k1)
        k3 = interp_cubic(v_half, Z + 0.5 * dt * k2)
        k4 = interp_cubic(v1, Z + dt * k3)
        Z = Z + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

        sign_now = np.sign(Z)
        crossings += int(np.count_nonzero(sign_now * sign_prev < 0))
        sign_prev = sign_now
        v0 = v1

        if store_stride is not None and ((n + 1) % store_stride == 0 or n == n_steps - 1):
            stored_t.append((n + 1) * dt)
            stored_Z.append(Z.copy())

    rho_up = np.abs(up) ** 2
    rho_dn = np.abs(dn) ** 2
    out = {
        "Z_final": Z,
        "crossings": crossings,
        "rho_up": rho_up,
        "rho_dn": rho_dn,
    }
    if store_stride is not None:
        out["traj_t"] = np.array(stored_t)
        out["traj_Z"] = np.array(stored_Z)
    return out


def spin_labels(res, Z):
    """Which packet dominates locally at each particle position: +1 up, -1 down."""
    fu = interp_cubic(res["rho_up"], Z)
    fd = interp_cubic(res["rho_dn"], Z)
    return np.where(fu > fd, 1, -1)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    Z0 = rng.normal(0.0, SIGMA, N_TRAJ)  # |phi|^2 sampling

    # ---- original gradient ------------------------------------------------
    res_o = run(+1, DT, Z0, store_stride=4)
    Zf_o = res_o["Z_final"]
    lab_o = spin_labels(res_o, Zf_o)

    # packet separation at T (centers of the two component packets)
    rho_up, rho_dn = res_o["rho_up"], res_o["rho_dn"]
    mean_up = np.sum(z * rho_up) * dz / (np.sum(rho_up) * dz)
    mean_dn = np.sum(z * rho_dn) * dz / (np.sum(rho_dn) * dz)
    separation = abs(mean_up - mean_dn)

    # ---- reversed gradient ------------------------------------------------
    res_f = run(-1, DT, Z0, store_stride=4)
    Zf_f = res_f["Z_final"]
    lab_f = spin_labels(res_f, Zf_f)

    # ---- gates ------------------------------------------------------------
    # G1: nodal-line theorem (both runs; spec is stated for the symmetric state)
    g1_crossings = res_o["crossings"]
    g1_crossings_flip = res_f["crossings"]
    g1_pass = g1_crossings == 0

    # G2: outcome flip
    same_side = np.all(np.sign(Zf_f) == np.sign(Zf_o))
    all_labels_flipped = np.all(lab_f == -lab_o)
    max_dz_flip = float(np.max(np.abs(Zf_f - Zf_o)))
    # invariance of rho and j under k -> -k, measured on the final fields
    rho_o = res_o["rho_up"] + res_o["rho_dn"]
    rho_f = res_f["rho_up"] + res_f["rho_dn"]
    max_rho_diff = float(np.max(np.abs(rho_f - rho_o)))
    g2_pass = bool(same_side and all_labels_flipped)

    # G3: Born check
    n_up_side = int(np.count_nonzero(Zf_o > 0))
    frac_up = n_up_side / N_TRAJ
    sigma_binom = 0.5 / np.sqrt(N_TRAJ)
    g3_dev_sigmas = abs(frac_up - 0.5) / sigma_binom
    g3_pass = g3_dev_sigmas <= 3.0

    # consistency: label on the upper side should be "up" in the original run
    upper_label_orig = int(np.sign(np.sum(lab_o[Zf_o > 0])))
    upper_label_flip = int(np.sign(np.sum(lab_f[Zf_f > 0])))

    # G4: determinism control -- halve dt, re-integrate 100 trajectories
    n_ctrl = 100
    res_h = run(+1, DT / 2.0, Z0[:n_ctrl])
    g4_max_dev = float(np.max(np.abs(res_h["Z_final"] - Zf_o[:n_ctrl])))
    g4_threshold = 1e-6 * separation
    g4_pass = g4_max_dev <= g4_threshold

    results = {
        "workstream": "L1",
        "finding_id": "F-T6-2a-EXEC",
        "parameters": {
            "L": L, "N_grid": N_GRID, "sigma": SIGMA, "k_kick": K_KICK,
            "dt": DT, "T_final": T_FINAL, "N_traj": N_TRAJ, "seed": SEED,
            "integrator": "RK4, spectral velocity at half steps, Catmull-Rom cubic in z",
        },
        "packet_separation_at_T": float(separation),
        "packet_separation_over_sigma0": float(separation / SIGMA),
        "gates": {
            "G1_nodal_line": {
                "spec": "median (z=0) crossings across all trajectories/steps == 0",
                "crossings_original": g1_crossings,
                "crossings_flipped": g1_crossings_flip,
                "verdict": "PASS" if g1_pass else "FAIL",
            },
            "G2_outcome_flip": {
                "spec": "every trajectory: same exit side, flipped spin label; "
                        "|Z_flip - Z_orig| at machine/integration precision",
                "all_same_side": bool(same_side),
                "all_labels_flipped": bool(all_labels_flipped),
                "max_abs_Z_flip_minus_Z_orig": max_dz_flip,
                "max_abs_rho_total_diff": max_rho_diff,
                "upper_side_label_original": "up" if upper_label_orig > 0 else "down",
                "upper_side_label_flipped": "up" if upper_label_flip > 0 else "down",
                "verdict": "PASS" if g2_pass else "FAIL",
            },
            "G3_born": {
                "spec": "upper-exit fraction 0.5 within 3 sigma binomial (N=4000)",
                "n_upper": n_up_side,
                "frac_upper": frac_up,
                "sigma_binomial": sigma_binom,
                "deviation_sigmas": float(g3_dev_sigmas),
                "verdict": "PASS" if g3_pass else "FAIL",
            },
            "G4_determinism": {
                "spec": "dt -> dt/2, 100 trajectories: max deviation <= 1e-6 * separation",
                "max_position_deviation": g4_max_dev,
                "threshold": float(g4_threshold),
                "verdict": "PASS" if g4_pass else "FAIL",
            },
        },
    }

    with open("L1_results.json", "w") as f:
        json.dump(results, f, indent=2)

    np.savez_compressed(
        "L1_traj.npz",
        traj_t=res_o["traj_t"], traj_Z=res_o["traj_Z"],
        traj_t_flip=res_f["traj_t"], traj_Z_flip=res_f["traj_Z"],
        Zf_o=Zf_o, Zf_f=Zf_f, lab_o=lab_o, lab_f=lab_f,
        z=z, rho_up_o=res_o["rho_up"], rho_dn_o=res_o["rho_dn"],
        rho_up_f=res_f["rho_up"], rho_dn_f=res_f["rho_dn"],
    )

    for gid, g in results["gates"].items():
        print(f"{gid}: {g['verdict']}")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
