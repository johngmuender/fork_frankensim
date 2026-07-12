//! Mass-matrix-aware symplectic Verlet and time-domain plane-wave
//! validation — the repo's first TIME-DOMAIN check of the micropolar
//! second-order system ρ₀ü = −K_uu u − K_uφ φ, J φ̈ = −K_φu u − K_φφ φ.
//!
//! fs-time's `verlet_step` assumes unit mass; here the kick–drift–kick is
//! written with q̇ = M⁻¹p directly (p += h/2·F(q); q += h·M⁻¹p; p += h/2·F(q)),
//! which is the same variational integrator after the M^(1/2) change of
//! variables. The complex 6-dof plane-wave amplitude at one wavevector is
//! evolved as its 12-real embedding x = (Re q̂; Im q̂) with the embedded
//! stiffness S = [[A, −B], [B, A]] and mass diag(M, M).
//!
//! For a branch eigenvector K v = ω² M v, the embedded state x₀ = (Re v; Im v)
//! is a generalized eigenvector of (S, M_e), so Verlet keeps the motion in a
//! 2-plane and the M-weighted projection s_n = ⟨x_n, x₀⟩_M / ⟨x₀, x₀⟩_M is
//! exactly cos(n·θ) with the DISCRETE frequency ω_h = (2/h)·asin(h·ω/2).
//! The frequency estimator uses the exact three-term identity
//! s_{n+1} − 2s_n + s_{n−1} = −4 sin²(θ/2) s_n averaged over the whole run
//! (well-conditioned at small θ, unlike acos of a near-1 cosine).

use fs_math::c64::C64;

/// One mass-aware Störmer–Verlet step for M q̈ = force(q):
/// kick–drift–kick with q̇ = M⁻¹p. `minv` holds 1/M_ii; `force` writes
/// F(q) into its output slice.
pub fn verlet_mass_step<F: Fn(&[f64], &mut [f64])>(
    q: &mut [f64],
    p: &mut [f64],
    h: f64,
    minv: &[f64],
    force: &F,
    scratch: &mut [f64],
) {
    let n = q.len();
    force(q, scratch);
    for i in 0..n {
        p[i] = (0.5 * h).mul_add(scratch[i], p[i]);
    }
    for i in 0..n {
        q[i] = (h * minv[i]).mul_add(p[i], q[i]);
    }
    force(q, scratch);
    for i in 0..n {
        p[i] = (0.5 * h).mul_add(scratch[i], p[i]);
    }
}

/// Result of evolving one plane-wave branch mode.
pub struct ModeRun {
    /// Continuum eigensolve frequency ω = √ω².
    pub omega: f64,
    /// Discrete Verlet frequency ω_h = (2/h)·asin(h·ω/2).
    pub omega_h: f64,
    /// Frequency measured from the trajectory.
    pub omega_meas: f64,
    /// Initial energy H = ½pᵀM⁻¹p + ½xᵀSx.
    pub e0: f64,
    /// max |H_n − H_0| / H_0 over the run.
    pub e_dev_max: f64,
    /// |mean(H) over last window − mean(H) over first window| / H_0
    /// (secular-drift witness; symplecticity keeps this at round-off).
    pub e_drift: f64,
    /// Number of steps taken.
    pub steps: usize,
    /// Step size.
    pub h: f64,
}

/// Evolve the branch eigenmode (w2, v) of the pencil (K, M) for `steps`
/// Verlet steps of size `h`, measuring the oscillation frequency and the
/// energy drift.
#[must_use]
pub fn evolve_mode(
    kmat: &[C64; 36],
    mdiag: &[f64; 6],
    v: &[C64],
    w2: f64,
    steps: usize,
    h: f64,
) -> ModeRun {
    // Embedded stiffness S (12×12) and mass (12).
    let mut s = [0.0f64; 144];
    for i in 0..6 {
        for j in 0..6 {
            let z = kmat[i * 6 + j];
            s[i * 12 + j] = z.re;
            s[i * 12 + 6 + j] = -z.im;
            s[(6 + i) * 12 + j] = z.im;
            s[(6 + i) * 12 + 6 + j] = z.re;
        }
    }
    let mut mass = [0.0f64; 12];
    let mut minv = [0.0f64; 12];
    for i in 0..6 {
        mass[i] = mdiag[i];
        mass[6 + i] = mdiag[i];
        minv[i] = 1.0 / mdiag[i];
        minv[6 + i] = 1.0 / mdiag[i];
    }
    // Initial state: x0 = (Re v; Im v), M-normalized; p0 = 0.
    let mut x0 = [0.0f64; 12];
    for i in 0..6 {
        x0[i] = v[i].re;
        x0[6 + i] = v[i].im;
    }
    let mut mnorm = 0.0f64;
    for i in 0..12 {
        mnorm += mass[i] * x0[i] * x0[i];
    }
    let inv = 1.0 / mnorm.sqrt();
    for xi in &mut x0 {
        *xi *= inv;
    }
    let force = |q: &[f64], out: &mut [f64]| {
        for i in 0..12 {
            let mut acc = 0.0f64;
            for j in 0..12 {
                acc = s[i * 12 + j].mul_add(q[j], acc);
            }
            out[i] = -acc;
        }
    };
    let energy = |q: &[f64], p: &[f64]| -> f64 {
        let mut e = 0.0f64;
        for i in 0..12 {
            e += 0.5 * minv[i] * p[i] * p[i];
            let mut sq = 0.0f64;
            for j in 0..12 {
                sq = s[i * 12 + j].mul_add(q[j], sq);
            }
            e += 0.5 * q[i] * sq;
        }
        e
    };
    let mut q = x0;
    let mut p = [0.0f64; 12];
    let mut scratch = [0.0f64; 12];
    let e0 = energy(&q, &p);
    let mut e_dev_max = 0.0f64;
    let window = (steps / 20).max(1);
    let mut e_first = 0.0f64;
    let mut e_last = 0.0f64;
    // M-weighted projection series s_n = ⟨x_n, x0⟩_M (x0 is M-unit).
    let mut proj = Vec::with_capacity(steps + 1);
    let project = |q: &[f64]| -> f64 {
        let mut a = 0.0f64;
        for i in 0..12 {
            a += mass[i] * q[i] * x0[i];
        }
        a
    };
    proj.push(project(&q));
    for n in 0..steps {
        verlet_mass_step(&mut q, &mut p, h, &minv, &force, &mut scratch);
        proj.push(project(&q));
        let e = energy(&q, &p);
        let dev = ((e - e0) / e0).abs();
        if dev > e_dev_max {
            e_dev_max = dev;
        }
        if n < window {
            e_first += e;
        }
        if n >= steps - window {
            e_last += e;
        }
    }
    let e_drift = ((e_last - e_first) / (window as f64) / e0).abs();
    // sin²(θ/2) = −Σ s_n (s_{n+1} − 2 s_n + s_{n−1}) / (4 Σ s_n²).
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for n in 1..steps {
        let r = proj[n + 1] - 2.0 * proj[n] + proj[n - 1];
        num -= proj[n] * r;
        den += proj[n] * proj[n];
    }
    let sin_half = (num / (4.0 * den)).max(0.0).sqrt();
    let omega_meas = 2.0 * sin_half.asin() / h;
    let omega = w2.max(0.0).sqrt();
    let omega_h = 2.0 * (0.5 * h * omega).asin() / h;
    ModeRun { omega, omega_h, omega_meas, e0, e_dev_max, e_drift, steps, h }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eigh::solve;
    use crate::moduli::Moduli;
    use crate::symbol::{mass_matrix, stiffness};

    #[test]
    fn mode_rings_at_discrete_frequency() {
        let mut m = Moduli::bench();
        m.j = 2.5; // non-unit mass matrix actually exercised
        let kmat = stiffness([0.0, 0.0, 1.0], &m);
        let md = mass_matrix(&m);
        let modes = solve(&kmat, &md, None);
        let n = 3; // a gapped mode
        let run = evolve_mode(&kmat, &md, &modes.vecs[n], modes.w2[n], 20_000, 1e-3 / modes.w2[n].sqrt());
        assert!(
            ((run.omega_meas - run.omega_h) / run.omega_h).abs() < 1e-9,
            "measured {} vs discrete {}",
            run.omega_meas,
            run.omega_h
        );
        assert!(run.e_dev_max < 1e-5, "energy dev {}", run.e_dev_max);
    }
}
