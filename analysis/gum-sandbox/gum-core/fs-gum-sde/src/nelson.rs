//! The Nelson (stochastic mechanics) sector — the physics glue.
//!
//! Units: ħ = m = 1 throughout (the born.py convention). Writing
//! W = ∇ψ/ψ:
//!
//! - current velocity  v = (ħ/m) Im W = Im W   (de Broglie guidance —
//!   born.py's `velocity`);
//! - osmotic velocity  u = (ħ/2m) ∇ln ρ = (ħ/m) Re W = Re W
//!   (ρ = |ψ|² ⇒ ∇ln ρ = 2 Re W; the factor ½ and the factor 2 cancel —
//!   the classic trap, stated explicitly here);
//! - forward Nelson SDE:  dX = (v + u) dt + √(ħ/m) dW = (v + u) dt + dW,
//!   i.e. constant ADDITIVE diffusion σ = √(ħ/m) = 1, noise variance
//!   2ν h per component with ν = ħ/2m = ½.
//!
//! Its Fokker–Planck equation ∂ρ = −∇·((v+u)ρ) + ν∇²ρ has ρ = |ψ|² as
//! solution: equilibrium ensembles stay equilibrium (gate G3 — a wrong
//! factor of 2 in σ or u fails it immediately).
//!
//! Boundaries: the box walls are REFLECTING (specular reflection of the
//! post-step position, repeated until inside), then a strict clip to
//! [ε, π−ε] with ε = 1e-4 (born.py's wall clip). born.py's plain clip is
//! not adequate once noise can carry a particle through a wall in one
//! substep — reflection preserves the Neumann (zero-flux) boundary of the
//! Fokker–Planck equation.
//!
//! Substepping / displacement cap (born.py's guard structure adapted to
//! noise): the macro step dt is subdivided so the DRIFT displacement per
//! substep obeys |b| h ≤ [`DISP_MAX`] (up to [`MAX_SUB`] substeps), and a
//! hard cap rescales the drift displacement if a mid-substep node spike
//! still exceeds it. Deviation from born.py, documented: the cap applies
//! to the drift contribution ONLY — capping the Wiener increment would
//! bias the diffusion coefficient itself, so the noise is added uncapped
//! and the walls handle overshoot by reflection.
//!
//! Wiener addressing: substep j of macro step `step` draws at slot
//! `step * MAX_SUB + j` (fixed capacity per macro step whether used or
//! not), so per-trajectory paths are order-independent and replayable in
//! isolation regardless of each trajectory's substep count.

use std::collections::HashMap;

use fs_math::det;

use crate::modes::{ModeCoeffs, ModeSystem, BOX_L};
use crate::sde::Sde;
use crate::wiener::WienerStream;

/// Keep particles strictly inside (ε, π − ε) — born.py's `EPS_WALL`.
pub const EPS_WALL: f64 = 1e-4;
/// Max allowed drift displacement per substep — born.py's `DISP_MAX`.
pub const DISP_MAX: f64 = 0.05;
/// Cap on per-particle substeps per macro step — born.py's `MAX_SUB`.
/// Also the slot stride of the Wiener addressing (see module docs).
pub const MAX_SUB: u64 = 32;

/// The Nelson SDE on the 2-D box as a generic [`Sde`] (per-trajectory
/// state [x, y]). Recomputes the coefficient matrix on every drift call —
/// bit-identical to the cached ensemble driver ([`nelson_macro_step`]),
/// which shares C(t) across trajectories; used by the isolation-replay
/// gate.
#[derive(Debug, Clone)]
pub struct NelsonBox2D {
    /// The wavefunction.
    pub modes: ModeSystem,
    /// Diffusion amplitude √(ħ/m) = 1 in campaign units.
    pub sigma: f64,
}

impl NelsonBox2D {
    /// Standard campaign units: σ = √(ħ/m) = 1.
    #[must_use]
    pub fn new(modes: ModeSystem) -> Self {
        NelsonBox2D { modes, sigma: 1.0 }
    }
}

impl Sde for NelsonBox2D {
    fn dim(&self) -> usize {
        2
    }

    fn drift(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let c = self.modes.coeffs(t);
        let ld = self.modes.log_derivative(&c, x[0], x[1]);
        out[0] = ld.wx_im + ld.wx_re; // v + u
        out[1] = ld.wy_im + ld.wy_re;
    }

    fn diffusion(&self, _t: f64, _x: &[f64], out: &mut [f64]) {
        out[0] = self.sigma;
        out[1] = self.sigma;
    }
}

/// Specular reflection of a post-step coordinate into [0, l].
#[must_use]
pub fn reflect(mut x: f64, l: f64) -> f64 {
    loop {
        if x < 0.0 {
            x = -x;
        } else if x > l {
            x = 2.0 * l - x;
        } else {
            return x;
        }
    }
}

/// Per-macro-step cache of coefficient matrices keyed by t bit-pattern.
/// C(t) is a pure function of t, so sharing it across trajectories cannot
/// change any trajectory's bits (the isolation-replay gate proves it).
#[derive(Debug, Default)]
pub struct CoeffCache {
    map: HashMap<u64, ModeCoeffs>,
}

impl CoeffCache {
    /// Empty cache.
    #[must_use]
    pub fn new() -> Self {
        CoeffCache::default()
    }

    /// Drop all entries (call once per macro step to keep it small).
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// Coefficients at t, computed once per distinct t.
    pub fn get(&mut self, modes: &ModeSystem, t: f64) -> &ModeCoeffs {
        self.map.entry(t.to_bits()).or_insert_with(|| modes.coeffs(t))
    }
}

/// One Nelson macro step of size `dt` for a single trajectory
/// (Euler–Maruyama per substep; see the module docs for the guard,
/// reflection, and addressing rules). `step` is the global macro-step
/// counter; noise comes from `ws` at slots `step*MAX_SUB + j`.
#[allow(clippy::too_many_arguments)]
pub fn nelson_macro_step(
    modes: &ModeSystem,
    ws: &WienerStream,
    cache: &mut CoeffCache,
    traj: u32,
    step: u64,
    t: f64,
    dt: f64,
    x: &mut f64,
    y: &mut f64,
) {
    // Drift at the macro-step start decides the substep count (born.py's
    // guard: |b| h_sub ≤ DISP_MAX) and is reused as the first substep's
    // drift (born.py reuses k1 the same way).
    let (mut bx, mut by);
    {
        let c0 = cache.get(modes, t);
        let ld = modes.log_derivative(c0, *x, *y);
        bx = ld.wx_im + ld.wx_re;
        by = ld.wy_im + ld.wy_re;
    }
    let speed = det::sqrt(bx * bx + by * by);
    let nsub = ((speed * dt / DISP_MAX).ceil() as u64).clamp(1, MAX_SUB);
    let h = dt / nsub as f64;
    let mut dw = [0.0f64; 2];
    for j in 0..nsub {
        if j > 0 {
            let cj = cache.get(modes, t + j as f64 * h);
            let ld = modes.log_derivative(cj, *x, *y);
            bx = ld.wx_im + ld.wx_re;
            by = ld.wy_im + ld.wy_re;
        }
        // Drift displacement, hard-capped (drift ONLY — see module docs).
        let mut ddx = bx * h;
        let mut ddy = by * h;
        let d = det::sqrt(ddx * ddx + ddy * ddy);
        if d > DISP_MAX {
            let f = DISP_MAX / d;
            ddx *= f;
            ddy *= f;
        }
        ws.dw(traj, step * MAX_SUB + j, h, &mut dw);
        *x = reflect(*x + ddx + dw[0], BOX_L).clamp(EPS_WALL, BOX_L - EPS_WALL);
        *y = reflect(*y + ddy + dw[1], BOX_L).clamp(EPS_WALL, BOX_L - EPS_WALL);
    }
}

/// Rejection sample from ρ₀ = |φ₁₁|² = (2/π)² sin²x sin²y (born.py
/// `sample_phi11`), one point per trajectory from that trajectory's own
/// stream — deterministic consumption (rejections advance the index like
/// any draw; fs-rand doctrine), so initialization is order-independent too.
#[must_use]
pub fn sample_phi11(seed: u64, kernel: u32, traj: u32) -> (f64, f64) {
    let mut s = fs_rand::StreamKey {
        seed,
        kernel,
        tile: traj,
    }
    .stream();
    loop {
        let px = s.next_f64() * BOX_L;
        let py = s.next_f64() * BOX_L;
        let u = s.next_f64();
        let sx = det::sin(px);
        let sy = det::sin(py);
        if u < sx * sx * sy * sy {
            return (
                px.clamp(EPS_WALL, BOX_L - EPS_WALL),
                py.clamp(EPS_WALL, BOX_L - EPS_WALL),
            );
        }
    }
}

/// Rejection sample from |ψ(0)|² (born.py `sample_psi2`, the equilibrium
/// control). `c0` = `modes.coeffs(0.0)`, `bound` =
/// [`ModeSystem::psi2_max_bound`] at t = 0.
#[must_use]
pub fn sample_psi2(
    modes: &ModeSystem,
    c0: &ModeCoeffs,
    bound: f64,
    seed: u64,
    kernel: u32,
    traj: u32,
) -> (f64, f64) {
    let mut s = fs_rand::StreamKey {
        seed,
        kernel,
        tile: traj,
    }
    .stream();
    loop {
        let px = s.next_f64() * BOX_L;
        let py = s.next_f64() * BOX_L;
        let u = s.next_f64() * bound;
        if u < modes.psi2_point(c0, px, py) {
            return (
                px.clamp(EPS_WALL, BOX_L - EPS_WALL),
                py.clamp(EPS_WALL, BOX_L - EPS_WALL),
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflect_folds_into_box() {
        let l = BOX_L;
        assert!((reflect(-0.3, l) - 0.3).abs() < 1e-15);
        assert!((reflect(l + 0.2, l) - (l - 0.2)).abs() < 1e-15);
        assert!((reflect(0.5, l) - 0.5).abs() == 0.0);
        // Double fold.
        let x = reflect(2.0 * l + 0.1, l);
        assert!((0.0..=l).contains(&x));
    }

    #[test]
    fn trajectory_is_replayable_in_isolation_and_order_independent() {
        // Run a tiny ensemble forward, then re-run one trajectory alone
        // with a fresh cache and in a different sweep order: bitwise equal.
        let modes = ModeSystem::born_m16();
        let ws = WienerStream::new(42, 99, 2);
        let n = 8usize;
        let dt = 2e-3;
        let steps = 25u64;
        let mut xs: Vec<f64> = (0..n).map(|i| 0.3 + 0.05 * i as f64).collect();
        let mut ys: Vec<f64> = (0..n).map(|i| 2.5 - 0.07 * i as f64).collect();
        let (x0, y0) = (xs.clone(), ys.clone());
        let mut cache = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * dt;
            cache.clear();
            for i in 0..n {
                nelson_macro_step(
                    &modes, &ws, &mut cache, i as u32, step, t, dt, &mut xs[i], &mut ys[i],
                );
            }
        }
        // Isolation replay of trajectory 5.
        let (mut rx, mut ry) = (x0[5], y0[5]);
        let mut fresh = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * dt;
            fresh.clear();
            nelson_macro_step(&modes, &ws, &mut fresh, 5, step, t, dt, &mut rx, &mut ry);
        }
        assert_eq!(rx.to_bits(), xs[5].to_bits());
        assert_eq!(ry.to_bits(), ys[5].to_bits());
        // Reverse-order sweep.
        let mut xr = x0.clone();
        let mut yr = y0.clone();
        let mut cache2 = CoeffCache::new();
        for step in 0..steps {
            let t = step as f64 * dt;
            cache2.clear();
            for i in (0..n).rev() {
                nelson_macro_step(
                    &modes, &ws, &mut cache2, i as u32, step, t, dt, &mut xr[i], &mut yr[i],
                );
            }
        }
        for i in 0..n {
            assert_eq!(xr[i].to_bits(), xs[i].to_bits(), "traj {i} x");
            assert_eq!(yr[i].to_bits(), ys[i].to_bits(), "traj {i} y");
        }
    }
}
