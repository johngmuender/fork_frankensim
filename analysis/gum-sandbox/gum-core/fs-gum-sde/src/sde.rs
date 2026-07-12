//! SDE problems and integrators.
//!
//! One problem trait ([`Sde`]: drift + diagonal diffusion), one integrator
//! trait ([`Stepper`]), two schemes:
//!
//! - [`EulerMaruyama`] — strong order 0.5 in general, 1.0 for ADDITIVE
//!   noise; weak order 1.
//! - [`Sra1`] — the additive-noise order-1.5 strong SRK (Rößler SRA1
//!   coefficients; Kloeden–Platen explicit order-1.5 class), using the
//!   correlated integral ΔZ ~ N(0, h³/3), cov(ΔW, ΔZ) = h²/2.
//!
//! # Why SRA1 and not Milstein
//!
//! The Milstein correction for diagonal noise is
//! `½ σ (∂σ/∂x) ((ΔW)² − h)`. Every SDE in this crate's scope — the OU
//! gates and the Nelson SDE with constant diffusion √(ħ/m) — has ADDITIVE
//! noise (`∂σ/∂x ≡ 0`), for which the correction vanishes identically and
//! Milstein degenerates to Euler–Maruyama exactly (survey dimension
//! "Nelson/stochastic-QM", missing-piece 6). The genuinely higher-order
//! scheme for additive noise is the SRA1-class SRK, which buys strong order
//! 1.5 at the price of one extra correlated Gaussian (ΔZ) per step per
//! component — two extra fixed draws, already reserved by the
//! [`crate::wiener`] addressing contract whether used or not.
//!
//! # SRA1 (additive noise dX = f(t,X) dt + g(t) dW, diagonal g)
//!
//! Rößler (2010) SRA1 tableau: c⁰ = (0, ¾), A⁰₂₁ = ¾, B⁰₂₁ = 3/2,
//! α = (⅓, ⅔), c¹ = (1, 0), β¹ = (1, 0), β² = (−1, 1):
//!
//! ```text
//! H2      = X + ¾ h f(t, X) + (3/2) g(t + h) (ΔZ/h)
//! X_{n+1} = X + h [⅓ f(t, X) + ⅔ f(t + ¾h, H2)]
//!             + g(t + h) (ΔW − ΔZ/h) + g(t) (ΔZ/h)
//! ```
//!
//! For constant g (OU, Nelson) the g-terms collapse to `g·ΔW`. SRA1's
//! order claims require additive noise: `diffusion` must not depend on `x`
//! (it MAY depend on `t`). Passing a state-dependent diffusion to [`Sra1`]
//! silently degrades the order — the convergence gates in the demo binary
//! are the guard.

use crate::wiener::WienerStream;

/// An Itô SDE with diagonal noise: dXᵢ = bᵢ(t, X) dt + σᵢ(t, X) dWᵢ.
pub trait Sde {
    /// State dimension.
    fn dim(&self) -> usize;
    /// Drift b(t, x) into `out` (len = dim).
    fn drift(&self, t: f64, x: &[f64], out: &mut [f64]);
    /// Diagonal diffusion σ(t, x) into `out` (len = dim). For [`Sra1`]
    /// validity σ must not depend on `x` (additive noise).
    fn diffusion(&self, t: f64, x: &[f64], out: &mut [f64]);
}

/// Reusable per-trajectory work buffers (no allocation inside step loops).
#[derive(Debug, Clone)]
pub struct Scratch {
    b1: Vec<f64>,
    b2: Vec<f64>,
    g0: Vec<f64>,
    g1: Vec<f64>,
    h2: Vec<f64>,
    dw: Vec<f64>,
    dz: Vec<f64>,
}

impl Scratch {
    /// Buffers for state dimension `dim`.
    #[must_use]
    pub fn new(dim: usize) -> Self {
        Scratch {
            b1: vec![0.0; dim],
            b2: vec![0.0; dim],
            g0: vec![0.0; dim],
            g1: vec![0.0; dim],
            h2: vec![0.0; dim],
            dw: vec![0.0; dim],
            dz: vec![0.0; dim],
        }
    }
}

/// One Euler–Maruyama step with EXPLICIT increments (ladder/gate entry
/// point): x += b(t,x) h + σ(t,x) ⊙ dw.
pub fn em_step_with(sde: &impl Sde, t: f64, h: f64, x: &mut [f64], dw: &[f64], s: &mut Scratch) {
    sde.drift(t, x, &mut s.b1);
    sde.diffusion(t, x, &mut s.g0);
    for i in 0..x.len() {
        x[i] += s.b1[i] * h + s.g0[i] * dw[i];
    }
}

/// One SRA1 step with EXPLICIT increments (see module docs for the tableau).
pub fn sra1_step_with(
    sde: &impl Sde,
    t: f64,
    h: f64,
    x: &mut [f64],
    dw: &[f64],
    dz: &[f64],
    s: &mut Scratch,
) {
    let n = x.len();
    sde.drift(t, x, &mut s.b1);
    sde.diffusion(t, x, &mut s.g0); // g(t)
    sde.diffusion(t + h, x, &mut s.g1); // g(t + h); x-independent by contract
    for i in 0..n {
        s.h2[i] = x[i] + 0.75 * h * s.b1[i] + 1.5 * s.g1[i] * (dz[i] / h);
    }
    sde.drift(t + 0.75 * h, &s.h2, &mut s.b2);
    for i in 0..n {
        let dzh = dz[i] / h;
        x[i] += h * (s.b1[i] / 3.0 + 2.0 * s.b2[i] / 3.0)
            + s.g1[i] * (dw[i] - dzh)
            + s.g0[i] * dzh;
    }
}

/// A one-step SDE integrator drawing its noise through the
/// [`WienerStream`] addressing contract.
pub trait Stepper {
    /// Advance `x` from `t` by `h`, using the increments addressed at
    /// `(traj, slot)`.
    fn step<S: Sde>(
        &self,
        sde: &S,
        ws: &WienerStream,
        traj: u32,
        slot: u64,
        t: f64,
        h: f64,
        x: &mut [f64],
        scratch: &mut Scratch,
    );
    /// Scheme name for gate tables.
    fn name(&self) -> &'static str;
}

/// Euler–Maruyama (strong 0.5 general / 1.0 additive; weak 1).
#[derive(Debug, Clone, Copy)]
pub struct EulerMaruyama;

impl Stepper for EulerMaruyama {
    fn step<S: Sde>(
        &self,
        sde: &S,
        ws: &WienerStream,
        traj: u32,
        slot: u64,
        t: f64,
        h: f64,
        x: &mut [f64],
        scratch: &mut Scratch,
    ) {
        ws.dw(traj, slot, h, &mut scratch.dw);
        // Split-borrow: move dw out of scratch view via raw indexing.
        let dw = std::mem::take(&mut scratch.dw);
        em_step_with(sde, t, h, x, &dw, scratch);
        scratch.dw = dw;
    }

    fn name(&self) -> &'static str {
        "euler-maruyama"
    }
}

/// Additive-noise SRA1-class order-1.5 SRK (see module docs; NOT Milstein —
/// Milstein degenerates to EM exactly for additive noise).
///
/// Order fine print (measured by the gate binary): the generic strong order
/// for additive noise is 1.5, limited by zero-mean residuals that carry
/// f''(x); on LINEAR drift (e.g. OU) those vanish and SRA1 superconverges
/// to strong order 2.0.
#[derive(Debug, Clone, Copy)]
pub struct Sra1;

impl Stepper for Sra1 {
    fn step<S: Sde>(
        &self,
        sde: &S,
        ws: &WienerStream,
        traj: u32,
        slot: u64,
        t: f64,
        h: f64,
        x: &mut [f64],
        scratch: &mut Scratch,
    ) {
        let mut dw = std::mem::take(&mut scratch.dw);
        let mut dz = std::mem::take(&mut scratch.dz);
        ws.dw_dz(traj, slot, h, &mut dw, &mut dz);
        sra1_step_with(sde, t, h, x, &dw, &dz, scratch);
        scratch.dw = dw;
        scratch.dz = dz;
    }

    fn name(&self) -> &'static str {
        "sra1"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// dX = -θX dt + σ dW.
    struct Ou {
        theta: f64,
        sigma: f64,
    }
    impl Sde for Ou {
        fn dim(&self) -> usize {
            1
        }
        fn drift(&self, _t: f64, x: &[f64], out: &mut [f64]) {
            out[0] = -self.theta * x[0];
        }
        fn diffusion(&self, _t: f64, _x: &[f64], out: &mut [f64]) {
            out[0] = self.sigma;
        }
    }

    #[test]
    fn zero_noise_reduces_to_deterministic_schemes() {
        let ou = Ou {
            theta: 1.0,
            sigma: 0.0,
        };
        let mut s = Scratch::new(1);
        let (dw, dz) = ([0.0], [0.0]);
        let h = 0.01;
        let mut x_em = [1.0];
        let mut x_sra = [1.0];
        let mut t = 0.0;
        for _ in 0..100 {
            em_step_with(&ou, t, h, &mut x_em, &dw, &mut s);
            sra1_step_with(&ou, t, h, &mut x_sra, &dw, &dz, &mut s);
            t += h;
        }
        let exact = (-1.0f64).exp();
        // EM = explicit Euler: O(h) global error; SRA1 drift tableau is a
        // 2-stage order-2 RK: O(h²).
        assert!((x_em[0] - exact).abs() < 2.5e-3);
        assert!((x_sra[0] - exact).abs() < 2.0e-5);
    }
}
