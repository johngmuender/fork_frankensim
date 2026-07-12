//! `WienerStream` — Brownian-increment addressing on fs-rand Philox streams.
//!
//! # The (trajectory, slot, component) → draw-index contract
//!
//! A Wiener increment is a pure function of `(seed, kernel, trajectory,
//! slot, component)` — never of execution order. The mapping is affine and
//! FIXED:
//!
//! ```text
//! stream key   = StreamKey { seed, kernel, tile: trajectory }
//! draw index   = (slot * dim + component) * DRAWS_PER_SLOT_COMPONENT
//! ```
//!
//! where `slot` is the integrator's step counter (or, for substepped
//! drivers, `step * MAX_SUB + substep` — the caller owns the slot
//! convention, this module owns the index arithmetic below it).
//!
//! Per (slot, component) exactly [`DRAWS_PER_SLOT_COMPONENT`] = 4 uniform
//! draws are RESERVED, covering two Box–Muller normals of 2 draws each
//! ([`DRAWS_PER_NORMAL`], the fs-rand `next_normal` contract):
//!
//! - `Z1` (draws 0–1): the ΔW normal, `ΔW = √h · Z1`;
//! - `Z2` (draws 2–3): the ΔZ partner for order-1.5 schemes,
//!   `ΔZ = ½ h^{3/2} (Z1 + Z2/√3)`, giving the exact joint law
//!   `ΔZ ~ N(0, h³/3)`, `cov(ΔW, ΔZ) = h²/2`.
//!
//! The reservation holds whether or not the integrator consumes `Z2`
//! (Euler–Maruyama reads only `Z1`; SRA1 reads both). Fixed consumption is
//! what makes the map affine, and the affine map is what buys the two
//! properties the spec demands:
//!
//! 1. **Order-independent ensembles**: trajectories can be stepped in any
//!    order (or in parallel, or grouped by substep count) and every
//!    trajectory sees bitwise-identical noise — fs-rand's shuffle-invariance
//!    property lifted to the Wiener layer.
//! 2. **Per-trajectory replay in isolation**: any single trajectory can be
//!    re-integrated alone (fresh caches, different ensemble, different
//!    machine) and reproduces its path bit-for-bit.
//!
//! Both properties are gated in the demo binary (gate G1).

use fs_math::det;
use fs_rand::{Stream, StreamKey};

/// Uniform draws consumed per Box–Muller normal (fs-rand `next_normal`
/// contract: exactly 2, always).
pub const DRAWS_PER_NORMAL: u64 = 2;

/// Normals reserved per (slot, component): Z1 (the ΔW normal) and Z2 (the
/// ΔZ partner). Both are RESERVED even when only Z1 is drawn — see the
/// module docs.
pub const NORMALS_PER_SLOT_COMPONENT: u64 = 2;

/// Uniform draws reserved per (slot, component) — the affine stride.
pub const DRAWS_PER_SLOT_COMPONENT: u64 = DRAWS_PER_NORMAL * NORMALS_PER_SLOT_COMPONENT;

/// Addressed Wiener-increment source for one ensemble: one Philox stream
/// per trajectory (`tile = trajectory`), increments random-access by
/// `(slot, component)`.
#[derive(Debug, Clone, Copy)]
pub struct WienerStream {
    seed: u64,
    kernel: u32,
    dim: u32,
}

impl WienerStream {
    /// New addressed source. `dim` is the per-trajectory state dimension
    /// (number of independent Wiener components per slot).
    #[must_use]
    pub fn new(seed: u64, kernel: u32, dim: u32) -> Self {
        assert!(dim > 0, "WienerStream needs at least one component");
        WienerStream { seed, kernel, dim }
    }

    /// State dimension (Wiener components per slot).
    #[must_use]
    pub fn dim(&self) -> u32 {
        self.dim
    }

    /// The affine draw-index map (exposed for contract tests).
    #[must_use]
    pub fn draw_index(&self, slot: u64, component: u32) -> u64 {
        debug_assert!(component < self.dim);
        (slot * u64::from(self.dim) + u64::from(component)) * DRAWS_PER_SLOT_COMPONENT
    }

    /// The raw normal pair (Z1, Z2) for one (trajectory, slot, component).
    #[must_use]
    pub fn normal_pair(&self, traj: u32, slot: u64, component: u32) -> (f64, f64) {
        let key = StreamKey {
            seed: self.seed,
            kernel: self.kernel,
            tile: traj,
        };
        let mut s = Stream::resume(key, self.draw_index(slot, component));
        let z1 = s.next_normal();
        let z2 = s.next_normal();
        (z1, z2)
    }

    /// Fill `dw[c] = √h · Z1(traj, slot, c)` — the Euler–Maruyama increment.
    /// Only Z1 is generated; Z2's index range stays reserved (unused), so EM
    /// and SRA1 runs share the same ΔW bits at every (traj, slot, c).
    pub fn dw(&self, traj: u32, slot: u64, h: f64, dw: &mut [f64]) {
        debug_assert_eq!(dw.len(), self.dim as usize);
        let sh = det::sqrt(h);
        for (c, w) in dw.iter_mut().enumerate() {
            let key = StreamKey {
                seed: self.seed,
                kernel: self.kernel,
                tile: traj,
            };
            let mut s = Stream::resume(key, self.draw_index(slot, c as u32));
            *w = sh * s.next_normal();
        }
    }

    /// Fill the coupled pair for order-1.5 schemes:
    /// `dw[c] = √h Z1`, `dz[c] = ½ h^{3/2} (Z1 + Z2/√3)`, so that
    /// `ΔZ ~ N(0, h³/3)` and `cov(ΔW, ΔZ) = h²/2` exactly.
    pub fn dw_dz(&self, traj: u32, slot: u64, h: f64, dw: &mut [f64], dz: &mut [f64]) {
        debug_assert_eq!(dw.len(), self.dim as usize);
        debug_assert_eq!(dz.len(), self.dim as usize);
        let sh = det::sqrt(h);
        let cz = 0.5 * h * sh; // ½ h^{3/2}
        let inv_sqrt3 = 1.0 / det::sqrt(3.0);
        for c in 0..self.dim as usize {
            let (z1, z2) = self.normal_pair(traj, slot, c as u32);
            dw[c] = sh * z1;
            dz[c] = cz * (z1 + z2 * inv_sqrt3);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn addressing_is_affine_and_reserved() {
        let ws = WienerStream::new(7, 3, 2);
        assert_eq!(ws.draw_index(0, 0), 0);
        assert_eq!(ws.draw_index(0, 1), 4);
        assert_eq!(ws.draw_index(1, 0), 8);
        assert_eq!(ws.draw_index(10, 1), 84);
    }

    #[test]
    fn em_and_sra1_share_dw_bits() {
        let ws = WienerStream::new(0xB4, 1, 3);
        let (mut a, mut w, mut z) = ([0.0; 3], [0.0; 3], [0.0; 3]);
        ws.dw(9, 1234, 0.01, &mut a);
        ws.dw_dz(9, 1234, 0.01, &mut w, &mut z);
        for c in 0..3 {
            assert_eq!(a[c].to_bits(), w[c].to_bits());
        }
    }

    #[test]
    fn dz_moments_match_the_joint_law() {
        // E[ΔZ²] = h³/3, E[ΔW ΔZ] = h²/2 (sampled, 3-sigma band).
        let ws = WienerStream::new(0x5EED, 2, 1);
        let h = 0.25;
        let n = 200_000u32;
        let (mut szz, mut swz, mut sww) = (0.0, 0.0, 0.0);
        let (mut w, mut z) = ([0.0], [0.0]);
        for traj in 0..n {
            ws.dw_dz(traj, 0, h, &mut w, &mut z);
            szz += z[0] * z[0];
            swz += w[0] * z[0];
            sww += w[0] * w[0];
        }
        let nf = f64::from(n);
        assert!((szz / nf / (h * h * h / 3.0) - 1.0).abs() < 0.02);
        assert!((swz / nf / (h * h / 2.0) - 1.0).abs() < 0.02);
        assert!((sww / nf / h - 1.0).abs() < 0.02);
    }
}
