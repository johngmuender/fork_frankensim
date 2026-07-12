//! The stored field type: cell-centred N^3 grid of 4-component quaternions,
//! SoA flat `Vec<f64>`, fixed-vacuum ghost rind (2 layers), unit-norm
//! renormalisation with a documented deterministic sweep order, and the
//! analytic/profile seeds (compacton, hedgehog).
//!
//! Layout (matches the numpy `(4, P, P, P)` C-order arrays of
//! `field3d_solve.py`, P = n + 2*GHOST): component-major, k fastest:
//! `idx(a, ip, jp, kp) = ((a*P + ip)*P + jp)*P + kp` with padded indices
//! `ip = i + GHOST` for real cell i.  Cell centres are
//! `x_i = -half + (i + 0.5) h`, `h = 2 half / n` (fs-cosserat-pilot
//! convention; HALF = 2.4 for the compacton, LBOX = 4.5 for eps = 0.05).

use fs_math::det;

/// Ghost-rind width: 2 layers, enough for the 4th-order central stencil
/// (and a fortiori for the 2nd-order and compact corner schemes).
pub const GHOST: usize = 2;

/// Cell-centred N^3 stored quaternion field with a fixed-vacuum ghost rind.
pub struct Field3 {
    n: usize,
    half: f64,
    h: f64,
    p: usize,
    data: Vec<f64>,
}

impl Field3 {
    /// New field over `[-half, half]^3` with `n` cells per axis, filled with
    /// the vacuum (1, 0, 0, 0) everywhere (real cells and ghosts).
    #[must_use]
    pub fn new(n: usize, half: f64) -> Self {
        let p = n + 2 * GHOST;
        let mut data = vec![0.0; 4 * p * p * p];
        for v in data.iter_mut().take(p * p * p) {
            *v = 1.0; // component 0 = 1: vacuum fill
        }
        Field3 { n, half, h: 2.0 * half / (n as f64), p, data }
    }

    /// Cells per axis.
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    /// Box half-width.
    #[must_use]
    pub fn half(&self) -> f64 {
        self.half
    }

    /// Grid spacing h = 2 half / n.
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    /// Padded extent P = n + 2*GHOST per axis.
    #[must_use]
    pub fn padded(&self) -> usize {
        self.p
    }

    /// Cell-centre coordinate of real index i (same expression as the
    /// pilot: `-half + (i + 0.5) h`).
    #[must_use]
    pub fn x(&self, i: usize) -> f64 {
        -self.half + (i as f64 + 0.5) * self.h
    }

    #[inline]
    fn off(&self, a: usize, ip: usize, jp: usize, kp: usize) -> usize {
        ((a * self.p + ip) * self.p + jp) * self.p + kp
    }

    /// Read-only view of the raw SoA storage: `4 * P^3` doubles,
    /// component-major, k fastest — `idx(a, ip, jp, kp) =
    /// ((a*P + ip)*P + jp)*P + kp` over PADDED indices (see module docs).
    /// Added for `fs-gum-kern`'s tiled sweeps (Phase G1); the layout is
    /// already a documented contract of this type.
    #[inline]
    #[must_use]
    pub fn data(&self) -> &[f64] {
        &self.data
    }

    /// Mutable view of the raw SoA storage (layout as [`Self::data`]).
    /// Callers own the invariants (unit norm after their sweep, vacuum
    /// ghosts untouched) — exactly the same responsibility `set` gives
    /// them per cell.  Added for `fs-gum-kern` (Phase G1).
    #[inline]
    #[must_use]
    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }

    /// Component `a` at padded indices (ghosts addressable: real cell i is
    /// at padded index i + GHOST).
    #[inline]
    #[must_use]
    pub fn get_pad(&self, a: usize, ip: usize, jp: usize, kp: usize) -> f64 {
        self.data[self.off(a, ip, jp, kp)]
    }

    /// All four components at real cell (i, j, k).
    #[inline]
    #[must_use]
    pub fn get(&self, i: usize, j: usize, k: usize) -> [f64; 4] {
        let (ip, jp, kp) = (i + GHOST, j + GHOST, k + GHOST);
        [
            self.get_pad(0, ip, jp, kp),
            self.get_pad(1, ip, jp, kp),
            self.get_pad(2, ip, jp, kp),
            self.get_pad(3, ip, jp, kp),
        ]
    }

    /// Set all four components at real cell (i, j, k).
    #[inline]
    pub fn set(&mut self, i: usize, j: usize, k: usize, q: [f64; 4]) {
        let (ip, jp, kp) = (i + GHOST, j + GHOST, k + GHOST);
        for (a, &v) in q.iter().enumerate() {
            let o = self.off(a, ip, jp, kp);
            self.data[o] = v;
        }
    }

    /// Real-cell value as an `fs_ga::Quat` (API-boundary interop).
    #[must_use]
    pub fn quat(&self, i: usize, j: usize, k: usize) -> fs_ga::Quat {
        let q = self.get(i, j, k);
        fs_ga::Quat { w: q[0], x: q[1], y: q[2], z: q[3] }
    }

    /// Set a real cell from an `fs_ga::Quat`.
    pub fn set_quat(&mut self, i: usize, j: usize, k: usize, q: fs_ga::Quat) {
        self.set(i, j, k, [q.w, q.x, q.y, q.z]);
    }

    /// Reset every ghost cell to the vacuum (1, 0, 0, 0) — the
    /// decay-to-vacuum boundary policy ("fixed-vacuum ghost rind, 2 layers"
    /// of field3d).  Fixed (a, ip, jp, kp) sweep order.
    pub fn fill_ghosts_vacuum(&mut self) {
        let (n, p) = (self.n, self.p);
        let real = GHOST..n + GHOST;
        for a in 0..4 {
            let vac = if a == 0 { 1.0 } else { 0.0 };
            for ip in 0..p {
                for jp in 0..p {
                    for kp in 0..p {
                        if !(real.contains(&ip) && real.contains(&jp) && real.contains(&kp)) {
                            let o = self.off(a, ip, jp, kp);
                            self.data[o] = vac;
                        }
                    }
                }
            }
        }
    }

    /// Pointwise unit-norm renormalisation of all REAL cells (ghosts are
    /// untouched — they stay vacuum, which is unit-norm already).
    ///
    /// Deterministic-iteration-order contract: cells are swept in flat
    /// ascending (i, j, k) order; per cell the 4 components are loaded,
    /// nrm = sqrt(q0^2 + q1^2 + q2^2 + q3^2) (IEEE sqrt), each component is
    /// divided by nrm.  Returns max |nrm - 1| BEFORE the fix.
    pub fn renormalize(&mut self) -> f64 {
        let n = self.n;
        let mut drift = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let q = self.get(i, j, k);
                    let nrm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
                    let d = (nrm - 1.0).abs();
                    if d > drift {
                        drift = d;
                    }
                    self.set(i, j, k, [q[0] / nrm, q[1] / nrm, q[2] / nrm, q[3] / nrm]);
                }
            }
        }
        drift
    }

    /// Boundary-tail diagnostic: max |q_vec| / pi over the six outermost
    /// REAL faces (field3d gate: < 1e-6 validates the decay-to-vacuum ghost
    /// policy).
    #[must_use]
    pub fn boundary_tail(&self) -> f64 {
        let n = self.n;
        let mut m = 0.0_f64;
        let mut upd = |q: [f64; 4]| {
            let fb = (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
            if fb > m {
                m = fb;
            }
        };
        for j in 0..n {
            for k in 0..n {
                upd(self.get(0, j, k));
                upd(self.get(n - 1, j, k));
            }
        }
        for i in 0..n {
            for k in 0..n {
                upd(self.get(i, 0, k));
                upd(self.get(i, n - 1, k));
            }
        }
        for i in 0..n {
            for j in 0..n {
                upd(self.get(i, j, 0));
                upd(self.get(i, j, n - 1));
            }
        }
        m / core::f64::consts::PI
    }

    /// Seed the exact eps = 0 BPS compacton (trig-free closed form, ported
    /// verbatim from fs-cosserat-pilot): f0 = 2 arccos(r/R*) gives
    /// q0 = 2 r^2/R*^2 - 1 and q_i = (2/R*^2) sqrt(R*^2 - r^2) x_i inside
    /// R*, vacuum outside.  Ghosts stay vacuum (the exact analytic value
    /// there whenever half - 2h >= R*).
    pub fn sample_compacton(&mut self) {
        let rs2 = crate::rstar_sq();
        let n = self.n;
        for i in 0..n {
            let x = self.x(i);
            for j in 0..n {
                let y = self.x(j);
                for k in 0..n {
                    let z = self.x(k);
                    self.set(i, j, k, q_compacton(x, y, z, rs2));
                }
            }
        }
    }

    /// Seed the hedgehog q = (cos f, sin f x/r, sin f y/r, sin f z/r) from a
    /// stored radial profile via monotone linear interpolation (np.interp
    /// semantics: clamped outside the node range), with the family-A
    /// dilation q_d(x) = q_base(x/d).  Trig from `fs_math::det`.  Ghosts
    /// stay vacuum.
    pub fn sample_hedgehog(&mut self, r_nodes: &[f64], f_nodes: &[f64], d: f64) {
        let n = self.n;
        for i in 0..n {
            let x = self.x(i);
            for j in 0..n {
                let y = self.x(j);
                for k in 0..n {
                    let z = self.x(k);
                    let r = (x * x + y * y + z * z).sqrt();
                    let f = crate::radial::interp(r / d, r_nodes, f_nodes);
                    let q0 = det::cos(f);
                    // cell-centred grid never hits r = 0, but guard anyway
                    let sfr = if r > 1.0e-300 { det::sin(f) / r } else { 0.0 };
                    self.set(i, j, k, [q0, sfr * x, sfr * y, sfr * z]);
                }
            }
        }
    }
}

/// The exact eps = 0 compacton hedgehog, trig-free closed form
/// (fs-cosserat-pilot, verbatim).
#[must_use]
pub fn q_compacton(x: f64, y: f64, z: f64, rs2: f64) -> [f64; 4] {
    let r2 = x * x + y * y + z * z;
    if r2 >= rs2 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    let k = 2.0 * (rs2 - r2).sqrt() / rs2;
    [2.0 * r2 / rs2 - 1.0, k * x, k * y, k * z]
}
