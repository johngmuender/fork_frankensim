//! Minimal stored-field 3-D grid container (the survey's `Field3`).
//!
//! Deliberately self-contained: a flat row-major `Vec<[f64; C]>` with a
//! fixed (i, j, k) order (the SAME loop order feeds every accumulation
//! and every certificate hash, per the pilot's determinism convention),
//! plus periodic or vacuum-padded ghost semantics. A sibling crate
//! (fs-gum-field) is building the production stored-field type
//! CONCURRENTLY; this container is intentionally minimal so a later
//! integration pass can unify the two without touching the algorithms.

/// Boundary condition of a stored field: periodic wrap, or a constant
/// vacuum ghost value outside the box (field3d_solve.py's `qp` padding —
/// vacuum ghost (1,0,0,0) for quaternions, z-hat for directors).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Bc<const C: usize> {
    /// Periodic wrap on all three axes.
    Periodic,
    /// Constant ghost value outside the stored box.
    Vacuum([f64; C]),
}

/// A C-component field stored on an n0 x n1 x n2 grid of nodes with
/// spacing `h`; node (i, j, k) sits at `origin + (i, j, k) * h`.
/// Row-major with k fastest: `idx = (i * n1 + j) * n2 + k`.
#[derive(Clone)]
pub struct Field3<const C: usize> {
    /// Nodes per axis.
    pub n: [usize; 3],
    /// Grid spacing (isotropic).
    pub h: f64,
    /// Coordinate of node (0, 0, 0).
    pub origin: [f64; 3],
    /// Boundary condition.
    pub bc: Bc<C>,
    /// Row-major node data, len = n0 * n1 * n2.
    pub data: Vec<[f64; C]>,
}

/// Unit-quaternion (S3) field.
pub type QField = Field3<4>;

/// Unit-vector (S2) director field. Fields meant as RP2 line fields are
/// simply consumed by the sign-gauge (disclination) diagnostics, which
/// never assume a global vector lift.
pub type DirectorField = Field3<3>;

impl<const C: usize> Field3<C> {
    /// Build from explicit parts (used to wrap derived data like the
    /// FFT-reconstructed vector potential).
    ///
    /// # Panics
    /// If `data.len() != n0 * n1 * n2`.
    #[must_use]
    pub fn from_parts(
        n: [usize; 3],
        h: f64,
        origin: [f64; 3],
        bc: Bc<C>,
        data: Vec<[f64; C]>,
    ) -> Self {
        assert_eq!(data.len(), n[0] * n[1] * n[2], "data length mismatch");
        Field3 { n, h, origin, bc, data }
    }

    /// Build by evaluating `f(x, y, z)` at every node, fixed (i, j, k)
    /// order.
    #[must_use]
    pub fn from_fn(
        n: [usize; 3],
        h: f64,
        origin: [f64; 3],
        bc: Bc<C>,
        f: impl Fn(f64, f64, f64) -> [f64; C],
    ) -> Self {
        let mut data = Vec::with_capacity(n[0] * n[1] * n[2]);
        for i in 0..n[0] {
            let x = origin[0] + i as f64 * h;
            for j in 0..n[1] {
                let y = origin[1] + j as f64 * h;
                for k in 0..n[2] {
                    let z = origin[2] + k as f64 * h;
                    data.push(f(x, y, z));
                }
            }
        }
        Field3 { n, h, origin, bc, data }
    }

    /// The campaign's standard cell-centred nn^3 grid over
    /// [-half, half]^3: h = 2 half / nn, node i at -half + (i + 0.5) h
    /// (the pilot's and field3d_solve.py's convention).
    #[must_use]
    pub fn cell_centered(
        nn: usize,
        half: f64,
        bc: Bc<C>,
        f: impl Fn(f64, f64, f64) -> [f64; C],
    ) -> Self {
        let h = 2.0 * half / nn as f64;
        let o = -half + 0.5 * h;
        Self::from_fn([nn; 3], h, [o; 3], bc, f)
    }

    /// Row-major index of node (i, j, k).
    #[inline]
    #[must_use]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.n[1] + j) * self.n[2] + k
    }

    /// Stored value at an in-range node.
    #[inline]
    #[must_use]
    pub fn at(&self, i: usize, j: usize, k: usize) -> [f64; C] {
        self.data[self.idx(i, j, k)]
    }

    /// BC-resolved value at any (possibly out-of-range) node index:
    /// periodic wrap, or the vacuum ghost outside the box.
    #[inline]
    #[must_use]
    pub fn get(&self, i: isize, j: isize, k: isize) -> [f64; C] {
        match self.bc {
            Bc::Periodic => {
                let w = |v: isize, n: usize| v.rem_euclid(n as isize) as usize;
                self.at(w(i, self.n[0]), w(j, self.n[1]), w(k, self.n[2]))
            }
            Bc::Vacuum(ghost) => {
                let inside = i >= 0
                    && (i as usize) < self.n[0]
                    && j >= 0
                    && (j as usize) < self.n[1]
                    && k >= 0
                    && (k as usize) < self.n[2];
                if inside {
                    self.at(i as usize, j as usize, k as usize)
                } else {
                    ghost
                }
            }
        }
    }

    /// Coordinate of node (i, j, k).
    #[inline]
    #[must_use]
    pub fn pos(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + i as f64 * self.h,
            self.origin[1] + j as f64 * self.h,
            self.origin[2] + k as f64 * self.h,
        ]
    }

    /// Trilinear interpolation at an arbitrary point, BC-resolved at the
    /// eight surrounding nodes.
    #[must_use]
    pub fn sample_trilinear(&self, p: [f64; 3]) -> [f64; C] {
        let g = [
            (p[0] - self.origin[0]) / self.h,
            (p[1] - self.origin[1]) / self.h,
            (p[2] - self.origin[2]) / self.h,
        ];
        let base = [g[0].floor(), g[1].floor(), g[2].floor()];
        let fr = [g[0] - base[0], g[1] - base[1], g[2] - base[2]];
        let (bi, bj, bk) = (base[0] as isize, base[1] as isize, base[2] as isize);
        let mut out = [0.0; C];
        for di in 0..2 {
            let wi = if di == 0 { 1.0 - fr[0] } else { fr[0] };
            for dj in 0..2 {
                let wj = if dj == 0 { 1.0 - fr[1] } else { fr[1] };
                for dk in 0..2 {
                    let wk = if dk == 0 { 1.0 - fr[2] } else { fr[2] };
                    let v = self.get(bi + di as isize, bj + dj as isize, bk + dk as isize);
                    let w = wi * wj * wk;
                    for c in 0..C {
                        out[c] += w * v[c];
                    }
                }
            }
        }
        out
    }
}

/// The Hopf map on one quaternion: n = R(q) z-hat with R the SO(3)
/// rotation represented by unit q (the survey's closed form).
#[inline]
#[must_use]
pub fn hopf_project_one(q: [f64; 4]) -> [f64; 3] {
    [
        2.0 * (q[1] * q[3] + q[0] * q[2]),
        2.0 * (q[2] * q[3] - q[0] * q[1]),
        1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]),
    ]
}

/// Director projection n = R(q) z-hat (S3 -> S2, the Hopf map) applied
/// node-wise to a stored quaternion field. Periodic BC carries over;
/// a vacuum ghost is itself Hopf-projected (the (1,0,0,0) vacuum maps
/// to the z-hat director vacuum).
#[must_use]
pub fn hopf_project(q: &QField) -> DirectorField {
    let bc = match q.bc {
        Bc::Periodic => Bc::Periodic,
        Bc::Vacuum(g) => Bc::Vacuum(hopf_project_one(g)),
    };
    let data = q.data.iter().map(|&v| hopf_project_one(v)).collect();
    Field3::from_parts(q.n, q.h, q.origin, bc, data)
}
