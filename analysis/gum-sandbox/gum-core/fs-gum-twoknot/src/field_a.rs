//! Anisotropic (nx, ny, nz) stored quaternion field — the minimal non-cubic
//! generalization of `fs_gum_field::Field3` (which is cubic-only), living in
//! THIS crate so the finished measurement/solver crates stay untouched.
//!
//! Layout is the exact anisotropic analogue of Field3: component-major SoA,
//! k fastest, 2-layer fixed-vacuum ghost rind:
//! `idx(a, ip, jp, kp) = ((a*PX + ip)*PY + jp)*PZ + kp`, padded `ip = i + GHOST`.
//! Cell centres `x_i = -halfx + (i + 0.5) h` with ONE shared spacing h per
//! axis (the frozen h = 0.09375 protocol); `halfx = nx h / 2` etc.
//!
//! Determinism contract inherited verbatim: all sweeps are plain sequential
//! f64 loops in fixed ascending (a-outer where stated, i, j, k) order; the
//! only intrinsic is IEEE sqrt.  On nx = ny = nz this type is loop-for-loop
//! identical to Field3 (asserted bitwise by the twin gate binary).

pub use fs_gum_field::GHOST;

/// Cell-centred nx*ny*nz stored quaternion field, fixed-vacuum ghost rind.
pub struct FieldA {
    nx: usize,
    ny: usize,
    nz: usize,
    h: f64,
    px: usize,
    py: usize,
    pz: usize,
    data: Vec<f64>,
}

impl FieldA {
    /// New vacuum-filled field with `n* ` cells per axis at shared spacing h.
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize, h: f64) -> Self {
        let (px, py, pz) = (nx + 2 * GHOST, ny + 2 * GHOST, nz + 2 * GHOST);
        let mut data = vec![0.0; 4 * px * py * pz];
        for v in data.iter_mut().take(px * py * pz) {
            *v = 1.0; // component 0 = 1: vacuum fill
        }
        FieldA { nx, ny, nz, h, px, py, pz, data }
    }

    #[must_use]
    pub fn nx(&self) -> usize {
        self.nx
    }
    #[must_use]
    pub fn ny(&self) -> usize {
        self.ny
    }
    #[must_use]
    pub fn nz(&self) -> usize {
        self.nz
    }
    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }
    #[must_use]
    pub fn halfx(&self) -> f64 {
        0.5 * self.h * self.nx as f64
    }
    #[must_use]
    pub fn halfy(&self) -> f64 {
        0.5 * self.h * self.ny as f64
    }
    #[must_use]
    pub fn halfz(&self) -> f64 {
        0.5 * self.h * self.nz as f64
    }

    /// Cell-centre coordinates (Field3 expression per axis).
    #[must_use]
    pub fn x(&self, i: usize) -> f64 {
        -self.halfx() + (i as f64 + 0.5) * self.h
    }
    #[must_use]
    pub fn y(&self, j: usize) -> f64 {
        -self.halfy() + (j as f64 + 0.5) * self.h
    }
    #[must_use]
    pub fn z(&self, k: usize) -> f64 {
        -self.halfz() + (k as f64 + 0.5) * self.h
    }

    #[inline]
    fn off(&self, a: usize, ip: usize, jp: usize, kp: usize) -> usize {
        ((a * self.px + ip) * self.py + jp) * self.pz + kp
    }

    /// Component `a` at padded indices.
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

    /// Pointwise unit-norm renormalisation of all REAL cells (Field3
    /// contract: fixed ascending (i, j, k) sweep, IEEE sqrt, ghosts
    /// untouched).  Returns max |nrm - 1| BEFORE the fix.
    pub fn renormalize(&mut self) -> f64 {
        let mut drift = 0.0_f64;
        for i in 0..self.nx {
            for j in 0..self.ny {
                for k in 0..self.nz {
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
    /// REAL faces (Field3::boundary_tail, anisotropic bounds).
    #[must_use]
    pub fn boundary_tail(&self) -> f64 {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let mut m = 0.0_f64;
        let mut upd = |q: [f64; 4]| {
            let fb = (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
            if fb > m {
                m = fb;
            }
        };
        for j in 0..ny {
            for k in 0..nz {
                upd(self.get(0, j, k));
                upd(self.get(nx - 1, j, k));
            }
        }
        for i in 0..nx {
            for k in 0..nz {
                upd(self.get(i, 0, k));
                upd(self.get(i, ny - 1, k));
            }
        }
        for i in 0..nx {
            for j in 0..ny {
                upd(self.get(i, j, 0));
                upd(self.get(i, j, nz - 1));
            }
        }
        m / core::f64::consts::PI
    }
}

/// Flat index into an auxiliary REAL-cell buffer (anisotropic Grad layout):
/// `((a nx + i) ny + j) nz + k`.
#[inline]
#[must_use]
pub fn aidx_a(nx: usize, ny: usize, nz: usize, a: usize, i: usize, j: usize, k: usize) -> usize {
    ((a * nx + i) * ny + j) * nz + k
}

/// 2nd/4th-order central and corner stencils, anisotropic bounds; exact
/// ports of `fs_gum_field::stencil` (same expression order).
#[must_use]
pub fn derivs4_a(f: &FieldA, i: usize, j: usize, k: usize) -> [[f64; 4]; 3] {
    let c1 = 8.0 / (12.0 * f.h());
    let c2 = 1.0 / (12.0 * f.h());
    let (ip, jp, kp) = (i + GHOST, j + GHOST, k + GHOST);
    let mut d = [[0.0; 4]; 3];
    for c in 0..4 {
        d[0][c] = (f.get_pad(c, ip + 1, jp, kp) - f.get_pad(c, ip - 1, jp, kp)) * c1
            - (f.get_pad(c, ip + 2, jp, kp) - f.get_pad(c, ip - 2, jp, kp)) * c2;
        d[1][c] = (f.get_pad(c, ip, jp + 1, kp) - f.get_pad(c, ip, jp - 1, kp)) * c1
            - (f.get_pad(c, ip, jp + 2, kp) - f.get_pad(c, ip, jp - 2, kp)) * c2;
        d[2][c] = (f.get_pad(c, ip, jp, kp + 1) - f.get_pad(c, ip, jp, kp - 1)) * c1
            - (f.get_pad(c, ip, jp, kp + 2) - f.get_pad(c, ip, jp, kp - 2)) * c2;
    }
    d
}

/// Compact corner scheme at corner (ci, cj, ck), 0..=n per axis (exact
/// nested-average order of `fs_gum_field::stencil::corner`).
#[must_use]
pub fn corner_a(f: &FieldA, ci: usize, cj: usize, ck: usize) -> ([f64; 4], [[f64; 4]; 3]) {
    let ih = 1.0 / f.h();
    let (i0, j0, k0) = (ci + GHOST - 1, cj + GHOST - 1, ck + GHOST - 1);
    let mut qe = [0.0; 4];
    let mut d = [[0.0; 4]; 3];
    for a in 0..4 {
        let v = |di: usize, dj: usize, dk: usize| f.get_pad(a, i0 + di, j0 + dj, k0 + dk);
        let ay = |di: usize, dk: usize| 0.5 * (v(di, 1, dk) + v(di, 0, dk));
        let ayz = |di: usize| 0.5 * (ay(di, 1) + ay(di, 0));
        d[0][a] = (ayz(1) - ayz(0)) * ih;
        qe[a] = 0.5 * (ayz(1) + ayz(0));
        let ax = |dj: usize, dk: usize| 0.5 * (v(1, dj, dk) + v(0, dj, dk));
        let axz = |dj: usize| 0.5 * (ax(dj, 1) + ax(dj, 0));
        d[1][a] = (axz(1) - axz(0)) * ih;
        let axy = |dk: usize| 0.5 * (ax(1, dk) + ax(0, dk));
        d[2][a] = (axy(1) - axy(0)) * ih;
    }
    (qe, d)
}
