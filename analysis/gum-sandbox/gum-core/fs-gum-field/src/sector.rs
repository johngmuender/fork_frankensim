//! Sector energies and topological degree on the stored field.
//!
//! All reductions are straight sequential f64 sums in fixed ascending
//! (i, j, k) order — the bit-determinism contract.  For `Scheme::Central2`
//! on a field whose stencil inputs equal the pilot's analytic evaluations
//! (e.g. the sampled compacton with a vacuum exterior), the accumulation
//! order below reproduces fs-cosserat-pilot's `grid_diag` bit-for-bit for
//! deg / E0 / E6 / I.

use crate::field::Field3;
use crate::stencil;

use core::f64::consts::PI;

/// Discretisation engine for the derivative sectors.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scheme {
    /// 2nd-order central differences at cell centres (pilot stencil).
    Central2,
    /// 4th-order central differences at cell centres (field3d measurement
    /// engine).
    Central4,
    /// Compact corner (midpoint) scheme at the (N+1)^3 corners (field3d
    /// descent engine; O(h^2), width-1 stencil).
    Corner,
}

/// Sector energies of one configuration (campaign conventions; see crate
/// docs).  `e0` and `i` are quadrature-exact cell sums, independent of the
/// scheme; `e2`, `e4`, `e6`, `deg` are evaluated at the scheme's points.
#[derive(Clone, Copy, Debug)]
pub struct Sectors {
    /// (1/4pi) INT sum_i |D_i q|^2
    pub e2: f64,
    /// (1/4pi) INT sum_{i<j} (|D_i|^2 |D_j|^2 - (D_i . D_j)^2)
    pub e4: f64,
    /// pi^3 INT b^2
    pub e6: f64,
    /// (1/4pi) INT (1 - q0)
    pub e0: f64,
    /// INT 2 (q1^2 + q2^2) (isorotation inertia, 3rd iso-axis)
    pub i: f64,
    /// INT b, b = -(1/2 pi^2) det[q, D_x q, D_y q, D_z q]
    pub deg: f64,
    /// t (E2 + E4) + E6 + E0 at the `t` passed to `measure`
    pub estat: f64,
}

pub(crate) fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

/// 4x4 determinant with rows (q, D_x q, D_y q, D_z q), first-row expansion
/// (fs-cosserat-pilot, verbatim; analytically equal to field3d's PAIRS
/// 6-term expansion).
#[must_use]
pub fn det4(r0: &[f64; 4], r1: &[f64; 4], r2: &[f64; 4], r3: &[f64; 4]) -> f64 {
    let minor = |r: &[f64; 4], skip: usize| -> [f64; 3] {
        let mut m = [0.0; 3];
        let mut j = 0;
        for (c, &v) in r.iter().enumerate() {
            if c != skip {
                m[j] = v;
                j += 1;
            }
        }
        m
    };
    r0[0] * det3(minor(r1, 0), minor(r2, 0), minor(r3, 0))
        - r0[1] * det3(minor(r1, 1), minor(r2, 1), minor(r3, 1))
        + r0[2] * det3(minor(r1, 2), minor(r2, 2), minor(r3, 2))
        - r0[3] * det3(minor(r1, 3), minor(r2, 3), minor(r3, 3))
}

struct Accum {
    s_e2: f64,
    s_e4: f64,
    s_deg: f64,
    s_b2: f64,
}

impl Accum {
    fn new() -> Self {
        Accum { s_e2: 0.0, s_e4: 0.0, s_deg: 0.0, s_b2: 0.0 }
    }

    /// One evaluation point: q (cell value or corner 8-average) and the
    /// three derivative quaternions.
    #[inline]
    fn point(&mut self, q: &[f64; 4], d: &[[f64; 4]; 3]) {
        let mut nn = [0.0_f64; 3];
        for ax in 0..3 {
            let v = &d[ax];
            nn[ax] = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
        }
        let dot = |i: usize, j: usize| -> f64 {
            d[i][0] * d[j][0] + d[i][1] * d[j][1] + d[i][2] * d[j][2] + d[i][3] * d[j][3]
        };
        let (d01, d02, d12) = (dot(0, 1), dot(0, 2), dot(1, 2));
        self.s_e2 += nn[0] + nn[1] + nn[2];
        self.s_e4 +=
            (nn[0] * nn[1] - d01 * d01) + (nn[0] * nn[2] - d02 * d02) + (nn[1] * nn[2] - d12 * d12);
        let b = -det4(q, &d[0], &d[1], &d[2]) / (2.0 * PI * PI);
        self.s_deg += b;
        self.s_b2 += b * b;
    }
}

/// Measure all sector energies of `f` with the given scheme at coupling `t`
/// (E_static = t (E2 + E4) + E6 + E0).  Plain h^3 sums, fixed (i, j, k)
/// ascending order for every accumulator.
#[must_use]
pub fn measure(f: &Field3, scheme: Scheme, t: f64) -> Sectors {
    let n = f.n();
    let h = f.h();
    let vol = h * h * h;
    let mut acc = Accum::new();
    match scheme {
        Scheme::Central2 => {
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let q = f.get(i, j, k);
                        let d = stencil::derivs2(f, i, j, k);
                        acc.point(&q, &d);
                    }
                }
            }
        }
        Scheme::Central4 => {
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let q = f.get(i, j, k);
                        let d = stencil::derivs4(f, i, j, k);
                        acc.point(&q, &d);
                    }
                }
            }
        }
        Scheme::Corner => {
            for ci in 0..=n {
                for cj in 0..=n {
                    for ck in 0..=n {
                        let (qe, d) = stencil::corner(f, ci, cj, ck);
                        acc.point(&qe, &d);
                    }
                }
            }
        }
    }
    // cell-density sectors (scheme-independent, quadrature-exact)
    let (mut s_e0, mut s_i) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                s_e0 += 1.0 - q[0];
                s_i += 2.0 * (q[1] * q[1] + q[2] * q[2]);
            }
        }
    }
    let fourpi = 4.0 * PI;
    let e2 = acc.s_e2 * vol / fourpi;
    let e4 = acc.s_e4 * vol / fourpi;
    let e6 = PI.powi(3) * acc.s_b2 * vol;
    let e0 = s_e0 * vol / fourpi;
    let deg = acc.s_deg * vol;
    let i_val = s_i * vol;
    Sectors { e2, e4, e6, e0, i: i_val, deg, estat: t * (e2 + e4) + e6 + e0 }
}

/// Weighted Derrick virial t E2 - t E4 - 3 E6 + 3 E0 (identically 0 on a
/// continuum stationary point).
#[must_use]
pub fn virial(s: &Sectors, t: f64) -> f64 {
    t * s.e2 - t * s.e4 - 3.0 * s.e6 + 3.0 * s.e0
}
