//! Complex-Hermitian eigensolve via the real-symmetric embedding, plus the
//! generalized reduction ω² M v = K v → M^(−1/2) K M^(−1/2).
//!
//! H = A + iB Hermitian (A symmetric, B antisymmetric) embeds as the real
//! symmetric 2n×2n matrix S = [[A, −B], [B, A]]. Its spectrum is that of H
//! with every eigenvalue EXACTLY doubled; any real eigenvector (x; y) of S
//! yields the complex eigenvector x + iy of H (nonzero automatically,
//! ‖x‖² + ‖y‖² = 1). Deduplication takes every second eigenvalue of the
//! ascending-sorted embedded spectrum and records the worst intra-pair
//! residual as a diagnostic. The solver is fs-la's deterministic cyclic
//! Jacobi (`jacobi_eigh`), so the whole path is replay-stable; fs-la's
//! complex QR `eig` serves as a values-only oracle.

use fs_la::eigen::jacobi_eigh;
use fs_la::eigen_complex::eig;
use fs_math::c64::C64;

/// Eigenmodes of a Hermitian pencil: ascending ω², eigenvectors, and the
/// worst intra-pair eigenvalue residual of the real embedding.
pub struct Modes {
    /// Eigenvalues ω², ascending.
    pub w2: Vec<f64>,
    /// Eigenvectors (one per eigenvalue). For the generalized solve these
    /// are q-space vectors M^(−1/2)·v (spectrum.py's `s[:, None] * V`).
    pub vecs: Vec<Vec<C64>>,
    /// max |λ_(2m+1) − λ_(2m)| over the deduplicated pairs.
    pub pair_resid: f64,
}

/// Standard Hermitian eigenproblem H v = λ v (row-major n×n) via the
/// 2n×2n real-symmetric embedding on `jacobi_eigh`.
#[must_use]
pub fn eigh_hermitian(h: &[C64], n: usize) -> Modes {
    assert_eq!(h.len(), n * n, "h must be n*n");
    let d = 2 * n;
    let mut s = vec![0.0f64; d * d];
    for i in 0..n {
        for j in 0..n {
            let z = h[i * n + j];
            s[i * d + j] = z.re;
            s[i * d + n + j] = -z.im;
            s[(n + i) * d + j] = z.im;
            s[(n + i) * d + n + j] = z.re;
        }
    }
    let (w, v) = jacobi_eigh(&s, d);
    let mut w2 = Vec::with_capacity(n);
    let mut vecs = Vec::with_capacity(n);
    let mut pair_resid = 0.0f64;
    for m in 0..n {
        let r = (w[2 * m + 1] - w[2 * m]).abs();
        if r > pair_resid {
            pair_resid = r;
        }
        w2.push(w[2 * m]);
        let col = 2 * m;
        let mut vc = Vec::with_capacity(n);
        for i in 0..n {
            vc.push(C64::new(v[i * d + col], v[(n + i) * d + col]));
        }
        vecs.push(vc);
    }
    Modes { w2, vecs, pair_resid }
}

/// Generalized eigenproblem ω² M v = K v with diagonal M, optionally on an
/// exactly-decoupled index sub-block (spectrum.py `solve`): reduce to
/// H = M^(−1/2) K M^(−1/2), solve, return q-space vectors M^(−1/2)·v.
#[must_use]
pub fn solve(kmat: &[C64; 36], mdiag: &[f64; 6], idx: Option<&[usize]>) -> Modes {
    let all: [usize; 6] = [0, 1, 2, 3, 4, 5];
    let idx: &[usize] = idx.unwrap_or(&all);
    let n = idx.len();
    let mut s = vec![0.0f64; n];
    for (a, &ia) in idx.iter().enumerate() {
        s[a] = 1.0 / mdiag[ia].sqrt();
    }
    let mut h = vec![C64::ZERO; n * n];
    for a in 0..n {
        for b in 0..n {
            h[a * n + b] = kmat[idx[a] * 6 + idx[b]].scale(s[a] * s[b]);
        }
    }
    // Enforce exact Hermiticity of the reduced block.
    for a in 0..n {
        for b in a..n {
            let z = (h[a * n + b] + h[b * n + a].conj()).scale(0.5);
            h[a * n + b] = z;
            h[b * n + a] = z.conj();
        }
    }
    let mut modes = eigh_hermitian(&h, n);
    for vc in &mut modes.vecs {
        for a in 0..n {
            vc[a] = vc[a].scale(s[a]);
        }
    }
    modes
}

/// Values-only oracle: run fs-la's complex QR `eig` on
/// H = M^(−1/2) K M^(−1/2) and compare against the embedding path.
/// Returns (max |Re λ_QR − ω²| over sorted spectra, max |Im λ_QR|).
#[must_use]
pub fn oracle_check(kmat: &[C64; 36], mdiag: &[f64; 6], w2: &[f64]) -> (f64, f64) {
    let mut s = [0.0f64; 6];
    for i in 0..6 {
        s[i] = 1.0 / mdiag[i].sqrt();
    }
    let mut h = vec![C64::ZERO; 36];
    for a in 0..6 {
        for b in 0..6 {
            h[a * 6 + b] = kmat[a * 6 + b].scale(s[a] * s[b]);
        }
    }
    let ev = eig(&h, 6).expect("complex QR oracle failed to converge");
    // eig sorts by (re, im); a Hermitian spectrum sorted by re matches our
    // ascending ω².
    let mut dre = 0.0f64;
    let mut dim = 0.0f64;
    for (i, z) in ev.iter().enumerate() {
        dre = dre.max((z.re - w2[i]).abs());
        dim = dim.max(z.im.abs());
    }
    (dre, dim)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::moduli::Moduli;
    use crate::symbol::{mass_matrix, stiffness};

    #[test]
    fn embedding_matches_oracle_and_residual() {
        let m = Moduli::bench();
        let k = stiffness([0.3, -0.2, 0.9], &m);
        let md = mass_matrix(&m);
        let modes = solve(&k, &md, None);
        assert_eq!(modes.w2.len(), 6);
        let scale = modes.w2[5].abs().max(1.0);
        assert!(modes.pair_resid <= 1e-12 * scale, "pair residual {}", modes.pair_resid);
        let (dre, dim) = oracle_check(&k, &md, &modes.w2);
        assert!(dre <= 1e-10 * scale && dim <= 1e-10 * scale, "oracle dev {dre} {dim}");
        // Eigenvector residual ‖K v − ω² M v‖ through the operator.
        for n in 0..6 {
            let v = &modes.vecs[n];
            let mut res = 0.0f64;
            for a in 0..6 {
                let mut kv = fs_math::c64::C64::ZERO;
                for b in 0..6 {
                    kv = kv + k[a * 6 + b] * v[b];
                }
                let r = kv - v[a].scale(modes.w2[n] * md[a]);
                res = res.max(r.abs());
            }
            assert!(res <= 1e-10 * scale, "mode {n} residual {res}");
        }
    }
}
