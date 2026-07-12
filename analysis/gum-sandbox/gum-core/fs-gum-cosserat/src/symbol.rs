//! Plane-wave symbol assembly for ARBITRARY k ∈ R³.
//!
//! spectrum.py hardcodes k ∥ z (its `build_ops` writes only the ∂_z row);
//! here every gradient is i·k_i, so K(k) is assembled for any direction.
//! The z-axis special case reproduces spectrum.py entry-for-entry, and the
//! gates binary verifies rotation invariance of the spectrum (random k̂ at
//! fixed |k| gives identical eigenvalues to 1e-12).
//!
//! Conventions (identical to spectrum.py):
//!   rows of a rank-2 operator are flattened C-order, row = 3·i + j;
//!   e_ij = i k_i u_j − ε_ijl φ_l;  Γ_ij = i k_i φ_j;
//!   ψ_a  = φ_a − (i/2) ε_abc k_b u_c  (ψ = φ − ½ curl u);
//!   each energy term c·Σ|T|² contributes K += 2c R†R  (W = ½ q†Kq).
//!
//! Phase E2 extension: the full W_χ chiral transduction of corpus II.C
//! eq (2.2) — χ₁ e_kkΓ_ll + χ₂ e_(ij)Γ_(ij) + χ₃ e_[ij]Γ_[ij] — assembled
//! as Hermitized cross terms K += χ(A + A†), with the (unpinned) χ₂
//! trace-vs-deviatoric convention switch on `Moduli::chi_deviatoric`.
//! Cross-validated against the real-space module [`crate::wchi`] (gate
//! G22): for a plane wave the grid energy equals (V/4)·v̂†K_χ(k_eff)v̂.

use crate::moduli::{MassCase, Moduli};
use fs_math::c64::C64;

/// Levi-Civita ε_ijk.
#[must_use]
pub const fn eps(i: usize, j: usize, k: usize) -> f64 {
    match (i, j, k) {
        (0, 1, 2) | (1, 2, 0) | (2, 0, 1) => 1.0,
        (0, 2, 1) | (2, 1, 0) | (1, 0, 2) => -1.0,
        _ => 0.0,
    }
}

/// Linear operators T(q̂) as matrices (rows = tensor components, 6 columns
/// = dof (u_x, u_y, u_z, φ_x, φ_y, φ_z)).
pub struct Ops {
    /// Relative strain e_ij, 9×6 row-major (row = 3i+j).
    pub e: [C64; 54],
    /// Wryness Γ_ij, 9×6.
    pub g: [C64; 54],
    /// ψ = φ − ½ curl u, 3×6.
    pub psi: [C64; 18],
    /// φ, 3×6.
    pub phi: [C64; 18],
}

/// Build the operator matrices at wavevector k (generalizes spectrum.py's
/// `build_ops(k)` from k·e_z to arbitrary k).
#[must_use]
pub fn build_ops(k: [f64; 3]) -> Ops {
    let mut e = [C64::ZERO; 54];
    let mut g = [C64::ZERO; 54];
    let mut psi = [C64::ZERO; 18];
    let mut phi = [C64::ZERO; 18];
    for i in 0..3 {
        for j in 0..3 {
            let r = 3 * i + j;
            // ∂_i u_j → i k_i u_j
            e[r * 6 + j] = e[r * 6 + j] + C64::new(0.0, k[i]);
            // − ε_ijl φ_l
            for l in 0..3 {
                let s = eps(i, j, l);
                if s != 0.0 {
                    e[r * 6 + 3 + l] = e[r * 6 + 3 + l] - C64::from_re(s);
                }
            }
            // Γ_ij = ∂_i φ_j → i k_i φ_j
            g[r * 6 + 3 + j] = g[r * 6 + 3 + j] + C64::new(0.0, k[i]);
        }
    }
    for a in 0..3 {
        psi[a * 6 + 3 + a] = C64::ONE;
        // ψ_a = φ_a − (i/2) ε_abc k_b u_c
        for c in 0..3 {
            let mut s = 0.0;
            for b in 0..3 {
                s += eps(a, b, c) * k[b];
            }
            if s != 0.0 {
                psi[a * 6 + c] = psi[a * 6 + c] - C64::new(0.0, 0.5 * s);
            }
        }
        phi[a * 6 + 3 + a] = C64::ONE;
    }
    Ops { e, g, psi, phi }
}

/// Trace (1×6), symmetric part (9×6), antisymmetric part (9×6) of a
/// rank-2 operator (9×6, row = 3i+j) — spectrum.py `tensor_parts`.
#[must_use]
pub fn tensor_parts(r: &[C64; 54]) -> ([C64; 6], [C64; 54], [C64; 54]) {
    let mut tr = [C64::ZERO; 6];
    let mut sym = [C64::ZERO; 54];
    let mut asym = [C64::ZERO; 54];
    for c in 0..6 {
        tr[c] = r[c] + r[4 * 6 + c] + r[8 * 6 + c];
    }
    for i in 0..3 {
        for j in 0..3 {
            let rij = 3 * i + j;
            let rji = 3 * j + i;
            for c in 0..6 {
                sym[rij * 6 + c] = (r[rij * 6 + c] + r[rji * 6 + c]).scale(0.5);
                asym[rij * 6 + c] = (r[rij * 6 + c] - r[rji * 6 + c]).scale(0.5);
            }
        }
    }
    (tr, sym, asym)
}

/// Deviatoric part of a symmetric rank-2 operator: dev(S)_ij = S_ij −
/// (1/3)δ_ij S_kk (rows 0, 4, 8 shifted by −tr/3). Used only by the
/// deviatoric χ₂ convention — see [`crate::moduli::Moduli::chi_deviatoric`].
#[must_use]
pub fn deviatoric(sym: &[C64; 54], tr: &[C64; 6]) -> [C64; 54] {
    let mut d = *sym;
    for diag in [0usize, 4, 8] {
        for c in 0..6 {
            d[diag * 6 + c] = d[diag * 6 + c] - tr[c].scale(1.0 / 3.0);
        }
    }
    d
}

/// K += coeff · (A + A†) with A = Ra†Rb over `rows` rows — the Hermitized
/// cross term of a real-space coupling coeff·Re(a_ij* b_ij), so that
/// ½ q†Kq picks up coeff·Re(Σ a* b). Loop order identical to the original
/// inline χ₃ block (bit-identical results for the χ₃ path).
fn add_cross(k: &mut [C64; 36], ra: &[C64], rb: &[C64], rows: usize, coeff: f64) {
    if coeff == 0.0 {
        return;
    }
    for a in 0..6 {
        for b in 0..6 {
            let mut acc = C64::ZERO;
            for row in 0..rows {
                acc = acc + ra[row * 6 + a].conj() * rb[row * 6 + b];
            }
            k[a * 6 + b] = k[a * 6 + b] + acc.scale(coeff);
            k[b * 6 + a] = k[b * 6 + a] + acc.conj().scale(coeff);
        }
    }
}

/// K += coeff · R†R over `rows` rows (fixed row-major accumulation order).
fn add_gram(k: &mut [C64; 36], r: &[C64], rows: usize, coeff: f64) {
    if coeff == 0.0 {
        return;
    }
    for a in 0..6 {
        for b in 0..6 {
            let mut acc = C64::ZERO;
            for row in 0..rows {
                acc = acc + r[row * 6 + a].conj() * r[row * 6 + b];
            }
            k[a * 6 + b] = k[a * 6 + b] + acc.scale(coeff);
        }
    }
}

/// 6×6 Hermitian stiffness symbol K(k) with W = ½ q†Kq
/// (spectrum.py `stiffness`, generalized to arbitrary k).
#[must_use]
pub fn stiffness(kvec: [f64; 3], m: &Moduli) -> [C64; 36] {
    let ops = build_ops(kvec);
    let (etr, esym, easym) = tensor_parts(&ops.e);
    let (gtr, gsym, gasym) = tensor_parts(&ops.g);

    let mut k = [C64::ZERO; 36];
    add_gram(&mut k, &etr, 1, m.lam); //          (λ/2)|e_kk|²
    add_gram(&mut k, &esym, 9, 2.0 * m.mu); //    μ e_(ij) e_(ij)
    add_gram(&mut k, &easym, 9, 2.0 * m.mu_c); // μ_c e_[ij] e_[ij]
    add_gram(&mut k, &gtr, 1, m.alpha); //        (α/2)|Γ_kk|²
    add_gram(&mut k, &gsym, 9, m.beta); //        (β/2)Γ_(ij)²
    add_gram(&mut k, &gasym, 9, m.gamma); //      (γ/2)Γ_[ij]²

    let mv2 = m.m_v * m.m_v;
    match m.case {
        MassCase::A => add_gram(&mut k, &ops.psi, 3, mv2), // (m_V²/2)|ψ|²
        MassCase::B => add_gram(&mut k, &ops.phi, 3, mv2), // (m_V²/2)|φ|²
    }

    // W_χ chiral transduction, corpus II.C eq (2.2):
    //   W_χ = χ₁ e_kk Γ_ll + χ₂ e_(ij)Γ_(ij) + χ₃ e_[ij]Γ_[ij]
    // (each coupling FREE, default 0; the (ij) contraction of χ₂ is
    // trace-included or deviatoric per m.chi_deviatoric — the text never
    // pins the choice). All three couplings enter as Hermitized cross
    // terms K += χ (A + A†) so that ½q†Kq = χ Re(Σ e-part* Γ-part).
    add_cross(&mut k, &etr, &gtr, 1, m.chi1); //   χ₁ e_kk Γ_ll
    if m.chi_deviatoric {
        let edev = deviatoric(&esym, &etr);
        let gdev = deviatoric(&gsym, &gtr);
        add_cross(&mut k, &edev, &gdev, 9, m.chi2); // χ₂ dev(e)_(ij) dev(Γ)_(ij)
    } else {
        add_cross(&mut k, &esym, &gsym, 9, m.chi2); // χ₂ e_(ij) Γ_(ij)
    }
    add_cross(&mut k, &easym, &gasym, 9, m.chi3); // χ₃ e_[ij] Γ_[ij]

    // Enforce exact Hermiticity against round-off.
    for a in 0..6 {
        for b in a..6 {
            let z = (k[a * 6 + b] + k[b * 6 + a].conj()).scale(0.5);
            k[a * 6 + b] = z;
            k[b * 6 + a] = z.conj();
        }
    }
    k
}

/// M = diag(ρ₀, ρ₀, ρ₀, J, J, J).
#[must_use]
pub fn mass_matrix(m: &Moduli) -> [f64; 6] {
    [m.rho0, m.rho0, m.rho0, m.j, m.j, m.j]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn max_abs_diff(a: &[C64; 36], b: &[C64; 36]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| (*x - *y).abs()).fold(0.0, f64::max)
    }

    #[test]
    fn z_axis_closed_forms() {
        // K along z: the u_z diagonal entry is (λ+2μ)k², the φ_z entry is
        // (α+β)k² + 4μ_c + m_V² (Case A gives ψ_z = φ_z).
        let m = Moduli::bench();
        let kv = 0.7;
        let k = stiffness([0.0, 0.0, kv], &m);
        let kuu = k[2 * 6 + 2];
        let kpp = k[5 * 6 + 5];
        let exp_u = (m.lam + 2.0 * m.mu) * kv * kv;
        let exp_p = (m.alpha + m.beta) * kv * kv + 4.0 * m.mu_c + m.m_v * m.m_v;
        assert!((kuu.re - exp_u).abs() < 1e-14 && kuu.im.abs() < 1e-16, "{kuu:?}");
        assert!((kpp.re - exp_p).abs() < 1e-14 && kpp.im.abs() < 1e-16, "{kpp:?}");
    }

    #[test]
    fn psi_identity_folds_mass_into_mu_c() {
        // e_[ij]e_[ij] = 2|ψ|² ⇒ Case A with (μ_c, m_V) equals Case A with
        // (μ_c + m_V²/4, 0) exactly.
        let mut m1 = Moduli::bench();
        m1.m_v = 2.0;
        let mut m2 = Moduli::bench();
        m2.mu_c += 1.0; // m_V²/4 = 1
        m2.m_v = 0.0;
        for kvec in [[0.3, -0.4, 1.1], [0.0, 0.0, 2.0], [1.0, 0.0, 0.0]] {
            let k1 = stiffness(kvec, &m1);
            let k2 = stiffness(kvec, &m2);
            assert!(max_abs_diff(&k1, &k2) < 1e-12, "identity broken at {kvec:?}");
        }
    }

    #[test]
    fn hermitian_by_construction() {
        let mut m = Moduli::bench();
        m.chi3 = 0.3;
        let k = stiffness([0.2, 0.5, 1.3], &m);
        for a in 0..6 {
            for b in 0..6 {
                let d = k[a * 6 + b] - k[b * 6 + a].conj();
                assert!(d.abs() == 0.0, "K not exactly Hermitian at ({a},{b})");
            }
        }
    }

    #[test]
    fn chi12_hermitian_and_nontrivial() {
        let mut m = Moduli::bench();
        m.chi1 = 0.4;
        m.chi2 = -0.25;
        m.chi3 = 0.3;
        let k = stiffness([0.7, -0.2, 1.1], &m);
        let k0 = stiffness([0.7, -0.2, 1.1], &Moduli::bench());
        let mut changed = false;
        for a in 0..6 {
            for b in 0..6 {
                let d = k[a * 6 + b] - k[b * 6 + a].conj();
                assert!(d.abs() == 0.0, "K not exactly Hermitian at ({a},{b})");
                changed |= (k[a * 6 + b] - k0[a * 6 + b]).abs() > 1e-12;
            }
        }
        assert!(changed, "chi1/chi2 must alter K");
    }

    #[test]
    fn deviatoric_is_chi1_shift() {
        // dev(e)_(ij) dev(Γ)_(ij) = e_(ij)Γ_(ij) − (1/3) e_kk Γ_ll, so
        // deviatoric(χ₁, χ₂) ≡ trace-included(χ₁ − χ₂/3, χ₂) exactly.
        let mut md = Moduli::bench();
        md.chi1 = 0.2;
        md.chi2 = 0.9;
        md.chi_deviatoric = true;
        let mut mt = Moduli::bench();
        mt.chi1 = 0.2 - 0.9 / 3.0;
        mt.chi2 = 0.9;
        for kvec in [[0.3, -0.4, 1.1], [0.0, 0.0, 2.0], [1.0, 0.5, 0.0]] {
            let kd = stiffness(kvec, &md);
            let kt = stiffness(kvec, &mt);
            let scale = kt.iter().map(|z| z.abs()).fold(1.0, f64::max);
            assert!(
                max_abs_diff(&kd, &kt) < 1e-13 * scale,
                "convention reparametrization broken at {kvec:?}"
            );
        }
    }
}
