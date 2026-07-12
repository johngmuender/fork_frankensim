//! Real-space W_χ chiral-coupling energy on a periodic grid (Phase E2).
//!
//! Implements the chiral-transduction term of the GUM constitutive action
//! for STORED micropolar fields (u, φ) on a periodic cell-centred N³ grid:
//!
//!   e_ij = ∂_i u_j − ε_ijk φ_k   (relative strain; parity tensor)
//!   Γ_ij = ∂_i φ_j               (wryness; pseudotensor)
//!   W_χ  = χ₁ e_kk Γ_ll + χ₂ e_(ij)Γ_(ij) + χ₃ e_[ij]Γ_[ij]
//!
//! # What the corpus text PINS (and where)
//!
//! - The three-coupling FORM of W_χ, with coefficient 1 on each contraction
//!   (no ½), built from e and Γ exactly as above: corpus
//!   `01-GUM-Omega-Paper-v2.0.1.md` Sec. II.C eq (2.2) ("chiral
//!   transduction, F2"). Because e is a tensor and Γ a pseudotensor, W_χ is
//!   parity-odd, and it lives on the MICROPOLAR (u, φ) fields — not on the
//!   quaternion knot texture.
//! - Units class [χ_i] = E·L⁻¹ (App. A per the survey entry).
//! - App. C.6 / Thm II.3: on-shell locking quenches χ₃ on the B2 branch at
//!   tree level; regeneration is (ka)²-suppressed (the Tier-1 k⁵ light
//!   splitting, gate G13, is the validated phenomenology).
//!
//! # What the text leaves FREE (and how this module handles it)
//!
//! - The VALUES of χ₁, χ₂, χ₃: free runtime parameters on [`Moduli`],
//!   **defaulting to 0**. At the default the density, total energy, and
//!   gradient below are all EXACTLY zero (hard guard, gate G23), so the
//!   campaign's χ-free knot-statics functional E_static = t(E2+E4)+E6+E0
//!   is untouched by this module's existence.
//! - Whether the (ij) contraction of the χ₂ term is trace-included or
//!   deviatoric: the text never says (survey missing-pieces item). Both
//!   conventions are implemented behind [`Moduli::chi_deviatoric`]
//!   (default: trace-included, matching this crate's W₂ convention). They
//!   are an exact reparametrization — dev·dev = sym·sym − (1/3)tr·tr, so
//!   deviatoric(χ₁, χ₂) ≡ trace-included(χ₁ − χ₂/3, χ₂) — gated at G23.
//! - Any coupling of W_χ to the quaternion knot sector: none appears in
//!   any solved statics functional anywhere in the corpus or campaign.
//!
//! # Discretization convention (this module's, documented — not pinned)
//!
//! Second-order central differences on the periodic grid:
//! (D_i g)(x) = [g(x + h ê_i) − g(x − h ê_i)]/(2h). The discrete plane
//! wave e^{ik·x} is an EXACT eigenfunction of D_i with symbol i·k_eff,
//! k_eff_i = sin(k_i h)/h, so the plane-wave consistency gate (G22)
//! compares the grid energy against (V/4)·v̂†K_χ(k_eff)v̂ built by
//! [`crate::symbol::stiffness`] from the same couplings — an exact
//! identity down to round-off, cross-validating both implementations.
//! (Factor: real fields u = Re(v̂ e^{ik·x}) average quadratics to half the
//! complex form; W = ½q̂†Kq̂ supplies the other half.)
//!
//! All reductions run in fixed (cell-major, component-minor) order — the
//! crate's bit-identical-replay discipline.

use crate::moduli::Moduli;
use crate::symbol::eps;
use fs_math::c64::C64;

/// Stored micropolar fields (u, φ) on a periodic cell-centred n³ grid of
/// box length `lbox` (spacing h = lbox/n). Data layout: `u[j][cell]`,
/// `phi[j][cell]` with cell = (ix·n + iy)·n + iz (row-major, fixed
/// reduction order).
pub struct MicroField {
    /// Cells per axis.
    pub n: usize,
    /// Grid spacing h = lbox/n.
    pub h: f64,
    /// Displacement components u_x, u_y, u_z.
    pub u: [Vec<f64>; 3],
    /// Micro-rotation components φ_x, φ_y, φ_z.
    pub phi: [Vec<f64>; 3],
}

impl MicroField {
    /// Zero fields on an n³ periodic grid of box length `lbox`.
    #[must_use]
    pub fn zeros(n: usize, lbox: f64) -> MicroField {
        let sz = n * n * n;
        MicroField {
            n,
            h: lbox / n as f64,
            u: [vec![0.0; sz], vec![0.0; sz], vec![0.0; sz]],
            phi: [vec![0.0; sz], vec![0.0; sz], vec![0.0; sz]],
        }
    }

    /// Box length L = n·h.
    #[must_use]
    pub fn lbox(&self) -> f64 {
        self.h * self.n as f64
    }

    /// Add one Fourier mode a·cos(k·x + p) to component `comp`
    /// (0–2 = u_x..u_z, 3–5 = φ_x..φ_z), k = 2π m/L. Used to build random
    /// SMOOTH periodic fields deterministically.
    pub fn add_mode(&mut self, comp: usize, m: [i64; 3], amp: f64, phase: f64) {
        let n = self.n as i64;
        let tau = 2.0 * std::f64::consts::PI;
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let dot = (m[0] * ix + m[1] * iy + m[2] * iz).rem_euclid(n);
                    let val = amp * (tau * dot as f64 / n as f64 + phase).cos();
                    let cell = ((ix * n + iy) * n + iz) as usize;
                    if comp < 3 {
                        self.u[comp][cell] += val;
                    } else {
                        self.phi[comp - 3][cell] += val;
                    }
                }
            }
        }
    }

    /// Overwrite the fields with the real plane wave
    /// (u_j, φ_j)(x) = Re(v̂_j e^{ik·x}), k = 2π m/L, with v̂ the complex
    /// 6-amplitude in the crate's (u_x, u_y, u_z, φ_x, φ_y, φ_z) order.
    pub fn set_plane_wave(&mut self, m: [i64; 3], v: &[C64; 6]) {
        let n = self.n as i64;
        let tau = 2.0 * std::f64::consts::PI;
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let dot = (m[0] * ix + m[1] * iy + m[2] * iz).rem_euclid(n);
                    let ang = tau * dot as f64 / n as f64;
                    let (s, c) = (ang.sin(), ang.cos());
                    let cell = ((ix * n + iy) * n + iz) as usize;
                    for j in 0..3 {
                        // Re(v e^{iθ}) = v.re cosθ − v.im sinθ
                        self.u[j][cell] = v[j].re * c - v[j].im * s;
                        self.phi[j][cell] = v[3 + j].re * c - v[3 + j].im * s;
                    }
                }
            }
        }
    }

    /// The exact symbol wavevector of the central-difference operator for
    /// integer mode m: k_eff_i = sin(2π m_i/n)/h (continuum k = 2π m/L).
    #[must_use]
    pub fn k_eff(&self, m: [i64; 3]) -> [f64; 3] {
        let n = self.n as f64;
        let tau = 2.0 * std::f64::consts::PI;
        [
            (tau * m[0] as f64 / n).sin() / self.h,
            (tau * m[1] as f64 / n).sin() / self.h,
            (tau * m[2] as f64 / n).sin() / self.h,
        ]
    }
}

/// Neighbour cell indices ±ê_d of `cell` on the periodic grid (n³).
#[inline]
fn wrap_pair(n: usize, ix: usize, iy: usize, iz: usize, d: usize) -> (usize, usize) {
    let (mut px, mut py, mut pz) = (ix, iy, iz);
    let (mut mx, mut my, mut mz) = (ix, iy, iz);
    match d {
        0 => {
            px = if ix + 1 == n { 0 } else { ix + 1 };
            mx = if ix == 0 { n - 1 } else { ix - 1 };
        }
        1 => {
            py = if iy + 1 == n { 0 } else { iy + 1 };
            my = if iy == 0 { n - 1 } else { iy - 1 };
        }
        _ => {
            pz = if iz + 1 == n { 0 } else { iz + 1 };
            mz = if iz == 0 { n - 1 } else { iz - 1 };
        }
    }
    ((px * n + py) * n + pz, (mx * n + my) * n + mz)
}

/// The nine (e_ij, Γ_ij) pairs at one cell (row = 3i + j), central
/// differences, periodic wrap.
#[inline]
fn tensors_at(f: &MicroField, ix: usize, iy: usize, iz: usize) -> ([f64; 9], [f64; 9]) {
    let n = f.n;
    let cell = (ix * n + iy) * n + iz;
    let inv2h = 1.0 / (2.0 * f.h);
    let mut e = [0.0f64; 9];
    let mut g = [0.0f64; 9];
    for i in 0..3 {
        let (p, m) = wrap_pair(n, ix, iy, iz, i);
        for j in 0..3 {
            let r = 3 * i + j;
            e[r] = (f.u[j][p] - f.u[j][m]) * inv2h;
            g[r] = (f.phi[j][p] - f.phi[j][m]) * inv2h;
            for l in 0..3 {
                let s = eps(i, j, l);
                if s != 0.0 {
                    e[r] -= s * f.phi[l][cell];
                }
            }
        }
    }
    (e, g)
}

/// Split a rank-2 value into (trace, sym[9], asym[9]); `dev` subtracts
/// tr/3 from the symmetric diagonal (the deviatoric χ₂ convention).
#[inline]
fn parts(t: &[f64; 9], dev: bool) -> (f64, [f64; 9], [f64; 9]) {
    let tr = t[0] + t[4] + t[8];
    let mut sym = [0.0f64; 9];
    let mut asym = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            let rij = 3 * i + j;
            let rji = 3 * j + i;
            sym[rij] = 0.5 * (t[rij] + t[rji]);
            asym[rij] = 0.5 * (t[rij] - t[rji]);
        }
    }
    if dev {
        for diag in [0usize, 4, 8] {
            sym[diag] -= tr / 3.0;
        }
    }
    (tr, sym, asym)
}

/// W_χ density at every cell plus the total energy E = h³ Σ w (fixed
/// cell-major summation order). Returns (density, total). All-zero χ is a
/// hard guard: returns exact zeros.
#[must_use]
pub fn wchi_energy(f: &MicroField, m: &Moduli) -> (Vec<f64>, f64) {
    let sz = f.n * f.n * f.n;
    if m.chi1 == 0.0 && m.chi2 == 0.0 && m.chi3 == 0.0 {
        return (vec![0.0; sz], 0.0);
    }
    let mut dens = vec![0.0f64; sz];
    let mut total = 0.0f64;
    let h3 = f.h * f.h * f.h;
    for ix in 0..f.n {
        for iy in 0..f.n {
            for iz in 0..f.n {
                let (e, g) = tensors_at(f, ix, iy, iz);
                let (etr, esym, easym) = parts(&e, m.chi_deviatoric);
                let (gtr, gsym, gasym) = parts(&g, m.chi_deviatoric);
                let mut w = m.chi1 * etr * gtr;
                for r in 0..9 {
                    w += m.chi2 * esym[r] * gsym[r] + m.chi3 * easym[r] * gasym[r];
                }
                dens[(ix * f.n + iy) * f.n + iz] = w;
                total += w * h3;
            }
        }
    }
    (dens, total)
}

/// Analytic gradient of the total W_χ energy with respect to the stored
/// fields: (∂E/∂u_j, ∂E/∂φ_j) per cell. Two-pass flux form:
///
///   P^e_ij = ∂w/∂e_ij = χ₁ δ_ij Γ_kk + χ₂ [dev]Γ_(ij) + χ₃ Γ_[ij]
///   P^Γ_ij = ∂w/∂Γ_ij = χ₁ δ_ij e_kk + χ₂ [dev]e_(ij) + χ₃ e_[ij]
///   ∂E/∂u_j = −h³ Σ_i D_i P^e_ij            (D self-adjoint-negative)
///   ∂E/∂φ_l = −h³ (Σ_ij ε_ijl P^e_ij + Σ_i D_i P^Γ_il)
///
/// (sym/asym/dev projectors are self-adjoint, so the projected Γ/e parts
/// ARE the derivatives). All-zero χ returns exact zeros.
#[must_use]
pub fn wchi_grad(f: &MicroField, m: &Moduli) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
    let sz = f.n * f.n * f.n;
    let zero3 = || [vec![0.0f64; sz], vec![0.0f64; sz], vec![0.0f64; sz]];
    let (mut gu, mut gphi) = (zero3(), zero3());
    if m.chi1 == 0.0 && m.chi2 == 0.0 && m.chi3 == 0.0 {
        return (gu, gphi);
    }
    // Pass 1: flux fields P^e, P^Γ (9 components each, row = 3i + j).
    let mut pe: Vec<[f64; 9]> = vec![[0.0; 9]; sz];
    let mut pg: Vec<[f64; 9]> = vec![[0.0; 9]; sz];
    for ix in 0..f.n {
        for iy in 0..f.n {
            for iz in 0..f.n {
                let (e, g) = tensors_at(f, ix, iy, iz);
                let (etr, esym, easym) = parts(&e, m.chi_deviatoric);
                let (gtr, gsym, gasym) = parts(&g, m.chi_deviatoric);
                let cell = (ix * f.n + iy) * f.n + iz;
                for i in 0..3 {
                    for j in 0..3 {
                        let r = 3 * i + j;
                        let dij = if i == j { 1.0 } else { 0.0 };
                        pe[cell][r] = m.chi1 * dij * gtr + m.chi2 * gsym[r] + m.chi3 * gasym[r];
                        pg[cell][r] = m.chi1 * dij * etr + m.chi2 * esym[r] + m.chi3 * easym[r];
                    }
                }
            }
        }
    }
    // Pass 2: adjoints. Central difference is anti-self-adjoint on the
    // periodic grid, so D† = −D and ∂E/∂field = −h³ D_i(P_i·).
    let h3 = f.h * f.h * f.h;
    let inv2h = 1.0 / (2.0 * f.h);
    for ix in 0..f.n {
        for iy in 0..f.n {
            for iz in 0..f.n {
                let cell = (ix * f.n + iy) * f.n + iz;
                for j in 0..3 {
                    let mut du = 0.0f64;
                    let mut dphi = 0.0f64;
                    for i in 0..3 {
                        let (p, mm) = wrap_pair(f.n, ix, iy, iz, i);
                        let r = 3 * i + j;
                        du -= (pe[p][r] - pe[mm][r]) * inv2h;
                        dphi -= (pg[p][r] - pg[mm][r]) * inv2h;
                    }
                    // pointwise −ε_ikl term of e: ∂e_ik/∂φ_j = −ε_ikj
                    for i in 0..3 {
                        for kx in 0..3 {
                            let s = eps(i, kx, j);
                            if s != 0.0 {
                                dphi -= s * pe[cell][3 * i + kx];
                            }
                        }
                    }
                    gu[j][cell] = h3 * du;
                    gphi[j][cell] = h3 * dphi;
                }
            }
        }
    }
    (gu, gphi)
}

/// Inner product Σ (gu·du + gphi·dphi) — fixed order, for the FD gate.
#[must_use]
pub fn grad_dot(
    gu: &[Vec<f64>; 3],
    gphi: &[Vec<f64>; 3],
    du: &MicroField,
) -> f64 {
    let mut acc = 0.0f64;
    for j in 0..3 {
        for c in 0..gu[j].len() {
            acc += gu[j][c] * du.u[j][c];
        }
    }
    for j in 0..3 {
        for c in 0..gphi[j].len() {
            acc += gphi[j][c] * du.phi[j][c];
        }
    }
    acc
}

/// E(f + t·d) for the central-difference FD gate (fields displaced along
/// the direction field d; f and d must share (n, h)).
#[must_use]
pub fn energy_displaced(f: &MicroField, d: &MicroField, t: f64, m: &Moduli) -> f64 {
    let mut shifted = MicroField::zeros(f.n, f.lbox());
    for j in 0..3 {
        for c in 0..f.u[j].len() {
            shifted.u[j][c] = f.u[j][c] + t * d.u[j][c];
            shifted.phi[j][c] = f.phi[j][c] + t * d.phi[j][c];
        }
    }
    wchi_energy(&shifted, m).1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::symbol::stiffness;

    fn chi_only(chi1: f64, chi2: f64, chi3: f64, dev: bool) -> Moduli {
        let mut m = Moduli::bench();
        m.lam = 0.0;
        m.mu = 0.0;
        m.mu_c = 0.0;
        m.alpha = 0.0;
        m.beta = 0.0;
        m.gamma = 0.0;
        m.m_v = 0.0;
        m.chi1 = chi1;
        m.chi2 = chi2;
        m.chi3 = chi3;
        m.chi_deviatoric = dev;
        m
    }

    fn smooth(seed: u64, n: usize) -> MicroField {
        let mut f = MicroField::zeros(n, 2.0 * std::f64::consts::PI);
        let mut s = seed;
        let mut next = move || {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            ((s >> 11) as f64) / ((1u64 << 53) as f64)
        };
        for comp in 0..6 {
            for _ in 0..4 {
                let m = [
                    (next() * 7.0) as i64 - 3,
                    (next() * 7.0) as i64 - 3,
                    (next() * 7.0) as i64 - 3,
                ];
                f.add_mode(comp, m, 2.0 * next() - 1.0, 2.0 * std::f64::consts::PI * next());
            }
        }
        f
    }

    #[test]
    fn zero_chi_is_exactly_zero() {
        let f = smooth(7, 6);
        let m = Moduli::bench(); // all χ default 0
        let (dens, tot) = wchi_energy(&f, &m);
        assert!(tot == 0.0 && dens.iter().all(|&d| d == 0.0));
        let (gu, gphi) = wchi_grad(&f, &m);
        assert!(gu.iter().chain(gphi.iter()).all(|v| v.iter().all(|&x| x == 0.0)));
    }

    #[test]
    fn fd_matches_analytic_gradient() {
        let f = smooth(11, 6);
        let d = smooth(23, 6);
        for m in [chi_only(0.4, -0.3, 0.25, false), chi_only(0.4, -0.3, 0.25, true)] {
            let (gu, gphi) = wchi_grad(&f, &m);
            let dot = grad_dot(&gu, &gphi, &d);
            let t = 1e-3; // E quadratic ⇒ central FD exact to round-off
            let fd = (energy_displaced(&f, &d, t, &m) - energy_displaced(&f, &d, -t, &m))
                / (2.0 * t);
            let rel = (fd - dot).abs() / dot.abs().max(fd.abs()).max(1e-30);
            assert!(rel < 1e-9, "FD {fd} vs analytic {dot} (rel {rel})");
        }
    }

    #[test]
    fn plane_wave_matches_symbol() {
        let mut f = MicroField::zeros(8, 2.0 * std::f64::consts::PI);
        let v = [
            C64::new(0.3, -0.1),
            C64::new(-0.7, 0.2),
            C64::new(0.1, 0.5),
            C64::new(0.4, 0.4),
            C64::new(-0.2, 0.6),
            C64::new(0.9, -0.3),
        ];
        let mode = [1i64, -2, 3];
        f.set_plane_wave(mode, &v);
        let m = chi_only(0.7, -0.4, 0.3, false);
        let (_, e_real) = wchi_energy(&f, &m);
        let keff = f.k_eff(mode);
        let kmat = stiffness(keff, &m);
        let vol = f.lbox().powi(3);
        let mut quad = 0.0f64;
        let mut scale = 0.0f64;
        for a in 0..6 {
            for b in 0..6 {
                let z = v[a].conj() * kmat[a * 6 + b] * v[b];
                quad += z.re;
                scale += v[a].abs() * kmat[a * 6 + b].abs() * v[b].abs();
            }
        }
        let e_sym = 0.25 * vol * quad;
        let dev = (e_real - e_sym).abs() / (0.25 * vol * scale).max(1e-30);
        assert!(dev < 1e-13, "real {e_real} vs symbol {e_sym} (scaled dev {dev})");
    }
}
