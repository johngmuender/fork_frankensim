//! Micropolar moduli and the Tier-1 benchmark point.

/// Which mass term closes the quadratic action.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MassCase {
    /// Case A: objective mass (m_V²/2)|ψ|², ψ = φ − ½ curl u.
    A,
    /// Case B: M-1 defect (m_V²/2)|φ|².
    B,
}

/// The full modulus set of the linearized micropolar continuum
/// (spectrum.py `BENCH` naming, dimensionless).
#[derive(Debug, Clone, Copy)]
pub struct Moduli {
    /// Mass density (translational inertia).
    pub rho0: f64,
    /// Rotational micro-inertia.
    pub j: f64,
    /// Lamé λ.
    pub lam: f64,
    /// Shear modulus μ.
    pub mu: f64,
    /// Cosserat couple modulus μ_c.
    pub mu_c: f64,
    /// Curvature (wryness) trace modulus α.
    pub alpha: f64,
    /// Curvature symmetric modulus β.
    pub beta: f64,
    /// Curvature antisymmetric modulus γ.
    pub gamma: f64,
    /// Mass parameter m_V.
    pub m_v: f64,
    /// Chiral coupling χ₃ on Re(e_[ij]* Γ_[ij]).
    pub chi3: f64,
    /// Mass-term case (A objective / B defect).
    pub case: MassCase,
}

impl Moduli {
    /// The Tier-1 benchmark: ρ₀ = J = λ = μ = 1, μ_c = 5,
    /// α = β = γ = 0.5, m_V = 1, χ₃ = 0, Case A.
    #[must_use]
    pub fn bench() -> Moduli {
        Moduli {
            rho0: 1.0,
            j: 1.0,
            lam: 1.0,
            mu: 1.0,
            mu_c: 5.0,
            alpha: 0.5,
            beta: 0.5,
            gamma: 0.5,
            m_v: 1.0,
            chi3: 0.0,
            case: MassCase::A,
        }
    }
}
