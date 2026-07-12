//! The Dzyaloshinskii soft-sector condensate functional (Phase E2 — the
//! E-H1 "mis-sited chiral condensate" audit machinery, corpus II.H).
//!
//! A conical texture of tilt θ and pitch q in a sector with twist
//! stiffness γ_s, chiral gain χ_s, and gap Δ_s. The corpus
//! (`01-GUM-Omega-Paper-v2.0.1.md` II.H(i); Substrate Course eq (6.1))
//! prints the density
//!
//!   f(θ, q) = sin²θ (−χ_s q + ½ γ_s q²) + ½ Δ_s² θ²          (II.H(i))
//!
//! and derives: pitch q\* = χ_s/γ_s (θ-independent), chiral gain at q\* =
//! −χ_s²/(2γ_s), and — from the SMALL-θ expansion — the Dzyaloshinskii
//! criterion: condensation iff χ_s² > γ_s Δ_s², margin
//! 𝔪 ≡ χ_s²/(γ_s Δ_s²), threshold exactly at 𝔪 = 1.
//!
//! # What the text PINS (and where)
//!
//! - The density above, TO QUADRATIC ORDER in θ (II.H(i) / Course 6.1 —
//!   the text only ever uses it in the small-θ expansion).
//! - q\* = χ_s/γ_s and the criterion χ_s² > γ_s Δ_s² (II.H(i)).
//! - Theorem H2-1 (ε-blindness): 𝔪 = χ̄²/(γ̄Δ̄²) is a pure number;
//!   ⟨r9⟩ pins **𝔪 = 1.9 ± 0.4** — the ONLY text-pinned numeric referee
//!   of the chiral sector (survey, Skyrme/chiral entry).
//! - Theorem H-3 (the blue fog): the condensate tilt
//!   **sinθ_c = √(1 − 1/𝔪)** = 0.69 ± 0.11 at 𝔪 = 1.9 (also App. K.5;
//!   Tier-5a check C7a validates the arithmetic: 0.688 → 0.69).
//!
//! # What the text leaves FREE — and the one finding this module records
//!
//! χ_s, γ_s, Δ_s individually are unpinned (only the combinations q\* and
//! 𝔪 are). And the FINITE-θ completion of the gap potential is unpinned:
//! the printed ½Δ_s²θ² taken literally to all orders is INCONSISTENT with
//! the pinned tilt — its stationarity condition is 𝔪 sinθ cosθ = θ, giving
//! sinθ_c ≈ 0.7924 at 𝔪 = 1.9, not 0.688. Requiring (a) the printed
//! chiral term −(χ_s²/2γ_s) sin²θ and (b) the pinned tilt
//! sin²θ_c = 1 − 1/𝔪 for EVERY 𝔪 > 1 forces the gap potential's
//! derivative to be G′(θ) = Δ_s² tanθ, i.e. determines it UNIQUELY:
//!
//!   G(θ) = −Δ_s² ln cosθ = ½Δ_s²θ² + Δ_s²θ⁴/12 + …
//!
//! whose quadratic truncation is exactly the printed gap term. This module
//! therefore implements BOTH readings: [`SoftSector::f_tilt_quadgap`] (the
//! literal printed form — correct criterion and threshold, tilt 0.7924 at
//! 𝔪 = 1.9) and [`SoftSector::f_tilt`] (the −Δ² ln cosθ completion — same
//! criterion, same threshold, and the Thm H-3 tilt exactly). The tilt
//! gates (G24/G25) run on the completed form; the quadratic-gap reading's
//! threshold is co-gated at G26.
//!
//! Everything here is deterministic (fixed-iteration golden section +
//! Newton polish; no RNG).

/// Soft-sector parameters (χ_s, γ_s, Δ_s) — individually FREE; only
/// 𝔪 = χ²/(γΔ²) and q\* = χ/γ are text-pinned combinations.
#[derive(Debug, Clone, Copy)]
pub struct SoftSector {
    /// Chiral gain χ_s (> 0).
    pub chi_s: f64,
    /// Twist stiffness γ_s (> 0).
    pub gamma_s: f64,
    /// Gap Δ_s (> 0).
    pub delta_s: f64,
}

impl SoftSector {
    /// Parameters realizing a given margin 𝔪 at chosen (γ_s, Δ_s):
    /// χ_s = √(𝔪 γ_s Δ_s²).
    #[must_use]
    pub fn with_margin(margin: f64, gamma_s: f64, delta_s: f64) -> SoftSector {
        SoftSector { chi_s: (margin * gamma_s * delta_s * delta_s).sqrt(), gamma_s, delta_s }
    }

    /// The condensation margin 𝔪 = χ_s²/(γ_s Δ_s²) (Thm H2-1; ⟨r9⟩ pins
    /// 𝔪 = 1.9 ± 0.4).
    #[must_use]
    pub fn margin(&self) -> f64 {
        self.chi_s * self.chi_s / (self.gamma_s * self.delta_s * self.delta_s)
    }

    /// Optimal pitch q\* = χ_s/γ_s (II.H(i); θ-independent).
    #[must_use]
    pub fn pitch_star(&self) -> f64 {
        self.chi_s / self.gamma_s
    }

    /// Dzyaloshinskii criterion: condensation iff χ_s² > γ_s Δ_s²
    /// (𝔪 > 1).
    #[must_use]
    pub fn condenses(&self) -> bool {
        self.margin() > 1.0
    }

    /// The density EXACTLY as printed (corpus II.H(i) / Course eq (6.1)):
    /// f(θ, q) = sin²θ(−χ_s q + ½γ_s q²) + ½Δ_s²θ².
    #[must_use]
    pub fn f_printed(&self, theta: f64, q: f64) -> f64 {
        let s2 = theta.sin() * theta.sin();
        s2 * (-self.chi_s * q + 0.5 * self.gamma_s * q * q)
            + 0.5 * self.delta_s * self.delta_s * theta * theta
    }

    /// q-optimized profile with the LITERAL quadratic gap:
    /// f(θ) = −(χ_s²/2γ_s) sin²θ + ½Δ_s²θ². Threshold at 𝔪 = 1; its
    /// stationarity is 𝔪 sinθcosθ = θ (NOT the Thm H-3 tilt — see the
    /// module docs).
    #[must_use]
    pub fn f_tilt_quadgap(&self, theta: f64) -> f64 {
        let s = theta.sin();
        -0.5 * self.chi_s * self.chi_s / self.gamma_s * s * s
            + 0.5 * self.delta_s * self.delta_s * theta * theta
    }

    /// q-optimized profile with the uniquely-determined saturating gap
    /// G(θ) = −Δ_s² ln cosθ (= ½Δ_s²θ² + O(θ⁴)):
    /// f(θ) = −(χ_s²/2γ_s) sin²θ − Δ_s² ln cosθ. Same criterion and
    /// threshold as the printed form; minimizer satisfies the pinned
    /// Thm H-3 tilt sin²θ_c = 1 − 1/𝔪 exactly, for every 𝔪 > 1.
    #[must_use]
    pub fn f_tilt(&self, theta: f64) -> f64 {
        let s = theta.sin();
        -0.5 * self.chi_s * self.chi_s / self.gamma_s * s * s
            - self.delta_s * self.delta_s * theta.cos().ln()
    }

    /// d f_tilt/dθ = sinθ (Δ_s²/cosθ − (χ_s²/γ_s) cosθ).
    #[must_use]
    pub fn df_tilt(&self, theta: f64) -> f64 {
        let (s, c) = (theta.sin(), theta.cos());
        s * (self.delta_s * self.delta_s / c - self.chi_s * self.chi_s / self.gamma_s * c)
    }

    /// d² f_tilt/dθ² = Δ_s²/cos²θ − (χ_s²/γ_s) cos 2θ; at θ = 0 this is
    /// Δ_s²(1 − 𝔪) — the sign flip AT 𝔪 = 1 is the threshold.
    #[must_use]
    pub fn d2f_tilt(&self, theta: f64) -> f64 {
        let c = theta.cos();
        self.delta_s * self.delta_s / (c * c)
            - self.chi_s * self.chi_s / self.gamma_s * (2.0 * theta).cos()
    }

    /// Analytic condensate tilt (Thm H-3): sinθ_c = √(1 − 1/𝔪) for
    /// 𝔪 > 1, else 0.
    #[must_use]
    pub fn sin_theta_c_analytic(&self) -> f64 {
        let m = self.margin();
        if m > 1.0 {
            (1.0 - 1.0 / m).sqrt()
        } else {
            0.0
        }
    }

    /// Numerical minimizer of [`SoftSector::f_tilt`] on [0, π/2):
    /// 120-iteration golden section (deterministic) plus up to 8 Newton
    /// polish steps on df/dθ (skipped when the golden minimum sits at the
    /// θ = 0 boundary, i.e. below threshold).
    #[must_use]
    pub fn theta_c_numeric(&self) -> f64 {
        let golden = 0.618_033_988_749_894_9_f64;
        let (mut a, mut b) = (0.0f64, std::f64::consts::FRAC_PI_2 * 0.999_999);
        let mut x1 = b - golden * (b - a);
        let mut x2 = a + golden * (b - a);
        let mut f1 = self.f_tilt(x1);
        let mut f2 = self.f_tilt(x2);
        for _ in 0..120 {
            if f1 <= f2 {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - golden * (b - a);
                f1 = self.f_tilt(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + golden * (b - a);
                f2 = self.f_tilt(x2);
            }
        }
        let mut th = 0.5 * (a + b);
        // Newton polish only when strictly interior (above threshold the
        // stationary point is unique and f'' > 0 there).
        if th > 1e-8 && self.margin() > 1.0 {
            for _ in 0..8 {
                let d2 = self.d2f_tilt(th);
                if d2 <= 0.0 {
                    break;
                }
                let step = self.df_tilt(th) / d2;
                th -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
        }
        th.max(0.0)
    }

    /// Numerical pitch: golden-section minimum of f_printed(θ_ref, q) over
    /// q ∈ [0, 4χ/γ] at fixed θ_ref (deterministic, 120 iterations).
    /// Referee: must equal q\* = χ_s/γ_s independent of θ_ref.
    #[must_use]
    pub fn pitch_numeric(&self, theta_ref: f64) -> f64 {
        let golden = 0.618_033_988_749_894_9_f64;
        let (mut a, mut b) = (0.0f64, 4.0 * self.chi_s / self.gamma_s);
        let mut x1 = b - golden * (b - a);
        let mut x2 = a + golden * (b - a);
        let mut f1 = self.f_printed(theta_ref, x1);
        let mut f2 = self.f_printed(theta_ref, x2);
        for _ in 0..120 {
            if f1 <= f2 {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - golden * (b - a);
                f1 = self.f_printed(theta_ref, x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + golden * (b - a);
                f2 = self.f_printed(theta_ref, x2);
            }
        }
        // Parabolic-vertex polish: f is EXACTLY quadratic in q, so the
        // three-point vertex formula is exact up to round-off (golden
        // section alone can only locate a flat quadratic minimum to
        // ~sqrt(eps) relative).
        let q0 = 0.5 * (a + b);
        let d = 1e-3 * self.chi_s / self.gamma_s;
        let (fp, f0, fm) = (
            self.f_printed(theta_ref, q0 + d),
            self.f_printed(theta_ref, q0),
            self.f_printed(theta_ref, q0 - d),
        );
        let denom = fp - 2.0 * f0 + fm;
        if denom > 0.0 {
            q0 - 0.5 * d * (fp - fm) / denom
        } else {
            q0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tilt_matches_thm_h3_at_r9_margin() {
        let s = SoftSector::with_margin(1.9, 0.7, 1.1);
        let th = s.theta_c_numeric();
        let want = (1.0 - 1.0 / 1.9_f64).sqrt(); // 0.6882472016116852
        assert!((th.sin() - want).abs() < 1e-13, "sin theta_c = {} vs {want}", th.sin());
        assert!(s.df_tilt(th).abs() < 1e-12);
        // pitch is theta-independent and equals chi/gamma
        for th_ref in [0.3, 0.7, 1.2] {
            let q = s.pitch_numeric(th_ref);
            assert!((q - s.pitch_star()).abs() < 1e-10 * s.pitch_star());
        }
    }

    #[test]
    fn threshold_at_margin_one() {
        for m in [0.3, 0.7, 0.999] {
            let s = SoftSector::with_margin(m, 1.3, 0.9);
            assert!(!s.condenses());
            assert!(s.theta_c_numeric() < 1e-6, "no tilt below threshold");
            assert!(s.d2f_tilt(0.0) > 0.0);
        }
        for m in [1.001, 1.1, 2.5] {
            let s = SoftSector::with_margin(m, 1.3, 0.9);
            assert!(s.condenses());
            let want = (1.0 - 1.0 / m).sqrt();
            assert!((s.theta_c_numeric().sin() - want).abs() < 1e-9);
            assert!(s.d2f_tilt(0.0) < 0.0);
        }
    }

    #[test]
    fn quadgap_reading_same_threshold_different_tilt() {
        // The literal printed gap: same threshold, but the stationarity is
        // m sin cos = theta, whose tilt at m = 1.9 is 0.7924, not 0.688.
        let s = SoftSector::with_margin(1.9, 1.0, 1.0);
        // stationarity residual of the quadgap form at the H-3 tilt is NOT 0
        let th_h3 = (1.0 - 1.0 / 1.9_f64).sqrt().asin();
        let m = 1.9;
        let resid = m * th_h3.sin() * th_h3.cos() - th_h3;
        assert!(resid.abs() > 1e-2, "quadgap stationarity should reject the H-3 tilt");
        // and below threshold quadgap is minimized at 0
        let s0 = SoftSector::with_margin(0.9, 1.0, 1.0);
        for th in [0.1, 0.5, 1.0, 1.4] {
            assert!(s0.f_tilt_quadgap(th) > 0.0);
        }
        assert!(s.f_tilt_quadgap(0.5) < 0.0, "above threshold tilting must pay");
    }
}
