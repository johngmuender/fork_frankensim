//! Fail-closed certification wrap (the fs-cosserat-pilot pattern made
//! reusable, domain "gum-core:topo:v1"): Check -> Certified<f64> ->
//! EvidencePackage claims -> BLAKE3 Merkle root, verified by fs-checker
//! in three modes (deny-all expect-fail, source-certificate capability
//! expect-pass, tampered-root expect-fail). Only deterministic math and
//! fixed loop orders feed the hashes; the gates binary runs the whole
//! pipeline twice and asserts a bit-identical root.

use fs_blake3::ContentHash;
use fs_evidence::{Evidence, ProvenanceHash};
use fs_ivl::Interval;
use fs_package::origin::{
    SourceCertificateRequest, SourceCertificateVerifier, VerificationCapabilities,
    VerificationDecision,
};
use fs_package::{Claim, EvidencePackage, Provenance};

/// The certification domain of this crate's claims.
pub const DOMAIN: &str = "gum-core:topo:v1";

/// One gate/referee check: a statement, its measured enclosure, the
/// target, the verdict, and a human note carrying the raw numbers.
pub struct Check {
    /// Short stable id (feeds the certificate hash; no '|').
    pub id: &'static str,
    /// The claim statement (no '|').
    pub statement: String,
    /// Measured value or enclosure.
    pub enclosure: Interval,
    /// The referee target.
    pub target: f64,
    /// Verdict.
    pub pass: bool,
    /// Raw numbers for the reader.
    pub note: String,
}

/// Convenience constructor.
#[must_use]
pub fn check(
    id: &'static str,
    statement: &str,
    enclosure: Interval,
    target: f64,
    pass: bool,
    note: String,
) -> Check {
    Check { id, statement: statement.to_string(), enclosure, target, pass, note }
}

/// Canonical certificate content for a claim: id | statement | exact
/// bit patterns of the interval endpoints, hashed in the crate domain.
fn certificate_hash(id: &str, statement: &str, lo: f64, hi: f64) -> ContentHash {
    let canon = format!("{id}|{statement}|{:016x}|{:016x}", lo.to_bits(), hi.to_bits());
    fs_blake3::hash_domain(DOMAIN, canon.as_bytes())
}

/// Recompute-and-compare certificate policy: accept a source-certificate
/// claim iff the declared artifact hash equals the recomputation from
/// the claim's own typed content.
struct TopoCertVerifier;

impl SourceCertificateVerifier for TopoCertVerifier {
    fn verify(&self, req: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let policy = fs_blake3::hash_domain(DOMAIN, b"policy:recompute-and-compare");
        let expect = certificate_hash(req.claim_id, req.statement, req.lo, req.hi);
        if req.producer == "fs-gum-topo" && req.certificate_hash == expect {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

/// Outcome of the certification pass.
pub struct CertOutcome {
    /// Claims that survived the fail-closed Certified<f64> gate.
    pub certified: usize,
    /// Total checks offered (failing checks are excluded by
    /// construction).
    pub total: usize,
    /// The package's BLAKE3 Merkle root (hex).
    pub root_hex: String,
    /// Mode 1, deny-all: MUST be false (anti-laundering).
    pub mode_deny: bool,
    /// Mode 2, certificate capability: MUST be true.
    pub mode_cert: bool,
    /// Mode 3, tampered root: MUST be false.
    pub mode_tamper: bool,
}

/// Bundle the passing checks into an EvidencePackage and verify it in
/// the three checker modes.
///
/// # Panics
/// If the package Merkle root cannot be assembled (empty package).
#[must_use]
pub fn certify(checks: &[Check]) -> CertOutcome {
    let mut pkg = EvidencePackage::new(Provenance::new(
        "fs-gum-topo-v0.1.0-gum-core-phase-b2-topological-diagnostics",
        "frankensim-fs-fft-fs-ivl-fs-evidence-fs-package-fs-checker-std-only",
    ));
    let mut certified = 0usize;
    for c in checks {
        if !c.pass || !c.enclosure.lo().is_finite() || !c.enclosure.hi().is_finite() {
            continue;
        }
        let prov = ProvenanceHash::of_bytes(
            format!(
                "{}|{}|{:016x}|{:016x}",
                c.id,
                c.statement,
                c.enclosure.lo().to_bits(),
                c.enclosure.hi().to_bits()
            )
            .as_bytes(),
        );
        let ev =
            Evidence::enclosed(c.enclosure.midpoint(), c.enclosure.lo(), c.enclosure.hi(), prov);
        if ev.certified().is_ok() {
            certified += 1;
            pkg = pkg.with_claim(Claim::from_certificate(
                c.id,
                &c.statement,
                c.enclosure.lo(),
                c.enclosure.hi(),
                "fs-gum-topo",
                certificate_hash(c.id, &c.statement, c.enclosure.lo(), c.enclosure.hi()).to_hex(),
            ));
        }
    }
    let root = pkg.try_merkle_root().expect("package root");
    let root_hex = root.to_hex();

    // Mode 1 — deny-all: Verified-origin claims must NOT be admitted.
    let mode_deny = fs_checker::check_against_root(&pkg, root).passed();

    // Mode 2 — recompute-and-compare certificate capability.
    let verifier = TopoCertVerifier;
    let caps = VerificationCapabilities::deny_all().with_source_certificates(&verifier);
    let mode_cert = fs_checker::check_with_capabilities(&pkg, Some(root), None, &caps).passed();

    // Mode 3 — a tampered expected root must fail.
    let mut wrong = root;
    wrong.0[0] ^= 0xff;
    let mode_tamper =
        fs_checker::check_with_capabilities(&pkg, Some(wrong), None, &caps).passed();

    CertOutcome {
        certified,
        total: checks.len(),
        root_hex,
        mode_deny,
        mode_cert,
        mode_tamper,
    }
}
