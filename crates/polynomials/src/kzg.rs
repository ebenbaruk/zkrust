use crate::DensePolynomial;
use zkrust_curves::{ate_pairing, msm_g1, G1Affine, G2Affine};
use zkrust_fields::{FieldElement, Fr};

/// KZG structured reference string (SRS) / public parameters.
///
/// Generated from a trusted setup ceremony with secret tau.
/// - `powers_g1[i]` = [tau^i]_1 for i = 0..max_degree
/// - `g2` = [1]_2 (G2 generator)
/// - `tau_g2` = [tau]_2
#[derive(Clone, Debug)]
pub struct KzgParams {
    pub powers_g1: Vec<G1Affine>,
    pub g2: G2Affine,
    pub tau_g2: G2Affine,
}

/// A KZG polynomial commitment (a single G1 point).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KzgCommitment(pub G1Affine);

/// A KZG opening proof (a single G1 point).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KzgProof(pub G1Affine);

impl KzgParams {
    /// Generate SRS for testing with a known tau (INSECURE — for tests only).
    pub fn setup(max_degree: usize, tau: Fr) -> Self {
        let g1 = G1Affine::generator();
        let g2 = G2Affine::generator();

        let mut powers_g1 = Vec::with_capacity(max_degree + 1);
        let mut tau_pow = Fr::ONE;
        for _ in 0..=max_degree {
            powers_g1.push((g1 * tau_pow).to_affine());
            tau_pow *= tau;
        }

        let tau_g2 = (g2 * tau).to_affine();

        Self {
            powers_g1,
            g2,
            tau_g2,
        }
    }

    /// Maximum polynomial degree this SRS supports.
    pub fn max_degree(&self) -> usize {
        self.powers_g1.len() - 1
    }
}

/// Commit to a polynomial: C = sum(coeffs[i] * [tau^i]_1).
pub fn kzg_commit(params: &KzgParams, poly: &DensePolynomial) -> KzgCommitment {
    if poly.is_zero() {
        return KzgCommitment(G1Affine::identity());
    }
    assert!(
        poly.coeffs.len() <= params.powers_g1.len(),
        "polynomial degree {} exceeds SRS max degree {}",
        poly.coeffs.len() - 1,
        params.max_degree()
    );

    let commitment = msm_g1(&params.powers_g1[..poly.coeffs.len()], &poly.coeffs);
    KzgCommitment(commitment.to_affine())
}

/// Create an opening proof for polynomial p at point z.
///
/// The quotient polynomial is q(x) = (p(x) - p(z)) / (x - z).
/// The proof is pi = commit(q).
pub fn kzg_open(params: &KzgParams, poly: &DensePolynomial, z: Fr) -> (Fr, KzgProof) {
    let value = poly.evaluate(&z);

    // q(x) = (p(x) - value) / (x - z)
    let numerator = poly.clone() - DensePolynomial::constant(value);
    let divisor = DensePolynomial::new(vec![-z, Fr::ONE]); // (x - z)
    let quotient = numerator.div_exact(&divisor);

    let proof = kzg_commit(params, &quotient);
    (value, KzgProof(proof.0))
}

/// Verify a KZG opening proof.
///
/// Check: e(C - [v]_1, [1]_2) == e(pi, [tau]_2 - [z]_2)
///
/// Equivalently: e(C - [v]_1, [1]_2) * e(-pi, [tau]_2 - [z]_2) == 1
pub fn kzg_verify(
    params: &KzgParams,
    commitment: &KzgCommitment,
    z: Fr,
    value: Fr,
    proof: &KzgProof,
) -> bool {
    let g1 = G1Affine::generator();

    // C - [v]_1
    let lhs_point = (commitment.0.to_projective() - (g1 * value)).to_affine();

    // [tau]_2 - [z]_2
    let rhs_g2 = (params.tau_g2.to_projective() - (params.g2 * z)).to_affine();

    // Check: e(lhs_point, [1]_2) == e(pi, rhs_g2)
    // i.e., e(lhs_point, [1]_2) * e(-pi, rhs_g2) == 1
    let e1 = ate_pairing(&lhs_point, &params.g2);
    let e2 = ate_pairing(&proof.0, &rhs_g2);

    e1 == e2
}

/// Batch verify multiple KZG openings at different points.
/// Uses random linear combination for efficiency.
pub fn kzg_batch_verify(
    params: &KzgParams,
    commitments: &[KzgCommitment],
    points: &[Fr],
    values: &[Fr],
    proofs: &[KzgProof],
) -> bool {
    assert_eq!(commitments.len(), points.len());
    assert_eq!(commitments.len(), values.len());
    assert_eq!(commitments.len(), proofs.len());

    // Simple: verify each individually
    // TODO: random linear combination batch verification
    for i in 0..commitments.len() {
        if !kzg_verify(params, &commitments[i], points[i], values[i], &proofs[i]) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_params() -> KzgParams {
        let tau = Fr::from(42u64);
        KzgParams::setup(32, tau)
    }

    #[test]
    fn test_commit_zero() {
        let params = test_params();
        let c = kzg_commit(&params, &DensePolynomial::zero());
        assert_eq!(c.0, G1Affine::identity());
    }

    #[test]
    fn test_commit_constant() {
        let params = test_params();
        let p = DensePolynomial::constant(Fr::from(7u64));
        let c = kzg_commit(&params, &p);
        let expected = (G1Affine::generator() * Fr::from(7u64)).to_affine();
        assert_eq!(c.0, expected);
    }

    #[test]
    fn test_kzg_open_verify() {
        let params = test_params();
        let poly = DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]); // 1 + 2x + 3x^2

        let commitment = kzg_commit(&params, &poly);
        let z = Fr::from(5u64);
        let (value, proof) = kzg_open(&params, &poly, z);

        // p(5) = 1 + 10 + 75 = 86
        assert_eq!(value, Fr::from(86u64));
        assert!(kzg_verify(&params, &commitment, z, value, &proof));
    }

    #[test]
    fn test_kzg_wrong_value_rejected() {
        let params = test_params();
        let poly = DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]);
        let commitment = kzg_commit(&params, &poly);
        let z = Fr::from(3u64);
        let (_value, proof) = kzg_open(&params, &poly, z);

        // Try verifying with wrong value
        let wrong_value = Fr::from(999u64);
        assert!(!kzg_verify(&params, &commitment, z, wrong_value, &proof));
    }

    #[test]
    fn test_kzg_wrong_point_rejected() {
        let params = test_params();
        let poly = DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]);
        let commitment = kzg_commit(&params, &poly);
        let z = Fr::from(5u64);
        let (value, proof) = kzg_open(&params, &poly, z);

        // Try verifying at a different point
        let wrong_z = Fr::from(6u64);
        assert!(!kzg_verify(&params, &commitment, wrong_z, value, &proof));
    }

    #[test]
    fn test_kzg_random_polynomial() {
        let params = test_params();
        let mut rng = rand::thread_rng();

        let coeffs: Vec<Fr> = (0..10).map(|_| Fr::random(&mut rng)).collect();
        let poly = DensePolynomial::new(coeffs);
        let commitment = kzg_commit(&params, &poly);

        let z = Fr::random(&mut rng);
        let (value, proof) = kzg_open(&params, &poly, z);

        assert_eq!(value, poly.evaluate(&z));
        assert!(kzg_verify(&params, &commitment, z, value, &proof));
    }

    #[test]
    fn test_kzg_batch_verify_all_valid() {
        let params = test_params();
        let polys = vec![
            DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]),
            DensePolynomial::new(vec![Fr::from(3u64), Fr::from(4u64), Fr::from(5u64)]),
        ];
        let points = vec![Fr::from(7u64), Fr::from(11u64)];

        let commitments: Vec<_> = polys.iter().map(|p| kzg_commit(&params, p)).collect();
        let mut values = Vec::new();
        let mut proofs = Vec::new();
        for (poly, z) in polys.iter().zip(points.iter()) {
            let (v, pi) = kzg_open(&params, poly, *z);
            values.push(v);
            proofs.push(pi);
        }

        assert!(kzg_batch_verify(
            &params,
            &commitments,
            &points,
            &values,
            &proofs,
        ));
    }
}
