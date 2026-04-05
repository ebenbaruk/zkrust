use crate::proof::{Proof, ProvingKey};
use rand::RngCore;
use zkrust_curves::{msm_g1, msm_g2};
use zkrust_fields::{FieldElement, Fr};
use zkrust_polynomials::{intt, ntt_mul};
use zkrust_r1cs::ConstraintSystem;

/// Generate a Groth16 proof for a satisfied constraint system.
pub fn prove(pk: &ProvingKey, cs: &ConstraintSystem, rng: &mut impl RngCore) -> Proof {
    let witness = cs.witness();
    let num_public_total = 1 + pk.num_public;

    // Random blinding factors
    let r = Fr::random(rng);
    let s = Fr::random(rng);

    // Compute h(x) polynomial
    let h_coeffs = compute_h(cs, pk.domain_size);

    // Proof element [A]_1 = [α]_1 + ∑ z_i·[u_i(τ)]_1 + r·[δ]_1
    let a_msm = msm_g1(&pk.a_g1, &witness);
    let a = (pk.alpha_g1.to_projective() + a_msm + (pk.delta_g1 * r)).to_affine();

    // Proof element [B]_2 = [β]_2 + ∑ z_i·[v_i(τ)]_2 + s·[δ]_2
    let b_msm = msm_g2(&pk.b_g2, &witness);
    let b = (pk.beta_g2.to_projective() + b_msm + (pk.delta_g2 * s)).to_affine();

    // [B]_1 = [β]_1 + ∑ z_i·[v_i(τ)]_1 + s·[δ]_1  (needed for [C]_1 computation)
    let b1_msm = msm_g1(&pk.b_g1, &witness);
    let b1 = (pk.beta_g1.to_projective() + b1_msm + (pk.delta_g1 * s)).to_affine();

    // Private witness part for [C]_1
    let private_witness: Vec<Fr> = witness[num_public_total..].to_vec();
    let l_msm = msm_g1(&pk.l_g1, &private_witness);

    // h(x) contribution: ∑ h_i · [τ^i·t(τ)/δ]_1
    let h_len = h_coeffs.len().min(pk.h_g1.len());
    let h_msm = if h_len > 0 {
        msm_g1(&pk.h_g1[..h_len], &h_coeffs[..h_len])
    } else {
        zkrust_curves::G1Projective::identity()
    };

    // [C]_1 = ∑ z_private·[L_i]_1 + ∑ h_i·[τ^i·t(τ)/δ]_1 + s·[A]_1 + r·[B]_1 - r·s·[δ]_1
    let c = (l_msm + h_msm + (a * s) + (b1 * r) - (pk.delta_g1 * (r * s))).to_affine();

    Proof { a, b, c }
}

/// Compute the h(x) polynomial such that A(x)·B(x) - C(x) = h(x)·t(x).
///
/// A(x) = ∑ z_i · u_i(x), evaluated at roots of unity = Az (matrix-vector product).
fn compute_h(cs: &ConstraintSystem, domain_size: usize) -> Vec<Fr> {
    let witness = cs.witness();
    let (a_mat, b_mat, c_mat) = cs.to_matrices();
    let _num_constraints = a_mat.len();

    // Compute Az, Bz, Cz evaluations at roots of unity
    let mut a_evals = vec![Fr::ZERO; domain_size];
    let mut b_evals = vec![Fr::ZERO; domain_size];
    let mut c_evals = vec![Fr::ZERO; domain_size];

    for (j, row) in a_mat.iter().enumerate() {
        for &(col, coeff) in row {
            a_evals[j] += coeff * witness[col];
        }
    }
    for (j, row) in b_mat.iter().enumerate() {
        for &(col, coeff) in row {
            b_evals[j] += coeff * witness[col];
        }
    }
    for (j, row) in c_mat.iter().enumerate() {
        for &(col, coeff) in row {
            c_evals[j] += coeff * witness[col];
        }
    }
    // Padded entries (j >= num_constraints) are already zero

    // Convert from evaluations to coefficient form via INTT
    let a_coeffs = intt(&a_evals);
    let b_coeffs = intt(&b_evals);
    let c_coeffs = intt(&c_evals);

    // Multiply A and B polynomials using NTT
    let ab_coeffs = ntt_mul(&a_coeffs, &b_coeffs);

    // Compute AB - C
    let mut numerator = ab_coeffs;
    for (i, c) in c_coeffs.iter().enumerate() {
        if i < numerator.len() {
            numerator[i] -= *c;
        }
    }

    // Divide by t(x) = x^domain_size - 1
    // t(x) divides the numerator exactly since the R1CS is satisfied
    divide_by_vanishing(&numerator, domain_size)
}

/// Divide a polynomial by x^n - 1.
///
/// Given f(x) such that (x^n - 1) | f(x),
/// compute q(x) = f(x) / (x^n - 1).
///
/// Uses the relation: f[i] = q[i-n] - q[i], so q[i-n] = f[i] + q[i].
/// Process from highest degree down.
fn divide_by_vanishing(coeffs: &[Fr], n: usize) -> Vec<Fr> {
    let deg = coeffs.len();
    let q_len = if deg > n { deg - n } else { return vec![] };

    let mut q = vec![Fr::ZERO; q_len];

    for i in (n..deg).rev() {
        let q_high = if i < q_len { q[i] } else { Fr::ZERO };
        q[i - n] = coeffs[i] + q_high;
    }

    while q.last().is_some_and(|x| x.is_zero()) {
        q.pop();
    }

    q
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_divide_by_vanishing() {
        // f(x) = x^4 - 1 = (x^2 - 1)(x^2 + 1)
        // Divide x^4 - 1 by x^2 - 1, should get x^2 + 1
        // But our function divides by x^n - 1 specifically
        // Let's test: (x^2 + 1)(x^2 - 1) = x^4 - 1
        // divide_by_vanishing(x^4 - 1, 2) should give x^2 + 1

        let f = vec![-Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ZERO, Fr::ONE]; // x^4 - 1
        let q = divide_by_vanishing(&f, 2);
        // q should be [1, 0, 1] = 1 + x^2
        assert_eq!(q.len(), 3);
        assert_eq!(q[0], Fr::ONE);
        assert!(q[1].is_zero());
        assert_eq!(q[2], Fr::ONE);
    }

    #[test]
    fn test_divide_by_vanishing_simple() {
        // f(x) = x^2 - 1 divided by x^2 - 1 = 1
        let f = vec![-Fr::ONE, Fr::ZERO, Fr::ONE];
        let q = divide_by_vanishing(&f, 2);
        assert_eq!(q.len(), 1);
        assert_eq!(q[0], Fr::ONE);
    }
}
