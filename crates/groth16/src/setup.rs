use crate::proof::{ProvingKey, VerifyingKey};
use zkrust_curves::{G1Affine, G2Affine};
use zkrust_fields::{FieldElement, Fr};
use zkrust_r1cs::ConstraintSystem;

/// Generate Groth16 proving and verifying keys from a constraint system.
///
/// WARNING: This uses the provided toxic waste directly. In production,
/// toxic waste must come from a secure multi-party computation ceremony.
pub fn setup(
    cs: &ConstraintSystem,
    tau: Fr,
    alpha: Fr,
    beta: Fr,
    gamma: Fr,
    delta: Fr,
) -> ProvingKey {
    let num_public = cs.num_public_inputs();
    let num_vars = cs.num_variables(); // includes the "1" variable
    let num_constraints = cs.num_constraints();

    // Domain must be a power of 2 for NTT
    let domain_size = num_constraints.next_power_of_two();

    let g1 = G1Affine::generator();
    let g2 = G2Affine::generator();

    // Precompute group elements for toxic waste
    let alpha_g1 = (g1 * alpha).to_affine();
    let beta_g1 = (g1 * beta).to_affine();
    let beta_g2 = (g2 * beta).to_affine();
    let gamma_g2 = (g2 * gamma).to_affine();
    let delta_g1 = (g1 * delta).to_affine();
    let delta_g2 = (g2 * delta).to_affine();

    let gamma_inv = gamma.inv().unwrap();
    let delta_inv = delta.inv().unwrap();

    // Compute Lagrange basis evaluations at tau:
    // L_j(τ) = (τ^m - 1) / (m · (τ - ω^j))
    // where m = domain_size, ω = primitive m-th root of unity
    let lagrange_at_tau = compute_lagrange_at_tau(tau, domain_size);

    // Extract sparse matrices
    let (a_mat, b_mat, c_mat) = cs.to_matrices();

    // Compute u_i(τ), v_i(τ), w_i(τ) for each variable i
    // u_i(τ) = ∑_j A[j][i] · L_j(τ)
    let mut u_at_tau = vec![Fr::ZERO; num_vars];
    let mut v_at_tau = vec![Fr::ZERO; num_vars];
    let mut w_at_tau = vec![Fr::ZERO; num_vars];

    for (j, row) in a_mat.iter().enumerate() {
        for &(col, coeff) in row {
            u_at_tau[col] += coeff * lagrange_at_tau[j];
        }
    }
    for (j, row) in b_mat.iter().enumerate() {
        for &(col, coeff) in row {
            v_at_tau[col] += coeff * lagrange_at_tau[j];
        }
    }
    for (j, row) in c_mat.iter().enumerate() {
        for &(col, coeff) in row {
            w_at_tau[col] += coeff * lagrange_at_tau[j];
        }
    }

    // [u_i(τ)]_1, [v_i(τ)]_1, [v_i(τ)]_2 for all variables
    let a_g1: Vec<G1Affine> = u_at_tau.iter().map(|u| (g1 * *u).to_affine()).collect();
    let b_g1: Vec<G1Affine> = v_at_tau.iter().map(|v| (g1 * *v).to_affine()).collect();
    let b_g2: Vec<G2Affine> = v_at_tau.iter().map(|v| (g2 * *v).to_affine()).collect();

    // Compute combined terms: β·u_i(τ) + α·v_i(τ) + w_i(τ)
    let combined: Vec<Fr> = (0..num_vars)
        .map(|i| beta * u_at_tau[i] + alpha * v_at_tau[i] + w_at_tau[i])
        .collect();

    // IC (for verifying key): public variables (indices 0..num_public+1, including "1" var)
    // [(β·u_i(τ) + α·v_i(τ) + w_i(τ))/γ]_1
    let num_public_total = 1 + num_public; // includes the "1" variable
    let ic: Vec<G1Affine> = (0..num_public_total)
        .map(|i| (g1 * (combined[i] * gamma_inv)).to_affine())
        .collect();

    // L (for proving key): private variables (indices num_public_total..num_vars)
    // [(β·u_i(τ) + α·v_i(τ) + w_i(τ))/δ]_1
    let l_g1: Vec<G1Affine> = (num_public_total..num_vars)
        .map(|i| (g1 * (combined[i] * delta_inv)).to_affine())
        .collect();

    // Compute [τ^i · t(τ) / δ]_1 for i = 0..domain_size-1
    // t(τ) = τ^domain_size - 1
    let t_at_tau = tau.pow(&[domain_size as u64, 0, 0, 0]) - Fr::ONE;
    let t_delta_inv = t_at_tau * delta_inv;

    let mut h_g1 = Vec::with_capacity(domain_size);
    let mut tau_pow = Fr::ONE;
    for _ in 0..domain_size {
        h_g1.push((g1 * (tau_pow * t_delta_inv)).to_affine());
        tau_pow *= tau;
    }

    let vk = VerifyingKey {
        alpha_g1,
        beta_g2,
        gamma_g2,
        delta_g2,
        ic,
    };

    ProvingKey {
        alpha_g1,
        beta_g1,
        beta_g2,
        delta_g1,
        delta_g2,
        a_g1,
        b_g1,
        b_g2,
        l_g1,
        h_g1,
        vk,
        num_public,
        domain_size,
    }
}

/// Compute Lagrange basis evaluations at tau for a domain of roots of unity.
///
/// L_j(τ) = (τ^m - 1) / (m · (τ - ω^j)) for j = 0..m-1
fn compute_lagrange_at_tau(tau: Fr, domain_size: usize) -> Vec<Fr> {
    let m = domain_size;
    let k = m.trailing_zeros();
    let omega = Fr::root_of_unity(k);

    let tau_m = tau.pow(&[m as u64, 0, 0, 0]);
    let numerator = tau_m - Fr::ONE;
    let m_inv = Fr::from(m as u64).inv().unwrap();

    let mut result = Vec::with_capacity(m);
    let mut omega_j = Fr::ONE;
    for _ in 0..m {
        let denom = tau - omega_j;
        if denom.is_zero() {
            // tau = omega^j, so L_j(tau) = 1, all others = 0
            // This is a degenerate case; shouldn't happen with random tau
            result.push(Fr::ONE);
        } else {
            result.push(omega_j * numerator * m_inv * denom.inv().unwrap());
        }
        omega_j *= omega;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_r1cs::{ConstraintSystem, LinearCombination};

    #[test]
    fn test_setup_produces_keys() {
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_private(Fr::from(9u64));
        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
        );

        let tau = Fr::from(17u64);
        let alpha = Fr::from(5u64);
        let beta = Fr::from(7u64);
        let gamma = Fr::from(11u64);
        let delta = Fr::from(13u64);

        let pk = setup(&cs, tau, alpha, beta, gamma, delta);

        assert_eq!(pk.num_public, 1);
        assert_eq!(pk.vk.ic.len(), 2); // "1" var + 1 public input
        assert!(!pk.l_g1.is_empty()); // private variable terms
        assert!(!pk.h_g1.is_empty()); // h polynomial terms
    }

    #[test]
    fn test_lagrange_at_tau() {
        // Lagrange basis should sum to 1 at any point
        let domain_size = 4;
        let tau = Fr::from(42u64);
        let lagrange = compute_lagrange_at_tau(tau, domain_size);

        let sum: Fr = lagrange.iter().fold(Fr::ZERO, |acc, x| acc + *x);
        assert_eq!(sum, Fr::ONE);
    }

    #[test]
    fn test_lagrange_at_root() {
        // L_j(ω^k) = delta_{j,k}
        let domain_size = 4;
        let omega = Fr::root_of_unity(2);
        let tau = omega; // evaluate at ω^0 ... wait, ω^0 = 1... tau = omega = ω^1

        let lagrange = compute_lagrange_at_tau(tau, domain_size);
        // L_1(ω) should be 1, others 0
        assert_eq!(lagrange[1], Fr::ONE);
        for (i, l) in lagrange.iter().enumerate() {
            if i != 1 {
                assert!(l.is_zero(), "L_{i}(omega) should be 0, got {:?}", l);
            }
        }
    }
}
