pub mod proof;
pub mod prover;
pub mod setup;
pub mod verifier;

pub use proof::{Proof, ProvingKey, VerifyingKey};
pub use prover::prove;
pub use setup::setup;
pub use verifier::verify;

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_fields::{FieldElement, Fr};
    use zkrust_r1cs::{ConstraintSystem, LinearCombination};

    /// Helper: build a simple "x * x = y" circuit.
    fn squaring_circuit(x_val: u64) -> ConstraintSystem {
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(x_val));
        let y = cs.alloc_public_input(Fr::from(x_val * x_val));
        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
        );
        cs
    }

    fn test_toxic_waste() -> (Fr, Fr, Fr, Fr, Fr) {
        (
            Fr::from(17u64), // tau
            Fr::from(5u64),  // alpha
            Fr::from(7u64),  // beta
            Fr::from(11u64), // gamma
            Fr::from(13u64), // delta
        )
    }

    #[test]
    fn test_prove_and_verify_squaring() {
        let cs = squaring_circuit(3);
        assert!(cs.is_satisfied());

        let (tau, alpha, beta, gamma, delta) = test_toxic_waste();
        let pk = setup(&cs, tau, alpha, beta, gamma, delta);

        let mut rng = rand::thread_rng();
        let proof = prove(&pk, &cs, &mut rng);

        let public_inputs = vec![Fr::from(3u64), Fr::from(9u64)];
        assert!(verify(&pk.vk, &public_inputs, &proof));
    }

    #[test]
    fn test_wrong_public_input_rejected() {
        let cs = squaring_circuit(3);
        let (tau, alpha, beta, gamma, delta) = test_toxic_waste();
        let pk = setup(&cs, tau, alpha, beta, gamma, delta);

        let mut rng = rand::thread_rng();
        let proof = prove(&pk, &cs, &mut rng);

        // Wrong public input
        let bad_inputs = vec![Fr::from(3u64), Fr::from(10u64)];
        assert!(!verify(&pk.vk, &bad_inputs, &proof));
    }

    #[test]
    fn test_different_valid_witnesses() {
        // Same circuit structure, different values
        let cs1 = squaring_circuit(3);
        let cs2 = squaring_circuit(7);

        let (tau, alpha, beta, gamma, delta) = test_toxic_waste();

        // Setup with cs1 (both have same structure)
        let pk = setup(&cs1, tau, alpha, beta, gamma, delta);

        let mut rng = rand::thread_rng();

        let proof1 = prove(&pk, &cs1, &mut rng);
        assert!(verify(&pk.vk, &[Fr::from(3u64), Fr::from(9u64)], &proof1));

        let proof2 = prove(&pk, &cs2, &mut rng);
        assert!(verify(&pk.vk, &[Fr::from(7u64), Fr::from(49u64)], &proof2));
    }

    #[test]
    fn test_quadratic_circuit() {
        // x^2 + x + 5 = y
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_public_input(Fr::from(17u64)); // 9 + 3 + 5 = 17
        let x_sq = cs.alloc_private(Fr::from(9u64));

        // x_sq = x * x
        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(x_sq),
        );

        // (x_sq + x + 5) * 1 = y
        let lhs = LinearCombination::from_variable(x_sq)
            + LinearCombination::from_variable(x)
            + LinearCombination::from_constant(Fr::from(5u64));
        cs.enforce(
            lhs,
            LinearCombination::from_constant(Fr::ONE),
            LinearCombination::from_variable(y),
        );

        assert!(cs.is_satisfied());

        let (tau, alpha, beta, gamma, delta) = test_toxic_waste();
        let pk = setup(&cs, tau, alpha, beta, gamma, delta);

        let mut rng = rand::thread_rng();
        let proof = prove(&pk, &cs, &mut rng);

        let public_inputs = vec![Fr::from(3u64), Fr::from(17u64)];
        assert!(verify(&pk.vk, &public_inputs, &proof));
    }
}
