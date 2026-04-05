//! Range proof example: prove that a secret value x is less than 2^32.
//!
//! Run: cargo run --example range_proof

use zkrust_fields::{FieldElement, Fr};
use zkrust_r1cs::{enforce_range, ConstraintSystem};

fn main() {
    let secret_value = 42u64;

    println!("=== Range Proof Example ===");
    println!("Proving that x = {} is in range [0, 2^32)", secret_value);
    println!();

    // Build the constraint system
    let mut cs = ConstraintSystem::new();
    let x = cs.alloc_public_input(Fr::from(secret_value));
    enforce_range(&mut cs, x, Fr::from(secret_value), 32);

    println!(
        "Circuit: {} constraints, {} variables",
        cs.num_constraints(),
        cs.num_variables()
    );
    assert!(cs.is_satisfied(), "constraints not satisfied");

    // Setup
    let mut rng = rand::thread_rng();
    let tau = Fr::random(&mut rng);
    let alpha = Fr::random(&mut rng);
    let beta = Fr::random(&mut rng);
    let gamma = Fr::random(&mut rng);
    let delta = Fr::random(&mut rng);

    let pk = zkrust_groth16::setup(&cs, tau, alpha, beta, gamma, delta);
    println!("Setup complete.");

    // Prove
    let proof = zkrust_groth16::prove(&pk, &cs, &mut rng);
    println!("Proof generated.");

    // Verify
    let public_inputs = vec![Fr::from(secret_value)];
    let valid = zkrust_groth16::verify(&pk.vk, &public_inputs, &proof);
    println!("Verification: {}", if valid { "VALID" } else { "INVALID" });
    assert!(valid);

    // Try with wrong input
    let wrong_inputs = vec![Fr::from(99u64)];
    let invalid = zkrust_groth16::verify(&pk.vk, &wrong_inputs, &proof);
    println!(
        "Wrong input verification: {} (expected INVALID)",
        if invalid { "VALID" } else { "INVALID" }
    );
    assert!(!invalid);

    println!("\nDone!");
}
