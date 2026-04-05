use crate::proof::{Proof, VerifyingKey};
use zkrust_curves::ate_pairing;
use zkrust_fields::Fr;

/// Verify a Groth16 proof given public inputs.
///
/// Checks the pairing equation:
///   e(A, B) = e(α, β) · e(∑ pub_i · IC_i, γ) · e(C, δ)
pub fn verify(vk: &VerifyingKey, public_inputs: &[Fr], proof: &Proof) -> bool {
    assert_eq!(
        public_inputs.len() + 1,
        vk.ic.len(),
        "wrong number of public inputs: expected {}, got {}",
        vk.ic.len() - 1,
        public_inputs.len()
    );

    // Compute the public input contribution: IC[0] + ∑ pub_i · IC[i+1]
    let mut acc = vk.ic[0].to_projective();
    for (i, &input) in public_inputs.iter().enumerate() {
        acc = acc + (vk.ic[i + 1] * input);
    }
    let public_g1 = acc.to_affine();

    // Pairing check:
    // e(A, B) == e(alpha, beta) * e(public_g1, gamma) * e(C, delta)
    //
    // Equivalently:
    // e(A, B) * e(-alpha, beta) * e(-public_g1, gamma) * e(-C, delta) == 1
    //
    // Or we can compute both sides:
    let lhs = ate_pairing(&proof.a, &proof.b);

    let e_alpha_beta = ate_pairing(&vk.alpha_g1, &vk.beta_g2);
    let e_public_gamma = ate_pairing(&public_g1, &vk.gamma_g2);
    let e_c_delta = ate_pairing(&proof.c, &vk.delta_g2);

    let rhs = e_alpha_beta * e_public_gamma * e_c_delta;

    lhs == rhs
}

#[cfg(test)]
mod tests {
    // Integration tests are in lib.rs since they need setup + prove + verify
}
