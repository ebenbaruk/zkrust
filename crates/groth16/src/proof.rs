use zkrust_curves::{G1Affine, G2Affine};

/// A Groth16 proof: three group elements.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Proof {
    pub a: G1Affine,
    pub b: G2Affine,
    pub c: G1Affine,
}

/// Proving key: everything the prover needs.
#[derive(Clone, Debug)]
pub struct ProvingKey {
    pub alpha_g1: G1Affine,
    pub beta_g1: G1Affine,
    pub beta_g2: G2Affine,
    pub delta_g1: G1Affine,
    pub delta_g2: G2Affine,

    /// [u_i(τ)]_1 for all variables i = 0..n
    pub a_g1: Vec<G1Affine>,
    /// [v_i(τ)]_1 for all variables i = 0..n
    pub b_g1: Vec<G1Affine>,
    /// [v_i(τ)]_2 for all variables i = 0..n
    pub b_g2: Vec<G2Affine>,

    /// [(β·u_i(τ) + α·v_i(τ) + w_i(τ))/δ]_1 for private variables
    pub l_g1: Vec<G1Affine>,

    /// [τ^i · t(τ) / δ]_1 for i = 0..domain_size-2
    pub h_g1: Vec<G1Affine>,

    /// The verifying key (included for convenience).
    pub vk: VerifyingKey,

    /// Number of public inputs (not counting the "1" variable).
    pub num_public: usize,
    /// Domain size (power of 2 ≥ num_constraints).
    pub domain_size: usize,
}

/// Verifying key: everything the verifier needs.
#[derive(Clone, Debug)]
pub struct VerifyingKey {
    pub alpha_g1: G1Affine,
    pub beta_g2: G2Affine,
    pub gamma_g2: G2Affine,
    pub delta_g2: G2Affine,

    /// IC[i] = [(β·u_i(τ) + α·v_i(τ) + w_i(τ))/γ]_1 for public variables i = 0..num_public+1
    /// (includes the "1" variable at index 0).
    pub ic: Vec<G1Affine>,
}
