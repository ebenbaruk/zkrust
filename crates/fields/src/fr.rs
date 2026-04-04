/// Scalar field element for BN254.
///
/// Stored in Montgomery form as 4 × u64 limbs.
/// Order r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Fr(pub(crate) [u64; 4]);

impl Fr {
    pub const ZERO: Self = Self([0; 4]);
}
