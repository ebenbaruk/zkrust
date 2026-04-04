/// Base field element for BN254.
///
/// Stored in Montgomery form as 4 × u64 limbs.
/// Modulus p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Fp(pub(crate) [u64; 4]);

impl Fp {
    pub const ZERO: Self = Self([0; 4]);
}
