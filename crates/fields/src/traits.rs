use rand::RngCore;
use std::fmt::Debug;

/// Trait for elements of a finite field.
///
/// All ZK cryptographic operations are built on top of field arithmetic.
/// Implementors must provide constant-time (where possible) operations
/// over a prime field or extension field.
pub trait FieldElement: Sized + Clone + Copy + PartialEq + Eq + Debug {
    /// The additive identity (zero element).
    const ZERO: Self;
    /// The multiplicative identity (one element).
    const ONE: Self;

    /// Field addition.
    fn add(&self, other: &Self) -> Self;
    /// Field subtraction.
    fn sub(&self, other: &Self) -> Self;
    /// Field multiplication.
    fn mul(&self, other: &Self) -> Self;
    /// Additive inverse.
    fn neg(&self) -> Self;
    /// Multiplicative inverse. Returns `None` for zero.
    fn inv(&self) -> Option<Self>;
    /// Squaring (often faster than generic multiplication).
    fn square(&self) -> Self;
    /// Exponentiation by a multi-limb exponent (little-endian u64 limbs).
    fn pow(&self, exp: &[u64]) -> Self;
    /// Check if this element is zero.
    fn is_zero(&self) -> bool;
    /// Sample a uniformly random field element.
    fn random(rng: &mut impl RngCore) -> Self;

    /// Double this element (add to itself).
    fn double(&self) -> Self {
        self.add(self)
    }
}
