use crate::FieldElement;
use rand::RngCore;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Base field element for BN254.
///
/// Internally stored in Montgomery form as 4 × u64 limbs (little-endian).
/// All arithmetic is performed modulo the BN254 base field prime:
///   p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
///     = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fp(pub(crate) [u64; 4]);

/// The BN254 base field modulus p in little-endian u64 limbs.
pub const MODULUS: [u64; 4] = [
    0x3C208C16D87CFD47,
    0x97816A916871CA8D,
    0xB85045B68181585D,
    0x30644E72E131A029,
];

/// R = 2^256 mod p (Montgomery form of 1).
const R: [u64; 4] = [
    0xD35D438DC58F0D9D,
    0x0A78EB28F5C70B3D,
    0x666EA36F7879462C,
    0x0E0A77C19A07DF2F,
];

/// R² = 2^512 mod p (used for converting to Montgomery form).
const R2: [u64; 4] = [
    0xF32CFC5B538AFA89,
    0xB5E71911D44501FB,
    0x47AB1EFF0A417FF6,
    0x06D89F71CAB8351F,
];

/// INV = -p^(-1) mod 2^64 (Montgomery reduction constant).
const INV: u64 = 0x87D20782E4866389;

impl Fp {
    /// Create an Fp element from a raw u64 array (NOT in Montgomery form).
    /// The input is interpreted as a little-endian representation.
    pub fn from_raw(val: [u64; 4]) -> Self {
        // Convert to Montgomery form by multiplying by R²
        Self(val).mul_mont(&Fp(R2))
    }

    /// Convert from Montgomery form back to standard representation.
    pub fn to_raw(&self) -> [u64; 4] {
        // Multiply by 1 (in standard form) to undo Montgomery
        self.mul_mont(&Fp([1, 0, 0, 0])).0
    }

    /// Montgomery multiplication: computes (a * b * R^{-1}) mod p
    /// using the CIOS (Coarsely Integrated Operand Scanning) method.
    ///
    /// This is the performance-critical inner loop of all field arithmetic.
    #[inline(always)]
    #[allow(clippy::needless_range_loop)]
    fn mul_mont(&self, other: &Self) -> Self {
        let a = &self.0;
        let b = &other.0;
        let mut t = [0u64; 5]; // 4 limbs + 1 carry

        for i in 0..4 {
            // Multiply a[i] through all of b, accumulating into t
            let mut carry: u64 = 0;
            for j in 0..4 {
                let product = (a[i] as u128) * (b[j] as u128) + (t[j] as u128) + (carry as u128);
                t[j] = product as u64;
                carry = (product >> 64) as u64;
            }
            t[4] = carry;

            // Montgomery reduction step
            let m = t[0].wrapping_mul(INV);
            // Add m * MODULUS to t
            let product = (m as u128) * (MODULUS[0] as u128) + (t[0] as u128);
            let mut carry: u64 = (product >> 64) as u64;

            for j in 1..4 {
                let product = (m as u128) * (MODULUS[j] as u128) + (t[j] as u128) + (carry as u128);
                t[j - 1] = product as u64;
                carry = (product >> 64) as u64;
            }

            let sum = (t[4] as u128) + (carry as u128);
            t[3] = sum as u64;
            t[4] = (sum >> 64) as u64;
        }

        // Final conditional subtraction
        let result = [t[0], t[1], t[2], t[3]];
        if t[4] != 0 || !lt(&result, &MODULUS) {
            sub_mod(&result, &MODULUS)
        } else {
            Fp(result)
        }
    }

    /// Compute self^exp where exp is given as a little-endian byte array.
    fn pow_vartime(&self, exp: &[u64]) -> Self {
        let mut result = Fp(R); // 1 in Montgomery form
        let mut base = *self;

        for &limb in exp {
            let mut e = limb;
            for _ in 0..64 {
                if e & 1 == 1 {
                    result = result.mul_mont(&base);
                }
                base = base.mul_mont(&base);
                e >>= 1;
            }
        }

        result
    }

    /// Compute the Legendre symbol: returns 1 if self is a QR, -1 if NQR, 0 if zero.
    pub fn legendre(&self) -> i8 {
        // Compute self^((p-1)/2) mod p
        let exp = [
            0x9E10460B6C3E7EA3,
            0xCBC0B548B438E546,
            0xDC2822DB40C0AC2E,
            0x183227397098D014,
        ];
        let s = self.pow_vartime(&exp);
        if s.is_zero() {
            0
        } else if s == Fp(R) {
            1
        } else {
            -1
        }
    }

    /// Compute the square root of self, if it exists.
    /// Since p ≡ 3 (mod 4), we can use sqrt(a) = a^((p+1)/4).
    pub fn sqrt(&self) -> Option<Self> {
        // (p + 1) / 4
        let exp = [
            0x4F082305B61F3F52,
            0x65E05AA45A1C72A3,
            0x6E14116DA0605617,
            0x0C19139CB84C680A,
        ];
        let root = self.pow_vartime(&exp);
        if root.mul_mont(&root) == *self {
            Some(root)
        } else {
            None
        }
    }
}

impl FieldElement for Fp {
    const ZERO: Self = Fp([0; 4]);
    const ONE: Self = Fp(R);

    #[inline]
    fn add(&self, other: &Self) -> Self {
        let (mut result, carry) = adc_array(&self.0, &other.0);
        if carry || !lt(&result, &MODULUS) {
            result = sub_mod(&result, &MODULUS).0;
        }
        Fp(result)
    }

    #[inline]
    fn sub(&self, other: &Self) -> Self {
        if lt(&self.0, &other.0) {
            // self < other, so add modulus: (self + p) - other
            let (sum, _) = adc_array(&self.0, &MODULUS);
            sub_mod(&sum, &other.0)
        } else {
            sub_mod(&self.0, &other.0)
        }
    }

    #[inline]
    fn mul(&self, other: &Self) -> Self {
        self.mul_mont(other)
    }

    #[inline]
    fn neg(&self) -> Self {
        if self.is_zero() {
            *self
        } else {
            sub_mod(&MODULUS, &self.0)
        }
    }

    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        // Fermat's little theorem: a^(-1) = a^(p-2) mod p
        let exp = [
            0x3C208C16D87CFD45, // MODULUS[0] - 2
            0x97816A916871CA8D,
            0xB85045B68181585D,
            0x30644E72E131A029,
        ];
        Some(self.pow_vartime(&exp))
    }

    #[inline]
    fn square(&self) -> Self {
        self.mul_mont(self)
    }

    fn pow(&self, exp: &[u64]) -> Self {
        self.pow_vartime(exp)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == [0u64; 4]
    }

    fn random(rng: &mut impl RngCore) -> Self {
        loop {
            let mut limbs = [0u64; 4];
            for limb in &mut limbs {
                *limb = rng.next_u64();
            }
            // Mask top limb to be at most ~modulus size
            limbs[3] &= 0x3FFFFFFFFFFFFFFF;

            if lt(&limbs, &MODULUS) {
                return Fp::from_raw(limbs);
            }
        }
    }

    fn double(&self) -> Self {
        self.add(self)
    }
}

impl From<u64> for Fp {
    fn from(val: u64) -> Self {
        Fp::from_raw([val, 0, 0, 0])
    }
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let raw = self.to_raw();
        write!(
            f,
            "Fp(0x{:016X}{:016X}{:016X}{:016X})",
            raw[3], raw[2], raw[1], raw[0]
        )
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let raw = self.to_raw();
        write!(
            f,
            "0x{:016X}{:016X}{:016X}{:016X}",
            raw[3], raw[2], raw[1], raw[0]
        )
    }
}

// Operator trait implementations

impl Add for Fp {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        FieldElement::add(&self, &rhs)
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        FieldElement::sub(&self, &rhs)
    }
}

impl Mul for Fp {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        FieldElement::mul(&self, &rhs)
    }
}

impl Neg for Fp {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        FieldElement::neg(&self)
    }
}

impl AddAssign for Fp {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = FieldElement::add(self, &rhs);
    }
}

impl SubAssign for Fp {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = FieldElement::sub(self, &rhs);
    }
}

impl MulAssign for Fp {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = FieldElement::mul(self, &rhs);
    }
}

// ---- Helper functions ----

/// Add two 4-limb numbers, returning the result and carry bit.
#[inline]
fn adc_array(a: &[u64; 4], b: &[u64; 4]) -> ([u64; 4], bool) {
    let mut result = [0u64; 4];
    let mut carry = false;

    for i in 0..4 {
        let (sum, c1) = a[i].overflowing_add(b[i]);
        let (sum, c2) = sum.overflowing_add(carry as u64);
        result[i] = sum;
        carry = c1 || c2;
    }

    (result, carry)
}

/// Subtract b from a (assumes a >= b), returning the Fp result.
#[inline]
fn sub_mod(a: &[u64; 4], b: &[u64; 4]) -> Fp {
    let mut result = [0u64; 4];
    let mut borrow = false;

    for i in 0..4 {
        let (diff, b1) = a[i].overflowing_sub(b[i]);
        let (diff, b2) = diff.overflowing_sub(borrow as u64);
        result[i] = diff;
        borrow = b1 || b2;
    }

    Fp(result)
}

/// Check if a < b (little-endian limb comparison).
#[inline]
fn lt(a: &[u64; 4], b: &[u64; 4]) -> bool {
    for i in (0..4).rev() {
        if a[i] < b[i] {
            return true;
        }
        if a[i] > b[i] {
            return false;
        }
    }
    false // equal
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_and_one() {
        assert!(Fp::ZERO.is_zero());
        assert!(!Fp::ONE.is_zero());
        assert_ne!(Fp::ZERO, Fp::ONE);
    }

    #[test]
    fn test_from_u64() {
        let a = Fp::from(0u64);
        assert!(a.is_zero());

        let one = Fp::from(1u64);
        assert_eq!(one, Fp::ONE);
    }

    #[test]
    fn test_montgomery_roundtrip() {
        let val = [42u64, 0, 0, 0];
        let fp = Fp::from_raw(val);
        let back = fp.to_raw();
        assert_eq!(back, val);
    }

    #[test]
    fn test_montgomery_roundtrip_large() {
        // A value smaller than the modulus
        let val = [
            0x1234567890ABCDEF,
            0xFEDCBA0987654321,
            0x1111111111111111,
            0x0000000000000001,
        ];
        let fp = Fp::from_raw(val);
        let back = fp.to_raw();
        assert_eq!(back, val);
    }

    #[test]
    fn test_addition_identity() {
        let a = Fp::from(42u64);
        assert_eq!(a + Fp::ZERO, a);
        assert_eq!(Fp::ZERO + a, a);
    }

    #[test]
    fn test_addition_commutativity() {
        let a = Fp::from(100u64);
        let b = Fp::from(200u64);
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn test_subtraction() {
        let a = Fp::from(100u64);
        let b = Fp::from(42u64);
        let c = a - b;
        assert_eq!(c + b, a);
    }

    #[test]
    fn test_subtraction_underflow() {
        // b > a in normal integers, but wraps in the field
        let a = Fp::from(10u64);
        let b = Fp::from(20u64);
        let c = a - b;
        // c + b should equal a
        assert_eq!(c + b, a);
    }

    #[test]
    fn test_multiplication_identity() {
        let a = Fp::from(42u64);
        assert_eq!(a * Fp::ONE, a);
        assert_eq!(Fp::ONE * a, a);
    }

    #[test]
    fn test_multiplication_zero() {
        let a = Fp::from(42u64);
        assert_eq!(a * Fp::ZERO, Fp::ZERO);
    }

    #[test]
    fn test_multiplication_commutativity() {
        let a = Fp::from(123u64);
        let b = Fp::from(456u64);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_multiplication_associativity() {
        let a = Fp::from(3u64);
        let b = Fp::from(5u64);
        let c = Fp::from(7u64);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_distributivity() {
        let a = Fp::from(3u64);
        let b = Fp::from(5u64);
        let c = Fp::from(7u64);
        assert_eq!(a * (b + c), a * b + a * c);
    }

    #[test]
    fn test_negation() {
        let a = Fp::from(42u64);
        let neg_a = -a;
        assert_eq!(a + neg_a, Fp::ZERO);
        assert_eq!(-Fp::ZERO, Fp::ZERO);
    }

    #[test]
    fn test_inversion() {
        let a = Fp::from(42u64);
        let a_inv = a.inv().unwrap();
        assert_eq!(a * a_inv, Fp::ONE);
    }

    #[test]
    fn test_inversion_of_one() {
        assert_eq!(Fp::ONE.inv().unwrap(), Fp::ONE);
    }

    #[test]
    fn test_inversion_of_zero() {
        assert!(Fp::ZERO.inv().is_none());
    }

    #[test]
    fn test_squaring() {
        let a = Fp::from(7u64);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_pow() {
        let a = Fp::from(2u64);
        // 2^10 = 1024
        let result = a.pow(&[10, 0, 0, 0]);
        assert_eq!(result, Fp::from(1024u64));
    }

    #[test]
    fn test_pow_zero() {
        let a = Fp::from(42u64);
        assert_eq!(a.pow(&[0, 0, 0, 0]), Fp::ONE);
    }

    #[test]
    fn test_sqrt() {
        // 4 has sqrt 2
        let four = Fp::from(4u64);
        let root = four.sqrt().unwrap();
        assert_eq!(root * root, four);
    }

    #[test]
    fn test_sqrt_one() {
        let one = Fp::ONE;
        let root = one.sqrt().unwrap();
        assert_eq!(root * root, one);
    }

    #[test]
    fn test_legendre() {
        assert_eq!(Fp::ZERO.legendre(), 0);
        assert_eq!(Fp::ONE.legendre(), 1);
        // 4 is a perfect square
        assert_eq!(Fp::from(4u64).legendre(), 1);
    }

    #[test]
    fn test_random_element_properties() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let a = Fp::random(&mut rng);
            let b = Fp::random(&mut rng);

            // Commutativity
            assert_eq!(a + b, b + a);
            assert_eq!(a * b, b * a);

            // Inverse
            if !a.is_zero() {
                let a_inv = a.inv().unwrap();
                assert_eq!(a * a_inv, Fp::ONE);
            }
        }
    }

    #[test]
    fn test_double() {
        let a = Fp::from(21u64);
        assert_eq!(a.double(), a + a);
        assert_eq!(a.double(), Fp::from(42u64));
    }

    #[test]
    fn test_assign_ops() {
        let mut a = Fp::from(10u64);
        let b = Fp::from(5u64);

        a += b;
        assert_eq!(a, Fp::from(15u64));

        a -= b;
        assert_eq!(a, Fp::from(10u64));

        a *= b;
        assert_eq!(a, Fp::from(50u64));
    }
}
