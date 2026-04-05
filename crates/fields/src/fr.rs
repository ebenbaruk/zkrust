use crate::FieldElement;
use rand::RngCore;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Scalar field element for BN254.
///
/// Internally stored in Montgomery form as 4 × u64 limbs (little-endian).
/// All arithmetic is performed modulo the BN254 scalar field order:
///   r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
///     = 0x30644E72E131A029B85045B68181585D2833E84879B9709143E1F593F0000001
///
/// The scalar field has two-adicity 28: r - 1 = 2^28 * t (where t is odd).
/// This is important for NTT-based polynomial multiplication.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fr(pub(crate) [u64; 4]);

/// The BN254 scalar field order r in little-endian u64 limbs.
pub const MODULUS: [u64; 4] = [
    0x43E1F593F0000001,
    0x2833E84879B97091,
    0xB85045B68181585D,
    0x30644E72E131A029,
];

/// R = 2^256 mod r (Montgomery form of 1).
const R: [u64; 4] = [
    0xAC96341C4FFFFFFB,
    0x36FC76959F60CD29,
    0x666EA36F7879462E,
    0x0E0A77C19A07DF2F,
];

/// R² = 2^512 mod r (used for converting to Montgomery form).
const R2: [u64; 4] = [
    0x1BB8E645AE216DA7,
    0x53FE3AB1E35C59E3,
    0x8C49833D53BB8085,
    0x0216D0B17F4E44A5,
];

/// INV = -r^(-1) mod 2^64 (Montgomery reduction constant).
const INV: u64 = 0xC2E1F593EFFFFFFF;

/// Two-adicity of the scalar field: r - 1 = 2^28 * t.
pub const TWO_ADICITY: u32 = 28;

/// Generator of the 2^28-th roots of unity subgroup.
/// This is used as the principal root of unity for NTT.
const TWO_ADIC_ROOT_OF_UNITY: [u64; 4] = [
    0x9BD61B6E725B19F0,
    0x402D111E41112ED4,
    0x00E0A7EB8EF62ABC,
    0x2A3C09F0A58A7E85,
];

impl Fr {
    /// Create an Fr element from a raw u64 array (NOT in Montgomery form).
    pub fn from_raw(val: [u64; 4]) -> Self {
        Self(val).mul_mont(&Fr(R2))
    }

    /// Convert from Montgomery form back to standard representation.
    pub fn to_raw(&self) -> [u64; 4] {
        self.mul_mont(&Fr([1, 0, 0, 0])).0
    }

    /// Montgomery multiplication using CIOS method.
    #[inline(always)]
    #[allow(clippy::needless_range_loop)]
    fn mul_mont(&self, other: &Self) -> Self {
        let a = &self.0;
        let b = &other.0;
        let mut t = [0u64; 5];

        for i in 0..4 {
            let mut carry: u64 = 0;
            for j in 0..4 {
                let product = (a[i] as u128) * (b[j] as u128) + (t[j] as u128) + (carry as u128);
                t[j] = product as u64;
                carry = (product >> 64) as u64;
            }
            t[4] = carry;

            let m = t[0].wrapping_mul(INV);
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

        let result = [t[0], t[1], t[2], t[3]];
        if t[4] != 0 || !lt(&result, &MODULUS) {
            sub_mod(&result, &MODULUS)
        } else {
            Fr(result)
        }
    }

    fn pow_vartime(&self, exp: &[u64]) -> Self {
        let mut result = Fr(R);
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

    /// Get a primitive 2^k-th root of unity in Fr.
    /// Panics if k > TWO_ADICITY (28).
    pub fn root_of_unity(k: u32) -> Self {
        assert!(k <= TWO_ADICITY, "k={k} exceeds two-adicity {TWO_ADICITY}");
        let mut root = Fr::from_raw(TWO_ADIC_ROOT_OF_UNITY);
        for _ in k..TWO_ADICITY {
            root = root.square();
        }
        root
    }

    /// Compute the square root using Tonelli-Shanks algorithm.
    /// The scalar field has two-adicity 28, so p ≡ 1 (mod 4) and we need
    /// the full Tonelli-Shanks rather than the simple (p+1)/4 exponent.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Fr::ZERO);
        }

        // r - 1 = 2^28 * t, where t is the odd part
        let t = [
            0x9B9709143E1F593F,
            0x181585D2833E8487,
            0x131A029B85045B68,
            0x000000030644E72E,
        ];

        // Find a non-residue: 5 is a non-residue in BN254 scalar field
        let non_residue = Fr::from(5u64);
        let mut c = non_residue.pow_vartime(&t); // c = g^t

        let mut r_val = self.pow_vartime(&[
            // (t + 1) / 2
            0xCDCB848A1F0FACA0,
            0x0C0AC2E9419F4243,
            0x098D014DC2822DB4,
            0x0000000183227397,
        ]);

        let mut t_val = self.pow_vartime(&t);
        let mut m = TWO_ADICITY;

        loop {
            if t_val == Fr(R) {
                // t_val == 1
                return Some(r_val);
            }

            // Find least i such that t_val^(2^i) == 1
            let mut i = 0u32;
            let mut tmp = t_val;
            while tmp != Fr(R) {
                tmp = tmp.square();
                i += 1;
                if i >= m {
                    return None; // Not a QR
                }
            }

            let b = {
                let mut b = c;
                for _ in 0..(m - i - 1) {
                    b = b.square();
                }
                b
            };

            r_val = r_val.mul_mont(&b);
            c = b.square();
            t_val = t_val.mul_mont(&c);
            m = i;
        }
    }
}

impl FieldElement for Fr {
    const ZERO: Self = Fr([0; 4]);
    const ONE: Self = Fr(R);

    #[inline]
    fn add(&self, other: &Self) -> Self {
        let (mut result, carry) = adc_array(&self.0, &other.0);
        if carry || !lt(&result, &MODULUS) {
            result = sub_mod(&result, &MODULUS).0;
        }
        Fr(result)
    }

    #[inline]
    fn sub(&self, other: &Self) -> Self {
        if lt(&self.0, &other.0) {
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
        let exp = [
            0x43E1F593EFFFFFFF, // MODULUS[0] - 2
            0x2833E84879B97091,
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
            limbs[3] &= 0x3FFFFFFFFFFFFFFF;

            if lt(&limbs, &MODULUS) {
                return Fr::from_raw(limbs);
            }
        }
    }

    fn double(&self) -> Self {
        self.add(self)
    }
}

impl From<u64> for Fr {
    fn from(val: u64) -> Self {
        Fr::from_raw([val, 0, 0, 0])
    }
}

impl fmt::Debug for Fr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let raw = self.to_raw();
        write!(
            f,
            "Fr(0x{:016X}{:016X}{:016X}{:016X})",
            raw[3], raw[2], raw[1], raw[0]
        )
    }
}

impl fmt::Display for Fr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let raw = self.to_raw();
        write!(
            f,
            "0x{:016X}{:016X}{:016X}{:016X}",
            raw[3], raw[2], raw[1], raw[0]
        )
    }
}

impl Add for Fr {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        FieldElement::add(&self, &rhs)
    }
}

impl Sub for Fr {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        FieldElement::sub(&self, &rhs)
    }
}

impl Mul for Fr {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        FieldElement::mul(&self, &rhs)
    }
}

impl Neg for Fr {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        FieldElement::neg(&self)
    }
}

impl AddAssign for Fr {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = FieldElement::add(self, &rhs);
    }
}

impl SubAssign for Fr {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = FieldElement::sub(self, &rhs);
    }
}

impl MulAssign for Fr {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = FieldElement::mul(self, &rhs);
    }
}

// ---- Helper functions ----

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

#[inline]
fn sub_mod(a: &[u64; 4], b: &[u64; 4]) -> Fr {
    let mut result = [0u64; 4];
    let mut borrow = false;
    for i in 0..4 {
        let (diff, b1) = a[i].overflowing_sub(b[i]);
        let (diff, b2) = diff.overflowing_sub(borrow as u64);
        result[i] = diff;
        borrow = b1 || b2;
    }
    Fr(result)
}

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
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_and_one() {
        assert!(Fr::ZERO.is_zero());
        assert!(!Fr::ONE.is_zero());
    }

    #[test]
    fn test_from_u64() {
        assert!(Fr::from(0u64).is_zero());
        assert_eq!(Fr::from(1u64), Fr::ONE);
    }

    #[test]
    fn test_montgomery_roundtrip() {
        let val = [42u64, 0, 0, 0];
        let fr = Fr::from_raw(val);
        assert_eq!(fr.to_raw(), val);
    }

    #[test]
    fn test_basic_arithmetic() {
        let a = Fr::from(100u64);
        let b = Fr::from(42u64);

        // Addition
        assert_eq!(a + b, Fr::from(142u64));

        // Subtraction
        assert_eq!(a - b, Fr::from(58u64));

        // Multiplication
        assert_eq!(a * b, Fr::from(4200u64));
    }

    #[test]
    fn test_subtraction_underflow() {
        let a = Fr::from(10u64);
        let b = Fr::from(20u64);
        let c = a - b;
        assert_eq!(c + b, a);
    }

    #[test]
    fn test_negation() {
        let a = Fr::from(42u64);
        assert_eq!(a + (-a), Fr::ZERO);
        assert_eq!(-Fr::ZERO, Fr::ZERO);
    }

    #[test]
    fn test_inversion() {
        let a = Fr::from(42u64);
        let a_inv = a.inv().unwrap();
        assert_eq!(a * a_inv, Fr::ONE);
        assert!(Fr::ZERO.inv().is_none());
    }

    #[test]
    fn test_pow() {
        let a = Fr::from(2u64);
        assert_eq!(a.pow(&[10, 0, 0, 0]), Fr::from(1024u64));
        assert_eq!(a.pow(&[0, 0, 0, 0]), Fr::ONE);
    }

    #[test]
    fn test_sqrt() {
        let four = Fr::from(4u64);
        let root = four.sqrt().unwrap();
        assert_eq!(root * root, four);

        // Test with a larger number
        let a = Fr::from(49u64);
        let root = a.sqrt().unwrap();
        assert_eq!(root * root, a);
    }

    #[test]
    fn test_root_of_unity() {
        // The 2^1-th root of unity squared should be 1
        let w = Fr::root_of_unity(1);
        assert_eq!(w.square(), Fr::ONE);
        assert_ne!(w, Fr::ONE); // It should be -1

        // The 2^k-th root of unity raised to 2^k should be 1
        let w = Fr::root_of_unity(4);
        let mut val = w;
        for _ in 1..16 {
            val = val * w;
        }
        assert_eq!(val, Fr::ONE);
    }

    #[test]
    fn test_random_algebraic_properties() {
        let mut rng = rand::thread_rng();
        for _ in 0..10 {
            let a = Fr::random(&mut rng);
            let b = Fr::random(&mut rng);
            let c = Fr::random(&mut rng);

            // Commutativity
            assert_eq!(a + b, b + a);
            assert_eq!(a * b, b * a);

            // Associativity
            assert_eq!((a + b) + c, a + (b + c));
            assert_eq!((a * b) * c, a * (b * c));

            // Distributivity
            assert_eq!(a * (b + c), a * b + a * c);

            // Inverse
            if !a.is_zero() {
                assert_eq!(a * a.inv().unwrap(), Fr::ONE);
            }
        }
    }
}
