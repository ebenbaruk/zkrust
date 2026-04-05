use crate::fp2::Fp2;
use crate::fp6::Fp6;
use crate::FieldElement;
use rand::RngCore;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Element of the dodecic extension field Fp12 = Fp6[w] / (w² - v).
///
/// Represented as c0 + c1*w where c0, c1 ∈ Fp6.
/// This is the target field of the BN254 pairing: e(G1, G2) → GT ⊂ Fp12*.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fp12 {
    pub c0: Fp6,
    pub c1: Fp6,
}

impl Fp12 {
    pub const fn new(c0: Fp6, c1: Fp6) -> Self {
        Self { c0, c1 }
    }

    /// Conjugate in Fp12: conjugate(c0 + c1*w) = c0 - c1*w.
    /// This is the unitary inverse for elements on the cyclotomic subgroup.
    pub fn conjugate(&self) -> Self {
        Self {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// Frobenius map: x -> x^p.
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let mut c1 = self.c1.frobenius_map();

        // Multiply c1 by the Frobenius constant for Fp12: ξ^((p-1)/6)
        let gamma = Fp2::new(
            crate::fp::Fp::from_raw([
                0xD60B35DADCC9E470,
                0x5C521E08292F2176,
                0xE8B99FDD76E68B60,
                0x1284B71C2865A7DF,
            ]),
            crate::fp::Fp::from_raw([
                0xCA5CF05F80F362AC,
                0x747992778EEEC7E5,
                0xA6327CFE12150B8E,
                0x246996F3B4FAE7E6,
            ]),
        );

        c1 = Fp6::new(c1.c0 * gamma, c1.c1 * gamma, c1.c2 * gamma);

        Self { c0, c1 }
    }

    /// Cyclotomic squaring for elements in the cyclotomic subgroup.
    pub fn cyclotomic_square(&self) -> Self {
        // TODO: implement Granger-Scott optimization
        self.square()
    }

    /// Exponentiation by the BN parameter x = 4965661367192848881.
    pub fn exp_by_x(&self) -> Self {
        let x: u64 = 4965661367192848881;
        let mut result = Self::ONE;
        let mut base = *self;
        let mut e = x;

        while e > 0 {
            if e & 1 == 1 {
                result *= base;
            }
            base = base.square();
            e >>= 1;
        }

        result
    }

    /// Multiply by a sparse Fp12 element from line functions.
    pub fn mul_by_034(&self, a0: &Fp2, a2: &Fp2, b1: &Fp2) -> Self {
        let sparse = Fp12::new(
            Fp6::new(*a0, Fp2::ZERO, *a2),
            Fp6::new(Fp2::ZERO, *b1, Fp2::ZERO),
        );
        *self * sparse
    }
}

impl FieldElement for Fp12 {
    const ZERO: Self = Fp12 {
        c0: Fp6::ZERO,
        c1: Fp6::ZERO,
    };
    const ONE: Self = Fp12 {
        c0: Fp6::ONE,
        c1: Fp6::ZERO,
    };

    fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }

    fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }

    /// Karatsuba multiplication for Fp12.
    /// w² = v in the tower.
    fn mul(&self, other: &Self) -> Self {
        let v0 = self.c0 * other.c0;
        let v1 = self.c1 * other.c1;

        let c0 = v0 + v1.mul_by_v();
        let c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - v0 - v1;

        Self { c0, c1 }
    }

    fn neg(&self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
        }
    }

    fn inv(&self) -> Option<Self> {
        let c0_sq = self.c0.square();
        let c1_sq_v = self.c1.square().mul_by_v();
        let denom = c0_sq - c1_sq_v;
        let denom_inv = denom.inv()?;

        Some(Self {
            c0: self.c0 * denom_inv,
            c1: (-self.c1) * denom_inv,
        })
    }

    fn square(&self) -> Self {
        let ab = self.c0 * self.c1;
        let c0 = (self.c0 + self.c1) * (self.c0 + self.c1.mul_by_v()) - ab - ab.mul_by_v();
        let c1 = ab.double();

        Self { c0, c1 }
    }

    fn pow(&self, exp: &[u64]) -> Self {
        let mut result = Self::ONE;
        let mut base = *self;
        for &limb in exp {
            let mut e = limb;
            for _ in 0..64 {
                if e & 1 == 1 {
                    result *= base;
                }
                base = base.square();
                e >>= 1;
            }
        }
        result
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn random(rng: &mut impl RngCore) -> Self {
        Self {
            c0: Fp6::random(rng),
            c1: Fp6::random(rng),
        }
    }

    fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }
}

impl fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp12({:?} + {:?}*w)", self.c0, self.c1)
    }
}

impl fmt::Display for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}*w)", self.c0, self.c1)
    }
}

impl Add for Fp12 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldElement::add(&self, &rhs)
    }
}

impl Sub for Fp12 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldElement::sub(&self, &rhs)
    }
}

impl Mul for Fp12 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldElement::mul(&self, &rhs)
    }
}

impl Neg for Fp12 {
    type Output = Self;
    fn neg(self) -> Self {
        FieldElement::neg(&self)
    }
}

impl AddAssign for Fp12 {
    fn add_assign(&mut self, rhs: Self) {
        *self = FieldElement::add(self, &rhs);
    }
}

impl SubAssign for Fp12 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = FieldElement::sub(self, &rhs);
    }
}

impl MulAssign for Fp12 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = FieldElement::mul(self, &rhs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_and_one() {
        assert!(Fp12::ZERO.is_zero());
        assert!(!Fp12::ONE.is_zero());
    }

    #[test]
    fn test_mul_identity() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        assert_eq!(a * Fp12::ONE, a);
    }

    #[test]
    fn test_mul_zero() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        assert_eq!(a * Fp12::ZERO, Fp12::ZERO);
    }

    #[test]
    fn test_squaring_matches_mul() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_inversion() {
        let mut rng = rand::thread_rng();
        for _ in 0..3 {
            let a = Fp12::random(&mut rng);
            if !a.is_zero() {
                let a_inv = a.inv().unwrap();
                assert_eq!(a * a_inv, Fp12::ONE);
            }
        }
    }

    #[test]
    fn test_negation() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        assert_eq!(a + (-a), Fp12::ZERO);
    }

    #[test]
    fn test_conjugate() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        let conj = a.conjugate();
        assert_eq!(conj.c0, a.c0);
        assert_eq!(conj.c1, -a.c1);
    }

    #[test]
    fn test_commutativity() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        let b = Fp12::random(&mut rng);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_associativity() {
        let mut rng = rand::thread_rng();
        let a = Fp12::random(&mut rng);
        let b = Fp12::random(&mut rng);
        let c = Fp12::random(&mut rng);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_w_squared_is_v() {
        let w = Fp12::new(Fp6::ZERO, Fp6::ONE);
        let w_sq = w.square();
        let v = Fp12::new(Fp6::new(Fp2::ZERO, Fp2::ONE, Fp2::ZERO), Fp6::ZERO);
        assert_eq!(w_sq, v);
    }
}
