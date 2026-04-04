use crate::fp::Fp;
use crate::FieldElement;
use rand::RngCore;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Element of the quadratic extension field Fp2 = Fp[u] / (u² + 1).
///
/// Represented as c0 + c1 * u where c0, c1 ∈ Fp.
/// This field is used for G2 curve points on BN254.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fp2 {
    pub c0: Fp,
    pub c1: Fp,
}

impl Fp2 {
    pub const fn new(c0: Fp, c1: Fp) -> Self {
        Self { c0, c1 }
    }

    /// Multiply by the non-residue u (the imaginary unit).
    /// u * (c0 + c1*u) = c0*u + c1*u² = -c1 + c0*u
    pub fn mul_by_nonresidue(&self) -> Self {
        Self {
            c0: -self.c1,
            c1: self.c0,
        }
    }

    /// Compute the norm: Norm(a) = a * conjugate(a) = c0² + c1² ∈ Fp.
    pub fn norm(&self) -> Fp {
        self.c0.square() + self.c1.square()
    }

    /// Complex conjugate: conjugate(c0 + c1*u) = c0 - c1*u.
    pub fn conjugate(&self) -> Self {
        Self {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    /// Frobenius map: x -> x^p.
    /// For Fp2, this is just conjugation since p ≡ 3 (mod 4).
    pub fn frobenius_map(&self) -> Self {
        self.conjugate()
    }

    /// Multiply this Fp2 element by an Fp scalar.
    pub fn scale(&self, s: &Fp) -> Self {
        Self {
            c0: self.c0 * *s,
            c1: self.c1 * *s,
        }
    }

    /// Square root in Fp2.
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Fp2::ZERO);
        }

        let norm = self.norm();
        let norm_sqrt = norm.sqrt()?;

        let two_inv = Fp::from(2u64).inv()?;
        let candidate_c0_sq = (self.c0 + norm_sqrt) * two_inv;
        if let Some(c0) = candidate_c0_sq.sqrt() {
            let c0_inv = c0.inv()?;
            let c1 = self.c1 * two_inv * c0_inv;
            let result = Fp2::new(c0, c1);
            if result.square() == *self {
                return Some(result);
            }
        }

        let candidate_c0_sq = (self.c0 - norm_sqrt) * two_inv;
        let c0 = candidate_c0_sq.sqrt()?;
        let c0_inv = c0.inv()?;
        let c1 = self.c1 * two_inv * c0_inv;
        let result = Fp2::new(c0, c1);
        if result.square() == *self {
            Some(result)
        } else {
            None
        }
    }
}

impl FieldElement for Fp2 {
    const ZERO: Self = Fp2 {
        c0: Fp::ZERO,
        c1: Fp::ZERO,
    };
    const ONE: Self = Fp2 {
        c0: Fp::ONE,
        c1: Fp::ZERO,
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

    /// Karatsuba multiplication for Fp2:
    /// (a0 + a1*u)(b0 + b1*u) = (a0*b0 - a1*b1) + (a0*b1 + a1*b0)*u
    /// Using Karatsuba: 3 Fp multiplications instead of 4.
    fn mul(&self, other: &Self) -> Self {
        let v0 = self.c0 * other.c0;
        let v1 = self.c1 * other.c1;

        let c0 = v0 - v1;
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
        let norm = self.norm();
        let norm_inv = norm.inv()?;
        Some(Self {
            c0: self.c0 * norm_inv,
            c1: (-self.c1) * norm_inv,
        })
    }

    fn square(&self) -> Self {
        let ab = self.c0 * self.c1;
        let c0 = (self.c0 + self.c1) * (self.c0 - self.c1);
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
            c0: Fp::random(rng),
            c1: Fp::random(rng),
        }
    }

    fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }
}

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp2({:?} + {:?}*u)", self.c0, self.c1)
    }
}

impl fmt::Display for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}*u)", self.c0, self.c1)
    }
}

impl Add for Fp2 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldElement::add(&self, &rhs)
    }
}

impl Sub for Fp2 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldElement::sub(&self, &rhs)
    }
}

impl Mul for Fp2 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldElement::mul(&self, &rhs)
    }
}

impl Neg for Fp2 {
    type Output = Self;
    fn neg(self) -> Self {
        FieldElement::neg(&self)
    }
}

impl AddAssign for Fp2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = FieldElement::add(self, &rhs);
    }
}

impl SubAssign for Fp2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = FieldElement::sub(self, &rhs);
    }
}

impl MulAssign for Fp2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = FieldElement::mul(self, &rhs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make(a: u64, b: u64) -> Fp2 {
        Fp2::new(Fp::from(a), Fp::from(b))
    }

    #[test]
    fn test_zero_and_one() {
        assert!(Fp2::ZERO.is_zero());
        assert!(!Fp2::ONE.is_zero());
    }

    #[test]
    fn test_addition() {
        let a = make(3, 4);
        let b = make(5, 7);
        let c = a + b;
        assert_eq!(c, make(8, 11));
    }

    #[test]
    fn test_subtraction() {
        let a = make(10, 20);
        let b = make(3, 7);
        let c = a - b;
        assert_eq!(c, make(7, 13));
    }

    #[test]
    fn test_multiplication() {
        let a = make(3, 4);
        let b = make(5, 7);
        let c = a * b;
        let expected = Fp2::new(
            Fp::from(3u64) * Fp::from(5u64) - Fp::from(4u64) * Fp::from(7u64),
            Fp::from(3u64) * Fp::from(7u64) + Fp::from(4u64) * Fp::from(5u64),
        );
        assert_eq!(c, expected);
    }

    #[test]
    fn test_u_squared_is_minus_one() {
        let u = Fp2::new(Fp::ZERO, Fp::ONE);
        let u_sq = u.square();
        let minus_one = Fp2::new(-Fp::ONE, Fp::ZERO);
        assert_eq!(u_sq, minus_one);
    }

    #[test]
    fn test_negation() {
        let a = make(3, 4);
        assert_eq!(a + (-a), Fp2::ZERO);
    }

    #[test]
    fn test_inversion() {
        let a = make(3, 4);
        let a_inv = a.inv().unwrap();
        assert_eq!(a * a_inv, Fp2::ONE);
    }

    #[test]
    fn test_squaring_matches_mul() {
        let a = make(7, 11);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_conjugate() {
        let a = make(3, 4);
        let conj = a.conjugate();
        assert_eq!(conj.c0, a.c0);
        assert_eq!(conj.c1, -a.c1);
    }

    #[test]
    fn test_frobenius_is_conjugation() {
        let a = make(3, 4);
        assert_eq!(a.frobenius_map(), a.conjugate());
    }

    #[test]
    fn test_norm() {
        let a = make(3, 4);
        assert_eq!(a.norm(), Fp::from(25u64));
    }

    #[test]
    fn test_random_properties() {
        let mut rng = rand::thread_rng();
        for _ in 0..5 {
            let a = Fp2::random(&mut rng);
            let b = Fp2::random(&mut rng);

            assert_eq!(a * b, b * a);
            if !a.is_zero() {
                assert_eq!(a * a.inv().unwrap(), Fp2::ONE);
            }
        }
    }
}
