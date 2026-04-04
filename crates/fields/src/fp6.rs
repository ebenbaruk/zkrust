use crate::fp::Fp;
use crate::fp2::Fp2;
use crate::FieldElement;
use rand::RngCore;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Element of the sextic extension field Fp6 = Fp2[v] / (v³ - ξ),
/// where ξ = 9 + u ∈ Fp2 is the cubic non-residue.
///
/// Represented as c0 + c1*v + c2*v² where c0, c1, c2 ∈ Fp2.
/// This is an intermediate tower used to build Fp12.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Fp6 {
    pub c0: Fp2,
    pub c1: Fp2,
    pub c2: Fp2,
}

/// The non-residue ξ = 9 + u ∈ Fp2 used in the tower construction.
/// v³ = ξ in Fp6.
pub fn xi() -> Fp2 {
    Fp2::new(Fp::from(9u64), Fp::from(1u64))
}

/// Multiply an Fp2 element by ξ = 9 + u.
pub fn mul_by_xi(a: &Fp2) -> Fp2 {
    let nine = Fp::from(9u64);
    Fp2::new(a.c0 * nine - a.c1, a.c0 + a.c1 * nine)
}

impl Fp6 {
    pub const fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Self {
        Self { c0, c1, c2 }
    }

    /// Multiply by the element v (shifts coefficients).
    /// v * (c0 + c1*v + c2*v²) = c0*v + c1*v² + c2*v³ = c2*ξ + c0*v + c1*v²
    pub fn mul_by_v(&self) -> Self {
        Self {
            c0: mul_by_xi(&self.c2),
            c1: self.c0,
            c2: self.c1,
        }
    }

    /// Multiply by an Fp2 element (scaling all coefficients).
    pub fn scale(&self, s: &Fp2) -> Self {
        Self {
            c0: self.c0 * *s,
            c1: self.c1 * *s,
            c2: self.c2 * *s,
        }
    }

    /// Frobenius map: x -> x^p.
    pub fn frobenius_map(&self) -> Self {
        let c0 = self.c0.frobenius_map();
        let c1 = self.c1.frobenius_map();
        let c2 = self.c2.frobenius_map();

        // Multiply c1 by γ₁ = ξ^((p-1)/3) and c2 by γ₂ = ξ^(2(p-1)/3)
        let gamma_1 = Fp2::new(
            Fp::from_raw([
                0x99E39557176F553D,
                0xB78CC310C2D3E10E,
                0xA6B0BD4AFC1B20F5,
                0x2615F9B0F50A0C0E,
            ]),
            Fp::from_raw([
                0x4F49FFFD8BFD0000,
                0xDC39573399E34607,
                0x1D1ADDE6DC65AE80,
                0x09EE3E9BC7B6C80C,
            ]),
        );
        let gamma_2 = Fp2::new(
            Fp::from_raw([
                0xCEE3DEBFAB2D7063,
                0x8C4944C69A2B3B39,
                0xDDB5D0F87CE6D3A0,
                0x2324B5BF2B1A6B5A,
            ]),
            Fp::from_raw([
                0x94C27CCBED40FE5D,
                0xB76C59614E0644CF,
                0xCD3B1E96BF3E0B6F,
                0x05C7ECDDAC1A4B7A,
            ]),
        );

        Self {
            c0,
            c1: c1 * gamma_1,
            c2: c2 * gamma_2,
        }
    }
}

impl FieldElement for Fp6 {
    const ZERO: Self = Fp6 {
        c0: Fp2::ZERO,
        c1: Fp2::ZERO,
        c2: Fp2::ZERO,
    };
    const ONE: Self = Fp6 {
        c0: Fp2::ONE,
        c1: Fp2::ZERO,
        c2: Fp2::ZERO,
    };

    fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
            c2: self.c2 + other.c2,
        }
    }

    fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
            c2: self.c2 - other.c2,
        }
    }

    /// Multiplication in Fp6 using Karatsuba-like method (6 Fp2 muls).
    fn mul(&self, other: &Self) -> Self {
        let a0b0 = self.c0 * other.c0;
        let a1b1 = self.c1 * other.c1;
        let a2b2 = self.c2 * other.c2;

        let t0 = (self.c1 + self.c2) * (other.c1 + other.c2) - a1b1 - a2b2;
        let c0 = a0b0 + mul_by_xi(&t0);

        let t1 = (self.c0 + self.c1) * (other.c0 + other.c1) - a0b0 - a1b1;
        let c1 = t1 + mul_by_xi(&a2b2);

        let c2 = (self.c0 + self.c2) * (other.c0 + other.c2) - a0b0 - a2b2 + a1b1;

        Self { c0, c1, c2 }
    }

    fn neg(&self) -> Self {
        Self {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }

    fn inv(&self) -> Option<Self> {
        let c0_sq = self.c0.square();
        let c1_sq = self.c1.square();
        let c2_sq = self.c2.square();

        let c0c1 = self.c0 * self.c1;
        let c0c2 = self.c0 * self.c2;
        let c1c2 = self.c1 * self.c2;

        let t0 = c0_sq - mul_by_xi(&c1c2);
        let t1 = mul_by_xi(&c2_sq) - c0c1;
        let t2 = c1_sq - c0c2;

        let factor = self.c0 * t0 + mul_by_xi(&(self.c2 * t1 + self.c1 * t2));
        let factor_inv = factor.inv()?;

        Some(Self {
            c0: t0 * factor_inv,
            c1: t1 * factor_inv,
            c2: t2 * factor_inv,
        })
    }

    fn square(&self) -> Self {
        // Chung-Hasan SQR3 squaring
        let s0 = self.c0.square();
        let ab = self.c0 * self.c1;
        let s1 = ab.double();
        let s2 = (self.c0 - self.c1 + self.c2).square();
        let bc = self.c1 * self.c2;
        let s3 = bc.double();
        let s4 = self.c2.square();

        let c0 = s0 + mul_by_xi(&s3);
        let c1 = s1 + mul_by_xi(&s4);
        let c2 = s1 + s2 + s3 - s0 - s4;

        Self { c0, c1, c2 }
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
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    fn random(rng: &mut impl RngCore) -> Self {
        Self {
            c0: Fp2::random(rng),
            c1: Fp2::random(rng),
            c2: Fp2::random(rng),
        }
    }

    fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
            c2: self.c2.double(),
        }
    }
}

impl fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Fp6({:?} + {:?}*v + {:?}*v²)",
            self.c0, self.c1, self.c2
        )
    }
}

impl fmt::Display for Fp6 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}*v + {}*v²)", self.c0, self.c1, self.c2)
    }
}

impl Add for Fp6 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        FieldElement::add(&self, &rhs)
    }
}

impl Sub for Fp6 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        FieldElement::sub(&self, &rhs)
    }
}

impl Mul for Fp6 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        FieldElement::mul(&self, &rhs)
    }
}

impl Neg for Fp6 {
    type Output = Self;
    fn neg(self) -> Self {
        FieldElement::neg(&self)
    }
}

impl AddAssign for Fp6 {
    fn add_assign(&mut self, rhs: Self) {
        *self = FieldElement::add(self, &rhs);
    }
}

impl SubAssign for Fp6 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = FieldElement::sub(self, &rhs);
    }
}

impl MulAssign for Fp6 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = FieldElement::mul(self, &rhs);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero_and_one() {
        assert!(Fp6::ZERO.is_zero());
        assert!(!Fp6::ONE.is_zero());
    }

    #[test]
    fn test_mul_identity() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        assert_eq!(a * Fp6::ONE, a);
        assert_eq!(Fp6::ONE * a, a);
    }

    #[test]
    fn test_mul_zero() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        assert_eq!(a * Fp6::ZERO, Fp6::ZERO);
    }

    #[test]
    fn test_squaring_matches_mul() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        assert_eq!(a.square(), a * a);
    }

    #[test]
    fn test_inversion() {
        let mut rng = rand::thread_rng();
        for _ in 0..5 {
            let a = Fp6::random(&mut rng);
            if !a.is_zero() {
                let a_inv = a.inv().unwrap();
                assert_eq!(a * a_inv, Fp6::ONE);
            }
        }
    }

    #[test]
    fn test_negation() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        assert_eq!(a + (-a), Fp6::ZERO);
    }

    #[test]
    fn test_distributivity() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        let b = Fp6::random(&mut rng);
        let c = Fp6::random(&mut rng);
        assert_eq!(a * (b + c), a * b + a * c);
    }

    #[test]
    fn test_associativity() {
        let mut rng = rand::thread_rng();
        let a = Fp6::random(&mut rng);
        let b = Fp6::random(&mut rng);
        let c = Fp6::random(&mut rng);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_v_cubed_is_xi() {
        let v = Fp6::new(Fp2::ZERO, Fp2::ONE, Fp2::ZERO);
        let v_cubed = v * v * v;
        let expected = Fp6::new(xi(), Fp2::ZERO, Fp2::ZERO);
        assert_eq!(v_cubed, expected);
    }
}
