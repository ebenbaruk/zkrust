use std::ops::{Add, Mul, Neg, Sub};
use zkrust_fields::{FieldElement, Fp, Fp2};

/// A point on the BN254 G2 twist curve in affine coordinates.
/// Curve equation: y² = x³ + b' where b' = 3/(9+u) over Fp2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct G2Affine {
    pub x: Fp2,
    pub y: Fp2,
    pub infinity: bool,
}

/// A point on the BN254 G2 twist curve in Jacobian projective coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct G2Projective {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
}

/// b' = 3/(9+u) = 3 * (9-u)/((9+u)(9-u)) = 3*(9-u)/(81+1) = 3*(9-u)/82
/// In BN254, the twist coefficient is b' = 3/(9+u) computed in Fp2.
fn twist_b() -> Fp2 {
    let xi = Fp2::new(Fp::from(9u64), Fp::from(1u64));
    let three = Fp2::new(Fp::from(3u64), Fp::ZERO);
    let xi_inv = xi.inv().unwrap();
    three * xi_inv
}

impl G2Affine {
    pub const fn identity() -> Self {
        Self {
            x: Fp2::ZERO,
            y: Fp2::ZERO,
            infinity: true,
        }
    }

    /// The standard BN254 G2 generator.
    pub fn generator() -> Self {
        Self {
            x: Fp2::new(
                Fp::from_raw([
                    0x46DEBD5CD992F6ED,
                    0x674322D4F75EDADD,
                    0x426A00665E5C4479,
                    0x1800DEEF121F1E76,
                ]),
                Fp::from_raw([
                    0x97E485B7AEF312C2,
                    0xF1AA493335A9E712,
                    0x7260BFB731FB5D25,
                    0x198E9393920D483A,
                ]),
            ),
            y: Fp2::new(
                Fp::from_raw([
                    0x4CE6CC0166FA7DAA,
                    0xE3D1E7690C43D37B,
                    0x4AAB71808DCB408F,
                    0x12C85EA5DB8C6DEB,
                ]),
                Fp::from_raw([
                    0x55ACDADCD122975B,
                    0xBC4B313370B38EF3,
                    0xEC9E99AD690C3395,
                    0x090689D0585FF075,
                ]),
            ),
            infinity: false,
        }
    }

    /// Check if this point lies on the G2 twist curve.
    pub fn is_on_curve(&self) -> bool {
        if self.infinity {
            return true;
        }
        let y2 = self.y.square();
        let x3_b = self.x.square() * self.x + twist_b();
        y2 == x3_b
    }

    pub fn is_identity(&self) -> bool {
        self.infinity
    }

    pub fn to_projective(&self) -> G2Projective {
        if self.infinity {
            G2Projective::identity()
        } else {
            G2Projective {
                x: self.x,
                y: self.y,
                z: Fp2::ONE,
            }
        }
    }
}

impl G2Projective {
    pub fn identity() -> Self {
        Self {
            x: Fp2::ZERO,
            y: Fp2::ONE,
            z: Fp2::ZERO,
        }
    }

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    pub fn to_affine(&self) -> G2Affine {
        if self.is_identity() {
            return G2Affine::identity();
        }

        let z_inv = self.z.inv().unwrap();
        let z_inv2 = z_inv.square();
        let z_inv3 = z_inv2 * z_inv;

        G2Affine {
            x: self.x * z_inv2,
            y: self.y * z_inv3,
            infinity: false,
        }
    }

    pub fn double(&self) -> Self {
        if self.is_identity() {
            return *self;
        }

        let a = self.x.square();
        let b = self.y.square();
        let c = b.square();

        let d = ((self.x + b).square() - a - c).double();
        let e = a.double() + a;
        let f = e.square();

        let x3 = f - d.double();
        let y3 = e * (d - x3) - c.double().double().double();
        let z3 = (self.y * self.z).double();

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn add_mixed(&self, rhs: &G2Affine) -> Self {
        if rhs.is_identity() {
            return *self;
        }
        if self.is_identity() {
            return rhs.to_projective();
        }

        let z1z1 = self.z.square();
        let u2 = rhs.x * z1z1;
        let s2 = rhs.y * self.z * z1z1;

        let h = u2 - self.x;
        let hh = h.square();
        let i = hh.double().double();
        let j = h * i;
        let r = (s2 - self.y).double();

        if h.is_zero() && r.is_zero() {
            return self.double();
        }

        let v = self.x * i;
        let x3 = r.square() - j - v.double();
        let y3 = r * (v - x3) - (self.y * j).double();
        let z3 = (self.z + h).square() - z1z1 - hh;

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn add_projective(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return *rhs;
        }
        if rhs.is_identity() {
            return *self;
        }

        let z1z1 = self.z.square();
        let z2z2 = rhs.z.square();
        let u1 = self.x * z2z2;
        let u2 = rhs.x * z1z1;
        let s1 = self.y * rhs.z * z2z2;
        let s2 = rhs.y * self.z * z1z1;

        let h = u2 - u1;
        let r = (s2 - s1).double();

        if h.is_zero() {
            if r.is_zero() {
                return self.double();
            }
            return Self::identity();
        }

        let i = h.double().square();
        let j = h * i;
        let v = u1 * i;
        let x3 = r.square() - j - v.double();
        let y3 = r * (v - x3) - (s1 * j).double();
        let z3 = ((self.z + rhs.z).square() - z1z1 - z2z2) * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    pub fn scalar_mul(&self, scalar: &[u64; 4]) -> Self {
        let mut result = Self::identity();
        let mut base = *self;

        for &limb in scalar.iter() {
            let mut s = limb;
            for _ in 0..64 {
                if s & 1 == 1 {
                    result = result.add_projective(&base);
                }
                base = base.double();
                s >>= 1;
            }
        }

        result
    }
}

impl Add for G2Projective {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.add_projective(&rhs)
    }
}

impl Sub for G2Projective {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.add_projective(&(-rhs))
    }
}

impl Neg for G2Projective {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_identity() {
            self
        } else {
            Self {
                x: self.x,
                y: -self.y,
                z: self.z,
            }
        }
    }
}

impl Neg for G2Affine {
    type Output = Self;
    fn neg(self) -> Self {
        if self.infinity {
            self
        } else {
            Self {
                x: self.x,
                y: -self.y,
                infinity: false,
            }
        }
    }
}

impl Add<G2Affine> for G2Projective {
    type Output = G2Projective;
    fn add(self, rhs: G2Affine) -> G2Projective {
        self.add_mixed(&rhs)
    }
}

impl Mul<zkrust_fields::Fr> for G2Projective {
    type Output = G2Projective;
    fn mul(self, rhs: zkrust_fields::Fr) -> G2Projective {
        self.scalar_mul(&rhs.to_raw())
    }
}

impl Mul<zkrust_fields::Fr> for G2Affine {
    type Output = G2Projective;
    fn mul(self, rhs: zkrust_fields::Fr) -> G2Projective {
        self.to_projective().scalar_mul(&rhs.to_raw())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_fields::Fr;

    #[test]
    fn test_generator_on_curve() {
        let g = G2Affine::generator();
        assert!(g.is_on_curve(), "G2 generator should be on curve");
    }

    #[test]
    fn test_identity_on_curve() {
        assert!(G2Affine::identity().is_on_curve());
    }

    #[test]
    fn test_add_identity() {
        let g = G2Affine::generator().to_projective();
        let id = G2Projective::identity();
        assert_eq!((g + id).to_affine(), g.to_affine());
        assert_eq!((id + g).to_affine(), g.to_affine());
    }

    #[test]
    fn test_double() {
        let g = G2Affine::generator().to_projective();
        let g2 = g.double();
        assert!(g2.to_affine().is_on_curve());
    }

    #[test]
    fn test_add_inverse_is_identity() {
        let g = G2Affine::generator().to_projective();
        let result = g + (-g);
        assert!(result.is_identity());
    }

    #[test]
    fn test_double_equals_add_self() {
        let g = G2Affine::generator().to_projective();
        assert_eq!(g.double().to_affine(), (g + g).to_affine());
    }

    #[test]
    fn test_order_of_generator() {
        let g = G2Affine::generator().to_projective();
        let r = zkrust_fields::fr::MODULUS;
        let result = g.scalar_mul(&r);
        assert!(result.is_identity());
    }

    #[test]
    fn test_scalar_mul_distributive() {
        let g = G2Affine::generator().to_projective();
        let a = Fr::from(77u64);
        let b = Fr::from(33u64);
        let lhs = g * (a + b);
        let rhs = (g * a) + (g * b);
        assert_eq!(lhs.to_affine(), rhs.to_affine());
    }
}
