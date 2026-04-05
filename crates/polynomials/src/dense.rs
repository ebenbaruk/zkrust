use std::ops::{Add, Mul, Neg, Sub};
use zkrust_fields::{FieldElement, Fr};

/// A polynomial represented by its coefficients in ascending degree order.
///
/// `coeffs[i]` is the coefficient of x^i.
/// Empty coeffs represents the zero polynomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DensePolynomial {
    pub coeffs: Vec<Fr>,
}

impl DensePolynomial {
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    pub fn new(coeffs: Vec<Fr>) -> Self {
        let mut p = Self { coeffs };
        p.truncate();
        p
    }

    /// Constant polynomial.
    pub fn constant(c: Fr) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self { coeffs: vec![c] }
        }
    }

    /// The polynomial coeff * x^n.
    pub fn monomial(n: usize, coeff: Fr) -> Self {
        if coeff.is_zero() {
            return Self::zero();
        }
        let mut coeffs = vec![Fr::ZERO; n + 1];
        coeffs[n] = coeff;
        Self { coeffs }
    }

    /// Construct from evaluation points using Lagrange interpolation.
    pub fn from_evaluations(points: &[(Fr, Fr)]) -> Self {
        if points.is_empty() {
            return Self::zero();
        }

        let n = points.len();
        let mut result = Self::zero();

        for i in 0..n {
            let (xi, yi) = points[i];
            if yi.is_zero() {
                continue;
            }

            let mut basis = Self::constant(Fr::ONE);
            let mut denom = Fr::ONE;

            for (j, &(xj, _)) in points.iter().enumerate() {
                if i == j {
                    continue;
                }
                basis = basis * Self::new(vec![-xj, Fr::ONE]);
                denom *= xi - xj;
            }

            let denom_inv = denom.inv().unwrap();
            basis = basis.scale(&(yi * denom_inv));
            result = result + basis;
        }

        result
    }

    /// The vanishing polynomial Z_H(x) = x^n - 1 for a domain of size n.
    pub fn vanishing(n: usize) -> Self {
        let mut coeffs = vec![Fr::ZERO; n + 1];
        coeffs[0] = -Fr::ONE;
        coeffs[n] = Fr::ONE;
        Self { coeffs }
    }

    /// Degree of the polynomial (-1 for zero polynomial, represented as None).
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// Evaluate polynomial at a point using Horner's method.
    pub fn evaluate(&self, x: &Fr) -> Fr {
        if self.coeffs.is_empty() {
            return Fr::ZERO;
        }
        let mut result = Fr::ZERO;
        for c in self.coeffs.iter().rev() {
            result = result * *x + *c;
        }
        result
    }

    /// Scale all coefficients by a scalar.
    pub fn scale(&self, s: &Fr) -> Self {
        if s.is_zero() {
            return Self::zero();
        }
        Self::new(self.coeffs.iter().map(|c| *c * *s).collect())
    }

    /// Naive O(n*m) polynomial multiplication.
    pub fn naive_mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut result = vec![Fr::ZERO; n];
        for (i, a) in self.coeffs.iter().enumerate() {
            for (j, b) in other.coeffs.iter().enumerate() {
                result[i + j] += *a * *b;
            }
        }
        Self::new(result)
    }

    /// Polynomial long division: returns (quotient, remainder) such that
    /// self = quotient * divisor + remainder.
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero(), "division by zero polynomial");

        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }

        let self_deg = self.degree().unwrap();
        let div_deg = divisor.degree().unwrap();

        if self_deg < div_deg {
            return (Self::zero(), self.clone());
        }

        let lead_inv = divisor.coeffs.last().unwrap().inv().unwrap();
        let mut remainder = self.coeffs.clone();
        let q_len = self_deg - div_deg + 1;
        let mut quotient = vec![Fr::ZERO; q_len];

        for i in (0..q_len).rev() {
            let coeff = remainder[i + div_deg] * lead_inv;
            quotient[i] = coeff;
            for (j, d) in divisor.coeffs.iter().enumerate() {
                remainder[i + j] -= coeff * *d;
            }
        }

        (Self::new(quotient), Self::new(remainder))
    }

    /// Exact division: panics if there is a nonzero remainder.
    pub fn div_exact(&self, divisor: &Self) -> Self {
        let (q, r) = self.div_rem(divisor);
        assert!(r.is_zero(), "polynomial division has nonzero remainder");
        q
    }

    /// Evaluate at all points in a domain {omega^0, omega^1, ..., omega^(n-1)}.
    /// This is just iterative evaluation (not FFT — see ntt module for FFT).
    pub fn evaluate_domain(&self, domain_size: usize) -> Vec<Fr> {
        let k = domain_size.trailing_zeros();
        let omega = Fr::root_of_unity(k);
        let mut results = Vec::with_capacity(domain_size);
        let mut x = Fr::ONE;
        for _ in 0..domain_size {
            results.push(self.evaluate(&x));
            x *= omega;
        }
        results
    }

    /// Remove trailing zero coefficients.
    fn truncate(&mut self) {
        while self.coeffs.last().is_some_and(|c| c.is_zero()) {
            self.coeffs.pop();
        }
    }

    /// Leading coefficient (None for zero polynomial).
    pub fn leading_coeff(&self) -> Option<Fr> {
        self.coeffs.last().copied()
    }
}

impl Add for DensePolynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        if self.is_zero() {
            return rhs;
        }
        if rhs.is_zero() {
            return self;
        }
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = vec![Fr::ZERO; n];
        for (i, c) in self.coeffs.iter().enumerate() {
            coeffs[i] += *c;
        }
        for (i, c) in rhs.coeffs.iter().enumerate() {
            coeffs[i] += *c;
        }
        Self::new(coeffs)
    }
}

impl Sub for DensePolynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if rhs.is_zero() {
            return self;
        }
        if self.is_zero() {
            return -rhs;
        }
        let n = self.coeffs.len().max(rhs.coeffs.len());
        let mut coeffs = vec![Fr::ZERO; n];
        for (i, c) in self.coeffs.iter().enumerate() {
            coeffs[i] += *c;
        }
        for (i, c) in rhs.coeffs.iter().enumerate() {
            coeffs[i] -= *c;
        }
        Self::new(coeffs)
    }
}

impl Mul for DensePolynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.naive_mul(&rhs)
    }
}

impl Neg for DensePolynomial {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(self.coeffs.iter().map(|c| -*c).collect())
    }
}

impl<'a> Add<&'a DensePolynomial> for &'a DensePolynomial {
    type Output = DensePolynomial;
    fn add(self, rhs: &'a DensePolynomial) -> DensePolynomial {
        self.clone() + rhs.clone()
    }
}

impl<'a> Sub<&'a DensePolynomial> for &'a DensePolynomial {
    type Output = DensePolynomial;
    fn sub(self, rhs: &'a DensePolynomial) -> DensePolynomial {
        self.clone() - rhs.clone()
    }
}

impl<'a> Mul<&'a DensePolynomial> for &'a DensePolynomial {
    type Output = DensePolynomial;
    fn mul(self, rhs: &'a DensePolynomial) -> DensePolynomial {
        self.naive_mul(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zero() {
        let z = DensePolynomial::zero();
        assert!(z.is_zero());
        assert_eq!(z.degree(), None);
        assert_eq!(z.evaluate(&Fr::from(5u64)), Fr::ZERO);
    }

    #[test]
    fn test_constant() {
        let c = DensePolynomial::constant(Fr::from(7u64));
        assert_eq!(c.degree(), Some(0));
        assert_eq!(c.evaluate(&Fr::from(100u64)), Fr::from(7u64));
    }

    #[test]
    fn test_evaluate_linear() {
        // p(x) = 2x + 3
        let p = DensePolynomial::new(vec![Fr::from(3u64), Fr::from(2u64)]);
        assert_eq!(p.evaluate(&Fr::from(5u64)), Fr::from(13u64));
    }

    #[test]
    fn test_add() {
        let a = DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]);
        let b = DensePolynomial::new(vec![Fr::from(3u64), Fr::from(4u64), Fr::from(5u64)]);
        let c = a + b;
        assert_eq!(
            c.coeffs,
            vec![Fr::from(4u64), Fr::from(6u64), Fr::from(5u64)]
        );
    }

    #[test]
    fn test_sub() {
        let a = DensePolynomial::new(vec![Fr::from(5u64), Fr::from(3u64)]);
        let b = DensePolynomial::new(vec![Fr::from(2u64), Fr::from(3u64)]);
        let c = a - b;
        assert_eq!(c.coeffs, vec![Fr::from(3u64)]);
    }

    #[test]
    fn test_mul() {
        // (1 + x) * (1 + x) = 1 + 2x + x^2
        let a = DensePolynomial::new(vec![Fr::ONE, Fr::ONE]);
        let b = a.clone();
        let c = a * b;
        assert_eq!(c.coeffs, vec![Fr::ONE, Fr::from(2u64), Fr::ONE]);
    }

    #[test]
    fn test_div_exact() {
        // (x^2 - 1) / (x - 1) = x + 1
        let num = DensePolynomial::new(vec![-Fr::ONE, Fr::ZERO, Fr::ONE]);
        let den = DensePolynomial::new(vec![-Fr::ONE, Fr::ONE]);
        let q = num.div_exact(&den);
        assert_eq!(q.coeffs, vec![Fr::ONE, Fr::ONE]);
    }

    #[test]
    fn test_div_rem() {
        // (x^2 + 1) / (x + 1) = (x - 1) remainder 2
        let num = DensePolynomial::new(vec![Fr::ONE, Fr::ZERO, Fr::ONE]);
        let den = DensePolynomial::new(vec![Fr::ONE, Fr::ONE]);
        let (q, r) = num.div_rem(&den);
        assert_eq!(q.coeffs, vec![-Fr::ONE, Fr::ONE]);
        assert_eq!(r.coeffs, vec![Fr::from(2u64)]);
    }

    #[test]
    fn test_lagrange_interpolation() {
        // Points: (0, 1), (1, 3), (2, 7) → p(x) = x^2 + x + 1
        let points = vec![
            (Fr::from(0u64), Fr::from(1u64)),
            (Fr::from(1u64), Fr::from(3u64)),
            (Fr::from(2u64), Fr::from(7u64)),
        ];
        let p = DensePolynomial::from_evaluations(&points);
        assert_eq!(p.degree(), Some(2));
        for (x, y) in &points {
            assert_eq!(p.evaluate(x), *y);
        }
    }

    #[test]
    fn test_vanishing_polynomial() {
        // Z_H(x) = x^4 - 1 should vanish at all 4th roots of unity
        let n = 4;
        let zh = DensePolynomial::vanishing(n);
        let omega = Fr::root_of_unity(2); // 4th root
        let mut x = Fr::ONE;
        for _ in 0..n {
            assert_eq!(zh.evaluate(&x), Fr::ZERO);
            x *= omega;
        }
    }

    #[test]
    fn test_monomial() {
        let p = DensePolynomial::monomial(3, Fr::from(5u64));
        assert_eq!(p.degree(), Some(3));
        assert_eq!(p.evaluate(&Fr::from(2u64)), Fr::from(40u64));
    }

    #[test]
    fn test_scale() {
        let p = DensePolynomial::new(vec![Fr::from(1u64), Fr::from(2u64)]);
        let s = p.scale(&Fr::from(3u64));
        assert_eq!(s.coeffs, vec![Fr::from(3u64), Fr::from(6u64)]);
    }

    #[test]
    fn test_div_then_mul() {
        let mut rng = rand::thread_rng();
        let a = DensePolynomial::new((0..5).map(|_| Fr::random(&mut rng)).collect());
        let b = DensePolynomial::new((0..3).map(|_| Fr::random(&mut rng)).collect());
        let product = a.clone() * b.clone();
        let recovered = product.div_exact(&b);
        assert_eq!(recovered, a);
    }
}
