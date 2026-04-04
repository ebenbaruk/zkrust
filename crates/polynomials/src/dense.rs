use zkrust_fields::Fr;

/// A polynomial represented by its coefficients in ascending degree order.
///
/// `coeffs[i]` is the coefficient of x^i.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DensePolynomial {
    pub coeffs: Vec<Fr>,
}

impl DensePolynomial {
    pub fn zero() -> Self {
        Self { coeffs: vec![] }
    }
}
