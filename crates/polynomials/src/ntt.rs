use zkrust_fields::{FieldElement, Fr};

/// Compute the NTT (Number Theoretic Transform) of a polynomial.
///
/// Input: coefficients of length n (must be a power of 2).
/// Output: evaluations at {omega^0, omega^1, ..., omega^(n-1)}.
///
/// Uses the radix-2 Cooley-Tukey DIT (decimation-in-time) algorithm.
pub fn ntt(coeffs: &[Fr]) -> Vec<Fr> {
    let n = coeffs.len();
    assert!(n.is_power_of_two(), "NTT size must be a power of 2");
    if n == 1 {
        return coeffs.to_vec();
    }

    let k = n.trailing_zeros();
    let omega = Fr::root_of_unity(k);

    let mut a = coeffs.to_vec();
    bit_reverse_permutation(&mut a);
    butterfly(&mut a, omega);
    a
}

/// Compute the inverse NTT: evaluations → coefficients.
///
/// Applies NTT with omega^(-1), then scales by n^(-1).
pub fn intt(evals: &[Fr]) -> Vec<Fr> {
    let n = evals.len();
    assert!(n.is_power_of_two(), "INTT size must be a power of 2");
    if n == 1 {
        return evals.to_vec();
    }

    let k = n.trailing_zeros();
    let omega_inv = Fr::root_of_unity(k).inv().unwrap();

    let mut a = evals.to_vec();
    bit_reverse_permutation(&mut a);
    butterfly(&mut a, omega_inv);

    // Scale by 1/n
    let n_inv = Fr::from(n as u64).inv().unwrap();
    for x in a.iter_mut() {
        *x *= n_inv;
    }
    a
}

/// Coset NTT: evaluate polynomial on coset {g, g*omega, g*omega^2, ...}
/// where g is a generator shift.
///
/// This is equivalent to NTT of [c0, c1*g, c2*g^2, ...].
pub fn coset_ntt(coeffs: &[Fr], generator: Fr) -> Vec<Fr> {
    let n = coeffs.len();
    let mut shifted = Vec::with_capacity(n);
    let mut g_pow = Fr::ONE;
    for c in coeffs {
        shifted.push(*c * g_pow);
        g_pow *= generator;
    }
    ntt(&shifted)
}

/// Inverse coset NTT: recover coefficients from coset evaluations.
pub fn coset_intt(evals: &[Fr], generator: Fr) -> Vec<Fr> {
    let mut coeffs = intt(evals);
    let g_inv = generator.inv().unwrap();
    let mut g_pow = Fr::ONE;
    for c in coeffs.iter_mut() {
        *c *= g_pow;
        g_pow *= g_inv;
    }
    coeffs
}

/// Multiply two polynomials using NTT.
///
/// Pads both to the next power of 2 that fits the product,
/// transforms, pointwise multiplies, and inverse-transforms.
pub fn ntt_mul(a: &[Fr], b: &[Fr]) -> Vec<Fr> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let result_len = a.len() + b.len() - 1;
    let n = result_len.next_power_of_two();

    let mut a_padded = vec![Fr::ZERO; n];
    let mut b_padded = vec![Fr::ZERO; n];
    a_padded[..a.len()].copy_from_slice(a);
    b_padded[..b.len()].copy_from_slice(b);

    let a_evals = ntt(&a_padded);
    let b_evals = ntt(&b_padded);

    let mut c_evals = vec![Fr::ZERO; n];
    for i in 0..n {
        c_evals[i] = a_evals[i] * b_evals[i];
    }

    let mut result = intt(&c_evals);
    result.truncate(result_len);
    result
}

/// Bit-reverse permutation of the array.
fn bit_reverse_permutation(a: &mut [Fr]) {
    let n = a.len();
    let log_n = n.trailing_zeros();

    for i in 0..n {
        let j = reverse_bits(i as u32, log_n) as usize;
        if i < j {
            a.swap(i, j);
        }
    }
}

/// Reverse the low `bits` bits of `x`.
fn reverse_bits(x: u32, bits: u32) -> u32 {
    let mut result = 0u32;
    let mut val = x;
    for _ in 0..bits {
        result = (result << 1) | (val & 1);
        val >>= 1;
    }
    result
}

/// In-place butterfly computation for NTT.
fn butterfly(a: &mut [Fr], omega: Fr) {
    let n = a.len();
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let step = n / len;
        // Precompute twiddle factors
        let w_base = omega.pow(&[step as u64, 0, 0, 0]);

        for start in (0..n).step_by(len) {
            let mut w = Fr::ONE;
            for j in 0..half {
                let u = a[start + j];
                let v = a[start + j + half] * w;
                a[start + j] = u + v;
                a[start + j + half] = u - v;
                w *= w_base;
            }
        }
        len <<= 1;
    }
}

/// A multiplicative coset generator for Fr.
/// Uses the multiplicative generator of Fr (a generator of the full Fr* group).
/// For coset FFT in Groth16, we use g = 5 (a common choice for BN254).
pub fn coset_generator() -> Fr {
    Fr::from(5u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ntt_intt_roundtrip() {
        let coeffs: Vec<Fr> = (0..8).map(|i| Fr::from(i as u64)).collect();
        let evals = ntt(&coeffs);
        let recovered = intt(&evals);
        assert_eq!(recovered, coeffs);
    }

    #[test]
    fn test_ntt_intt_roundtrip_large() {
        let mut rng = rand::thread_rng();
        let n = 256;
        let coeffs: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
        let evals = ntt(&coeffs);
        let recovered = intt(&evals);
        assert_eq!(recovered, coeffs);
    }

    #[test]
    fn test_ntt_evaluates_correctly() {
        // NTT should produce evaluations at powers of omega
        let coeffs = vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
        ];
        let evals = ntt(&coeffs);

        let omega = Fr::root_of_unity(2); // 4th root of unity
        let mut x = Fr::ONE;
        for (i, eval) in evals.iter().enumerate() {
            let mut expected = Fr::ZERO;
            let mut x_pow = Fr::ONE;
            for c in &coeffs {
                expected = expected + *c * x_pow;
                x_pow = x_pow * x;
            }
            assert_eq!(*eval, expected, "mismatch at index {}", i);
            x = x * omega;
        }
    }

    #[test]
    fn test_ntt_mul_matches_naive() {
        let a = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let b = vec![Fr::from(4u64), Fr::from(5u64)];

        let ntt_result = ntt_mul(&a, &b);

        // Naive: (1 + 2x + 3x^2)(4 + 5x) = 4 + 13x + 22x^2 + 15x^3
        let expected = vec![
            Fr::from(4u64),
            Fr::from(13u64),
            Fr::from(22u64),
            Fr::from(15u64),
        ];
        assert_eq!(ntt_result, expected);
    }

    #[test]
    fn test_ntt_mul_random() {
        use crate::DensePolynomial;
        let mut rng = rand::thread_rng();

        let a_coeffs: Vec<Fr> = (0..32).map(|_| Fr::random(&mut rng)).collect();
        let b_coeffs: Vec<Fr> = (0..32).map(|_| Fr::random(&mut rng)).collect();

        let a_poly = DensePolynomial::new(a_coeffs.clone());
        let b_poly = DensePolynomial::new(b_coeffs.clone());
        let naive = a_poly.naive_mul(&b_poly);

        let fast = ntt_mul(&a_coeffs, &b_coeffs);
        let fast_poly = DensePolynomial::new(fast);

        assert_eq!(naive, fast_poly);
    }

    #[test]
    fn test_coset_ntt_roundtrip() {
        let mut rng = rand::thread_rng();
        let n = 16;
        let coeffs: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
        let g = coset_generator();
        let evals = coset_ntt(&coeffs, g);
        let recovered = coset_intt(&evals, g);
        assert_eq!(recovered, coeffs);
    }

    #[test]
    fn test_ntt_single_element() {
        let coeffs = vec![Fr::from(42u64)];
        let evals = ntt(&coeffs);
        assert_eq!(evals, coeffs);
        let recovered = intt(&evals);
        assert_eq!(recovered, coeffs);
    }
}
