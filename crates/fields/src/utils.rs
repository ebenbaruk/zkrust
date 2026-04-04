use crate::FieldElement;

/// Batch inversion using Montgomery's trick.
///
/// Given `[a₀, a₁, ..., aₙ₋₁]`, computes `[a₀⁻¹, a₁⁻¹, ..., aₙ₋₁⁻¹]`
/// using only a single field inversion and 3(n-1) field multiplications.
///
/// This is a massive optimization when many inversions are needed (e.g.,
/// converting projective points to affine, computing Lagrange coefficients).
///
/// Panics if any element is zero.
pub fn batch_inversion<F: FieldElement>(elements: &[F]) -> Vec<F> {
    if elements.is_empty() {
        return vec![];
    }

    let n = elements.len();

    // Phase 1: Compute prefix products
    // prefix[i] = a₀ * a₁ * ... * aᵢ
    let mut prefix = Vec::with_capacity(n);
    prefix.push(elements[0]);
    for i in 1..n {
        prefix.push(prefix[i - 1].mul(&elements[i]));
    }

    // Phase 2: Invert the total product
    let mut inv = prefix[n - 1]
        .inv()
        .expect("batch_inversion: encountered zero element");

    // Phase 3: Unwind to get individual inverses
    let mut result = vec![F::ZERO; n];
    for i in (1..n).rev() {
        // result[i] = inv * prefix[i-1] = (a₀...aₙ₋₁)⁻¹ * (a₀...aᵢ₋₁) = aᵢ⁻¹ * ...
        result[i] = inv.mul(&prefix[i - 1]);
        inv = inv.mul(&elements[i]);
    }
    result[0] = inv;

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fp;

    #[test]
    fn test_batch_inversion_empty() {
        let result: Vec<Fp> = batch_inversion(&[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_batch_inversion_single() {
        let a = Fp::from(42u64);
        let result = batch_inversion(&[a]);
        assert_eq!(result[0], a.inv().unwrap());
    }

    #[test]
    fn test_batch_inversion_matches_individual() {
        let elements: Vec<Fp> = (1..=10).map(|i| Fp::from(i as u64)).collect();
        let batch_result = batch_inversion(&elements);

        for (elem, inv) in elements.iter().zip(batch_result.iter()) {
            assert_eq!(*elem * *inv, Fp::ONE);
            assert_eq!(*inv, elem.inv().unwrap());
        }
    }

    #[test]
    fn test_batch_inversion_random() {
        use crate::Fr;
        let mut rng = rand::thread_rng();

        let elements: Vec<Fr> = (0..20)
            .map(|_| {
                let mut elem = Fr::random(&mut rng);
                // Ensure non-zero
                while elem.is_zero() {
                    elem = Fr::random(&mut rng);
                }
                elem
            })
            .collect();

        let batch_result = batch_inversion(&elements);

        for (elem, inv) in elements.iter().zip(batch_result.iter()) {
            assert_eq!(*elem * *inv, Fr::ONE);
        }
    }

    #[test]
    #[should_panic(expected = "zero")]
    fn test_batch_inversion_panics_on_zero() {
        let elements = vec![Fp::from(1u64), Fp::ZERO, Fp::from(3u64)];
        batch_inversion(&elements);
    }
}
