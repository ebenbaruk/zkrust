use crate::g1::{G1Affine, G1Projective};
use crate::g2::{G2Affine, G2Projective};
use zkrust_fields::Fr;

/// Multi-scalar multiplication for G1: compute sum(scalars[i] * bases[i]).
///
/// Uses the Pippenger bucket method for large inputs,
/// falling back to naive summation for small inputs.
pub fn msm_g1(bases: &[G1Affine], scalars: &[Fr]) -> G1Projective {
    assert_eq!(bases.len(), scalars.len());
    let n = bases.len();

    if n == 0 {
        return G1Projective::identity();
    }
    if n < 8 {
        return naive_msm_g1(bases, scalars);
    }

    pippenger_g1(bases, scalars)
}

/// Multi-scalar multiplication for G2.
pub fn msm_g2(bases: &[G2Affine], scalars: &[Fr]) -> G2Projective {
    assert_eq!(bases.len(), scalars.len());
    let n = bases.len();

    if n == 0 {
        return G2Projective::identity();
    }
    if n < 8 {
        return naive_msm_g2(bases, scalars);
    }

    pippenger_g2(bases, scalars)
}

fn naive_msm_g1(bases: &[G1Affine], scalars: &[Fr]) -> G1Projective {
    let mut result = G1Projective::identity();
    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        result = result + (*base * *scalar);
    }
    result
}

fn naive_msm_g2(bases: &[G2Affine], scalars: &[Fr]) -> G2Projective {
    let mut result = G2Projective::identity();
    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        result = result + (*base * *scalar);
    }
    result
}

/// Pippenger's bucket method for G1.
/// Window size chosen based on input size for optimal performance.
fn pippenger_g1(bases: &[G1Affine], scalars: &[Fr]) -> G1Projective {
    let n = bases.len();
    let window = optimal_window(n);
    let num_windows = 256_usize.div_ceil(window);
    let num_buckets = (1 << window) - 1;

    let scalar_limbs: Vec<[u64; 4]> = scalars.iter().map(|s| s.to_raw()).collect();

    let mut result = G1Projective::identity();

    for w in (0..num_windows).rev() {
        // Shift result by window bits
        for _ in 0..window {
            result = result.double();
        }

        // Fill buckets
        let mut buckets = vec![G1Projective::identity(); num_buckets];

        for (i, base) in bases.iter().enumerate() {
            let scalar_bits = get_window_bits(&scalar_limbs[i], w, window);
            if scalar_bits > 0 {
                buckets[scalar_bits - 1] = buckets[scalar_bits - 1] + *base;
            }
        }

        // Accumulate buckets: bucket[k] contributes (k+1) times
        let mut running_sum = G1Projective::identity();
        let mut window_sum = G1Projective::identity();
        for bucket in buckets.iter().rev() {
            running_sum = running_sum + *bucket;
            window_sum = window_sum + running_sum;
        }

        result = result + window_sum;
    }

    result
}

/// Pippenger's bucket method for G2.
fn pippenger_g2(bases: &[G2Affine], scalars: &[Fr]) -> G2Projective {
    let n = bases.len();
    let window = optimal_window(n);
    let num_windows = 256_usize.div_ceil(window);
    let num_buckets = (1 << window) - 1;

    let scalar_limbs: Vec<[u64; 4]> = scalars.iter().map(|s| s.to_raw()).collect();

    let mut result = G2Projective::identity();

    for w in (0..num_windows).rev() {
        for _ in 0..window {
            result = result.double();
        }

        let mut buckets = vec![G2Projective::identity(); num_buckets];

        for (i, base) in bases.iter().enumerate() {
            let scalar_bits = get_window_bits(&scalar_limbs[i], w, window);
            if scalar_bits > 0 {
                buckets[scalar_bits - 1] = buckets[scalar_bits - 1] + *base;
            }
        }

        let mut running_sum = G2Projective::identity();
        let mut window_sum = G2Projective::identity();
        for bucket in buckets.iter().rev() {
            running_sum = running_sum + *bucket;
            window_sum = window_sum + running_sum;
        }

        result = result + window_sum;
    }

    result
}

/// Extract a window of bits from a scalar at position `window_idx`.
fn get_window_bits(scalar: &[u64; 4], window_idx: usize, window_size: usize) -> usize {
    let bit_offset = window_idx * window_size;
    let mask = (1u64 << window_size) - 1;

    let limb_idx = bit_offset / 64;
    let bit_idx = bit_offset % 64;

    if limb_idx >= 4 {
        return 0;
    }

    let mut bits = (scalar[limb_idx] >> bit_idx) & mask;

    // Handle cross-limb boundary
    if bit_idx + window_size > 64 && limb_idx + 1 < 4 {
        let overflow = bit_idx + window_size - 64;
        bits |= (scalar[limb_idx + 1] & ((1u64 << overflow) - 1)) << (64 - bit_idx);
    }

    bits as usize
}

/// Choose optimal window size based on input size.
fn optimal_window(n: usize) -> usize {
    if n < 32 {
        3
    } else if n < 256 {
        5
    } else if n < 1024 {
        7
    } else if n < 4096 {
        9
    } else {
        11
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_fields::FieldElement;

    #[test]
    fn test_msm_g1_empty() {
        let result = msm_g1(&[], &[]);
        assert!(result.is_identity());
    }

    #[test]
    fn test_msm_g1_single() {
        let g = G1Affine::generator();
        let s = Fr::from(42u64);
        let result = msm_g1(&[g], &[s]);
        let expected = g * s;
        assert_eq!(result.to_affine(), expected.to_affine());
    }

    #[test]
    fn test_msm_g1_matches_naive() {
        let g = G1Affine::generator();
        let mut rng = rand::thread_rng();

        let n = 16;
        let scalars: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
        let bases: Vec<G1Affine> = (0..n)
            .map(|i| (g * Fr::from((i + 1) as u64)).to_affine())
            .collect();

        let result = msm_g1(&bases, &scalars);
        let expected = naive_msm_g1(&bases, &scalars);
        assert_eq!(result.to_affine(), expected.to_affine());
    }

    #[test]
    fn test_msm_g1_linearity() {
        // sum(2*s_i * g_i) = 2 * sum(s_i * g_i)
        let g = G1Affine::generator();

        let scalars: Vec<Fr> = (1..=5).map(|i| Fr::from(i as u64)).collect();
        let double_scalars: Vec<Fr> = scalars.iter().map(|s| s.double()).collect();
        let bases: Vec<G1Affine> = (1..=5)
            .map(|i| (g * Fr::from(i as u64)).to_affine())
            .collect();

        let result1 = msm_g1(&bases, &double_scalars);
        let result2 = msm_g1(&bases, &scalars).double();
        assert_eq!(result1.to_affine(), result2.to_affine());
    }

    #[test]
    fn test_msm_g2_single() {
        let g = G2Affine::generator();
        let s = Fr::from(42u64);
        let result = msm_g2(&[g], &[s]);
        let expected = g * s;
        assert_eq!(result.to_affine(), expected.to_affine());
    }
}
