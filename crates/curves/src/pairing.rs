use crate::g1::G1Affine;
use crate::g2::G2Affine;
use zkrust_fields::{FieldElement, Fp, Fp12, Fp2, Fp6};

/// BN parameter x = 0x44e992b44a6909f1.
/// This is the generator of the BN254 curve family: p = 36x⁴ + 36x³ + 24x² + 6x + 1.
const BN_X: u64 = 4965661367192848881;

/// (p-1)/3 in little-endian u64 limbs, for Frobenius on G2.
const P_MINUS_1_OVER_3: [u64; 4] = [
    0x69602EB24829A9C2,
    0xDD2B2385CD7B4384,
    0xE81AC1E7808072C9,
    0x10216F7BA065E00D,
];

/// (p-1)/2 in little-endian u64 limbs, for Frobenius on G2.
const P_MINUS_1_OVER_2: [u64; 4] = [
    0x9E10460B6C3E7EA3,
    0xCBC0B548B438E546,
    0xDC2822DB40C0AC2E,
    0x183227397098D014,
];

/// Compute the optimal Ate pairing e(P, Q) for P ∈ G1 and Q ∈ G2.
pub fn ate_pairing(p: &G1Affine, q: &G2Affine) -> Fp12 {
    if p.is_identity() || q.is_identity() {
        return Fp12::ONE;
    }
    let f = miller_loop(p, q);
    final_exponentiation(&f)
}

/// Multi-pairing: shares final exponentiation.
pub fn multi_pairing(pairs: &[(G1Affine, G2Affine)]) -> Fp12 {
    let mut f = Fp12::ONE;
    for (p, q) in pairs {
        if !p.is_identity() && !q.is_identity() {
            f *= miller_loop(p, q);
        }
    }
    final_exponentiation(&f)
}

/// Line evaluation result: sparse Fp12 at positions 0, 3, 4.
///
/// For a D-type sextic twist with tower Fp12 = Fp6[w]/(w²-v), Fp6 = Fp2[v]/(v³-ξ):
///   - Position 0 (= 1):  r0 → Fp12.c0.c0
///   - Position 3 (= w):  r3 → Fp12.c1.c0
///   - Position 4 (= vw): r4 → Fp12.c1.c1
///
/// The untwisting map ψ⁻¹(x', y') = (x'·ω², y'·ω³) where ω² = v, ω³ = vw
/// produces line evaluations at P with this sparse structure.
struct LineEval {
    r0: Fp2,
    r3: Fp2,
    r4: Fp2,
}

impl LineEval {
    fn to_fp12(&self) -> Fp12 {
        Fp12::new(
            Fp6::new(self.r0, Fp2::ZERO, Fp2::ZERO),
            Fp6::new(self.r3, self.r4, Fp2::ZERO),
        )
    }
}

/// Miller loop for the optimal Ate pairing.
/// Uses affine coordinates for the running point T on E'.
fn miller_loop(p: &G1Affine, q: &G2Affine) -> Fp12 {
    let six_x_plus_2: u128 = 6 * (BN_X as u128) + 2;

    let mut f = Fp12::ONE;
    let mut tx = q.x;
    let mut ty = q.y;

    let bits = {
        let mut v = Vec::new();
        let mut val = six_x_plus_2;
        while val > 0 {
            v.push((val & 1) as u8);
            val >>= 1;
        }
        v.reverse();
        v
    };

    for &bit in bits[1..].iter() {
        f = f.square();

        let line = line_double(&mut tx, &mut ty, p);
        f *= line.to_fp12();

        if bit == 1 {
            let line = line_add(&mut tx, &mut ty, &q.x, &q.y, p);
            f *= line.to_fp12();
        }
    }

    // Frobenius corrections: add π(Q) and -π²(Q)
    let (q1x, q1y) = frobenius_g2(q, 1);
    let (q2x, q2y_neg) = {
        let (x, y) = frobenius_g2(q, 2);
        (x, -y) // negate y for -π²(Q)
    };

    let line = line_add(&mut tx, &mut ty, &q1x, &q1y, p);
    f *= line.to_fp12();

    let line = line_add(&mut tx, &mut ty, &q2x, &q2y_neg, p);
    f *= line.to_fp12();

    f
}

/// Doubling step: tangent line at T = (tx, ty) on E', evaluated at P ∈ G1.
/// Updates T := 2T.
///
/// The tangent slope on E' is λ' = 3·x²/(2·y).
/// After untwisting, the line at P = (xP, yP) is:
///   l(P) = yP + (-λ'·xP)·ω + (λ'·xT - yT)·ω³
fn line_double(tx: &mut Fp2, ty: &mut Fp2, p: &G1Affine) -> LineEval {
    let xx = tx.square();
    let three_xx = xx.double() + xx;
    let two_y = ty.double();
    let lambda = three_xx * two_y.inv().unwrap();

    // 2T in affine
    let new_x = lambda.square() - tx.double();
    let new_y = lambda * (*tx - new_x) - *ty;

    // Line evaluation
    let r0 = Fp2::new(p.y, Fp::ZERO);
    let r3 = Fp2::new(-(lambda.c0 * p.x), -(lambda.c1 * p.x));
    let r4 = lambda * *tx - *ty;

    *tx = new_x;
    *ty = new_y;

    LineEval { r0, r3, r4 }
}

/// Addition step: chord through T and Q on E', evaluated at P ∈ G1.
/// Updates T := T + Q.
///
/// The chord slope on E' is λ' = (yT - yQ)/(xT - xQ).
/// After untwisting, the line at P is:
///   l(P) = yP + (-λ'·xP)·ω + (λ'·xQ - yQ)·ω³
fn line_add(tx: &mut Fp2, ty: &mut Fp2, qx: &Fp2, qy: &Fp2, p: &G1Affine) -> LineEval {
    let lambda = (*ty - *qy) * (*tx - *qx).inv().unwrap();

    let new_x = lambda.square() - *tx - *qx;
    let new_y = lambda * (*qx - new_x) - *qy;

    let r0 = Fp2::new(p.y, Fp::ZERO);
    let r3 = Fp2::new(-(lambda.c0 * p.x), -(lambda.c1 * p.x));
    let r4 = lambda * *qx - *qy;

    *tx = new_x;
    *ty = new_y;

    LineEval { r0, r3, r4 }
}

/// Frobenius endomorphism on G2 twist coordinates.
/// π^k(x, y) = (x^(p^k) · ξ^((p^k-1)/3), y^(p^k) · ξ^((p^k-1)/2))
///
/// Applied iteratively: each step conjugates (= x^p for Fp2) and
/// multiplies by ξ^((p-1)/3) or ξ^((p-1)/2).
fn frobenius_g2(q: &G2Affine, power: usize) -> (Fp2, Fp2) {
    if power == 0 || q.is_identity() {
        return (q.x, q.y);
    }

    let gamma1 = xi_pow_p_minus_1_over_3();
    let gamma2 = xi_pow_p_minus_1_over_2();

    let mut x = q.x;
    let mut y = q.y;

    for _ in 0..power {
        x = x.frobenius_map() * gamma1;
        y = y.frobenius_map() * gamma2;
    }

    (x, y)
}

/// ξ^((p-1)/3) where ξ = 9 + u.
fn xi_pow_p_minus_1_over_3() -> Fp2 {
    let xi = Fp2::new(Fp::from(9u64), Fp::from(1u64));
    xi.pow(&P_MINUS_1_OVER_3)
}

/// ξ^((p-1)/2) where ξ = 9 + u.
fn xi_pow_p_minus_1_over_2() -> Fp2 {
    let xi = Fp2::new(Fp::from(9u64), Fp::from(1u64));
    xi.pow(&P_MINUS_1_OVER_2)
}

/// Final exponentiation: f^((p^12 - 1) / r).
/// Split into easy part and hard part.
fn final_exponentiation(f: &Fp12) -> Fp12 {
    // Easy part: f^((p^6 - 1)(p^2 + 1))
    let f_conj = f.conjugate();
    let f_inv = f.inv().unwrap();
    let r0 = f_conj * f_inv; // f^(p^6 - 1)
    let r1 = r0.frobenius_map().frobenius_map() * r0; // r0^(p^2 + 1)
    hard_part(&r1)
}

/// Hard part of final exponentiation.
/// Computes f^((p^4 - p^2 + 1) / r).
fn hard_part(f: &Fp12) -> Fp12 {
    // (p^4 - p^2 + 1) / r expressed as little-endian u64 limbs.
    #[rustfmt::skip]
    const HARD_EXP: [u64; 12] = [
        0xE81BB482CCDF42B1, 0x5ABF5CC4F49C36D4,
        0xF1154E7E1DA014FD, 0xDCC7B44C87CDBACF,
        0xAAA441E3954BCF8A, 0x6B887D56D5095F23,
        0x79581E16F3FD90C6, 0x3B1B1355D189227D,
        0x4E529A5861876F6B, 0x6C0EB522D5B12278,
        0x331EC15183177FAF, 0x01BAAA710B0759AD,
    ];
    f.pow(&HARD_EXP)
}

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_fields::Fr;

    #[test]
    fn test_pairing_identity_g1() {
        assert_eq!(
            ate_pairing(&G1Affine::identity(), &G2Affine::generator()),
            Fp12::ONE
        );
    }

    #[test]
    fn test_pairing_identity_g2() {
        assert_eq!(
            ate_pairing(&G1Affine::generator(), &G2Affine::identity()),
            Fp12::ONE
        );
    }

    #[test]
    fn test_pairing_non_degeneracy() {
        let result = ate_pairing(&G1Affine::generator(), &G2Affine::generator());
        assert_ne!(result, Fp12::ONE);
        assert_ne!(result, Fp12::ZERO);
    }

    #[test]
    fn test_pairing_inverse() {
        let p = G1Affine::generator();
        let q = G2Affine::generator();
        let e1 = ate_pairing(&p, &q);
        let e2 = ate_pairing(&(-p), &q);
        assert_eq!(e1 * e2, Fp12::ONE, "e(P,Q)*e(-P,Q) should be 1");
    }

    #[test]
    fn test_pairing_bilinearity_lhs() {
        let p = G1Affine::generator();
        let q = G2Affine::generator();
        let a = Fr::from(7u64);
        let ap = (p * a).to_affine();
        let lhs = ate_pairing(&ap, &q);
        let rhs = ate_pairing(&p, &q).pow(&[7, 0, 0, 0]);
        assert_eq!(lhs, rhs, "e(aP, Q) should equal e(P, Q)^a");
    }

    #[test]
    fn test_pairing_bilinearity_rhs() {
        let p = G1Affine::generator();
        let q = G2Affine::generator();
        let b = Fr::from(5u64);
        let bq = (q * b).to_affine();
        let lhs = ate_pairing(&p, &bq);
        let rhs = ate_pairing(&p, &q).pow(&[5, 0, 0, 0]);
        assert_eq!(lhs, rhs, "e(P, bQ) should equal e(P, Q)^b");
    }

    #[test]
    fn test_pairing_bilinearity_both() {
        let p = G1Affine::generator();
        let q = G2Affine::generator();
        let a = Fr::from(3u64);
        let b = Fr::from(5u64);
        let lhs = ate_pairing(&(p * a).to_affine(), &(q * b).to_affine());
        let rhs = ate_pairing(&p, &q).pow(&[15, 0, 0, 0]);
        assert_eq!(lhs, rhs, "e(aP, bQ) should equal e(P, Q)^(ab)");
    }
}
