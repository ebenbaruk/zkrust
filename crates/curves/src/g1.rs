use zkrust_fields::Fp;

/// A point on the BN254 G1 curve in affine coordinates.
///
/// Curve equation: y² = x³ + 3 over Fp.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct G1Affine {
    pub x: Fp,
    pub y: Fp,
    pub infinity: bool,
}

impl G1Affine {
    pub const fn identity() -> Self {
        Self {
            x: Fp::ZERO,
            y: Fp::ZERO,
            infinity: true,
        }
    }
}
