use zkrust_fields::{FieldElement, Fr};
use zkrust_r1cs::{ConstraintSystem, LinearCombination};

/// Build a named circuit with the given witness values.
/// For setup, values can be empty (uses dummy values).
pub fn build_circuit(name: &str, values: &[u64]) -> Option<ConstraintSystem> {
    match name {
        "square" => Some(build_square(values)),
        "range32" => Some(build_range32(values)),
        "cubic" => Some(build_cubic(values)),
        _ => None,
    }
}

/// Square circuit: prove x * x = y.
/// Witness: [x, y] or [] for setup.
fn build_square(values: &[u64]) -> ConstraintSystem {
    let (x_val, y_val) = if values.len() >= 2 {
        (values[0], values[1])
    } else {
        (1, 1) // dummy for setup
    };

    let mut cs = ConstraintSystem::new();
    let x = cs.alloc_public_input(Fr::from(x_val));
    let y = cs.alloc_public_input(Fr::from(y_val));
    cs.enforce(
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(y),
    );
    cs
}

/// Range proof: prove x < 2^32 by bit decomposition.
/// Witness: [x] or [] for setup.
fn build_range32(values: &[u64]) -> ConstraintSystem {
    let x_val = if values.is_empty() { 1 } else { values[0] };

    let mut cs = ConstraintSystem::new();
    let x = cs.alloc_public_input(Fr::from(x_val));
    zkrust_r1cs::enforce_range(&mut cs, x, Fr::from(x_val), 32);
    cs
}

/// Cubic circuit: prove x^3 + x + 5 = y.
/// Witness: [x, y] or [] for setup.
fn build_cubic(values: &[u64]) -> ConstraintSystem {
    let (x_val, y_val) = if values.len() >= 2 {
        (values[0], values[1])
    } else {
        (1, 7) // 1 + 1 + 5 = 7
    };

    let mut cs = ConstraintSystem::new();
    let x = cs.alloc_public_input(Fr::from(x_val));
    let y = cs.alloc_public_input(Fr::from(y_val));

    // x_sq = x * x
    let x_sq_val = x_val.wrapping_mul(x_val);
    let x_sq = cs.alloc_private(Fr::from(x_sq_val));
    cs.enforce(
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(x_sq),
    );

    // x_cu = x_sq * x
    let x_cu_val = x_sq_val.wrapping_mul(x_val);
    let x_cu = cs.alloc_private(Fr::from(x_cu_val));
    cs.enforce(
        LinearCombination::from_variable(x_sq),
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(x_cu),
    );

    // (x_cu + x + 5) * 1 = y
    let lhs = LinearCombination::from_variable(x_cu)
        + LinearCombination::from_variable(x)
        + LinearCombination::from_constant(Fr::from(5u64));
    cs.enforce(
        lhs,
        LinearCombination::from_constant(Fr::ONE),
        LinearCombination::from_variable(y),
    );

    cs
}
