use crate::constraint::{LinearCombination, Variable};
use crate::system::ConstraintSystem;
use zkrust_fields::{FieldElement, Fr};

/// Allocate a private boolean variable and enforce it is 0 or 1.
pub fn alloc_boolean(cs: &mut ConstraintSystem, value: bool) -> Variable {
    let v = cs.alloc_private(if value { Fr::ONE } else { Fr::ZERO });
    cs.enforce_boolean(v);
    v
}

/// Decompose a field element into `num_bits` boolean variables (little-endian).
/// Enforces that the bits reconstruct the original value.
pub fn alloc_bits(cs: &mut ConstraintSystem, value: Fr, num_bits: usize) -> Vec<Variable> {
    let raw = value.to_raw();
    let mut bits = Vec::with_capacity(num_bits);

    for i in 0..num_bits {
        let limb = i / 64;
        let bit_idx = i % 64;
        let bit_val = if limb < 4 {
            (raw[limb] >> bit_idx) & 1 == 1
        } else {
            false
        };
        bits.push(alloc_boolean(cs, bit_val));
    }

    // Enforce: sum(bits[i] * 2^i) = value
    let mut lc = LinearCombination::zero();
    let mut power = Fr::ONE;
    let two = Fr::from(2u64);
    for &bit in &bits {
        lc = lc + LinearCombination::from_variable(bit) * power;
        power *= two;
    }

    let value_var = cs.alloc_public_input(value);
    cs.enforce(
        lc,
        LinearCombination::from_constant(Fr::ONE),
        LinearCombination::from_variable(value_var),
    );

    bits
}

/// Enforce a < 2^num_bits by decomposing into bits.
/// Returns the bit decomposition.
pub fn enforce_range(
    cs: &mut ConstraintSystem,
    var: Variable,
    value: Fr,
    num_bits: usize,
) -> Vec<Variable> {
    let raw = value.to_raw();
    let mut bits = Vec::with_capacity(num_bits);

    for i in 0..num_bits {
        let limb = i / 64;
        let bit_idx = i % 64;
        let bit_val = if limb < 4 {
            (raw[limb] >> bit_idx) & 1 == 1
        } else {
            false
        };
        bits.push(alloc_boolean(cs, bit_val));
    }

    // Enforce: sum(bits[i] * 2^i) = var
    let mut lc = LinearCombination::zero();
    let mut power = Fr::ONE;
    let two = Fr::from(2u64);
    for &bit in &bits {
        lc = lc + LinearCombination::from_variable(bit) * power;
        power *= two;
    }

    cs.enforce(
        lc,
        LinearCombination::from_constant(Fr::ONE),
        LinearCombination::from_variable(var),
    );

    bits
}

/// Conditionally select: if condition then a else b.
/// Returns a variable equal to condition * a + (1 - condition) * b.
pub fn conditional_select(
    cs: &mut ConstraintSystem,
    condition: Variable,
    a: Variable,
    b: Variable,
    cond_val: bool,
    a_val: Fr,
    b_val: Fr,
) -> Variable {
    let result_val = if cond_val { a_val } else { b_val };
    let result = cs.alloc_private(result_val);

    // result = condition * (a - b) + b
    // Rearranged: condition * (a - b) = result - b
    cs.enforce(
        LinearCombination::from_variable(condition),
        LinearCombination::from_variable(a) - LinearCombination::from_variable(b),
        LinearCombination::from_variable(result) - LinearCombination::from_variable(b),
    );

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alloc_boolean() {
        let mut cs = ConstraintSystem::new();
        let _b0 = alloc_boolean(&mut cs, false);
        let _b1 = alloc_boolean(&mut cs, true);
        assert!(cs.is_satisfied());
        assert_eq!(cs.num_constraints(), 2);
    }

    #[test]
    fn test_alloc_bits() {
        let mut cs = ConstraintSystem::new();
        let value = Fr::from(42u64); // 101010 in binary
        let bits = alloc_bits(&mut cs, value, 8);
        assert_eq!(bits.len(), 8);
        assert!(cs.is_satisfied());
        // 8 boolean constraints + 1 reconstruction constraint
        assert_eq!(cs.num_constraints(), 9);
    }

    #[test]
    fn test_range_proof() {
        let mut cs = ConstraintSystem::new();
        let value = Fr::from(255u64);
        let var = cs.alloc_public_input(value);
        let _bits = enforce_range(&mut cs, var, value, 8);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn test_range_proof_overflow() {
        let mut cs = ConstraintSystem::new();
        let value = Fr::from(256u64); // doesn't fit in 8 bits
        let var = cs.alloc_public_input(value);
        let _bits = enforce_range(&mut cs, var, value, 8);
        // This should fail because 256 can't be decomposed into 8 bits
        assert!(!cs.is_satisfied());
    }

    #[test]
    fn test_conditional_select() {
        let mut cs = ConstraintSystem::new();
        let a = cs.alloc_private(Fr::from(10u64));
        let b = cs.alloc_private(Fr::from(20u64));

        // Select a (condition = true)
        let cond = alloc_boolean(&mut cs, true);
        let result =
            conditional_select(&mut cs, cond, a, b, true, Fr::from(10u64), Fr::from(20u64));

        // Check result = 10
        cs.enforce_equal(result, a);
        assert!(cs.is_satisfied());
    }

    #[test]
    fn test_conditional_select_false() {
        let mut cs = ConstraintSystem::new();
        let a = cs.alloc_private(Fr::from(10u64));
        let b = cs.alloc_private(Fr::from(20u64));

        let cond = alloc_boolean(&mut cs, false);
        let result =
            conditional_select(&mut cs, cond, a, b, false, Fr::from(10u64), Fr::from(20u64));

        cs.enforce_equal(result, b);
        assert!(cs.is_satisfied());
    }
}
