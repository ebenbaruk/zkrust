use crate::constraint::{Constraint, LinearCombination, Variable};
use zkrust_fields::{FieldElement, Fr};

/// A sparse matrix: each row is a list of (column_index, coefficient) pairs.
pub type SparseMatrix = Vec<Vec<(usize, Fr)>>;

/// Error types for constraint system operations.
#[derive(Debug, thiserror::Error)]
pub enum SynthesisError {
    #[error("constraint not satisfied at index {0}")]
    UnsatisfiedConstraint(usize),
    #[error("witness vector has wrong length: expected {expected}, got {got}")]
    WrongWitnessLength { expected: usize, got: usize },
    #[error("circuit synthesis error: {0}")]
    Other(String),
}

/// Builder for an R1CS constraint system.
///
/// Witness vector layout: z = (1, public_inputs..., private_inputs...).
pub struct ConstraintSystem {
    /// Number of public input variables (not counting the "1" variable).
    num_public: usize,
    /// Number of private (auxiliary) variables.
    num_private: usize,
    /// The constraints: A_i * B_i = C_i.
    constraints: Vec<Constraint>,
    /// Values assigned to public inputs.
    public_values: Vec<Fr>,
    /// Values assigned to private inputs.
    private_values: Vec<Fr>,
}

impl ConstraintSystem {
    pub fn new() -> Self {
        Self {
            num_public: 0,
            num_private: 0,
            constraints: Vec::new(),
            public_values: Vec::new(),
            private_values: Vec::new(),
        }
    }

    /// Allocate a public input variable with the given value.
    pub fn alloc_public_input(&mut self, value: Fr) -> Variable {
        let idx = self.num_public;
        self.num_public += 1;
        self.public_values.push(value);
        Variable::Public(idx)
    }

    /// Allocate a private (auxiliary/witness) variable with the given value.
    pub fn alloc_private(&mut self, value: Fr) -> Variable {
        let idx = self.num_private;
        self.num_private += 1;
        self.private_values.push(value);
        Variable::Private(idx)
    }

    /// Add a constraint: a * b = c.
    pub fn enforce(&mut self, a: LinearCombination, b: LinearCombination, c: LinearCombination) {
        self.constraints.push(Constraint::new(a, b, c));
    }

    /// Enforce that a variable is boolean: v * (1 - v) = 0.
    pub fn enforce_boolean(&mut self, v: Variable) {
        let one_minus_v =
            LinearCombination::from_constant(Fr::ONE) - LinearCombination::from_variable(v);
        self.enforce(
            LinearCombination::from_variable(v),
            one_minus_v,
            LinearCombination::zero(),
        );
    }

    /// Enforce that two variables are equal: (a - b) * 1 = 0.
    pub fn enforce_equal(&mut self, a: Variable, b: Variable) {
        let diff = LinearCombination::from_variable(a) - LinearCombination::from_variable(b);
        self.enforce(
            diff,
            LinearCombination::from_constant(Fr::ONE),
            LinearCombination::zero(),
        );
    }

    /// Number of public inputs (not counting the "1" variable).
    pub fn num_public_inputs(&self) -> usize {
        self.num_public
    }

    /// Number of private (auxiliary) variables.
    pub fn num_private_inputs(&self) -> usize {
        self.num_private
    }

    /// Total number of variables including the "1" variable.
    pub fn num_variables(&self) -> usize {
        1 + self.num_public + self.num_private
    }

    /// Number of constraints.
    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }

    /// Build the full witness vector: z = (1, public..., private...).
    pub fn witness(&self) -> Vec<Fr> {
        let mut w = Vec::with_capacity(self.num_variables());
        w.push(Fr::ONE);
        w.extend_from_slice(&self.public_values);
        w.extend_from_slice(&self.private_values);
        w
    }

    /// Check all constraints are satisfied by the current witness.
    pub fn is_satisfied(&self) -> bool {
        let witness = self.witness();
        self.constraints
            .iter()
            .all(|c| c.is_satisfied(&witness, self.num_public))
    }

    /// Check all constraints, returning the index of the first unsatisfied one.
    pub fn verify(&self) -> Result<(), SynthesisError> {
        let witness = self.witness();
        for (i, constraint) in self.constraints.iter().enumerate() {
            if !constraint.is_satisfied(&witness, self.num_public) {
                return Err(SynthesisError::UnsatisfiedConstraint(i));
            }
        }
        Ok(())
    }

    /// Extract the sparse A, B, C matrices for Groth16.
    ///
    /// Each matrix is `Vec<Vec<(usize, Fr)>>` where:
    /// - Outer vec: one entry per constraint (row)
    /// - Inner vec: (column_index, coefficient) pairs
    /// - Column indices correspond to the witness vector z = (1, public..., private...)
    pub fn to_matrices(&self) -> (SparseMatrix, SparseMatrix, SparseMatrix) {
        let m = self.constraints.len();
        let mut a_matrix = Vec::with_capacity(m);
        let mut b_matrix = Vec::with_capacity(m);
        let mut c_matrix = Vec::with_capacity(m);

        for constraint in &self.constraints {
            a_matrix.push(lc_to_sparse(&constraint.a, self.num_public));
            b_matrix.push(lc_to_sparse(&constraint.b, self.num_public));
            c_matrix.push(lc_to_sparse(&constraint.c, self.num_public));
        }

        (a_matrix, b_matrix, c_matrix)
    }

    /// Access the constraints.
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// Access the public input values.
    pub fn public_inputs(&self) -> &[Fr] {
        &self.public_values
    }

    /// Access the private input values.
    pub fn private_inputs(&self) -> &[Fr] {
        &self.private_values
    }
}

impl Default for ConstraintSystem {
    fn default() -> Self {
        Self::new()
    }
}

/// Convert a linear combination to sparse matrix row format.
fn lc_to_sparse(lc: &LinearCombination, num_public: usize) -> Vec<(usize, Fr)> {
    lc.terms
        .iter()
        .map(|&(coeff, var)| (var.index(num_public), coeff))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_multiply() {
        // Prove: x * y = z where x=3, y=4, z=12
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_private(Fr::from(4u64));
        let z = cs.alloc_public_input(Fr::from(12u64));

        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
            LinearCombination::from_variable(z),
        );

        assert!(cs.is_satisfied());
        assert_eq!(cs.num_constraints(), 1);
        assert_eq!(cs.num_public_inputs(), 2);
        assert_eq!(cs.num_private_inputs(), 1);
        assert_eq!(cs.num_variables(), 4);
    }

    #[test]
    fn test_unsatisfied_constraint() {
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_private(Fr::from(4u64));
        let z = cs.alloc_public_input(Fr::from(11u64)); // wrong!

        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
            LinearCombination::from_variable(z),
        );

        assert!(!cs.is_satisfied());
        assert!(cs.verify().is_err());
    }

    #[test]
    fn test_boolean_constraint() {
        let mut cs = ConstraintSystem::new();
        let b = cs.alloc_private(Fr::ONE); // boolean = 1
        cs.enforce_boolean(b);
        assert!(cs.is_satisfied());

        let mut cs2 = ConstraintSystem::new();
        let b2 = cs2.alloc_private(Fr::ZERO); // boolean = 0
        cs2.enforce_boolean(b2);
        assert!(cs2.is_satisfied());

        let mut cs3 = ConstraintSystem::new();
        let b3 = cs3.alloc_private(Fr::from(2u64)); // not boolean!
        cs3.enforce_boolean(b3);
        assert!(!cs3.is_satisfied());
    }

    #[test]
    fn test_equality_constraint() {
        let mut cs = ConstraintSystem::new();
        let a = cs.alloc_public_input(Fr::from(42u64));
        let b = cs.alloc_private(Fr::from(42u64));
        cs.enforce_equal(a, b);
        assert!(cs.is_satisfied());

        let mut cs2 = ConstraintSystem::new();
        let a2 = cs2.alloc_public_input(Fr::from(42u64));
        let b2 = cs2.alloc_private(Fr::from(43u64));
        cs2.enforce_equal(a2, b2);
        assert!(!cs2.is_satisfied());
    }

    #[test]
    fn test_matrix_extraction() {
        // x * y = z
        let mut cs = ConstraintSystem::new();
        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_private(Fr::from(4u64));
        let z = cs.alloc_public_input(Fr::from(12u64));

        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
            LinearCombination::from_variable(z),
        );

        let (a, b, c) = cs.to_matrices();
        assert_eq!(a.len(), 1);
        assert_eq!(b.len(), 1);
        assert_eq!(c.len(), 1);

        // A row: x = Public(0) → index 1
        assert_eq!(a[0], vec![(1, Fr::ONE)]);
        // B row: y = Private(0) → index 3
        assert_eq!(b[0], vec![(3, Fr::ONE)]);
        // C row: z = Public(1) → index 2
        assert_eq!(c[0], vec![(2, Fr::ONE)]);
    }

    #[test]
    fn test_witness_vector() {
        let mut cs = ConstraintSystem::new();
        cs.alloc_public_input(Fr::from(5u64));
        cs.alloc_public_input(Fr::from(7u64));
        cs.alloc_private(Fr::from(35u64));

        let w = cs.witness();
        assert_eq!(w.len(), 4);
        assert_eq!(w[0], Fr::ONE);
        assert_eq!(w[1], Fr::from(5u64));
        assert_eq!(w[2], Fr::from(7u64));
        assert_eq!(w[3], Fr::from(35u64));
    }

    #[test]
    fn test_quadratic_circuit() {
        // Prove: x^2 + x + 5 = y
        // Constraints:
        //   sym_1 = x * x        (intermediate squaring)
        //   y = sym_1 + x + 5    (rearranged as: (sym_1 + x + 5) * 1 = y)
        let mut cs = ConstraintSystem::new();

        let x = cs.alloc_public_input(Fr::from(3u64));
        let y = cs.alloc_public_input(Fr::from(17u64)); // 9 + 3 + 5 = 17
        let sym_1 = cs.alloc_private(Fr::from(9u64)); // x^2

        // sym_1 = x * x
        cs.enforce(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(sym_1),
        );

        // (sym_1 + x + 5) * 1 = y
        let lhs = LinearCombination::from_variable(sym_1)
            + LinearCombination::from_variable(x)
            + LinearCombination::from_constant(Fr::from(5u64));
        cs.enforce(
            lhs,
            LinearCombination::from_constant(Fr::ONE),
            LinearCombination::from_variable(y),
        );

        assert!(cs.is_satisfied());
        assert_eq!(cs.num_constraints(), 2);
    }
}
