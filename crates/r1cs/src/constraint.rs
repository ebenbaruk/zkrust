use std::ops::{Add, Mul, Neg, Sub};
use zkrust_fields::{FieldElement, Fr};

/// A variable in the constraint system.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Variable {
    /// The constant "1" variable (index 0 in the witness).
    One,
    /// A public input variable.
    Public(usize),
    /// A private (auxiliary/witness) variable.
    Private(usize),
}

impl Variable {
    /// Map this variable to its index in the full witness vector z = (1, public..., private...).
    pub fn index(&self, num_public: usize) -> usize {
        match self {
            Variable::One => 0,
            Variable::Public(i) => 1 + i,
            Variable::Private(i) => 1 + num_public + i,
        }
    }
}

/// A linear combination of variables: sum(coeff_i * var_i).
#[derive(Clone, Debug)]
pub struct LinearCombination {
    pub terms: Vec<(Fr, Variable)>,
}

impl LinearCombination {
    pub fn zero() -> Self {
        Self { terms: vec![] }
    }

    pub fn from_variable(var: Variable) -> Self {
        Self {
            terms: vec![(Fr::ONE, var)],
        }
    }

    pub fn from_constant(c: Fr) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Self {
                terms: vec![(c, Variable::One)],
            }
        }
    }

    /// Evaluate this linear combination given the full witness vector and num_public.
    pub fn evaluate_indexed(&self, witness: &[Fr], num_public: usize) -> Fr {
        let mut result = Fr::ZERO;
        for &(coeff, var) in &self.terms {
            let idx = var.index(num_public);
            result += coeff * witness[idx];
        }
        result
    }
}

impl Add for LinearCombination {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        self.terms.extend(rhs.terms);
        self
    }
}

impl Sub for LinearCombination {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        for (coeff, var) in rhs.terms {
            self.terms.push((-coeff, var));
        }
        self
    }
}

impl Neg for LinearCombination {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            terms: self.terms.into_iter().map(|(c, v)| (-c, v)).collect(),
        }
    }
}

impl Mul<Fr> for LinearCombination {
    type Output = Self;
    fn mul(self, scalar: Fr) -> Self {
        Self {
            terms: self
                .terms
                .into_iter()
                .map(|(c, v)| (c * scalar, v))
                .collect(),
        }
    }
}

impl From<Variable> for LinearCombination {
    fn from(var: Variable) -> Self {
        Self::from_variable(var)
    }
}

/// An R1CS constraint: A * B = C where A, B, C are linear combinations.
#[derive(Clone, Debug)]
pub struct Constraint {
    pub a: LinearCombination,
    pub b: LinearCombination,
    pub c: LinearCombination,
}

impl Constraint {
    pub fn new(a: LinearCombination, b: LinearCombination, c: LinearCombination) -> Self {
        Self { a, b, c }
    }

    /// Check if this constraint is satisfied by the given witness.
    pub fn is_satisfied(&self, witness: &[Fr], num_public: usize) -> bool {
        let a_val = self.a.evaluate_indexed(witness, num_public);
        let b_val = self.b.evaluate_indexed(witness, num_public);
        let c_val = self.c.evaluate_indexed(witness, num_public);
        a_val * b_val == c_val
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use zkrust_fields::FieldElement;

    #[test]
    fn test_variable_index() {
        let num_public = 3;
        assert_eq!(Variable::One.index(num_public), 0);
        assert_eq!(Variable::Public(0).index(num_public), 1);
        assert_eq!(Variable::Public(2).index(num_public), 3);
        assert_eq!(Variable::Private(0).index(num_public), 4);
        assert_eq!(Variable::Private(1).index(num_public), 5);
    }

    #[test]
    fn test_lc_evaluate() {
        // LC: 3*x + 2*y where x=Public(0)=5, y=Private(0)=7
        let lc = LinearCombination {
            terms: vec![
                (Fr::from(3u64), Variable::Public(0)),
                (Fr::from(2u64), Variable::Private(0)),
            ],
        };
        // witness = [1, 5, 7] (one, public(0), private(0))
        let witness = vec![Fr::ONE, Fr::from(5u64), Fr::from(7u64)];
        let result = lc.evaluate_indexed(&witness, 1);
        assert_eq!(result, Fr::from(29u64)); // 3*5 + 2*7 = 29
    }

    #[test]
    fn test_constraint_satisfied() {
        // x * y = z where x=3, y=4, z=12
        let x = Variable::Public(0);
        let y = Variable::Private(0);
        let z = Variable::Private(1);

        let constraint = Constraint::new(
            LinearCombination::from_variable(x),
            LinearCombination::from_variable(y),
            LinearCombination::from_variable(z),
        );

        // witness = [1, 3, 4, 12]
        let witness = vec![Fr::ONE, Fr::from(3u64), Fr::from(4u64), Fr::from(12u64)];
        assert!(constraint.is_satisfied(&witness, 1));

        // wrong witness
        let bad_witness = vec![Fr::ONE, Fr::from(3u64), Fr::from(4u64), Fr::from(11u64)];
        assert!(!constraint.is_satisfied(&bad_witness, 1));
    }

    #[test]
    fn test_lc_add_sub() {
        let a = LinearCombination::from_variable(Variable::Public(0));
        let b = LinearCombination::from_variable(Variable::Private(0));
        let sum = a + b;
        assert_eq!(sum.terms.len(), 2);

        let c = LinearCombination::from_variable(Variable::Public(0));
        let d = LinearCombination::from_variable(Variable::Private(0));
        let diff = c - d;
        assert_eq!(diff.terms.len(), 2);
    }
}
