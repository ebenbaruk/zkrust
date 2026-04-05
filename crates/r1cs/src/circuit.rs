use crate::system::{ConstraintSystem, SynthesisError};

/// Trait for circuits that can be synthesized into an R1CS constraint system.
///
/// A circuit defines the constraint structure and populates the witness.
pub trait Circuit {
    /// Synthesize the circuit into the given constraint system.
    ///
    /// This method should:
    /// 1. Allocate public and private variables
    /// 2. Add constraints relating those variables
    fn synthesize(&self, cs: &mut ConstraintSystem) -> Result<(), SynthesisError>;
}
