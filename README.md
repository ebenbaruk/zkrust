# zkrust

A from-scratch implementation of the **Groth16 zk-SNARK** proving system in Rust, built on hand-rolled **BN254** elliptic curve and finite field arithmetic. No external cryptographic dependencies — every algorithm from Montgomery multiplication to the Ate pairing is implemented in this repository.

## Overview

zkrust implements the full Groth16 zero-knowledge proof pipeline: finite field arithmetic, elliptic curve operations, polynomial commitments, constraint system compilation, and proof generation/verification. The library is structured as a Cargo workspace with six purpose-built crates totaling ~6,500 lines of Rust and 159 tests.

## Architecture

```
crates/
  fields/        Fp, Fr, Fp2, Fp6, Fp12  —  Montgomery-form arithmetic, extension field tower
  curves/        G1, G2, Ate pairing, Pippenger MSM
  polynomials/   Dense polynomials, radix-2 NTT, KZG commitments
  r1cs/          Variables, linear combinations, constraint system builder, gadgets
  groth16/       Trusted setup, prover, verifier
  cli/           Command-line interface for setup/prove/verify
```

## Quick Start

```bash
# Build and test
cargo build --workspace
cargo test --workspace

# End-to-end: prove that 5^2 = 25
cargo run --bin zkrust -- setup -c square -o ./keys
cargo run --bin zkrust -- prove -c square -w 5,25 -k ./keys -o proof.bin
cargo run --bin zkrust -- verify -p 5,25 -k ./keys -f proof.bin
# => Proof is VALID

# Run the range proof example
cargo run --example range_proof -p zkrust-groth16

# Benchmarks
cargo bench -p zkrust-fields --bench field_ops
cargo bench -p zkrust-groth16 --bench groth16_bench
```

## Built-in Circuits

| Circuit   | Description                     | Public Inputs | Constraints |
|-----------|---------------------------------|---------------|-------------|
| `square`  | Prove x * x = y                | x, y          | 1           |
| `cubic`   | Prove x^3 + x + 5 = y         | x, y          | 3           |
| `range32` | Prove x < 2^32 (bit decomposition) | x         | 33          |

## Implementation Details

### Finite Fields

- **Fp / Fr**: 254-bit prime fields in Montgomery form with CIOS multiplication
- **Extension tower**: Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[w]/(w^2-v)
- Fermat's little theorem inversion, Tonelli-Shanks square roots, Frobenius endomorphisms

### Elliptic Curves

- **G1**: y^2 = x^3 + 3 over Fp, Jacobian projective coordinates
- **G2**: D-type sextic twist over Fp2
- **Pairing**: Optimal Ate pairing with Miller loop over |6x+2| and final exponentiation
- **MSM**: Pippenger bucket method with adaptive window sizing

### Polynomials

- Dense polynomial arithmetic with Horner evaluation and long division
- Radix-2 Cooley-Tukey NTT using Fr's 2^28-th root of unity
- Coset FFT/iFFT for Groth16 quotient polynomial computation
- KZG polynomial commitment scheme (commit, open, verify via pairings)

### Groth16

- Trusted setup via Lagrange basis evaluation at secret tau
- Prover computes h(x) = (Az * Bz - Cz) / Z_H via NTT, builds proof with MSMs and random blinding
- Verifier checks pairing equation: e(A, B) = e(alpha, beta) * e(pub, gamma) * e(C, delta)
- Proof size: 3 group elements (1 G1 + 1 G2 + 1 G1)

## Benchmark Results

Measured on Apple Silicon (single-threaded, `--release`):

| Operation          | Time     |
|--------------------|----------|
| Fp multiply        | ~55 ns   |
| Fp12 multiply      | ~1.7 us  |
| Groth16 setup      | ~2.0 ms  |
| Groth16 prove      | ~1.0 ms  |
| Groth16 verify     | ~13.5 ms |

## Technical Specifications

| Property        | Value                          |
|-----------------|--------------------------------|
| Curve           | BN254 (alt-bn128)              |
| Base field      | p = 21888242871839275222246405745257275088696311157297823662689037894645226208583 |
| Scalar field    | r = 21888242871839275222246405745257275088548364400416034343698204186575808495617 |
| BN parameter    | x = 4965661367192848881         |
| Two-adicity (Fr)| 28                             |
| Proof system    | Groth16                        |
| MSRV            | Rust 1.75+                     |

## References

- Groth, J. (2016). [On the Size of Pairing-Based Non-interactive Arguments](https://eprint.iacr.org/2016/260.pdf). EUROCRYPT 2016.
- Barreto, P. & Naehrig, M. (2005). [Pairing-Friendly Elliptic Curves of Prime Order](https://eprint.iacr.org/2005/133.pdf). SAC 2005.
- [EIP-196](https://eips.ethereum.org/EIPS/eip-196) / [EIP-197](https://eips.ethereum.org/EIPS/eip-197) — Ethereum precompiles for BN254.
- Least Authority. [The MoonMath Manual](https://leastauthority.com/community-matters/moonmath-manual/).

## License

MIT

---

*This is an educational and portfolio project. The trusted setup is simulated locally and should not be used to secure real assets.*
