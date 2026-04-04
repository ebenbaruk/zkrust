# zkrust

A from-scratch implementation of the **Groth16 zk-SNARK** proving system in Rust, built on hand-rolled **BN254** elliptic curve and finite field arithmetic.

> *Prove you know a secret — without revealing it.*

## Architecture

```
crates/
├── fields/          # Finite field arithmetic (Fp, Fr, extension fields)
├── curves/          # Elliptic curve groups (BN254 G1, G2, pairing)
├── polynomials/     # Polynomial arithmetic, NTT, KZG commitments
├── r1cs/            # Rank-1 Constraint System compiler
├── groth16/         # Prover, Verifier, Trusted Setup
└── cli/             # Developer-facing CLI tool
```

## Features

- **Montgomery-form** finite field arithmetic with optimized CIOS multiplication
- **BN254** curve (Ethereum-compatible) with Ate pairing
- **NTT-based** polynomial multiplication
- **KZG** polynomial commitment scheme
- **Groth16** proving system with full trusted setup
- **Circuit DSL** for defining R1CS constraints in idiomatic Rust

## Example Circuits

- **Sudoku** — Prove you solved a Sudoku puzzle without revealing the solution
- **Range Proof** — Prove a value is within a range without revealing it
- **Merkle Proof** — Prove membership in a Merkle tree without revealing the leaf

## Quick Start

```bash
# Run trusted setup for a circuit
cargo run --bin zkrust -- setup --circuit range_proof --output ./keys

# Generate a proof
cargo run --bin zkrust -- prove --circuit range_proof --witness witness.json --keys ./keys

# Verify a proof
cargo run --bin zkrust -- verify --proof proof.bin --keys ./keys --public-inputs inputs.json
```

## Building

```bash
cargo build --workspace
cargo test --workspace
cargo bench
```

## Technical Details

| Property | Value |
|---|---|
| Curve | BN254 (alt-bn128) |
| Proving scheme | Groth16 |
| Field size | 254-bit prime |
| Proof size | 128 bytes (compressed) |
| MSRV | Rust 1.75+ |

## References

- [Groth16](https://eprint.iacr.org/2016/260.pdf) — Jens Groth, 2016
- [BN254 (EIP-196/197)](https://eips.ethereum.org/EIPS/eip-196) — Ethereum curve specification
- [The MoonMath Manual](https://leastauthority.com/community-matters/moonmath-manual/) — ZK math for engineers

## License

MIT

---

*This is an educational/portfolio project. The trusted setup is simulated locally and should not be used to secure real assets.*
