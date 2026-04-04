# 📄 Product Requirements Document
## Zero-Knowledge Proof System in Rust
**Version:** 1.0
**Author:** TBD
**Status:** Draft

---

## 1. Overview

### 1.1 Summary
A from-scratch implementation of a Zero-Knowledge Proof system in Rust, exposing a clean developer-facing library and CLI. The system will allow a **prover** to convince a **verifier** that they know a secret value satisfying some constraint — without revealing the secret itself.

The project targets Groth16 (zk-SNARK) as the primary proving scheme, built on top of hand-rolled elliptic curve and finite field arithmetic.

### 1.2 Goals
- Demonstrate deep cryptographic engineering competence in Rust
- Build a working, benchmarked, documented ZK library publishable on `crates.io`
- Serve as a portfolio centerpiece for roles at ZK-focused companies (StarkWare, Aztec, Polygon, zkSync, Risc0)
- Produce an accompanying technical blog post

### 1.3 Non-Goals
- Production security hardening (no side-channel resistance)
- A fully general zkEVM
- Competing with `arkworks` or `bellman` in performance

---

## 2. Background & Motivation

Zero-Knowledge Proofs are the most mathematically demanding and in-demand skill in the blockchain industry. Most engineers copy-paste `arkworks`. Building one from scratch signals:

- You understand the **math** (elliptic curves, pairings, polynomial IOP)
- You can write **safe and unsafe Rust** at a systems level
- You are not afraid of **hard problems**

Companies hiring for this: **StarkWare, Aztec Network, Polygon, zkSync (Matter Labs), Risc Zero, Succinct Labs, Nil Foundation, Penumbra.**

---

## 3. System Architecture

```
zk-rust/
├── crates/
│   ├── fields/          # Finite field arithmetic (Fp, Fr)
│   ├── curves/          # Elliptic curve groups (BN254)
│   ├── polynomials/     # Polynomial arithmetic & commitments
│   ├── r1cs/            # Rank-1 Constraint System compiler
│   ├── groth16/         # Prover, Verifier, Trusted Setup
│   └── cli/             # Developer-facing CLI tool
├── examples/
│   ├── sudoku/          # Prove you solved a Sudoku
│   ├── range_proof/     # Prove a number is in [0, 2^n)
│   └── merkle/          # Prove membership in a Merkle tree
├── benches/             # Criterion benchmarks
└── docs/                # Architecture writeups
```

---

## 4. Core Components

### 4.1 Finite Field Arithmetic (`crates/fields`)
The bedrock of everything. All ZK math lives in finite fields.

| Feature | Description |
|---|---|
| `Fp` — Base field | Prime field over BN254's base prime |
| `Fr` — Scalar field | Prime field over BN254's scalar prime |
| Montgomery multiplication | Fast modular mul using Montgomery form |
| Field traits | `Add`, `Sub`, `Mul`, `Inv`, `Pow` impls |
| Batch inversion | Montgomery batch inversion trick |

**Acceptance Criteria:**
- All field operations pass reference test vectors
- Montgomery multiplication benchmarks within 2x of `arkworks`

---

### 4.2 Elliptic Curve Groups (`crates/curves`)
BN254 (also called alt-bn128) — the curve used by Ethereum's precompiles.

| Feature | Description |
|---|---|
| G1 group | Short Weierstrass curve over `Fp` |
| G2 group | Twist curve over `Fp2` |
| Point addition / doubling | Both affine and projective coordinates |
| Scalar multiplication | Double-and-add, then windowed NAF |
| Ate pairing | Bilinear pairing `e: G1 × G2 → GT` |
| Subgroup checks | Reject points not in the correct group |

**Acceptance Criteria:**
- Pairing outputs match Ethereum's `bn256` precompile on test vectors
- Scalar mul handles edge cases (point at infinity, zero scalar)

---

### 4.3 Polynomial Arithmetic (`crates/polynomials`)

| Feature | Description |
|---|---|
| Dense polynomial | Coefficient-vector representation |
| Evaluation / interpolation | Lagrange interpolation |
| FFT/NTT | Fast Number Theoretic Transform for multiplication |
| KZG commitments | Kate-Zaverucha-Goldberg polynomial commitments |

**Acceptance Criteria:**
- KZG open/verify round-trips correctly
- FFT-based polynomial multiplication matches naive O(n²) for all degrees ≤ 2^16

---

### 4.4 R1CS Compiler (`crates/r1cs`)
Rank-1 Constraint Systems are the "circuit" format Groth16 consumes.

| Feature | Description |
|---|---|
| Constraint builder DSL | Ergonomic Rust API to define circuits |
| Witness generation | Populate all wire values given private inputs |
| R1CS validation | Check `Az ∘ Bz = Cz` over the scalar field |
| Circuit examples | Included: hash preimage, range check, Merkle proof |

**Acceptance Criteria:**
- A developer can define a new circuit in < 50 lines of Rust
- Constraint system correctly rejects invalid witnesses

---

### 4.5 Groth16 Prover & Verifier (`crates/groth16`)

| Feature | Description |
|---|---|
| Trusted Setup (Powers of Tau) | Phase 1 & Phase 2 ceremony simulation (local) |
| Proving Key / Verifying Key | Generated from R1CS + toxic waste |
| Prover | Generates a proof `π = (A, B, C)` in G1/G2 |
| Verifier | Checks pairing equation on-chain style |
| Proof serialization | Bincode + hex encoding for interop |

**Acceptance Criteria:**
- Honest proofs always verify
- Tampered proofs always fail
- Prove + verify round-trip on all 3 example circuits

---

### 4.6 CLI (`crates/cli`)

```bash
# Run trusted setup for a circuit
zk setup --circuit sudoku --output ./keys

# Generate a proof
zk prove --circuit sudoku --input witness.json --keys ./keys

# Verify a proof
zk verify --proof proof.bin --keys ./keys --public-inputs inputs.json
```

---

## 5. Example Circuits (Demos)

These are the "wow" demos that go in the README.

### 5.1 Sudoku Proof
> *"I know the solution to this Sudoku — but I won't show you."*

- Public input: the puzzle (given cells)
- Private input: the full solution
- Constraints: rows, columns, boxes each contain 1–9 exactly once

### 5.2 Range Proof
> *"I know a number x such that 0 ≤ x < 2^32 — without revealing x."*

- Classic ZK primitive used in confidential transactions
- Demonstrates bit decomposition circuit pattern

### 5.3 Merkle Membership Proof
> *"I know a leaf in this Merkle tree — without revealing which one."*

- Public input: Merkle root
- Private input: leaf + sibling path
- Demonstrates hash circuit composition (Poseidon hash in-circuit)

---

## 6. Phased Build Plan

### Phase 1 — Math Foundations (Weeks 1–3)
- [ ] Implement `Fp` and `Fr` with Montgomery form
- [ ] Implement G1 and G2 point arithmetic
- [ ] Implement and verify Ate pairing against Ethereum test vectors
- [ ] Write 100% unit test coverage for all field/curve ops

### Phase 2 — Polynomial Layer (Weeks 4–5)
- [ ] Dense polynomial arithmetic
- [ ] NTT-based fast multiplication
- [ ] KZG commitment scheme

### Phase 3 — R1CS & Circuit DSL (Weeks 6–7)
- [ ] R1CS struct and constraint builder
- [ ] Witness generation engine
- [ ] Implement range proof circuit as first test

### Phase 4 — Groth16 (Weeks 8–10)
- [ ] Trusted setup (Powers of Tau, local simulation)
- [ ] Prover implementation
- [ ] Verifier implementation
- [ ] Full round-trip tests

### Phase 5 — Polish & Ship (Weeks 11–12)
- [ ] Sudoku + Merkle example circuits
- [ ] CLI tool
- [ ] Criterion benchmarks
- [ ] README with diagrams
- [ ] Blog post: *"Building a ZK Proof System from Scratch in Rust"*
- [ ] Publish to `crates.io`

---

## 7. Technical Constraints

| Constraint | Decision |
|---|---|
| Curve | BN254 (Ethereum-compatible, well-tested vectors) |
| Proving scheme | Groth16 (most common, most understood) |
| Async | None — all CPU-bound, use Rayon for parallelism |
| Unsafe | Allowed only in `fields/` for Montgomery mul hot path |
| MSRV | Rust 1.75+ (stable) |
| Dependencies | Minimal — `rayon`, `criterion`, `serde`, `thiserror` |

---

## 8. Success Metrics

| Metric | Target |
|---|---|
| Groth16 proof generation (Sudoku circuit) | < 5 seconds on M2 MacBook |
| Verification time | < 10ms |
| GitHub stars (6 months post-launch) | 200+ |
| crates.io downloads | 500+ |
| Blog post reach | Featured on r/rust or Hacker News |
| Recruiter inbounds after launch | At least 3 from target companies |

---

## 9. Risks

| Risk | Mitigation |
|---|---|
| Math complexity causes project stall | Start with finite fields only — ship incrementally |
| BN254 pairing is hard to implement correctly | Validate against Ethereum precompile at every step |
| Scope creep into a full zkEVM | Strictly enforce non-goals — Groth16 only |
| Trusted setup "toxic waste" misunderstood | Add big disclaimers: "for educational use, not production" |

---

## 10. References

- [Groth16 original paper](https://eprint.iacr.org/2016/260.pdf) — Jens Groth, 2016
- [arkworks-rs](https://github.com/arkworks-rs) — Reference Rust ZK library
- [BN254 EIP-196/197](https://eips.ethereum.org/EIPS/eip-196) — Ethereum curve spec
- [The MoonMath Manual](https://leastauthority.com/community-matters/moonmath-manual/) — ZK math for engineers
- [ZKProof Standards](https://zkproof.org/) — Community reference

---

*This is an educational/portfolio project. The trusted setup is simulated locally and should never be used to secure real assets.*

---