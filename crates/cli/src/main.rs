use clap::{Parser, Subcommand};
use std::fs;
use std::path::{Path, PathBuf};
use zkrust_fields::{FieldElement, Fr};

mod circuits;

#[derive(Parser)]
#[command(name = "zkrust", about = "Groth16 zk-SNARK proof system over BN254")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate proving and verifying keys for a circuit
    Setup {
        /// Circuit name (e.g., "square", "range32")
        #[arg(short, long)]
        circuit: String,

        /// Output directory for keys
        #[arg(short, long, default_value = "./keys")]
        output: PathBuf,
    },

    /// Generate a proof
    Prove {
        /// Circuit name
        #[arg(short, long)]
        circuit: String,

        /// Witness values as comma-separated integers (public first, then private)
        #[arg(short, long)]
        witness: String,

        /// Directory containing proving key
        #[arg(short, long, default_value = "./keys")]
        keys: PathBuf,

        /// Output file for proof
        #[arg(short, long, default_value = "proof.bin")]
        output: PathBuf,
    },

    /// Verify a proof
    Verify {
        /// Public inputs as comma-separated integers
        #[arg(short, long)]
        public_inputs: String,

        /// Directory containing verifying key
        #[arg(short, long, default_value = "./keys")]
        keys: PathBuf,

        /// Proof file
        #[arg(short = 'f', long, default_value = "proof.bin")]
        proof: PathBuf,
    },

    /// List available circuits
    List,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Setup { circuit, output } => cmd_setup(&circuit, &output),
        Commands::Prove {
            circuit,
            witness,
            keys,
            output,
        } => cmd_prove(&circuit, &witness, &keys, &output),
        Commands::Verify {
            public_inputs,
            keys,
            proof,
        } => cmd_verify(&public_inputs, &keys, &proof),
        Commands::List => cmd_list(),
    }
}

fn cmd_setup(circuit_name: &str, output_dir: &Path) {
    let cs = circuits::build_circuit(circuit_name, &[]);
    if cs.is_none() {
        eprintln!("Unknown circuit: {circuit_name}. Use 'zkrust list' to see available circuits.");
        std::process::exit(1);
    }
    let cs = cs.unwrap();

    println!(
        "Circuit '{circuit_name}': {} constraints, {} variables ({} public, {} private)",
        cs.num_constraints(),
        cs.num_variables(),
        cs.num_public_inputs(),
        cs.num_private_inputs()
    );

    // Generate toxic waste (INSECURE — for demonstration only)
    let mut rng = rand::thread_rng();
    let tau = Fr::random(&mut rng);
    let alpha = Fr::random(&mut rng);
    let beta = Fr::random(&mut rng);
    let gamma = Fr::random(&mut rng);
    let delta = Fr::random(&mut rng);

    let pk = zkrust_groth16::setup(&cs, tau, alpha, beta, gamma, delta);

    fs::create_dir_all(output_dir).expect("failed to create output directory");

    let pk_bytes = bincode::serialize(&pk_to_serializable(&pk)).expect("failed to serialize pk");
    let vk_bytes = bincode::serialize(&vk_to_serializable(&pk.vk)).expect("failed to serialize vk");

    fs::write(output_dir.join("proving_key.bin"), &pk_bytes).expect("failed to write proving key");
    fs::write(output_dir.join("verifying_key.bin"), &vk_bytes)
        .expect("failed to write verifying key");

    println!(
        "Keys written to {}/ ({} bytes PK, {} bytes VK)",
        output_dir.display(),
        pk_bytes.len(),
        vk_bytes.len()
    );
}

fn cmd_prove(circuit_name: &str, witness_str: &str, keys_dir: &Path, output: &Path) {
    let values: Vec<u64> = witness_str
        .split(',')
        .map(|s| s.trim().parse().expect("invalid witness value"))
        .collect();

    let cs = circuits::build_circuit(circuit_name, &values);
    if cs.is_none() {
        eprintln!("Unknown circuit: {circuit_name}");
        std::process::exit(1);
    }
    let cs = cs.unwrap();

    if !cs.is_satisfied() {
        eprintln!("Error: witness does not satisfy the circuit constraints");
        std::process::exit(1);
    }

    let pk_bytes = fs::read(keys_dir.join("proving_key.bin")).expect("failed to read proving key");
    let pk_ser: SerializableProvingKey =
        bincode::deserialize(&pk_bytes).expect("failed to deserialize proving key");
    let pk = pk_from_serializable(&pk_ser);

    let mut rng = rand::thread_rng();
    let proof = zkrust_groth16::prove(&pk, &cs, &mut rng);

    let proof_bytes =
        bincode::serialize(&proof_to_serializable(&proof)).expect("failed to serialize proof");
    fs::write(output, &proof_bytes).expect("failed to write proof");

    println!(
        "Proof written to {} ({} bytes)",
        output.display(),
        proof_bytes.len()
    );
}

fn cmd_verify(inputs_str: &str, keys_dir: &Path, proof_path: &Path) {
    let public_inputs: Vec<Fr> = inputs_str
        .split(',')
        .map(|s| Fr::from(s.trim().parse::<u64>().expect("invalid public input")))
        .collect();

    let vk_bytes =
        fs::read(keys_dir.join("verifying_key.bin")).expect("failed to read verifying key");
    let vk_ser: SerializableVerifyingKey =
        bincode::deserialize(&vk_bytes).expect("failed to deserialize verifying key");
    let vk = vk_from_serializable(&vk_ser);

    let proof_bytes = fs::read(proof_path).expect("failed to read proof");
    let proof_ser: SerializableProof =
        bincode::deserialize(&proof_bytes).expect("failed to deserialize proof");
    let proof = proof_from_serializable(&proof_ser);

    if zkrust_groth16::verify(&vk, &public_inputs, &proof) {
        println!("Proof is VALID");
    } else {
        println!("Proof is INVALID");
        std::process::exit(1);
    }
}

fn cmd_list() {
    println!("Available circuits:");
    println!("  square    - Prove x * x = y (2 public inputs: x, y)");
    println!("  range32   - Prove x < 2^32 (1 public input: x)");
    println!("  cubic     - Prove x^3 + x + 5 = y (2 public inputs: x, y)");
}

// --- Serialization helpers ---
// We serialize field/curve elements as raw u64 arrays since serde
// isn't derived on the core types.

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
struct SerializableG1([u64; 4], [u64; 4], bool);

#[derive(Serialize, Deserialize)]
struct SerializableG2([u64; 4], [u64; 4], [u64; 4], [u64; 4], bool);

#[derive(Serialize, Deserialize)]
struct SerializableProof {
    a: SerializableG1,
    b: SerializableG2,
    c: SerializableG1,
}

#[derive(Serialize, Deserialize)]
struct SerializableVerifyingKey {
    alpha_g1: SerializableG1,
    beta_g2: SerializableG2,
    gamma_g2: SerializableG2,
    delta_g2: SerializableG2,
    ic: Vec<SerializableG1>,
}

#[derive(Serialize, Deserialize)]
struct SerializableProvingKey {
    alpha_g1: SerializableG1,
    beta_g1: SerializableG1,
    beta_g2: SerializableG2,
    delta_g1: SerializableG1,
    delta_g2: SerializableG2,
    a_g1: Vec<SerializableG1>,
    b_g1: Vec<SerializableG1>,
    b_g2: Vec<SerializableG2>,
    l_g1: Vec<SerializableG1>,
    h_g1: Vec<SerializableG1>,
    vk: SerializableVerifyingKey,
    num_public: usize,
    domain_size: usize,
}

use zkrust_curves::{G1Affine, G2Affine};
use zkrust_fields::{Fp, Fp2};

fn g1_to_ser(p: &G1Affine) -> SerializableG1 {
    SerializableG1(p.x.to_raw(), p.y.to_raw(), p.infinity)
}

fn g1_from_ser(s: &SerializableG1) -> G1Affine {
    G1Affine {
        x: Fp::from_raw(s.0),
        y: Fp::from_raw(s.1),
        infinity: s.2,
    }
}

fn g2_to_ser(p: &G2Affine) -> SerializableG2 {
    SerializableG2(
        p.x.c0.to_raw(),
        p.x.c1.to_raw(),
        p.y.c0.to_raw(),
        p.y.c1.to_raw(),
        p.infinity,
    )
}

fn g2_from_ser(s: &SerializableG2) -> G2Affine {
    G2Affine {
        x: Fp2::new(Fp::from_raw(s.0), Fp::from_raw(s.1)),
        y: Fp2::new(Fp::from_raw(s.2), Fp::from_raw(s.3)),
        infinity: s.4,
    }
}

fn proof_to_serializable(p: &zkrust_groth16::Proof) -> SerializableProof {
    SerializableProof {
        a: g1_to_ser(&p.a),
        b: g2_to_ser(&p.b),
        c: g1_to_ser(&p.c),
    }
}

fn proof_from_serializable(s: &SerializableProof) -> zkrust_groth16::Proof {
    zkrust_groth16::Proof {
        a: g1_from_ser(&s.a),
        b: g2_from_ser(&s.b),
        c: g1_from_ser(&s.c),
    }
}

fn vk_to_serializable(vk: &zkrust_groth16::VerifyingKey) -> SerializableVerifyingKey {
    SerializableVerifyingKey {
        alpha_g1: g1_to_ser(&vk.alpha_g1),
        beta_g2: g2_to_ser(&vk.beta_g2),
        gamma_g2: g2_to_ser(&vk.gamma_g2),
        delta_g2: g2_to_ser(&vk.delta_g2),
        ic: vk.ic.iter().map(g1_to_ser).collect(),
    }
}

fn vk_from_serializable(s: &SerializableVerifyingKey) -> zkrust_groth16::VerifyingKey {
    zkrust_groth16::VerifyingKey {
        alpha_g1: g1_from_ser(&s.alpha_g1),
        beta_g2: g2_from_ser(&s.beta_g2),
        gamma_g2: g2_from_ser(&s.gamma_g2),
        delta_g2: g2_from_ser(&s.delta_g2),
        ic: s.ic.iter().map(g1_from_ser).collect(),
    }
}

fn pk_to_serializable(pk: &zkrust_groth16::ProvingKey) -> SerializableProvingKey {
    SerializableProvingKey {
        alpha_g1: g1_to_ser(&pk.alpha_g1),
        beta_g1: g1_to_ser(&pk.beta_g1),
        beta_g2: g2_to_ser(&pk.beta_g2),
        delta_g1: g1_to_ser(&pk.delta_g1),
        delta_g2: g2_to_ser(&pk.delta_g2),
        a_g1: pk.a_g1.iter().map(g1_to_ser).collect(),
        b_g1: pk.b_g1.iter().map(g1_to_ser).collect(),
        b_g2: pk.b_g2.iter().map(g2_to_ser).collect(),
        l_g1: pk.l_g1.iter().map(g1_to_ser).collect(),
        h_g1: pk.h_g1.iter().map(g1_to_ser).collect(),
        vk: vk_to_serializable(&pk.vk),
        num_public: pk.num_public,
        domain_size: pk.domain_size,
    }
}

fn pk_from_serializable(s: &SerializableProvingKey) -> zkrust_groth16::ProvingKey {
    zkrust_groth16::ProvingKey {
        alpha_g1: g1_from_ser(&s.alpha_g1),
        beta_g1: g1_from_ser(&s.beta_g1),
        beta_g2: g2_from_ser(&s.beta_g2),
        delta_g1: g1_from_ser(&s.delta_g1),
        delta_g2: g2_from_ser(&s.delta_g2),
        a_g1: s.a_g1.iter().map(g1_from_ser).collect(),
        b_g1: s.b_g1.iter().map(g1_from_ser).collect(),
        b_g2: s.b_g2.iter().map(g2_from_ser).collect(),
        l_g1: s.l_g1.iter().map(g1_from_ser).collect(),
        h_g1: s.h_g1.iter().map(g1_from_ser).collect(),
        vk: vk_from_serializable(&s.vk),
        num_public: s.num_public,
        domain_size: s.domain_size,
    }
}
