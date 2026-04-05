use criterion::{criterion_group, criterion_main, Criterion};
use zkrust_fields::{FieldElement, Fr};
use zkrust_r1cs::{ConstraintSystem, LinearCombination};

fn squaring_circuit() -> ConstraintSystem {
    let mut cs = ConstraintSystem::new();
    let x = cs.alloc_public_input(Fr::from(3u64));
    let y = cs.alloc_public_input(Fr::from(9u64));
    cs.enforce(
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(x),
        LinearCombination::from_variable(y),
    );
    cs
}

fn bench_setup(c: &mut Criterion) {
    let cs = squaring_circuit();
    let mut rng = rand::thread_rng();

    c.bench_function("groth16_setup_square", |b| {
        b.iter(|| {
            let tau = Fr::random(&mut rng);
            let alpha = Fr::random(&mut rng);
            let beta = Fr::random(&mut rng);
            let gamma = Fr::random(&mut rng);
            let delta = Fr::random(&mut rng);
            zkrust_groth16::setup(&cs, tau, alpha, beta, gamma, delta)
        })
    });
}

fn bench_prove(c: &mut Criterion) {
    let cs = squaring_circuit();
    let mut rng = rand::thread_rng();
    let tau = Fr::random(&mut rng);
    let alpha = Fr::random(&mut rng);
    let beta = Fr::random(&mut rng);
    let gamma = Fr::random(&mut rng);
    let delta = Fr::random(&mut rng);
    let pk = zkrust_groth16::setup(&cs, tau, alpha, beta, gamma, delta);

    c.bench_function("groth16_prove_square", |b| {
        b.iter(|| zkrust_groth16::prove(&pk, &cs, &mut rng))
    });
}

fn bench_verify(c: &mut Criterion) {
    let cs = squaring_circuit();
    let mut rng = rand::thread_rng();
    let tau = Fr::random(&mut rng);
    let alpha = Fr::random(&mut rng);
    let beta = Fr::random(&mut rng);
    let gamma = Fr::random(&mut rng);
    let delta = Fr::random(&mut rng);
    let pk = zkrust_groth16::setup(&cs, tau, alpha, beta, gamma, delta);
    let proof = zkrust_groth16::prove(&pk, &cs, &mut rng);
    let public_inputs = vec![Fr::from(3u64), Fr::from(9u64)];

    c.bench_function("groth16_verify_square", |b| {
        b.iter(|| zkrust_groth16::verify(&pk.vk, &public_inputs, &proof))
    });
}

criterion_group!(benches, bench_setup, bench_prove, bench_verify);
criterion_main!(benches);
