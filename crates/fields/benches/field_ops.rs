use criterion::{criterion_group, criterion_main, Criterion};
use zkrust_fields::{FieldElement, Fp, Fp12, Fp2, Fr};

fn bench_fp_mul(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a = Fp::random(&mut rng);
    let b = Fp::random(&mut rng);
    c.bench_function("fp_mul", |bench| bench.iter(|| a * b));
}

fn bench_fp_inv(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a = Fp::random(&mut rng);
    c.bench_function("fp_inv", |bench| bench.iter(|| a.inv()));
}

fn bench_fr_mul(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a = Fr::random(&mut rng);
    let b = Fr::random(&mut rng);
    c.bench_function("fr_mul", |bench| bench.iter(|| a * b));
}

fn bench_fp2_mul(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a = Fp2::random(&mut rng);
    let b = Fp2::random(&mut rng);
    c.bench_function("fp2_mul", |bench| bench.iter(|| a * b));
}

fn bench_fp12_mul(c: &mut Criterion) {
    let mut rng = rand::thread_rng();
    let a = Fp12::random(&mut rng);
    let b = Fp12::random(&mut rng);
    c.bench_function("fp12_mul", |bench| bench.iter(|| a * b));
}

criterion_group!(
    benches,
    bench_fp_mul,
    bench_fp_inv,
    bench_fr_mul,
    bench_fp2_mul,
    bench_fp12_mul,
);
criterion_main!(benches);
