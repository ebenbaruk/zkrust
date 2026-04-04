use criterion::{criterion_group, criterion_main, Criterion};

fn curve_benchmarks(_c: &mut Criterion) {
    // TODO: Add curve operation benchmarks
}

criterion_group!(benches, curve_benchmarks);
criterion_main!(benches);
