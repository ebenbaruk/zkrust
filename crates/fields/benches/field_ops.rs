use criterion::{criterion_group, criterion_main, Criterion};

fn field_benchmarks(_c: &mut Criterion) {
    // TODO: Add field operation benchmarks
}

criterion_group!(benches, field_benchmarks);
criterion_main!(benches);
