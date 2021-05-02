use criterion::{black_box, criterion_group, criterion_main, Criterion};
use reservoirs::prelude::*;

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("sim_73", |b| {
        b.iter(|| {
            Reservoir::new()
                .input(&0.73)?
                .output(&0.73)?
                .sim(black_box(&1000.0))
        })
    });
}

pub fn transit_benchmark(c: &mut Criterion) {
    c.bench_function("transit_times", |b| {
        b.iter(|| {
            Model::new(
                Reservoir::new()
                    .input(&0.73)
                    .unwrap()
                    .output(&0.73)
                    .unwrap(),
            )
            .period(&1000.0)
            .runs(100)
            .transit_times()
        })
    });
}

criterion_group!(benches, criterion_benchmark, transit_benchmark);
criterion_main!(benches);
