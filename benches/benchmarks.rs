use criterion::{Criterion, criterion_group, criterion_main};
use std::{hint::black_box, time::Duration};

#[inline]
pub fn murmur3_old(value: u64) -> u64 {
    murmurhash3::murmurhash3_x64_128(&value.to_le_bytes(), 42).0
}

#[inline]
pub fn murmur3_new(value: u64) -> u64 {
    fastmurmur3::murmur3_x64_128(&value.to_le_bytes(), 42) as u64
}

#[inline]
pub fn xxhash3(value: u64) -> u64 {
    xxhash_rust::xxh3::xxh3_64(&value.to_le_bytes())
}

#[inline]
pub fn jamhash(value: u64) -> u64 {
    jamhash::jamhash_u64(value)
}

fn bench_hash_functions(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_hash");
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

    // Pre-generate enough values for the benchmark
    let mut values = (0..100000u64).cycle();

    group.bench_function("xxhash3", |b| {
        b.iter(|| xxhash3(black_box(values.next().unwrap())))
    });

    group.bench_function("murmur3_old", |b| {
        b.iter(|| murmur3_old(black_box(values.next().unwrap())))
    });

    group.bench_function("murmur3_new", |b| {
        b.iter(|| murmur3_new(black_box(values.next().unwrap())))
    });

    group.bench_function("jamhash", |b| {
        b.iter(|| jamhash(black_box(values.next().unwrap())))
    });

    group.finish();
}

criterion_group!(benches, bench_hash_functions);
criterion_main!(benches);
