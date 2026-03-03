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

#[inline]
pub fn wang64(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key.wrapping_add(key << 31)
}

fn bench_hash_functions(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_hash");
    group.warm_up_time(Duration::from_millis(500));
    group.measurement_time(Duration::from_secs(2));

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

    group.bench_function("wang64", |b| {
        b.iter(|| wang64(black_box(values.next().unwrap())))
    });

    group.finish();
}

criterion_group!(benches, bench_hash_functions);
criterion_main!(benches);
