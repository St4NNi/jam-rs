use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;

#[inline]
pub fn murmur3_old(kmer: &[u8]) -> u64 {
    murmurhash3::murmurhash3_x64_128(kmer, 42).0
}

#[inline]
pub fn murmur3_new(kmer: &[u8]) -> u64 {
    fastmurmur3::murmur3_x64_128(kmer, 42) as u64
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Hashes");
    group.warm_up_time(Duration::from_millis(100));
    group.measurement_time(Duration::from_millis(100));

    for x in u64::MAX - 20..u64::MAX {
        group.bench_with_input(format!("xxhash_{}", x), &x, |b, &x| {
            b.iter(|| jam_rs::hash_functions::xxhash3(&x.to_be_bytes()));
        });
        group.bench_with_input(format!("ahash_{}", x), &x, |b, &x| {
            b.iter(|| jam_rs::hash_functions::ahash(x));
        });
        group.bench_with_input(format!("murmur3_old_{}", x), &x, |b, &x| {
            b.iter(|| murmur3_old(&x.to_be_bytes()));
        });
        group.bench_with_input(format!("murmur3_new_{}", x), &x, |b, &x| {
            b.iter(|| murmur3_new(&x.to_be_bytes()));
        });
    }
    group.finish();
}
criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
