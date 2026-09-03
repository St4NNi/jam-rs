use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use jam_rs::alignment::{AlignmentConfig, AlignmentWorkspace};
use std::hint::black_box;

fn sequence(length: usize) -> Vec<u8> {
    let mut state = 7u64;
    (0..length)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn affine_alignment(criterion: &mut Criterion) {
    let query = sequence(5_000);
    let mut target = query.clone();
    for position in (50..target.len()).step_by(100) {
        target[position] = match target[position] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let config = AlignmentConfig::default();
    let mut workspace = AlignmentWorkspace::default();
    assert!(workspace.align(&query, &target, config).unwrap().identity() > 0.98);

    let mut group = criterion.benchmark_group("trace_alignment");
    group.throughput(Throughput::Bytes(query.len() as u64));
    group.bench_function("5kb_99pct", |bencher| {
        bencher.iter(|| {
            workspace
                .align(black_box(&query), black_box(&target), config)
                .unwrap()
        })
    });
    group.finish();
}

criterion_group!(benches, affine_alignment);
criterion_main!(benches);
