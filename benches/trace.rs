use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use jam_rs::alignment::{AlignmentConfig, AlignmentWorkspace};
use jam_rs::owner_postings;
use std::hint::black_box;
use std::time::Duration;

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

fn owner_postings(criterion: &mut Criterion) {
    use owner_postings::{OwnerKey, OwnerMember, OwnerOccurrence};

    let mut key = 0u64;
    let keys = (0..256)
        .map(|ordinal| {
            key += 1 + (ordinal * 65_537) % 1_000_003;
            OwnerKey {
                key,
                members: (0..if ordinal % 64 == 0 { 128 } else { 1 })
                    .map(|document_id| OwnerMember {
                        document_id,
                        occurrences: (0..if ordinal % 64 == 0 { 16 } else { 1 })
                            .map(|position| OwnerOccurrence {
                                local_contig: position / 4,
                                position: u64::from(position % 4) * 17 + 31,
                                canonical_orientation: position % 2 == 0,
                            })
                            .collect(),
                    })
                    .collect(),
            }
        })
        .collect::<Vec<_>>();
    let encoded = owner_postings::encode_block(&keys).unwrap();
    assert_eq!(
        owner_postings::decode_block(&encoded.hot, &encoded.cold).unwrap(),
        keys
    );
    let mut group = criterion.benchmark_group("owner_postings");
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    group.throughput(Throughput::Bytes(
        (encoded.hot.len() + encoded.cold.len()) as u64,
    ));
    group.bench_function("encode_mixed", |bencher| {
        bencher.iter(|| owner_postings::encode_block(black_box(&keys)).unwrap())
    });
    group.bench_function("decode_mixed", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_block(black_box(&encoded.hot), black_box(&encoded.cold)).unwrap()
        })
    });
    group.throughput(Throughput::Bytes(encoded.hot.len() as u64));
    group.bench_function("hot_directory", |bencher| {
        bencher.iter(|| owner_postings::parse_hot(black_box(&encoded.hot)).unwrap())
    });
    group.bench_function("absent_key", |bencher| {
        bencher.iter(|| {
            owner_postings::lookup_hot(black_box(&encoded.hot), black_box(key + 1)).unwrap()
        })
    });
    group.finish();
}

criterion_group!(benches, affine_alignment, owner_postings);
criterion_main!(benches);
