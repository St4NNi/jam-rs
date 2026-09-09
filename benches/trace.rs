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
    let widths = [32; 128];
    let encoded = owner_postings::encode_block(&keys, &widths, |_, occurrence| {
        Ok(u64::from(occurrence.local_contig) * 1024 + occurrence.position)
    })
    .unwrap();
    assert_eq!(
        owner_postings::decode_block(
            &encoded.hot,
            &encoded.cold,
            &widths,
            64 * 1024 * 1024,
            |_, locus, canonical_orientation| {
                Ok(OwnerOccurrence {
                    local_contig: (locus / 1024) as u32,
                    position: locus % 1024,
                    canonical_orientation,
                })
            },
        )
        .unwrap(),
        keys
    );
    let hot = owner_postings::parse_hot(&encoded.hot).unwrap();
    let selected = owner_postings::find_key(&hot, keys[0].key).unwrap();
    let member = owner_postings::key_members(&hot, selected).unwrap()[0];
    let window = owner_postings::locate_member(&hot, member.member_ordinal).unwrap();
    let start = window.start_bit / 8;
    let end = window.end_bit.div_ceil(8);
    let tiny_key = owner_postings::find_key(&hot, keys[1].key).unwrap();
    let tiny_member = owner_postings::key_members(&hot, tiny_key).unwrap()[0];
    let tiny_window = owner_postings::locate_member(&hot, tiny_member.member_ordinal).unwrap();
    let tiny_start = tiny_window.start_bit / 8;
    let tiny_end = tiny_window.end_bit.div_ceil(8);
    let mut group = criterion.benchmark_group("owner_postings");
    group.sample_size(20);
    group.warm_up_time(Duration::from_millis(250));
    group.measurement_time(Duration::from_secs(1));
    group.throughput(Throughput::Bytes(
        (encoded.hot.len() + encoded.cold.len()) as u64,
    ));
    group.bench_function("encode_mixed", |bencher| {
        bencher.iter(|| {
            owner_postings::encode_block(black_box(&keys), black_box(&widths), |_, occurrence| {
                Ok(u64::from(occurrence.local_contig) * 1024 + occurrence.position)
            })
            .unwrap()
        })
    });
    group.bench_function("decode_mixed", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_block(
                black_box(&encoded.hot),
                black_box(&encoded.cold),
                black_box(&widths),
                64 * 1024 * 1024,
                |_, locus, canonical_orientation| {
                    Ok(OwnerOccurrence {
                        local_contig: (locus / 1024) as u32,
                        position: locus % 1024,
                        canonical_orientation,
                    })
                },
            )
            .unwrap()
        })
    });
    group.throughput(Throughput::Bytes(encoded.hot.len() as u64));
    group.bench_function("hot_directory", |bencher| {
        bencher.iter(|| owner_postings::parse_hot(black_box(&encoded.hot)).unwrap())
    });
    group.bench_function("absent_parse_and_find", |bencher| {
        bencher.iter(|| {
            let parsed = owner_postings::parse_hot(black_box(&encoded.hot)).unwrap();
            owner_postings::find_key(&parsed, black_box(key + 1))
                .map(|entry| entry.document_frequency)
        })
    });
    group.bench_function("absent_count", |bencher| {
        bencher.iter(|| owner_postings::find_key(black_box(&hot), black_box(key + 1)))
    });
    group.throughput(Throughput::Elements(member.occurrence_count));
    group.bench_function("selected_member", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_member_window(
                black_box(&encoded.cold[start as usize..end as usize]),
                start,
                window,
                black_box(&hot),
                black_box(&widths),
                64 * 1024 * 1024,
            )
            .unwrap()
        })
    });
    group.throughput(Throughput::Elements(tiny_member.occurrence_count));
    group.bench_function("selected_tiny_member", |bencher| {
        bencher.iter(|| {
            owner_postings::decode_member_window(
                black_box(&encoded.cold[tiny_start as usize..tiny_end as usize]),
                tiny_start,
                tiny_window,
                black_box(&hot),
                black_box(&widths),
                64 * 1024 * 1024,
            )
            .unwrap()
        })
    });
    group.finish();
}

criterion_group!(benches, affine_alignment, owner_postings);
criterion_main!(benches);
