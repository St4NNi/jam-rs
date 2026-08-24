//! Bounded native sequence/archive policy measurements.
//!
//! The benchmark keeps policy selection separate from end-to-end query
//! timing.  It covers the required fixed 16/64/256 KiB and 1 MiB blocks,
//! independently compressed two-bit payloads, and the three DNA Gear table
//! families.  Archive build/size, RSS, and remote range measurements are
//! collected by the manifest runners under `evaluation/archive-seeds/`.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use jam_rs::jma::sequence_builder::{IndexedSequenceBuildConfig, build_indexed_sequence_blocks};
use jam_rs::sequence::{
    BlockCodec, GearTable, SequenceBlockPolicy, decode_block_range, decode_stored_block,
};
use std::hint::black_box;

fn synthetic_dna(length: usize) -> Vec<u8> {
    let mut state = 0x8f31_4a72_d0ce_119bu64;
    (0..length)
        .map(|index| {
            state = state
                .wrapping_mul(2_862_933_555_777_941_757)
                .wrapping_add(3_039)
                .rotate_left(17);
            if index % 9_973 == 0 {
                b'N'
            } else {
                b"ACGT"[((state >> 32) & 3) as usize]
            }
        })
        .collect()
}

fn build_case(
    sequence: &[u8],
    policy: SequenceBlockPolicy,
    codec: BlockCodec,
) -> (usize, usize, usize) {
    let blocks =
        build_indexed_sequence_blocks(0, sequence, IndexedSequenceBuildConfig { policy, codec })
            .expect("benchmark policy must encode");
    let stored_bytes = blocks.iter().map(|block| block.payload.len()).sum();
    let ambiguity_bytes = blocks
        .iter()
        .map(|block| block.ambiguity_payload.len())
        .sum();
    black_box((blocks.len(), stored_bytes, ambiguity_bytes))
}

fn bench_build_policies(c: &mut Criterion) {
    let sequence = synthetic_dna(1_200_000);
    let mut group = c.benchmark_group("archive_sequence_policy_build");
    for block_bases in [16 * 1024, 64 * 1024, 256 * 1024, 1024 * 1024] {
        for codec in [BlockCodec::Raw2Bit, BlockCodec::Zstd2Bit] {
            let codec_name = match codec {
                BlockCodec::Raw2Bit => "raw2bit",
                BlockCodec::Zstd2Bit => "zstd2bit",
            };
            let name = format!("fixed-{}-{}", block_bases / 1024, codec_name);
            group.bench_with_input(
                BenchmarkId::from_parameter(name),
                &(block_bases, codec),
                |bench, &(block_bases, codec)| {
                    bench.iter(|| {
                        build_case(
                            black_box(&sequence),
                            SequenceBlockPolicy::Fixed { block_bases },
                            codec,
                        )
                    })
                },
            );
        }
    }
    for table in [
        GearTable::SingleBase,
        GearTable::Dinucleotide,
        GearTable::PackedFourBase,
    ] {
        let table_name = match table {
            GearTable::SingleBase => "single-base",
            GearTable::Dinucleotide => "dinucleotide",
            GearTable::PackedFourBase => "packed-four-base",
        };
        for &(min_bases, target_bases, max_bases) in &[
            (16 * 1024, 64 * 1024, 256 * 1024),
            (32 * 1024, 128 * 1024, 512 * 1024),
        ] {
            let name = format!(
                "gear-{}-{}-{}-{}",
                table_name,
                min_bases / 1024,
                target_bases / 1024,
                max_bases / 1024
            );
            group.bench_with_input(
                BenchmarkId::from_parameter(name),
                &(min_bases, target_bases, max_bases, table),
                |bench, &(min_bases, target_bases, max_bases, table)| {
                    bench.iter(|| {
                        build_case(
                            black_box(&sequence),
                            SequenceBlockPolicy::Gear {
                                min_bases,
                                target_bases,
                                max_bases,
                                table,
                            },
                            BlockCodec::Zstd2Bit,
                        )
                    })
                },
            );
        }
    }
    group.finish();
}

fn bench_selected_range_reads(c: &mut Criterion) {
    let sequence = synthetic_dna(1_200_000);
    let mut group = c.benchmark_group("archive_sequence_selective_reads");
    for (name, policy, codec) in [
        (
            "raw-fixed-64k",
            SequenceBlockPolicy::Fixed {
                block_bases: 64 * 1024,
            },
            BlockCodec::Raw2Bit,
        ),
        (
            "zstd-fixed-256k",
            SequenceBlockPolicy::Fixed {
                block_bases: 256 * 1024,
            },
            BlockCodec::Zstd2Bit,
        ),
        (
            "zstd-gear-64k",
            SequenceBlockPolicy::Gear {
                min_bases: 16 * 1024,
                target_bases: 64 * 1024,
                max_bases: 256 * 1024,
                table: GearTable::SingleBase,
            },
            BlockCodec::Zstd2Bit,
        ),
    ] {
        let blocks = build_indexed_sequence_blocks(
            0,
            &sequence,
            IndexedSequenceBuildConfig { policy, codec },
        )
        .expect("benchmark policy must encode");
        let request = 487_321u64..532_901u64;
        group.bench_function(name, |bench| {
            bench.iter(|| {
                let decoded = blocks
                    .iter()
                    .filter_map(|block| {
                        let start = request.start.max(block.record.base_start);
                        let end = request.end.min(block.record.base_end()?);
                        (start < end).then(|| {
                            decode_block_range(
                                &block.record,
                                &block.payload,
                                &block.ambiguity_payload,
                                start..end,
                            )
                            .expect("selected range decodes")
                        })
                    })
                    .flatten()
                    .collect::<Vec<_>>();
                black_box(decoded)
            })
        });
    }
    group.finish();
}

fn bench_full_block_decode(c: &mut Criterion) {
    let sequence = synthetic_dna(512 * 1024);
    let mut group = c.benchmark_group("archive_sequence_block_decode");
    for codec in [BlockCodec::Raw2Bit, BlockCodec::Zstd2Bit] {
        let blocks = build_indexed_sequence_blocks(
            0,
            &sequence,
            IndexedSequenceBuildConfig {
                policy: SequenceBlockPolicy::Fixed {
                    block_bases: 64 * 1024,
                },
                codec,
            },
        )
        .expect("benchmark policy must encode");
        let name = match codec {
            BlockCodec::Raw2Bit => "raw2bit",
            BlockCodec::Zstd2Bit => "zstd2bit",
        };
        group.bench_function(name, |bench| {
            bench.iter(|| {
                let bases = blocks
                    .iter()
                    .flat_map(|block| {
                        decode_stored_block(&block.record, &block.payload, &block.ambiguity_payload)
                            .expect("block decodes")
                    })
                    .collect::<Vec<_>>();
                black_box(bases)
            })
        });
    }
    group.finish();
}

criterion_group!(
    archive_benches,
    bench_build_policies,
    bench_selected_range_reads,
    bench_full_block_decode
);
criterion_main!(archive_benches);
