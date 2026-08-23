//! Component measurements for the positional trace workflow.
//!
//! This target measures deterministic in-process operations. End-to-end and
//! resource measurements live under `evaluation/trace/`.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use jam_rs::jma::sequence::{decode_range, decode_sequence_blocks, encode_sequence_blocks};
use jam_rs::jma::sequence_builder::{SequenceBuildConfig, build_sequence_blocks};
use jam_rs::jma::writer::{decode_seed_section, encode_seed_section};
use jam_rs::jma::{ContigMetadata, SeedOccurrence, SequenceRange};
use jam_rs::resource::ResourceMetrics;
use jam_rs::trace::alignment::{AlignmentOptions, AlignmentWorkspace};
use jam_rs::trace::anchors::{Anchor, SeedOccurrenceGroup, generate_anchors};
use jam_rs::trace::chain::{ChainConfig, chain_anchors};
use jam_rs::trace::config::{AlignmentScoring, SeedSensitivity, SensitivityConfig};
use jam_rs::trace::model::{
    BaseAlignment, BaseInterval, CandidateResult, CoverageSummary, EditOperation, EditRun,
    InputResource, Strand, TraceMetagenomeResult, TraceRunFooter, TraceRunHeader, TraceStatus,
};
use jam_rs::trace::output::TraceJsonlWriter;
use jam_rs::trace::seeds::extract_seed_level;

fn dna(length: usize) -> Vec<u8> {
    let mut state = 0xD1B54A32D192ED03u64;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[((state >> 62) & 3) as usize]
        })
        .collect()
}

fn bench_seed_generation(c: &mut Criterion) {
    let sequence = dna(50_000);
    let config = SeedSensitivity {
        k: 31,
        scale: 200,
        max_occurrences: 128,
    };
    let mut group = c.benchmark_group("trace_seed_generation");
    group.bench_function("k31_scale200_50kb", |bench| {
        bench.iter(|| extract_seed_level(black_box(&sequence), config).unwrap())
    });
    group.finish();
}

fn bench_bucket_decoding(c: &mut Criterion) {
    let sequence = dna(50_000);
    let input = jam_rs::jma::seed_builder::SeedInput {
        contig_id: 0,
        sequence: &sequence,
    };
    let section = jam_rs::jma::seed_builder::build_seed_section(&[input], 31, 1, None)
        .unwrap()
        .section;
    let encoded = encode_seed_section(&section).unwrap();
    let mut group = c.benchmark_group("trace_bucket_decoding");
    group.bench_function("k31_scale1", |bench| {
        bench.iter(|| decode_seed_section(black_box(&encoded)).unwrap())
    });
    group.finish();
}

fn bench_chaining(c: &mut Criterion) {
    let anchors = (0..500u64)
        .map(|index| Anchor {
            query_position: index * 90,
            target_position: index * 92,
            contig_id: 0,
            strand: Strand::Forward,
            k: 31,
            hash: index + 1,
            canonical_kmer: index + 1001,
            query_reverse: false,
            target_reverse: false,
        })
        .collect::<Vec<_>>();
    let config = ChainConfig {
        max_chains: 8,
        min_anchors: 3,
        max_predecessors: 256,
        max_query_gap: 1_000_000,
        max_target_gap: 1_000_000,
        gap_penalty: 1,
    };
    let mut group = c.benchmark_group("trace_anchor_chaining");
    group.bench_function("500_anchors", |bench| {
        bench.iter(|| chain_anchors(black_box(&anchors), 50_000, config).unwrap())
    });
    group.finish();
}

fn bench_alignment(c: &mut Criterion) {
    let query = dna(4_000);
    let mut target = query.clone();
    for (index, base) in target.iter_mut().enumerate() {
        if index % 173 == 0 {
            *base = match *base {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }
    let options = AlignmentOptions::new(AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 128,
    })
    .with_max_cells(4_000_000);
    let mut workspace = AlignmentWorkspace::new();
    let mut group = c.benchmark_group("trace_alignment");
    group.bench_function("4kb_band128_reused_workspace", |bench| {
        bench.iter(|| {
            workspace
                .align(black_box(&query), black_box(&target), options)
                .unwrap()
        })
    });
    group.finish();
}

fn bench_sequence_decoding(c: &mut Criterion) {
    let sequence = dna(100_000);
    let blocks =
        build_sequence_blocks(0, &sequence, SequenceBuildConfig { block_bases: 4096 }).unwrap();
    let encoded = encode_sequence_blocks(&blocks).unwrap();
    let contigs = [ContigMetadata {
        id: 0,
        name: "contig-0".to_string(),
        length: sequence.len() as u64,
    }];
    let range = SequenceRange::new(12_345, 67_890).unwrap();
    let mut group = c.benchmark_group("trace_sequence_decoding");
    group.bench_function("decode_blocks_and_55kb_range", |bench| {
        bench.iter(|| {
            let decoded = decode_sequence_blocks(black_box(&encoded)).unwrap();
            decode_range(&decoded, &contigs, 0, range).unwrap()
        })
    });
    group.finish();
}

fn example_header() -> TraceRunHeader {
    TraceRunHeader {
        schema_version: jam_rs::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: "bench-run".to_string(),
        jam_rs_version: env!("CARGO_PKG_VERSION").to_string(),
        source_commit: Some("bench".to_string()),
        started_at_utc: "2026-01-01T00:00:00Z".to_string(),
        command: vec!["jam".to_string(), "trace".to_string()],
        plasmid_id: "query".to_string(),
        plasmid_length: 4_000,
        sensitivity: SensitivityConfig::default(),
        inputs: vec![InputResource {
            role: "plasmid".to_string(),
            redacted_locator: "query.fasta".to_string(),
            sha256: None,
        }],
    }
}

fn example_result() -> TraceMetagenomeResult {
    TraceMetagenomeResult {
        schema_version: jam_rs::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: "bench-run".to_string(),
        plasmid_id: "query".to_string(),
        metagenome_id: "sample".to_string(),
        status: TraceStatus::Complete,
        candidate: Some(CandidateResult {
            metagenome_id: "sample".to_string(),
            shared_hashes: 100,
            plasmid_hashes: 120,
            metagenome_hashes: 300,
            plasmid_containment: 0.83,
            metagenome_containment: 0.33,
            rank: 1,
            score_mode: "uniform".to_string(),
            bias_weighted_plasmid_containment: None,
            uniform_hash_e_value: Some(0.01),
        }),
        alignments: vec![BaseAlignment {
            plasmid_id: "query".to_string(),
            metagenome_id: "sample".to_string(),
            contig_id: "contig".to_string(),
            strand: Strand::Forward,
            query_segments: vec![BaseInterval {
                start: 0,
                end: 4_000,
            }],
            target_interval: BaseInterval {
                start: 10,
                end: 4_010,
            },
            query_length: 4_000,
            target_length: 4_000,
            origin_crossing: false,
            score: 7_900,
            matches: 3_950,
            substitutions: 50,
            insertions: 0,
            deletions: 0,
            cigar: "3950=50X".to_string(),
            edit_script: vec![
                EditRun {
                    operation: EditOperation::Equal,
                    length: 3_950,
                },
                EditRun {
                    operation: EditOperation::Substitution,
                    length: 50,
                },
            ],
            chain_score: 100,
            primary: true,
        }],
        coverage: Some(CoverageSummary {
            plasmid_length: 4_000,
            supported_bases: 4_000,
            supported_fraction: 1.0,
            primary_intervals: vec![BaseInterval {
                start: 0,
                end: 4_000,
            }],
            secondary_intervals: Vec::new(),
            gaps: Vec::new(),
            largest_gap: 0,
        }),
        failures: Vec::new(),
        resource_metrics: ResourceMetrics::default(),
    }
}

fn bench_json_writing(c: &mut Criterion) {
    let header = example_header();
    let result = example_result();
    let footer = TraceRunFooter {
        schema_version: jam_rs::trace::TRACE_JSON_SCHEMA_VERSION.to_string(),
        run_id: "bench-run".to_string(),
        completed_at_utc: "2026-01-01T00:00:01Z".to_string(),
        metagenomes_total: 1,
        metagenomes_with_candidates: 1,
        metagenomes_aligned: 1,
        metagenomes_failed: 0,
        alignments_total: 1,
        resource_metrics: ResourceMetrics::default(),
    };
    let mut group = c.benchmark_group("trace_json_writing");
    group.bench_function("header_result_footer_jsonl", |bench| {
        bench.iter(|| {
            let mut writer = TraceJsonlWriter::new(Vec::with_capacity(8 * 1024));
            writer.write_header(black_box(&header)).unwrap();
            writer.write_metagenome_result(black_box(&result)).unwrap();
            writer.write_footer(black_box(&footer)).unwrap();
            black_box(writer.finish().unwrap())
        })
    });
    group.finish();
}

fn bench_anchor_generation(c: &mut Criterion) {
    let groups = (0..250u64)
        .map(|index| SeedOccurrenceGroup {
            seed: jam_rs::trace::seeds::QuerySeed {
                position: index * 100,
                hash: index + 1,
                canonical_kmer: index + 10_001,
                reverse: false,
            },
            k: 31,
            occurrences: vec![SeedOccurrence {
                contig_id: 0,
                position: index * 101,
                reverse: false,
            }],
        })
        .collect::<Vec<_>>();
    let mut group = c.benchmark_group("trace_bucket_anchor_generation");
    group.bench_function("250_exact_occurrence_groups", |bench| {
        bench.iter(|| generate_anchors(black_box(&groups), 128, 100_000))
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_seed_generation,
    bench_bucket_decoding,
    bench_chaining,
    bench_alignment,
    bench_sequence_decoding,
    bench_json_writing,
    bench_anchor_generation
);
criterion_main!(benches);
