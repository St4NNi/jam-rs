//! Bounded component smoke measurements for the positional trace path.
//!
//! This target registers one labelled representative case per component so
//! `cargo test --all-targets` and explicit Criterion runs remain bounded.
//! The complete thread/profile/query/candidate/block matrix is executed only
//! by the Python harnesses under `evaluation/trace/`; process RSS, child CPU
//! time, and remote request counters are collected there.

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use jam_rs::jma::seed_builder::{SeedInput, build_seed_section};
use jam_rs::jma::sequence::{decode_range, decode_sequence_blocks, encode_sequence_blocks};
use jam_rs::jma::sequence_builder::{SequenceBuildConfig, build_sequence_blocks};
use jam_rs::jma::writer::{decode_seed_section, encode_seed_section};
use jam_rs::jma::{ContigMetadata, SeedOccurrence, SequenceRange};
use jam_rs::trace::alignment::{AlignmentOptions, AlignmentWorkspace};
use jam_rs::trace::anchors::{SeedOccurrenceGroup, generate_anchors};
use jam_rs::trace::chain::{ChainConfig, chain_anchors};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use jam_rs::trace::seeds::extract_seed_level;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
enum QueryShape {
    Small,
}

#[derive(Clone, Copy, Debug)]
struct MemoryCase {
    threads: usize,
    profile: SensitivityProfile,
    shape: QueryShape,
    candidates: usize,
    cache_block_bytes: usize,
}

#[derive(Clone, Copy)]
enum Component {
    Seed,
    Bucket,
    Anchor,
    Chain,
    Sequence,
    Alignment,
    Output,
}

const COMPONENTS: [Component; 7] = [
    Component::Seed,
    Component::Bucket,
    Component::Anchor,
    Component::Chain,
    Component::Sequence,
    Component::Alignment,
    Component::Output,
];

fn shape_name(shape: QueryShape) -> &'static str {
    match shape {
        QueryShape::Small => "small",
    }
}

fn profile_name(profile: SensitivityProfile) -> &'static str {
    match profile {
        SensitivityProfile::Fast => "fast",
        SensitivityProfile::Balanced => "balanced",
        SensitivityProfile::Sensitive => "sensitive",
    }
}

fn component_name(component: Component) -> &'static str {
    match component {
        Component::Seed => "seed",
        Component::Bucket => "bucket",
        Component::Anchor => "anchor",
        Component::Chain => "chain",
        Component::Sequence => "sequence",
        Component::Alignment => "alignment",
        Component::Output => "output",
    }
}

fn query_sequence(shape: QueryShape) -> Vec<u8> {
    let length = match shape {
        QueryShape::Small => 4_096,
    };
    let mut state = 0xD1B5_4A32_D192_ED03_u64;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[((state >> 62) & 3) as usize]
        })
        .collect()
}

fn for_case(case: MemoryCase, component: Component) -> usize {
    let query = query_sequence(case.shape);
    let target = query.clone();
    let sensitivity = SensitivityConfig::for_profile(case.profile);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(case.threads)
        .build()
        .expect("measurement thread pool must build");
    pool.install(|| {
        (0..case.candidates)
            .into_par_iter()
            .map(|candidate| {
                component_work(
                    &query,
                    &target,
                    sensitivity.clone(),
                    case,
                    component,
                    candidate,
                )
            })
            .sum()
    })
}

fn component_work(
    query: &[u8],
    target: &[u8],
    sensitivity: SensitivityConfig,
    case: MemoryCase,
    component: Component,
    candidate: usize,
) -> usize {
    let seed_level = extract_seed_level(query, sensitivity.primary).expect("fixed seed config");
    if matches!(component, Component::Seed) {
        return black_box(seed_level.seeds.len());
    }

    let input = [SeedInput {
        contig_id: 0,
        sequence: target,
    }];
    if matches!(component, Component::Bucket) {
        let section = build_seed_section(
            &input,
            sensitivity.primary.k,
            sensitivity.primary.scale,
            None,
        )
        .expect("fixed bucket config")
        .section;
        let encoded = encode_seed_section(&section).expect("seed section encodes");
        return black_box(
            decode_seed_section(&encoded)
                .expect("seed section decodes")
                .levels
                .len(),
        );
    }

    let groups = seed_level
        .seeds
        .iter()
        .take(512)
        .map(|seed| SeedOccurrenceGroup {
            seed: *seed,
            k: seed_level.k,
            occurrences: vec![SeedOccurrence {
                contig_id: 0,
                position: seed.position,
                reverse: seed.reverse,
            }],
        })
        .collect::<Vec<_>>();
    let anchors = generate_anchors(
        &groups,
        sensitivity.primary.max_occurrences,
        sensitivity.max_anchors_per_candidate,
    );
    if matches!(component, Component::Anchor) {
        return black_box(anchors.anchors.len());
    }

    let chains = chain_anchors(
        &anchors.anchors,
        query.len() as u64,
        ChainConfig::from_sensitivity(&sensitivity),
    )
    .expect("fixed chain config");
    if matches!(component, Component::Chain) {
        return black_box(chains.len());
    }

    let block_bases = case.cache_block_bytes.min(target.len().max(1));
    let blocks = build_sequence_blocks(0, target, SequenceBuildConfig { block_bases })
        .expect("fixed sequence block config");
    let encoded_blocks = encode_sequence_blocks(&blocks).expect("sequence blocks encode");
    let decoded_blocks = decode_sequence_blocks(&encoded_blocks).expect("sequence blocks decode");
    if matches!(component, Component::Sequence) {
        let contigs = [ContigMetadata {
            id: 0,
            name: "target".to_string(),
            length: target.len() as u64,
        }];
        let range = SequenceRange::new(0, target.len().min(8_192) as u64).expect("range");
        let decoded = decode_range(&decoded_blocks, &contigs, 0, range).expect("range decodes");
        return black_box(decoded.len() + decoded_blocks.len());
    }

    let alignment_length = query.len().min(target.len()).min(8_192);
    let options = AlignmentOptions::new(sensitivity.alignment)
        .with_max_cells((alignment_length + 1).saturating_mul(513));
    if matches!(component, Component::Alignment) {
        let mut workspace = AlignmentWorkspace::new();
        let result = workspace
            .align(
                &query[..alignment_length],
                &target[..alignment_length],
                options,
            )
            .expect("identical measurement sequences align");
        return black_box(result.matches as usize + workspace.capacity_cells());
    }

    let record = serde_json::json!({
        "component": "output",
        "candidate": candidate,
        "profile": profile_name(case.profile),
        "threads": case.threads,
        "shape": shape_name(case.shape),
        "seed_count": seed_level.seeds.len(),
        "anchor_count": anchors.anchors.len(),
        "chain_count": chains.len(),
        "sequence_blocks": decoded_blocks.len(),
    });
    black_box(
        serde_json::to_vec(&record)
            .expect("measurement record serializes")
            .len(),
    )
}

fn bench_memory_matrix(c: &mut Criterion) {
    let mut group = c.benchmark_group("trace_memory_smoke");
    group.sample_size(10);
    let case = MemoryCase {
        threads: 1,
        profile: SensitivityProfile::Balanced,
        shape: QueryShape::Small,
        candidates: 1,
        cache_block_bytes: 64 * 1024,
    };
    for component in COMPONENTS {
        let name = format!(
            "{}__{}__{}__t{}__c{}__b{}",
            component_name(component),
            profile_name(case.profile),
            shape_name(case.shape),
            case.threads,
            case.candidates,
            case.cache_block_bytes
        );
        group.bench_function(name, |bench| {
            bench.iter(|| black_box(for_case(case, component)))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_memory_matrix);
criterion_main!(benches);
