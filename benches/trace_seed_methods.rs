//! Bounded positional-seed ablations used by `archive-seeds` manifests.
//!
//! A digest is measured as a candidate-generation operation only.  The
//! production lookup path still verifies packed selected bases or sequence
//! bytes before creating an anchor.  This target covers the conventional
//! FracMinHash levels, all three Gear variants, syncmer/spaced/strobemer
//! rescue, and the selected gap-directed cascade.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use jam_rs::trace::anchors::Anchor;
use jam_rs::trace::chain::{
    AnchorClass, ChainConfig, WeightedAnchor, chain_anchors, chain_weighted_anchors,
};
use jam_rs::trace::config::SeedSensitivity;
use jam_rs::trace::model::{CoordinateModel, Strand};
use jam_rs::trace::seeds::extract_seed_level;
use jam_rs::trace::seeds::gear::{
    FragmentationMode, GearAnchorConfig, GearConfig, GearTableKind, MiniSketchConfig,
    build_fragment_mini_sketches, fragment_sequence, gear_selected_anchors,
};
use jam_rs::trace::seeds::spaced::{SpacedSeedConfig, extract_spaced_seeds};
use jam_rs::trace::seeds::strobemer::{StrobemerConfig, extract_strobemers};
use jam_rs::trace::seeds::syncmer::{SyncmerConfig, SyncmerMode, extract_syncmers};
use std::hint::black_box;

fn synthetic_dna(length: usize) -> Vec<u8> {
    let mut state = 0x4b1d_93c5_afe0_7811u64;
    (0..length)
        .map(|index| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            if index % 7_777 == 0 {
                b'N'
            } else {
                b"ACGT"[((state >> 62) & 3) as usize]
            }
        })
        .collect()
}

fn seed_count(sequence: &[u8], k: u8, scale: u64) -> usize {
    extract_seed_level(
        sequence,
        SeedSensitivity {
            k,
            scale,
            max_occurrences: 128,
        },
    )
    .expect("fixed FracMinHash configuration")
    .seeds
    .len()
}

fn bench_fracminhash(c: &mut Criterion) {
    let sequence = synthetic_dna(120_000);
    let mut group = c.benchmark_group("trace_seed_fracminhash");
    for (name, k, scale) in [
        ("fracminhash-k31", 31, 200),
        ("fracminhash-k21", 21, 500),
        ("denser-k21-gap", 21, 100),
    ] {
        group.bench_with_input(
            BenchmarkId::from_parameter(name),
            &(k, scale),
            |bench, &(k, scale)| {
                bench.iter(|| black_box(seed_count(black_box(&sequence), k, scale)))
            },
        );
    }
    group.finish();
}

fn bench_gear_variants(c: &mut Criterion) {
    let sequence = synthetic_dna(120_000);
    let config = GearConfig::fine(GearTableKind::SingleBase, 0x5eed);
    let mut group = c.benchmark_group("trace_seed_gear");
    group.bench_function("gear-fragment-exact", |bench| {
        bench.iter(|| {
            black_box(
                fragment_sequence(
                    black_box(&sequence),
                    0,
                    config,
                    FragmentationMode::StrandSymmetric,
                )
                .expect("Gear fragmentation")
                .len(),
            )
        })
    });
    for (name, anchor_config) in [
        ("gear-anchor-k31", GearAnchorConfig::k31()),
        ("gear-anchor-k21-rescue", GearAnchorConfig::k21_rescue()),
    ] {
        group.bench_function(name, |bench| {
            bench.iter(|| {
                black_box(
                    gear_selected_anchors(
                        black_box(&sequence),
                        config,
                        FragmentationMode::StrandSymmetric,
                        anchor_config,
                    )
                    .expect("Gear anchor selection")
                    .len(),
                )
            })
        });
    }
    group.bench_function("gear-fragment-mini-sketch", |bench| {
        bench.iter(|| {
            black_box(
                build_fragment_mini_sketches(
                    black_box(&sequence),
                    0,
                    config,
                    FragmentationMode::StrandSymmetric,
                    MiniSketchConfig::k21(8),
                )
                .expect("Gear mini-sketch")
                .len(),
            )
        })
    });
    group.finish();
}

fn bench_alternative_seeds(c: &mut Criterion) {
    let sequence = synthetic_dna(60_000);
    let mut group = c.benchmark_group("trace_seed_alternatives");
    for (name, config) in [
        ("syncmer-k31", SyncmerConfig::k31(SyncmerMode::Closed)),
        (
            "syncmer-k21-rescue",
            SyncmerConfig::k21_rescue(SyncmerMode::Closed),
        ),
    ] {
        group.bench_function(name, |bench| {
            bench.iter(|| {
                black_box(
                    extract_syncmers(black_box(&sequence), config)
                        .expect("syncmer configuration")
                        .len(),
                )
            })
        });
    }
    group.bench_function("spaced-31-weight21", |bench| {
        bench.iter(|| {
            black_box(
                extract_spaced_seeds(black_box(&sequence), SpacedSeedConfig::span31_weight21())
                    .expect("spaced configuration")
                    .len(),
            )
        })
    });
    for (name, config) in [
        ("strobemer-21x2", StrobemerConfig::paired_k21()),
        ("strobemer-21-minhash", StrobemerConfig::strobemer_k21()),
    ] {
        group.bench_function(name, |bench| {
            bench.iter(|| {
                black_box(
                    extract_strobemers(black_box(&sequence), config)
                        .expect("strobemer configuration")
                        .len(),
                )
            })
        });
    }
    group.finish();
}

fn bench_gap_directed_cascade(c: &mut Criterion) {
    let sequence = synthetic_dna(120_000);
    let config = GearConfig::fine(GearTableKind::SingleBase, 0x5eed);
    let mut group = c.benchmark_group("trace_seed_cascade");
    group.bench_function("selected-mixed-fracminhash", |bench| {
        bench.iter(|| {
            let primary = seed_count(&sequence, 31, 200);
            let rescue = seed_count(&sequence, 21, 100);
            black_box(primary + rescue)
        })
    });
    group.bench_function("experimental-gear-spaced-strobemer", |bench| {
        bench.iter(|| {
            let exact = fragment_sequence(
                black_box(&sequence),
                0,
                config,
                FragmentationMode::StrandSymmetric,
            )
            .expect("Gear fragments");
            let primary = seed_count(&sequence, 31, 200);
            let rescue = seed_count(&sequence, 21, 100);
            let spaced = extract_spaced_seeds(&sequence, SpacedSeedConfig::span31_weight21())
                .expect("spaced rescue")
                .len();
            let strobemers = extract_strobemers(&sequence, StrobemerConfig::paired_k21())
                .expect("paired rescue")
                .len();
            black_box(exact.len() + primary + rescue + spaced + strobemers)
        })
    });
    group.finish();
}

fn bench_mixed_chains(c: &mut Criterion) {
    let k31 = (0..2_000u64)
        .map(|index| Anchor {
            query_position: index * 50,
            target_position: index * 50 + 1_000,
            contig_id: 0,
            strand: Strand::Forward,
            k: 31,
            hash: index + 1,
            canonical_kmer: index + 11,
            query_reverse: false,
            target_reverse: false,
        })
        .collect::<Vec<_>>();
    let k21 = (0..2_000u64)
        .map(|index| Anchor {
            query_position: index * 50 + 25,
            target_position: index * 50 + 1_025,
            contig_id: 0,
            strand: Strand::Forward,
            k: 21,
            hash: index + 10_001,
            canonical_kmer: index + 20_001,
            query_reverse: false,
            target_reverse: false,
        })
        .collect::<Vec<_>>();
    let mixed = k31
        .iter()
        .copied()
        .map(|anchor| WeightedAnchor::new(anchor, AnchorClass::SpecificK31))
        .chain(
            k21.iter()
                .copied()
                .map(|anchor| WeightedAnchor::new(anchor, AnchorClass::SpecificK21)),
        )
        .collect::<Vec<_>>();
    let config = ChainConfig {
        max_chains: 16,
        min_anchors: 2,
        max_predecessors: 256,
        max_query_gap: 200_000,
        max_target_gap: 200_000,
        gap_penalty: 1,
        coordinate_model: CoordinateModel::Linear,
    };
    let mut group = c.benchmark_group("trace_seed_chaining");
    group.bench_function("separate-k31-k21", |bench| {
        bench.iter(|| {
            let primary = chain_anchors(black_box(&k31), 100_000, config).unwrap();
            let rescue = chain_anchors(black_box(&k21), 100_000, config).unwrap();
            black_box((primary, rescue))
        })
    });
    group.bench_function("mixed-k31-k21", |bench| {
        bench.iter(|| black_box(chain_weighted_anchors(black_box(&mixed), 100_000, config)))
    });
    group.finish();
}

criterion_group!(
    seed_benches,
    bench_fracminhash,
    bench_gear_variants,
    bench_alternative_seeds,
    bench_gap_directed_cascade,
    bench_mixed_chains
);
criterion_main!(seed_benches);
