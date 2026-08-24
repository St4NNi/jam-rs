//! Bounded Gear-fragment and mini-sketch benchmark helpers.
//!
//! The fine/coarse stream sizes and independent table variants are recovered
//! from the historical Gear proposal.  Positional anchors and mini-sketches
//! are new search variants; this bench intentionally measures candidate
//! generation only and does not treat a digest as alignment evidence.

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use jam_rs::trace::seeds::gear::{
    FragmentationMode, GearAnchorConfig, GearConfig, GearTableKind, MiniSketchConfig,
    build_fragment_mini_sketches, fragment_sequence, gear_selected_anchors,
};

fn synthetic_dna(length: usize) -> Vec<u8> {
    let mut state = 0x1234_5678_9abc_def0u64;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(2_862_933_555_777_941_757)
                .wrapping_add(3_039)
                .rotate_left(17);
            b"ACGT"[((state >> 32) & 3) as usize]
        })
        .collect()
}

fn bench_fragments(c: &mut Criterion) {
    let sequence = synthetic_dna(250_000);
    let mut group = c.benchmark_group("gear_fragments");
    for kind in [
        GearTableKind::SingleBase,
        GearTableKind::Dinucleotide,
        GearTableKind::PackedFourBase,
    ] {
        let id = kind.identifier();
        group.bench_with_input(BenchmarkId::new("fragment", id), &kind, |bench, kind| {
            let config = GearConfig::fine(*kind, 0x5eed);
            bench.iter(|| {
                fragment_sequence(&sequence, 0, config, FragmentationMode::Forward).unwrap()
            });
        });
    }
    group.finish();
}

fn bench_anchors_and_minisketch(c: &mut Criterion) {
    let sequence = synthetic_dna(250_000);
    let config = GearConfig::fine(GearTableKind::SingleBase, 0x5eed);
    let mut group = c.benchmark_group("gear_positional_variants");
    group.bench_function("selected-k31", |bench| {
        bench.iter(|| {
            gear_selected_anchors(
                &sequence,
                config,
                FragmentationMode::Forward,
                GearAnchorConfig::k31(),
            )
            .unwrap()
        });
    });
    group.bench_function("selected-k21-rescue", |bench| {
        bench.iter(|| {
            gear_selected_anchors(
                &sequence,
                config,
                FragmentationMode::Forward,
                GearAnchorConfig::k21_rescue(),
            )
            .unwrap()
        });
    });
    group.bench_function("fragment-mini-sketch", |bench| {
        bench.iter(|| {
            build_fragment_mini_sketches(
                &sequence,
                0,
                config,
                FragmentationMode::Forward,
                MiniSketchConfig::k21(8),
            )
            .unwrap()
        });
    });
    group.finish();
}

criterion_group!(gear_benches, bench_fragments, bench_anchors_and_minisketch);
criterion_main!(gear_benches);
