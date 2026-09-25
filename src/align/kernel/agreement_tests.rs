//! Tests that the kernel paths agree: AVX2 and scalar, resident and checkpointed traceback,
//! i16 and i32 scores.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::rebased::{
    RebaseLimits, SENTINEL_BASE, TEST_HOOKS, TestHooks, rebased_local_scores, update_hooks,
};
use super::tests::mutate;
use super::*;
use crate::test_support::{random_bases, rng};

/// Whether the host runs the AVX2 kernels; tests that compare them say so when it does not.
fn avx2() -> bool {
    let vector = Backend::Detected.vector();
    if !vector {
        eprintln!("skipped: this host has no AVX2");
    }
    vector
}

fn align(
    mut workspace: AlignmentWorkspace,
    query: &[u8],
    contig: &[u8],
    strand: Strand,
    config: AlignmentConfig,
) -> Result<Alignment, AlignmentError> {
    let input = AlignmentInput {
        query,
        query_offset: 3,
        contig,
        contig_offset: 17,
    };
    workspace.align(input, strand, config)
}

/// The scalar result, checked against the AVX2 path with rebased and with i32 scores.
fn assert_paths_agree(
    query: &[u8],
    contig: &[u8],
    strand: Strand,
    config: AlignmentConfig,
) -> Result<Alignment, AlignmentError> {
    let scalar = align(AlignmentWorkspace::scalar(), query, contig, strand, config);
    let case = format!("{config:?} {strand:?} {} {}", query.len(), contig.len());
    let vector = align(AlignmentWorkspace::default(), query, contig, strand, config);
    assert_eq!(vector, scalar, "{case}");
    update_hooks(|hooks| hooks.force_i32 = true);
    let wide = align(AlignmentWorkspace::default(), query, contig, strand, config);
    update_hooks(|hooks| hooks.force_i32 = false);
    assert_eq!(wide, scalar, "{case}");
    scalar
}

#[test]
fn avx2_matches_scalar() {
    if !avx2() {
        return;
    }
    let mut rng = rng(0x9e37_79b9_7f4a_7c15);
    let mut aligned = 0;
    for case in 0..160u64 {
        let query = random_bases(&mut rng, 24 + (case as usize * 7) % 300);
        let mut contig = mutate(&mut rng, &query, 20 + case % 9 * 40, case % 2 == 0);
        if case % 5 == 0 {
            contig.splice(
                contig.len() / 2..contig.len() / 2,
                random_bases(&mut rng, 30),
            );
        }
        if case % 7 == 0 {
            contig = contig.to_ascii_lowercase();
        }
        let config = AlignmentConfig {
            band_width: [1, 5, 16, 40, 128][case as usize % 5],
            diagonal_offset: case as i64 % 11 - 5,
            ..AlignmentConfig::default()
        };
        for strand in [Strand::Forward, Strand::Reverse] {
            aligned += u32::from(assert_paths_agree(&query, &contig, strand, config).is_ok());
        }
    }
    assert!(aligned > 150, "{aligned}");
}

#[test]
fn unknown_bases_unmatched() {
    // N and IUPAC codes on both sides, in either case, are substitution columns on every path.
    let mut rng = rng(0x5eed_0f00_0000_0015);
    for case in 0..40u64 {
        let mut query = random_bases(&mut rng, 80 + case as usize * 9);
        let codes = b"NnRy-";
        for at in (case as usize % 7..query.len()).step_by(11) {
            query[at] = codes[at % codes.len()];
        }
        let contig = mutate(&mut rng, &query, 20 + case % 4 * 30, case % 2 == 0);
        let config = AlignmentConfig {
            band_width: [3, 16, 128][case as usize % 3],
            ..AlignmentConfig::default()
        };
        for strand in [Strand::Forward, Strand::Reverse] {
            let contig = match strand {
                Strand::Forward => contig.clone(),
                Strand::Reverse => crate::seeds::reverse_complement(&contig),
            };
            if let Ok(alignment) = assert_paths_agree(&query, &contig, strand, config) {
                let query_bases = &query[alignment.query_interval.start as usize - 3..]
                    [..alignment.query_interval.len() as usize];
                let known = query_bases
                    .iter()
                    .filter(|base| b"ACGTacgt".contains(base))
                    .count() as u64;
                assert!(alignment.counts().matches <= known, "case {case}");
            }
        }
    }
    let unknown = vec![b'N'; 200];
    let result = assert_paths_agree(
        &unknown,
        &unknown,
        Strand::Forward,
        AlignmentConfig::default(),
    );
    assert_eq!(result, Err(AlignmentError::NoAlignment));
}

#[test]
fn checkpoints_match_resident() {
    let mut rng = rng(0x1234_5678_9abc_def1);
    for case in 0..48u64 {
        let query = random_bases(&mut rng, 60 + (case as usize * 13) % 400);
        let contig = mutate(&mut rng, &query, 80, true);
        let base = AlignmentConfig {
            band_width: 2 + case as u32 % 15,
            diagonal_offset: case as i64 % 7 - 3,
            max_cells: usize::MAX,
            ..AlignmentConfig::default()
        };
        let expected = align(
            AlignmentWorkspace::scalar(),
            &query,
            &contig,
            Strand::Forward,
            base,
        );
        let smallest = base.band_width as usize + 1;
        for max_cells in [smallest, smallest + 7, 500, 2_000] {
            let config = AlignmentConfig { max_cells, ..base };
            let chunked = assert_paths_agree(&query, &contig, Strand::Forward, config);
            assert_eq!(chunked, expected, "case {case} max_cells {max_cells}");
        }
    }
}

#[test]
fn scores_filled_first() {
    // Above `max_cells` the vector path fills all waves without traceback first.
    let mut rng = rng(0x0bad_5eed_1234_5678);
    let query = random_bases(&mut rng, 6_000);
    let contig = mutate(&mut rng, &query, 30, true);
    let config = AlignmentConfig {
        band_width: 64,
        max_cells: 100_000,
        ..AlignmentConfig::default()
    };
    let alignment = assert_paths_agree(&query, &contig, Strand::Forward, config).unwrap();
    assert!(alignment.query_interval.len() > 5_500);
}

#[test]
fn base_shifts_exact() {
    if !avx2() {
        return;
    }
    let mut rng = rng(0x7777_1111_9999_3333);
    // Thresholds this low move the base every few waves, so short inputs cover shifts up,
    // shifts down, zero states stored below the base, and chunk restarts from a saved base.
    TEST_HOOKS.set(TestHooks {
        limits: Some(RebaseLimits {
            high: 150,
            settle: 60,
            low: 30,
        }),
        ..TestHooks::NONE
    });
    for round in 0..14u64 {
        let length = 200 + 37 * round as usize;
        let query = random_bases(&mut rng, length);
        let mut contig = mutate(&mut rng, &query, 10 + round % 5 * 20, round % 2 == 0);
        if round % 3 == 0 {
            contig.splice(length / 2..length / 2, random_bases(&mut rng, 40));
        }
        for band_width in [5, 16, 40] {
            for max_cells in [usize::MAX, 900] {
                let config = AlignmentConfig {
                    band_width,
                    max_cells,
                    diagonal_offset: [0, 7, -9][round as usize % 3],
                    ..AlignmentConfig::default()
                };
                let _ = assert_paths_agree(&query, &contig, Strand::Forward, config);
            }
        }
    }
    let hooks = TEST_HOOKS.replace(TestHooks::NONE);
    let (up, down) = (hooks.shifts_up, hooks.shifts_down);
    assert!(up > 100 && down > 20, "shifts up {up} and down {down}");
}

#[test]
fn reaches_sentinel_bases() {
    if !avx2() {
        return;
    }
    // Scores above 40,000 need a base beyond the i16 range of a zero state.
    let mut rng = rng(0x5555_aaaa_1234_4321);
    let query = random_bases(&mut rng, 30_000);
    let contig = mutate(&mut rng, &query, 12, false);
    let config = AlignmentConfig {
        band_width: 16,
        ..AlignmentConfig::default()
    };
    TEST_HOOKS.set(TestHooks::NONE);
    let alignment = assert_paths_agree(&query, &contig, Strand::Forward, config).unwrap();
    let hooks = TEST_HOOKS.get();
    let (highest, up) = (hooks.highest_base, hooks.shifts_up);
    assert!(highest > SENTINEL_BASE && up > 0, "{highest} {up}");
    assert_eq!(hooks.reruns, 0);
    assert!(alignment.score > 40_000);
}

#[test]
fn violation_reruns_i32() {
    if !avx2() {
        return;
    }
    let mut rng = rng(0x3141_5926_5358_9793);
    let config = AlignmentConfig {
        band_width: 16,
        ..AlignmentConfig::default()
    };
    let query = random_bases(&mut rng, 4_000);
    let contig = mutate(&mut rng, &query, 15, true);
    let expected = align(
        AlignmentWorkspace::scalar(),
        &query,
        &contig,
        Strand::Forward,
        config,
    );
    TEST_HOOKS.set(TestHooks {
        violation_wave: Some(500),
        ..TestHooks::NONE
    });
    let rebased = align(
        AlignmentWorkspace::default(),
        &query,
        &contig,
        Strand::Forward,
        config,
    );
    let hooks = TEST_HOOKS.replace(TestHooks::NONE);
    assert_eq!(rebased, expected);
    assert_eq!(hooks.reruns, 1);
}

#[test]
fn rebased_spread_bound() {
    let default = AlignmentConfig::default();
    assert!(rebased_local_scores(default));
    assert!(rebased_local_scores(AlignmentConfig {
        band_width: 1_088,
        ..default
    }));
    // A very wide band or steep gap costs can spread the scores of one wave too far.
    assert!(!rebased_local_scores(AlignmentConfig {
        band_width: 4_000,
        ..default
    }));
    assert!(!rebased_local_scores(AlignmentConfig {
        gap_extend_score: -20_000,
        ..default
    }));
    assert!(!rebased_local_scores(AlignmentConfig {
        match_score: 40_000,
        ..default
    }));
}
