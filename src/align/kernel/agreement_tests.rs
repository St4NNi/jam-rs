//! The kernel paths agree: AVX2 and scalar, resident and checkpointed traceback, rebased i16
//! and i32 scores.

use super::driver::WIDE_LOCAL_KERNEL;
use super::rebased::{
    REBASE_FORCED_VIOLATION, REBASE_REACHED, REBASE_THRESHOLDS, SENTINEL_BASE, rebased_local_scores,
};
use super::tests::{test_dna, test_mutate};
use super::*;

fn avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
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
    WIDE_LOCAL_KERNEL.set(true);
    let wide = align(AlignmentWorkspace::default(), query, contig, strand, config);
    WIDE_LOCAL_KERNEL.set(false);
    assert_eq!(wide, scalar, "{case}");
    scalar
}

#[test]
fn avx2_and_scalar_agree_on_random_inputs() {
    if !avx2() {
        return;
    }
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut aligned = 0;
    for case in 0..160u64 {
        let query = test_dna(&mut state, 24 + (case as usize * 7) % 300);
        let mut contig = test_mutate(&mut state, &query, 20 + case % 9 * 40, case % 2 == 0);
        if case % 5 == 0 {
            contig.splice(contig.len() / 2..contig.len() / 2, test_dna(&mut state, 30));
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
fn checkpointed_traceback_matches_resident() {
    let mut state = 0x1234_5678_9abc_def1u64;
    for case in 0..48u64 {
        let query = test_dna(&mut state, 60 + (case as usize * 13) % 400);
        let contig = test_mutate(&mut state, &query, 80, true);
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
fn long_passes_fill_scores_first() {
    // Above `max_cells` the vector path fills all waves without traceback first.
    let mut state = 0x0bad_5eed_1234_5678u64;
    let query = test_dna(&mut state, 6_000);
    let contig = test_mutate(&mut state, &query, 30, true);
    let config = AlignmentConfig {
        band_width: 64,
        max_cells: 100_000,
        ..AlignmentConfig::default()
    };
    let alignment = assert_paths_agree(&query, &contig, Strand::Forward, config).unwrap();
    assert!(alignment.query_interval.len() > 5_500);
}

#[test]
fn rebased_base_shifts_match_i32() {
    if !avx2() {
        return;
    }
    let mut state = 0x7777_1111_9999_3333u64;
    // Thresholds this low move the base every few waves, so short inputs cover shifts up,
    // shifts down, zero states stored below the base, and chunk restarts from a saved base.
    REBASE_THRESHOLDS.set(Some((150, 60, 30)));
    REBASE_REACHED.set((0, 0, 0, 0));
    for round in 0..14u64 {
        let length = 200 + 37 * round as usize;
        let query = test_dna(&mut state, length);
        let mut contig = test_mutate(&mut state, &query, 10 + round % 5 * 20, round % 2 == 0);
        if round % 3 == 0 {
            contig.splice(length / 2..length / 2, test_dna(&mut state, 40));
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
    REBASE_THRESHOLDS.set(None);
    let (_, up, down, _) = REBASE_REACHED.get();
    assert!(up > 100 && down > 20, "shifts up {up} and down {down}");
}

#[test]
fn rebased_reaches_sentinel_bases() {
    if !avx2() {
        return;
    }
    // Scores above 40,000 need a base beyond the i16 range of a zero state.
    let mut state = 0x5555_aaaa_1234_4321u64;
    let query = test_dna(&mut state, 30_000);
    let contig = test_mutate(&mut state, &query, 12, false);
    let config = AlignmentConfig {
        band_width: 16,
        ..AlignmentConfig::default()
    };
    REBASE_REACHED.set((0, 0, 0, 0));
    let alignment = assert_paths_agree(&query, &contig, Strand::Forward, config).unwrap();
    let (highest, up, _, fallbacks) = REBASE_REACHED.get();
    assert!(highest > SENTINEL_BASE && up > 0, "{highest} {up}");
    assert_eq!(fallbacks, 0);
    assert!(alignment.score > 40_000);
}

#[test]
fn violation_reruns_on_i32() {
    if !avx2() {
        return;
    }
    let mut state = 0x3141_5926_5358_9793u64;
    let config = AlignmentConfig {
        band_width: 16,
        ..AlignmentConfig::default()
    };
    let query = test_dna(&mut state, 4_000);
    let contig = test_mutate(&mut state, &query, 15, true);
    let expected = align(
        AlignmentWorkspace::scalar(),
        &query,
        &contig,
        Strand::Forward,
        config,
    );
    REBASE_REACHED.set((0, 0, 0, 0));
    REBASE_FORCED_VIOLATION.set(Some(500));
    let rebased = align(
        AlignmentWorkspace::default(),
        &query,
        &contig,
        Strand::Forward,
        config,
    );
    REBASE_FORCED_VIOLATION.set(None);
    assert_eq!(rebased, expected);
    assert_eq!(REBASE_REACHED.get().3, 1);
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
