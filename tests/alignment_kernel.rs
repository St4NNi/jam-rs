use jam_rs::trace::alignment::{
    AlignmentError, AlignmentOptions, AlignmentWorkspace, align, align_circular, align_oriented,
};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::model::{BaseInterval, EditOperation, Strand};

fn scoring() -> AlignmentScoring {
    AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 32,
    }
}

#[test]
fn exact_match_reports_cigar_and_reconstructs() {
    let query = b"ACGTACGT";
    let target = b"TTACGTACGTAA";
    let result = align(query, target, scoring()).unwrap();

    assert_eq!(result.cigar, "8=");
    assert_eq!(result.matches, 8);
    assert_eq!(result.substitutions, 0);
    assert_eq!(result.query_interval, BaseInterval::new(0, 8).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(2, 10).unwrap());
    assert_eq!(
        result.query_segments,
        vec![BaseInterval::new(0, 8).unwrap()]
    );

    let (aligned_query, aligned_target) = result.reconstruct(query, target).unwrap();
    assert_eq!(aligned_query, query);
    assert_eq!(aligned_target, b"ACGTACGT");
}

#[test]
fn reverse_complement_uses_forward_target_coordinates() {
    let target_forward = b"GGGAAGTCCTTAAA";
    let query = b"GACTT";
    let result = align_oriented(query, target_forward, 100, Strand::Reverse, scoring()).unwrap();

    assert_eq!(result.strand, Strand::Reverse);
    assert_eq!(result.target_interval, BaseInterval::new(103, 108).unwrap());
    assert_eq!(result.cigar, "5=");
    assert!(!result.origin_crossing);
}

#[test]
fn partial_trace_and_edit_script_keep_direct_evidence() {
    let query = b"GGGACGTCCC";
    let target = b"TTTACGACCCAAA";
    let result = align(query, target, scoring()).unwrap();

    assert_eq!(result.query_interval, BaseInterval::new(3, 10).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(3, 10).unwrap());
    assert_eq!(result.cigar, "3=1X3=");
    assert_eq!(result.matches, 6);
    assert_eq!(result.substitutions, 1);
    assert_eq!(result.insertions, 0);
    assert_eq!(result.deletions, 0);
    assert_eq!(
        result.edit_script,
        vec![
            jam_rs::trace::model::EditRun {
                operation: EditOperation::Equal,
                length: 3,
            },
            jam_rs::trace::model::EditRun {
                operation: EditOperation::Substitution,
                length: 1,
            },
            jam_rs::trace::model::EditRun {
                operation: EditOperation::Equal,
                length: 3,
            },
        ]
    );
}

#[test]
fn affine_insertion_and_deletion_reconstruct() {
    let insertion = align(b"ACGTACGT", b"ACGTGACGT", scoring()).unwrap();
    assert_eq!(insertion.cigar, "4=1I4=");
    assert_eq!(insertion.insertions, 1);
    assert_eq!(
        insertion.reconstruct(b"ACGTACGT", b"ACGTGACGT").unwrap(),
        (b"ACGT-ACGT".to_vec(), b"ACGTGACGT".to_vec())
    );

    let mut deletion_scoring = scoring();
    deletion_scoring.gap_open_score = -4;
    let deletion = align(b"ACGTACGT", b"ACGACGT", deletion_scoring).unwrap();
    assert_eq!(deletion.cigar, "3=1D4=");
    assert_eq!(deletion.deletions, 1);
    assert_eq!(
        deletion.reconstruct(b"ACGTACGT", b"ACGACGT").unwrap(),
        (b"ACGTACGT".to_vec(), b"ACG-ACGT".to_vec())
    );
}

#[test]
fn long_gaps_and_mixed_edits_reconstruct_and_score() {
    let mut gap_scoring = scoring();
    gap_scoring.gap_open_score = -1;
    gap_scoring.gap_extend_score = -1;

    let insertion_target = b"ACGTNNNNNACGT";
    let insertion = align(b"ACGTACGT", insertion_target, gap_scoring).unwrap();
    assert_eq!(insertion.insertions, 5);
    assert_alignment_invariants(b"ACGTACGT", insertion_target, &insertion, gap_scoring);

    let deletion_query = b"ACGTNNNNNNACGT";
    let deletion = align(deletion_query, b"ACGTACGT", gap_scoring).unwrap();
    assert_eq!(deletion.deletions, 6);
    assert_alignment_invariants(deletion_query, b"ACGTACGT", &deletion, gap_scoring);

    let mixed_query = b"ACGTACGTGACGTACGT";
    let mixed_target = b"ACGTACGTACGTTCGT";
    let mixed = align(mixed_query, mixed_target, gap_scoring).unwrap();
    assert!(mixed.substitutions > 0 || mixed.insertions > 0 || mixed.deletions > 0);
    assert_alignment_invariants(mixed_query, mixed_target, &mixed, gap_scoring);
}

#[test]
fn full_affine_oracle_matches_wide_band_and_result_invariants() {
    let cases = [
        (b"ACGTACGT".as_slice(), b"ACGTACGT".as_slice()),
        (b"ACGTACGT".as_slice(), b"TTACGTACGTAA".as_slice()),
        (b"ACGTTACGT".as_slice(), b"ACGTACGT".as_slice()),
        (b"GGACGTTA".as_slice(), b"TTACCTAA".as_slice()),
        (b"AACCGGTTAACCGG".as_slice(), b"AACCTTAACCGG".as_slice()),
        (b"GATTACAGATTACA".as_slice(), b"GCATGCUGATTACA".as_slice()),
        (
            b"AAAAACCCCCGGGGG".as_slice(),
            b"TTAAAAACCCCCAAAGGGGGCC".as_slice(),
        ),
        (b"acNRY".as_slice(), b"TTACNRYAA".as_slice()),
    ];
    let mut wide = scoring();
    wide.band_width = 64;
    for (query, target) in cases {
        let result = align(query, target, wide).unwrap();
        assert_eq!(
            result.score,
            oracle_local_affine_score(query, target, wide),
            "full-band score mismatch for query={query:?}, target={target:?}"
        );
        assert_alignment_invariants(query, target, &result, wide);
    }
}

#[test]
fn alignment_covers_query_and_target_clipping() {
    let query_clipped = align(b"TTACGTAA", b"ACGT", scoring()).unwrap();
    assert_eq!(
        query_clipped.query_interval,
        BaseInterval::new(2, 6).unwrap()
    );
    assert_eq!(
        query_clipped.target_interval,
        BaseInterval::new(0, 4).unwrap()
    );
    assert_alignment_invariants(b"TTACGTAA", b"ACGT", &query_clipped, scoring());

    let target_clipped = align(b"ACGT", b"TTACGTAA", scoring()).unwrap();
    assert_eq!(
        target_clipped.query_interval,
        BaseInterval::new(0, 4).unwrap()
    );
    assert_eq!(
        target_clipped.target_interval,
        BaseInterval::new(2, 6).unwrap()
    );
    assert_alignment_invariants(b"ACGT", b"TTACGTAA", &target_clipped, scoring());
}

#[test]
fn local_alignment_selects_the_best_of_several_regions() {
    let query = b"AAAAACCCCCGGGGG";
    let target = b"TTAAAAACCCCCAAAAAAAAAAGGGGGCC";
    let result = align(query, target, scoring()).unwrap();

    assert_eq!(result.score, 20);
    assert_eq!(result.matches, 10);
    assert_eq!(result.query_interval, BaseInterval::new(0, 10).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(2, 12).unwrap());
    assert_alignment_invariants(query, target, &result, scoring());
}

#[test]
fn band_edge_is_accepted_and_excluded_band_is_explicit() {
    let query = b"ACGT";
    let target = b"TTACGTAA";
    let mut edge_scoring = scoring();
    edge_scoring.band_width = 0;
    let edge_options = AlignmentOptions::new(edge_scoring).with_diagonal_offset(2);
    let edge = align(query, target, edge_options).unwrap();
    assert_eq!(edge.cigar, "4=");
    assert_alignment_invariants(query, target, &edge, edge_scoring);

    let excluded = AlignmentOptions::new(edge_scoring).with_diagonal_offset(100);
    assert!(matches!(
        align(query, target, excluded),
        Err(AlignmentError::BandExcludesInput)
    ));

    let too_narrow = AlignmentOptions::new(edge_scoring).with_diagonal_offset(0);
    assert!(matches!(
        align(query, target, too_narrow),
        Err(AlignmentError::NoAlignment)
    ));
}

#[test]
fn ambiguous_bases_are_case_insensitive_but_not_wildcards() {
    let result = align(b"acNRY", b"TTACNRYAA", scoring()).unwrap();
    assert_eq!(result.cigar, "5=");
    assert_eq!(result.matches, 5);
    assert_alignment_invariants(b"acNRY", b"TTACNRYAA", &result, scoring());

    let mismatch = align(b"AAAAAAA", b"AAACAAA", scoring()).unwrap();
    assert_eq!(mismatch.substitutions, 1);
    assert_eq!(mismatch.matches, 6);
    assert_alignment_invariants(b"AAAAAAA", b"AAACAAA", &mismatch, scoring());
}

#[test]
fn circular_query_crossing_origin_has_two_ordered_segments() {
    let result = align_circular(b"ACGTAC", 4, 4, b"ACAC", 50, Strand::Forward, scoring()).unwrap();

    assert!(result.origin_crossing);
    assert_eq!(
        result.query_segments,
        vec![
            BaseInterval::new(4, 6).unwrap(),
            BaseInterval::new(0, 2).unwrap()
        ]
    );
    assert_eq!(result.target_interval, BaseInterval::new(50, 54).unwrap());
}

#[test]
fn circular_alignment_covering_full_query_keeps_origin_order() {
    let query = b"ACGTAC";
    let result = align_circular(query, 4, 6, b"ACACGT", 0, Strand::Forward, scoring()).unwrap();

    assert_eq!(result.cigar, "6=");
    assert!(result.origin_crossing);
    assert_eq!(
        result.query_segments,
        vec![
            BaseInterval::new(4, 6).unwrap(),
            BaseInterval::new(0, 4).unwrap()
        ]
    );
    assert_eq!(
        result
            .query_segments
            .iter()
            .map(|segment| segment.len())
            .sum::<u64>(),
        6
    );
    assert_alignment_invariants(b"ACACGT", b"ACACGT", &result, scoring());
}

#[test]
fn reverse_complement_alignment_maps_to_same_forward_interval() {
    let query = b"ACGTTAG";
    let target_forward = b"TTGGGCTAACGTAAA";
    let target_reverse = reverse_complement(target_forward);
    let forward = align(query, target_reverse.as_slice(), scoring()).unwrap();
    let reverse = align_oriented(query, target_forward, 100, Strand::Reverse, scoring()).unwrap();

    let mapped_start = target_forward.len() as u64 - forward.target_interval.end + 100;
    let mapped_end = target_forward.len() as u64 - forward.target_interval.start + 100;
    assert_eq!(
        reverse.target_interval,
        BaseInterval::new(mapped_start, mapped_end).unwrap()
    );
    assert_eq!(forward.cigar, reverse.cigar);
    assert_eq!(forward.score, reverse.score);
    assert_alignment_invariants(query, target_reverse.as_slice(), &forward, scoring());
}

#[test]
fn workspace_reuses_allocations_and_limits_matrix() {
    let mut workspace = AlignmentWorkspace::new();
    let options = AlignmentOptions::new(scoring());
    let _ = workspace.align(b"ACGTACGT", b"ACGTACGT", options).unwrap();
    let capacity = workspace.capacity_cells();
    let _ = workspace.align(b"ACGT", b"ACGT", options).unwrap();
    assert_eq!(workspace.capacity_cells(), capacity);

    let limited = AlignmentOptions::new(scoring()).with_max_cells(1);
    assert!(matches!(
        workspace.align(b"ACGT", b"ACGT", limited),
        Err(AlignmentError::MatrixTooLarge { .. })
    ));
}

#[test]
fn empty_and_nonmatching_inputs_have_explicit_errors() {
    assert!(matches!(
        align(b"", b"ACGT", scoring()),
        Err(AlignmentError::EmptyQuery)
    ));
    assert!(matches!(
        align(b"ACGT", b"", scoring()),
        Err(AlignmentError::EmptyTarget)
    ));
    assert!(matches!(
        align(b"AAAA", b"TTTT", scoring()),
        Err(AlignmentError::NoAlignment)
    ));
}

#[test]
fn local_affine_score_matches_test_only_oracle() {
    let cases = [
        (b"ACGT".as_slice(), b"ACGT".as_slice()),
        (b"ACGTT".as_slice(), b"ACGT".as_slice()),
        (b"GGACGTA".as_slice(), b"TTACCTA".as_slice()),
        (b"AACCGGTT".as_slice(), b"AATCGGTT".as_slice()),
        (b"GATTACA".as_slice(), b"GCATGCU".as_slice()),
    ];
    let mut oracle_scoring = scoring();
    oracle_scoring.gap_open_score = 0;
    oracle_scoring.gap_extend_score = -1;
    for (query, target) in cases {
        let mut options = AlignmentOptions::new(oracle_scoring);
        options.scoring.band_width = 32;
        let result = align(query, target, options).unwrap();
        assert_eq!(
            result.score,
            oracle_local_affine_score(query, target, oracle_scoring)
        );
    }
}

fn oracle_local_affine_score(query: &[u8], target: &[u8], scoring: AlignmentScoring) -> i32 {
    let mut match_state = vec![vec![0_i32; target.len() + 1]; query.len() + 1];
    let mut insertion = match_state.clone();
    let mut deletion = match_state.clone();
    let mut best = 0_i32;
    for i in 1..=query.len() {
        for j in 1..=target.len() {
            let substitution = if query[i - 1].eq_ignore_ascii_case(&target[j - 1]) {
                scoring.match_score
            } else {
                scoring.mismatch_score
            };
            let previous = match_state[i - 1][j - 1]
                .max(insertion[i - 1][j - 1])
                .max(deletion[i - 1][j - 1]);
            match_state[i][j] = (previous + substitution).max(0);
            insertion[i][j] = (insertion[i][j - 1] + scoring.gap_extend_score)
                .max(match_state[i][j - 1] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(deletion[i][j - 1] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(0);
            deletion[i][j] = (deletion[i - 1][j] + scoring.gap_extend_score)
                .max(match_state[i - 1][j] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(insertion[i - 1][j] + scoring.gap_open_score + scoring.gap_extend_score)
                .max(0);
            best = best
                .max(match_state[i][j])
                .max(insertion[i][j])
                .max(deletion[i][j]);
        }
    }
    best
}

fn assert_alignment_invariants(
    query: &[u8],
    target: &[u8],
    result: &jam_rs::trace::alignment::AlignmentResult,
    scoring: AlignmentScoring,
) {
    let query_start = result.query_interval.start as usize;
    let query_end = result.query_interval.end as usize;
    let target_start = result.target_interval.start as usize;
    let target_end = result.target_interval.end as usize;
    assert!(query_start <= query_end && query_end <= query.len());
    assert!(target_start <= target_end && target_end <= target.len());

    let (aligned_query, aligned_target) = result.reconstruct(query, target).unwrap();
    let mut query_consumed = 0_u64;
    let mut target_consumed = 0_u64;
    let mut matches = 0_u64;
    let mut substitutions = 0_u64;
    let mut insertions = 0_u64;
    let mut deletions = 0_u64;
    let mut recomputed_score = 0_i32;
    let mut cigar = String::new();
    for run in &result.edit_script {
        assert!(run.length > 0);
        let code = match run.operation {
            EditOperation::Equal => {
                query_consumed += u64::from(run.length);
                target_consumed += u64::from(run.length);
                matches += u64::from(run.length);
                recomputed_score += scoring.match_score * run.length as i32;
                '='
            }
            EditOperation::Substitution => {
                query_consumed += u64::from(run.length);
                target_consumed += u64::from(run.length);
                substitutions += u64::from(run.length);
                recomputed_score += scoring.mismatch_score * run.length as i32;
                'X'
            }
            EditOperation::Insertion => {
                target_consumed += u64::from(run.length);
                insertions += u64::from(run.length);
                recomputed_score +=
                    scoring.gap_open_score + scoring.gap_extend_score * run.length as i32;
                'I'
            }
            EditOperation::Deletion => {
                query_consumed += u64::from(run.length);
                deletions += u64::from(run.length);
                recomputed_score +=
                    scoring.gap_open_score + scoring.gap_extend_score * run.length as i32;
                'D'
            }
            EditOperation::SoftClip => panic!("local alignments do not emit soft clips"),
        };
        use std::fmt::Write as _;
        write!(&mut cigar, "{}{}", run.length, code).unwrap();
    }

    assert_eq!(query_consumed, result.query_interval.len());
    assert_eq!(target_consumed, result.target_interval.len());
    assert_eq!(query_consumed, result.query_length);
    assert_eq!(target_consumed, result.target_length);
    assert_eq!(matches, result.matches);
    assert_eq!(substitutions, result.substitutions);
    assert_eq!(insertions, result.insertions);
    assert_eq!(deletions, result.deletions);
    assert_eq!(cigar, result.cigar);
    assert_eq!(aligned_query.len(), aligned_target.len());
    assert_eq!(
        aligned_query.iter().filter(|base| **base != b'-').count() as u64,
        query_consumed
    );
    assert_eq!(
        aligned_target.iter().filter(|base| **base != b'-').count() as u64,
        target_consumed
    );
    let mut column = 0usize;
    for run in &result.edit_script {
        for _ in 0..run.length {
            let query_base = aligned_query[column];
            let target_base = aligned_target[column];
            match run.operation {
                EditOperation::Equal => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Substitution => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(!query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Insertion => {
                    assert_eq!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                }
                EditOperation::Deletion => {
                    assert_ne!(query_base, b'-');
                    assert_eq!(target_base, b'-');
                }
                EditOperation::SoftClip => panic!("local alignments do not emit soft clips"),
            }
            column += 1;
        }
    }
    assert_eq!(column, aligned_query.len());
    assert_eq!(recomputed_score, result.score);
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => b'N',
        })
        .collect()
}
