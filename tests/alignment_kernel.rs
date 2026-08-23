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
