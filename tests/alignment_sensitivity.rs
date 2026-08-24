use jam_rs::trace::alignment::{
    AlignmentCorridor, AlignmentError, AlignmentOptions, AlignmentWorkspace, CorridorAnchor,
};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::model::{BaseInterval, Strand};

fn scoring() -> AlignmentScoring {
    AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 32,
    }
}

fn corridor(anchors: impl Into<Vec<CorridorAnchor>>, safety: u32) -> AlignmentCorridor {
    AlignmentCorridor::new(anchors, safety).unwrap()
}

#[test]
fn exact_corridor_alignment_keeps_fast_local_result() {
    let query = b"ACGTACGT";
    let target = b"TTACGTACGTAA";
    let anchors = corridor(
        [CorridorAnchor::new(0, 2, 4), CorridorAnchor::new(4, 6, 4)],
        4,
    );
    let result = AlignmentWorkspace::new()
        .align_with_retries(query, target, &anchors, scoring())
        .unwrap();

    assert_eq!(result.cigar, "8=");
    assert_eq!(result.query_interval, BaseInterval::new(0, 8).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(2, 10).unwrap());
    assert_eq!(result.retry_metadata.attempted_widths, vec![64]);
    assert_eq!(result.retry_metadata.retries, 0);
    assert!(!result.retry_metadata.band_edge_touched);
}

#[test]
fn corridor_preserves_substitutions_and_oracle_score() {
    let query = b"ACGTACGT";
    let target = b"TTACGTTCGTAA";
    let anchors = corridor(
        [CorridorAnchor::new(0, 2, 4), CorridorAnchor::new(4, 6, 4)],
        8,
    );
    let mut workspace = AlignmentWorkspace::new();
    let result = workspace
        .align_with_retries(query, target, &anchors, scoring())
        .unwrap();
    let oracle = workspace.align(query, target, scoring()).unwrap();

    assert_eq!(result.score, oracle.score);
    assert_eq!(result.cigar, "4=1X3=");
    assert_eq!(result.substitutions, 1);
    assert_eq!(result.reconstruct(query, target).unwrap().0.len(), 8);
}

#[test]
fn piecewise_corridor_tracks_a_long_indel() {
    let mut query = b"AAAACCCC".to_vec();
    query.extend(std::iter::repeat_n(b'G', 8));
    query.extend_from_slice(b"TTTT");
    let mut target = b"AAAACCCC".to_vec();
    target.extend(std::iter::repeat_n(b'N', 100));
    target.extend(std::iter::repeat_n(b'G', 8));
    target.extend_from_slice(b"TTTT");
    let anchors = corridor(
        [
            CorridorAnchor::new(0, 0, 4),
            CorridorAnchor::new(12, 112, 4),
        ],
        8,
    );
    let mut indel_scoring = scoring();
    indel_scoring.gap_open_score = -1;
    indel_scoring.gap_extend_score = 0;
    let result = AlignmentWorkspace::new()
        .align_with_corridor(&query, &target, &anchors, indel_scoring)
        .unwrap();

    assert_eq!(result.query_interval, BaseInterval::new(0, 20).unwrap());
    assert_eq!(result.insertions, 100);
    assert!(result.cigar.contains("100I"));
    assert!(result.retry_metadata.predicted_drift);
}

#[test]
fn automatic_retries_widen_for_predicted_offset_drift() {
    let query = vec![b'A'; 32];
    let mut target = vec![b'T'; 200];
    target.extend_from_slice(&query);
    let anchors = corridor([CorridorAnchor::new(0, 0, 4)], 0);
    let result = AlignmentWorkspace::new()
        .align_with_retries(&query, &target, &anchors, scoring())
        .unwrap();

    assert_eq!(result.cigar, "32=");
    assert_eq!(result.target_interval, BaseInterval::new(200, 232).unwrap());
    assert_eq!(result.retry_metadata.attempted_widths, vec![64, 128, 256]);
    assert_eq!(result.retry_metadata.retries, 2);
}

#[test]
fn retry_cap_is_reported_when_the_final_edge_is_still_used() {
    let query = vec![b'A'; 32];
    let mut target = vec![b'T'; 2048];
    target.extend_from_slice(&query);
    let anchors = corridor([CorridorAnchor::new(0, 0, 4)], 0);
    let result = AlignmentWorkspace::new()
        .align_with_retries(&query, &target, &anchors, scoring())
        .unwrap();

    assert_eq!(result.retry_metadata.attempted_widths.last(), Some(&2048));
    assert!(result.retry_metadata.retry_cap_reached);
    assert!(result.retry_metadata.band_edge_touched);
}

#[test]
fn memory_bound_stops_retry_without_allocating_unbounded_state() {
    let query = b"ACGTACGT";
    let target = b"ACGTACGT";
    let anchors = corridor([CorridorAnchor::new(0, 0, 8)], 0);
    let options = AlignmentOptions::new(scoring()).with_max_cells(1);
    let error = AlignmentWorkspace::new().align_with_retries(query, target, &anchors, options);

    assert!(matches!(
        error,
        Err(AlignmentError::MatrixTooLarge { max_cells: 1, .. })
    ));
}

#[test]
fn anchored_semiglobal_keeps_true_fragment_ends() {
    let target = b"TTACGTAA";
    let anchors = corridor([CorridorAnchor::new(0, 2, 4)], 8);
    let result = AlignmentWorkspace::new()
        .align_anchored_semiglobal(b"ACGT", target, &anchors, scoring())
        .unwrap();

    assert!(result.retry_metadata.semiglobal_attempted);
    assert!(result.retry_metadata.semiglobal_accepted);
    assert_eq!(result.query_interval, BaseInterval::new(0, 4).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(2, 6).unwrap());
}

#[test]
fn anchored_semiglobal_allows_contig_end_clipping() {
    let target = b"ACGT";
    let anchors = corridor([CorridorAnchor::new(2, 0, 4)], 8);
    let result = AlignmentWorkspace::new()
        .align_anchored_semiglobal(b"TTACGTAA", target, &anchors, scoring())
        .unwrap();

    assert!(result.retry_metadata.semiglobal_attempted);
    assert!(!result.retry_metadata.semiglobal_accepted);
    assert_eq!(result.query_interval, BaseInterval::new(2, 6).unwrap());
    assert_eq!(result.target_interval, BaseInterval::new(0, 4).unwrap());
}

#[test]
fn reverse_complement_corridor_maps_back_to_forward_coordinates() {
    let query = b"GACTT";
    let target_forward = b"GGGAAGTCCTTAAA";
    // The query occurs at oriented [6, 11), corresponding to forward [3, 8).
    let anchors = corridor([CorridorAnchor::new(0, 6, 5)], 8);
    let result = AlignmentWorkspace::new()
        .align_with_retries_oriented(
            query,
            target_forward,
            100,
            Strand::Reverse,
            &anchors,
            scoring(),
        )
        .unwrap();

    assert_eq!(result.strand, Strand::Reverse);
    assert_eq!(result.target_interval, BaseInterval::new(103, 108).unwrap());
    assert_eq!(result.cigar, "5=");
}

#[test]
fn corridor_retry_result_matches_wide_band_oracle() {
    let query = b"ACGTGACCTAGT";
    let target = b"TTACGTGACCTAGTAA";
    let anchors = corridor(
        [CorridorAnchor::new(0, 2, 6), CorridorAnchor::new(6, 8, 6)],
        8,
    );
    let mut workspace = AlignmentWorkspace::new();
    let corridor_result = workspace
        .align_with_retries(query, target, &anchors, scoring())
        .unwrap();
    let mut wide = scoring();
    wide.band_width = 128;
    let oracle = workspace.align(query, target, wide).unwrap();

    assert_eq!(corridor_result.score, oracle.score);
    assert_eq!(corridor_result.cigar, oracle.cigar);
    assert_eq!(corridor_result.matches, oracle.matches);
}
