use jam_rs::trace::alignment::refine::{
    BridgeConfig, ChainEndpoint, ScoreOnlyWorkspace, ScoreProbe, bridge_chains,
    refine_with_score_only, score_only_corridor_with_scratch, select_retry_width,
};
use jam_rs::trace::alignment::{
    AlignmentCorridor, AlignmentError, AlignmentOptions, AlignmentResult, AlignmentRetryMetadata,
    AlignmentWorkspace, CorridorAnchor,
};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::model::{BaseInterval, Strand};

fn interval(start: u64, end: u64) -> BaseInterval {
    BaseInterval::new(start, end).unwrap()
}

fn endpoint(
    contig_id: u32,
    strand: Strand,
    query: (u64, u64),
    target: (u64, u64),
) -> ChainEndpoint {
    ChainEndpoint {
        contig_id,
        strand,
        query_interval: interval(query.0, query.1),
        target_interval: interval(target.0, target.1),
        support_bases: 31,
    }
}

fn dummy_result() -> AlignmentResult {
    AlignmentResult {
        score: 42,
        strand: Strand::Forward,
        query_interval: interval(0, 4),
        query_segments: vec![interval(0, 4)],
        target_interval: interval(0, 4),
        origin_crossing: false,
        query_length: 4,
        target_length: 4,
        matches: 4,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: "4=".to_string(),
        edit_script: vec![],
        chain_score: 0,
        retry_metadata: AlignmentRetryMetadata::default(),
    }
}

fn scoring() -> AlignmentScoring {
    AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 32,
    }
}

fn corridor(anchors: impl Into<Vec<CorridorAnchor>>, safety_margin: u32) -> AlignmentCorridor {
    AlignmentCorridor::new(anchors, safety_margin).unwrap()
}

#[test]
fn score_only_exact_corridor_matches_full_kernel() {
    let query = b"ACGTACGTACGT";
    let target = b"TTACGTACGTACGTAA";
    let exact_corridor = corridor(
        [CorridorAnchor::new(0, 2, 6), CorridorAnchor::new(6, 8, 6)],
        8,
    );
    let options = AlignmentOptions::new(scoring());
    let mut scratch = ScoreOnlyWorkspace::new();
    let probe =
        score_only_corridor_with_scratch(&mut scratch, query, target, &exact_corridor, options)
            .unwrap();
    let full = AlignmentWorkspace::new()
        .align_with_corridor(query, target, &exact_corridor, options)
        .unwrap();
    assert_eq!(probe.score, full.score);
    assert_eq!(
        probe.band_edge_touched,
        full.retry_metadata.band_edge_touched
    );
}

#[test]
fn score_only_substitution_and_indel_match_full_kernel() {
    let query = b"ACGTACGT";
    let substituted = b"TTACGTTCGTAA";
    let substitution_corridor = corridor(
        [CorridorAnchor::new(0, 2, 4), CorridorAnchor::new(4, 6, 4)],
        8,
    );
    let options = AlignmentOptions::new(scoring());
    let substitution = score_only_corridor_with_scratch(
        &mut ScoreOnlyWorkspace::new(),
        query,
        substituted,
        &substitution_corridor,
        options,
    )
    .unwrap();
    let full_substitution = AlignmentWorkspace::new()
        .align_with_corridor(query, substituted, &substitution_corridor, options)
        .unwrap();
    assert_eq!(substitution.score, full_substitution.score);

    let query_with_gap = b"AAAACCCCGGGGTTTT";
    let mut target_with_gap = b"AAAACCCC".to_vec();
    target_with_gap.extend(std::iter::repeat_n(b'N', 24));
    target_with_gap.extend_from_slice(b"GGGGTTTT");
    let gap_corridor = corridor(
        [CorridorAnchor::new(0, 0, 4), CorridorAnchor::new(12, 36, 4)],
        8,
    );
    let mut gap_scoring = scoring();
    gap_scoring.gap_open_score = -1;
    gap_scoring.gap_extend_score = 0;
    let gap_options = AlignmentOptions::new(gap_scoring);
    let gap = score_only_corridor_with_scratch(
        &mut ScoreOnlyWorkspace::new(),
        query_with_gap,
        &target_with_gap,
        &gap_corridor,
        gap_options,
    )
    .unwrap();
    let full_gap = AlignmentWorkspace::new()
        .align_with_corridor(query_with_gap, &target_with_gap, &gap_corridor, gap_options)
        .unwrap();
    assert_eq!(gap.score, full_gap.score);
}

#[test]
fn score_only_scratch_is_bounded_by_widest_row() {
    let query = b"ACGT".repeat(2_500);
    let corridor = corridor([CorridorAnchor::new(0, 0, 21)], 64);
    let options = AlignmentOptions::new(scoring());
    let mut scratch = ScoreOnlyWorkspace::new();
    let probe =
        score_only_corridor_with_scratch(&mut scratch, &query, &query, &corridor, options).unwrap();
    assert_eq!(probe.score, query.len() as i32 * scoring().match_score);
    assert!(scratch.max_row_width() <= 2 * scoring().band_width as usize + 1 + 128);
    assert!(scratch.capacity_cells() <= scratch.max_row_width() * 2);
    assert_eq!(scratch.rows_processed(), query.len() + 1);
}

#[test]
fn score_only_retries_traceback_once_at_first_safe_width() {
    let mut probes = Vec::new();
    let outcome = refine_with_score_only(
        &[64, 128, 128, 256],
        |width| {
            probes.push(width);
            Ok(ScoreProbe {
                score: width as i32,
                band_edge_touched: width < 128,
                under_explained: false,
                predicted_drift: false,
            })
        },
        |width| {
            assert_eq!(width, 128);
            Ok(dummy_result())
        },
    )
    .unwrap();

    assert_eq!(probes, vec![64, 128]);
    assert_eq!(outcome.selection.attempted_widths, vec![64, 128]);
    assert_eq!(outcome.selection.retries, 1);
    assert_eq!(outcome.result.retry_metadata.selected_width, 128);
    assert!(!outcome.selection.retry_cap_reached);
}

#[test]
fn retry_selection_marks_final_edge_when_no_width_explains_chain() {
    let result = select_retry_width(&[32, 64, 128], |width| {
        Ok(ScoreProbe {
            score: width as i32,
            band_edge_touched: true,
            under_explained: false,
            predicted_drift: false,
        })
    })
    .unwrap();

    assert_eq!(result.selected_width, 128);
    assert!(result.retry_cap_reached);
    assert!(result.band_edge_touched);
}

#[test]
fn empty_retry_widths_are_explicitly_rejected() {
    assert!(matches!(
        select_retry_width(&[0, 0], |_| {
            Ok(ScoreProbe {
                score: 0,
                band_edge_touched: false,
                under_explained: false,
                predicted_drift: false,
            })
        }),
        Err(AlignmentError::NoAlignment)
    ));
}

#[test]
fn bridge_requires_same_contig_and_strand_and_similar_gaps() {
    let config = BridgeConfig {
        max_query_gap: 100,
        max_target_gap: 100,
        max_gap_difference: 3,
        min_flank_support: 20,
    };
    let left = endpoint(1, Strand::Forward, (10, 40), (100, 130));
    let right = endpoint(1, Strand::Forward, (60, 90), (153, 183));
    let bridge = bridge_chains(left, right, 300, config).unwrap();
    assert_eq!(bridge.query_gap_bases, 20);
    assert_eq!(bridge.target_gap_bases, 23);
    assert_eq!(bridge.query_gap, interval(40, 60));
    assert_eq!(bridge.target_gap, interval(130, 153));

    assert!(
        bridge_chains(
            left,
            endpoint(2, Strand::Forward, (60, 90), (153, 183)),
            300,
            config
        )
        .is_none()
    );
    assert!(
        bridge_chains(
            left,
            endpoint(1, Strand::Reverse, (60, 90), (117, 147)),
            300,
            config
        )
        .is_none()
    );
    assert!(
        bridge_chains(
            left,
            endpoint(1, Strand::Forward, (60, 90), (190, 220)),
            300,
            config
        )
        .is_none()
    );
}

#[test]
fn reverse_bridge_is_compared_in_oriented_order() {
    let config = BridgeConfig {
        max_query_gap: 100,
        max_target_gap: 100,
        max_gap_difference: 0,
        min_flank_support: 21,
    };
    // Forward target intervals descend as the query advances on the reverse
    // strand: [200,230] then [150,180] in a 300-base contig.
    let left = endpoint(1, Strand::Reverse, (10, 40), (200, 230));
    let right = endpoint(1, Strand::Reverse, (60, 90), (150, 180));
    let bridge = bridge_chains(left, right, 300, config).unwrap();
    assert_eq!(bridge.target_gap, interval(180, 200));
    assert_eq!(bridge.target_gap_bases, 20);
}
