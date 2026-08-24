use jam_rs::trace::alignment::exact_blocks::ExactBlock;
use jam_rs::trace::alignment::refine::{
    BridgeConfig, ChainEndpoint, ExactBlockRefinementConfig, ScoreOnlyWorkspace, ScoreProbe,
    bridge_chains, refine_with_exact_blocks, refine_with_score_only,
    score_only_corridor_with_scratch, select_retry_width,
};
use jam_rs::trace::alignment::{
    AlignmentCorridor, AlignmentError, AlignmentOptions, AlignmentResult, AlignmentRetryMetadata,
    AlignmentWorkspace, CorridorAnchor,
};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::model::{BaseInterval, EditOperation, Strand};

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

fn exact_block(
    query_start: u64,
    query_end: u64,
    oriented_target_start: u64,
    oriented_target_end: u64,
    target_length: u64,
    contig_id: u32,
    strand: Strand,
) -> ExactBlock {
    let oriented_target_interval = interval(oriented_target_start, oriented_target_end);
    let target_interval = match strand {
        Strand::Forward => oriented_target_interval,
        Strand::Reverse => interval(
            target_length - oriented_target_end,
            target_length - oriented_target_start,
        ),
    };
    ExactBlock {
        query_interval: interval(query_start, query_end),
        target_interval,
        oriented_target_interval,
        contig_id,
        strand,
        anchor_count: 1,
        minimum_anchor_k: 21,
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[test]
fn exact_blocks_bypass_dp_and_emit_direct_equal_run() {
    let query = b"ACGT".repeat(8);
    let block = exact_block(
        0,
        query.len() as u64,
        0,
        query.len() as u64,
        query.len() as u64,
        4,
        Strand::Forward,
    );
    // A zero matrix budget proves that no affine-gap allocation is attempted
    // when the exact block covers the selected result.
    let options = AlignmentOptions::new(scoring()).with_max_cells(0);
    let result = refine_with_exact_blocks(
        &query,
        &query,
        &[block],
        options,
        ExactBlockRefinementConfig {
            max_gap_bases: 0,
            max_terminal_flank_bases: 0,
        },
    )
    .unwrap();
    assert_eq!(result.cigar, "32=");
    assert_eq!(result.edit_script[0].operation, EditOperation::Equal);
    assert_eq!(result.matches, 32);
}

#[test]
fn exact_block_stitch_aligns_only_bounded_indel_gap() {
    let query = b"ACGT".repeat(16);
    let mut target = query[..21].to_vec();
    target.extend(std::iter::repeat_n(b'N', 4));
    target.extend_from_slice(&query[21..]);
    let blocks = [
        exact_block(0, 21, 0, 21, target.len() as u64, 4, Strand::Forward),
        exact_block(43, 64, 47, 68, target.len() as u64, 4, Strand::Forward),
    ];
    let result = refine_with_exact_blocks(
        &query,
        &target,
        &blocks,
        scoring(),
        ExactBlockRefinementConfig {
            max_gap_bases: 64,
            max_terminal_flank_bases: 0,
        },
    )
    .unwrap();
    assert_eq!(result.query_interval, interval(0, 64));
    assert_eq!(result.target_interval, interval(0, 68));
    assert_eq!(result.cigar, "21=4I43=");
    assert_eq!(result.insertions, 4);
    let (aligned_query, aligned_target) = result.reconstruct(&query, &target).unwrap();
    assert_eq!(
        aligned_query
            .into_iter()
            .filter(|base| *base != b'-')
            .collect::<Vec<_>>(),
        query
    );
    assert_eq!(
        aligned_target
            .into_iter()
            .filter(|base| *base != b'-')
            .collect::<Vec<_>>(),
        target
    );
}

#[test]
fn reverse_exact_block_stitch_maps_target_interval_without_dp() {
    let query = b"ACGT".repeat(8);
    let target_forward = reverse_complement(&query);
    let block = exact_block(
        0,
        query.len() as u64,
        0,
        query.len() as u64,
        target_forward.len() as u64,
        4,
        Strand::Reverse,
    );
    let result = refine_with_exact_blocks(
        &query,
        &target_forward,
        &[block],
        AlignmentOptions::new(scoring()).with_max_cells(0),
        ExactBlockRefinementConfig {
            max_gap_bases: 0,
            max_terminal_flank_bases: 0,
        },
    )
    .unwrap();
    assert_eq!(result.strand, Strand::Reverse);
    assert_eq!(result.target_interval, interval(0, 32));
    assert_eq!(result.cigar, "32=");
}

#[test]
fn exact_block_stitch_rejects_contig_or_strand_changes() {
    let query = b"ACGT".repeat(8);
    let target = query.clone();
    let first = exact_block(0, 21, 0, 21, 32, 4, Strand::Forward);
    let other_contig = exact_block(21, 32, 21, 32, 32, 5, Strand::Forward);
    let contig_error = refine_with_exact_blocks(
        &query,
        &target,
        &[first.clone(), other_contig],
        scoring(),
        ExactBlockRefinementConfig::default(),
    )
    .unwrap_err();
    assert!(matches!(
        contig_error,
        jam_rs::trace::alignment::refine::ExactBlockRefinementError::ContigMismatch { .. }
    ));

    let other_strand = exact_block(21, 32, 0, 11, 32, 4, Strand::Reverse);
    let strand_error = refine_with_exact_blocks(
        &query,
        &target,
        &[first, other_strand],
        scoring(),
        ExactBlockRefinementConfig::default(),
    )
    .unwrap_err();
    assert!(matches!(
        strand_error,
        jam_rs::trace::alignment::refine::ExactBlockRefinementError::StrandMismatch
    ));
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
