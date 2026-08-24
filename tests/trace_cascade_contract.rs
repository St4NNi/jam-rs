use jam_rs::archive::SeedSchemeDescriptor;
use jam_rs::trace::alignment::{
    AlignmentCorridor, AlignmentOptions, AlignmentWorkspace, CorridorAnchor,
};
use jam_rs::trace::chain::{AnchorClass, ChainConfig, WeightedAnchor, chain_weighted_anchors};
use jam_rs::trace::model::{BaseInterval, Strand, TraceStageMetrics};
use jam_rs::trace::runner::cascade::{
    CascadeStageKind, CascadeState, CascadeStopReason, GEAR_FRAGMENT_ALGORITHM_ID,
    advertised_exact_gear_scheme,
};
use jam_rs::trace::seeds::gear::{
    ExactFragment, FragmentOrientation, FragmentationMode, GearConfig, GearTableKind,
    VerifiedFragment, fragment_sequence, merge_exact_runs, verify_exact_fragment,
};

#[test]
fn mixed_k31_k21_chain_is_retained_as_one_path() {
    let anchors = [
        WeightedAnchor::new(
            jam_rs::trace::anchors::Anchor {
                query_position: 10,
                target_position: 110,
                contig_id: 1,
                strand: Strand::Forward,
                k: 31,
                hash: 1,
                canonical_kmer: 2,
                query_reverse: false,
                target_reverse: false,
            },
            AnchorClass::SpecificK31,
        ),
        WeightedAnchor::new(
            jam_rs::trace::anchors::Anchor {
                query_position: 80,
                target_position: 180,
                contig_id: 1,
                strand: Strand::Forward,
                k: 21,
                hash: 3,
                canonical_kmer: 4,
                query_reverse: false,
                target_reverse: false,
            },
            AnchorClass::SpecificK21,
        ),
    ];
    let chains = chain_weighted_anchors(
        &anchors,
        256,
        ChainConfig {
            max_chains: 4,
            min_anchors: 2,
            max_predecessors: 8,
            max_query_gap: 256,
            max_target_gap: 256,
            gap_penalty: 1,
            coordinate_model: jam_rs::trace::model::CoordinateModel::Linear,
        },
    )
    .unwrap();
    assert_eq!(chains.len(), 1);
    assert_eq!(chains[0].anchors.len(), 2);
    assert!(chains[0].has_high_specificity_anchor());
}

#[test]
fn corridor_retry_keeps_retry_metadata() {
    let query = b"ACGT".repeat(80);
    let mut target = query.clone();
    target.splice(80..80, b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".iter().copied());
    let corridor = AlignmentCorridor::new(
        [
            CorridorAnchor::new(0, 0, 31),
            CorridorAnchor::new(100, 132, 31),
        ],
        16,
    )
    .unwrap();
    let mut workspace = AlignmentWorkspace::new();
    let result = workspace
        .align_with_retries(
            &query,
            &target,
            &corridor,
            AlignmentOptions::new(jam_rs::trace::SensitivityConfig::default().alignment),
        )
        .unwrap();
    assert!(result.retry_metadata.corridor_used);
    assert!(!result.retry_metadata.attempted_widths.is_empty());
}

#[test]
fn exact_gear_run_emits_direct_equal_cigar() {
    let sequence = b"ACGT".repeat(160);
    let fragments = fragment_sequence(
        &sequence,
        0,
        GearConfig::fine(GearTableKind::SingleBase, 0),
        FragmentationMode::StrandSymmetric,
    )
    .unwrap();
    let first = fragments[0];
    let run = merge_exact_runs(&[VerifiedFragment {
        query_start: first.start,
        target_axis_start: first.start,
        length: u64::from(first.length),
        contig_id: 1,
        orientation: FragmentOrientation::Forward,
    }])
    .unwrap();
    assert_eq!(run.len(), 1);
    assert_eq!(run[0].direct_cigar(), Some(format!("{}=", first.length)));
}

#[test]
fn exact_gear_reverse_occurrence_is_verified_in_query_orientation() {
    let query = b"ACGTTGCATGTCAGTAGGCATCAGTACCGATG".to_vec();
    let target_forward = query
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect::<Vec<_>>();
    let query_fragment = ExactFragment {
        contig_id: 0,
        start: 0,
        length: u32::try_from(query.len()).unwrap(),
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&query),
    };
    let target_fragment = ExactFragment {
        contig_id: 1,
        start: 0,
        length: u32::try_from(target_forward.len()).unwrap(),
        orientation: FragmentOrientation::Reverse,
        digest: query_fragment.digest,
    };
    assert!(verify_exact_fragment(
        &query,
        query_fragment,
        &target_forward,
        target_fragment
    ));
}

#[test]
fn failed_stage_does_not_discard_previous_support() {
    let mut state = CascadeState::new(100);
    let first = TraceStageMetrics {
        stage: 1,
        name: CascadeStageKind::SparseK31.name().to_string(),
        ..TraceStageMetrics::default()
    };
    assert_eq!(
        state.record_stage(first, [BaseInterval { start: 5, end: 35 }]),
        30
    );
    let failed = TraceStageMetrics {
        stage: 2,
        name: CascadeStageKind::DenseGapK31.name().to_string(),
        ..TraceStageMetrics::default()
    };
    assert_eq!(
        state.record_stage(failed, [BaseInterval { start: 12, end: 30 }]),
        0
    );
    assert_eq!(
        state.stop_after_stage(0, false, true, true),
        CascadeStopReason::NoGain
    );
    assert_eq!(
        state.accepted_intervals(),
        &[BaseInterval { start: 5, end: 35 }]
    );
}

#[test]
fn exact_gear_descriptor_is_the_only_fast_path_gate() {
    let descriptor = SeedSchemeDescriptor {
        scheme_id: 99,
        algorithm_id: GEAR_FRAGMENT_ALGORITHM_ID,
        span: 0,
        informative_bases: 0,
        density_parameter: 0,
        bucket_bits: 8,
        key_encoding: 2,
        occurrence_encoding: 2,
        flags: 0,
    };
    assert_eq!(
        advertised_exact_gear_scheme(&[descriptor]),
        Some(descriptor)
    );
    assert!(
        advertised_exact_gear_scheme(&[SeedSchemeDescriptor {
            algorithm_id: GEAR_FRAGMENT_ALGORITHM_ID,
            key_encoding: 1,
            ..descriptor
        }])
        .is_none()
    );
}
