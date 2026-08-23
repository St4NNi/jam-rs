use jam_rs::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, CoordinateModel, EditOperation, EditRun,
    SeedEvidence, Strand, TopologyEvidence, TopologyRequested,
};
use jam_rs::trace::mosaic::{assess_topology, assign_alignment_ids, summarize_fragment_mosaic};

fn alignment(
    contig: &str,
    query_segments: Vec<BaseInterval>,
    edit_script: Vec<EditRun>,
    target_length: u64,
    score: i64,
) -> BaseAlignment {
    let query_length = query_segments.iter().map(|segment| segment.len()).sum();
    let matches = edit_script
        .iter()
        .filter(|run| run.operation == EditOperation::Equal)
        .map(|run| u64::from(run.length))
        .sum();
    let substitutions = edit_script
        .iter()
        .filter(|run| run.operation == EditOperation::Substitution)
        .map(|run| u64::from(run.length))
        .sum();
    let insertions = edit_script
        .iter()
        .filter(|run| run.operation == EditOperation::Insertion)
        .map(|run| u64::from(run.length))
        .sum();
    let deletions = edit_script
        .iter()
        .filter(|run| run.operation == EditOperation::Deletion)
        .map(|run| u64::from(run.length))
        .sum();
    BaseAlignment {
        alignment_id: format!("{contig}:0-{target_length}"),
        plasmid_id: "query".to_string(),
        metagenome_id: "metagenome".to_string(),
        contig_id: contig.to_string(),
        strand: Strand::Forward,
        query_segments,
        target_interval: BaseInterval {
            start: 0,
            end: target_length,
        },
        query_length,
        target_length,
        origin_crossing: false,
        score,
        matches,
        substitutions,
        insertions,
        deletions,
        cigar: String::new(),
        edit_script,
        chain_score: score,
        identity: if query_length == 0 {
            0.0
        } else {
            matches as f64 / query_length as f64
        },
        seed_evidence: SeedEvidence::default(),
        primary_supported_bases: 0,
        secondary_supported_bases: 0,
        newly_supported_bases: 0,
        role: AlignmentRole::AlternativeMapping,
        primary: false,
    }
}

fn exact(contig: &str, start: u64, end: u64) -> BaseAlignment {
    let length = end - start;
    alignment(
        contig,
        vec![BaseInterval { start, end }],
        vec![EditRun {
            operation: EditOperation::Equal,
            length: u32::try_from(length).unwrap(),
        }],
        length,
        length as i64 * 2,
    )
}

fn with_seed_evidence(
    mut alignment: BaseAlignment,
    nonrepetitive_anchor_count: u64,
    common_anchor_count: u64,
    repetitive_seed_count: u64,
) -> BaseAlignment {
    alignment.seed_evidence = SeedEvidence {
        nonrepetitive_anchor_count,
        common_anchor_count,
        repetitive_seed_count,
        ..SeedEvidence::default()
    };
    alignment
}

fn wrapped(contig: &str, start: u64, end: u64, query_length: u64) -> BaseAlignment {
    let first = BaseInterval {
        start,
        end: query_length,
    };
    let second = BaseInterval { start: 0, end };
    let length = first.len() + second.len();
    let mut result = exact(contig, 0, length);
    result.query_segments = vec![first, second];
    result.origin_crossing = true;
    result
}

#[test]
fn atomic_mosaic_unions_overlap_and_splits_at_boundaries() {
    let summary = summarize_fragment_mosaic(
        10,
        CoordinateModel::Linear,
        &[exact("contig_a", 0, 6), exact("contig_b", 4, 10)],
    )
    .unwrap();
    assert_eq!(summary.base_covered_bases, 10);
    assert_eq!(summary.aligned_span_bases, 10);
    assert_eq!(summary.unsupported_gaps, Vec::new());
    assert_eq!(
        summary
            .atomic_intervals
            .iter()
            .map(|atomic| atomic.interval)
            .collect::<Vec<_>>(),
        vec![
            BaseInterval { start: 0, end: 4 },
            BaseInterval { start: 4, end: 6 },
            BaseInterval { start: 6, end: 10 },
        ]
    );
    assert_eq!(summary.selection_components.redundant_overlap_bases, 2);
}

#[test]
fn duplicate_support_is_retained_as_alternative_without_inflation() {
    let summary = summarize_fragment_mosaic(
        12,
        CoordinateModel::Linear,
        &[exact("repeat_a", 2, 8), exact("repeat_b", 2, 8)],
    )
    .unwrap();
    assert_eq!(summary.base_covered_bases, 6);
    assert_eq!(summary.alternative_alignment_count, 1);
    assert_eq!(summary.selection_components.fragment_count, 1);
    assert_eq!(summary.atomic_intervals.len(), 3);
    assert!(
        summary
            .atomic_intervals
            .iter()
            .any(|atomic| atomic.alternative_alignment_ids.len() == 1)
    );
}

#[test]
fn seed_evidence_populates_category_unions_and_roles() {
    let nonrepetitive = with_seed_evidence(exact("nonrepeat", 0, 4), 3, 0, 0);
    let common = with_seed_evidence(exact("common", 4, 8), 0, 2, 0);
    let repeat = with_seed_evidence(exact("repeat", 8, 10), 0, 0, 4);
    let summary = summarize_fragment_mosaic(
        10,
        CoordinateModel::Linear,
        &[nonrepetitive, common, repeat],
    )
    .unwrap();
    assert_eq!(summary.nonrepetitive_supported_bases, 4);
    assert_eq!(summary.common_sequence_supported_bases, 4);
    assert_eq!(summary.repeat_only_supported_bases, 2);
    assert_eq!(
        summary.common_sequence_intervals,
        vec![BaseInterval { start: 4, end: 8 }]
    );
    assert_eq!(
        summary.repeat_only_intervals,
        vec![BaseInterval { start: 8, end: 10 }]
    );
    assert_eq!(
        summary.selection_components.nonrepetitive_anchor_evidence,
        3
    );
    assert!(
        summary
            .alignment_evidence
            .iter()
            .any(|evidence| evidence.role == AlignmentRole::CommonSequence)
    );
    assert!(
        summary
            .alignment_evidence
            .iter()
            .any(|evidence| evidence.role == AlignmentRole::RepeatOnly)
    );
}

#[test]
fn deletion_is_aligned_span_but_not_base_support_or_gap_semantics() {
    let summary = summarize_fragment_mosaic(
        10,
        CoordinateModel::Linear,
        &[alignment(
            "deletion",
            vec![BaseInterval { start: 0, end: 10 }],
            vec![
                EditRun {
                    operation: EditOperation::Equal,
                    length: 4,
                },
                EditRun {
                    operation: EditOperation::Deletion,
                    length: 2,
                },
                EditRun {
                    operation: EditOperation::Equal,
                    length: 4,
                },
            ],
            8,
            10,
        )],
    )
    .unwrap();
    assert_eq!(summary.base_covered_bases, 8);
    assert_eq!(summary.aligned_span_bases, 10);
    assert_eq!(
        summary.alignment_deletions,
        vec![BaseInterval { start: 4, end: 6 }]
    );
    assert_eq!(
        summary.unsupported_gaps[0].interval,
        BaseInterval { start: 4, end: 6 }
    );
    assert_eq!(summary.unsupported_gaps[0].length, 2);
}

#[test]
fn linear_model_never_wraps_an_origin_crossing_alignment() {
    let crossing = wrapped("origin", 8, 2, 10);
    let summary = summarize_fragment_mosaic(10, CoordinateModel::Linear, &[crossing]).unwrap();
    assert_eq!(summary.coordinate_model, CoordinateModel::Linear);
    assert_eq!(summary.base_covered_bases, 0);
    assert!(summary.covered_intervals.is_empty());
    assert_eq!(
        summary.unsupported_gaps[0].interval,
        BaseInterval { start: 0, end: 10 }
    );
}

#[test]
fn wrap_model_is_origin_invariant_in_coverage_lengths() {
    let first =
        summarize_fragment_mosaic(10, CoordinateModel::Wrap, &[wrapped("origin", 8, 2, 10)])
            .unwrap();
    let rotated =
        summarize_fragment_mosaic(10, CoordinateModel::Wrap, &[exact("rotated", 1, 5)]).unwrap();
    assert_eq!(first.base_covered_bases, 4);
    assert_eq!(rotated.base_covered_bases, 4);
    assert_eq!(first.unsupported_gaps.len(), 1);
    assert!(!first.unsupported_gaps[0].wraps_origin);
    assert_eq!(
        first.unsupported_gaps[0].segments,
        vec![BaseInterval { start: 2, end: 8 }]
    );
    assert_eq!(first.unsupported_gaps[0].length, 6);
    assert_eq!(rotated.unsupported_gaps.len(), 1);
    assert!(rotated.unsupported_gaps[0].wraps_origin);
    let terminal =
        summarize_fragment_mosaic(10, CoordinateModel::Wrap, &[exact("terminal", 4, 6)]).unwrap();
    assert_eq!(terminal.unsupported_gaps.len(), 1);
    assert!(terminal.unsupported_gaps[0].wraps_origin);
    assert_eq!(
        terminal.unsupported_gaps[0].segments,
        vec![
            BaseInterval { start: 0, end: 4 },
            BaseInterval { start: 6, end: 10 },
        ]
    );
    assert_eq!(terminal.unsupported_gaps[0].length, 8);
    assert_eq!(
        first
            .unsupported_gaps
            .iter()
            .map(|gap| gap.length)
            .sum::<u64>(),
        rotated
            .unsupported_gaps
            .iter()
            .map(|gap| gap.length)
            .sum::<u64>()
    );
}

#[test]
fn auto_and_unknown_do_not_turn_coordinate_models_into_biology_claims() {
    let crossing = wrapped("origin", 8, 2, 10);
    let auto = assess_topology(
        10,
        TopologyRequested::Auto,
        1,
        std::slice::from_ref(&crossing),
    )
    .unwrap();
    assert_eq!(auto.coordinate_model, CoordinateModel::Wrap);
    assert_eq!(auto.topology_evidence, TopologyEvidence::WrapSupported);

    let unknown = assess_topology(10, TopologyRequested::Unknown, 1, &[crossing]).unwrap();
    assert_eq!(unknown.coordinate_model, CoordinateModel::Undetermined);
    assert_eq!(unknown.topology_evidence, TopologyEvidence::Undetermined);
    assert!(unknown.wrap_model.is_some());
}

#[test]
fn auto_returns_both_compatible_for_equal_nonwrapping_evidence() {
    let auto = assess_topology(10, TopologyRequested::Auto, 2, &[exact("linear", 2, 8)]).unwrap();
    assert_eq!(auto.coordinate_model, CoordinateModel::Undetermined);
    assert_eq!(auto.topology_evidence, TopologyEvidence::BothCompatible);
}

#[test]
fn reordered_input_has_identical_generated_ids_and_selection() {
    let mut first = vec![exact("b", 4, 10), exact("a", 0, 6)];
    let mut second = vec![first[1].clone(), first[0].clone()];
    assign_alignment_ids(&mut first);
    assign_alignment_ids(&mut second);
    assert_eq!(
        summarize_fragment_mosaic(10, CoordinateModel::Linear, &first).unwrap(),
        summarize_fragment_mosaic(10, CoordinateModel::Linear, &second).unwrap()
    );
}
