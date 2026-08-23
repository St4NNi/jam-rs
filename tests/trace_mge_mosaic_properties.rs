use jam_rs::trace::intervals::circular_gap_complement;
use jam_rs::trace::model::{
    AlignmentRole, BaseAlignment, BaseInterval, CoordinateModel, EditOperation, EditRun,
    SeedEvidence, Strand, TopologyEvidence, TopologyRequested,
};
use jam_rs::trace::mosaic::{assess_topology, assign_alignment_ids, summarize_fragment_mosaic};

fn exact(contig: &str, start: u64, end: u64) -> BaseAlignment {
    let length = end - start;
    BaseAlignment {
        alignment_id: String::new(),
        plasmid_id: "query".to_string(),
        metagenome_id: "metagenome".to_string(),
        contig_id: contig.to_string(),
        strand: Strand::Forward,
        query_segments: vec![BaseInterval { start, end }],
        target_interval: BaseInterval {
            start: 0,
            end: length,
        },
        query_length: length,
        target_length: length,
        origin_crossing: false,
        score: (length * 2) as i64,
        matches: length,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: format!("{length}="),
        edit_script: vec![EditRun {
            operation: EditOperation::Equal,
            length: u32::try_from(length).unwrap(),
        }],
        chain_score: (length * 2) as i64,
        identity: 1.0,
        seed_evidence: SeedEvidence::default(),
        primary_supported_bases: 0,
        secondary_supported_bases: 0,
        newly_supported_bases: 0,
        role: AlignmentRole::AlternativeMapping,
        primary: false,
    }
}

fn wrapped(contig: &str, start: u64, end: u64, query_length: u64) -> BaseAlignment {
    let first = BaseInterval {
        start,
        end: query_length,
    };
    let second = BaseInterval { start: 0, end };
    let length = first.len() + second.len();
    let mut alignment = exact(contig, 0, length);
    alignment.query_segments = vec![first, second];
    alignment.origin_crossing = true;
    alignment
}

#[test]
fn canonical_ids_make_atomic_selection_order_independent() {
    let mut forward = vec![exact("b", 4, 10), exact("a", 0, 6)];
    let mut reverse = vec![forward[1].clone(), forward[0].clone()];
    assign_alignment_ids(&mut forward);
    assign_alignment_ids(&mut reverse);

    let left = summarize_fragment_mosaic(10, CoordinateModel::Linear, &forward).unwrap();
    let right = summarize_fragment_mosaic(10, CoordinateModel::Linear, &reverse).unwrap();
    assert_eq!(left, right);
    assert_eq!(left.base_covered_bases, 10);
    assert_eq!(
        left.atomic_intervals
            .iter()
            .map(|atomic| atomic.interval)
            .collect::<Vec<_>>(),
        vec![
            BaseInterval { start: 0, end: 4 },
            BaseInterval { start: 4, end: 6 },
            BaseInterval { start: 6, end: 10 },
        ]
    );
}

#[test]
fn duplicate_and_overlapping_fragments_never_inflate_support() {
    let mut alignments = vec![
        exact("a", 0, 8),
        exact("b", 4, 12),
        exact("duplicate", 0, 8),
    ];
    assign_alignment_ids(&mut alignments);
    let summary = summarize_fragment_mosaic(12, CoordinateModel::Linear, &alignments).unwrap();
    assert_eq!(summary.base_covered_bases, 12);
    assert_eq!(summary.selection_components.redundant_overlap_bases, 12);
    assert_eq!(summary.selection_components.fragment_count, 2);
    assert!(
        summary
            .alignment_evidence
            .iter()
            .any(|evidence| evidence.secondary_supported_bases > 0)
    );
}

#[test]
fn terminal_repeat_evidence_is_labeled_without_double_counting() {
    let mut left = exact("left_repeat", 0, 3);
    left.seed_evidence.repetitive_seed_count = 4;
    let mut right = exact("right_repeat", 9, 12);
    right.seed_evidence.repetitive_seed_count = 4;
    let mut duplicate = left.clone();
    duplicate.contig_id = "duplicate_repeat".to_string();
    duplicate.seed_evidence.repetitive_seed_count = 4;

    let summary =
        summarize_fragment_mosaic(12, CoordinateModel::Linear, &[left, right, duplicate]).unwrap();
    assert_eq!(summary.base_covered_bases, 6);
    assert_eq!(summary.repeat_only_supported_bases, 6);
    assert_eq!(summary.repeat_only_intervals.len(), 2);
    assert!(
        summary
            .alignment_evidence
            .iter()
            .all(|evidence| evidence.role == AlignmentRole::RepeatOnly)
    );
}

#[test]
fn linear_terminal_gaps_remain_distinct_and_nonwrapping() {
    let summary =
        summarize_fragment_mosaic(10, CoordinateModel::Linear, &[exact("middle", 3, 7)]).unwrap();
    assert_eq!(summary.unsupported_gaps.len(), 2);
    assert!(
        summary
            .unsupported_gaps
            .iter()
            .all(|gap| !gap.wraps_origin && gap.segments.len() == 1)
    );
    assert_eq!(
        summary.unsupported_gaps[0].segments[0],
        BaseInterval { start: 0, end: 3 }
    );
    assert_eq!(
        summary.unsupported_gaps[1].segments[0],
        BaseInterval { start: 7, end: 10 }
    );
}

#[test]
fn circular_rotation_preserves_support_and_gap_lengths() {
    let original =
        summarize_fragment_mosaic(12, CoordinateModel::Wrap, &[wrapped("origin", 10, 3, 12)])
            .unwrap();
    let rotated =
        summarize_fragment_mosaic(12, CoordinateModel::Wrap, &[exact("rotated", 1, 6)]).unwrap();
    assert_eq!(original.base_covered_bases, 5);
    assert_eq!(rotated.base_covered_bases, 5);
    assert_eq!(original.aligned_span_bases, rotated.aligned_span_bases);
    assert_eq!(
        original
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
fn circular_gap_complement_preserves_terminal_pieces() {
    let gaps = circular_gap_complement(
        &[
            BaseInterval { start: 2, end: 4 },
            BaseInterval { start: 8, end: 10 },
        ],
        12,
    )
    .unwrap();
    assert_eq!(
        gaps,
        vec![
            BaseInterval { start: 0, end: 2 },
            BaseInterval { start: 4, end: 8 },
            BaseInterval { start: 10, end: 12 },
        ]
    );
    assert_eq!(gaps.iter().map(|gap| gap.len()).sum::<u64>(), 8);
}

#[test]
fn unknown_topology_keeps_linear_display_and_stays_undetermined() {
    let crossing = wrapped("crossing", 10, 2, 12);
    let assessment = assess_topology(
        12,
        TopologyRequested::Unknown,
        2,
        std::slice::from_ref(&crossing),
    )
    .unwrap();
    assert_eq!(assessment.coordinate_model, CoordinateModel::Undetermined);
    assert_eq!(assessment.topology_evidence, TopologyEvidence::Undetermined);
    assert_eq!(assessment.linear_model.mosaic.base_covered_bases, 0);
    assert_eq!(assessment.linear_model.mosaic.unsupported_gaps.len(), 1);
    assert!(assessment.wrap_model.is_some());
    assert_eq!(
        assessment
            .wrap_model
            .as_ref()
            .unwrap()
            .mosaic
            .base_covered_bases,
        4
    );
}
