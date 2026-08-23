use jam_rs::trace::model::{BaseAlignment, BaseInterval, EditOperation, EditRun, Strand};
use jam_rs::trace::mosaic::{OverlapClass, SupportClass, classify_overlap, select_primary};

fn alignment(contig: &str, start: u64, end: u64) -> BaseAlignment {
    alignment_segments(contig, vec![BaseInterval { start, end }])
}

fn alignment_segments(contig: &str, query_segments: Vec<BaseInterval>) -> BaseAlignment {
    let length = query_segments
        .iter()
        .map(|segment| segment.len())
        .sum::<u64>();
    BaseAlignment {
        plasmid_id: "p".to_string(),
        metagenome_id: "m".to_string(),
        contig_id: contig.to_string(),
        strand: Strand::Forward,
        query_segments,
        target_interval: BaseInterval {
            start: 0,
            end: length,
        },
        query_length: length,
        target_length: length,
        origin_crossing: false,
        score: length as i64,
        matches: length,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: format!("{length}="),
        edit_script: vec![EditRun {
            operation: EditOperation::Equal,
            length: u32::try_from(length).unwrap(),
        }],
        chain_score: length as i64,
        primary: false,
    }
}

#[test]
fn overlapping_contigs_are_union_deduplicated() {
    let alignments = vec![
        alignment("contig_a", 0, 100),
        alignment("contig_b", 80, 160),
    ];
    let result = select_primary(176, &alignments).unwrap();
    assert_eq!(result.coverage.supported_bases, 160);
    assert_eq!(
        result.coverage.primary_intervals,
        vec![BaseInterval { start: 0, end: 160 }]
    );
    assert_eq!(result.primary_indices, vec![0, 1]);
}

#[test]
fn identical_alignment_is_retained_as_secondary_overlap() {
    let alignments = vec![alignment("contig_a", 0, 100), alignment("contig_b", 0, 100)];
    let result = select_primary(176, &alignments).unwrap();
    assert_eq!(result.primary_indices, vec![0]);
    assert_eq!(result.secondary_indices, vec![1]);
    assert_eq!(
        result.classifications[1].support_class,
        SupportClass::SecondaryOverlap
    );
    assert_eq!(
        result.classifications[1].overlap_class,
        OverlapClass::Identical
    );
}

#[test]
fn overlap_class_is_stable_for_reordered_input() {
    let left = alignment("contig_a", 0, 100);
    let right = alignment("contig_b", 80, 160);
    assert_eq!(
        classify_overlap(176, &left, &right).unwrap(),
        OverlapClass::Partial
    );
    let first = select_primary(176, &[left.clone(), right.clone()])
        .unwrap()
        .coverage;
    let second = select_primary(176, &[right, left]).unwrap().coverage;
    assert_eq!(first, second);
}

#[test]
fn zero_and_origin_crossing_coverage_are_bounded_and_gap_aware() {
    let zero = select_primary(12, &[]).unwrap();
    assert_eq!(zero.coverage.supported_bases, 0);
    assert_eq!(zero.coverage.largest_gap, 12);
    assert_eq!(zero.coverage.gaps.len(), 1);

    let wrapped = alignment_segments(
        "origin",
        vec![
            BaseInterval { start: 10, end: 12 },
            BaseInterval { start: 0, end: 4 },
        ],
    );
    let result = select_primary(12, &[wrapped]).unwrap();
    assert_eq!(result.coverage.supported_bases, 6);
    assert_eq!(
        result.coverage.supported_bases,
        12 - result.coverage.largest_gap
    );
    assert_eq!(
        result.coverage.primary_intervals,
        vec![
            BaseInterval { start: 0, end: 4 },
            BaseInterval { start: 10, end: 12 },
        ]
    );
    assert_eq!(
        result.coverage.gaps[0].interval,
        BaseInterval { start: 4, end: 10 }
    );
}

#[test]
fn repeated_support_is_secondary_without_inflating_any_summary() {
    let result = select_primary(
        120,
        &[
            alignment("repeat_a", 20, 80),
            alignment("repeat_b", 20, 80),
            alignment("repeat_c", 20, 80),
        ],
    )
    .unwrap();
    assert_eq!(result.coverage.supported_bases, 60);
    assert_eq!(result.coverage.primary_intervals.len(), 1);
    assert_eq!(result.coverage.secondary_intervals.len(), 1);
    assert_eq!(result.primary_indices, vec![0]);
    assert_eq!(result.secondary_indices, vec![1, 2]);
    assert!(
        result
            .coverage
            .primary_intervals
            .iter()
            .all(|interval| interval.len() <= 120)
    );
}

#[test]
fn query_segment_tie_break_is_independent_of_input_order() {
    let split = alignment_segments(
        "same_contig",
        vec![
            BaseInterval { start: 0, end: 5 },
            BaseInterval { start: 5, end: 10 },
        ],
    );
    let joined = alignment_segments("same_contig", vec![BaseInterval { start: 0, end: 10 }]);

    let forward = vec![joined.clone(), split.clone()];
    let reverse = vec![split.clone(), joined.clone()];
    let first = select_primary(40, &forward).unwrap();
    let second = select_primary(40, &reverse).unwrap();
    assert_eq!(
        forward[first.primary_indices[0]].query_segments,
        split.query_segments
    );
    assert_eq!(
        reverse[second.primary_indices[0]].query_segments,
        split.query_segments
    );
    assert_eq!(first.coverage, second.coverage);
}
