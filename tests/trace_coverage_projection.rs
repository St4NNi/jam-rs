use jam_rs::trace::coverage::{
    CoverageError, parse_cigar, project_alignment, project_cigar, project_edit_script,
    summarize_coverage,
};
use jam_rs::trace::intervals::circular_gap_complement;
use jam_rs::trace::model::{BaseAlignment, BaseInterval, EditOperation, EditRun, Strand};

#[test]
fn cigar_projection_excludes_internal_deletion_from_supported_query() {
    let projection = project_cigar(&[BaseInterval { start: 10, end: 19 }], "3=2I2D4X").unwrap();
    assert_eq!(projection.query_consumed, 9);
    assert_eq!(projection.target_consumed, 9);
    assert_eq!(projection.supported_bases, 7);
    assert_eq!(projection.insertions, 2);
    assert_eq!(projection.deletions, 2);
    assert_eq!(
        projection.deletion_intervals,
        vec![BaseInterval { start: 13, end: 15 }]
    );
    assert_eq!(
        projection.supported_intervals,
        vec![
            BaseInterval { start: 10, end: 13 },
            BaseInterval { start: 15, end: 19 }
        ]
    );
}

#[test]
fn circular_projection_preserves_two_ordered_query_segments() {
    let projection = project_edit_script(
        &[
            BaseInterval { start: 8, end: 10 },
            BaseInterval { start: 0, end: 4 },
        ],
        &[
            EditRun {
                operation: EditOperation::Equal,
                length: 2,
            },
            EditRun {
                operation: EditOperation::Deletion,
                length: 1,
            },
            EditRun {
                operation: EditOperation::Equal,
                length: 3,
            },
        ],
    )
    .unwrap();
    assert_eq!(projection.query_consumed, 6);
    assert_eq!(
        projection.supported_intervals,
        vec![
            BaseInterval { start: 8, end: 10 },
            BaseInterval { start: 1, end: 4 }
        ]
    );
    assert_eq!(
        projection.deletion_intervals,
        vec![BaseInterval { start: 0, end: 1 }]
    );
}

#[test]
fn malformed_cigar_and_length_mismatch_are_explicit() {
    assert!(matches!(
        parse_cigar("3"),
        Err(CoverageError::MissingLength { .. })
    ));
    assert!(matches!(
        parse_cigar("0="),
        Err(CoverageError::ZeroLengthRun { .. })
    ));

    let alignment = BaseAlignment {
        plasmid_id: "p".to_string(),
        metagenome_id: "m".to_string(),
        contig_id: "c".to_string(),
        strand: Strand::Forward,
        query_segments: vec![BaseInterval { start: 0, end: 4 }],
        target_interval: BaseInterval { start: 0, end: 3 },
        query_length: 4,
        target_length: 3,
        origin_crossing: false,
        score: 1,
        matches: 3,
        substitutions: 0,
        insertions: 0,
        deletions: 1,
        cigar: "3=".to_string(),
        edit_script: vec![EditRun {
            operation: EditOperation::Equal,
            length: 3,
        }],
        chain_score: 1,
        primary: true,
    };
    let error = summarize_coverage(10, &[alignment], &[]).unwrap_err();
    assert!(matches!(error, CoverageError::QueryLengthMismatch { .. }));
}

#[test]
fn insertion_and_deletion_projection_keeps_supported_union_bounded() {
    let projection = project_cigar(&[BaseInterval { start: 0, end: 9 }], "2=2I2D5=").unwrap();
    assert_eq!(projection.query_consumed, 9);
    assert_eq!(projection.target_consumed, 9);
    assert_eq!(projection.supported_bases, 7);
    assert_eq!(projection.aligned_query_bases, 9);
    assert_eq!(
        projection.deletion_intervals,
        vec![BaseInterval { start: 2, end: 4 }]
    );
    assert_eq!(
        projection.supported_intervals,
        vec![
            BaseInterval { start: 0, end: 2 },
            BaseInterval { start: 4, end: 9 }
        ]
    );
}

#[test]
fn alignment_target_coordinates_are_checked_before_coverage() {
    let mut alignment = BaseAlignment {
        plasmid_id: "p".to_string(),
        metagenome_id: "m".to_string(),
        contig_id: "c".to_string(),
        strand: Strand::Forward,
        query_segments: vec![BaseInterval { start: 0, end: 3 }],
        target_interval: BaseInterval { start: 8, end: 4 },
        query_length: 3,
        target_length: 0,
        origin_crossing: false,
        score: 3,
        matches: 3,
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        cigar: "3=".to_string(),
        edit_script: vec![EditRun {
            operation: EditOperation::Equal,
            length: 3,
        }],
        chain_score: 3,
        primary: true,
    };
    assert!(matches!(
        project_alignment(&alignment),
        Err(CoverageError::ReversedTargetInterval { .. })
    ));

    alignment.target_interval = BaseInterval { start: 8, end: 12 };
    assert!(matches!(
        project_alignment(&alignment),
        Err(CoverageError::TargetIntervalLengthMismatch { .. })
    ));
}

#[test]
fn overlapping_query_segments_are_rejected_before_projection() {
    let error = project_cigar(
        &[
            BaseInterval { start: 0, end: 3 },
            BaseInterval { start: 1, end: 4 },
        ],
        "6=",
    )
    .unwrap_err();
    assert!(matches!(
        error,
        CoverageError::OverlappingQuerySegments { .. }
    ));
}

#[test]
fn zero_complete_and_separate_coverage_have_consistent_gaps() {
    let zero = jam_rs::trace::coverage::summary_from_intervals(10, &[], &[]).unwrap();
    assert_eq!(zero.supported_bases, 0);
    assert_eq!(zero.supported_fraction, 0.0);
    assert_eq!(zero.gaps.len(), 1);
    assert_eq!(zero.gaps[0].interval, BaseInterval { start: 0, end: 10 });

    let complete = jam_rs::trace::coverage::summary_from_intervals(
        10,
        &[
            BaseInterval { start: 0, end: 6 },
            BaseInterval { start: 6, end: 10 },
        ],
        &[],
    )
    .unwrap();
    assert_eq!(complete.supported_bases, 10);
    assert!(complete.gaps.is_empty());

    let separate = jam_rs::trace::coverage::summary_from_intervals(
        10,
        &[
            BaseInterval { start: 0, end: 2 },
            BaseInterval { start: 6, end: 8 },
        ],
        &[],
    )
    .unwrap();
    assert_eq!(separate.supported_bases, 4);
    assert_eq!(
        separate
            .gaps
            .iter()
            .map(|gap| gap.length)
            .collect::<Vec<_>>(),
        vec![4, 2]
    );
}

#[test]
fn circular_gap_and_origin_rotation_preserve_total_evidence() {
    let wrapped = jam_rs::trace::coverage::summary_from_intervals(
        10,
        &[
            BaseInterval { start: 8, end: 10 },
            BaseInterval { start: 0, end: 2 },
        ],
        &[],
    )
    .unwrap();
    assert_eq!(wrapped.supported_bases, 4);
    assert_eq!(wrapped.largest_gap, 6);
    assert_eq!(
        circular_gap_complement(
            &[
                BaseInterval { start: 6, end: 10 },
                BaseInterval { start: 0, end: 1 },
            ],
            10,
        )
        .unwrap(),
        vec![BaseInterval { start: 1, end: 6 }]
    );

    // Rotate the stored circular origin by three bases.  The coordinates move,
    // but the supported and unsupported lengths must not.
    let rotated = jam_rs::trace::coverage::summary_from_intervals(
        10,
        &[
            BaseInterval { start: 1, end: 3 },
            BaseInterval { start: 3, end: 5 },
        ],
        &[],
    )
    .unwrap();
    assert_eq!(rotated.supported_bases, wrapped.supported_bases);
    assert_eq!(rotated.supported_fraction, wrapped.supported_fraction);
    assert_eq!(rotated.largest_gap, wrapped.largest_gap);
    assert_eq!(
        rotated.gaps.iter().map(|gap| gap.length).sum::<u64>(),
        wrapped.gaps.iter().map(|gap| gap.length).sum::<u64>()
    );
}

#[test]
fn largest_gap_joins_serialized_pieces_across_origin() {
    let summary = jam_rs::trace::coverage::summary_from_intervals(
        10,
        &[BaseInterval { start: 2, end: 6 }],
        &[],
    )
    .unwrap();
    assert_eq!(
        summary
            .gaps
            .iter()
            .map(|gap| gap.interval)
            .collect::<Vec<_>>(),
        vec![
            BaseInterval { start: 0, end: 2 },
            BaseInterval { start: 6, end: 10 },
        ]
    );
    assert_eq!(summary.largest_gap, 6);
}
