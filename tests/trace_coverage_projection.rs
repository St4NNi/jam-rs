use jam_rs::trace::coverage::{
    CoverageError, parse_cigar, project_cigar, project_edit_script, summarize_coverage,
};
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
