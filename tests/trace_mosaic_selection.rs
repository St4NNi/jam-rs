use jam_rs::trace::model::{BaseAlignment, BaseInterval, EditOperation, EditRun, Strand};
use jam_rs::trace::mosaic::{OverlapClass, SupportClass, classify_overlap, select_primary};

fn alignment(contig: &str, start: u64, end: u64) -> BaseAlignment {
    let length = end - start;
    BaseAlignment {
        plasmid_id: "p".to_string(),
        metagenome_id: "m".to_string(),
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
