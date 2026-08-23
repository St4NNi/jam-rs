use jam_rs::trace::alignment::{AlignmentResult, align, align_circular, align_oriented};
use jam_rs::trace::config::AlignmentScoring;
use jam_rs::trace::coverage::cigar_from_edit_script;
use jam_rs::trace::model::{BaseInterval, EditOperation, Strand};

fn scoring() -> AlignmentScoring {
    AlignmentScoring {
        match_score: 2,
        mismatch_score: -3,
        gap_open_score: -5,
        gap_extend_score: -1,
        band_width: 32,
    }
}

#[test]
fn plasmid_phage_and_synthetic_queries_use_the_same_evidence_contract() {
    let cases = [
        ("plasmid", b"ACGTTGCACTGA".as_slice()),
        ("phage", b"GATTACAGGCCT".as_slice()),
        ("synthetic", b"TCCGATGACGTA".as_slice()),
    ];

    for (query_class, query) in cases {
        let mut target = b"TT".to_vec();
        target.extend_from_slice(query);
        target.extend_from_slice(b"AA");
        let result = align(query, &target, scoring()).unwrap();

        assert_eq!(result.matches, query.len() as u64, "{query_class} query");
        assert_eq!(result.query_interval.len(), query.len() as u64);
        assert!(!result.origin_crossing, "{query_class} query");
        assert_eq!(result.query_segments.len(), 1);
        assert_reconstructible(query, &target, &result);
    }
}

#[test]
fn forward_and_reverse_alignments_reconstruct_cigar_and_edit_script() {
    let query = b"GACTTAGC";
    let forward_target = reverse_complement(query);
    let reverse_target = reverse_complement(&forward_target);

    let forward = align(query, &reverse_target, scoring()).unwrap();
    assert_eq!(forward.strand, Strand::Forward);
    assert_reconstructible(query, &reverse_target, &forward);

    let reverse = align_oriented(query, &forward_target, 0, Strand::Reverse, scoring()).unwrap();
    assert_eq!(reverse.strand, Strand::Reverse);
    assert_eq!(reverse.cigar, forward.cigar);
    assert_eq!(reverse.score, forward.score);
    // `align_oriented` uses the reverse-complement window internally.  This
    // exact-window case keeps its reported forward interval at [0, L), so the
    // same reconstruction contract can be checked for both strands.
    assert_reconstructible(query, &reverse_target, &reverse);
}

#[test]
fn linear_calls_never_emit_wrap_coordinates() {
    let query = b"ACGTTGCA";
    let forward = align(query, query, scoring()).unwrap();
    let reverse_target = reverse_complement(query);
    let reverse = align_oriented(query, &reverse_target, 17, Strand::Reverse, scoring()).unwrap();

    for result in [&forward, &reverse] {
        assert!(!result.origin_crossing);
        assert_eq!(result.query_segments, vec![result.query_interval]);
        assert!(result.query_interval.end <= query.len() as u64);
        assert_eq!(
            result
                .query_segments
                .iter()
                .map(|segment| segment.len())
                .sum::<u64>(),
            result.query_length
        );
    }
}

#[test]
fn explicit_query_rotation_has_equivalent_normalized_support() {
    let canonical = b"ACGTTGCA";
    let origin = 6_u64;
    let target = b"CAACGT";
    let wrapped = align_circular(
        canonical,
        origin,
        target.len() as u64,
        target,
        0,
        Strand::Forward,
        scoring(),
    )
    .unwrap();

    let mut rotated = canonical[origin as usize..].to_vec();
    rotated.extend_from_slice(&canonical[..origin as usize]);
    let linearized = align(&rotated, target, scoring()).unwrap();

    assert!(wrapped.origin_crossing);
    assert_eq!(wrapped.cigar, linearized.cigar);
    assert_eq!(wrapped.score, linearized.score);
    assert_eq!(wrapped.query_length, linearized.query_length);
    assert_eq!(
        wrapped.query_segments,
        vec![
            BaseInterval { start: 6, end: 8 },
            BaseInterval { start: 0, end: 4 },
        ]
    );

    let canonical_support = support_bitmap(&wrapped.query_segments, canonical.len());
    let rotated_support =
        rotated_support_bitmap(linearized.query_interval, origin as usize, canonical.len());
    assert_eq!(canonical_support, rotated_support);
    assert_eq!(
        canonical_support
            .iter()
            .filter(|supported| **supported)
            .count(),
        6
    );
}

#[test]
fn terminal_direct_and_inverted_repeats_do_not_expand_one_alignment() {
    let cases = [
        ("direct", b"GACCTAGGGGGACCTA".as_slice()),
        ("inverted", b"GACCTAGGGGGTAGGTC".as_slice()),
    ];
    let repeat = b"GACCTA";

    for (repeat_kind, query) in cases {
        let result = align(query, repeat, scoring()).unwrap();
        assert_eq!(result.matches, repeat.len() as u64, "{repeat_kind} repeat");
        assert!(result.query_interval.len() <= repeat.len() as u64);
        assert_eq!(result.query_segments.len(), 1);
        assert!(!result.origin_crossing);
        assert!(
            result
                .query_segments
                .iter()
                .all(|segment| segment.end <= query.len() as u64)
        );
    }
}

#[test]
fn every_accepted_span_is_bounded_by_its_query_length() {
    let linear = align(b"ACGTACGT", b"TTACGTACGTAA", scoring()).unwrap();
    assert_linear_span_bound(&linear, 8);

    let circular = align_circular(b"ACGTAC", 4, 4, b"ACAC", 0, Strand::Forward, scoring()).unwrap();
    assert!(
        circular
            .query_segments
            .iter()
            .map(|segment| segment.len())
            .sum::<u64>()
            <= 6
    );
    assert!(
        circular
            .query_segments
            .iter()
            .all(|segment| segment.end <= 6)
    );
}

fn assert_reconstructible(query: &[u8], target: &[u8], result: &AlignmentResult) {
    let (aligned_query, aligned_target) = result.reconstruct(query, target).unwrap();
    assert_eq!(cigar_from_edit_script(&result.edit_script), result.cigar);
    assert_eq!(aligned_query.len(), aligned_target.len());
    assert_eq!(
        aligned_query.len() as u64,
        result.query_length + result.insertions
    );
    assert_eq!(
        aligned_target.len() as u64,
        result.target_length + result.deletions
    );

    let mut column = 0;
    for run in &result.edit_script {
        assert!(run.length > 0);
        for _ in 0..run.length {
            let query_base = aligned_query[column];
            let target_base = aligned_target[column];
            match run.operation {
                EditOperation::Equal => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Substitution => {
                    assert_ne!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                    assert!(!query_base.eq_ignore_ascii_case(&target_base));
                }
                EditOperation::Insertion => {
                    assert_eq!(query_base, b'-');
                    assert_ne!(target_base, b'-');
                }
                EditOperation::Deletion => {
                    assert_ne!(query_base, b'-');
                    assert_eq!(target_base, b'-');
                }
                EditOperation::SoftClip => panic!("local alignments do not emit soft clips"),
            }
            column += 1;
        }
    }
    assert_eq!(column, aligned_query.len());
}

fn assert_linear_span_bound(result: &AlignmentResult, query_length: u64) {
    assert!(!result.origin_crossing);
    assert_eq!(result.query_segments, vec![result.query_interval]);
    assert!(result.query_interval.end <= query_length);
    assert!(
        result
            .query_segments
            .iter()
            .map(|segment| segment.len())
            .sum::<u64>()
            <= query_length
    );
}

fn support_bitmap(intervals: &[BaseInterval], length: usize) -> Vec<bool> {
    let mut support = vec![false; length];
    for interval in intervals {
        for slot in support
            .iter_mut()
            .take(interval.end as usize)
            .skip(interval.start as usize)
        {
            *slot = true;
        }
    }
    support
}

fn rotated_support_bitmap(interval: BaseInterval, origin: usize, length: usize) -> Vec<bool> {
    let mut support = vec![false; length];
    for offset in interval.start as usize..interval.end as usize {
        support[(origin + offset) % length] = true;
    }
    support
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}
