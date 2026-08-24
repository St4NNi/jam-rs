use jam_rs::jamhash_u64_v1;
use jam_rs::trace::alignment::exact_blocks::{
    ExactBlockError, PackedTwoBit, compare_oriented_range, extend_anchor, first_mismatch,
    first_mismatch_reverse_complement,
};
use jam_rs::trace::model::{EditOperation, Strand};
use needletail::Sequence;

fn canonical(sequence: &[u8]) -> u64 {
    sequence
        .normalize(false)
        .bit_kmers(sequence.len() as u8, true)
        .next()
        .unwrap()
        .1
        .0
}

fn anchor(
    query: &[u8],
    query_position: usize,
    target_position: usize,
    contig_id: u32,
    strand: Strand,
    k: usize,
) -> jam_rs::trace::anchors::Anchor {
    let query_kmer = &query[query_position..query_position + k];
    let packed = canonical(query_kmer);
    jam_rs::trace::anchors::Anchor {
        query_position: query_position as u64,
        target_position: target_position as u64,
        contig_id,
        strand,
        k: k as u8,
        hash: jamhash_u64_v1(packed),
        canonical_kmer: packed,
        query_reverse: false,
        target_reverse: strand == Strand::Reverse,
    }
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

#[test]
fn packed_comparison_uses_word_path_and_reports_first_mismatch() {
    let mut left = vec![b'A'; 64];
    let right = left.clone();
    assert_eq!(first_mismatch(&left, &right), 64);
    left[37] = b'C';
    assert_eq!(first_mismatch(&left, &right), 37);

    let packed = PackedTwoBit::from_bytes(b"ACGTacgt");
    let same = PackedTwoBit::from_bytes(b"ACGTACGT");
    assert!(packed.compare_range(0, &same, 0, 8));
}

#[test]
fn exact_extension_reconstructs_full_forward_run_and_emits_only_equal() {
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let target = b"TTTACGTACGTACGTACGTACGTACGTACGTACGTAAA";
    let key = anchor(query, 5, 8, 11, Strand::Forward, 21);
    let block = extend_anchor(query, target, key).unwrap();

    assert_eq!(block.query_interval.start, 0);
    assert_eq!(block.query_interval.end, query.len() as u64);
    assert_eq!(block.target_interval.start, 3);
    assert_eq!(block.target_interval.end, 3 + query.len() as u64);
    assert_eq!(
        block.edit_runs_checked().unwrap()[0].operation,
        EditOperation::Equal
    );
    assert_eq!(block.cigar(), "32=");
}

#[test]
fn reverse_complement_extension_maps_forward_target_coordinates() {
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let reverse = reverse_complement(query);
    let mut target = b"GGG".to_vec();
    target.extend_from_slice(&reverse);
    target.extend_from_slice(b"CCC");
    let key = anchor(query, 2, 12, 7, Strand::Reverse, 21);
    let block = extend_anchor(query, &target, key).unwrap();

    assert_eq!(
        block.query_interval,
        jam_rs::trace::model::BaseInterval::new(0, 32).unwrap()
    );
    assert_eq!(
        block.target_interval,
        jam_rs::trace::model::BaseInterval::new(3, 35).unwrap()
    );
    assert_eq!(block.oriented_target_interval.start, 3);
    assert_eq!(first_mismatch_reverse_complement(query, &target[3..35]), 32);
    assert!(compare_oriented_range(
        query,
        &target,
        0,
        3,
        32,
        Strand::Reverse
    ));
}

#[test]
fn extension_stops_at_indel_and_sequence_edges() {
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let target = b"ACGTTACGTACGTACGTACGTACGTACGTACGT";
    // The inserted T at target position 4 matches query position 3, so the
    // maximal byte-equal run extends through that coincidental base and then
    // stops at the preceding mismatch.
    let key = anchor(query, 5, 6, 1, Strand::Forward, 21);
    let block = extend_anchor(query, target, key).unwrap();
    assert_eq!(
        block.query_interval,
        jam_rs::trace::model::BaseInterval::new(3, 32).unwrap()
    );
    assert_eq!(
        block.target_interval,
        jam_rs::trace::model::BaseInterval::new(4, 33).unwrap()
    );

    let edge_key = anchor(query, 1, 1, 1, Strand::Forward, 21);
    let edge = extend_anchor(query, query, edge_key).unwrap();
    assert_eq!(edge.query_interval.start, 0);
    assert_eq!(edge.query_interval.end, query.len() as u64);
}

#[test]
fn malformed_hash_or_packed_key_never_becomes_an_exact_block() {
    let query = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
    let mut key = anchor(query, 0, 0, 0, Strand::Forward, 21);
    key.hash ^= 1;
    assert!(matches!(
        extend_anchor(query, query, key),
        Err(ExactBlockError::HashMismatch { .. })
    ));
}
