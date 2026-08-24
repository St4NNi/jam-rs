use jam_rs::jamhash_u64_v1;
use jam_rs::trace::alignment::exact_blocks::extend_and_merge_anchors;
use jam_rs::trace::anchors::Anchor;
use jam_rs::trace::model::Strand;
use needletail::Sequence;

fn packed(sequence: &[u8]) -> u64 {
    sequence
        .normalize(false)
        .bit_kmers(sequence.len() as u8, true)
        .next()
        .unwrap()
        .1
        .0
}

fn anchor(query: &[u8], position: usize, target_position: usize, contig_id: u32) -> Anchor {
    let key = packed(&query[position..position + 21]);
    Anchor {
        query_position: position as u64,
        target_position: target_position as u64,
        contig_id,
        strand: Strand::Forward,
        k: 21,
        hash: jamhash_u64_v1(key),
        canonical_kmer: key,
        query_reverse: false,
        target_reverse: false,
    }
}

#[test]
fn overlapping_k21_and_k31_style_exact_anchors_merge_once() {
    let query = b"ACGT".repeat(16);
    let target = [b"TT".as_slice(), query.as_slice(), b"AA".as_slice()].concat();
    let anchors = [anchor(&query, 1, 3, 4), anchor(&query, 20, 22, 4)];
    let blocks = extend_and_merge_anchors(&query, &target, &anchors).unwrap();
    assert_eq!(blocks.len(), 1);
    assert_eq!(blocks[0].query_interval.start, 0);
    assert_eq!(blocks[0].query_interval.end, query.len() as u64);
    assert_eq!(blocks[0].anchor_count, 2);
}

#[test]
fn blocks_never_merge_across_contigs_or_unverified_gaps() {
    let query = b"ACGT".repeat(16);
    let target = [b"TT".as_slice(), query.as_slice(), b"AA".as_slice()].concat();
    let different_contig = [anchor(&query, 1, 3, 4), anchor(&query, 20, 22, 5)];
    let blocks = extend_and_merge_anchors(&query, &target, &different_contig).unwrap();
    assert_eq!(blocks.len(), 2);

    let mut gapped_target = b"TT".to_vec();
    gapped_target.extend_from_slice(&query[..32]);
    gapped_target.push(b'N');
    gapped_target.extend_from_slice(&query[32..]);
    gapped_target.extend_from_slice(b"AA");
    let left = anchor(&query, 1, 3, 4);
    let right = anchor(&query, 40, 43, 4);
    let blocks = extend_and_merge_anchors(&query, &gapped_target, &[left, right]).unwrap();
    assert_eq!(blocks.len(), 2);
}
