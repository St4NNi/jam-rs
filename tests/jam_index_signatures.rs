use jam_rs::jam_index::{MetagenomeSignatureBuilder, ScreenSelectionPolicy};

fn sequence(length: usize) -> Vec<u8> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            b"ACGT"[(state as usize) & 3]
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[test]
fn standalone_160_base_contig_receives_the_minimum_budget() {
    let policy = ScreenSelectionPolicy::default_signatures();
    let mut builder = MetagenomeSignatureBuilder::new(policy).unwrap();
    let signature = builder.add_contig(&sequence(160)).unwrap();
    assert_eq!(signature.requested_budget, 16);
    assert_eq!(signature.hashes.len(), 16);
    assert!(signature.hashes.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(signature.hashes.iter().all(|hash| *hash != 0));
}

#[test]
fn contig_minhash_is_reverse_complement_invariant_and_bounded() {
    let sequence = sequence(20_000);
    let reverse = reverse_complement(&sequence);
    let policy = ScreenSelectionPolicy::default_signatures();
    let mut forward_builder = MetagenomeSignatureBuilder::new(policy.clone()).unwrap();
    let mut reverse_builder = MetagenomeSignatureBuilder::new(policy).unwrap();
    let forward = forward_builder.add_contig(&sequence).unwrap();
    let reverse = reverse_builder.add_contig(&reverse).unwrap();
    assert_eq!(forward.hashes, reverse.hashes);
    assert!(forward.hashes.len() <= 256);
}

#[test]
fn whole_metagenome_sketch_is_fixed_and_union_is_deduplicated() {
    let policy = ScreenSelectionPolicy::smaller_signatures();
    let mut builder = MetagenomeSignatureBuilder::new(policy).unwrap();
    let first = builder.add_contig(&sequence(30_000)).unwrap();
    let second = builder.add_contig(&sequence(40_000)).unwrap();
    let result = builder.finish();
    assert_eq!(result.contig_count, 2);
    assert_eq!(result.total_bases, 70_000);
    assert!(result.whole_metagenome_hashes.len() <= 256);
    assert!(result.union_hashes.len() <= first.hashes.len() + second.hashes.len() + 256);
    assert!(result.union_hashes.windows(2).all(|pair| pair[0] < pair[1]));
}
