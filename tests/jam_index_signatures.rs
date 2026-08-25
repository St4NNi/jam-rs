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

#[test]
fn spatial_policies_have_no_fixed_sixteen_signature_floor() {
    let sequence = sequence(1_000);
    let mut spatial_256 =
        MetagenomeSignatureBuilder::new(ScreenSelectionPolicy::spatial_256(512)).unwrap();
    let mut spatial_256_two =
        MetagenomeSignatureBuilder::new(ScreenSelectionPolicy::spatial_256_two(512)).unwrap();
    let selected_256 = spatial_256.add_contig(&sequence).unwrap();
    let selected_256_two = spatial_256_two.add_contig(&sequence).unwrap();
    assert_eq!(selected_256.requested_budget, 4);
    assert_eq!(selected_256.hashes.len(), 4);
    assert_eq!(selected_256_two.requested_budget, 8);
    assert_eq!(selected_256_two.hashes.len(), 8);

    let mut short =
        MetagenomeSignatureBuilder::new(ScreenSelectionPolicy::spatial_256(512)).unwrap();
    let selected_short = short.add_contig(&sequence[..160]).unwrap();
    assert_eq!(selected_short.requested_budget, 1);
    assert_eq!(selected_short.hashes.len(), 1);
}

#[test]
fn spatial_signatures_remain_sorted_deduplicated_and_nonzero() {
    let repeated = b"ACGT".repeat(1_000);
    let mut builder =
        MetagenomeSignatureBuilder::new(ScreenSelectionPolicy::spatial_256_two(1_024)).unwrap();
    let selected = builder.add_contig(&repeated).unwrap();
    assert!(selected.hashes.len() < selected.requested_budget as usize);
    assert!(selected.hashes.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(selected.hashes.iter().all(|hash| *hash != 0));
}

#[test]
fn adaptive_spatial_policy_adds_second_minima_only_at_threshold() {
    for threshold in [512, 768, 1_024] {
        let policy = ScreenSelectionPolicy::spatial_256_adaptive(threshold, 512).unwrap();
        let mut below = MetagenomeSignatureBuilder::new(policy.clone()).unwrap();
        let below = below
            .add_contig(&sequence(usize::try_from(threshold - 1).unwrap()))
            .unwrap();
        assert_eq!(
            below.requested_budget,
            u32::try_from((threshold - 1).div_ceil(256)).unwrap()
        );

        let mut at = MetagenomeSignatureBuilder::new(policy).unwrap();
        let at = at
            .add_contig(&sequence(usize::try_from(threshold).unwrap()))
            .unwrap();
        assert_eq!(
            at.requested_budget,
            u32::try_from(threshold.div_ceil(256) * 2).unwrap()
        );
    }
}
