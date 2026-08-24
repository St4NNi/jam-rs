use jam_rs::trace::seeds::gear::{
    ExactFragment, FragmentOrientation, FragmentationMode, GearConfig, GearStream, GearTableKind,
    INDEPENDENT_TABLE_SEEDS, VerifiedFragment, boundary_distribution, fragment_boundaries,
    fragment_bytes, fragment_sequence, merge_exact_runs, strand_symmetric_distribution,
    verify_exact_fragment,
};

fn dna(length: usize, mut state: u64) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(length);
    for _ in 0..length {
        state ^= state << 7;
        state ^= state >> 9;
        state ^= state << 8;
        sequence.push(b"ACGT"[(state & 3) as usize]);
    }
    sequence
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|byte| match byte {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[test]
fn recovered_fine_and_coarse_streams_are_exactly_declared() {
    assert_eq!(GearStream::fine().min_size, 64);
    assert_eq!(GearStream::fine().target_size, 192);
    assert_eq!(GearStream::fine().max_size, 384);
    assert_eq!(GearStream::coarse().min_size, 256);
    assert_eq!(GearStream::coarse().target_size, 768);
    assert_eq!(GearStream::coarse().max_size, 1_536);
}

#[test]
fn fragment_boundaries_are_deterministic_for_all_table_families() {
    let sequence = dna(20_000, 0x1234_5678);
    for kind in [
        GearTableKind::SingleBase,
        GearTableKind::Dinucleotide,
        GearTableKind::PackedFourBase,
    ] {
        let config = GearConfig::fine(kind, 0xfeed_face);
        let first = fragment_boundaries(&sequence, config).unwrap();
        let second = fragment_boundaries(&sequence, config).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.first(), Some(&0));
        assert_eq!(first.last(), Some(&(sequence.len() as u64)));
        assert!(first.windows(2).all(|window| window[1] > window[0]));
        let fragments =
            fragment_sequence(&sequence, 7, config, FragmentationMode::Forward).unwrap();
        let stats = boundary_distribution(&fragments);
        assert!(stats.count > 0);
        assert!(stats.max <= u64::from(config.max_size));
        assert!(
            fragments
                .iter()
                .take(fragments.len().saturating_sub(1))
                .all(|fragment| fragment.length >= config.min_size)
        );
    }
}

#[test]
fn ordinary_random_dna_has_target_calibrated_distribution_for_each_table() {
    let sequence = dna(200_000, 0x1234_5678);
    let target = 192.0;
    let target_u64 = target as u64;
    for kind in [
        GearTableKind::SingleBase,
        GearTableKind::Dinucleotide,
        GearTableKind::PackedFourBase,
    ] {
        let config = GearConfig::fine(kind, 0xfeed_face);
        let fragments =
            fragment_sequence(&sequence, 7, config, FragmentationMode::Forward).unwrap();
        let stats = boundary_distribution(&fragments);
        eprintln!(
            "Gear {:?}: mean={:.3}, p50={}, p95={}, count={}, min={}, max={}",
            kind, stats.mean, stats.p50, stats.p95, stats.count, stats.min, stats.max
        );
        // The last fragment may be short.  The calibrated stream should not
        // exhibit the old target+minimum bias on the ordinary random case.
        assert!(stats.mean > target * 0.80 && stats.mean < target * 1.25);
        assert!(stats.p50 > config.min_size as u64 && stats.p50 < target_u64);
        assert!(stats.p95 <= config.max_size as u64);
    }
}

#[test]
fn independent_table_seeds_are_stable_and_distinct() {
    assert_eq!(INDEPENDENT_TABLE_SEEDS.len(), 8);
    let mut seeds = INDEPENDENT_TABLE_SEEDS.to_vec();
    seeds.sort_unstable();
    seeds.dedup();
    assert_eq!(seeds.len(), INDEPENDENT_TABLE_SEEDS.len());
}

#[test]
fn boundary_policy_handles_sequence_classes_and_ambiguity() {
    let cases = [
        vec![b'A'; 3_000],
        b"ACGT".repeat(800),
        b"G".repeat(3_000),
        b"AT".repeat(1_500),
        b"GC".repeat(1_500),
        {
            let mut value = dna(3_000, 99);
            value[500..560].fill(b'N');
            value
        },
    ];
    for sequence in cases {
        for kind in [
            GearTableKind::SingleBase,
            GearTableKind::Dinucleotide,
            GearTableKind::PackedFourBase,
        ] {
            let config = GearConfig::fine(kind, 42);
            let fragments =
                fragment_sequence(&sequence, 1, config, FragmentationMode::Forward).unwrap();
            assert!(!fragments.is_empty());
            assert!(fragments.iter().all(|fragment| fragment.length > 0));
            assert!(fragments.iter().all(|fragment| {
                fragment
                    .end()
                    .is_some_and(|end| end <= sequence.len() as u64)
            }));
            assert!(fragments.iter().all(|fragment| {
                fragment_bytes(&sequence, *fragment)
                    .map(|bytes| bytes.len() == fragment.length as usize)
                    .unwrap_or(false)
            }));
        }
    }
}

#[test]
fn both_strands_keep_forward_coordinates_and_exact_oriented_bytes() {
    let sequence = dna(2_500, 0xabcdef);
    let config = GearConfig::fine(GearTableKind::SingleBase, 8);
    let records = fragment_sequence(&sequence, 9, config, FragmentationMode::BothStrands).unwrap();
    assert!(
        records
            .iter()
            .any(|fragment| fragment.orientation == FragmentOrientation::Forward)
    );
    assert!(
        records
            .iter()
            .any(|fragment| fragment.orientation == FragmentOrientation::Reverse)
    );
    for fragment in records {
        let bytes = fragment_bytes(&sequence, fragment).unwrap();
        assert_eq!(bytes.len(), fragment.length as usize);
        assert_eq!(
            jam_rs::trace::seeds::gear::digest_fragment(&bytes),
            fragment.digest
        );
    }
}

#[test]
fn strand_symmetric_boundaries_survive_reverse_complement() {
    let sequence = dna(4_000, 0x1020_3040);
    let reverse = reverse_complement(&sequence);
    let config = GearConfig::fine(GearTableKind::Dinucleotide, 77);
    let forward =
        fragment_sequence(&sequence, 1, config, FragmentationMode::StrandSymmetric).unwrap();
    let reversed =
        fragment_sequence(&reverse, 1, config, FragmentationMode::StrandSymmetric).unwrap();
    let mut forward_intervals = forward
        .iter()
        .map(|fragment| (fragment.start, fragment.length))
        .collect::<Vec<_>>();
    let mut reversed_intervals = reversed
        .iter()
        .map(|fragment| {
            (
                sequence.len() as u64 - fragment.start - u64::from(fragment.length),
                fragment.length,
            )
        })
        .collect::<Vec<_>>();
    forward_intervals.sort_unstable();
    reversed_intervals.sort_unstable();
    assert_eq!(forward_intervals, reversed_intervals);
}

#[test]
fn strand_symmetric_union_reports_tiny_fragments_instead_of_hiding_them() {
    let sequence = dna(20_000, 0x1020_3040);
    let config = GearConfig::fine(GearTableKind::SingleBase, 77);
    let distribution = strand_symmetric_distribution(&sequence, config).unwrap();
    assert_eq!(distribution.configured_min, u64::from(config.min_size));
    assert!(distribution.stats.count > 0);
    assert!(distribution.below_min_count <= distribution.stats.count);
    assert!(distribution.below_min_bases <= sequence.len() as u64);
    assert!((0.0..=1.0).contains(&distribution.below_min_fraction));
}

#[test]
fn strand_symmetric_palindrome_uses_forward_orientation_on_exact_tie() {
    let sequence = b"ACGT".repeat(16);
    let config = GearConfig {
        min_size: 64,
        target_size: 64,
        max_size: 64,
        table_kind: GearTableKind::SingleBase,
        table_seed: 1,
    };
    let fragments =
        fragment_sequence(&sequence, 0, config, FragmentationMode::StrandSymmetric).unwrap();
    assert_eq!(fragments.len(), 1);
    assert_eq!(fragments[0].orientation, FragmentOrientation::Forward);
    assert_eq!(
        fragments[0].digest,
        jam_rs::trace::seeds::gear::digest_fragment(&sequence)
    );
}

#[test]
fn digest_collision_never_becomes_exact_fragment_evidence() {
    let query = b"ACGTACGTACGTACGT";
    let target = b"ACGTACGTACGTACGA";
    let query_fragment = ExactFragment {
        contig_id: 1,
        start: 0,
        length: query.len() as u32,
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(query),
    };
    let mut collision = query_fragment;
    collision.start = 0;
    collision.digest = query_fragment.digest;
    assert!(!verify_exact_fragment(
        query,
        query_fragment,
        target,
        collision
    ));
}

#[test]
fn coordinate_overflow_is_reported_instead_of_saturating() {
    let fragment = ExactFragment {
        contig_id: 0,
        start: u64::MAX,
        length: 1,
        orientation: FragmentOrientation::Forward,
        digest: 0,
    };
    assert_eq!(fragment.end(), None);
    let matches = [VerifiedFragment {
        query_start: u64::MAX,
        target_axis_start: 0,
        length: 1,
        contig_id: 0,
        orientation: FragmentOrientation::Forward,
    }];
    assert!(merge_exact_runs(&matches).is_err());
}

#[test]
fn adjacent_verified_fragments_merge_and_gap_lengths_are_preserved() {
    let matches = [
        VerifiedFragment {
            query_start: 10,
            target_axis_start: 100,
            length: 64,
            contig_id: 4,
            orientation: FragmentOrientation::Forward,
        },
        VerifiedFragment {
            query_start: 74,
            target_axis_start: 164,
            length: 64,
            contig_id: 4,
            orientation: FragmentOrientation::Forward,
        },
    ];
    let runs = merge_exact_runs(&matches).unwrap();
    assert_eq!(runs.len(), 1);
    assert_eq!(runs[0].fragment_count, 2);
    assert_eq!(runs[0].verified_bases, 128);
    assert_eq!(runs[0].direct_cigar().as_deref(), Some("128="));
}
