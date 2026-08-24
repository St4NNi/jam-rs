use jam_rs::trace::seeds::gear::{
    ExactFragment, FragmentMiniSketchIndex, FragmentOrientation, MiniSketchConfig,
    build_fragment_mini_sketch, lookup_related_fragments,
};

fn dna(length: usize) -> Vec<u8> {
    let mut state = 0xfeed_face_u64;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[((state >> 32) & 3) as usize]
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
fn overlapping_related_fragments_share_internal_exact_seed_votes() {
    let sequence = dna(500);
    let target = ExactFragment {
        contig_id: 3,
        start: 50,
        length: 220,
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&sequence[50..270]),
    };
    let query = ExactFragment {
        contig_id: 0,
        start: 0,
        length: 220,
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&sequence[..220]),
    };
    assert_ne!(query.digest, target.digest);
    let config = MiniSketchConfig::k21(0);
    let target_sketch = build_fragment_mini_sketch(&sequence, target, config).unwrap();
    let query_sketch = build_fragment_mini_sketch(&sequence, query, config).unwrap();
    assert!(!target_sketch.seeds.is_empty());
    assert!(!query_sketch.seeds.is_empty());
    let mut index = FragmentMiniSketchIndex::default();
    index.insert(target_sketch);
    let related = lookup_related_fragments(&index, &query_sketch, 2);
    assert_eq!(related.len(), 1);
    assert!(related[0].votes >= 2);
    assert!(!related[0].local_pairs.is_empty());
}

#[test]
fn related_fragment_pairs_preserve_forward_source_orientation() {
    let sequence = dna(500);
    let target_bytes = reverse_complement(&sequence[50..270]);
    let target = ExactFragment {
        contig_id: 3,
        start: 50,
        length: 220,
        orientation: FragmentOrientation::Reverse,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&target_bytes),
    };
    let query = ExactFragment {
        contig_id: 0,
        start: 0,
        length: 220,
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&sequence[..220]),
    };
    let config = MiniSketchConfig::k21(0);
    let target_sketch = build_fragment_mini_sketch(&sequence, target, config).unwrap();
    let query_sketch = build_fragment_mini_sketch(&sequence, query, config).unwrap();
    let mut index = FragmentMiniSketchIndex::default();
    index.insert(target_sketch);
    let related = lookup_related_fragments(&index, &query_sketch, 2);
    assert_eq!(related.len(), 1);
    assert!(
        related[0]
            .local_pairs
            .iter()
            .all(|pair| { pair.query_reverse == pair.target_reverse && !pair.relative_reverse })
    );
}

#[test]
fn mini_sketch_rejects_k_below_21_and_zero_scale() {
    let sequence = dna(200);
    let fragment = ExactFragment {
        contig_id: 0,
        start: 0,
        length: 200,
        orientation: FragmentOrientation::Forward,
        digest: jam_rs::trace::seeds::gear::digest_fragment(&sequence),
    };
    assert!(matches!(
        build_fragment_mini_sketch(
            &sequence,
            fragment,
            MiniSketchConfig {
                k: 20,
                max_seeds: 0,
                scale: 1,
            },
        ),
        Err(jam_rs::trace::seeds::gear::MiniSketchError::UnsupportedK(
            20
        ))
    ));
    assert!(matches!(
        build_fragment_mini_sketch(
            &sequence,
            fragment,
            MiniSketchConfig {
                k: 21,
                max_seeds: 0,
                scale: 0,
            },
        ),
        Err(jam_rs::trace::seeds::gear::MiniSketchError::ZeroScale)
    ));
}
