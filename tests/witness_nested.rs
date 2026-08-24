use std::collections::BTreeMap;

use jam_rs::router::witness::{DEFAULT_QUERY_WINDOW_SIZES, extract_nested_witnesses};
use jam_rs::router::{HashAlgorithmId, WITNESS_K, WitnessKey, WitnessScheme};
use needletail::Sequence;

fn scheme() -> WitnessScheme {
    WitnessScheme {
        scheme_id: 1,
        k: WITNESS_K,
        base_scale: 20,
        available_scales: vec![20, 50, 100, 200, 500],
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    }
}

fn sequence(length: usize) -> Vec<u8> {
    let alphabet = b"ACGT";
    let mut state = 0x1234_5678_9abc_def0_u64;
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            alphabet[(state as usize) % alphabet.len()]
        })
        .collect()
}

#[test]
fn nested_tiers_reuse_one_exact_key_and_label_query_windows() {
    let result = extract_nested_witnesses(&sequence(1_200), &scheme(), 128).unwrap();
    assert!(!result.witnesses.is_empty());
    assert_eq!(result.query_window_count(), 10);
    assert_eq!(result.eligible_window_count(), 10);
    assert!(result.witnesses.iter().all(|witness| {
        !witness.query_window_ids.is_empty()
            && witness
                .retained_scales
                .windows(2)
                .all(|pair| pair[0] < pair[1])
    }));
    assert!(result.at_scale(500).unwrap().len() <= result.at_scale(20).unwrap().len());
    assert_eq!(result.unique_keys().len(), {
        let mut keys = result
            .witnesses
            .iter()
            .map(|witness| witness.key)
            .collect::<Vec<_>>();
        keys.sort_unstable();
        keys.dedup();
        keys.len()
    });
    assert_eq!(DEFAULT_QUERY_WINDOW_SIZES, [128, 256, 512, 1024]);
}

#[test]
fn reverse_complement_preserves_canonical_key_set() {
    let forward = sequence(700);
    let reverse = forward
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect::<Vec<_>>();
    let mut left = extract_nested_witnesses(&forward, &scheme(), 256)
        .unwrap()
        .unique_keys();
    let mut right = extract_nested_witnesses(&reverse, &scheme(), 256)
        .unwrap()
        .unique_keys();
    left.sort_unstable();
    right.sort_unstable();
    assert_eq!(left, right);
}

#[test]
fn nested_occurrences_preserve_canonical_reverse_orientation() {
    let sequence = sequence(700);
    let normalized = sequence.normalize(false);
    let expected = normalized
        .bit_kmers(WITNESS_K, true)
        .filter_map(|(position, kmer, reverse)| {
            WitnessKey::from_packed(kmer.0)
                .ok()
                .map(|key| ((position as u64, key.packed), reverse))
        })
        .collect::<BTreeMap<_, _>>();
    let result = extract_nested_witnesses(&sequence, &scheme(), 256).unwrap();
    assert!(result.witnesses.iter().any(|witness| witness.query_reverse));
    for witness in result.witnesses {
        assert_eq!(
            expected.get(&(witness.query_position, witness.key.packed)),
            Some(&witness.query_reverse),
            "orientation at query position {}",
            witness.query_position
        );
    }
}

#[test]
fn ambiguous_and_zero_hash_windows_are_excluded() {
    let mut data = vec![b'A'; 21];
    data.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNN");
    data.extend_from_slice(&sequence(300));
    let result = extract_nested_witnesses(&data, &scheme(), 128).unwrap();
    assert!(
        result
            .witnesses
            .iter()
            .all(|witness| { witness.key.packed != 0 && witness.key.jamhash != 0 })
    );
}
