use std::collections::BTreeSet;

use jam_rs::router::witness::{DEFAULT_QUERY_WINDOW_SIZES, local_bottom_r_witnesses};
use jam_rs::router::{WITNESS_K, WitnessKey};
use needletail::Sequence;

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
fn local_bottom_r_supports_128_256_512_and_1024_windows() {
    let sequence = sequence(2_500);
    for window_size in DEFAULT_QUERY_WINDOW_SIZES {
        for r in [1, 2, 4] {
            let result = local_bottom_r_witnesses(&sequence, window_size, r).unwrap();
            assert_eq!(result.window_size, window_size);
            assert_eq!(result.r, r);
            assert!(result.witnesses.iter().all(|witness| {
                !witness.query_window_ids.is_empty()
                    && witness
                        .query_window_ids
                        .windows(2)
                        .all(|pair| pair[0] < pair[1])
            }));
            assert!(result.witnesses.len() <= result.eligible_window_count as usize * r);
        }
    }
}

#[test]
fn local_bottom_r_skips_ambiguous_positions_without_panicking() {
    let mut sequence = sequence(700);
    sequence[220..260].fill(b'N');
    let result = local_bottom_r_witnesses(&sequence, 128, 2).unwrap();
    assert!(
        result
            .witnesses
            .iter()
            .all(|witness| { witness.key.packed != 0 && witness.key.jamhash != 0 })
    );
}

#[test]
fn local_bottom_r_matches_jamhash_minima_in_every_window() {
    let sequence = sequence(2_500);
    let window_size = 128_u32;
    let r = 4_usize;
    let result = local_bottom_r_witnesses(&sequence, window_size, r).unwrap();
    let normalized = sequence.normalize(false);
    let observations = normalized
        .bit_kmers(WITNESS_K, true)
        .filter_map(|(position, kmer, _reverse)| {
            WitnessKey::from_packed(kmer.0)
                .ok()
                .map(|key| (position as u64, key))
        })
        .collect::<Vec<_>>();

    for window_id in 0..result.eligible_window_count {
        let start = u64::from(window_id) * u64::from(window_size);
        let end = (start + u64::from(window_size)).min(sequence.len() as u64);
        let end_exclusive = end - u64::from(WITNESS_K - 1);
        let mut expected = observations
            .iter()
            .filter(|(position, _)| *position >= start && *position < end_exclusive)
            .map(|(position, key)| (key.jamhash, key.packed, *position))
            .collect::<Vec<_>>();
        expected.sort_unstable();
        let expected = expected
            .into_iter()
            .take(r)
            .map(|(_hash, packed, position)| (packed, position))
            .collect::<BTreeSet<_>>();
        let actual = result
            .witnesses
            .iter()
            .filter(|witness| witness.query_window_ids.contains(&window_id))
            .map(|witness| (witness.key.packed, witness.query_position))
            .collect::<BTreeSet<_>>();
        assert_eq!(actual, expected, "window {window_id}");
    }
}
