use jam_rs::router::witness::fixed_bottom_k_witnesses;

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
fn fixed_bottom_k_is_bounded_and_deterministic() {
    let sequence = sequence(4_000);
    let first = fixed_bottom_k_witnesses(&sequence, 16, 256).unwrap();
    let second = fixed_bottom_k_witnesses(&sequence, 16, 256).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.len(), 16);
    let ranks = first
        .iter()
        .map(|witness| (witness.key.jamhash, witness.key.packed))
        .collect::<Vec<_>>();
    assert!(ranks.windows(2).all(|pair| pair[0] < pair[1]));
}
