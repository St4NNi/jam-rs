use jam_rs::trace::seeds::spaced::{SpacedSeedConfig, extract_spaced_seeds};
use jam_rs::trace::seeds::strobemer::{StrobemerConfig, extract_strobemers};
use jam_rs::trace::seeds::syncmer::{SyncmerConfig, SyncmerMode, extract_syncmers};

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

#[test]
fn selected_schemes_are_reverse_complement_symmetric() {
    let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
    let reverse = reverse_complement(query);
    for mode in [SyncmerMode::Closed, SyncmerMode::Open] {
        let config = SyncmerConfig {
            mode,
            ..SyncmerConfig::k31(mode)
        };
        assert_eq!(
            extract_syncmers(query, config).unwrap().len(),
            extract_syncmers(&reverse, config).unwrap().len()
        );
    }
    let spaced = SpacedSeedConfig::span31_weight21();
    assert_eq!(
        extract_spaced_seeds(query, spaced).unwrap().len(),
        extract_spaced_seeds(&reverse, spaced).unwrap().len()
    );
    let paired = StrobemerConfig::paired_k21();
    assert_eq!(
        extract_strobemers(query, paired).unwrap().len(),
        extract_strobemers(&reverse, paired).unwrap().len()
    );
}

#[test]
fn bounded_extraction_and_exact_collision_material_are_explicit() {
    let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCATGCATGCA";
    let syncmers = extract_syncmers(
        query,
        SyncmerConfig {
            max_seeds: 3,
            ..SyncmerConfig::k31(SyncmerMode::Closed)
        },
    )
    .unwrap();
    assert!(syncmers.len() <= 3);
    if let Some(seed) = syncmers.first() {
        let mut collision = seed.seed_key();
        collision.verification[0] ^= 1;
        assert_ne!(collision, seed.seed_key());
    }

    let spaced = extract_spaced_seeds(
        query,
        SpacedSeedConfig {
            max_seeds: 2,
            ..SpacedSeedConfig::span31_weight21()
        },
    )
    .unwrap();
    assert!(spaced.len() <= 2);
    let paired = extract_strobemers(
        query,
        StrobemerConfig {
            max_seeds: 2,
            ..StrobemerConfig::paired_k21()
        },
    )
    .unwrap();
    assert!(paired.len() <= 2);
    assert!(
        paired
            .iter()
            .all(|seed| seed.verification_bytes().len() == 14)
    );
}

#[test]
fn ambiguity_never_becomes_exact_alternative_seed() {
    let query = b"ACGTCAGTACCGATGCTAGCTAGGCTAACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    assert!(
        extract_spaced_seeds(query, SpacedSeedConfig::span31_weight21())
            .unwrap()
            .is_empty()
    );
    assert!(
        extract_strobemers(query, StrobemerConfig::paired_k21())
            .unwrap()
            .is_empty()
    );
}
