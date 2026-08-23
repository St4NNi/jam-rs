use jam_rs::trace::config::SeedSensitivity;
use jam_rs::trace::seeds::{
    HASH_ID, PRIMARY_K, RESCUE_K, SeedError, extract_for_sensitivity, extract_seed_level,
    extract_seed_levels, retained_hash,
};

#[test]
fn trace_seed_identity_and_zero_rule_are_explicit() {
    assert_eq!(HASH_ID, "jamhash_u64_v1");
    assert_eq!(PRIMARY_K, 31);
    assert_eq!(RESCUE_K, 21);
    assert!(!retained_hash(0, 100));
}

#[test]
fn all_zero_hash_windows_are_skipped_at_both_supported_kmers() {
    for (k, windows) in [(31, 34u64), (21, 44u64)] {
        let sequence = vec![b'A'; 64];
        let level = extract_seed_level(
            &sequence,
            SeedSensitivity {
                k,
                scale: 1,
                max_occurrences: 1,
            },
        )
        .unwrap();
        assert!(level.seeds.is_empty());
        assert_eq!(level.skipped_hash_zero, windows);
    }
}

#[test]
fn exact_packed_seeds_are_positioned_and_levels_are_nested() {
    let sequence = b"ACGTGCACTGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    let sketch = extract_seed_levels(
        sequence,
        &[
            SeedSensitivity {
                k: 31,
                scale: 500,
                max_occurrences: 32,
            },
            SeedSensitivity {
                k: 31,
                scale: 100,
                max_occurrences: 32,
            },
        ],
    )
    .unwrap();
    assert_eq!(sketch.sequence_length, sequence.len() as u64);
    assert_eq!(sketch.levels.len(), 2);
    assert_eq!(sketch.levels[0].scale, 100);
    assert!(
        sketch.levels[1]
            .seeds
            .iter()
            .all(|seed| sketch.levels[0].seeds.contains(seed))
    );
    assert!(
        sketch.levels[0]
            .seeds
            .iter()
            .all(|seed| seed.position + u64::from(PRIMARY_K) <= sketch.sequence_length)
    );
    assert!(sketch.levels[0].seeds.iter().all(|seed| seed.hash != 0));
}

#[test]
fn sensitivity_extracts_primary_and_optional_rescue_levels() {
    let sketch = extract_for_sensitivity(
        b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        SeedSensitivity {
            k: 31,
            scale: 200,
            max_occurrences: 64,
        },
        Some(SeedSensitivity {
            k: 21,
            scale: 500,
            max_occurrences: 64,
        }),
    )
    .unwrap();
    assert_eq!(sketch.primary().map(|level| level.k), Some(31));
    assert_eq!(sketch.rescue().map(|level| level.k), Some(21));
}

#[test]
fn unsupported_k_and_zero_scale_fail_before_hashing() {
    let invalid_k = extract_seed_level(
        b"ACGT",
        SeedSensitivity {
            k: 11,
            scale: 1,
            max_occurrences: 1,
        },
    );
    assert_eq!(invalid_k, Err(SeedError::UnsupportedK(11)));
    let invalid_scale = extract_seed_level(
        b"ACGT",
        SeedSensitivity {
            k: 21,
            scale: 0,
            max_occurrences: 1,
        },
    );
    assert_eq!(invalid_scale, Err(SeedError::ZeroScale));
}

#[test]
fn duplicate_seed_levels_are_rejected() {
    let config = SeedSensitivity {
        k: 31,
        scale: 1,
        max_occurrences: 1,
    };
    assert_eq!(
        extract_seed_levels(b"ACGTACGTACGTACGTACGTACGTACGTACGT", &[config, config]),
        Err(SeedError::DuplicateLevel)
    );
}
