use jam_rs::jma::SeedOccurrence;
use jam_rs::trace::anchors::{SeedOccurrenceGroup, summarize_occurrence_evidence};
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

#[test]
fn occurrence_evidence_is_exact_keyed_and_marks_common_repetitive_seeds() {
    let repeated = jam_rs::trace::seeds::QuerySeed {
        position: 5,
        hash: 11,
        canonical_kmer: 22,
        reverse: false,
    };
    let second_query_occurrence = jam_rs::trace::seeds::QuerySeed {
        position: 37,
        ..repeated
    };
    let groups = [
        SeedOccurrenceGroup {
            seed: repeated,
            k: 31,
            occurrences: vec![
                SeedOccurrence {
                    contig_id: 0,
                    position: 8,
                    reverse: false,
                },
                SeedOccurrence {
                    contig_id: 1,
                    position: 12,
                    reverse: true,
                },
                SeedOccurrence {
                    contig_id: 1,
                    position: 44,
                    reverse: false,
                },
            ],
        },
        SeedOccurrenceGroup {
            seed: second_query_occurrence,
            k: 31,
            occurrences: vec![SeedOccurrence {
                contig_id: 0,
                position: 64,
                reverse: false,
            }],
        },
        SeedOccurrenceGroup {
            seed: jam_rs::trace::seeds::QuerySeed {
                position: 10,
                hash: 12,
                canonical_kmer: 23,
                reverse: false,
            },
            k: 31,
            occurrences: vec![],
        },
        SeedOccurrenceGroup {
            seed: jam_rs::trace::seeds::QuerySeed {
                position: 1,
                hash: 0,
                canonical_kmer: 99,
                reverse: false,
            },
            k: 31,
            occurrences: vec![SeedOccurrence {
                contig_id: 0,
                position: 0,
                reverse: false,
            }],
        },
    ];
    let evidence = summarize_occurrence_evidence(&groups, 2, 3);
    assert_eq!(evidence.len(), 2);
    assert_eq!(evidence[0].key.canonical_kmer, 22);
    assert_eq!(evidence[0].query_occurrence_count, 2);
    assert_eq!(evidence[0].candidate_occurrence_count, 4);
    assert_eq!(evidence[0].candidate_occurrence_group_count, 2);
    assert!(evidence[0].is_common);
    assert!(evidence[0].is_repetitive);
    assert_eq!(evidence[0].collection_document_frequency, None);
    assert_eq!(evidence[1].candidate_occurrence_count, 0);
    assert!(!evidence[1].is_common);
    assert!(!evidence[1].is_repetitive);
}
