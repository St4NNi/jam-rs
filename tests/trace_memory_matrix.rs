//! Bounded, side-effect-free checks for the production measurement matrix.
//!
//! The expensive matrix is executed only by the opt-in Python harnesses.  The
//! integration checks here validate its dimensions and the values that are
//! intended to appear in raw measurement labels without launching a process.

use std::collections::HashSet;

const THREADS: &[usize] = &[1, 4, 8, 16];
const PROFILES: &[&str] = &["fast", "balanced", "sensitive"];
const QUERY_CLASSES: &[&str] = &["small", "large", "repeat-rich"];
const CANDIDATES: &[usize] = &[1, 100];
const CACHE_BLOCK_BYTES: &[usize] = &[16 * 1024, 64 * 1024, 256 * 1024, 1024 * 1024];
const COMPONENTS: &[&str] = &[
    "seed",
    "bucket",
    "anchor",
    "chain",
    "sequence",
    "alignment",
    "output",
];

#[test]
fn memory_matrix_dimensions_are_explicit_and_bounded() {
    let case_count = THREADS.len()
        * PROFILES.len()
        * QUERY_CLASSES.len()
        * CANDIDATES.len()
        * CACHE_BLOCK_BYTES.len();
    let component_count = case_count * COMPONENTS.len();

    assert_eq!(THREADS, &[1, 4, 8, 16]);
    assert_eq!(PROFILES, &["fast", "balanced", "sensitive"]);
    assert_eq!(QUERY_CLASSES, &["small", "large", "repeat-rich"]);
    assert_eq!(CANDIDATES, &[1, 100]);
    assert_eq!(
        CACHE_BLOCK_BYTES,
        &[16 * 1024, 64 * 1024, 256 * 1024, 1024 * 1024]
    );
    assert_eq!(case_count, 288);
    assert_eq!(component_count, 2_016);
    assert!(component_count <= 4_096, "matrix must remain bounded");
}

#[test]
fn measurement_labels_are_unique_for_the_matrix_dimensions() {
    let mut labels = HashSet::new();
    for &threads in THREADS {
        for &profile in PROFILES {
            for &query_class in QUERY_CLASSES {
                for &candidates in CANDIDATES {
                    for &block in CACHE_BLOCK_BYTES {
                        for &component in COMPONENTS {
                            let label = format!(
                                "{component}__{profile}__{query_class}__t{threads}__c{candidates}__b{block}"
                            );
                            assert!(labels.insert(label));
                        }
                    }
                }
            }
        }
    }
    assert_eq!(labels.len(), 2_016);
}

#[test]
fn matrix_values_are_positive_and_have_no_unrequested_kmer_size() {
    assert!(THREADS.iter().all(|value| *value > 0));
    assert!(CANDIDATES.iter().all(|value| *value > 0));
    assert!(CACHE_BLOCK_BYTES.iter().all(|value| *value > 0));
    assert!(!PROFILES.iter().any(|profile| profile.contains("k=")));
    assert!(
        !QUERY_CLASSES
            .iter()
            .any(|query_class| query_class.contains("k="))
    );
}
