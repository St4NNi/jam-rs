use jam_rs::core_utils::{passes_entropy_filter, shannon_entropy};

#[test]
fn test_shannon_entropy() {
    // Test entropy calculation
    let entropy1 = shannon_entropy(0, 21); // All A's encoded as 0s
    let entropy2 = shannon_entropy(u64::MAX >> 22, 21); // Mixed sequence

    // Entropy should be non-negative
    assert!(entropy1 >= 0.0);
    assert!(entropy2 >= 0.0);

    // Test different k-mer lengths
    for kmer_len in [15, 21, 31] {
        let entropy = shannon_entropy(0x12345678, kmer_len);
        assert!(
            entropy >= 0.0,
            "Entropy should be non-negative for k-mer length {}",
            kmer_len
        );
    }
}

#[test]
fn test_passes_entropy_filter() {
    let kmer_length = 21;

    // Test with different entropy thresholds
    assert!(passes_entropy_filter(0x12345678, kmer_length, 0.0)); // Should always pass with 0 threshold

    // Test that low-complexity k-mers (like all A's) fail high entropy thresholds
    let low_complexity_kmer = 0; // All A's
    assert!(passes_entropy_filter(low_complexity_kmer, kmer_length, 0.0));

    // High entropy threshold should filter out low-complexity k-mers
    assert!(!passes_entropy_filter(low_complexity_kmer, kmer_length, 1.5));
}

#[test]
fn test_entropy_edge_cases() {
    // Test with k-mer length 1
    let entropy = shannon_entropy(1, 1); // Single nucleotide
    assert!(entropy >= 0.0);

    // Test with maximum practical k-mer length
    let entropy_max = shannon_entropy(0, 31);
    assert!(entropy_max >= 0.0);

    // Test passes_entropy_filter with very high threshold
    assert!(!passes_entropy_filter(0, 21, 10.0)); // Should fail with impossible threshold

    // Test with very low threshold
    assert!(passes_entropy_filter(0, 21, -1.0)); // Should pass with negative threshold
}

#[test]
fn test_homopolymer_entropy() {
    // Homopolymer should have zero entropy
    let homopolymer = 0; // All A's
    assert_eq!(shannon_entropy(homopolymer, 8), 0.0);
}

#[test]
fn test_max_entropy() {
    // Maximum entropy for 4-mer with all different bases (ATCG = 2.0 bits)
    let all_diff = 0b11100100; // TCGA
    let max_entropy = shannon_entropy(all_diff, 4);
    assert!((max_entropy - 2.0).abs() < 0.01);
}
