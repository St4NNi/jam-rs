use jam_rs::core_utils::{passes_entropy_filter, shannon_entropy};

#[test]
fn test_shannon_entropy() {
    let entropy1 = shannon_entropy(0, 21);
    let entropy2 = shannon_entropy(u64::MAX >> 22, 21);

    assert!(entropy1 >= 0.0);
    assert!(entropy2 >= 0.0);

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

    assert!(passes_entropy_filter(0x12345678, kmer_length, 0.0));

    let low_complexity_kmer = 0;
    assert!(passes_entropy_filter(low_complexity_kmer, kmer_length, 0.0));

    assert!(!passes_entropy_filter(
        low_complexity_kmer,
        kmer_length,
        1.5
    ));
}

#[test]
fn test_entropy_edge_cases() {
    let entropy = shannon_entropy(1, 1);
    assert!(entropy >= 0.0);

    let entropy_max = shannon_entropy(0, 31);
    assert!(entropy_max >= 0.0);

    assert!(!passes_entropy_filter(0, 21, 10.0));

    assert!(passes_entropy_filter(0, 21, -1.0));
}

#[test]
fn test_homopolymer_entropy() {
    let homopolymer = 0;
    assert_eq!(shannon_entropy(homopolymer, 8), 0.0);
}

#[test]
fn test_max_entropy() {
    let all_diff = 0b11100100;
    let max_entropy = shannon_entropy(all_diff, 4);
    assert!((max_entropy - 2.0).abs() < 0.01);
}
