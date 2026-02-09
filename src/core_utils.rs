pub fn shannon_entropy(kmer: u64, kmer_length: u8) -> f64 {
    let mut counts = [0u8; 4];

    for i in 0..kmer_length {
        let nucleotide = (kmer >> (2 * i)) & 0b11;
        counts[nucleotide as usize] += 1;
    }

    let length = kmer_length as f64;
    let mut entropy = 0.0;

    for count in counts.iter() {
        if *count > 0 {
            let p = (*count as f64) / length;
            entropy -= p * p.log2();
        }
    }

    entropy
}

pub fn passes_entropy_filter(kmer: u64, kmer_length: u8, min_entropy: f64) -> bool {
    shannon_entropy(kmer, kmer_length) >= min_entropy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shannon_entropy() {
        let homopolymer = 0; // All A's
        assert!(shannon_entropy(homopolymer, 8) < 1.0);

        let mixed = 0b10011100; // ATGC pattern
        assert!(shannon_entropy(mixed, 4) > 1.5);
    }

    #[test]
    fn test_passes_entropy_filter() {
        let homopolymer = 0;
        assert!(!passes_entropy_filter(homopolymer, 8, 1.5));

        let mixed = 0b10011100;
        assert!(passes_entropy_filter(mixed, 4, 1.5));
    }

    #[test]
    fn test_entropy_edge_cases() {
        assert_eq!(shannon_entropy(0, 1), 0.0);

        let all_diff = 0b11100100; // TCGA
        let max_entropy = shannon_entropy(all_diff, 4);
        assert!((max_entropy - 2.0).abs() < 0.01);
    }
}
