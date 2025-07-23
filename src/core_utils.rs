use std::collections::BinaryHeap;

use serde::{Deserialize, Serialize};

/// Enhanced metadata for each hash
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord)]
pub struct HashMetadata {
    /// Index pointing to file/sequence metadata in secondary database
    pub file_index: u32,
    /// GC percentage category (0-255, with focus on 30-70%)
    pub gc_category: u16,
    /// Length category (0-255, with focus on 1kB-500kB)
    pub length_category: u16,
}

impl HashMetadata {
    pub fn new(file_index: u32, gc_percent: f64, length: usize) -> Self {
        Self {
            file_index,
            gc_category: categorize_gc_percent(gc_percent),
            length_category: categorize_length(length),
        }
    }

    /// Pack metadata into u64 for storage
    pub fn pack(&self) -> u64 {
        ((self.file_index as u64) << 32)
            | ((self.gc_category as u64) << 16)
            | (self.length_category as u64)
    }

    /// Unpack metadata from u64
    pub fn unpack(packed: u64) -> Self {
        Self {
            file_index: (packed >> 32) as u32,
            gc_category: ((packed >> 16) & 0xFFFF) as u16,
            length_category: (packed & 0xFFFF) as u16,
        }
    }
}

/// File/sequence metadata stored in secondary database
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileMetadata {
    pub filename: String,
    pub file_size: u64,
    pub sequence_name: String,
    pub sequence_length: usize,
    pub total_sequences: usize,
    pub total_hashes: usize,
}

/// Categorize GC percentage using 0.01% precision
/// Uses categories 0-10000 (0.00% to 100.00%) within u16 range
pub fn categorize_gc_percent(gc_percent: f64) -> u16 {
    let gc = gc_percent.clamp(0.0, 100.0);
    (gc * 100.0).round() as u16
}

/// Categorize sequence length with 200bp per category
pub fn categorize_length(length: usize) -> u16 {
    (length / 200).try_into().unwrap_or(u16::MAX)
}

/// Shannon entropy calculation for k-mers
pub fn shannon_entropy(kmer: u64, kmer_length: u8) -> f64 {
    let mut counts = [0u8; 4]; // A, C, G, T

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

/// Filter k-mers based on Shannon entropy
pub fn passes_entropy_filter(kmer: u64, kmer_length: u8, min_entropy: f64) -> bool {
    shannon_entropy(kmer, kmer_length) >= min_entropy
}

/// Calculate GC content for a sequence
pub fn calculate_gc_content(sequence: &[u8]) -> f64 {
    let mut gc_count = 0;
    let mut total_count = 0;

    for &base in sequence {
        match base.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc_count += 1;
                total_count += 1;
            }
            b'A' | b'T' => {
                total_count += 1;
            }
            _ => {} // Skip ambiguous bases
        }
    }

    if total_count == 0 {
        0.0
    } else {
        (gc_count as f64 / total_count as f64) * 100.0
    }
}

/// Heap for managing top N smallest hashes per sequence
pub struct TopNHashCollector {
    heap: BinaryHeap<(u64, HashMetadata)>, // Max heap, so largest is at top
    max_size: usize,
}

impl TopNHashCollector {
    pub fn new(max_size: usize) -> Self {
        Self {
            heap: BinaryHeap::with_capacity(max_size + 1),
            max_size,
        }
    }

    pub fn add_hash(&mut self, hash: u64, metadata: HashMetadata) {
        self.heap.push((hash, metadata));

        if self.heap.len() > self.max_size {
            self.heap.pop(); // Remove largest hash
        }
    }

    pub fn into_sorted_vec(self) -> Vec<(u64, HashMetadata)> {
        let mut results: Vec<_> = self.heap.into_iter().collect();
        results.sort_by_key(|&(hash, _)| hash); // Sort by hash value
        results
    }

    pub fn len(&self) -> usize {
        self.heap.len()
    }

    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gc_categorization() {
        assert_eq!(categorize_gc_percent(25.0), 2500); // Below 30%
        assert_eq!(categorize_gc_percent(35.0), 3500); // 30-70% range
        assert_eq!(categorize_gc_percent(50.0), 5000); // Mid range
        assert_eq!(categorize_gc_percent(75.0), 7500); // Above 70%
    }

    #[test]
    fn test_length_categorization() {
        assert_eq!(categorize_length(500), 2); // Below 1kB
        assert_eq!(categorize_length(11_000), 11000 / 200); // 1kB-500kB range
        assert_eq!(categorize_length(250_000), (250000usize / 200) as u16); // Mid range
        assert_eq!(categorize_length(15_000_000), u16::MAX); // Above 10MB -> Overflow u16*200
    }

    #[test]
    fn test_metadata_packing() {
        let metadata = HashMetadata::new(12345, 45.5, 150_000);
        let packed = metadata.pack();
        let unpacked = HashMetadata::unpack(packed);

        assert_eq!(metadata.file_index, unpacked.file_index);
        assert_eq!(metadata.gc_category, unpacked.gc_category);
        assert_eq!(metadata.length_category, unpacked.length_category);
    }

    #[test]
    fn test_shannon_entropy() {
        // Homopolymer should have low entropy
        let homopolymer = 0; // All A's
        assert!(shannon_entropy(homopolymer, 8) < 1.0);

        // Mixed sequence should have higher entropy
        let mixed = 0b10011100; // ATGC pattern
        assert!(shannon_entropy(mixed, 4) > 1.5);
    }

    #[test]
    fn test_top_n_collector() {
        let mut collector = TopNHashCollector::new(3);
        let metadata = HashMetadata::new(1, 50.0, 100_000);

        collector.add_hash(5, metadata);
        collector.add_hash(2, metadata);
        collector.add_hash(8, metadata);
        collector.add_hash(1, metadata);

        let results = collector.into_sorted_vec();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].0, 1); // Smallest
        assert_eq!(results[1].0, 2);
        assert_eq!(results[2].0, 5);
    }
}
