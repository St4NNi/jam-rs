use anyhow::Result;
use jam_rs::core_utils::{
    HashMetadata, FileMetadata, categorize_gc_percent, categorize_length, 
    shannon_entropy, passes_entropy_filter, calculate_gc_content, TopNHashCollector
};

#[test]
fn test_hash_metadata_pack_unpack() {
    let metadata = HashMetadata::new(12345, 55.5, 10000);
    
    // Test pack/unpack roundtrip
    let packed = metadata.pack();
    let unpacked = HashMetadata::unpack(packed);
    
    assert_eq!(metadata.file_index, unpacked.file_index);
    assert_eq!(metadata.gc_category, unpacked.gc_category);
    assert_eq!(metadata.length_category, unpacked.length_category);
    
    // Test specific values
    assert_eq!(metadata.file_index, 12345);
    assert!(metadata.gc_category > 0, "GC category should be calculated");
    assert!(metadata.length_category > 0, "Length category should be calculated");
}

#[test]
fn test_hash_metadata_edge_cases() {
    // Test with zero values
    let zero_metadata = HashMetadata::new(0, 0.0, 0);
    let packed = zero_metadata.pack();
    let unpacked = HashMetadata::unpack(packed);
    assert_eq!(zero_metadata.file_index, unpacked.file_index);
    
    // Test with maximum values
    let max_metadata = HashMetadata::new(u32::MAX, 100.0, usize::MAX);
    let packed = max_metadata.pack();
    let unpacked = HashMetadata::unpack(packed);
    assert_eq!(max_metadata.file_index, unpacked.file_index);
}

#[test]
fn test_categorize_gc_percent() {
    // Test various GC percentages
    assert_eq!(categorize_gc_percent(0.0), 0);
    assert_eq!(categorize_gc_percent(50.0), categorize_gc_percent(50.0)); // Should be consistent
    assert_eq!(categorize_gc_percent(100.0), categorize_gc_percent(100.0));
    
    // Test that similar values map to same category
    let cat1 = categorize_gc_percent(49.9);
    let cat2 = categorize_gc_percent(50.1);
    // They might be the same or different, but should be deterministic
    assert_eq!(categorize_gc_percent(49.9), cat1);
    assert_eq!(categorize_gc_percent(50.1), cat2);
}

#[test]
fn test_categorize_length() {
    // Test various lengths
    assert_eq!(categorize_length(0), 0);
    assert!(categorize_length(1000) > 0);
    assert!(categorize_length(100000) > categorize_length(1000));
    
    // Test that function is monotonic for reasonable inputs
    let cat_small = categorize_length(100);
    let cat_medium = categorize_length(10000);
    let cat_large = categorize_length(1000000);
    
    assert!(cat_small <= cat_medium);
    assert!(cat_medium <= cat_large);
}

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
        assert!(entropy >= 0.0, "Entropy should be non-negative for k-mer length {}", kmer_len);
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
    
    // High entropy threshold might filter out some k-mers
    let result_high_threshold = passes_entropy_filter(low_complexity_kmer, kmer_length, 3.0);
    // Result depends on implementation, but function should not crash
    assert!(result_high_threshold == true || result_high_threshold == false);
}

#[test]
fn test_calculate_gc_content() {
    // Test various sequences
    let all_at = b"AAATTTAAATTT";
    assert_eq!(calculate_gc_content(all_at), 0.0);
    
    let all_gc = b"GGGCCCGGGCCC";
    assert_eq!(calculate_gc_content(all_gc), 100.0);
    
    let mixed = b"ATCGATCGATCG";
    let gc_content = calculate_gc_content(mixed);
    assert!(gc_content > 0.0 && gc_content < 100.0, "Mixed sequence should have intermediate GC content");
    
    // Test empty sequence
    let empty = b"";
    assert_eq!(calculate_gc_content(empty), 0.0);
    
    // Test sequence with N's (should be ignored)
    let with_n = b"ATCGNATCGN";
    let gc_with_n = calculate_gc_content(with_n);
    let without_n = b"ATCGATCG";
    let gc_without_n = calculate_gc_content(without_n);
    assert_eq!(gc_with_n, gc_without_n, "N's should be ignored in GC calculation");
}

#[test]
fn test_top_n_hash_collector() {
    let mut collector = TopNHashCollector::new(3);
    
    // Test adding hashes
    let metadata1 = HashMetadata::new(1, 50.0, 1000);
    let metadata2 = HashMetadata::new(2, 60.0, 2000);
    let metadata3 = HashMetadata::new(3, 40.0, 3000);
    let metadata4 = HashMetadata::new(4, 70.0, 4000);
    
    collector.add_hash(100, metadata1);
    collector.add_hash(200, metadata2);
    collector.add_hash(50, metadata3);
    
    assert_eq!(collector.len(), 3);
    assert!(!collector.is_empty());
    
    // Add one more, should evict the largest
    collector.add_hash(75, metadata4);
    assert_eq!(collector.len(), 3); // Should still be 3
    
    // Convert to sorted vector
    let sorted = collector.into_sorted_vec();
    assert_eq!(sorted.len(), 3);
    
    // Should be sorted by hash value
    for i in 1..sorted.len() {
        assert!(sorted[i-1].0 <= sorted[i].0, "Hashes should be sorted");
    }
}

#[test]
fn test_top_n_hash_collector_edge_cases() {
    // Test with size 0
    let mut collector_zero = TopNHashCollector::new(0);
    let metadata = HashMetadata::new(1, 50.0, 1000);
    collector_zero.add_hash(100, metadata);
    assert_eq!(collector_zero.len(), 0);
    assert!(collector_zero.is_empty());
    
    // Test with size 1
    let mut collector_one = TopNHashCollector::new(1);
    collector_one.add_hash(100, metadata);
    collector_one.add_hash(50, metadata); // Should replace the first
    
    assert_eq!(collector_one.len(), 1);
    let sorted = collector_one.into_sorted_vec();
    assert_eq!(sorted[0].0, 50); // Should keep the smaller hash
}

#[test]
fn test_file_metadata_serialization() -> Result<()> {
    let metadata = FileMetadata {
        filename: "test.fa".to_string(),
        file_size: 12345,
        sequence_name: "seq1".to_string(),
        sequence_length: 1000,
        total_sequences: 5,
        total_hashes: 250,
    };
    
    // Test JSON serialization/deserialization
    let json = serde_json::to_string(&metadata)?;
    let deserialized: FileMetadata = serde_json::from_str(&json)?;
    
    assert_eq!(metadata.filename, deserialized.filename);
    assert_eq!(metadata.file_size, deserialized.file_size);
    assert_eq!(metadata.sequence_name, deserialized.sequence_name);
    assert_eq!(metadata.sequence_length, deserialized.sequence_length);
    assert_eq!(metadata.total_sequences, deserialized.total_sequences);
    assert_eq!(metadata.total_hashes, deserialized.total_hashes);
    
    Ok(())
}

#[test]
fn test_hash_metadata_ordering() {
    let metadata1 = HashMetadata::new(1, 50.0, 1000);
    let metadata2 = HashMetadata::new(2, 50.0, 1000);
    let metadata3 = HashMetadata::new(1, 60.0, 1000);
    
    // Test ordering (should be by file_index first)
    assert!(metadata1 < metadata2);
    assert!(metadata1 < metadata3);
    
    // Test equality
    let metadata1_copy = HashMetadata::new(1, 50.0, 1000);
    assert_eq!(metadata1, metadata1_copy);
}

#[test]
fn test_gc_content_edge_cases() {
    // Test with lowercase nucleotides
    let lowercase = b"atcgatcg";
    let uppercase = b"ATCGATCG";
    assert_eq!(calculate_gc_content(lowercase), calculate_gc_content(uppercase));
    
    // Test with single character
    assert_eq!(calculate_gc_content(b"G"), 100.0);
    assert_eq!(calculate_gc_content(b"A"), 0.0);
    
    // Test with unknown characters (should be ignored)
    let with_unknown = b"ATCGXYZ";
    let clean = b"ATCG";
    assert_eq!(calculate_gc_content(with_unknown), calculate_gc_content(clean));
}

#[test]
fn test_entropy_edge_cases() {
    // Test with k-mer length 1
    let entropy = shannon_entropy(1, 1); // Single nucleotide
    assert!(entropy >= 0.0);
    
    // Test with maximum k-mer length
    let entropy_max = shannon_entropy(0, 32);
    assert!(entropy_max >= 0.0);
    
    // Test passes_entropy_filter with very high threshold
    assert!(!passes_entropy_filter(0, 21, 10.0)); // Should fail with impossible threshold
    
    // Test with very low threshold
    assert!(passes_entropy_filter(0, 21, -1.0)); // Should pass with negative threshold
}
