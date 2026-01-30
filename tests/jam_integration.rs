use jam_rs::format::{BUCKET_COUNT, bucket_id};
use jam_rs::query::QueryEngine;
use jam_rs::reader::JamReader;
use jam_rs::writer::{BuildConfig, build};
use std::io::Write;
use tempfile::NamedTempFile;

fn make_fasta(seqs: &[(&str, &str)]) -> NamedTempFile {
    let mut f = NamedTempFile::with_suffix(".fa").unwrap();
    for (name, seq) in seqs {
        writeln!(f, ">{name}").unwrap();
        writeln!(f, "{seq}").unwrap();
    }
    f
}

fn extract_hashes_from_db(reader: &JamReader) -> Vec<u64> {
    let mut hashes = Vec::new();
    for bucket_idx in 0..BUCKET_COUNT {
        for entry in reader.bucket_entries(bucket_idx) {
            hashes.push(entry.hash);
        }
    }
    hashes.sort_unstable();
    hashes.dedup();
    hashes
}

#[test]
fn test_full_roundtrip() {
    let seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    let input = make_fasta(&[("seq1", seq)]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1,
        num_threads: 2,
        memory: 1,
        ..Default::default()
    };

    let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();
    assert!(stats.total_entries > 0);
    assert_eq!(stats.sample_count, 1);
    assert_eq!(stats.kmer_size, 11);

    // Query with hashes from the database itself
    let reader = JamReader::open(&output_path).unwrap();
    let query_hashes = extract_hashes_from_db(&reader);
    assert!(!query_hashes.is_empty(), "Database should have hashes");

    let engine = QueryEngine::open(&output_path).unwrap();
    let result = engine.query(&query_hashes);
    assert!(result.has_matches());
    assert!(result.hashes_found > 0);

    // Should have 100% containment since we're querying with the exact same hashes
    let top = result.top(1);
    assert!(!top.is_empty());
    assert!(
        top[0].containment >= 1.0,
        "Expected 100% containment, got {}",
        top[0].containment
    );
}

#[test]
fn test_singleton_mode() {
    let input = make_fasta(&[
        ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
    ]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1,
        singleton: true,
        num_threads: 1,
        memory: 1,
        ..Default::default()
    };

    let stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();
    assert_eq!(
        stats.sample_count, 2,
        "Each sequence should be a separate sample"
    );

    let reader = JamReader::open(&output_path).unwrap();
    let reader_stats = reader.stats();
    assert_eq!(reader_stats.sample_count, 2);
}

#[test]
fn test_multiple_samples_shared_hashes() {
    // Two identical sequences should share all hashes
    let seq = "ATCGATCGATCGATCGATCGATCGATCGATCG";
    let input = make_fasta(&[("seq1", seq), ("seq2", seq)]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1,
        singleton: true,
        num_threads: 1,
        memory: 1,
        ..Default::default()
    };

    build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

    let reader = JamReader::open(&output_path).unwrap();
    let query_hashes = extract_hashes_from_db(&reader);

    let engine = QueryEngine::open(&output_path).unwrap();
    let result = engine.query(&query_hashes);
    // Both samples should have matches
    assert!(result.matches.len() >= 2 || result.matches.iter().any(|m| m.hit_count > 0));
}

#[test]
fn test_empty_buckets() {
    // Use very restrictive fscale to ensure most buckets are empty
    let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1_000_000,
        num_threads: 1,
        memory: 1,
        ..Default::default()
    };

    let result = build(&[input.path().to_path_buf()], &output_path, &config);
    assert!(result.is_ok());

    let reader = JamReader::open(&output_path).unwrap();
    let stats = reader.stats();

    // Most buckets should be empty
    let empty_count = stats
        .bucket_entry_counts
        .iter()
        .filter(|&&c| c == 0)
        .count();
    assert!(
        empty_count > 200,
        "Most buckets should be empty with high fscale"
    );
}

#[test]
fn test_reader_stats() {
    let input = make_fasta(&[
        ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
    ]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 21,
        fscale: 10,
        num_threads: 2,
        memory: 1,
        ..Default::default()
    };

    let build_stats = build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

    let reader = JamReader::open(&output_path).unwrap();
    let reader_stats = reader.stats();

    assert_eq!(reader_stats.entry_count, build_stats.total_entries);
    assert_eq!(reader_stats.unique_hash_count, build_stats.unique_hashes);
    assert_eq!(reader_stats.sample_count, build_stats.sample_count);
    assert_eq!(reader_stats.kmer_size, 21);
    assert!(reader_stats.file_size > 0);
}

#[test]
fn test_bucket_distribution() {
    let input = make_fasta(&[
        ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"),
        ("seq2", "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"),
        ("seq3", "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG"),
    ]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1,
        num_threads: 2,
        memory: 1,
        ..Default::default()
    };

    build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

    let reader = JamReader::open(&output_path).unwrap();

    // Verify entries are correctly bucketed
    for bucket_idx in 0..BUCKET_COUNT {
        let entries = reader.bucket_entries(bucket_idx);
        for entry in entries {
            assert_eq!(bucket_id(entry.hash), bucket_idx, "Entry in wrong bucket");
        }

        // Verify entries are sorted within bucket
        for window in entries.windows(2) {
            assert!(window[0] <= window[1], "Entries not sorted");
        }
    }
}

#[test]
fn test_query_batch() {
    let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")]);
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("test.jam");

    let config = BuildConfig {
        kmer_size: 11,
        fscale: 1,
        num_threads: 1,
        memory: 1,
        ..Default::default()
    };

    build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

    let engine = QueryEngine::open(&output_path).unwrap();

    // Create multiple queries
    let reader = JamReader::open(&output_path).unwrap();
    let mut all_hashes = Vec::new();
    for bucket_idx in 0..BUCKET_COUNT {
        for entry in reader.bucket_entries(bucket_idx).iter().take(3) {
            all_hashes.push(entry.hash);
        }
    }

    // Split into batches
    let batch1: Vec<u64> = all_hashes.iter().step_by(2).copied().collect();
    let batch2: Vec<u64> = all_hashes.iter().skip(1).step_by(2).copied().collect();

    let results = engine.query_batch(&[batch1, batch2]);
    assert_eq!(results.len(), 2);
}
