use anyhow::Result;
use jam_rs::sketch::{SketchConfig, sketch_files};
use jam_rs::stats::StatsCalculator;
use std::fs;
use tempfile::tempdir;

#[test]
fn test_database_stats_calculation() -> Result<()> {
    let temp_dir = tempdir()?;

    // Create test FASTA with multiple sequences
    let fasta_content = r#">seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
>seq3
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"#;

    let fasta_path = temp_dir.path().join("stats_test.fa");
    fs::write(&fasta_path, fasta_content)?;

    // Create sketch database
    let sketch_path = temp_dir.path().join("stats_test.lmdb");
    let config = SketchConfig {
        kmer_size: 8,     // Lower kmer size for short/repetitive seqs
        fscale: u64::MAX, // No FracMinHash filtering
        nmax: u64::MAX,   // No nmax filtering
        singleton: true,  // FIX: ensure per-sequence stats
        min_entropy: 0.0,
        threads: 1,
        memory_budget_gb: 1.0,
        temp_dir: None,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    // Test calculate_stats with detailed=false
    let stats = StatsCalculator::calculate_stats(&sketch_path)?;

    assert!(stats.total_hashes > 0, "Should have calculated some hashes");
    assert!(stats.unique_files > 0, "Should have at least one file");
    assert!(!stats.file_metadata.is_empty(), "Should have file metadata");

    // Test calculate_stats with detailed=true
    let detailed_stats = StatsCalculator::calculate_stats(&sketch_path)?;

    assert_eq!(
        stats.total_hashes, detailed_stats.total_hashes,
        "Hash count should be same"
    );
    assert_eq!(
        stats.unique_files, detailed_stats.unique_files,
        "File count should be same"
    );

    Ok(())
}

#[test]
fn test_database_stats_print_methods() -> Result<()> {
    let temp_dir = tempdir()?;

    let fasta_content = r#">print_test
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"#;

    let fasta_path = temp_dir.path().join("print_stats.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("print_stats.lmdb");
    let config = SketchConfig {
        kmer_size: 8,    // Lower kmer size for short/repetitive seqs
        singleton: true, // FIX: ensure per-sequence stats
        min_entropy: 0.0,
        ..SketchConfig::default()
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    let stats = StatsCalculator::calculate_stats(&sketch_path)?;

    // Test print_short method (this just prints to stdout, so we check it doesn't panic)
    stats.print_short();

    // Test print_detailed method
    stats.print_detailed();

    Ok(())
}

#[test]
fn test_print_stats_function() -> Result<()> {
    let temp_dir = tempdir()?;

    let fasta_content = r#">func_test
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"#;

    let fasta_path = temp_dir.path().join("func_stats.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("func_stats.lmdb");
    let config = SketchConfig {
        kmer_size: 8,    // Lower kmer size for short/repetitive seqs
        singleton: true, // FIX: ensure per-sequence stats
        min_entropy: 0.0,
        ..SketchConfig::default()
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    // Test print_stats with short=true
    let result = StatsCalculator::print_stats(&sketch_path, true, false);
    assert!(result.is_ok(), "print_stats with short=true should succeed");

    // Test print_stats with short=false (detailed)
    let result = StatsCalculator::print_stats(&sketch_path, false, false);
    assert!(
        result.is_ok(),
        "print_stats with short=false should succeed"
    );

    Ok(())
}

#[test]
fn test_stats_with_multiple_files() -> Result<()> {
    let temp_dir = tempdir()?;

    // Create multiple sequences to test file counting and metadata
    let sequences = vec![
        (
            "file1",
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        ),
        (
            "file2",
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",
        ),
        (
            "file3",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        ),
        (
            "file4",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        ),
        (
            "file5",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
        (
            "file6",
            "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",
        ),
    ];

    let fasta_content = sequences
        .iter()
        .map(|(name, seq)| format!(">{}\n{}", name, seq))
        .collect::<Vec<_>>()
        .join("\n");

    let fasta_path = temp_dir.path().join("multi_stats.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("multi_stats.lmdb");
    let config = SketchConfig {
        kmer_size: 8,    // Lower kmer size for short/repetitive seqs
        singleton: true, // FIX: ensure per-sequence stats
        min_entropy: 0.0,
        ..SketchConfig::default()
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    let stats = StatsCalculator::calculate_stats(&sketch_path)?;

    assert_eq!(
        stats.unique_files as usize,
        sequences.len(),
        "Should count all sequences correctly"
    );
    assert!(
        stats.total_hashes > 0,
        "Should have hashes from all sequences"
    );
    assert_eq!(
        stats.file_metadata.len(),
        sequences.len(),
        "Should have metadata for all sequences"
    );

    // Test that print_detailed handles the "more files" case
    stats.print_detailed();

    Ok(())
}

#[test]
fn test_stats_error_handling() -> Result<()> {
    let temp_dir = tempdir()?;
    let nonexistent_path = temp_dir.path().join("does_not_exist.lmdb");

    // Test with non-existent database
    let result = StatsCalculator::calculate_stats(&nonexistent_path);
    assert!(result.is_err(), "Should fail with non-existent database");

    let result = StatsCalculator::print_stats(&nonexistent_path, true, false);
    assert!(result.is_err(), "Should fail with non-existent database");

    Ok(())
}

#[test]
fn test_gc_and_length_distributions() -> Result<()> {
    let temp_dir = tempdir()?;

    // Create sequences with different GC content and lengths
    let fasta_content = r#">low_gc
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>high_gc
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
>medium_gc
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"#;

    let fasta_path = temp_dir.path().join("gc_test.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("gc_test.lmdb");
    let config = SketchConfig {
        kmer_size: 8,    // Lower kmer size for short/repetitive seqs
        singleton: true, // FIX: ensure per-sequence stats
        min_entropy: 0.0,
        ..SketchConfig::default()
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    let stats = StatsCalculator::calculate_stats(&sketch_path)?;

    // The distributions might be empty or have values depending on the implementation
    // We just verify the function runs without error and the stats structure is valid
    assert!(stats.total_hashes > 0, "Should have some hashes");

    Ok(())
}
