use anyhow::Result;
use jam_rs::distance::{
    DistanceConfig, LengthCategoryMode, OutputFormat, StreamingDistanceCalculator,
};
use jam_rs::sketch::{SketchConfig, sketch_files};
use std::fs;

/// Tests specifically targeting the write_json_results function and print_results function
/// to improve coverage for the JSON output functionality that was mentioned as missing

#[test]
fn test_json_output_direct_coverage() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;

    // Create a simple FASTA file with content that will generate sketches
    let fasta_content = r#">seq1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>seq2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"#;

    let fasta_path = temp_dir.path().join("test.fa");
    fs::write(&fasta_path, fasta_content)?;

    // Create sketch database
    let sketch_path = temp_dir.path().join("test.lmdb");
    let config = SketchConfig {
        kmer_size: 15, // Smaller k-mer to ensure we get some hashes
        fscale: 100,
        nmax: 1000,
        singleton: false,
        min_entropy: 0.0, // Very low entropy to accept any k-mer
        threads: 1,
        memory_budget_gb: 1.0,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config)?;

    // Test JSON output with very permissive settings
    let distance_config = DistanceConfig {
        cutoff: 0.0, // Accept all results
        output_format: OutputFormat::Json,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };

    // Create calculator and test self-comparison to ensure we get results
    let calculator = StreamingDistanceCalculator::new(distance_config, false, &sketch_path)?;
    let results = calculator.calculate_distances_streaming(&sketch_path, None)?;

    if !results.is_empty() {
        // Test JSON file output
        let json_output = temp_dir.path().join("results.json");
        let _results_with_file =
            calculator.calculate_distances_streaming(&sketch_path, Some(&json_output))?;

        // Verify JSON file was created and contains valid JSON
        assert!(json_output.exists(), "JSON output file should be created");
        let json_content = fs::read_to_string(&json_output)?;
        let _parsed: serde_json::Value = serde_json::from_str(&json_content)?;
    }

    Ok(())
}

#[test]
fn test_streaming_coverage_edge_cases() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;

    // Create test with very long sequences to test different code paths
    let long_seq = "A".repeat(1000);
    let fasta_content = format!(">long_seq\n{}", long_seq);

    let fasta_path = temp_dir.path().join("long.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("long.lmdb");
    let config = SketchConfig {
        kmer_size: 21,
        fscale: 1000,
        nmax: 1000,
        singleton: false,
        min_entropy: 0.0,
        threads: 1,
        memory_budget_gb: 1.0,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config)?;

    // Test with different length categories
    let configs = vec![
        DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::QueryAndBelow,
        },
        DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::Range(1),
        },
        DistanceConfig {
            cutoff: 0.0,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::Range(5),
        },
    ];

    for config in configs {
        let calculator = StreamingDistanceCalculator::new(config, false, &sketch_path)?;
        let _results = calculator.calculate_distances_streaming(&sketch_path, None)?;
    }

    Ok(())
}

#[test]
fn test_output_file_handling() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;

    // Create minimal test data
    let fasta_content =
        ">test\nATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    let fasta_path = temp_dir.path().join("minimal.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("minimal.lmdb");
    let config = SketchConfig {
        kmer_size: 15,
        fscale: 1000,
        nmax: 1000,
        singleton: false,
        min_entropy: 0.0,
        threads: 1,
        memory_budget_gb: 1.0,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config)?;

    // Test TSV output
    let tsv_config = DistanceConfig {
        cutoff: 0.0,
        output_format: OutputFormat::Tsv,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };

    let calculator = StreamingDistanceCalculator::new(tsv_config, false, &sketch_path)?;

    // Test with output file
    let tsv_output = temp_dir.path().join("results.tsv");
    let _results = calculator.calculate_distances_streaming(&sketch_path, Some(&tsv_output))?;

    // Drop the calculator to ensure LMDB environment is closed
    drop(calculator);

    // Test JSON output
    let json_config = DistanceConfig {
        cutoff: 0.0,
        output_format: OutputFormat::Json,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };

    let calculator = StreamingDistanceCalculator::new(json_config, false, &sketch_path)?;
    let json_output = temp_dir.path().join("results.json");
    let _results = calculator.calculate_distances_streaming(&sketch_path, Some(&json_output))?;
    Ok(())
}

#[test]
fn test_singleton_mode_coverage() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;

    let fasta_content = r#">singleton_test
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"#;

    let fasta_path = temp_dir.path().join("singleton.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("singleton.lmdb");
    let config = SketchConfig {
        kmer_size: 15,
        fscale: 1000,
        nmax: 1000,
        singleton: true, // Test singleton mode in sketching
        min_entropy: 0.0,
        threads: 1,
        memory_budget_gb: 1.0,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config)?;

    // Test with singleton mode in distance calculation
    let distance_config = DistanceConfig {
        cutoff: 0.0,
        output_format: OutputFormat::Tsv,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };

    let calculator = StreamingDistanceCalculator::new(distance_config, true, &sketch_path)?;
    let _results = calculator.calculate_distances_streaming(&sketch_path, None)?;

    Ok(())
}

#[test]
fn test_different_cutoff_values() -> Result<()> {
    let temp_dir = tempfile::tempdir()?;

    let fasta_content = r#">cutoff_test
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"#;

    let fasta_path = temp_dir.path().join("cutoff.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("cutoff.lmdb");
    let config = SketchConfig {
        kmer_size: 15,
        fscale: 1000,
        nmax: 1000,
        singleton: false,
        min_entropy: 0.0,
        threads: 1,
        memory_budget_gb: 1.0,
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config)?;

    // Test with different cutoff values to exercise filtering logic
    let cutoffs = vec![0.0, 0.1, 0.5, 0.9, 1.0];

    for cutoff in cutoffs {
        let distance_config = DistanceConfig {
            cutoff,
            output_format: OutputFormat::Tsv,
            length_category_mode: LengthCategoryMode::QueryAndBelow,
        };

        let calculator = StreamingDistanceCalculator::new(distance_config, false, &sketch_path)?;
        let _results = calculator.calculate_distances_streaming(&sketch_path, None)?;
    }

    Ok(())
}
