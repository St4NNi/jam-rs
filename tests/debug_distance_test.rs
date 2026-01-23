use anyhow::Result;
use jam_rs::distance::{
    DistanceConfig, LengthCategoryMode, OutputFormat, calculate_distances_streaming,
};
use jam_rs::sketch::{SketchConfig, sketch_files};
use std::fs;
use tempfile::TempDir;

#[test]
fn test_debug_distance_simple() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let fasta_content = r#">seq1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>seq2
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"#;

    let fasta_path = temp_dir.path().join("test.fa");
    fs::write(&fasta_path, fasta_content)?;

    let sketch_path = temp_dir.path().join("test.lmdb");
    let config = SketchConfig {
        kmer_size: 15,
        fscale: 100,
        nmax: 1000,
        singleton: false,
        min_entropy: 0.0,
        ..SketchConfig::default()
    };

    sketch_files(&[fasta_path], sketch_path.clone(), config, false)?;

    let distance_config = DistanceConfig {
        cutoff: 0.0,
        output_format: OutputFormat::Json,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };

    let _results =
        calculate_distances_streaming(&sketch_path, &sketch_path, None, distance_config, false, None)?;

    // Test passed if we reach this point without panicking

    Ok(())
}
