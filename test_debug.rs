use anyhow::Result;
use jam_rs::distance::{StreamingDistanceCalculator, DistanceConfig, LengthCategoryMode, OutputFormat};
use jam_rs::sketch::{sketch_files, SketchConfig};
use std::io::Write;
use tempfile::NamedTempFile;

fn main() -> Result<()> {
    println!("Creating input file...");
    let mut input_file = NamedTempFile::new()?;
    writeln!(input_file, ">test_seq")?;
    writeln!(input_file, "ATCGATCGATCGATCGATCGATCGATCGATCG")?;
    input_file.flush()?;
    let input_path = input_file.into_temp_path().to_path_buf();
    
    println!("Creating output path...");
    let output_file = NamedTempFile::new()?;
    let output_path = output_file.into_temp_path().to_path_buf();
    
    println!("Input: {:?}", input_path);
    println!("Output: {:?}", output_path);
    
    let config = SketchConfig {
        kmer_size: 21,
        min_entropy: 1.5,
        ..Default::default()
    };
    
    println!("Sketching files...");
    sketch_files(&[input_path], output_path.clone(), config)?;
    
    println!("Output file exists: {}", output_path.exists());
    
    println!("Attempting to open distance calculator...");
    let distance_config = DistanceConfig {
        cutoff: 0.0,
        output_format: OutputFormat::Tsv,
        length_category_mode: LengthCategoryMode::QueryAndBelow,
    };
    
    match StreamingDistanceCalculator::new(distance_config, false, &output_path) {
        Ok(_) => println!("SUCCESS: Distance calculator created!"),
        Err(e) => println!("ERROR: Failed to create distance calculator: {}", e),
    }
    
    Ok(())
}
