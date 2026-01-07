use anyhow::Result;
use std::path::{Path, PathBuf};

pub fn expand_input_paths(input_paths: &[PathBuf]) -> Result<Vec<PathBuf>> {
    let mut expanded_paths = Vec::new();

    for path in input_paths {
        if path.is_dir() {
            for entry in std::fs::read_dir(path)? {
                let entry = entry?;
                let file_path = entry.path();

                if is_sequence_file(&file_path) {
                    expanded_paths.push(file_path);
                }
            }
        } else if path.is_file() {
            if is_sequence_file(path) {
                expanded_paths.push(path.clone());
            } else {
                let content = std::fs::read_to_string(path)?;
                for line in content.lines() {
                    let file_path = PathBuf::from(line.trim());
                    if file_path.exists() && is_sequence_file(&file_path) {
                        expanded_paths.push(file_path);
                    }
                }
            }
        }
    }

    if expanded_paths.is_empty() {
        return Err(anyhow::anyhow!(
            "No valid sequence files found in input paths"
        ));
    }

    expanded_paths.sort();
    Ok(expanded_paths)
}

pub fn is_sequence_file(path: &Path) -> bool {
    if let Some(ext) = path.extension().map(|e| e.to_string_lossy().to_lowercase()) {
        if ext == "gz"
            && let Some(stem_ext) = path.file_stem().and_then(|s| Path::new(s).extension())
        {
            let stem_ext = stem_ext.to_string_lossy().to_lowercase();
            return matches!(
                stem_ext.as_str(),
                "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
            );
        }
        return matches!(
            ext.as_str(),
            "fasta" | "fa" | "fas" | "fna" | "fastq" | "fq"
        );
    }
    false
}
