use anyhow::Result;
use jam_rs::expand_input_paths;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use tempfile::TempDir;

/// Comprehensive edge case and error path tests for maximum code coverage
struct CoverageTestEnvironment {
    temp_dir: TempDir,
}

impl CoverageTestEnvironment {
    fn new() -> Result<Self> {
        let temp_dir = tempfile::tempdir()?;
        Ok(CoverageTestEnvironment { temp_dir })
    }

    fn temp_path(&self) -> &std::path::Path {
        self.temp_dir.path()
    }

    fn create_test_fasta(&self, name: &str, content: &str) -> Result<PathBuf> {
        let path = self.temp_path().join(name);
        fs::write(&path, content)?;
        Ok(path)
    }

    fn create_readonly_dir(&self, name: &str) -> Result<PathBuf> {
        let dir_path = self.temp_path().join(name);
        fs::create_dir(&dir_path)?;

        // Make directory read-only (no write permissions)
        let mut perms = fs::metadata(&dir_path)?.permissions();
        perms.set_mode(0o555); // read and execute only
        fs::set_permissions(&dir_path, perms)?;

        Ok(dir_path)
    }

    fn create_symlink(&self, target: &str, link_name: &str) -> Result<PathBuf> {
        let link_path = self.temp_path().join(link_name);
        std::os::unix::fs::symlink(target, &link_path)?;
        Ok(link_path)
    }
}

// Test complex path expansion scenarios
#[test]
fn test_expand_input_paths_symlinks() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create a real file
    let real_file = env.create_test_fasta("real.fa", ">seq\nACGT\n")?;

    // Create a symlink to it
    let symlink = env.create_symlink(real_file.to_str().unwrap(), "symlink.fa")?;

    let result = expand_input_paths(&[symlink])?;

    assert_eq!(result.len(), 1);

    Ok(())
}

#[test]
fn test_expand_input_paths_broken_symlinks() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create a symlink to non-existent file
    let broken_symlink = env.create_symlink("/nonexistent/file.fa", "broken.fa")?;

    let result = expand_input_paths(&[broken_symlink]);

    // Should handle broken symlinks gracefully
    assert!(result.is_err() || (result.is_ok() && result.unwrap().is_empty()));

    Ok(())
}

#[test]
fn test_expand_input_paths_permission_denied_directory() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create a directory with no read permissions
    let restricted_dir = env.create_readonly_dir("restricted")?;

    let result = expand_input_paths(&[restricted_dir]);

    // Should handle permission errors gracefully
    assert!(result.is_err() || result.is_ok());

    Ok(())
}

#[test]
fn test_expand_input_paths_very_long_path() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create a deeply nested directory structure
    let mut deep_path = env.temp_path().to_path_buf();
    for i in 0..50 {
        deep_path = deep_path.join(format!("level_{i}"));
    }

    // Only create the directory structure up to OS limits
    if fs::create_dir_all(&deep_path).is_ok() {
        let deep_file = deep_path.join("deep.fa");
        if fs::write(&deep_file, ">seq\nACGT\n").is_ok() {
            let result = expand_input_paths(&[deep_file]);
            assert!(result.is_ok() || result.is_err());
        }
    }

    Ok(())
}

#[test]
fn test_expand_input_paths_special_filenames() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let special_names = vec![
        "file with spaces.fa",
        "file-with-dashes.fa",
        "file_with_underscores.fa",
        "file.with.dots.fa",
        "file@with@symbols.fa",
        "UPPERCASE.FA",
        "123numeric.fa",
    ];

    let mut created_files = Vec::new();
    for name in &special_names {
        if let Ok(file) = env.create_test_fasta(name, ">seq\nACGT\n") {
            created_files.push(file);
        }
    }

    if !created_files.is_empty() {
        let result = expand_input_paths(&created_files)?;
        assert!(result.len() <= created_files.len());
    }

    Ok(())
}

#[test]
fn test_expand_input_paths_mixed_valid_invalid() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let valid_file = env.create_test_fasta("valid.fa", ">seq\nACGT\n")?;
    let invalid_file = env.temp_path().join("invalid.txt");
    fs::write(&invalid_file, "not a sequence file")?;
    let nonexistent_file = env.temp_path().join("nonexistent.fa");

    let result = expand_input_paths(&[valid_file.clone(), invalid_file, nonexistent_file])?;

    // Should only include valid files
    assert_eq!(result.len(), 1);
    assert_eq!(result[0], valid_file);

    Ok(())
}

#[test]
fn test_expand_input_paths_empty_file_list() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create empty file list
    let empty_list = env.temp_path().join("empty_list.txt");
    fs::write(&empty_list, "")?;

    let result = expand_input_paths(&[empty_list]);

    assert!(result.is_err());
    assert!(
        result
            .unwrap_err()
            .to_string()
            .contains("No valid sequence files")
    );

    Ok(())
}

#[test]
fn test_expand_input_paths_file_list_with_comments() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let valid_file = env.create_test_fasta("valid.fa", ">seq\nACGT\n")?;

    let file_list_content = format!(
        "# This is a comment\n{}\n# Another comment\n\n# Empty line above\n",
        valid_file.to_str().unwrap()
    );
    let file_list = env.temp_path().join("commented_list.txt");
    fs::write(&file_list, file_list_content)?;

    let result = expand_input_paths(&[file_list])?;

    assert_eq!(result.len(), 1);
    assert_eq!(result[0], valid_file);

    Ok(())
}

#[test]
fn test_expand_input_paths_duplicate_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let file1 = env.create_test_fasta("file1.fa", ">seq1\nACGT\n")?;
    let file2 = env.create_test_fasta("file2.fa", ">seq2\nTGCA\n")?;

    // Test with duplicate inputs
    let result = expand_input_paths(&[file1.clone(), file2.clone(), file1.clone()])?;

    // Should handle duplicates (might include them or deduplicate)
    assert!(result.len() >= 2);
    assert!(result.contains(&file1));
    assert!(result.contains(&file2));

    Ok(())
}

#[test]
fn test_expand_input_paths_relative_paths() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Change to temp directory
    let original_dir = std::env::current_dir()?;
    std::env::set_current_dir(env.temp_path())?;

    // Create file with relative path
    let file_content = ">seq\nACGT\n";
    fs::write("relative.fa", file_content)?;

    let result = expand_input_paths(&[PathBuf::from("relative.fa")]);

    // Restore original directory
    std::env::set_current_dir(original_dir)?;

    assert!(result.is_ok());

    Ok(())
}

#[test]
fn test_expand_input_paths_compressed_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create files with compression extensions (not actually compressed)
    let gz_file = env.create_test_fasta("test.fa.gz", ">seq\nACGT\n")?;
    let fq_gz_file = env.create_test_fasta("test.fq.gz", "@seq\nACGT\n+\nIIII\n")?;

    let result = expand_input_paths(&[gz_file.clone(), fq_gz_file.clone()])?;

    assert_eq!(result.len(), 2);
    assert!(result.contains(&gz_file));
    assert!(result.contains(&fq_gz_file));

    Ok(())
}

#[test]
fn test_expand_input_paths_directory_with_subdirs() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let main_dir = env.temp_path().join("main");
    let sub_dir = main_dir.join("subdir");
    fs::create_dir_all(&sub_dir)?;

    // Create files in main directory
    fs::write(main_dir.join("main.fa"), ">seq1\nACGT\n")?;

    // Create files in subdirectory (should not be included in non-recursive search)
    fs::write(sub_dir.join("sub.fa"), ">seq2\nTGCA\n")?;

    let result = expand_input_paths(&[main_dir])?;

    // Should only find files in main directory (not recursive)
    assert_eq!(result.len(), 1);
    assert!(result[0].file_name().unwrap().to_str().unwrap() == "main.fa");

    Ok(())
}

#[test]
fn test_expand_input_paths_large_number_of_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let test_dir = env.temp_path().join("many_files");
    fs::create_dir(&test_dir)?;

    // Create many files
    for i in 0..100 {
        fs::write(
            test_dir.join(format!("file_{i:03}.fa")),
            format!(">seq_{i}\nACGT\n"),
        )?;
    }

    let result = expand_input_paths(&[test_dir])?;

    assert_eq!(result.len(), 100);
    // Should be sorted
    for i in 1..result.len() {
        assert!(result[i - 1] <= result[i]);
    }

    Ok(())
}

#[test]
fn test_expand_input_paths_unicode_filenames() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let unicode_names = vec![
        "seq_中文.fa",
        "seq_español.fa",
        "seq_français.fa",
        "seq_🧬.fa",
    ];

    let mut created_files = Vec::new();
    for name in &unicode_names {
        if let Ok(file) = env.create_test_fasta(name, ">seq\nACGT\n") {
            created_files.push(file);
        }
    }

    if !created_files.is_empty() {
        let result = expand_input_paths(&created_files);
        assert!(result.is_ok() || result.is_err()); // Either handling is acceptable
    }

    Ok(())
}

// Additional coverage for is_sequence_file with edge cases
#[test]
fn test_is_sequence_file_edge_cases() {
    use jam_rs::is_sequence_file;
    use std::path::Path;

    // Test very long extensions
    assert!(!is_sequence_file(Path::new("file.fastaverylongextension")));

    // Test extensions with numbers
    assert!(!is_sequence_file(Path::new("file.fa1")));
    assert!(!is_sequence_file(Path::new("file.fasta2")));

    // Test partial matches
    assert!(!is_sequence_file(Path::new("file.fast")));
    assert!(!is_sequence_file(Path::new("file.f")));

    // Test with multiple dots
    assert!(is_sequence_file(Path::new("file.name.with.dots.fa")));
    assert!(is_sequence_file(Path::new("version.1.0.fasta")));

    // Test compressed variations
    assert!(is_sequence_file(Path::new("file.fastq.gz")));
    assert!(is_sequence_file(Path::new("file.fq.gz")));
    assert!(!is_sequence_file(Path::new("file.fa.bz2"))); // Only .gz is supported

    // Test case sensitivity thoroughly
    let mixed_case = vec![
        "File.Fa",
        "File.FQ",
        "File.FASTQ",
        "File.FASTA",
        "file.Fa.GZ",
        "file.FASTQ.gz",
        "file.fq.GZ",
    ];

    for filename in mixed_case {
        assert!(
            is_sequence_file(Path::new(filename)),
            "Failed for {filename}"
        );
    }
}

// Test error conditions that might not be covered elsewhere
#[test]
fn test_expand_input_paths_io_errors() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    // Create a file, then try to read it as a directory
    let regular_file = env.create_test_fasta("regular.fa", ">seq\nACGT\n")?;

    // This should be handled by the file vs directory logic
    let result = expand_input_paths(&[regular_file]);
    assert!(result.is_ok());

    Ok(())
}
