use anyhow::Result;
use jam_rs::expand_input_paths;
use std::fs;
use std::path::PathBuf;
use tempfile::TempDir;

/// Unit tests for lib.rs internal functions to maximize code coverage
struct LibTestEnvironment {
    temp_dir: TempDir,
}

impl LibTestEnvironment {
    fn new() -> Result<Self> {
        let temp_dir = tempfile::tempdir()?;
        Ok(LibTestEnvironment { temp_dir })
    }

    fn temp_path(&self) -> &std::path::Path {
        self.temp_dir.path()
    }

    fn create_test_fasta(&self, name: &str, content: &str) -> Result<PathBuf> {
        let path = self.temp_path().join(name);
        fs::write(&path, content)?;
        Ok(path)
    }

    fn create_directory_with_files(
        &self,
        dir_name: &str,
        files: &[(&str, &str)],
    ) -> Result<PathBuf> {
        let dir_path = self.temp_path().join(dir_name);
        fs::create_dir_all(&dir_path)?;

        for (filename, content) in files {
            let file_path = dir_path.join(filename);
            fs::write(file_path, content)?;
        }

        Ok(dir_path)
    }

    fn _get_test_file(&self, filename: &str) -> PathBuf {
        PathBuf::from("tests/testfiles").join(filename)
    }
}

#[test]
fn test_expand_input_paths_single_file() -> Result<()> {
    let env = LibTestEnvironment::new()?;
    let test_file = env.create_test_fasta("test.fa", ">seq\nACGT\n")?;

    let result = expand_input_paths(std::slice::from_ref(&test_file))?;

    assert_eq!(result.len(), 1);
    assert_eq!(result[0], test_file);

    Ok(())
}

#[test]
fn test_expand_input_paths_directory() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let test_files = vec![
        ("seq1.fa", ">seq1\nACGT\n"),
        ("seq2.fasta", ">seq2\nTGCA\n"),
        ("seq3.fastq", "@seq3\nGGCC\n+\nIIII\n"),
        ("readme.txt", "Not a sequence file"),
    ];

    let test_dir = env.create_directory_with_files("test_dir", &test_files)?;
    let result = expand_input_paths(&[test_dir])?;

    // Should find 3 sequence files, sorted alphabetically
    assert_eq!(result.len(), 3);
    assert!(
        result[0]
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("seq1")
    );
    assert!(
        result[1]
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("seq2")
    );
    assert!(
        result[2]
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .starts_with("seq3")
    );

    Ok(())
}

#[test]
fn test_expand_input_paths_file_list() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let file1 = env.create_test_fasta("file1.fa", ">seq1\nACGT\n")?;
    let file2 = env.create_test_fasta("file2.fa", ">seq2\nTGCA\n")?;

    let file_list_content = format!("{}\n{}\n", file1.to_str().unwrap(), file2.to_str().unwrap());
    let file_list = env.temp_path().join("files.txt");
    fs::write(&file_list, file_list_content)?;

    let result = expand_input_paths(&[file_list])?;

    assert_eq!(result.len(), 2);
    assert!(result.contains(&file1));
    assert!(result.contains(&file2));

    Ok(())
}

#[test]
fn test_expand_input_paths_mixed_extensions() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let test_files = vec![
        ("test.fa", ">seq\nACGT\n"),
        ("test.fasta", ">seq\nACGT\n"),
        ("test.fas", ">seq\nACGT\n"),
        ("test.fna", ">seq\nACGT\n"),
        ("test.fastq", "@seq\nACGT\n+\nIIII\n"),
        ("test.fq", "@seq\nACGT\n+\nIIII\n"),
        ("test.fa.gz", ">seq\nACGT\n"), // Would be compressed in real usage
        ("test.fastq.gz", "@seq\nACGT\n+\nIIII\n"), // Would be compressed in real usage
    ];

    let test_dir = env.create_directory_with_files("extensions", &test_files)?;
    let result = expand_input_paths(&[test_dir])?;

    assert_eq!(result.len(), 8);

    Ok(())
}

#[test]
fn test_expand_input_paths_empty_directory() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let empty_dir = env.temp_path().join("empty");
    fs::create_dir(&empty_dir)?;

    let result = expand_input_paths(&[empty_dir]);

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
fn test_expand_input_paths_nonexistent_file() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let nonexistent = env.temp_path().join("nonexistent.fa");
    let result = expand_input_paths(&[nonexistent]);

    assert!(result.is_err());

    Ok(())
}

#[test]
fn test_expand_input_paths_multiple_inputs() -> Result<()> {
    let env = LibTestEnvironment::new()?;

    let file1 = env.create_test_fasta("single.fa", ">seq1\nACGT\n")?;
    let test_dir = env.create_directory_with_files("dir", &[("dir_file.fa", ">seq2\nTGCA\n")])?;

    let result = expand_input_paths(&[file1.clone(), test_dir])?;

    assert_eq!(result.len(), 2);
    // Results should be sorted
    assert!(result.contains(&file1));

    Ok(())
}

// Now let's test the handle functions directly by importing them from lib.rs
// We need to make these functions public in lib.rs first, or create a test module

#[cfg(test)]
mod handle_function_tests {
    use super::*;
    use jam_rs::expand_input_paths;

    // Test expand_input_paths edge cases
    #[test]
    fn test_expand_paths_file_list_with_invalid_entries() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        let valid_file = env.create_test_fasta("valid.fa", ">seq\nACGT\n")?;
        let invalid_path = env.temp_path().join("nonexistent.fa");

        let file_list_content = format!(
            "{}\n{}\n# comment line\n\n",
            valid_file.to_str().unwrap(),
            invalid_path.to_str().unwrap()
        );
        let file_list = env.temp_path().join("mixed_list.txt");
        fs::write(&file_list, file_list_content)?;

        let result = expand_input_paths(&[file_list])?;

        // Should only include valid files
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], valid_file);

        Ok(())
    }

    #[test]
    fn test_expand_paths_nested_directories() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        // Create nested structure
        let level1 = env.temp_path().join("level1");
        let level2 = level1.join("level2");
        fs::create_dir_all(&level2)?;

        // Only files in the direct directory should be found (not recursive)
        fs::write(level1.join("file1.fa"), ">seq1\nACGT\n")?;
        fs::write(level2.join("file2.fa"), ">seq2\nTGCA\n")?;

        let result = expand_input_paths(&[level1])?;

        // Should only find file1.fa, not file2.fa (no recursive search)
        assert_eq!(result.len(), 1);
        assert!(result[0].file_name().unwrap().to_str().unwrap() == "file1.fa");

        Ok(())
    }

    #[test]
    fn test_expand_paths_case_insensitive_extensions() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        let test_files = vec![
            ("test.FA", ">seq\nACGT\n"),
            ("test.FASTA", ">seq\nACGT\n"),
            ("test.FQ", "@seq\nACGT\n+\nIIII\n"),
        ];

        let test_dir = env.create_directory_with_files("case_test", &test_files)?;
        let result = expand_input_paths(&[test_dir])?;

        assert_eq!(result.len(), 3);

        Ok(())
    }

    #[test]
    fn test_expand_paths_sorting() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        let test_files = vec![
            ("z_last.fa", ">seq\nACGT\n"),
            ("a_first.fa", ">seq\nACGT\n"),
            ("m_middle.fa", ">seq\nACGT\n"),
        ];

        let test_dir = env.create_directory_with_files("sort_test", &test_files)?;
        let result = expand_input_paths(&[test_dir])?;

        assert_eq!(result.len(), 3);
        // Should be sorted alphabetically
        assert!(
            result[0]
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("a_first")
        );
        assert!(
            result[1]
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("m_middle")
        );
        assert!(
            result[2]
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .starts_with("z_last")
        );

        Ok(())
    }

    #[test]
    fn test_expand_paths_whitespace_in_file_list() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        let file1 = env.create_test_fasta("file1.fa", ">seq1\nACGT\n")?;
        let file2 = env.create_test_fasta("file2.fa", ">seq2\nTGCA\n")?;

        let file_list_content = format!(
            "  {}\t\n\n  {}  \n",
            file1.to_str().unwrap(),
            file2.to_str().unwrap()
        );
        let file_list = env.temp_path().join("whitespace_list.txt");
        fs::write(&file_list, file_list_content)?;

        let result = expand_input_paths(&[file_list])?;

        assert_eq!(result.len(), 2);
        assert!(result.contains(&file1));
        assert!(result.contains(&file2));

        Ok(())
    }
}

// Test sequence file detection function
#[cfg(test)]
mod sequence_file_tests {
    use super::*;

    // We need to access the is_sequence_file function - it should be made public or tested indirectly

    #[test]
    fn test_sequence_file_detection_through_expand() -> Result<()> {
        let env = LibTestEnvironment::new()?;

        let test_files = vec![
            ("valid.fa", ">seq\nACGT\n", true),
            ("valid.fasta", ">seq\nACGT\n", true),
            ("valid.fastq", "@seq\nACGT\n+\nIIII\n", true),
            ("valid.fq", "@seq\nACGT\n+\nIIII\n", true),
            ("valid.fa.gz", ">seq\nACGT\n", true),
            ("invalid.txt", "not a sequence", false),
            ("invalid.csv", "col1,col2\n1,2\n", false),
            ("no_extension", ">seq\nACGT\n", false),
        ];

        let mut valid_count = 0;
        let test_dir = env.temp_path().join("file_detection");
        fs::create_dir(&test_dir)?;

        for (filename, content, should_be_valid) in &test_files {
            fs::write(test_dir.join(filename), content)?;
            if *should_be_valid {
                valid_count += 1;
            }
        }

        let result = expand_input_paths(&[test_dir])?;
        assert_eq!(result.len(), valid_count);

        Ok(())
    }
}
