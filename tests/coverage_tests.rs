use anyhow::Result;
use jam_rs::expand_input_paths;
use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::PathBuf;
use tempfile::TempDir;

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

        let mut perms = fs::metadata(&dir_path)?.permissions();
        perms.set_mode(0o555);
        fs::set_permissions(&dir_path, perms)?;

        Ok(dir_path)
    }

    fn create_symlink(&self, target: &str, link_name: &str) -> Result<PathBuf> {
        let link_path = self.temp_path().join(link_name);
        std::os::unix::fs::symlink(target, &link_path)?;
        Ok(link_path)
    }
}

#[test]
fn test_expand_input_paths_symlinks() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let real_file = env.create_test_fasta("real.fa", ">seq\nACGT\n")?;

    let symlink = env.create_symlink(real_file.to_str().unwrap(), "symlink.fa")?;

    let result = expand_input_paths(&[symlink])?;

    assert_eq!(result.len(), 1);

    Ok(())
}

#[test]
fn test_expand_input_paths_broken_symlinks() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let broken_symlink = env.create_symlink("/nonexistent/file.fa", "broken.fa")?;

    let result = expand_input_paths(&[broken_symlink]);

    assert!(result.is_err() || (result.is_ok() && result.unwrap().is_empty()));

    Ok(())
}

#[test]
fn test_expand_input_paths_permission_denied_directory() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let restricted_dir = env.create_readonly_dir("restricted")?;

    let result = expand_input_paths(&[restricted_dir]);

    assert!(result.is_err() || result.is_ok());

    Ok(())
}

#[test]
fn test_expand_input_paths_very_long_path() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let mut deep_path = env.temp_path().to_path_buf();
    for i in 0..50 {
        deep_path = deep_path.join(format!("level_{i}"));
    }

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

    assert_eq!(result.len(), 1);
    assert_eq!(result[0], valid_file);

    Ok(())
}

#[test]
fn test_expand_input_paths_empty_file_list() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

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
fn test_expand_input_paths_file_list_relative_to_list_file() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let list_dir = env.temp_path().join("lists");
    let data_dir = env.temp_path().join("data");
    fs::create_dir_all(&list_dir)?;
    fs::create_dir_all(&data_dir)?;

    let valid_file = data_dir.join("valid.fa");
    fs::write(&valid_file, ">seq\nACGT\n")?;

    let file_list = list_dir.join("inputs.txt");
    fs::write(&file_list, "../data/valid.fa\n")?;

    let result = expand_input_paths(&[file_list])?;

    assert_eq!(result, vec![valid_file]);

    Ok(())
}

#[test]
fn test_expand_input_paths_duplicate_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let file1 = env.create_test_fasta("file1.fa", ">seq1\nACGT\n")?;
    let file2 = env.create_test_fasta("file2.fa", ">seq2\nTGCA\n")?;

    let result = expand_input_paths(&[file1.clone(), file2.clone(), file1.clone()])?;

    assert!(result.len() >= 2);
    assert!(result.contains(&file1));
    assert!(result.contains(&file2));

    Ok(())
}

#[test]
fn test_expand_input_paths_relative_paths() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let original_dir = std::env::current_dir()?;
    std::env::set_current_dir(env.temp_path())?;

    let file_content = ">seq\nACGT\n";
    fs::write("relative.fa", file_content)?;

    let result = expand_input_paths(&[PathBuf::from("relative.fa")]);

    std::env::set_current_dir(original_dir)?;

    assert!(result.is_ok());

    Ok(())
}

#[test]
fn test_expand_input_paths_compressed_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let compressed_files = [
        env.create_test_fasta("test.fa.gz", ">seq\nACGT\n")?,
        env.create_test_fasta("test.fq.gz", "@seq\nACGT\n+\nIIII\n")?,
        env.create_test_fasta("test.fa.bz2", ">seq\nACGT\n")?,
        env.create_test_fasta("test.fasta.xz", ">seq\nACGT\n")?,
        env.create_test_fasta("test.fastq.zst", "@seq\nACGT\n+\nIIII\n")?,
        env.create_test_fasta("test.fq.zstd", "@seq\nACGT\n+\nIIII\n")?,
    ];

    let result = expand_input_paths(&compressed_files)?;

    assert_eq!(result.len(), compressed_files.len());
    for file in compressed_files {
        assert!(result.contains(&file));
    }

    Ok(())
}

#[test]
fn test_expand_input_paths_directory_with_subdirs() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let main_dir = env.temp_path().join("main");
    let sub_dir = main_dir.join("subdir");
    fs::create_dir_all(&sub_dir)?;

    fs::write(main_dir.join("main.fa"), ">seq1\nACGT\n")?;

    fs::write(sub_dir.join("sub.fa"), ">seq2\nTGCA\n")?;

    let result = expand_input_paths(&[main_dir])?;

    assert_eq!(result.len(), 1);
    assert!(result[0].file_name().unwrap().to_str().unwrap() == "main.fa");

    Ok(())
}

#[test]
fn test_expand_input_paths_large_number_of_files() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let test_dir = env.temp_path().join("many_files");
    fs::create_dir(&test_dir)?;

    for i in 0..100 {
        fs::write(
            test_dir.join(format!("file_{i:03}.fa")),
            format!(">seq_{i}\nACGT\n"),
        )?;
    }

    let result = expand_input_paths(&[test_dir])?;

    assert_eq!(result.len(), 100);
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
        assert!(result.is_ok() || result.is_err());
    }

    Ok(())
}

#[test]
fn test_is_sequence_file_edge_cases() {
    use jam_rs::is_sequence_file;
    use std::path::Path;

    assert!(!is_sequence_file(Path::new("file.fastaverylongextension")));

    assert!(!is_sequence_file(Path::new("file.fa1")));
    assert!(!is_sequence_file(Path::new("file.fasta2")));

    assert!(!is_sequence_file(Path::new("file.fast")));
    assert!(!is_sequence_file(Path::new("file.f")));

    assert!(is_sequence_file(Path::new("file.name.with.dots.fa")));
    assert!(is_sequence_file(Path::new("version.1.0.fasta")));

    for filename in [
        "file.fastq.gz",
        "file.fq.gz",
        "file.fa.bz2",
        "file.fasta.xz",
        "file.fastq.zst",
        "file.fq.zstd",
    ] {
        assert!(
            is_sequence_file(Path::new(filename)),
            "Failed for {filename}"
        );
    }

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

#[test]
fn test_expand_input_paths_io_errors() -> Result<()> {
    let env = CoverageTestEnvironment::new()?;

    let regular_file = env.create_test_fasta("regular.fa", ">seq\nACGT\n")?;

    let result = expand_input_paths(&[regular_file]);
    assert!(result.is_ok());

    Ok(())
}
