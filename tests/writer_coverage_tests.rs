use anyhow::Result;
use jam_rs::writer::MergeIterator;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use tempfile::TempDir;

/// Test environment for writer tests
struct WriterTestEnvironment {
    temp_dir: TempDir,
}

impl WriterTestEnvironment {
    fn new() -> Result<Self> {
        let temp_dir = tempfile::tempdir()?;
        Ok(Self { temp_dir })
    }

    fn temp_path(&self) -> &std::path::Path {
        self.temp_dir.path()
    }

    /// Create a test binary file with (u64, u64) pairs in big-endian format
    fn create_binary_file(&self, name: &str, data: &[(u64, u64)]) -> Result<PathBuf> {
        let path = self.temp_path().join(name);
        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);

        for &(first, second) in data {
            writer.write_all(&first.to_be_bytes())?;
            writer.write_all(&second.to_be_bytes())?;
        }

        writer.flush()?;
        Ok(path)
    }
}

#[test]
fn test_merge_iterator_basic() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Create test files with sorted data
    let file1_data = vec![(1, 100), (3, 300), (5, 500)];
    let file2_data = vec![(2, 200), (4, 400), (6, 600)];
    let file3_data = vec![(7, 700), (8, 800), (9, 900)];

    let file1 = env.create_binary_file("file1.bin", &file1_data)?;
    let file2 = env.create_binary_file("file2.bin", &file2_data)?;
    let file3 = env.create_binary_file("file3.bin", &file3_data)?;

    // Test MergeIterator::new and iteration
    let file_paths = vec![file1, file2, file3];
    let memory_limit_gb = 1.0;
    let iterator = MergeIterator::new(file_paths, memory_limit_gb);

    // Collect all results
    let results: Vec<(u64, u64)> = iterator.collect();

    // Verify merged results are sorted
    let expected = vec![
        (1, 100),
        (2, 200),
        (3, 300),
        (4, 400),
        (5, 500),
        (6, 600),
        (7, 700),
        (8, 800),
        (9, 900),
    ];

    assert_eq!(results, expected, "Merged results should be sorted");

    Ok(())
}

#[test]
fn test_merge_iterator_empty_files() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Create one empty file and one with data
    let empty_file = env.create_binary_file("empty.bin", &[])?;
    let data_file = env.create_binary_file("data.bin", &[(1, 100), (2, 200)])?;

    let file_paths = vec![empty_file, data_file];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();
    assert_eq!(
        results,
        vec![(1, 100), (2, 200)],
        "Should handle empty files gracefully"
    );

    Ok(())
}

#[test]
fn test_merge_iterator_single_file() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    let data = vec![(10, 1000), (20, 2000), (30, 3000)];
    let file = env.create_binary_file("single.bin", &data)?;

    let file_paths = vec![file];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();
    assert_eq!(results, data, "Single file should work correctly");

    Ok(())
}

#[test]
fn test_merge_iterator_duplicates() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Create files with duplicate keys
    let file1_data = vec![(1, 100), (2, 200), (3, 300)];
    let file2_data = vec![(1, 101), (2, 201), (4, 400)]; // Same keys with different values

    let file1 = env.create_binary_file("dup1.bin", &file1_data)?;
    let file2 = env.create_binary_file("dup2.bin", &file2_data)?;

    let file_paths = vec![file1, file2];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();

    // Should preserve order and include duplicates
    assert_eq!(
        results.len(),
        6,
        "Should include all entries including duplicates"
    );
    assert!(
        results.iter().any(|&(k, v)| k == 1 && v == 100),
        "Should have first duplicate"
    );
    assert!(
        results.iter().any(|&(k, v)| k == 1 && v == 101),
        "Should have second duplicate"
    );

    Ok(())
}

#[test]
fn test_merge_iterator_large_files() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Create larger files to test batch loading
    let mut file1_data = Vec::new();
    let mut file2_data = Vec::new();

    // Generate 5000 entries per file to exceed typical batch sizes
    for i in 0..5000 {
        file1_data.push((i * 2, i * 2 + 1000)); // Even numbers
        file2_data.push((i * 2 + 1, i * 2 + 1001)); // Odd numbers
    }

    let file1 = env.create_binary_file("large1.bin", &file1_data)?;
    let file2 = env.create_binary_file("large2.bin", &file2_data)?;

    let file_paths = vec![file1, file2];
    let memory_limit_gb = 0.001; // Very small memory limit to force frequent batch reloading
    let iterator = MergeIterator::new(file_paths, memory_limit_gb);

    let results: Vec<(u64, u64)> = iterator.collect();

    assert_eq!(
        results.len(),
        10000,
        "Should merge all entries from large files"
    );

    // Verify that results are still sorted
    for i in 1..results.len() {
        assert!(
            results[i - 1].0 <= results[i].0,
            "Results should remain sorted even with batching"
        );
    }

    Ok(())
}

#[test]
fn test_merge_iterator_memory_limits() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    let data1 = vec![(1, 10), (3, 30), (5, 50)];
    let data2 = vec![(2, 20), (4, 40), (6, 60)];

    let file1 = env.create_binary_file("mem1.bin", &data1)?;
    let file2 = env.create_binary_file("mem2.bin", &data2)?;

    let file_paths = vec![file1, file2];

    // Test with different memory limits
    for memory_limit in [0.001, 0.01, 0.1, 1.0] {
        let iterator = MergeIterator::new(file_paths.clone(), memory_limit);
        let results: Vec<(u64, u64)> = iterator.collect();

        let expected = vec![(1, 10), (2, 20), (3, 30), (4, 40), (5, 50), (6, 60)];
        assert_eq!(
            results, expected,
            "Results should be consistent regardless of memory limit"
        );
    }

    Ok(())
}

#[test]
fn test_merge_iterator_reload_functionality() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Create files with enough data to trigger reloading
    let mut data1 = Vec::new();
    let mut data2 = Vec::new();

    for i in 0..2000 {
        data1.push((i * 3, i * 3 + 100)); // Every 3rd number
        data2.push((i * 3 + 1, i * 3 + 101)); // Every 3rd number + 1
    }

    let file1 = env.create_binary_file("reload1.bin", &data1)?;
    let file2 = env.create_binary_file("reload2.bin", &data2)?;

    let file_paths = vec![file1, file2];
    let memory_limit_gb = 0.01; // Small limit to force reloading
    let mut iterator = MergeIterator::new(file_paths, memory_limit_gb);

    // Manually iterate to test reload functionality
    let mut count = 0;
    let mut last_value = 0;

    while let Some((key, _value)) = iterator.next() {
        assert!(key >= last_value, "Keys should be non-decreasing");
        last_value = key;
        count += 1;

        // Stop early to avoid long test execution
        if count >= 1000 {
            break;
        }
    }

    assert!(
        count == 1000,
        "Should successfully iterate through many entries with reloading"
    );

    Ok(())
}

#[test]
fn test_merge_iterator_edge_cases() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Test with all empty files
    let empty1 = env.create_binary_file("empty1.bin", &[])?;
    let empty2 = env.create_binary_file("empty2.bin", &[])?;

    let file_paths = vec![empty1, empty2];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();
    assert_eq!(
        results.len(),
        0,
        "All empty files should produce no results"
    );

    // Test with files containing identical data
    let identical_data = vec![(1, 100), (2, 200)];
    let file1 = env.create_binary_file("identical1.bin", &identical_data)?;
    let file2 = env.create_binary_file("identical2.bin", &identical_data)?;

    let file_paths = vec![file1, file2];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();
    assert_eq!(
        results.len(),
        4,
        "Identical files should produce duplicate entries"
    );

    Ok(())
}

#[test]
fn test_merge_iterator_boundary_values() -> Result<()> {
    let env = WriterTestEnvironment::new()?;

    // Test with extreme values
    let data1 = vec![(0, 0), (u64::MAX / 2, 1000)];
    let data2 = vec![(1, 1), (u64::MAX, 2000)];

    let file1 = env.create_binary_file("boundary1.bin", &data1)?;
    let file2 = env.create_binary_file("boundary2.bin", &data2)?;

    let file_paths = vec![file1, file2];
    let iterator = MergeIterator::new(file_paths, 1.0);

    let results: Vec<(u64, u64)> = iterator.collect();

    assert_eq!(results.len(), 4, "Should handle boundary values correctly");
    assert_eq!(results[0], (0, 0), "Should handle zero correctly");
    assert_eq!(
        results[3],
        (u64::MAX, 2000),
        "Should handle MAX value correctly"
    );

    Ok(())
}
