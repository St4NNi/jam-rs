use crate::bias::BiasTable;
use crate::format::{
    bucket_id, BucketMeta, Entry, FormatError, Header, BUCKET_COUNT, BUCKET_TABLE_SIZE,
    ENTRY_SIZE, FLAG_HAS_BIAS_TABLE, HEADER_SIZE,
};
use crate::writer::FILTER_DESCRIPTOR_SIZE;
use memmap2::Mmap;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use xorf::{BinaryFuse8Ref, Filter, FilterRef};

#[derive(Debug, thiserror::Error)]
pub enum ReaderError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("Format error: {0}")]
    Format(#[from] FormatError),

    #[error("Invalid filter data for bucket {bucket}: {message}")]
    InvalidFilter { bucket: usize, message: String },

    #[error("File too small: expected at least {expected} bytes, got {actual}")]
    FileTooSmall { expected: usize, actual: usize },
}

#[derive(Debug, Clone)]
pub struct ReaderStats {
    pub entry_count: u64,
    pub unique_hash_count: u64,
    pub sample_count: u32,
    pub file_size: u64,
    pub kmer_size: u8,
    pub hash_threshold: u64,
    pub bucket_entry_counts: [u64; BUCKET_COUNT],
    pub has_bias_table: bool,
}

/// Heap-owned filter data that can be used to construct BinaryFuse8Ref on demand.
/// The descriptor is small (20 bytes) and stays hot in L1 cache.
struct LoadedFilter {
    descriptor: [u8; FILTER_DESCRIPTOR_SIZE],
    fingerprints: Vec<u8>,
}

impl LoadedFilter {
    /// Check if the filter might contain the given hash.
    #[inline]
    fn contains(&self, hash: &u64) -> bool {
        BinaryFuse8Ref::from_dma(&self.descriptor, &self.fingerprints).contains(hash)
    }
}

pub struct JamReader {
    mmap: Mmap,
    header: Header,
    bucket_table: Vec<BucketMeta>,
    filters: Vec<Option<LoadedFilter>>, // None for empty buckets
    bias_table: Option<Arc<BiasTable>>,
}

impl JamReader {
    /// Open a .jam file for reading.
    ///
    /// This loads the header, bucket table, and all filters into memory.
    /// Filters are copied to heap for fast access; entries remain mmap'd.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, ReaderError> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        // Validate minimum size for header
        if mmap.len() < HEADER_SIZE {
            return Err(ReaderError::FileTooSmall {
                expected: HEADER_SIZE,
                actual: mmap.len(),
            });
        }

        // Parse and validate header
        let header: Header = *bytemuck::from_bytes(&mmap[..HEADER_SIZE]);
        header.validate()?;

        // Validate size for bucket table
        let table_end = HEADER_SIZE + BUCKET_TABLE_SIZE;
        if mmap.len() < table_end {
            return Err(ReaderError::FileTooSmall {
                expected: table_end,
                actual: mmap.len(),
            });
        }

        // Parse bucket table
        let bucket_table: Vec<BucketMeta> =
            bytemuck::cast_slice(&mmap[HEADER_SIZE..table_end]).to_vec();

        // Load all filters into heap memory
        let mut filters = Vec::with_capacity(BUCKET_COUNT);
        for (i, meta) in bucket_table.iter().enumerate() {
            if meta.filter_size == 0 {
                filters.push(None);
                continue;
            }

            let filter = parse_filter(&mmap, meta, i)?;
            filters.push(Some(filter));
        }

        // Load embedded bias table if present
        let bias_table = if header.flags & FLAG_HAS_BIAS_TABLE != 0
            && header.bias_table_offset > 0
            && header.bias_table_size > 0
        {
            let offset = header.bias_table_offset as usize;
            let size = header.bias_table_size as usize;
            if offset + size > mmap.len() {
                return Err(ReaderError::FileTooSmall {
                    expected: offset + size,
                    actual: mmap.len(),
                });
            }
            let bias_data = &mmap[offset..offset + size];
            let table = BiasTable::from_bytes(bias_data).map_err(|e| {
                ReaderError::InvalidFilter {
                    bucket: 0,
                    message: format!("Failed to parse embedded bias table: {}", e),
                }
            })?;
            Some(Arc::new(table))
        } else {
            None
        };

        Ok(Self {
            mmap,
            header,
            bucket_table,
            filters,
            bias_table,
        })
    }

    /// Get the hash threshold used for FracMinHash filtering.
    #[inline]
    pub fn threshold(&self) -> u64 {
        self.header.hash_threshold
    }

    /// Get the k-mer size used when building this database.
    #[inline]
    pub fn kmer_size(&self) -> u8 {
        self.header.kmer_size
    }

    /// Get the embedded bias table, if one was stored during database creation.
    #[inline]
    pub fn bias_table(&self) -> Option<Arc<BiasTable>> {
        self.bias_table.clone()
    }

    /// Check if the database has an embedded bias table.
    #[inline]
    pub fn has_bias_table(&self) -> bool {
        self.bias_table.is_some()
    }

    /// Get statistics about this database.
    pub fn stats(&self) -> ReaderStats {
        let mut bucket_entry_counts = [0u64; BUCKET_COUNT];
        for (i, meta) in self.bucket_table.iter().enumerate() {
            bucket_entry_counts[i] = meta.entry_count;
        }

        ReaderStats {
            entry_count: self.header.entry_count,
            unique_hash_count: self.header.unique_hash_count,
            sample_count: self.header.sample_count,
            file_size: self.mmap.len() as u64,
            kmer_size: self.header.kmer_size,
            hash_threshold: self.header.hash_threshold,
            bucket_entry_counts,
            has_bias_table: self.bias_table.is_some(),
        }
    }

    /// Get entries for a specific bucket as a slice.
    #[inline]
    pub fn bucket_entries(&self, bucket_idx: usize) -> &[Entry] {
        let meta = &self.bucket_table[bucket_idx];
        if meta.entry_count == 0 {
            return &[];
        }

        let start = meta.entry_offset as usize;
        let end = start + (meta.entry_count as usize) * ENTRY_SIZE;
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    /// Check if a hash exists in the database (filter check + verification).
    #[inline]
    pub fn contains(&self, hash: u64) -> bool {
        let bucket_idx = bucket_id(hash);

        // Filter check first
        if let Some(ref filter) = self.filters[bucket_idx] {
            if !filter.contains(&hash) {
                return false; // Definitely not present
            }
        } else {
            return false; // Empty bucket
        }

        // Filter says "maybe" - verify with actual lookup
        let entries = self.bucket_entries(bucket_idx);
        self.interpolation_search(entries, hash).is_some()
    }

    /// Search for all sample IDs containing a hash.
    ///
    /// Returns an iterator over sample IDs. Uses filter check + interpolation search.
    #[inline]
    pub fn search(&self, hash: u64) -> impl Iterator<Item = u32> + '_ {
        let bucket_idx = bucket_id(hash);

        // Filter check first
        let dominated = if let Some(ref filter) = self.filters[bucket_idx] {
            filter.contains(&hash)
        } else {
            false
        };

        let entries = if dominated {
            self.bucket_entries(bucket_idx)
        } else {
            &[]
        };

        // Find start position using interpolation, then iterate matches
        let start = if entries.is_empty() {
            0
        } else {
            self.interpolation_find_start(entries, hash)
        };

        entries[start..]
            .iter()
            .skip_while(move |e| e.hash < hash)
            .take_while(move |e| e.hash == hash)
            .map(|e| e.sample_id)
    }

    /// Interpolation search to find if hash exists.
    fn interpolation_search(&self, entries: &[Entry], key: u64) -> Option<usize> {
        if entries.is_empty() {
            return None;
        }

        let start = self.interpolation_find_start(entries, key);

        // Forward scan to find exact match
        for (i, entry) in entries[start..].iter().enumerate() {
            if entry.hash == key {
                return Some(start + i);
            }
            if entry.hash > key {
                break;
            }
        }

        None
    }

    /// Find starting position for a hash using interpolation.
    ///
    /// Uses the known uniform distribution of FracMinHash to estimate position.
    /// Returns a position that's at or before the target hash.
    #[inline]
    fn interpolation_find_start(&self, entries: &[Entry], key: u64) -> usize {
        let count = entries.len();
        let threshold = self.threshold();

        // Hashes are uniform in [0, threshold), so position ≈ key/threshold * count
        let est = ((key as u128 * count as u128) / threshold as u128) as usize;

        // Bias low to almost always land before target (avoids backtracking)
        let est = est.saturating_sub(16).min(count - 1);

        // Rare backtrack if we overshot
        if entries[est].hash > key {
            let mut i = est;
            while i > 0 && entries[i - 1].hash >= key {
                i -= 1;
            }
            i
        } else {
            est
        }
    }
}

/// Parse filter data from mmap into heap-owned LoadedFilter.
fn parse_filter(mmap: &Mmap, meta: &BucketMeta, bucket_idx: usize) -> Result<LoadedFilter, ReaderError> {
    let start = meta.filter_offset as usize;
    let end = start + meta.filter_size as usize;

    if end > mmap.len() {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!("filter extends beyond file: {}..{} > {}", start, end, mmap.len()),
        });
    }

    let data = &mmap[start..end];

    // Format: descriptor_size(u32) + fingerprints_size(u32) + descriptor + fingerprints
    if data.len() < 8 {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: "filter data too small for header".to_string(),
        });
    }

    let descriptor_size = u32::from_le_bytes(data[0..4].try_into().unwrap()) as usize;
    let fingerprints_size = u32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;

    if descriptor_size != FILTER_DESCRIPTOR_SIZE {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "unexpected descriptor size: {} (expected {})",
                descriptor_size, FILTER_DESCRIPTOR_SIZE
            ),
        });
    }

    let expected_size = 8 + descriptor_size + fingerprints_size;
    if data.len() < expected_size {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "filter data too small: {} < {}",
                data.len(),
                expected_size
            ),
        });
    }

    let descriptor_data = &data[8..8 + descriptor_size];
    let fingerprints_data = &data[8 + descriptor_size..8 + descriptor_size + fingerprints_size];

    // Copy to heap-owned storage
    let mut descriptor = [0u8; FILTER_DESCRIPTOR_SIZE];
    descriptor.copy_from_slice(descriptor_data);

    Ok(LoadedFilter {
        descriptor,
        fingerprints: fingerprints_data.to_vec(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::{build, BuildConfig};
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

    #[test]
    fn test_reader_open() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
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
        let stats = reader.stats();

        assert!(stats.entry_count > 0);
        assert_eq!(stats.sample_count, 1);
        assert_eq!(stats.kmer_size, 11);
    }

    #[test]
    fn test_reader_search() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1, // Keep all hashes
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();

        // Get a hash that we know exists
        let entries = reader.bucket_entries(0);
        if !entries.is_empty() {
            let test_hash = entries[0].hash;
            assert!(reader.contains(test_hash));

            let samples: Vec<_> = reader.search(test_hash).collect();
            assert!(!samples.is_empty());
        }
    }

    #[test]
    fn test_reader_nonexistent_hash() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1000, // Keep only ~0.1% of hashes
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();

        // Search for a hash that almost certainly doesn't exist
        // (threshold is u64::MAX / 1000, so most hashes are above it)
        let fake_hash = u64::MAX - 1;
        assert!(!reader.contains(fake_hash));

        let samples: Vec<_> = reader.search(fake_hash).collect();
        assert!(samples.is_empty());
    }

    #[test]
    fn test_reader_multiple_samples() {
        let input = make_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "ATCGATCGATCGATCGATCGATCGATCGATCG"), // Same sequence
        ]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            singleton: true, // Each sequence is a separate sample
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();
        assert_eq!(reader.stats().sample_count, 2);

        // Find a hash that should be in both samples
        for bucket_idx in 0..BUCKET_COUNT {
            let entries = reader.bucket_entries(bucket_idx);
            if entries.len() >= 2 {
                // Find a hash with multiple sample hits
                let test_hash = entries[0].hash;
                let samples: Vec<_> = reader.search(test_hash).collect();
                if samples.len() == 2 {
                    assert!(samples.contains(&0) || samples.contains(&1));
                    return; // Test passed
                }
            }
        }
    }

    #[test]
    fn test_reader_bucket_entries() {
        let input = make_fasta(&[("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG")]);
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

        let reader = JamReader::open(&output_path).unwrap();

        // Verify entries are sorted within each bucket
        for bucket_idx in 0..BUCKET_COUNT {
            let entries = reader.bucket_entries(bucket_idx);
            for window in entries.windows(2) {
                assert!(window[0] <= window[1], "Entries not sorted in bucket {}", bucket_idx);
            }

            // Verify all entries belong to this bucket
            for entry in entries {
                assert_eq!(bucket_id(entry.hash), bucket_idx);
            }
        }
    }
}
