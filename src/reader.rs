use crate::bias::HashBiasTable;
use crate::format::{
    BUCKET_COUNT, BUCKET_TABLE_SIZE, BucketMeta, ENTRY_SIZE, Entry, FLAG_HAS_BIAS_TABLE,
    FormatError, HEADER_SIZE, Header, PAGE_SIZE, bucket_id,
};
use crate::writer::FILTER_DESCRIPTOR_SIZE;
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use xorf::{BinaryFuse8Ref, Filter, FilterRef};

#[cfg(unix)]
use memmap2::{Advice, UncheckedAdvice};

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

    #[error("Invalid sample data: {message}")]
    InvalidSampleData { message: String },
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

struct FilterMeta {
    descriptor_offset: usize,
    fingerprints_offset: usize,
    fingerprints_size: usize,
}

pub struct BucketFilter<'a> {
    mmap: &'a [u8],
    meta: &'a FilterMeta,
}

impl BucketFilter<'_> {
    #[inline]
    pub fn contains(&self, hash: &u64) -> bool {
        let descriptor = &self.mmap
            [self.meta.descriptor_offset..self.meta.descriptor_offset + FILTER_DESCRIPTOR_SIZE];
        let fingerprints = &self.mmap[self.meta.fingerprints_offset
            ..self.meta.fingerprints_offset + self.meta.fingerprints_size];
        BinaryFuse8Ref::from_dma(descriptor, fingerprints).contains(hash)
    }
}

struct CachedFilterMeta {
    descriptor_start: usize,
    descriptor_size: usize,
    fingerprints_start: usize,
    fingerprints_size: usize,
}

pub struct BucketRegion {
    mmap: Mmap,
    data_offset: usize,
    filter_size: usize,
    entry_count: usize,
    filter_meta: Option<CachedFilterMeta>,
}

impl BucketRegion {
    #[inline]
    pub fn filter_contains(&self, hash: &u64) -> bool {
        let meta = match &self.filter_meta {
            Some(m) => m,
            None => return false,
        };

        let descriptor =
            &self.mmap[meta.descriptor_start..meta.descriptor_start + meta.descriptor_size];
        let fingerprints =
            &self.mmap[meta.fingerprints_start..meta.fingerprints_start + meta.fingerprints_size];
        BinaryFuse8Ref::from_dma(descriptor, fingerprints).contains(hash)
    }

    #[inline]
    pub fn entries(&self) -> &[Entry] {
        if self.entry_count == 0 {
            return &[];
        }
        let start = self.data_offset + self.filter_size;
        let end = start + self.entry_count * ENTRY_SIZE;
        bytemuck::cast_slice(&self.mmap[start..end])
    }

    #[inline]
    pub fn entry_count(&self) -> usize {
        self.entry_count
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.filter_size == 0 && self.entry_count == 0
    }
}

pub struct JamReader {
    file: Arc<File>,
    mmap: Mmap,
    header: Header,
    bucket_table: Vec<BucketMeta>,
    filters: Vec<Option<FilterMeta>>,
    bias_table: Option<Arc<HashBiasTable>>,
    sample_names: Vec<String>,
    sample_sizes: Vec<u64>,
}

impl JamReader {
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, ReaderError> {
        let file = Arc::new(File::open(path.as_ref())?);
        let mmap = unsafe { Mmap::map(file.as_ref())? };

        if mmap.len() < HEADER_SIZE {
            return Err(ReaderError::FileTooSmall {
                expected: HEADER_SIZE,
                actual: mmap.len(),
            });
        }

        let header: Header = *bytemuck::from_bytes(&mmap[..HEADER_SIZE]);
        header.validate()?;

        let table_end = HEADER_SIZE + BUCKET_TABLE_SIZE;
        if mmap.len() < table_end {
            return Err(ReaderError::FileTooSmall {
                expected: table_end,
                actual: mmap.len(),
            });
        }

        let bucket_table: Vec<BucketMeta> =
            bytemuck::cast_slice(&mmap[HEADER_SIZE..table_end]).to_vec();

        let mut filters = Vec::with_capacity(BUCKET_COUNT);
        for (i, meta) in bucket_table.iter().enumerate() {
            if meta.filter_size == 0 {
                filters.push(None);
                continue;
            }

            let filter_meta = parse_filter_meta(&mmap, meta, i)?;
            filters.push(Some(filter_meta));
        }

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
            let table =
                HashBiasTable::from_bytes(bias_data).map_err(|e| ReaderError::InvalidFilter {
                    bucket: 0,
                    message: format!("Failed to parse embedded bias table: {}", e),
                })?;
            Some(Arc::new(table))
        } else {
            None
        };

        let sample_names = if header.sample_names_offset > 0 && header.sample_names_size > 0 {
            let offset = header.sample_names_offset as usize;
            let size = header.sample_names_size as usize;
            if offset + size > mmap.len() {
                return Err(ReaderError::FileTooSmall {
                    expected: offset + size,
                    actual: mmap.len(),
                });
            }
            let names = parse_sample_names(&mmap[offset..offset + size], header.sample_count)?;
            if names.len() != header.sample_count as usize {
                return Err(ReaderError::InvalidSampleData {
                    message: format!(
                        "sample names count mismatch: got {}, expected {}",
                        names.len(),
                        header.sample_count
                    ),
                });
            }
            names
        } else {
            (0..header.sample_count)
                .map(|i| format!("sample_{}", i))
                .collect()
        };

        let sample_sizes = if header.sample_sizes_offset > 0 && header.sample_sizes_size > 0 {
            let offset = header.sample_sizes_offset as usize;
            let size = header.sample_sizes_size as usize;
            let expected_size = header.sample_count as usize * 8;
            if size != expected_size {
                return Err(ReaderError::InvalidSampleData {
                    message: format!(
                        "sample sizes section size mismatch: got {}, expected {}",
                        size, expected_size
                    ),
                });
            }
            if offset + size > mmap.len() {
                return Err(ReaderError::FileTooSmall {
                    expected: offset + size,
                    actual: mmap.len(),
                });
            }
            parse_sample_sizes(&mmap[offset..offset + size])
        } else {
            vec![0u64; header.sample_count as usize]
        };

        Ok(Self {
            file,
            mmap,
            header,
            bucket_table,
            filters,
            bias_table,
            sample_names,
            sample_sizes,
        })
    }

    pub fn open_bucket_region(&self, bucket_idx: usize) -> Result<BucketRegion, ReaderError> {
        let meta = &self.bucket_table[bucket_idx];

        if meta.filter_size == 0 && meta.entry_count == 0 {
            let empty_mmap = MmapOptions::new().len(1).map_anon()?.make_read_only()?;
            return Ok(BucketRegion {
                mmap: empty_mmap,
                data_offset: 0,
                filter_size: 0,
                entry_count: 0,
                filter_meta: None,
            });
        }

        let region_start = meta.filter_offset as usize;
        let data_size = meta.filter_size as usize + (meta.entry_count as usize) * ENTRY_SIZE;

        let page_start = region_start & !(PAGE_SIZE - 1);
        let data_offset = region_start - page_start;
        let mmap_len = data_offset + data_size;

        let mmap = unsafe {
            MmapOptions::new()
                .offset(page_start as u64)
                .len(mmap_len)
                .map(self.file.as_ref())?
        };

        #[cfg(unix)]
        {
            let _ = mmap.advise(Advice::Sequential);
        }

        let filter_meta = if meta.filter_size > 0 {
            let filter_data_start = data_offset;
            let filter_data =
                &mmap[filter_data_start..filter_data_start + meta.filter_size as usize];

            if filter_data.len() >= 8 {
                let descriptor_size =
                    u32::from_le_bytes(filter_data[0..4].try_into().unwrap()) as usize;
                let fingerprints_size =
                    u32::from_le_bytes(filter_data[4..8].try_into().unwrap()) as usize;

                if descriptor_size != FILTER_DESCRIPTOR_SIZE {
                    return Err(ReaderError::InvalidFilter {
                        bucket: bucket_idx,
                        message: format!(
                            "unexpected descriptor size in bucket region: {} (expected {})",
                            descriptor_size, FILTER_DESCRIPTOR_SIZE
                        ),
                    });
                }

                if filter_data.len() >= 8 + descriptor_size + fingerprints_size {
                    Some(CachedFilterMeta {
                        descriptor_start: filter_data_start + 8,
                        descriptor_size,
                        fingerprints_start: filter_data_start + 8 + descriptor_size,
                        fingerprints_size,
                    })
                } else {
                    return Err(ReaderError::InvalidFilter {
                        bucket: bucket_idx,
                        message: format!(
                            "filter data truncated: need {} bytes, have {}",
                            8 + descriptor_size + fingerprints_size,
                            filter_data.len()
                        ),
                    });
                }
            } else {
                return Err(ReaderError::InvalidFilter {
                    bucket: bucket_idx,
                    message: format!(
                        "filter header too small: need 8 bytes, have {}",
                        filter_data.len()
                    ),
                });
            }
        } else {
            None
        };

        Ok(BucketRegion {
            mmap,
            data_offset,
            filter_size: meta.filter_size as usize,
            entry_count: meta.entry_count as usize,
            filter_meta,
        })
    }

    #[inline]
    pub fn bucket_meta(&self, bucket_idx: usize) -> &BucketMeta {
        &self.bucket_table[bucket_idx]
    }

    #[inline]
    pub fn threshold(&self) -> u64 {
        self.header.hash_threshold
    }

    #[inline]
    pub fn kmer_size(&self) -> u8 {
        self.header.kmer_size
    }

    #[inline]
    pub fn min_entropy(&self) -> f64 {
        self.header.min_entropy as f64
    }

    #[inline]
    pub fn bias_table(&self) -> Option<Arc<HashBiasTable>> {
        self.bias_table.clone()
    }

    #[inline]
    pub fn has_bias_table(&self) -> bool {
        self.bias_table.is_some()
    }

    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    pub fn sample_name(&self, id: u32) -> Option<&str> {
        self.sample_names.get(id as usize).map(|s| s.as_str())
    }

    pub fn sample_sizes(&self) -> &[u64] {
        &self.sample_sizes
    }

    pub fn sample_size(&self, id: u32) -> Option<u64> {
        self.sample_sizes.get(id as usize).copied()
    }

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

    #[inline]
    pub fn bucket_entry_byte_range(&self, bucket_idx: usize) -> (usize, usize) {
        let meta = &self.bucket_table[bucket_idx];
        let start = meta.entry_offset as usize;
        let end = start + (meta.entry_count as usize) * ENTRY_SIZE;
        (start, end)
    }

    #[inline]
    pub fn bucket_filter_byte_range(&self, bucket_idx: usize) -> (usize, usize) {
        let meta = &self.bucket_table[bucket_idx];
        let start = meta.filter_offset as usize;
        let end = start + meta.filter_size as usize;
        (start, end)
    }

    #[cfg(unix)]
    pub fn release_pages(&self, start: usize, end: usize) {
        if start >= end {
            return;
        }
        let page_start = start & !(PAGE_SIZE - 1);
        let page_end = (end + PAGE_SIZE - 1) & !(PAGE_SIZE - 1);
        let len = page_end.saturating_sub(page_start);
        if len > 0 && page_end <= self.mmap.len() {
            let _ = unsafe {
                self.mmap
                    .unchecked_advise_range(UncheckedAdvice::DontNeed, page_start, len)
            };
        }
    }

    #[cfg(not(unix))]
    pub fn release_pages(&self, _start: usize, _end: usize) {}

    pub fn release_bucket(&self, bucket_idx: usize) {
        let (filter_start, filter_end) = self.bucket_filter_byte_range(bucket_idx);
        let (entry_start, entry_end) = self.bucket_entry_byte_range(bucket_idx);
        self.release_pages(filter_start, filter_end);
        self.release_pages(entry_start, entry_end);
    }

    #[cfg(unix)]
    pub fn advise_random(&self) {
        let _ = self.mmap.advise(Advice::Random);
    }

    #[cfg(not(unix))]
    pub fn advise_random(&self) {}

    #[inline]
    pub fn bucket_filter(&self, bucket_idx: usize) -> Option<BucketFilter<'_>> {
        self.filters[bucket_idx].as_ref().map(|meta| BucketFilter {
            mmap: &self.mmap,
            meta,
        })
    }

    #[inline]
    pub fn contains(&self, hash: u64) -> bool {
        let bucket_idx = bucket_id(hash);

        if let Some(filter) = self.bucket_filter(bucket_idx) {
            if !filter.contains(&hash) {
                return false;
            }
        } else {
            return false;
        }

        let entries = self.bucket_entries(bucket_idx);
        self.interpolation_search(entries, hash).is_some()
    }

    #[inline]
    pub fn search(&self, hash: u64) -> impl Iterator<Item = u32> + '_ {
        let bucket_idx = bucket_id(hash);

        let dominated = self
            .bucket_filter(bucket_idx)
            .is_some_and(|f| f.contains(&hash));

        let entries = if dominated {
            self.bucket_entries(bucket_idx)
        } else {
            &[]
        };

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

    fn interpolation_search(&self, entries: &[Entry], key: u64) -> Option<usize> {
        if entries.is_empty() {
            return None;
        }

        let start = self.interpolation_find_start(entries, key);

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

    #[inline]
    fn interpolation_find_start(&self, entries: &[Entry], key: u64) -> usize {
        let count = entries.len();
        let threshold = self.threshold();

        let est = ((key as u128 * count as u128) / threshold as u128) as usize;

        let est = est.saturating_sub(16).min(count - 1);

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

fn parse_sample_names(data: &[u8], count: u32) -> Result<Vec<String>, ReaderError> {
    let mut names = Vec::with_capacity(count as usize);
    let mut offset = 0;

    for i in 0..count {
        if offset + 2 > data.len() {
            return Err(ReaderError::InvalidSampleData {
                message: format!(
                    "truncated sample names section: cannot read length for sample {}",
                    i
                ),
            });
        }
        let len = u16::from_le_bytes(data[offset..offset + 2].try_into().unwrap()) as usize;
        offset += 2;
        if offset + len > data.len() {
            return Err(ReaderError::InvalidSampleData {
                message: format!(
                    "truncated sample names section: cannot read name for sample {} (need {} bytes, have {})",
                    i,
                    len,
                    data.len() - offset
                ),
            });
        }
        names.push(String::from_utf8_lossy(&data[offset..offset + len]).to_string());
        offset += len;
    }

    Ok(names)
}

fn parse_sample_sizes(data: &[u8]) -> Vec<u64> {
    data.chunks_exact(8)
        .map(|chunk| u64::from_le_bytes(chunk.try_into().unwrap()))
        .collect()
}

fn parse_filter_meta(
    mmap: &Mmap,
    meta: &BucketMeta,
    bucket_idx: usize,
) -> Result<FilterMeta, ReaderError> {
    let start = meta.filter_offset as usize;
    let end = start + meta.filter_size as usize;

    if end > mmap.len() {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "filter extends beyond file: {}..{} > {}",
                start,
                end,
                mmap.len()
            ),
        });
    }

    let data = &mmap[start..end];

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
            message: format!("filter data too small: {} < {}", data.len(), expected_size),
        });
    }

    Ok(FilterMeta {
        descriptor_offset: start + 8,
        fingerprints_offset: start + 8 + descriptor_size,
        fingerprints_size,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::{BuildConfig, build};
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
            fscale: 1,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();

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
            fscale: 1000,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();

        let fake_hash = u64::MAX - 1;
        assert!(!reader.contains(fake_hash));

        let samples: Vec<_> = reader.search(fake_hash).collect();
        assert!(samples.is_empty());
    }

    #[test]
    fn test_reader_multiple_samples() {
        let input = make_fasta(&[
            ("seq1", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
            ("seq2", "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        ]);
        let output_dir = tempfile::tempdir().unwrap();
        let output_path = output_dir.path().join("test.jam");

        let config = BuildConfig {
            kmer_size: 11,
            fscale: 1,
            singleton: true,
            num_threads: 1,
            memory: 1,
            ..Default::default()
        };

        build(&[input.path().to_path_buf()], &output_path, &config).unwrap();

        let reader = JamReader::open(&output_path).unwrap();
        assert_eq!(reader.stats().sample_count, 2);

        for bucket_idx in 0..BUCKET_COUNT {
            let entries = reader.bucket_entries(bucket_idx);
            if entries.len() >= 2 {
                let test_hash = entries[0].hash;
                let samples: Vec<_> = reader.search(test_hash).collect();
                if samples.len() == 2 {
                    assert!(samples.contains(&0) || samples.contains(&1));
                    return;
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

        for bucket_idx in 0..BUCKET_COUNT {
            let entries = reader.bucket_entries(bucket_idx);
            for window in entries.windows(2) {
                assert!(
                    window[0] <= window[1],
                    "Entries not sorted in bucket {}",
                    bucket_idx
                );
            }

            for entry in entries {
                assert_eq!(bucket_id(entry.hash), bucket_idx);
            }
        }
    }
}
