use crate::bias::HashBiasTable;
use crate::format::{
    BUCKET_BITS, BUCKET_COUNT, BUCKET_TABLE_SIZE, BucketMeta, DATA_START, ENTRY_SIZE, Entry,
    FLAG_HAS_BIAS_TABLE, FormatError, HEADER_SIZE, Header, PAGE_SIZE, bucket_id,
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

    #[error("Invalid metadata: {message}")]
    InvalidMetadata { message: String },
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
        validate_header_layout(&header)?;

        let table_end = HEADER_SIZE + BUCKET_TABLE_SIZE;
        if mmap.len() < table_end {
            return Err(ReaderError::FileTooSmall {
                expected: table_end,
                actual: mmap.len(),
            });
        }

        let bucket_table: Vec<BucketMeta> =
            bytemuck::cast_slice(&mmap[HEADER_SIZE..table_end]).to_vec();

        for (i, meta) in bucket_table.iter().enumerate() {
            validate_bucket_meta(meta, i, mmap.len())?;
        }
        validate_bucket_totals(&header, &bucket_table)?;

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
            let (offset, end) = checked_file_range(
                header.bias_table_offset,
                header.bias_table_size,
                mmap.len(),
                "bias table",
            )?;
            let bias_data = &mmap[offset..end];
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
            let (offset, end) = checked_file_range(
                header.sample_names_offset,
                header.sample_names_size,
                mmap.len(),
                "sample names",
            )?;
            if end > mmap.len() {
                return Err(ReaderError::FileTooSmall {
                    expected: end,
                    actual: mmap.len(),
                });
            }
            let names = parse_sample_names(&mmap[offset..end], header.sample_count)?;
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
            let (offset, end) = checked_file_range(
                header.sample_sizes_offset,
                header.sample_sizes_size,
                mmap.len(),
                "sample sizes",
            )?;
            let size = end - offset;
            let expected_size = (header.sample_count as usize)
                .checked_mul(8)
                .ok_or_else(|| ReaderError::InvalidSampleData {
                    message: "sample sizes section size overflow".to_string(),
                })?;
            if size != expected_size {
                return Err(ReaderError::InvalidSampleData {
                    message: format!(
                        "sample sizes section size mismatch: got {}, expected {}",
                        size, expected_size
                    ),
                });
            }
            if end > mmap.len() {
                return Err(ReaderError::FileTooSmall {
                    expected: end,
                    actual: mmap.len(),
                });
            }
            parse_sample_sizes(&mmap[offset..end])
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

        let region_start =
            usize::try_from(meta.filter_offset).map_err(|_| ReaderError::InvalidMetadata {
                message: format!("bucket {bucket_idx} filter offset does not fit usize"),
            })?;
        let filter_size =
            usize::try_from(meta.filter_size).map_err(|_| ReaderError::InvalidMetadata {
                message: format!("bucket {bucket_idx} filter size does not fit usize"),
            })?;
        let entry_count =
            usize::try_from(meta.entry_count).map_err(|_| ReaderError::InvalidMetadata {
                message: format!("bucket {bucket_idx} entry count does not fit usize"),
            })?;
        let entry_size =
            entry_count
                .checked_mul(ENTRY_SIZE)
                .ok_or_else(|| ReaderError::InvalidMetadata {
                    message: format!("bucket {bucket_idx} entry byte size overflows"),
                })?;
        let data_size =
            filter_size
                .checked_add(entry_size)
                .ok_or_else(|| ReaderError::InvalidMetadata {
                    message: format!("bucket {bucket_idx} mapped region size overflows"),
                })?;

        let page_start = region_start & !(PAGE_SIZE - 1);
        let data_offset = region_start - page_start;
        let mmap_len =
            data_offset
                .checked_add(data_size)
                .ok_or_else(|| ReaderError::InvalidMetadata {
                    message: format!("bucket {bucket_idx} mmap length overflows"),
                })?;

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
            let filter_data = &mmap[filter_data_start..filter_data_start + filter_size];

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
            filter_size,
            entry_count,
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
        let idx = lower_bound_hash(entries, hash);
        entries.get(idx).is_some_and(|entry| entry.hash == hash)
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

        let start = lower_bound_hash(entries, hash);

        entries[start..]
            .iter()
            .take_while(move |e| e.hash == hash)
            .map(|e| e.sample_id)
    }
}

#[inline]
pub(crate) fn lower_bound_hash(entries: &[Entry], key: u64) -> usize {
    const MAX_INTERPOLATION_PROBES: usize = 8;
    const BINARY_SEARCH_CUTOFF: usize = 64;

    let mut lo = 0usize;
    let mut hi = entries.len();
    if hi == 0 {
        return 0;
    }

    if entries[lo].hash >= key {
        return lo;
    }
    if entries[hi - 1].hash < key {
        return hi;
    }

    for _ in 0..MAX_INTERPOLATION_PROBES {
        if hi - lo <= BINARY_SEARCH_CUTOFF {
            break;
        }

        let lo_hash = entries[lo].hash;
        let hi_hash = entries[hi - 1].hash;
        if lo_hash >= key {
            return lo;
        }
        if hi_hash < key {
            return hi;
        }
        if lo_hash == hi_hash {
            break;
        }

        let width = hi - lo - 1;
        let span = hi_hash - lo_hash;
        let rel = key.saturating_sub(lo_hash);
        let probe = lo + ((rel as u128 * width as u128) / span as u128) as usize;
        let probe_hash = entries[probe].hash;

        if probe_hash < key {
            lo = probe + 1;
        } else {
            hi = probe + 1;
        }
    }

    lower_bound_hash_binary(entries, key, lo, hi)
}

#[inline]
fn lower_bound_hash_binary(entries: &[Entry], key: u64, mut lo: usize, mut hi: usize) -> usize {
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if entries[mid].hash < key {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

fn checked_file_range(
    offset: u64,
    size: u64,
    file_len: usize,
    label: &str,
) -> Result<(usize, usize), ReaderError> {
    let start = usize::try_from(offset).map_err(|_| ReaderError::InvalidMetadata {
        message: format!("{label} offset does not fit usize: {offset}"),
    })?;
    let size = usize::try_from(size).map_err(|_| ReaderError::InvalidMetadata {
        message: format!("{label} size does not fit usize: {size}"),
    })?;
    let end = start
        .checked_add(size)
        .ok_or_else(|| ReaderError::InvalidMetadata {
            message: format!("{label} range overflows: offset={offset}, size={size}"),
        })?;

    if end > file_len {
        return Err(ReaderError::FileTooSmall {
            expected: end,
            actual: file_len,
        });
    }

    Ok((start, end))
}

fn invalid_metadata(message: impl Into<String>) -> ReaderError {
    ReaderError::InvalidMetadata {
        message: message.into(),
    }
}

fn validate_header_layout(header: &Header) -> Result<(), ReaderError> {
    if header.bucket_bits != BUCKET_BITS {
        return Err(invalid_metadata(format!(
            "invalid bucket_bits: got {}, expected {}",
            header.bucket_bits, BUCKET_BITS
        )));
    }
    if header.bucket_table_offset != HEADER_SIZE as u64 {
        return Err(invalid_metadata(format!(
            "invalid bucket_table_offset: got {}, expected {}",
            header.bucket_table_offset, HEADER_SIZE
        )));
    }
    Ok(())
}

fn validate_bucket_totals(header: &Header, bucket_table: &[BucketMeta]) -> Result<(), ReaderError> {
    let mut entry_count = 0u64;
    let mut entries_size = 0u64;
    let mut filters_size = 0u64;

    for (bucket_idx, meta) in bucket_table.iter().enumerate() {
        entry_count = entry_count.checked_add(meta.entry_count).ok_or_else(|| {
            invalid_metadata(format!("bucket {bucket_idx} entry count total overflows"))
        })?;
        let bucket_entry_size =
            meta.entry_count
                .checked_mul(ENTRY_SIZE as u64)
                .ok_or_else(|| {
                    invalid_metadata(format!("bucket {bucket_idx} entry byte size overflows"))
                })?;
        entries_size = entries_size.checked_add(bucket_entry_size).ok_or_else(|| {
            invalid_metadata(format!("bucket {bucket_idx} entry byte total overflows"))
        })?;
        filters_size = filters_size.checked_add(meta.filter_size).ok_or_else(|| {
            invalid_metadata(format!("bucket {bucket_idx} filter byte total overflows"))
        })?;
    }

    if entry_count != header.entry_count {
        return Err(invalid_metadata(format!(
            "entry_count mismatch: bucket table sums to {}, header reports {}",
            entry_count, header.entry_count
        )));
    }
    if entries_size != header.entries_size {
        return Err(invalid_metadata(format!(
            "entries_size mismatch: bucket table sums to {}, header reports {}",
            entries_size, header.entries_size
        )));
    }
    if filters_size != header.filters_size {
        return Err(invalid_metadata(format!(
            "filters_size mismatch: bucket table sums to {}, header reports {}",
            filters_size, header.filters_size
        )));
    }

    Ok(())
}

fn validate_bucket_meta(
    meta: &BucketMeta,
    bucket_idx: usize,
    file_len: usize,
) -> Result<(), ReaderError> {
    let entry_count =
        usize::try_from(meta.entry_count).map_err(|_| ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!("entry count does not fit usize: {}", meta.entry_count),
        })?;
    let entry_size =
        entry_count
            .checked_mul(ENTRY_SIZE)
            .ok_or_else(|| ReaderError::InvalidFilter {
                bucket: bucket_idx,
                message: format!("entry range size overflows: count={}", meta.entry_count),
            })?;
    if entry_size > 0 && meta.entry_offset < DATA_START as u64 {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "entry range starts before data section: {} < {}",
                meta.entry_offset, DATA_START
            ),
        });
    }
    if meta.filter_size > 0 && meta.filter_offset < DATA_START as u64 {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "filter range starts before data section: {} < {}",
                meta.filter_offset, DATA_START
            ),
        });
    }
    checked_file_range(
        meta.entry_offset,
        entry_size as u64,
        file_len,
        &format!("bucket {bucket_idx} entries"),
    )?;
    checked_file_range(
        meta.filter_offset,
        meta.filter_size,
        file_len,
        &format!("bucket {bucket_idx} filter"),
    )?;

    if meta.entry_count > 0 && meta.filter_size == 0 {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: "non-empty bucket is missing its filter".to_string(),
        });
    }

    if meta.filter_size > 0 && meta.entry_count > 0 {
        let expected_entry_offset = meta
            .filter_offset
            .checked_add(meta.filter_size)
            .ok_or_else(|| ReaderError::InvalidFilter {
                bucket: bucket_idx,
                message: "filter range end overflows".to_string(),
            })?;
        if meta.entry_offset != expected_entry_offset {
            return Err(ReaderError::InvalidFilter {
                bucket: bucket_idx,
                message: format!(
                    "entry range must start immediately after filter: got {}, expected {}",
                    meta.entry_offset, expected_entry_offset
                ),
            });
        }
    }
    Ok(())
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

    if offset != data.len() {
        return Err(ReaderError::InvalidSampleData {
            message: format!(
                "sample names section has {} trailing bytes",
                data.len() - offset
            ),
        });
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
    let (start, end) = checked_file_range(
        meta.filter_offset,
        meta.filter_size,
        mmap.len(),
        &format!("bucket {bucket_idx} filter"),
    )?;

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

    let expected_size = 8usize
        .checked_add(descriptor_size)
        .and_then(|size| size.checked_add(fingerprints_size))
        .ok_or_else(|| ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: "filter data size overflows".to_string(),
        })?;
    if data.len() != expected_size {
        return Err(ReaderError::InvalidFilter {
            bucket: bucket_idx,
            message: format!(
                "filter data size mismatch: got {}, expected {}",
                data.len(),
                expected_size
            ),
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
    use crate::format::{MAGIC, VERSION};
    use crate::writer::{BuildConfig, build, serialize_filter};
    use std::io::Write;
    use tempfile::NamedTempFile;
    use xorf::BinaryFuse8;

    fn make_fasta(seqs: &[(&str, &str)]) -> NamedTempFile {
        let mut f = NamedTempFile::with_suffix(".fa").unwrap();
        for (name, seq) in seqs {
            writeln!(f, ">{name}").unwrap();
            writeln!(f, "{seq}").unwrap();
        }
        f
    }

    fn valid_test_header(entry_count: u64, unique_hash_count: u64, sample_count: u32) -> Header {
        Header {
            magic: MAGIC,
            version: VERSION,
            flags: 0,
            entry_count,
            unique_hash_count,
            sample_count,
            bucket_count: BUCKET_COUNT as u16,
            bucket_bits: 8,
            entry_size: ENTRY_SIZE as u8,
            hash_threshold: 1024,
            kmer_size: 11,
            _param_reserved: [0; 3],
            min_entropy: 0.0,
            bucket_table_offset: HEADER_SIZE as u64,
            entries_offset: DATA_START as u64,
            filters_offset: DATA_START as u64,
            bias_table_offset: 0,
            entries_size: entry_count * ENTRY_SIZE as u64,
            filters_size: 0,
            bias_table_size: 0,
            sample_names_offset: 0,
            sample_names_size: 0,
            sample_sizes_offset: 0,
            sample_sizes_size: 0,
            _padding: [0; 16],
        }
    }

    fn write_minimal_jam(bucket_table: &[BucketMeta], body_len: usize) -> NamedTempFile {
        let mut file = NamedTempFile::with_suffix(".jam").unwrap();
        let header = valid_test_header(0, 0, 0);
        let mut bytes = vec![0u8; DATA_START + body_len];
        bytes[..HEADER_SIZE].copy_from_slice(bytemuck::bytes_of(&header));
        bytes[HEADER_SIZE..HEADER_SIZE + BUCKET_TABLE_SIZE]
            .copy_from_slice(bytemuck::cast_slice(bucket_table));
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_reader_search_duplicate_run_returns_all_samples() {
        let hash = 512u64;
        let entry_count = 64usize;
        let entries: Vec<Entry> = (0..entry_count)
            .map(|sample_id| Entry::new(hash, sample_id as u32))
            .collect();
        let filter = BinaryFuse8::try_from(&[hash][..]).unwrap();
        let filter_bytes = serialize_filter(&filter);

        let entry_offset = DATA_START + filter_bytes.len();
        let file_len = entry_offset + entries.len() * ENTRY_SIZE;
        let mut bucket_table = vec![BucketMeta::default(); BUCKET_COUNT];
        bucket_table[0] = BucketMeta {
            entry_offset: entry_offset as u64,
            entry_count: entries.len() as u64,
            filter_offset: DATA_START as u64,
            filter_size: filter_bytes.len() as u64,
        };

        let mut header = valid_test_header(entries.len() as u64, 1, entries.len() as u32);
        header.filters_size = filter_bytes.len() as u64;
        let mut bytes = vec![0u8; file_len];
        bytes[..HEADER_SIZE].copy_from_slice(bytemuck::bytes_of(&header));
        bytes[HEADER_SIZE..HEADER_SIZE + BUCKET_TABLE_SIZE]
            .copy_from_slice(bytemuck::cast_slice(&bucket_table));
        bytes[DATA_START..entry_offset].copy_from_slice(&filter_bytes);
        bytes[entry_offset..file_len].copy_from_slice(bytemuck::cast_slice(&entries));

        let mut file = NamedTempFile::with_suffix(".jam").unwrap();
        file.write_all(&bytes).unwrap();
        file.flush().unwrap();

        let reader = JamReader::open(file.path()).unwrap();
        let samples: Vec<_> = reader.search(hash).collect();
        assert_eq!(samples, (0..entry_count as u32).collect::<Vec<_>>());
    }

    #[test]
    fn test_lower_bound_hash_handles_skewed_distribution() {
        let mut hashes: Vec<u64> = (0..9000).map(|i| i * 2).collect();
        hashes.extend((0..1000).map(|i| 1_000_000 + i * 2));
        hashes.splice(4500..4500, [9000, 9000, 9000]);

        let entries: Vec<Entry> = hashes
            .iter()
            .enumerate()
            .map(|(sample_id, &hash)| Entry::new(hash, sample_id as u32))
            .collect();

        for key in [
            0, 1, 8999, 9000, 9001, 20_000, 1_000_000, 1_001_999, 2_000_000,
        ] {
            let expected = entries.partition_point(|entry| entry.hash < key);
            assert_eq!(lower_bound_hash(&entries, key), expected, "key {key}");
        }
    }

    #[test]
    fn test_reader_open_rejects_malformed_bucket_entry_range() {
        let mut bucket_table = vec![BucketMeta::default(); BUCKET_COUNT];
        bucket_table[0] = BucketMeta {
            entry_offset: u64::MAX,
            entry_count: 1,
            filter_offset: 0,
            filter_size: 0,
        };
        let file = write_minimal_jam(&bucket_table, 0);

        assert!(JamReader::open(file.path()).is_err());
    }

    #[test]
    fn test_reader_open_rejects_malformed_bucket_filter_range() {
        let mut bucket_table = vec![BucketMeta::default(); BUCKET_COUNT];
        bucket_table[0] = BucketMeta {
            entry_offset: 0,
            entry_count: 0,
            filter_offset: DATA_START as u64,
            filter_size: 1,
        };
        let file = write_minimal_jam(&bucket_table, 0);

        assert!(JamReader::open(file.path()).is_err());
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
