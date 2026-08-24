//! Range-backed access to compact JMA seed sections.
//!
//! The reader fetches the fixed header and directory prefix first, then reads
//! only key and position pages selected by exact query hashes.

use crate::jma::seeds::{
    COLLECTION_HEADER_SIZE, LookupOptions, POSITION_PAGE_RECORD_SIZE, POSTING_HEADER_RECORD_SIZE,
    SeedBinding, SeedCollectionEntry, SeedCollectionHeader, SeedHeader, SeedKeyPage,
    SeedLookupResult, SeedPositionPage, SeedPostingHeader, SeedPrefix, decode_collection_directory,
    decode_collection_header, decode_key_lookup_prefix, decode_key_page_records,
    decode_position_page_record, decode_position_posting, decode_posting_header_record,
    decode_seed_header, membership_contains, validate_collection_entries,
};
use crate::jma::{JmaError, JmaResult, SeedOccurrence, SeedQuery};
use crate::resource::{ByteRange, RangeReader};
use std::collections::BTreeMap;

/// A range-backed seed collection directory. Opening the collection reads
/// only the fixed header and directory; scheme blobs remain unopened.
pub struct RemoteSeedCollection<R> {
    resource: R,
    section_offset: u64,
    header: SeedCollectionHeader,
    entries: Vec<SeedCollectionEntry>,
}

impl<R: RangeReader + Clone> RemoteSeedCollection<R> {
    pub fn open(resource: R, section_offset: u64, section_length: u64) -> JmaResult<Self> {
        let header_bytes = read(
            &resource,
            section_offset,
            u64::try_from(COLLECTION_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?,
        )?;
        let header = decode_collection_header(&header_bytes, section_length)?;
        let directory_end = header
            .directory_offset
            .checked_add(header.directory_length)
            .ok_or(JmaError::OffsetOverflow)?;
        let directory = read(&resource, section_offset, directory_end)?;
        let entries = decode_collection_directory(&directory, &header)?;
        validate_collection_entries(&header, &entries)?;
        Ok(Self {
            resource,
            section_offset,
            header,
            entries,
        })
    }

    #[must_use]
    pub const fn header(&self) -> &SeedCollectionHeader {
        &self.header
    }

    #[must_use]
    pub fn entries(&self) -> &[SeedCollectionEntry] {
        &self.entries
    }

    pub fn entry(&self, scheme_id: u32) -> JmaResult<SeedCollectionEntry> {
        self.entries
            .binary_search_by_key(&scheme_id, |entry| entry.scheme_id)
            .map(|index| self.entries[index])
            .map_err(|_| {
                JmaError::CorruptSection(format!("compact seed scheme {scheme_id} is unavailable"))
            })
    }

    pub fn open_scheme(
        &self,
        scheme_id: u32,
        expected: Option<SeedBinding>,
    ) -> JmaResult<RemoteSeedIndex<R>> {
        let entry = self.entry(scheme_id)?;
        if let Some(expected) = expected
            && (expected.scheme_id != entry.scheme_id
                || expected.k != entry.k
                || expected.scale != entry.scale)
        {
            return Err(JmaError::CorruptSection(
                "compact seed scheme binding does not match collection entry".to_string(),
            ));
        }
        let offset = self
            .section_offset
            .checked_add(entry.offset)
            .ok_or(JmaError::OffsetOverflow)?;
        let index = RemoteSeedIndex::open(self.resource.clone(), offset, entry.length, expected)?;
        if index.header().scheme_id != entry.scheme_id
            || index.header().k != entry.k
            || index.header().scale != entry.scale
            || index.header().section_checksum != entry.checksum
        {
            return Err(JmaError::CorruptSection(
                "compact seed collection entry identity mismatch".to_string(),
            ));
        }
        Ok(index)
    }
}

/// A compact seed section opened against a local or remote range resource.
pub struct RemoteSeedIndex<R> {
    resource: R,
    section_offset: u64,
    header: SeedHeader,
    filter: Vec<u8>,
    prefixes: Vec<SeedPrefix>,
    key_pages: Vec<SeedKeyPage>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
#[allow(dead_code)]
pub struct ExportedSeedRecord {
    pub query: SeedQuery,
    pub class: crate::jma::seeds::PostingClass,
    pub occurrences: Vec<SeedOccurrence>,
}

impl<R: RangeReader> RemoteSeedIndex<R> {
    pub fn open(
        resource: R,
        section_offset: u64,
        section_length: u64,
        expected: Option<SeedBinding>,
    ) -> JmaResult<Self> {
        let header_bytes = read(&resource, section_offset, 320)?;
        let header = decode_seed_header(&header_bytes, section_length)?;
        if let Some(expected) = expected {
            expected.matches(&header)?;
        }
        let directory_length =
            usize::try_from(header.posting_header_offset).map_err(|_| JmaError::OffsetOverflow)?;
        let directory = read(
            &resource,
            section_offset,
            u64::try_from(directory_length).map_err(|_| JmaError::OffsetOverflow)?,
        )?;
        let (header, prefixes, key_pages) = decode_key_lookup_prefix(&directory, section_length)?;
        let filter_range = range(header.filter_offset, header.filter_length)?;
        let filter = directory
            .get(filter_range)
            .ok_or_else(|| {
                JmaError::CorruptSection("compact seed filter is unavailable".to_string())
            })?
            .to_vec();
        Ok(Self {
            resource,
            section_offset,
            header,
            filter,
            prefixes,
            key_pages,
        })
    }

    #[must_use]
    pub const fn header(&self) -> &SeedHeader {
        &self.header
    }

    pub fn lookup_batch(
        &self,
        queries: &[SeedQuery],
        options: LookupOptions,
    ) -> JmaResult<SeedLookupResult> {
        let mut result = SeedLookupResult::default();
        let mut page_queries: BTreeMap<u32, Vec<(usize, SeedQuery)>> = BTreeMap::new();
        for &query in queries {
            crate::jma::seeds::validate_query(query, self.header.k)?;
            if !membership_contains(&self.filter, query.hash) {
                result.metrics.filter_misses = result.metrics.filter_misses.saturating_add(1);
                continue;
            }
            let Some(prefix) = self.find_prefix(query.hash)? else {
                result.metrics.filter_misses = result.metrics.filter_misses.saturating_add(1);
                continue;
            };
            let page_end = prefix
                .first_page
                .checked_add(prefix.page_count)
                .ok_or(JmaError::OffsetOverflow)?;
            let match_index = result.matches.len();
            result.matches.push(crate::jma::seeds::SeedMatch {
                query,
                occurrence_count: 0,
                class: crate::jma::seeds::PostingClass::Suppressed,
                occurrences: Vec::new(),
            });
            for page_id in prefix.first_page..page_end {
                page_queries
                    .entry(page_id)
                    .or_default()
                    .push((match_index, query));
            }
        }

        let mut page_cache = BTreeMap::new();
        let mut pending = Vec::new();
        for (page_id, requests) in page_queries {
            let page = *self
                .key_pages
                .get(usize::try_from(page_id).map_err(|_| JmaError::OffsetOverflow)?)
                .ok_or_else(|| {
                    JmaError::CorruptSection("compact seed key page is unavailable".to_string())
                })?;
            let bytes = read(
                &self.resource,
                self.section_offset
                    .checked_add(page.data_offset)
                    .ok_or(JmaError::OffsetOverflow)?,
                page.data_length,
            )?;
            result.metrics.key_pages_read = result.metrics.key_pages_read.saturating_add(1);
            result.metrics.bytes_read = result
                .metrics
                .bytes_read
                .saturating_add(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
            let keys =
                decode_key_page_records(&bytes, page, self.header.k, self.header.key_encoding)?;
            page_cache.insert(page_id, keys);
            let mut page_matches = Vec::new();
            for (match_index, query) in requests {
                let keys = page_cache.get(&page_id).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed key page cache miss".to_string())
                })?;
                result.metrics.keys_tested = result
                    .metrics
                    .keys_tested
                    .saturating_add(u64::try_from(keys.len()).unwrap_or(u64::MAX));
                let Some(local) = keys
                    .binary_search_by_key(&(query.hash, query.canonical_kmer), |key| *key)
                    .ok()
                else {
                    continue;
                };
                let global = page
                    .first_key
                    .checked_add(u64::try_from(local).map_err(|_| JmaError::OffsetOverflow)?)
                    .ok_or(JmaError::OffsetOverflow)?;
                page_matches.push((match_index, global));
            }
            let header_indices = page_matches
                .iter()
                .map(|(_, key_index)| *key_index)
                .collect::<Vec<_>>();
            let (posting_headers, header_bytes) = self.read_posting_headers(&header_indices)?;
            result.metrics.bytes_read = result.metrics.bytes_read.saturating_add(header_bytes);
            for (match_index, global) in page_matches {
                let posting = *posting_headers.get(&global).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed posting header cache miss".to_string())
                })?;
                let target = result.matches.get_mut(match_index).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed match index is invalid".to_string())
                })?;
                target.occurrence_count = posting.occurrence_count;
                target.class = posting.class;
                result.metrics.postings_seen = result.metrics.postings_seen.saturating_add(1);
                if !options.include_positions
                    || options
                        .max_occurrences
                        .is_some_and(|limit| posting.occurrence_count > limit)
                    || posting.class == crate::jma::seeds::PostingClass::Suppressed
                {
                    result.metrics.postings_skipped =
                        result.metrics.postings_skipped.saturating_add(1);
                } else {
                    pending.push((match_index, posting));
                }
            }
        }

        let position_page_ids = pending
            .iter()
            .map(|(_, posting)| posting.position_page_id)
            .collect::<Vec<_>>();
        let (position_descriptors, descriptor_bytes) =
            self.read_position_pages(&position_page_ids)?;
        result.metrics.bytes_read = result.metrics.bytes_read.saturating_add(descriptor_bytes);
        let mut position_cache = BTreeMap::new();
        let mut decoded: BTreeMap<u64, Vec<SeedOccurrence>> = BTreeMap::new();
        for (match_index, posting) in pending {
            let page = *position_descriptors
                .get(&posting.position_page_id)
                .ok_or_else(|| {
                    JmaError::CorruptSection("compact seed position page cache miss".to_string())
                })?;
            if let std::collections::btree_map::Entry::Vacant(entry) =
                position_cache.entry(page.page_id)
            {
                let bytes = read(
                    &self.resource,
                    self.section_offset
                        .checked_add(page.data_offset)
                        .ok_or(JmaError::OffsetOverflow)?,
                    page.data_length,
                )?;
                if crate::jma::format::checksum(&bytes) != page.checksum {
                    return Err(JmaError::ChecksumMismatch(format!(
                        "compact seed position page {}",
                        page.page_id
                    )));
                }
                result.metrics.position_pages_read =
                    result.metrics.position_pages_read.saturating_add(1);
                result.metrics.bytes_read = result
                    .metrics
                    .bytes_read
                    .saturating_add(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
                entry.insert(bytes);
            }
            let bytes = position_cache.get(&page.page_id).ok_or_else(|| {
                JmaError::CorruptSection("compact seed position page cache miss".to_string())
            })?;
            let occurrences = if let Some(occurrences) = decoded.get(&posting.key_index) {
                occurrences.clone()
            } else {
                result.metrics.position_payload_reads =
                    result.metrics.position_payload_reads.saturating_add(1);
                let start = usize::try_from(posting.position_offset)
                    .map_err(|_| JmaError::OffsetOverflow)?;
                let end = start
                    .checked_add(
                        usize::try_from(posting.position_length)
                            .map_err(|_| JmaError::OffsetOverflow)?,
                    )
                    .ok_or(JmaError::OffsetOverflow)?;
                let payload = bytes.get(start..end).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed posting range is invalid".to_string())
                })?;
                let occurrences = decode_position_posting(
                    payload,
                    posting,
                    self.header.k,
                    options.max_occurrences.map(|value| value as usize),
                )?;
                decoded.insert(posting.key_index, occurrences.clone());
                occurrences
            };
            if let Some(target) = result.matches.get_mut(match_index) {
                target.occurrences = occurrences;
            }
        }
        result.matches.retain(|entry| entry.occurrence_count != 0);
        Ok(result)
    }

    pub fn lookup(&self, query: SeedQuery, options: LookupOptions) -> JmaResult<SeedLookupResult> {
        self.lookup_batch(&[query], options)
    }

    /// Decode every exact key and occurrence for immutable router builds.
    #[allow(dead_code)]
    pub fn export_records(&self) -> JmaResult<Vec<ExportedSeedRecord>> {
        let mut records = Vec::with_capacity(
            usize::try_from(self.header.key_count).map_err(|_| JmaError::OffsetOverflow)?,
        );
        let mut position_descriptors = BTreeMap::new();
        let mut position_pages = BTreeMap::new();
        for page in &self.key_pages {
            let bytes = read(
                &self.resource,
                self.section_offset
                    .checked_add(page.data_offset)
                    .ok_or(JmaError::OffsetOverflow)?,
                page.data_length,
            )?;
            let keys =
                decode_key_page_records(&bytes, *page, self.header.k, self.header.key_encoding)?;
            for (local_index, (hash, packed)) in keys.into_iter().enumerate() {
                let key_index = page
                    .first_key
                    .checked_add(u64::try_from(local_index).map_err(|_| JmaError::OffsetOverflow)?)
                    .ok_or(JmaError::OffsetOverflow)?;
                let posting = self.read_posting_header(key_index)?;
                let descriptor =
                    if let Some(descriptor) = position_descriptors.get(&posting.position_page_id) {
                        *descriptor
                    } else {
                        let descriptor = self.read_position_page(posting.position_page_id)?;
                        position_descriptors.insert(posting.position_page_id, descriptor);
                        descriptor
                    };
                if let std::collections::btree_map::Entry::Vacant(entry) =
                    position_pages.entry(descriptor.page_id)
                {
                    let bytes = read(
                        &self.resource,
                        self.section_offset
                            .checked_add(descriptor.data_offset)
                            .ok_or(JmaError::OffsetOverflow)?,
                        descriptor.data_length,
                    )?;
                    if crate::jma::format::checksum(&bytes) != descriptor.checksum {
                        return Err(JmaError::ChecksumMismatch(format!(
                            "compact seed position page {}",
                            descriptor.page_id
                        )));
                    }
                    entry.insert(bytes);
                }
                let page_bytes = position_pages.get(&descriptor.page_id).ok_or_else(|| {
                    JmaError::CorruptSection(
                        "compact seed export position page is unavailable".to_string(),
                    )
                })?;
                let start = usize::try_from(posting.position_offset)
                    .map_err(|_| JmaError::OffsetOverflow)?;
                let end = start
                    .checked_add(
                        usize::try_from(posting.position_length)
                            .map_err(|_| JmaError::OffsetOverflow)?,
                    )
                    .ok_or(JmaError::OffsetOverflow)?;
                let payload = page_bytes.get(start..end).ok_or_else(|| {
                    JmaError::CorruptSection(
                        "compact seed export posting range is invalid".to_string(),
                    )
                })?;
                records.push(ExportedSeedRecord {
                    query: SeedQuery {
                        k: self.header.k,
                        hash,
                        canonical_kmer: packed,
                    },
                    class: posting.class,
                    occurrences: decode_position_posting(payload, posting, self.header.k, None)?,
                });
            }
        }
        Ok(records)
    }

    fn find_prefix(&self, hash: u64) -> JmaResult<Option<SeedPrefix>> {
        let prefix = crate::jma::format::hash_prefix(hash, self.header.bucket_bits)?;
        Ok(self
            .prefixes
            .binary_search_by_key(&prefix, |entry| entry.hash_prefix)
            .ok()
            .map(|index| self.prefixes[index]))
    }

    fn read_posting_header(&self, key_index: u64) -> JmaResult<SeedPostingHeader> {
        if key_index >= self.header.posting_count {
            return Err(JmaError::CorruptSection(
                "compact seed posting header index is out of bounds".to_string(),
            ));
        }
        let relative = key_index
            .checked_mul(POSTING_HEADER_RECORD_SIZE as u64)
            .and_then(|offset| self.header.posting_header_offset.checked_add(offset))
            .ok_or(JmaError::OffsetOverflow)?;
        let bytes = read(
            &self.resource,
            self.section_offset
                .checked_add(relative)
                .ok_or(JmaError::OffsetOverflow)?,
            POSTING_HEADER_RECORD_SIZE as u64,
        )?;
        decode_posting_header_record(&bytes, key_index, &self.header)
    }

    fn read_posting_headers(
        &self,
        key_indices: &[u64],
    ) -> JmaResult<(BTreeMap<u64, SeedPostingHeader>, u64)> {
        let mut indices = key_indices.to_vec();
        indices.sort_unstable();
        indices.dedup();
        let Some(&first) = indices.first() else {
            return Ok((BTreeMap::new(), 0));
        };
        let last = *indices.last().ok_or(JmaError::OffsetOverflow)?;
        if last >= self.header.posting_count {
            return Err(JmaError::CorruptSection(
                "compact seed posting header index is out of bounds".to_string(),
            ));
        }
        let record_size = POSTING_HEADER_RECORD_SIZE as u64;
        let relative = first
            .checked_mul(record_size)
            .and_then(|offset| self.header.posting_header_offset.checked_add(offset))
            .ok_or(JmaError::OffsetOverflow)?;
        let length = last
            .checked_sub(first)
            .and_then(|span| span.checked_add(1))
            .and_then(|count| count.checked_mul(record_size))
            .ok_or(JmaError::OffsetOverflow)?;
        let bytes = read(
            &self.resource,
            self.section_offset
                .checked_add(relative)
                .ok_or(JmaError::OffsetOverflow)?,
            length,
        )?;
        let mut headers = BTreeMap::new();
        for key_index in indices {
            let start = usize::try_from(
                key_index
                    .checked_sub(first)
                    .and_then(|index| index.checked_mul(record_size))
                    .ok_or(JmaError::OffsetOverflow)?,
            )
            .map_err(|_| JmaError::OffsetOverflow)?;
            let end = start
                .checked_add(POSTING_HEADER_RECORD_SIZE)
                .ok_or(JmaError::OffsetOverflow)?;
            let record = bytes.get(start..end).ok_or_else(|| {
                JmaError::CorruptSection(
                    "compact seed posting header batch is truncated".to_string(),
                )
            })?;
            headers.insert(
                key_index,
                decode_posting_header_record(record, key_index, &self.header)?,
            );
        }
        Ok((headers, length))
    }

    fn read_position_page(&self, page_id: u32) -> JmaResult<SeedPositionPage> {
        if page_id >= self.header.position_page_count {
            return Err(JmaError::CorruptSection(
                "compact seed position page index is out of bounds".to_string(),
            ));
        }
        let relative = u64::from(page_id)
            .checked_mul(POSITION_PAGE_RECORD_SIZE as u64)
            .and_then(|offset| self.header.position_page_dir_offset.checked_add(offset))
            .ok_or(JmaError::OffsetOverflow)?;
        let bytes = read(
            &self.resource,
            self.section_offset
                .checked_add(relative)
                .ok_or(JmaError::OffsetOverflow)?,
            POSITION_PAGE_RECORD_SIZE as u64,
        )?;
        decode_position_page_record(&bytes, page_id, &self.header)
    }

    fn read_position_pages(
        &self,
        page_ids: &[u32],
    ) -> JmaResult<(BTreeMap<u32, SeedPositionPage>, u64)> {
        const MAX_GAP_RECORDS: u32 = 64;
        const MAX_SPAN_RECORDS: u32 = 1_024;

        let mut ids = page_ids.to_vec();
        ids.sort_unstable();
        ids.dedup();
        if ids
            .last()
            .is_some_and(|page_id| *page_id >= self.header.position_page_count)
        {
            return Err(JmaError::CorruptSection(
                "compact seed position page index is out of bounds".to_string(),
            ));
        }

        let mut pages = BTreeMap::new();
        let mut bytes_read = 0_u64;
        let mut start_index = 0_usize;
        while start_index < ids.len() {
            let first = ids[start_index];
            let mut end_index = start_index + 1;
            while end_index < ids.len() {
                let previous = ids[end_index - 1];
                let next = ids[end_index];
                if next.saturating_sub(previous) > MAX_GAP_RECORDS
                    || next.saturating_sub(first) >= MAX_SPAN_RECORDS
                {
                    break;
                }
                end_index += 1;
            }
            let last = ids[end_index - 1];
            let record_size = POSITION_PAGE_RECORD_SIZE as u64;
            let relative = u64::from(first)
                .checked_mul(record_size)
                .and_then(|offset| self.header.position_page_dir_offset.checked_add(offset))
                .ok_or(JmaError::OffsetOverflow)?;
            let length = u64::from(last - first + 1)
                .checked_mul(record_size)
                .ok_or(JmaError::OffsetOverflow)?;
            let bytes = read(
                &self.resource,
                self.section_offset
                    .checked_add(relative)
                    .ok_or(JmaError::OffsetOverflow)?,
                length,
            )?;
            bytes_read = bytes_read.saturating_add(length);
            for &page_id in &ids[start_index..end_index] {
                let start = usize::try_from(u64::from(page_id - first) * record_size)
                    .map_err(|_| JmaError::OffsetOverflow)?;
                let end = start
                    .checked_add(POSITION_PAGE_RECORD_SIZE)
                    .ok_or(JmaError::OffsetOverflow)?;
                let record = bytes.get(start..end).ok_or_else(|| {
                    JmaError::CorruptSection(
                        "compact seed position page batch is truncated".to_string(),
                    )
                })?;
                pages.insert(
                    page_id,
                    decode_position_page_record(record, page_id, &self.header)?,
                );
            }
            start_index = end_index;
        }
        Ok((pages, bytes_read))
    }
}

fn range(offset: u64, length: u64) -> JmaResult<std::ops::Range<usize>> {
    let end = offset.checked_add(length).ok_or(JmaError::OffsetOverflow)?;
    Ok(
        usize::try_from(offset).map_err(|_| JmaError::OffsetOverflow)?
            ..usize::try_from(end).map_err(|_| JmaError::OffsetOverflow)?,
    )
}

fn read<R: RangeReader>(resource: &R, offset: u64, length: u64) -> JmaResult<Vec<u8>> {
    resource
        .read_range_bytes(ByteRange::new(offset, length)?)
        .map(|bytes| bytes.as_ref().to_vec())
        .map_err(JmaError::Resource)
}
