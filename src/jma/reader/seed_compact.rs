//! Range-backed access to compact JMA seed sections.
//!
//! The reader fetches the fixed header and directory prefix first, then reads
//! only key and position pages selected by exact query hashes.

use crate::jma::seeds::{
    COLLECTION_HEADER_SIZE, LookupOptions, SeedBinding, SeedCollectionEntry, SeedCollectionHeader,
    SeedHeader, SeedKeyPage, SeedLookupResult, SeedPositionPage, SeedPostingHeader, SeedPrefix,
    decode_collection_directory, decode_collection_header, decode_directory_prefix,
    decode_key_page_records, decode_position_posting, decode_seed_header, membership_contains,
    validate_collection_entries,
};
use crate::jma::{JmaError, JmaResult, SeedOccurrence, SeedQuery};
use crate::resource::{ByteRange, RangeReader};
use std::collections::{BTreeMap, BTreeSet};

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
    posting_headers: Vec<SeedPostingHeader>,
    position_pages: Vec<SeedPositionPage>,
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
            usize::try_from(header.key_data_offset).map_err(|_| JmaError::OffsetOverflow)?;
        let directory = read(
            &resource,
            section_offset,
            u64::try_from(directory_length).map_err(|_| JmaError::OffsetOverflow)?,
        )?;
        let (header, prefixes, key_pages, posting_headers, position_pages) =
            decode_directory_prefix(&directory, section_length)?;
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
            posting_headers,
            position_pages,
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
                let posting = *self
                    .posting_headers
                    .get(usize::try_from(global).map_err(|_| JmaError::OffsetOverflow)?)
                    .ok_or_else(|| {
                        JmaError::CorruptSection(
                            "compact seed posting header is unavailable".to_string(),
                        )
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

        let mut position_cache = BTreeMap::new();
        let mut checked_pages = BTreeSet::new();
        let mut decoded: BTreeMap<u64, Vec<SeedOccurrence>> = BTreeMap::new();
        for (match_index, posting) in pending {
            let page = *self
                .position_pages
                .get(
                    usize::try_from(posting.position_page_id)
                        .map_err(|_| JmaError::OffsetOverflow)?,
                )
                .ok_or_else(|| {
                    JmaError::CorruptSection(
                        "compact seed position page is unavailable".to_string(),
                    )
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
                if checked_pages.insert(page.page_id) {
                    result.metrics.position_pages_read =
                        result.metrics.position_pages_read.saturating_add(1);
                    result.metrics.bytes_read = result
                        .metrics
                        .bytes_read
                        .saturating_add(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
                }
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

    fn find_prefix(&self, hash: u64) -> JmaResult<Option<SeedPrefix>> {
        let prefix = crate::jma::format::hash_prefix(hash, self.header.bucket_bits)?;
        Ok(self
            .prefixes
            .binary_search_by_key(&prefix, |entry| entry.hash_prefix)
            .ok()
            .map(|index| self.prefixes[index]))
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
        .read_range(ByteRange::new(offset, length)?)
        .map_err(JmaError::Resource)
}
