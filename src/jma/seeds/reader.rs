use super::builder::filter_bits;
use super::format::{
    COMPACT_HEADER_SIZE, KEY_PAGE_RECORD_SIZE, POSITION_PAGE_RECORD_SIZE,
    POSTING_FLAG_POSITION_BEARING, POSTING_HEADER_RECORD_SIZE, PREFIX_RECORD_SIZE, SeedHeader,
    SeedKeyPage, SeedPositionPage, SeedPostingHeader, SeedPrefix, decode_header, get_array, get_u8,
    get_u16, get_u32, get_u64, section_range,
};
use super::{OccurrenceEncoding, PostingClass, SeedBinding, validate_occurrence, validate_query};
use crate::jma::{JmaError, JmaResult, SeedOccurrence, SeedQuery};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LookupOptions {
    pub include_positions: bool,
    pub max_occurrences: Option<u32>,
}

impl Default for LookupOptions {
    fn default() -> Self {
        Self {
            include_positions: true,
            max_occurrences: None,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct SeedLookupMetrics {
    pub filter_misses: u64,
    pub key_pages_read: u64,
    pub position_pages_read: u64,
    pub keys_tested: u64,
    pub postings_seen: u64,
    pub postings_skipped: u64,
    pub position_payload_reads: u64,
    pub bytes_read: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedMatch {
    pub query: SeedQuery,
    pub occurrence_count: u32,
    pub class: PostingClass,
    pub occurrences: Vec<SeedOccurrence>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SeedLookupResult {
    pub matches: Vec<SeedMatch>,
    pub metrics: SeedLookupMetrics,
}

pub type SeedDirectoryPrefix = (
    SeedHeader,
    Vec<SeedPrefix>,
    Vec<SeedKeyPage>,
    Vec<SeedPostingHeader>,
    Vec<SeedPositionPage>,
);

pub type SeedKeyLookupPrefix = (SeedHeader, Vec<SeedPrefix>, Vec<SeedKeyPage>);

/// An opened compact section backed by immutable bytes. Remote readers can
/// use the same directory decoders and fetch only the selected ranges.
#[derive(Clone, Debug)]
pub struct SeedIndex {
    bytes: Vec<u8>,
    header: SeedHeader,
    filter: std::ops::Range<usize>,
    prefixes: Vec<SeedPrefix>,
    key_pages: Vec<SeedKeyPage>,
    posting_headers: Vec<SeedPostingHeader>,
    position_pages: Vec<SeedPositionPage>,
}

impl SeedIndex {
    pub fn open(bytes: Vec<u8>, expected: Option<SeedBinding>) -> JmaResult<Self> {
        let section_length = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let header = decode_header(&bytes, section_length)?;
        if let Some(expected) = expected {
            expected.matches(&header)?;
        }
        let mut checksum_input = bytes.clone();
        checksum_input[288..320].fill(0);
        if crate::jma::format::checksum(&checksum_input) != header.section_checksum {
            return Err(JmaError::ChecksumMismatch(
                "compact seed section".to_string(),
            ));
        }
        let expected_filter =
            u64::try_from(COMPACT_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
        if header.filter_offset != expected_filter
            || header.prefix_offset != add(header.filter_offset, header.filter_length)?
            || header.key_page_dir_offset
                != add(
                    header.prefix_offset,
                    u64::from(header.prefix_count)
                        .checked_mul(PREFIX_RECORD_SIZE as u64)
                        .ok_or(JmaError::OffsetOverflow)?,
                )?
            || header.posting_header_offset
                != add(
                    header.key_page_dir_offset,
                    u64::from(header.key_page_count)
                        .checked_mul(KEY_PAGE_RECORD_SIZE as u64)
                        .ok_or(JmaError::OffsetOverflow)?,
                )?
            || header.position_page_dir_offset
                != add(
                    header.posting_header_offset,
                    u64::from(header.posting_header_count)
                        .checked_mul(POSTING_HEADER_RECORD_SIZE as u64)
                        .ok_or(JmaError::OffsetOverflow)?,
                )?
            || header.key_data_offset
                != add(
                    header.position_page_dir_offset,
                    u64::from(header.position_page_count)
                        .checked_mul(POSITION_PAGE_RECORD_SIZE as u64)
                        .ok_or(JmaError::OffsetOverflow)?,
                )?
            || header.position_data_offset != add(header.key_data_offset, header.key_data_length)?
            || header.section_length
                != add(header.position_data_offset, header.position_data_length)?
        {
            return Err(JmaError::CorruptSection(
                "compact seed section ranges are inconsistent".to_string(),
            ));
        }
        let filter = section_range(
            header.filter_offset,
            header.filter_length,
            header.section_length,
            bytes.len(),
            "membership filter",
        )?;
        let prefixes = decode_prefixes(&bytes, &header)?;
        let key_pages = decode_key_pages(&bytes, &header)?;
        let posting_headers = decode_posting_headers(&bytes, &header)?;
        let position_pages = decode_position_pages(&bytes, &header)?;
        validate_directory(
            &header,
            &prefixes,
            &key_pages,
            &posting_headers,
            &position_pages,
        )?;
        Ok(Self {
            bytes,
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

    #[must_use]
    pub fn key_pages(&self) -> &[SeedKeyPage] {
        &self.key_pages
    }

    #[must_use]
    pub fn posting_headers(&self) -> &[SeedPostingHeader] {
        &self.posting_headers
    }

    /// Lookup exact packed keys, grouping all requests that use a key page or
    /// position page so each selected page is decoded once.
    pub fn lookup_batch(
        &self,
        queries: &[SeedQuery],
        options: LookupOptions,
    ) -> JmaResult<SeedLookupResult> {
        let mut result = SeedLookupResult::default();
        let mut page_queries: BTreeMap<u32, Vec<(usize, SeedQuery)>> = BTreeMap::new();
        for &query in queries {
            validate_query(query, self.header.k)?;
            if !self.filter_contains(query.hash) {
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
            for page_id in prefix.first_page..page_end {
                page_queries
                    .entry(page_id)
                    .or_default()
                    .push((result.matches.len(), query));
            }
            result.matches.push(SeedMatch {
                query,
                occurrence_count: 0,
                class: PostingClass::Suppressed,
                occurrences: Vec::new(),
            });
        }

        let mut pending = Vec::new();
        let mut page_cache = BTreeMap::new();
        for (page_id, requests) in page_queries {
            let page = self
                .key_pages
                .get(usize::try_from(page_id).map_err(|_| JmaError::OffsetOverflow)?)
                .ok_or_else(|| {
                    JmaError::CorruptSection("compact seed page is unavailable".to_string())
                })?;
            let page_bytes = self.key_page_bytes(*page)?;
            result.metrics.key_pages_read = result.metrics.key_pages_read.saturating_add(1);
            result.metrics.bytes_read = result
                .metrics
                .bytes_read
                .saturating_add(u64::try_from(page_bytes.len()).unwrap_or(u64::MAX));
            let keys = decode_key_page(page_bytes, *page, self.header.k, self.header.key_encoding)?;
            page_cache.insert(page_id, keys);
            for (match_index, query) in requests {
                let keys = page_cache.get(&page_id).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed page cache miss".to_string())
                })?;
                result.metrics.keys_tested = result
                    .metrics
                    .keys_tested
                    .saturating_add(u64::try_from(keys.len()).unwrap_or(u64::MAX));
                let Some(local) = keys
                    .binary_search_by_key(&(query.hash, query.canonical_kmer), |key| (key.0, key.1))
                    .ok()
                else {
                    continue;
                };
                let global = page
                    .first_key
                    .checked_add(u64::try_from(local).map_err(|_| JmaError::OffsetOverflow)?)
                    .ok_or(JmaError::OffsetOverflow)?;
                let posting = self
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
                let skip = !options.include_positions
                    || options
                        .max_occurrences
                        .is_some_and(|limit| posting.occurrence_count > limit)
                    || posting.class == PostingClass::Suppressed;
                if skip {
                    result.metrics.postings_skipped =
                        result.metrics.postings_skipped.saturating_add(1);
                } else {
                    pending.push((match_index, *posting));
                }
            }
        }

        let mut position_cache: BTreeMap<u32, Vec<u8>> = BTreeMap::new();
        let mut decoded_cache: BTreeMap<u64, Vec<SeedOccurrence>> = BTreeMap::new();
        let mut checked_pages = BTreeSet::new();
        for (match_index, posting) in pending {
            let page = self
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
            let page_bytes = if let Some(bytes) = position_cache.get(&page.page_id) {
                bytes.as_slice()
            } else {
                let bytes = self.position_page_bytes(*page)?;
                if checked_pages.insert(page.page_id) {
                    result.metrics.position_pages_read =
                        result.metrics.position_pages_read.saturating_add(1);
                    result.metrics.bytes_read = result
                        .metrics
                        .bytes_read
                        .saturating_add(u64::try_from(bytes.len()).unwrap_or(u64::MAX));
                }
                position_cache.insert(page.page_id, bytes);
                position_cache
                    .get(&page.page_id)
                    .expect("inserted position page")
                    .as_slice()
            };
            let occurrences = if let Some(occurrences) = decoded_cache.get(&posting.key_index) {
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
                let payload = page_bytes.get(start..end).ok_or_else(|| {
                    JmaError::CorruptSection(
                        "compact seed posting range is outside its page".to_string(),
                    )
                })?;
                let limit = options.max_occurrences.map(|value| value as usize);
                let occurrences = decode_position_posting(payload, posting, self.header.k, limit)?;
                decoded_cache.insert(posting.key_index, occurrences.clone());
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

    fn filter_contains(&self, hash: u64) -> bool {
        let bytes = &self.bytes[self.filter.clone()];
        let bit_count = bytes.len().saturating_mul(8);
        bit_count != 0
            && filter_bits(hash, bit_count)
                .into_iter()
                .all(|bit| bytes[bit / 8] & (1u8 << (bit % 8)) != 0)
    }

    fn find_prefix(&self, hash: u64) -> JmaResult<Option<SeedPrefix>> {
        let prefix = crate::jma::format::hash_prefix(hash, self.header.bucket_bits)?;
        Ok(self
            .prefixes
            .binary_search_by_key(&prefix, |entry| entry.hash_prefix)
            .ok()
            .map(|index| self.prefixes[index]))
    }

    fn key_page_bytes(&self, page: SeedKeyPage) -> JmaResult<&[u8]> {
        let range = section_range(
            page.data_offset,
            page.data_length,
            self.header.section_length,
            self.bytes.len(),
            "key page",
        )?;
        Ok(&self.bytes[range])
    }

    fn position_page_bytes(&self, page: SeedPositionPage) -> JmaResult<Vec<u8>> {
        let range = section_range(
            page.data_offset,
            page.data_length,
            self.header.section_length,
            self.bytes.len(),
            "position page",
        )?;
        let bytes = self.bytes[range].to_vec();
        if crate::jma::format::checksum(&bytes) != page.checksum {
            return Err(JmaError::ChecksumMismatch(format!(
                "compact seed position page {}",
                page.page_id
            )));
        }
        Ok(bytes)
    }
}

/// Decode the header and all non-payload directories from a prefix that may
/// omit key and position data ranges. This is the remote-reader entry point.
pub fn decode_directory_prefix(
    bytes: &[u8],
    section_length: u64,
) -> JmaResult<SeedDirectoryPrefix> {
    let header = decode_header(bytes, section_length)?;
    validate_header_offsets(&header)?;
    let prefixes = decode_prefixes(bytes, &header)?;
    let key_pages = decode_key_pages(bytes, &header)?;
    let posting_headers = decode_posting_headers(bytes, &header)?;
    let position_pages = decode_position_pages(bytes, &header)?;
    validate_directory(
        &header,
        &prefixes,
        &key_pages,
        &posting_headers,
        &position_pages,
    )?;
    Ok((header, prefixes, key_pages, posting_headers, position_pages))
}

/// Decode only the filter/prefix/key-page directory needed before an exact
/// key match. Posting and position descriptors remain lazy range reads.
pub fn decode_key_lookup_prefix(
    bytes: &[u8],
    section_length: u64,
) -> JmaResult<SeedKeyLookupPrefix> {
    let header = decode_header(bytes, section_length)?;
    validate_header_offsets(&header)?;
    let prefixes = decode_prefixes(bytes, &header)?;
    let key_pages = decode_key_pages(bytes, &header)?;
    validate_key_directory(&header, &prefixes, &key_pages)?;
    Ok((header, prefixes, key_pages))
}

pub fn decode_posting_header_record(
    bytes: &[u8],
    key_index: u64,
    header: &SeedHeader,
) -> JmaResult<SeedPostingHeader> {
    if bytes.len() != POSTING_HEADER_RECORD_SIZE {
        return Err(JmaError::CorruptSection(
            "compact seed posting header record is truncated".to_string(),
        ));
    }
    let mut decoded = decode_posting_headers(
        bytes,
        &SeedHeader {
            posting_header_offset: 0,
            posting_header_count: 1,
            section_length: POSTING_HEADER_RECORD_SIZE as u64,
            ..*header
        },
    )?;
    let mut posting = decoded.pop().ok_or_else(|| {
        JmaError::CorruptSection("compact seed posting header is unavailable".to_string())
    })?;
    if posting.key_index != key_index {
        return Err(JmaError::CorruptSection(
            "compact seed posting header key index mismatch".to_string(),
        ));
    }
    posting.key_index = key_index;
    Ok(posting)
}

pub fn decode_position_page_record(
    bytes: &[u8],
    page_id: u32,
    header: &SeedHeader,
) -> JmaResult<SeedPositionPage> {
    if bytes.len() != POSITION_PAGE_RECORD_SIZE {
        return Err(JmaError::CorruptSection(
            "compact seed position page record is truncated".to_string(),
        ));
    }
    let mut decoded = decode_position_pages(
        bytes,
        &SeedHeader {
            position_page_dir_offset: 0,
            position_page_count: 1,
            section_length: POSITION_PAGE_RECORD_SIZE as u64,
            ..*header
        },
    )?;
    let page = decoded.pop().ok_or_else(|| {
        JmaError::CorruptSection("compact seed position page is unavailable".to_string())
    })?;
    if page.page_id != page_id {
        return Err(JmaError::CorruptSection(
            "compact seed position page id mismatch".to_string(),
        ));
    }
    Ok(page)
}

fn validate_header_offsets(header: &SeedHeader) -> JmaResult<()> {
    let expected_filter =
        u64::try_from(COMPACT_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    if header.filter_offset != expected_filter
        || header.prefix_offset != add(header.filter_offset, header.filter_length)?
        || header.key_page_dir_offset
            != add(
                header.prefix_offset,
                u64::from(header.prefix_count)
                    .checked_mul(PREFIX_RECORD_SIZE as u64)
                    .ok_or(JmaError::OffsetOverflow)?,
            )?
        || header.posting_header_offset
            != add(
                header.key_page_dir_offset,
                u64::from(header.key_page_count)
                    .checked_mul(KEY_PAGE_RECORD_SIZE as u64)
                    .ok_or(JmaError::OffsetOverflow)?,
            )?
        || header.position_page_dir_offset
            != add(
                header.posting_header_offset,
                u64::from(header.posting_header_count)
                    .checked_mul(POSTING_HEADER_RECORD_SIZE as u64)
                    .ok_or(JmaError::OffsetOverflow)?,
            )?
        || header.key_data_offset
            != add(
                header.position_page_dir_offset,
                u64::from(header.position_page_count)
                    .checked_mul(POSITION_PAGE_RECORD_SIZE as u64)
                    .ok_or(JmaError::OffsetOverflow)?,
            )?
        || header.position_data_offset != add(header.key_data_offset, header.key_data_length)?
        || header.section_length != add(header.position_data_offset, header.position_data_length)?
    {
        return Err(JmaError::CorruptSection(
            "compact seed section ranges are inconsistent".to_string(),
        ));
    }
    Ok(())
}

/// Decode only the fixed header before a remote reader requests directories.
pub fn decode_seed_header(bytes: &[u8], section_length: u64) -> JmaResult<SeedHeader> {
    decode_header(bytes, section_length)
}

/// Decode one exact key page and return `(jamhash, packed_kmer)` tuples in
/// sorted order. The checksum remains mandatory for every decoded page.
pub fn decode_key_page_records(
    bytes: &[u8],
    page: SeedKeyPage,
    k: u8,
    encoding: super::KeyEncoding,
) -> JmaResult<Vec<(u64, u64)>> {
    decode_key_page(bytes, page, k, encoding)
}

fn add(left: u64, right: u64) -> JmaResult<u64> {
    left.checked_add(right).ok_or(JmaError::OffsetOverflow)
}

fn decode_prefixes(bytes: &[u8], header: &SeedHeader) -> JmaResult<Vec<SeedPrefix>> {
    let range = section_range(
        header.prefix_offset,
        u64::from(header.prefix_count)
            .checked_mul(PREFIX_RECORD_SIZE as u64)
            .ok_or(JmaError::OffsetOverflow)?,
        header.section_length,
        bytes.len(),
        "prefix directory",
    )?;
    let mut result = Vec::with_capacity(header.prefix_count as usize);
    for offset in (range.start..range.end).step_by(PREFIX_RECORD_SIZE) {
        result.push(SeedPrefix {
            hash_prefix: get_u32(bytes, offset)?,
            first_key: get_u64(bytes, offset + 4)?,
            key_count: get_u32(bytes, offset + 12)?,
            first_page: get_u32(bytes, offset + 16)?,
            page_count: get_u32(bytes, offset + 20)?,
            first_hash: get_u64(bytes, offset + 24)?,
            last_hash: get_u64(bytes, offset + 32)?,
        });
    }
    Ok(result)
}

fn decode_key_pages(bytes: &[u8], header: &SeedHeader) -> JmaResult<Vec<SeedKeyPage>> {
    let range = section_range(
        header.key_page_dir_offset,
        u64::from(header.key_page_count)
            .checked_mul(KEY_PAGE_RECORD_SIZE as u64)
            .ok_or(JmaError::OffsetOverflow)?,
        header.section_length,
        bytes.len(),
        "key page directory",
    )?;
    let mut result = Vec::with_capacity(header.key_page_count as usize);
    for offset in (range.start..range.end).step_by(KEY_PAGE_RECORD_SIZE) {
        if bytes[offset + 16..offset + 20]
            .iter()
            .any(|byte| *byte != 0)
        {
            return Err(JmaError::CorruptSection(
                "compact seed key page reserved bytes are non-zero".to_string(),
            ));
        }
        result.push(SeedKeyPage {
            page_id: get_u32(bytes, offset)?,
            first_key: get_u64(bytes, offset + 4)?,
            key_count: get_u32(bytes, offset + 12)?,
            first_hash: get_u64(bytes, offset + 20)?,
            last_hash: get_u64(bytes, offset + 28)?,
            data_offset: get_u64(bytes, offset + 36)?,
            data_length: get_u64(bytes, offset + 44)?,
            checksum: get_array(bytes, offset + 52)?,
        });
    }
    Ok(result)
}

fn decode_posting_headers(bytes: &[u8], header: &SeedHeader) -> JmaResult<Vec<SeedPostingHeader>> {
    let range = section_range(
        header.posting_header_offset,
        u64::from(header.posting_header_count)
            .checked_mul(POSTING_HEADER_RECORD_SIZE as u64)
            .ok_or(JmaError::OffsetOverflow)?,
        header.section_length,
        bytes.len(),
        "posting header directory",
    )?;
    let mut result = Vec::with_capacity(header.posting_header_count as usize);
    for offset in (range.start..range.end).step_by(POSTING_HEADER_RECORD_SIZE) {
        let class = PostingClass::from_byte(get_u8(bytes, offset + 12)?).ok_or_else(|| {
            JmaError::CorruptSection("compact seed posting class is invalid".to_string())
        })?;
        let encoding =
            OccurrenceEncoding::from_byte(get_u8(bytes, offset + 13)?).ok_or_else(|| {
                JmaError::CorruptSection("compact seed posting encoding is invalid".to_string())
            })?;
        if get_u32(bytes, offset + 28)? != 0 {
            return Err(JmaError::CorruptSection(
                "compact seed posting reserved bytes are non-zero".to_string(),
            ));
        }
        result.push(SeedPostingHeader {
            key_index: get_u64(bytes, offset)?,
            occurrence_count: get_u32(bytes, offset + 8)?,
            class,
            encoding,
            flags: get_u16(bytes, offset + 14)?,
            position_page_id: get_u32(bytes, offset + 16)?,
            position_offset: get_u32(bytes, offset + 20)?,
            position_length: get_u32(bytes, offset + 24)?,
            checksum: get_array(bytes, offset + 32)?,
        });
    }
    Ok(result)
}

fn decode_position_pages(bytes: &[u8], header: &SeedHeader) -> JmaResult<Vec<SeedPositionPage>> {
    let range = section_range(
        header.position_page_dir_offset,
        u64::from(header.position_page_count)
            .checked_mul(POSITION_PAGE_RECORD_SIZE as u64)
            .ok_or(JmaError::OffsetOverflow)?,
        header.section_length,
        bytes.len(),
        "position page directory",
    )?;
    let mut result = Vec::with_capacity(header.position_page_count as usize);
    for offset in (range.start..range.end).step_by(POSITION_PAGE_RECORD_SIZE) {
        if get_u32(bytes, offset + 12)? != 0 {
            return Err(JmaError::CorruptSection(
                "compact seed position page reserved bytes are non-zero".to_string(),
            ));
        }
        result.push(SeedPositionPage {
            page_id: get_u32(bytes, offset)?,
            first_header: get_u32(bytes, offset + 4)?,
            header_count: get_u32(bytes, offset + 8)?,
            data_offset: get_u64(bytes, offset + 16)?,
            data_length: get_u64(bytes, offset + 24)?,
            checksum: get_array(bytes, offset + 32)?,
        });
    }
    Ok(result)
}

fn validate_directory(
    header: &SeedHeader,
    prefixes: &[SeedPrefix],
    key_pages: &[SeedKeyPage],
    postings: &[SeedPostingHeader],
    position_pages: &[SeedPositionPage],
) -> JmaResult<()> {
    validate_key_directory(header, prefixes, key_pages)?;
    if position_pages
        .iter()
        .enumerate()
        .any(|(index, page)| page.page_id != u32::try_from(index).unwrap_or(u32::MAX))
    {
        return Err(JmaError::CorruptSection(
            "compact seed directories are not sorted".to_string(),
        ));
    }
    for (index, posting) in postings.iter().enumerate() {
        if posting.key_index != u64::try_from(index).map_err(|_| JmaError::OffsetOverflow)?
            || posting.occurrence_count == 0
            || posting.encoding != header.occurrence_encoding
            || posting.flags != POSTING_FLAG_POSITION_BEARING
            || usize::try_from(posting.position_page_id)
                .ok()
                .is_none_or(|page| page >= position_pages.len())
        {
            return Err(JmaError::CorruptSection(
                "compact seed posting header is invalid".to_string(),
            ));
        }
    }
    let mut expected_header = 0u32;
    for page in position_pages {
        if page.first_header != expected_header || page.header_count == 0 {
            return Err(JmaError::CorruptSection(
                "compact seed position page header range is invalid".to_string(),
            ));
        }
        expected_header = expected_header
            .checked_add(page.header_count)
            .ok_or(JmaError::OffsetOverflow)?;
        let page_end = page
            .data_offset
            .checked_add(page.data_length)
            .ok_or(JmaError::OffsetOverflow)?;
        if page.data_offset < header.position_data_offset
            || page_end > header.position_data_offset + header.position_data_length
        {
            return Err(JmaError::CorruptSection(
                "compact seed position page exceeds payload".to_string(),
            ));
        }
        for posting in postings
            .iter()
            .skip(page.first_header as usize)
            .take(page.header_count as usize)
        {
            if posting.position_page_id != page.page_id {
                return Err(JmaError::CorruptSection(
                    "compact seed posting points to the wrong position page".to_string(),
                ));
            }
            let end = u64::from(posting.position_offset)
                .checked_add(u64::from(posting.position_length))
                .ok_or(JmaError::OffsetOverflow)?;
            if end > page.data_length {
                return Err(JmaError::CorruptSection(
                    "compact seed posting exceeds position page".to_string(),
                ));
            }
        }
    }
    if expected_header != u32::try_from(postings.len()).map_err(|_| JmaError::OffsetOverflow)? {
        return Err(JmaError::CorruptSection(
            "compact seed position pages do not cover postings".to_string(),
        ));
    }
    Ok(())
}

fn validate_key_directory(
    header: &SeedHeader,
    prefixes: &[SeedPrefix],
    key_pages: &[SeedKeyPage],
) -> JmaResult<()> {
    if prefixes
        .windows(2)
        .any(|pair| pair[0].hash_prefix >= pair[1].hash_prefix)
        || key_pages
            .iter()
            .enumerate()
            .any(|(index, page)| page.page_id != u32::try_from(index).unwrap_or(u32::MAX))
    {
        return Err(JmaError::CorruptSection(
            "compact seed directories are not sorted".to_string(),
        ));
    }
    let mut key_total = 0u64;
    for prefix in prefixes {
        if prefix.key_count == 0
            || prefix.first_key != key_total
            || prefix.first_page >= header.key_page_count
            || prefix.page_count == 0
            || prefix
                .first_page
                .checked_add(prefix.page_count)
                .is_none_or(|end| end > header.key_page_count)
            || prefix.first_hash > prefix.last_hash
        {
            return Err(JmaError::CorruptSection(
                "compact seed prefix range is invalid".to_string(),
            ));
        }
        key_total = key_total
            .checked_add(u64::from(prefix.key_count))
            .ok_or(JmaError::OffsetOverflow)?;
    }
    if key_total != header.key_count {
        return Err(JmaError::CorruptSection(
            "compact seed prefix count does not match keys".to_string(),
        ));
    }
    let key_record_size =
        u64::try_from(header.key_encoding.record_size()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut expected_key = 0u64;
    for page in key_pages {
        let page_end = page
            .data_offset
            .checked_add(page.data_length)
            .ok_or(JmaError::OffsetOverflow)?;
        if page.first_key != expected_key
            || page.key_count == 0
            || page.data_length != u64::from(page.key_count) * key_record_size
            || page.data_offset < header.key_data_offset
            || page_end > header.key_data_offset + header.key_data_length
        {
            return Err(JmaError::CorruptSection(
                "compact seed key page range is invalid".to_string(),
            ));
        }
        expected_key = expected_key
            .checked_add(u64::from(page.key_count))
            .ok_or(JmaError::OffsetOverflow)?;
    }
    if expected_key != header.key_count {
        return Err(JmaError::CorruptSection(
            "compact seed key pages do not cover all keys".to_string(),
        ));
    }
    Ok(())
}

fn decode_key_page(
    bytes: &[u8],
    page: SeedKeyPage,
    k: u8,
    encoding: super::KeyEncoding,
) -> JmaResult<Vec<(u64, u64)>> {
    if crate::jma::format::checksum(bytes) != page.checksum {
        return Err(JmaError::ChecksumMismatch(format!(
            "compact seed key page {}",
            page.page_id
        )));
    }
    let size = encoding.record_size();
    if bytes.len() != page.key_count as usize * size {
        return Err(JmaError::CorruptSection(
            "compact seed key page length is invalid".to_string(),
        ));
    }
    let mut keys = Vec::with_capacity(page.key_count as usize);
    for index in 0..page.key_count as usize {
        let offset = index * size;
        let packed = match encoding {
            super::KeyEncoding::Packed8 => get_u64(bytes, offset)?,
            super::KeyEncoding::Packed6 => {
                let mut value = [0u8; 8];
                value[..6].copy_from_slice(&bytes[offset..offset + 6]);
                u64::from_le_bytes(value)
            }
        };
        if !crate::jma::format::valid_canonical_kmer(k, packed) {
            return Err(JmaError::CorruptSection(
                "compact seed key has invalid packed k-mer".to_string(),
            ));
        }
        let hash = crate::jamhash_u64_v1(packed);
        if hash == 0 || hash < page.first_hash || hash > page.last_hash {
            return Err(JmaError::CorruptSection(
                "compact seed key is outside page hash range".to_string(),
            ));
        }
        keys.push((hash, packed));
    }
    if keys.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(JmaError::CorruptSection(
            "compact seed key page is not sorted".to_string(),
        ));
    }
    Ok(keys)
}

/// Decode one posting without allocating beyond `limit`. The complete byte
/// payload is checksum-verified even when the caller caps returned positions.
pub fn decode_position_posting(
    bytes: &[u8],
    header: SeedPostingHeader,
    k: u8,
    limit: Option<usize>,
) -> JmaResult<Vec<SeedOccurrence>> {
    if bytes.len()
        != usize::try_from(header.position_length).map_err(|_| JmaError::OffsetOverflow)?
        || crate::jma::format::checksum(bytes) != header.checksum
    {
        return Err(JmaError::ChecksumMismatch(
            "compact seed posting".to_string(),
        ));
    }
    let output_limit = limit.unwrap_or(header.occurrence_count as usize);
    let mut occurrences: Vec<SeedOccurrence> =
        Vec::with_capacity(output_limit.min(header.occurrence_count as usize));
    match header.encoding {
        OccurrenceEncoding::Fixed16 => {
            let expected = (header.occurrence_count as usize)
                .checked_mul(16)
                .ok_or(JmaError::OffsetOverflow)?;
            if bytes.len() != expected {
                return Err(JmaError::CorruptSection(
                    "compact fixed posting length is invalid".to_string(),
                ));
            }
            let mut previous = None;
            for index in 0..header.occurrence_count as usize {
                let offset = index * 16;
                let contig_id = get_u32(bytes, offset)?;
                let reverse = match get_u8(bytes, offset + 4)? {
                    0 => false,
                    1 => true,
                    _ => {
                        return Err(JmaError::CorruptSection(
                            "compact seed orientation is invalid".to_string(),
                        ));
                    }
                };
                if bytes[offset + 5..offset + 8].iter().any(|byte| *byte != 0) {
                    return Err(JmaError::CorruptSection(
                        "compact fixed posting reserved bytes are non-zero".to_string(),
                    ));
                }
                let occurrence = SeedOccurrence {
                    contig_id,
                    position: get_u64(bytes, offset + 8)?,
                    reverse,
                };
                validate_occurrence(occurrence, k)?;
                if previous.is_some_and(|previous: SeedOccurrence| {
                    (previous.contig_id, previous.position, previous.reverse)
                        > (
                            occurrence.contig_id,
                            occurrence.position,
                            occurrence.reverse,
                        )
                }) {
                    return Err(JmaError::CorruptSection(
                        "compact seed posting is not sorted".to_string(),
                    ));
                }
                previous = Some(occurrence);
                if occurrences.len() < output_limit {
                    occurrences.push(occurrence);
                }
            }
        }
        OccurrenceEncoding::DeltaVarint => {
            let mut cursor = 0usize;
            let mut previous_contig = 0u32;
            let mut previous_position = 0u64;
            for index in 0..header.occurrence_count as usize {
                let contig_delta = read_varint(bytes, &mut cursor)?;
                let contig = previous_contig
                    .checked_add(u32::try_from(contig_delta).map_err(|_| JmaError::OffsetOverflow)?)
                    .ok_or(JmaError::OffsetOverflow)?;
                let position_value = read_varint(bytes, &mut cursor)?;
                let position = if contig == previous_contig {
                    previous_position
                        .checked_add(position_value)
                        .ok_or(JmaError::OffsetOverflow)?
                } else {
                    position_value
                };
                let reverse = match get_u8(bytes, cursor)? {
                    0 => false,
                    1 => true,
                    _ => {
                        return Err(JmaError::CorruptSection(
                            "compact seed orientation is invalid".to_string(),
                        ));
                    }
                };
                cursor = cursor.checked_add(1).ok_or(JmaError::OffsetOverflow)?;
                let occurrence = SeedOccurrence {
                    contig_id: contig,
                    position,
                    reverse,
                };
                validate_occurrence(occurrence, k)?;
                if index > 0
                    && (contig < previous_contig
                        || (contig == previous_contig && position < previous_position))
                {
                    return Err(JmaError::CorruptSection(
                        "compact seed posting is not sorted".to_string(),
                    ));
                }
                if occurrences.len() < output_limit {
                    occurrences.push(occurrence);
                }
                previous_contig = contig;
                previous_position = position;
            }
            if cursor != bytes.len() {
                return Err(JmaError::CorruptSection(
                    "compact delta posting has trailing bytes".to_string(),
                ));
            }
        }
    }
    Ok(occurrences)
}

fn read_varint(bytes: &[u8], cursor: &mut usize) -> JmaResult<u64> {
    let mut value = 0u64;
    for shift in (0..64).step_by(7) {
        let byte = *bytes.get(*cursor).ok_or_else(|| {
            JmaError::CorruptSection("compact delta posting is truncated".to_string())
        })?;
        *cursor = cursor.checked_add(1).ok_or(JmaError::OffsetOverflow)?;
        let payload = u64::from(byte & 0x7f);
        value |= payload.checked_shl(shift).ok_or(JmaError::OffsetOverflow)?;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(JmaError::CorruptSection(
        "compact delta posting varint is too long".to_string(),
    ))
}
