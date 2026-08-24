use super::format::{
    COMPACT_HEADER_SIZE, KEY_PAGE_RECORD_SIZE, POSITION_PAGE_RECORD_SIZE,
    POSTING_FLAG_POSITION_BEARING, POSTING_HEADER_RECORD_SIZE, PREFIX_RECORD_SIZE, SeedHeader,
    SeedKeyPage, SeedPositionPage, SeedPostingHeader, SeedPrefix, encode_header,
};
use super::{
    JAMHASH_ALGORITHM_ID, KeyEncoding, OccurrenceEncoding, PostingClass, SeedBinding,
    validate_occurrence, validate_query,
};
use crate::jma::format::hash_prefix;
use crate::jma::writer::SeedRecord;
use crate::jma::{JmaError, JmaResult};
use std::collections::{BTreeMap, BTreeSet};

/// Build controls for the compact JMA seed section.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedBuildConfig {
    pub key_encoding: KeyEncoding,
    pub occurrence_encoding: OccurrenceEncoding,
    pub bucket_bits: u8,
    pub key_page_bytes: u32,
    pub position_page_bytes: u32,
    pub filter_bits_per_key: u8,
    pub rare_occurrence_limit: u32,
    pub common_occurrence_limit: u32,
    pub suppressed_occurrence_limit: u32,
}

impl Default for SeedBuildConfig {
    fn default() -> Self {
        Self {
            key_encoding: KeyEncoding::Packed8,
            occurrence_encoding: OccurrenceEncoding::DeltaVarint,
            bucket_bits: 12,
            key_page_bytes: 4096,
            position_page_bytes: 64 * 1024,
            filter_bits_per_key: 10,
            rare_occurrence_limit: 3,
            common_occurrence_limit: 64,
            suppressed_occurrence_limit: 4096,
        }
    }
}

impl SeedBuildConfig {
    fn validate(self, binding: SeedBinding) -> JmaResult<()> {
        binding.validate()?;
        if !(1..=24).contains(&self.bucket_bits) {
            return Err(JmaError::CorruptSection(
                "compact seed bucket width must be between 1 and 24".to_string(),
            ));
        }
        if self.key_page_bytes < self.key_encoding.record_size() as u32
            || self.position_page_bytes == 0
            || self.filter_bits_per_key < 2
            || self.rare_occurrence_limit > self.common_occurrence_limit
            || self.common_occurrence_limit > self.suppressed_occurrence_limit
        {
            return Err(JmaError::CorruptSection(
                "compact seed page or occurrence limits are invalid".to_string(),
            ));
        }
        if self.key_encoding == KeyEncoding::Packed6 && binding.k > 24 {
            return Err(JmaError::CorruptSection(
                "six-byte packed keys support k at most 24".to_string(),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct KeyData {
    hash: u64,
    packed: u64,
    occurrences: Vec<(u32, u64, bool)>,
}

type ExactKey = (u64, u64);
type OccurrenceKey = (u32, u64, bool);

#[derive(Clone, Debug)]
struct PositionPageData {
    first_header: u32,
    bytes: Vec<u8>,
}

/// Encode one immutable compact seed section. Records are deduplicated and
/// sorted by exact key, then by contig, position, and orientation.
pub fn encode_seed_section(
    records: &[SeedRecord],
    binding: SeedBinding,
    config: SeedBuildConfig,
) -> JmaResult<Vec<u8>> {
    config.validate(binding)?;
    let mut grouped: BTreeMap<ExactKey, BTreeSet<OccurrenceKey>> = BTreeMap::new();
    for record in records {
        validate_query(record.query, binding.k)?;
        validate_occurrence(record.occurrence, binding.k)?;
        grouped
            .entry((record.query.hash, record.query.canonical_kmer))
            .or_default()
            .insert((
                record.occurrence.contig_id,
                record.occurrence.position,
                record.occurrence.reverse,
            ));
    }
    let keys = grouped
        .into_iter()
        .map(|((hash, packed), occurrences)| KeyData {
            hash,
            packed,
            occurrences: occurrences.into_iter().collect(),
        })
        .collect::<Vec<_>>();
    let key_count = u64::try_from(keys.len()).map_err(|_| JmaError::OffsetOverflow)?;
    if key_count > u64::from(u32::MAX) {
        return Err(JmaError::OffsetOverflow);
    }

    let filter = build_filter(&keys, config.filter_bits_per_key)?;
    let mut key_data = Vec::new();
    let mut key_pages = Vec::new();
    let key_record_size = config.key_encoding.record_size();
    let page_key_limit = usize::try_from(config.key_page_bytes)
        .map_err(|_| JmaError::OffsetOverflow)?
        / key_record_size;
    let page_key_limit = page_key_limit.max(1);
    for (index, chunk) in keys.chunks(page_key_limit).enumerate() {
        let first_key = u64::try_from(index)
            .map_err(|_| JmaError::OffsetOverflow)?
            .checked_mul(u64::try_from(page_key_limit).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let local_offset = u64::try_from(key_data.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let mut page_bytes = Vec::with_capacity(chunk.len() * key_record_size);
        for key in chunk {
            encode_packed_key(&mut page_bytes, key.packed, config.key_encoding);
        }
        let first_hash = chunk.first().map_or(0, |key| key.hash);
        let last_hash = chunk.last().map_or(0, |key| key.hash);
        let checksum = crate::jma::format::checksum(&page_bytes);
        key_data.extend_from_slice(&page_bytes);
        key_pages.push(SeedKeyPage {
            page_id: u32::try_from(index).map_err(|_| JmaError::OffsetOverflow)?,
            first_key,
            key_count: u32::try_from(chunk.len()).map_err(|_| JmaError::OffsetOverflow)?,
            first_hash,
            last_hash,
            data_offset: local_offset,
            data_length: u64::try_from(page_bytes.len()).map_err(|_| JmaError::OffsetOverflow)?,
            checksum,
        });
    }
    let prefixes = build_prefixes(&keys, &key_pages, config.bucket_bits)?;

    let mut posting_headers = Vec::with_capacity(keys.len());
    let mut position_pages = Vec::new();
    let position_page_limit =
        usize::try_from(config.position_page_bytes).map_err(|_| JmaError::OffsetOverflow)?;
    let mut current_page = Vec::new();
    let mut current_first_header = 0u32;
    for (key_index, key) in keys.iter().enumerate() {
        let payload = encode_occurrences(&key.occurrences, config.occurrence_encoding)?;
        if !current_page.is_empty()
            && current_page.len().saturating_add(payload.len()) > position_page_limit
        {
            position_pages.push(PositionPageData {
                first_header: current_first_header,
                bytes: std::mem::take(&mut current_page),
            });
            current_first_header =
                u32::try_from(key_index).map_err(|_| JmaError::OffsetOverflow)?;
        }
        let page_id = u32::try_from(position_pages.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let position_offset =
            u32::try_from(current_page.len()).map_err(|_| JmaError::OffsetOverflow)?;
        current_page.extend_from_slice(&payload);
        let class = classify(key.occurrences.len(), config);
        posting_headers.push(SeedPostingHeader {
            key_index: u64::try_from(key_index).map_err(|_| JmaError::OffsetOverflow)?,
            occurrence_count: u32::try_from(key.occurrences.len())
                .map_err(|_| JmaError::OffsetOverflow)?,
            class,
            encoding: config.occurrence_encoding,
            flags: POSTING_FLAG_POSITION_BEARING,
            position_page_id: page_id,
            position_offset,
            position_length: u32::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?,
            checksum: crate::jma::format::checksum(&payload),
        });
        if current_page.len() >= position_page_limit {
            position_pages.push(PositionPageData {
                first_header: current_first_header,
                bytes: std::mem::take(&mut current_page),
            });
            current_first_header =
                u32::try_from(key_index + 1).map_err(|_| JmaError::OffsetOverflow)?;
        }
    }
    if !current_page.is_empty() {
        position_pages.push(PositionPageData {
            first_header: current_first_header,
            bytes: current_page,
        });
    }

    let filter_offset = u64::try_from(COMPACT_HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    let prefix_offset = add_offset(filter_offset, filter.len())?;
    let key_page_dir_offset = add_offset(
        prefix_offset,
        prefixes
            .len()
            .checked_mul(PREFIX_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let posting_header_offset = add_offset(
        key_page_dir_offset,
        key_pages
            .len()
            .checked_mul(KEY_PAGE_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let position_page_dir_offset = add_offset(
        posting_header_offset,
        posting_headers
            .len()
            .checked_mul(POSTING_HEADER_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let key_data_offset = add_offset(
        position_page_dir_offset,
        position_pages
            .len()
            .checked_mul(POSITION_PAGE_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let position_data_offset = add_offset(key_data_offset, key_data.len())?;
    let section_length = add_offset(
        position_data_offset,
        position_pages
            .iter()
            .map(|page| page.bytes.len())
            .sum::<usize>(),
    )?;

    for page in &mut key_pages {
        page.data_offset = key_data_offset
            .checked_add(page.data_offset)
            .ok_or(JmaError::OffsetOverflow)?;
    }
    let mut position_cursor = position_data_offset;
    let mut position_page_records = Vec::with_capacity(position_pages.len());
    for (page_id, page) in position_pages.iter().enumerate() {
        let data_offset = position_cursor;
        let data_length = u64::try_from(page.bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
        position_page_records.push(SeedPositionPage {
            page_id: u32::try_from(page_id).map_err(|_| JmaError::OffsetOverflow)?,
            first_header: page.first_header,
            header_count: position_pages
                .get(page_id + 1)
                .map_or(
                    u32::try_from(posting_headers.len()).unwrap_or(u32::MAX),
                    |next| next.first_header,
                )
                .saturating_sub(page.first_header),
            data_offset,
            data_length,
            checksum: crate::jma::format::checksum(&page.bytes),
        });
        position_cursor = position_cursor
            .checked_add(data_length)
            .ok_or(JmaError::OffsetOverflow)?;
    }
    let mut header = SeedHeader {
        scheme_id: binding.scheme_id,
        k: binding.k,
        key_encoding: config.key_encoding,
        occurrence_encoding: config.occurrence_encoding,
        bucket_bits: config.bucket_bits,
        scale: binding.scale,
        hash_algorithm_id: JAMHASH_ALGORITHM_ID,
        key_count,
        posting_count: key_count,
        filter_offset,
        filter_length: u64::try_from(filter.len()).map_err(|_| JmaError::OffsetOverflow)?,
        prefix_offset,
        prefix_count: u32::try_from(prefixes.len()).map_err(|_| JmaError::OffsetOverflow)?,
        key_page_dir_offset,
        key_page_count: u32::try_from(key_pages.len()).map_err(|_| JmaError::OffsetOverflow)?,
        posting_header_offset,
        posting_header_count: u32::try_from(posting_headers.len())
            .map_err(|_| JmaError::OffsetOverflow)?,
        position_page_dir_offset,
        position_page_count: u32::try_from(position_page_records.len())
            .map_err(|_| JmaError::OffsetOverflow)?,
        key_data_offset,
        key_data_length: u64::try_from(key_data.len()).map_err(|_| JmaError::OffsetOverflow)?,
        position_data_offset,
        position_data_length: position_cursor - position_data_offset,
        section_length,
        flags: 0,
        catalog_checksum: binding.catalog_checksum,
        archive_checksum: binding.archive_checksum,
        scheme_checksum: binding.scheme_checksum,
        section_checksum: [0; 32],
    };

    let mut output = encode_header(header);
    output.extend_from_slice(&filter);
    for prefix in &prefixes {
        encode_prefix(&mut output, *prefix);
    }
    for page in &key_pages {
        encode_key_page(&mut output, *page);
    }
    for posting in &posting_headers {
        encode_posting_header(&mut output, *posting);
    }
    for page in &position_page_records {
        encode_position_page(&mut output, *page);
    }
    output.extend_from_slice(&key_data);
    for page in &position_pages {
        output.extend_from_slice(&page.bytes);
    }
    if u64::try_from(output.len()).map_err(|_| JmaError::OffsetOverflow)? != section_length {
        return Err(JmaError::CorruptSection(
            "compact seed section length is inconsistent".to_string(),
        ));
    }
    header.section_checksum = crate::jma::format::checksum(&output);
    let encoded_header = encode_header(header);
    output[..COMPACT_HEADER_SIZE].copy_from_slice(&encoded_header);
    Ok(output)
}

fn add_offset(offset: u64, length: usize) -> JmaResult<u64> {
    offset
        .checked_add(u64::try_from(length).map_err(|_| JmaError::OffsetOverflow)?)
        .ok_or(JmaError::OffsetOverflow)
}

fn encode_packed_key(bytes: &mut Vec<u8>, packed: u64, encoding: KeyEncoding) {
    match encoding {
        KeyEncoding::Packed8 => bytes.extend_from_slice(&packed.to_le_bytes()),
        KeyEncoding::Packed6 => bytes.extend_from_slice(&packed.to_le_bytes()[..6]),
    }
}

fn build_filter(keys: &[KeyData], bits_per_key: u8) -> JmaResult<Vec<u8>> {
    let requested = keys.len().saturating_mul(usize::from(bits_per_key)).max(64);
    let bit_count = requested
        .checked_next_power_of_two()
        .ok_or(JmaError::OffsetOverflow)?;
    let byte_count = bit_count.checked_div(8).ok_or(JmaError::OffsetOverflow)?;
    let mut filter = vec![0u8; byte_count];
    for key in keys {
        for bit in filter_bits(key.hash, bit_count) {
            filter[bit / 8] |= 1u8 << (bit % 8);
        }
    }
    Ok(filter)
}

pub(crate) fn filter_bits(hash: u64, bit_count: usize) -> [usize; 3] {
    let mixed = hash.wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(29) ^ hash.rotate_right(17);
    [
        usize::try_from(hash).unwrap_or(0) % bit_count,
        usize::try_from(mixed).unwrap_or(0) % bit_count,
        usize::try_from(mixed.rotate_left(23)).unwrap_or(0) % bit_count,
    ]
}

fn build_prefixes(
    keys: &[KeyData],
    pages: &[SeedKeyPage],
    bucket_bits: u8,
) -> JmaResult<Vec<SeedPrefix>> {
    let mut result = Vec::new();
    let mut start = 0usize;
    while start < keys.len() {
        let prefix = hash_prefix(keys[start].hash, bucket_bits)?;
        let mut end = start + 1;
        while end < keys.len() && hash_prefix(keys[end].hash, bucket_bits)? == prefix {
            end += 1;
        }
        let first_page = page_for_key(pages, start)?;
        let last_page = page_for_key(pages, end - 1)?;
        result.push(SeedPrefix {
            hash_prefix: prefix,
            first_key: u64::try_from(start).map_err(|_| JmaError::OffsetOverflow)?,
            key_count: u32::try_from(end - start).map_err(|_| JmaError::OffsetOverflow)?,
            first_page,
            page_count: last_page
                .checked_sub(first_page)
                .and_then(|value| value.checked_add(1))
                .ok_or(JmaError::OffsetOverflow)?,
            first_hash: keys[start].hash,
            last_hash: keys[end - 1].hash,
        });
        start = end;
    }
    Ok(result)
}

fn page_for_key(pages: &[SeedKeyPage], key: usize) -> JmaResult<u32> {
    pages
        .iter()
        .find(|page| {
            let first = usize::try_from(page.first_key).unwrap_or(usize::MAX);
            let end = first.saturating_add(page.key_count as usize);
            (first..end).contains(&key)
        })
        .map(|page| page.page_id)
        .ok_or_else(|| JmaError::CorruptSection("compact seed key has no page".to_string()))
}

fn classify(count: usize, config: SeedBuildConfig) -> PostingClass {
    let count = u32::try_from(count).unwrap_or(u32::MAX);
    if count > config.suppressed_occurrence_limit {
        PostingClass::Suppressed
    } else if count > config.common_occurrence_limit {
        PostingClass::CollectionCommon
    } else if count > config.rare_occurrence_limit {
        PostingClass::ModeratelyCommon
    } else {
        PostingClass::Rare
    }
}

fn encode_occurrences(
    occurrences: &[(u32, u64, bool)],
    encoding: OccurrenceEncoding,
) -> JmaResult<Vec<u8>> {
    if occurrences.is_empty() {
        return Err(JmaError::CorruptSection(
            "compact seed posting has no occurrences".to_string(),
        ));
    }
    for pair in occurrences.windows(2) {
        if pair[1] < pair[0] {
            return Err(JmaError::CorruptSection(
                "compact seed occurrences are not sorted".to_string(),
            ));
        }
    }
    let mut bytes = Vec::new();
    match encoding {
        OccurrenceEncoding::Fixed16 => {
            for &(contig, position, reverse) in occurrences {
                bytes.extend_from_slice(&contig.to_le_bytes());
                bytes.push(u8::from(reverse));
                bytes.extend_from_slice(&[0, 0, 0]);
                bytes.extend_from_slice(&position.to_le_bytes());
            }
        }
        OccurrenceEncoding::DeltaVarint => {
            let mut previous_contig = 0u32;
            let mut previous_position = 0u64;
            for &(contig, position, reverse) in occurrences {
                let contig_delta = contig.checked_sub(previous_contig).ok_or_else(|| {
                    JmaError::CorruptSection("compact seed contig order is invalid".to_string())
                })?;
                put_varint(&mut bytes, u64::from(contig_delta));
                let position_delta = if contig == previous_contig {
                    position.checked_sub(previous_position).ok_or_else(|| {
                        JmaError::CorruptSection(
                            "compact seed position order is invalid".to_string(),
                        )
                    })?
                } else {
                    position
                };
                put_varint(&mut bytes, position_delta);
                bytes.push(u8::from(reverse));
                previous_contig = contig;
                previous_position = position;
            }
        }
    }
    Ok(bytes)
}

pub(crate) fn put_varint(bytes: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        bytes.push((value as u8) | 0x80);
        value >>= 7;
    }
    bytes.push(value as u8);
}

fn encode_prefix(bytes: &mut Vec<u8>, prefix: SeedPrefix) {
    bytes.extend_from_slice(&prefix.hash_prefix.to_le_bytes());
    bytes.extend_from_slice(&prefix.first_key.to_le_bytes());
    bytes.extend_from_slice(&prefix.key_count.to_le_bytes());
    bytes.extend_from_slice(&prefix.first_page.to_le_bytes());
    bytes.extend_from_slice(&prefix.page_count.to_le_bytes());
    bytes.extend_from_slice(&prefix.first_hash.to_le_bytes());
    bytes.extend_from_slice(&prefix.last_hash.to_le_bytes());
}

fn encode_key_page(bytes: &mut Vec<u8>, page: SeedKeyPage) {
    bytes.extend_from_slice(&page.page_id.to_le_bytes());
    bytes.extend_from_slice(&page.first_key.to_le_bytes());
    bytes.extend_from_slice(&page.key_count.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&page.first_hash.to_le_bytes());
    bytes.extend_from_slice(&page.last_hash.to_le_bytes());
    bytes.extend_from_slice(&page.data_offset.to_le_bytes());
    bytes.extend_from_slice(&page.data_length.to_le_bytes());
    bytes.extend_from_slice(&page.checksum);
}

fn encode_posting_header(bytes: &mut Vec<u8>, header: SeedPostingHeader) {
    bytes.extend_from_slice(&header.key_index.to_le_bytes());
    bytes.extend_from_slice(&header.occurrence_count.to_le_bytes());
    bytes.push(header.class as u8);
    bytes.push(header.encoding as u8);
    bytes.extend_from_slice(&header.flags.to_le_bytes());
    bytes.extend_from_slice(&header.position_page_id.to_le_bytes());
    bytes.extend_from_slice(&header.position_offset.to_le_bytes());
    bytes.extend_from_slice(&header.position_length.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&header.checksum);
}

fn encode_position_page(bytes: &mut Vec<u8>, page: SeedPositionPage) {
    bytes.extend_from_slice(&page.page_id.to_le_bytes());
    bytes.extend_from_slice(&page.first_header.to_le_bytes());
    bytes.extend_from_slice(&page.header_count.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&page.data_offset.to_le_bytes());
    bytes.extend_from_slice(&page.data_length.to_le_bytes());
    bytes.extend_from_slice(&page.checksum);
}
