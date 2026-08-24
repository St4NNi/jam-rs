//! Embedded binary lookup directories for JMA format 1.
//!
//! This module deliberately contains no text or external-index format. The
//! sequence directory and seed page directories are sections inside the one
//! JMA object and point at independently checksummed payload ranges.

use crate::archive::SeedSchemeDescriptor;
use crate::jma::format::{
    SECTION_VERSION, checked_end, checksum, get_u32, get_u64, hash_prefix, valid_canonical_kmer,
};
use crate::jma::writer::{SeedRecord, SeedSection};
use crate::jma::{ContigMetadata, JmaError, JmaResult};
use crate::sequence::{
    DEFAULT_MAX_DECODED_BLOCK_BYTES, EncodedSequenceBlock, SequenceBlockDirectory,
    SequenceBlockRecord, decode_stored_block_bounded,
};

pub const SEED_SCHEME_RECORD_SIZE: usize = 64;
pub const SEED_PAGE_RECORD_SIZE: usize = 104;
pub const SEED_SCHEME_RECORD_SIZE_U64: u64 = 64;
pub const SEED_PAGE_RECORD_SIZE_U64: u64 = 104;
pub const SEED_BUCKET_BITS: u8 = 12;
pub const MAX_SEED_PAGE_RECORDS: usize = 4096;
pub const SEED_DIRECTORY_HEADER_SIZE: usize = 40;
pub const SEED_KEY_RECORD_SIZE: usize = 24;
pub const SEED_OCCURRENCE_RECORD_SIZE: usize = 24;
pub const SEED_DIRECTORY_HEADER_SIZE_U64: u64 = 40;
pub const SEED_KEY_RECORD_SIZE_U64: u64 = 24;
pub const SEED_OCCURRENCE_RECORD_SIZE_U64: u64 = 24;

/// A seed page has separate fixed-width key and occurrence ranges. A digest
/// match remains only a candidate until the packed canonical k-mer is checked.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedPageIndex {
    pub scheme_id: u32,
    pub page_id: u32,
    pub hash_prefix: u32,
    pub key_count: u32,
    pub occurrence_count: u32,
    pub key_offset: u64,
    pub key_length: u64,
    pub occurrence_offset: u64,
    pub occurrence_length: u64,
    pub first_hash: u64,
    pub last_hash: u64,
    pub checksum: [u8; 32],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedSchemeIndex {
    pub descriptor: SeedSchemeDescriptor,
    pub pages: Vec<SeedPageIndex>,
    pub record_count: u64,
    declared_page_count: u32,
    declared_page_first: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedIndexDirectory {
    pub schemes: Vec<SeedSchemeIndex>,
}

/// Returns the byte length of the seed directory prefix (header, scheme
/// records, and page records), excluding key and occurrence payload pages.
pub fn seed_directory_prefix_length(bytes: &[u8]) -> JmaResult<usize> {
    if bytes.len() < SEED_DIRECTORY_HEADER_SIZE {
        return Err(JmaError::CorruptSection(
            "seed directory header is truncated".to_string(),
        ));
    }
    if get_u32(bytes, 0)? != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported seed directory version".to_string(),
        ));
    }
    let scheme_count = get_u32(bytes, 4)?;
    let page_count = get_u32(bytes, 8)?;
    let scheme_end = checked_end(
        get_u64(bytes, 16)?,
        u64::from(scheme_count)
            .checked_mul(SEED_SCHEME_RECORD_SIZE_U64)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let page_end = checked_end(
        get_u64(bytes, 24)?,
        u64::from(page_count)
            .checked_mul(SEED_PAGE_RECORD_SIZE_U64)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let data_offset = get_u64(bytes, 32)?;
    let prefix = usize::try_from(data_offset).map_err(|_| JmaError::OffsetOverflow)?;
    if scheme_end > data_offset || page_end > data_offset || prefix > bytes.len() {
        return Err(JmaError::CorruptSection(
            "seed directory prefix exceeds section".to_string(),
        ));
    }
    Ok(prefix)
}

/// Builds the final sequence directory from the shared codec records. Each
/// block's encoded two-bit payload and ambiguity payload remain independent;
/// they are laid out adjacent to permit one coalesced range request.
pub fn encode_shared_sequence_directory(
    blocks: &[EncodedSequenceBlock],
    contigs: &[ContigMetadata],
) -> JmaResult<(Vec<u8>, Vec<u8>)> {
    let mut ordered = blocks.to_vec();
    ordered.sort_by_key(|block| {
        (
            block.record.contig_id,
            block.record.base_start,
            block.record.payload_offset,
        )
    });
    let mut records = Vec::with_capacity(ordered.len());
    let mut payload = Vec::new();
    for block in &ordered {
        let record = block.record;
        if record.flags != 0 || record.base_count == 0 {
            return Err(JmaError::CorruptSection(
                "sequence block has unsupported flags or an empty range".to_string(),
            ));
        }
        let contig = contigs
            .iter()
            .find(|contig| contig.id == record.contig_id)
            .ok_or(JmaError::UnknownContig(record.contig_id))?;
        let base_end = record.base_end().ok_or(JmaError::OffsetOverflow)?;
        if base_end > contig.length {
            return Err(JmaError::CorruptSection(
                "sequence block exceeds contig length".to_string(),
            ));
        }
        if record.stored_length
            != u64::try_from(block.payload.len()).map_err(|_| JmaError::OffsetOverflow)?
            || record.ambiguity_length
                != u64::try_from(block.ambiguity_payload.len())
                    .map_err(|_| JmaError::OffsetOverflow)?
        {
            return Err(JmaError::CorruptSection(
                "sequence block payload lengths do not match its record".to_string(),
            ));
        }
        decode_stored_block_bounded(
            &record,
            &block.payload,
            &block.ambiguity_payload,
            DEFAULT_MAX_DECODED_BLOCK_BYTES,
        )
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
        let payload_offset = u64::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
        payload.extend_from_slice(&block.payload);
        let ambiguity_offset =
            u64::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
        payload.extend_from_slice(&block.ambiguity_payload);
        records.push(record.with_offsets(payload_offset, ambiguity_offset));
    }
    let directory = SequenceBlockDirectory::new(records)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?
        .encode()
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
    Ok((directory, payload))
}

/// Rebinds shared sequence directory offsets from section-relative to
/// absolute object coordinates after the writer has laid out all sections.
pub fn rebase_shared_sequence_directory(
    bytes: &[u8],
    payload_section_offset: u64,
) -> JmaResult<Vec<u8>> {
    let directory = SequenceBlockDirectory::decode(bytes)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
    let records = directory
        .records
        .into_iter()
        .map(|record| {
            let payload_offset = record
                .payload_offset
                .checked_add(payload_section_offset)
                .ok_or(JmaError::OffsetOverflow)?;
            let ambiguity_offset = record
                .ambiguity_offset
                .checked_add(payload_section_offset)
                .ok_or(JmaError::OffsetOverflow)?;
            Ok(record.with_offsets(payload_offset, ambiguity_offset))
        })
        .collect::<JmaResult<Vec<_>>>()?;
    SequenceBlockDirectory::new(records)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?
        .encode()
        .map_err(|error| JmaError::CorruptSection(error.to_string()))
}

/// Decodes and bounds-checks shared sequence records while retaining absolute
/// object offsets for selective range reads.
pub fn decode_shared_sequence_directory(
    bytes: &[u8],
    payload_size: u64,
    payload_section_offset: u64,
    contigs: &[ContigMetadata],
    object_size: u64,
) -> JmaResult<Vec<SequenceBlockRecord>> {
    let directory = SequenceBlockDirectory::decode(bytes)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
    directory
        .validate_object_size(object_size)
        .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
    let section_end = checked_end(payload_section_offset, payload_size)?;
    let mut ranges = directory
        .records
        .iter()
        .flat_map(|record| {
            [
                (record.payload_offset, record.stored_length),
                (record.ambiguity_offset, record.ambiguity_length),
            ]
        })
        .collect::<Vec<_>>();
    ranges.sort_unstable_by_key(|range| range.0);
    let mut previous_end = payload_section_offset;
    for record in &directory.records {
        if record.payload_offset < payload_section_offset
            || record.payload_end().ok_or(JmaError::OffsetOverflow)? > section_end
            || record.ambiguity_offset < payload_section_offset
            || record.ambiguity_end().ok_or(JmaError::OffsetOverflow)? > section_end
        {
            return Err(JmaError::CorruptSection(
                "sequence block range is outside the payload section".to_string(),
            ));
        }
        let contig = contigs
            .iter()
            .find(|contig| contig.id == record.contig_id)
            .ok_or(JmaError::UnknownContig(record.contig_id))?;
        if record.base_end().ok_or(JmaError::OffsetOverflow)? > contig.length {
            return Err(JmaError::CorruptSection(
                "sequence block exceeds contig length".to_string(),
            ));
        }
    }
    for (offset, length) in ranges {
        if offset < previous_end {
            return Err(JmaError::CorruptSection(
                "sequence payload and ambiguity ranges overlap".to_string(),
            ));
        }
        previous_end = checked_end(offset, length)?;
    }
    Ok(directory.records)
}

/// Builds all seed descriptors and pages. Each scheme is independent and can
/// be opened without decoding another scheme's records.
pub fn encode_seed_index(sections: &[SeedSection]) -> JmaResult<Vec<u8>> {
    let mut schemes = Vec::new();
    let mut identities = std::collections::BTreeSet::new();
    let mut pages = Vec::new();
    let mut page_payload = Vec::new();
    let mut next_scheme_id = 0x4a4d_0000u32;
    let mut ordered_sections = sections.to_vec();
    ordered_sections.sort_by_key(|section| section.k);
    for section in &ordered_sections {
        let mut levels = section.levels.clone();
        levels.sort_by_key(|level| level.level.scale);
        for level in levels {
            if !identities.insert((section.k, level.level.scale)) {
                return Err(JmaError::CorruptSection(
                    "duplicate seed scheme identity".to_string(),
                ));
            }
            let scheme_id = next_scheme_id;
            next_scheme_id = next_scheme_id
                .checked_add(1)
                .ok_or(JmaError::OffsetOverflow)?;
            let density_parameter = u32::try_from(level.level.scale).map_err(|_| {
                JmaError::CorruptSection("seed scale does not fit scheme descriptor".to_string())
            })?;
            let descriptor = SeedSchemeDescriptor {
                scheme_id,
                algorithm_id: 1,
                span: u16::from(section.k),
                informative_bases: u16::from(section.k),
                density_parameter,
                bucket_bits: SEED_BUCKET_BITS,
                key_encoding: 1,
                occurrence_encoding: 1,
                flags: 0,
            };
            let mut records = level.records;
            records.sort_by_key(|record| {
                (
                    record.query.hash,
                    record.query.canonical_kmer,
                    record.occurrence.contig_id,
                    record.occurrence.position,
                    record.occurrence.reverse,
                )
            });
            let first_page = pages.len();
            let mut record_offset = 0usize;
            while record_offset < records.len() {
                let prefix = hash_prefix(records[record_offset].query.hash, SEED_BUCKET_BITS)?;
                let max_page_end = (record_offset + MAX_SEED_PAGE_RECORDS).min(records.len());
                let mut page_records_end = record_offset + 1;
                while page_records_end < max_page_end
                    && hash_prefix(records[page_records_end].query.hash, SEED_BUCKET_BITS)?
                        == prefix
                {
                    page_records_end += 1;
                }
                let page_records = &records[record_offset..page_records_end];
                let key_offset =
                    u64::try_from(page_payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
                let mut keys = Vec::with_capacity(page_records.len() * SEED_KEY_RECORD_SIZE);
                let mut occurrences =
                    Vec::with_capacity(page_records.len() * SEED_OCCURRENCE_RECORD_SIZE);
                for (local, record) in page_records.iter().enumerate() {
                    keys.extend_from_slice(&record.query.hash.to_le_bytes());
                    keys.extend_from_slice(&record.query.canonical_kmer.to_le_bytes());
                    keys.extend_from_slice(
                        &u64::try_from(local)
                            .map_err(|_| JmaError::OffsetOverflow)?
                            .to_le_bytes(),
                    );
                    let occurrence = record.occurrence;
                    occurrences.extend_from_slice(&occurrence.contig_id.to_le_bytes());
                    occurrences.push(u8::from(occurrence.reverse));
                    occurrences.extend_from_slice(&[0u8; 3]);
                    occurrences.extend_from_slice(&occurrence.position.to_le_bytes());
                    occurrences.extend_from_slice(&u16::from(section.k).to_le_bytes());
                    occurrences.extend_from_slice(&0u16.to_le_bytes());
                    occurrences.extend_from_slice(&0u32.to_le_bytes());
                }
                page_payload.extend_from_slice(&keys);
                let occurrence_offset =
                    u64::try_from(page_payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
                page_payload.extend_from_slice(&occurrences);
                let mut digest_input = Vec::with_capacity(keys.len() + occurrences.len());
                digest_input.extend_from_slice(&keys);
                digest_input.extend_from_slice(&occurrences);
                pages.push(SeedPageIndex {
                    scheme_id,
                    page_id: u32::try_from(pages.len() - first_page)
                        .map_err(|_| JmaError::OffsetOverflow)?,
                    hash_prefix: prefix,
                    key_count: u32::try_from(page_records.len())
                        .map_err(|_| JmaError::OffsetOverflow)?,
                    occurrence_count: u32::try_from(page_records.len())
                        .map_err(|_| JmaError::OffsetOverflow)?,
                    key_offset,
                    key_length: u64::try_from(keys.len()).map_err(|_| JmaError::OffsetOverflow)?,
                    occurrence_offset,
                    occurrence_length: u64::try_from(occurrences.len())
                        .map_err(|_| JmaError::OffsetOverflow)?,
                    first_hash: page_records[0].query.hash,
                    last_hash: page_records[page_records.len() - 1].query.hash,
                    checksum: checksum(&digest_input),
                });
                record_offset = page_records_end;
            }
            let page_count = pages.len() - first_page;
            schemes.push(SeedSchemeIndex {
                descriptor,
                pages: pages[first_page..].to_vec(),
                record_count: u64::try_from(records.len()).map_err(|_| JmaError::OffsetOverflow)?,
                declared_page_count: u32::try_from(page_count)
                    .map_err(|_| JmaError::OffsetOverflow)?,
                declared_page_first: 0,
            });
            // `pages` is retained globally only to avoid a second allocation
            // while calculating deterministic page IDs; scheme copies above
            // own the public directory returned to readers.
            let _ = page_count;
        }
    }
    let scheme_count = schemes.len();
    let all_pages = schemes
        .iter()
        .flat_map(|scheme| scheme.pages.iter())
        .collect::<Vec<_>>();
    let scheme_offset = SEED_DIRECTORY_HEADER_SIZE_U64;
    let page_offset = scheme_offset
        .checked_add(
            u64::try_from(scheme_count)
                .map_err(|_| JmaError::OffsetOverflow)?
                .checked_mul(SEED_SCHEME_RECORD_SIZE_U64)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .ok_or(JmaError::OffsetOverflow)?;
    let page_data_offset = page_offset
        .checked_add(
            u64::try_from(all_pages.len())
                .map_err(|_| JmaError::OffsetOverflow)?
                .checked_mul(SEED_PAGE_RECORD_SIZE_U64)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .ok_or(JmaError::OffsetOverflow)?;
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.extend_from_slice(
        &u32::try_from(scheme_count)
            .map_err(|_| JmaError::OffsetOverflow)?
            .to_le_bytes(),
    );
    bytes.extend_from_slice(
        &u32::try_from(all_pages.len())
            .map_err(|_| JmaError::OffsetOverflow)?
            .to_le_bytes(),
    );
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&scheme_offset.to_le_bytes());
    bytes.extend_from_slice(&page_offset.to_le_bytes());
    bytes.extend_from_slice(&page_data_offset.to_le_bytes());
    let mut page_first = 0u32;
    for scheme in &schemes {
        encode_scheme_record(&mut bytes, scheme, page_first)?;
        page_first = page_first
            .checked_add(u32::try_from(scheme.pages.len()).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
    }
    for page in all_pages {
        let mut page = *page;
        page.key_offset = page
            .key_offset
            .checked_add(page_data_offset)
            .ok_or(JmaError::OffsetOverflow)?;
        page.occurrence_offset = page
            .occurrence_offset
            .checked_add(page_data_offset)
            .ok_or(JmaError::OffsetOverflow)?;
        encode_page_record(&mut bytes, page)?;
    }
    bytes.extend_from_slice(&page_payload);
    Ok(bytes)
}

pub fn decode_seed_index_directory(bytes: &[u8]) -> JmaResult<SeedIndexDirectory> {
    let section_length = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
    decode_seed_index_directory_with_length(bytes, section_length)
}

/// Decodes the seed directory prefix while validating page ranges against the
/// complete seed-section length. The payload pages themselves may be omitted
/// from `bytes`, which is what the remote reader uses at archive open.
pub fn decode_seed_index_directory_with_length(
    bytes: &[u8],
    section_length: u64,
) -> JmaResult<SeedIndexDirectory> {
    if bytes.len() < SEED_DIRECTORY_HEADER_SIZE {
        return Err(JmaError::CorruptSection(
            "seed directory header is truncated".to_string(),
        ));
    }
    if get_u32(bytes, 0)? != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported seed directory version".to_string(),
        ));
    }
    let scheme_count = get_u32(bytes, 4)?;
    let page_count = get_u32(bytes, 8)?;
    let scheme_offset = get_u64(bytes, 16)?;
    let page_offset = get_u64(bytes, 24)?;
    let data_offset = get_u64(bytes, 32)?;
    let scheme_bytes = u64::from(scheme_count)
        .checked_mul(SEED_SCHEME_RECORD_SIZE_U64)
        .ok_or(JmaError::OffsetOverflow)?;
    let page_bytes = u64::from(page_count)
        .checked_mul(SEED_PAGE_RECORD_SIZE_U64)
        .ok_or(JmaError::OffsetOverflow)?;
    let scheme_end = checked_end(scheme_offset, scheme_bytes)?;
    let page_end = checked_end(page_offset, page_bytes)?;
    let directory_bytes_length =
        u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
    if scheme_offset < SEED_DIRECTORY_HEADER_SIZE_U64
        || scheme_end > directory_bytes_length
        || page_offset < scheme_end
        || page_end > directory_bytes_length
        || data_offset < page_end
        || data_offset > section_length
        || data_offset > directory_bytes_length
    {
        return Err(JmaError::CorruptSection(
            "seed directory offsets are outside the section".to_string(),
        ));
    }
    let mut schemes = Vec::with_capacity(usize::try_from(scheme_count).unwrap_or(0));
    let mut scheme_ids = std::collections::BTreeSet::new();
    for index in 0..usize::try_from(scheme_count).map_err(|_| JmaError::OffsetOverflow)? {
        let offset = usize::try_from(scheme_offset)
            .map_err(|_| JmaError::OffsetOverflow)?
            .checked_add(
                index
                    .checked_mul(SEED_SCHEME_RECORD_SIZE)
                    .ok_or(JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        let scheme = decode_scheme_record(bytes, offset)?;
        if !scheme_ids.insert(scheme.descriptor.scheme_id) {
            return Err(JmaError::CorruptSection(
                "duplicate seed scheme identifier".to_string(),
            ));
        }
        schemes.push(scheme);
    }
    let mut pages = Vec::with_capacity(usize::try_from(page_count).unwrap_or(0));
    for index in 0..usize::try_from(page_count).map_err(|_| JmaError::OffsetOverflow)? {
        let offset = usize::try_from(page_offset)
            .map_err(|_| JmaError::OffsetOverflow)?
            .checked_add(
                index
                    .checked_mul(SEED_PAGE_RECORD_SIZE)
                    .ok_or(JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        let page = decode_page_record(bytes, offset)?;
        if page.key_offset < data_offset
            || checked_end(page.key_offset, page.key_length)? > section_length
            || page.occurrence_offset < data_offset
            || checked_end(page.occurrence_offset, page.occurrence_length)? > section_length
            || page.key_length
                != u64::from(page.key_count)
                    .checked_mul(SEED_KEY_RECORD_SIZE_U64)
                    .ok_or(JmaError::OffsetOverflow)?
            || page.occurrence_length
                != u64::from(page.occurrence_count)
                    .checked_mul(SEED_OCCURRENCE_RECORD_SIZE_U64)
                    .ok_or(JmaError::OffsetOverflow)?
            || page.key_count != page.occurrence_count
        {
            return Err(JmaError::CorruptSection(
                "seed page range or count is invalid".to_string(),
            ));
        }
        if page.first_hash > page.last_hash
            || hash_prefix(page.first_hash, SEED_BUCKET_BITS)? != page.hash_prefix
            || hash_prefix(page.last_hash, SEED_BUCKET_BITS)? != page.hash_prefix
        {
            return Err(JmaError::CorruptSection(
                "seed page hash range is invalid".to_string(),
            ));
        }
        pages.push(page);
    }
    if pages
        .iter()
        .any(|page| !scheme_ids.contains(&page.scheme_id))
    {
        return Err(JmaError::CorruptSection(
            "seed page references an unknown scheme identifier".to_string(),
        ));
    }
    for occurrence_kind in [true, false] {
        let mut ranges = pages
            .iter()
            .map(|page| {
                if occurrence_kind {
                    (page.key_offset, page.key_length)
                } else {
                    (page.occurrence_offset, page.occurrence_length)
                }
            })
            .collect::<Vec<_>>();
        ranges.sort_unstable_by_key(|range| range.0);
        let mut previous_end = data_offset;
        for (offset, length) in ranges {
            if offset < previous_end {
                return Err(JmaError::CorruptSection(
                    "seed page ranges overlap".to_string(),
                ));
            }
            previous_end = checked_end(offset, length)?;
        }
    }
    let mut all_ranges = pages
        .iter()
        .flat_map(|page| {
            [
                (page.key_offset, page.key_length),
                (page.occurrence_offset, page.occurrence_length),
            ]
        })
        .collect::<Vec<_>>();
    all_ranges.sort_unstable_by_key(|range| range.0);
    let mut previous_end = data_offset;
    for (offset, length) in all_ranges {
        if offset < previous_end {
            return Err(JmaError::CorruptSection(
                "seed key and occurrence ranges overlap".to_string(),
            ));
        }
        previous_end = checked_end(offset, length)?;
    }
    for scheme in &mut schemes {
        let mut selected = pages
            .iter()
            .filter(|page| page.scheme_id == scheme.descriptor.scheme_id)
            .copied()
            .collect::<Vec<_>>();
        selected.sort_by_key(|page| (page.hash_prefix, page.page_id));
        if u32::try_from(selected.len()).map_err(|_| JmaError::OffsetOverflow)?
            != scheme.declared_page_count
        {
            return Err(JmaError::CorruptSection(
                "seed scheme page count does not match its page directory".to_string(),
            ));
        }
        let first =
            usize::try_from(scheme.declared_page_first).map_err(|_| JmaError::OffsetOverflow)?;
        let end = first
            .checked_add(
                usize::try_from(scheme.declared_page_count)
                    .map_err(|_| JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        if end > pages.len()
            || pages[first..end]
                .iter()
                .any(|page| page.scheme_id != scheme.descriptor.scheme_id)
        {
            return Err(JmaError::CorruptSection(
                "seed scheme page ownership range is invalid".to_string(),
            ));
        }
        for pair in selected.windows(2) {
            if (pair[1].hash_prefix, pair[1].page_id) <= (pair[0].hash_prefix, pair[0].page_id)
                || pair[1].first_hash < pair[0].first_hash
                || pair[1].page_id != pair[0].page_id.saturating_add(1)
            {
                return Err(JmaError::CorruptSection(
                    "seed pages are not sorted deterministically".to_string(),
                ));
            }
        }
        if selected.first().is_some_and(|page| page.page_id != 0) {
            return Err(JmaError::CorruptSection(
                "seed page IDs do not start at zero".to_string(),
            ));
        }
        let selected_records = selected.iter().try_fold(0u64, |total, page| {
            total
                .checked_add(u64::from(page.key_count))
                .ok_or(JmaError::OffsetOverflow)
        })?;
        if selected_records != scheme.record_count {
            return Err(JmaError::CorruptSection(
                "seed scheme record count does not match its pages".to_string(),
            ));
        }
        let mut key_ranges = selected
            .iter()
            .map(|page| (page.key_offset, page.key_length))
            .collect::<Vec<_>>();
        key_ranges.sort_unstable_by_key(|range| range.0);
        let mut occurrence_ranges = selected
            .iter()
            .map(|page| (page.occurrence_offset, page.occurrence_length))
            .collect::<Vec<_>>();
        occurrence_ranges.sort_unstable_by_key(|range| range.0);
        for ranges in [key_ranges, occurrence_ranges] {
            let mut previous_end = data_offset;
            for (offset, length) in ranges {
                if offset < previous_end {
                    return Err(JmaError::CorruptSection(
                        "seed page ranges overlap".to_string(),
                    ));
                }
                previous_end = checked_end(offset, length)?;
            }
        }
        scheme.pages = selected;
        if scheme
            .pages
            .iter()
            .any(|page| page.scheme_id != scheme.descriptor.scheme_id)
        {
            return Err(JmaError::CorruptSection(
                "seed page references an unknown scheme".to_string(),
            ));
        }
    }
    Ok(SeedIndexDirectory { schemes })
}

/// Decodes one pair of fixed key/occurrence ranges after checksum validation.
pub fn decode_seed_page(
    key_bytes: &[u8],
    occurrence_bytes: &[u8],
    page: &SeedPageIndex,
    k: u8,
) -> JmaResult<Vec<SeedRecord>> {
    if key_bytes.len() != usize::try_from(page.key_length).map_err(|_| JmaError::OffsetOverflow)?
        || occurrence_bytes.len()
            != usize::try_from(page.occurrence_length).map_err(|_| JmaError::OffsetOverflow)?
    {
        return Err(JmaError::CorruptSection(
            "seed page range length does not match directory".to_string(),
        ));
    }
    let mut digest_input = Vec::with_capacity(key_bytes.len() + occurrence_bytes.len());
    digest_input.extend_from_slice(key_bytes);
    digest_input.extend_from_slice(occurrence_bytes);
    if checksum(&digest_input) != page.checksum {
        return Err(JmaError::ChecksumMismatch(format!(
            "seed page scheme={} page={}",
            page.scheme_id, page.page_id
        )));
    }
    let mut records = Vec::with_capacity(usize::try_from(page.key_count).unwrap_or(0));
    for index in 0..usize::try_from(page.key_count).map_err(|_| JmaError::OffsetOverflow)? {
        let key_offset = index
            .checked_mul(SEED_KEY_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        let occurrence_index = get_u64(key_bytes, key_offset + 16)?;
        let occurrence_offset = usize::try_from(occurrence_index)
            .ok()
            .and_then(|index| index.checked_mul(SEED_OCCURRENCE_RECORD_SIZE))
            .ok_or(JmaError::OffsetOverflow)?;
        let hash = get_u64(key_bytes, key_offset)?;
        let canonical_kmer = get_u64(key_bytes, key_offset + 8)?;
        if hash_prefix(hash, SEED_BUCKET_BITS)? != page.hash_prefix
            || hash < page.first_hash
            || hash > page.last_hash
        {
            return Err(JmaError::CorruptSection(
                "seed key is outside its page hash range".to_string(),
            ));
        }
        if !valid_canonical_kmer(k, canonical_kmer)
            || occurrence_offset + SEED_OCCURRENCE_RECORD_SIZE > occurrence_bytes.len()
        {
            return Err(JmaError::CorruptSection(
                "seed page key verification or occurrence index is invalid".to_string(),
            ));
        }
        let contig_id = get_u32(occurrence_bytes, occurrence_offset)?;
        let reverse = match occurrence_bytes[occurrence_offset + 4] {
            0 => false,
            1 => true,
            value => {
                return Err(JmaError::CorruptSection(format!(
                    "invalid seed orientation {value}"
                )));
            }
        };
        if occurrence_bytes[occurrence_offset + 5..occurrence_offset + 8]
            .iter()
            .any(|byte| *byte != 0)
            || occurrence_bytes[occurrence_offset + 18..occurrence_offset + 24]
                .iter()
                .any(|byte| *byte != 0)
        {
            return Err(JmaError::CorruptSection(
                "seed occurrence reserved bytes are non-zero".to_string(),
            ));
        }
        let position = get_u64(occurrence_bytes, occurrence_offset + 8)?;
        let span = u16::from_le_bytes([
            occurrence_bytes[occurrence_offset + 16],
            occurrence_bytes[occurrence_offset + 17],
        ]);
        if span != u16::from(k) {
            return Err(JmaError::CorruptSection(
                "seed occurrence span does not match scheme".to_string(),
            ));
        }
        records.push(SeedRecord {
            query: crate::jma::SeedQuery {
                k,
                hash,
                canonical_kmer,
            },
            occurrence: crate::jma::SeedOccurrence {
                contig_id,
                position,
                reverse,
            },
        });
    }
    Ok(records)
}

fn encode_scheme_record(
    bytes: &mut Vec<u8>,
    scheme: &SeedSchemeIndex,
    first_page: u32,
) -> JmaResult<()> {
    let descriptor = scheme.descriptor;
    let page_count = u32::try_from(scheme.pages.len()).map_err(|_| JmaError::OffsetOverflow)?;
    bytes.extend_from_slice(&descriptor.scheme_id.to_le_bytes());
    bytes.extend_from_slice(&descriptor.algorithm_id.to_le_bytes());
    bytes.extend_from_slice(&descriptor.span.to_le_bytes());
    bytes.extend_from_slice(&descriptor.informative_bases.to_le_bytes());
    bytes.extend_from_slice(&descriptor.density_parameter.to_le_bytes());
    bytes.push(descriptor.bucket_bits);
    bytes.push(descriptor.key_encoding);
    bytes.push(descriptor.occurrence_encoding);
    bytes.push(0);
    bytes.extend_from_slice(&descriptor.flags.to_le_bytes());
    bytes.extend_from_slice(&first_page.to_le_bytes());
    bytes.extend_from_slice(&page_count.to_le_bytes());
    bytes.extend_from_slice(&scheme.record_count.to_le_bytes());
    bytes.resize(bytes.len() + (SEED_SCHEME_RECORD_SIZE - 40), 0);
    Ok(())
}

fn decode_scheme_record(bytes: &[u8], offset: usize) -> JmaResult<SeedSchemeIndex> {
    if offset + SEED_SCHEME_RECORD_SIZE > bytes.len() {
        return Err(JmaError::CorruptSection(
            "seed scheme record is truncated".to_string(),
        ));
    }
    let descriptor = SeedSchemeDescriptor {
        scheme_id: get_u32(bytes, offset)?,
        algorithm_id: get_u32(bytes, offset + 4)?,
        span: u16::from_le_bytes([bytes[offset + 8], bytes[offset + 9]]),
        informative_bases: u16::from_le_bytes([bytes[offset + 10], bytes[offset + 11]]),
        density_parameter: get_u32(bytes, offset + 12)?,
        bucket_bits: bytes[offset + 16],
        key_encoding: bytes[offset + 17],
        occurrence_encoding: bytes[offset + 18],
        flags: get_u32(bytes, offset + 20)?,
    };
    if descriptor.bucket_bits != SEED_BUCKET_BITS
        || descriptor.key_encoding != 1
        || descriptor.occurrence_encoding != 1
        || !(21..=32).contains(&descriptor.span)
        || descriptor.density_parameter == 0
    {
        return Err(JmaError::CorruptSection(
            "unsupported or invalid seed scheme descriptor".to_string(),
        ));
    }
    let page_count = get_u32(bytes, offset + 28)?;
    let record_count = get_u64(bytes, offset + 32)?;
    if bytes[offset + 19] != 0
        || bytes[offset + 40..offset + SEED_SCHEME_RECORD_SIZE]
            .iter()
            .any(|byte| *byte != 0)
    {
        return Err(JmaError::CorruptSection(
            "seed scheme record reserved bytes are non-zero".to_string(),
        ));
    }
    Ok(SeedSchemeIndex {
        descriptor,
        pages: Vec::with_capacity(usize::try_from(page_count).unwrap_or(0)),
        record_count,
        declared_page_count: page_count,
        declared_page_first: get_u32(bytes, offset + 24)?,
    })
}

fn encode_page_record(bytes: &mut Vec<u8>, page: SeedPageIndex) -> JmaResult<()> {
    bytes.extend_from_slice(&page.scheme_id.to_le_bytes());
    bytes.extend_from_slice(&page.page_id.to_le_bytes());
    bytes.extend_from_slice(&page.hash_prefix.to_le_bytes());
    bytes.extend_from_slice(&page.key_count.to_le_bytes());
    bytes.extend_from_slice(&page.occurrence_count.to_le_bytes());
    bytes.extend_from_slice(&0u32.to_le_bytes());
    bytes.extend_from_slice(&page.key_offset.to_le_bytes());
    bytes.extend_from_slice(&page.key_length.to_le_bytes());
    bytes.extend_from_slice(&page.occurrence_offset.to_le_bytes());
    bytes.extend_from_slice(&page.occurrence_length.to_le_bytes());
    bytes.extend_from_slice(&page.first_hash.to_le_bytes());
    bytes.extend_from_slice(&page.last_hash.to_le_bytes());
    bytes.extend_from_slice(&page.checksum);
    Ok(())
}

fn decode_page_record(bytes: &[u8], offset: usize) -> JmaResult<SeedPageIndex> {
    if offset + SEED_PAGE_RECORD_SIZE > bytes.len() {
        return Err(JmaError::CorruptSection(
            "seed page record is truncated".to_string(),
        ));
    }
    if bytes[offset + 20..offset + 24]
        .iter()
        .any(|byte| *byte != 0)
    {
        return Err(JmaError::CorruptSection(
            "seed page reserved bytes are non-zero".to_string(),
        ));
    }
    Ok(SeedPageIndex {
        scheme_id: get_u32(bytes, offset)?,
        page_id: get_u32(bytes, offset + 4)?,
        hash_prefix: get_u32(bytes, offset + 8)?,
        key_count: get_u32(bytes, offset + 12)?,
        occurrence_count: get_u32(bytes, offset + 16)?,
        key_offset: get_u64(bytes, offset + 24)?,
        key_length: get_u64(bytes, offset + 32)?,
        occurrence_offset: get_u64(bytes, offset + 40)?,
        occurrence_length: get_u64(bytes, offset + 48)?,
        first_hash: get_u64(bytes, offset + 56)?,
        last_hash: get_u64(bytes, offset + 64)?,
        checksum: crate::jma::format::array_32_at(bytes, offset + 72)?,
    })
}
