//! Checked local-mmap and range readers for JAM Witness Router format 1.

use crate::resource::{ByteRange, RangeReader};
use crate::router::WitnessKey;
use crate::router::format::{
    KeyBlockDescriptor, MembershipFilter, MetagenomeEntry, ObjectRange, PostingHeader, PrefixEntry,
    RouterFormatError, RouterHeader, SECTION_FILTER, SECTION_KEYS, SECTION_METAGENOMES,
    SECTION_POSITION_PAYLOAD, SECTION_POSTING_HEADERS, SECTION_POSTING_PAYLOAD, SECTION_PREFIX,
    SECTION_SCHEMES, SUPERBLOCK_SIZE, SectionDescriptor, SectionDirectory, array_32, checked_end,
    checked_usize, checksum, decode_filter, decode_header, decode_metagenomes, decode_schemes,
    decode_section_directory, get_u32, get_u64, hash_prefix, range_in_object,
    verify_object_checksum, verify_section_checksum,
};
use memmap2::Mmap;
use std::collections::BTreeSet;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use thiserror::Error;

const KEY_HEADER_SIZE: usize = 40;
const KEY_BLOCK_RECORD_SIZE: usize = 88;
const PREFIX_HEADER_SIZE: usize = 24;
const PREFIX_RECORD_SIZE: usize = 56;
const POSTING_HEADER_SIZE: usize = 24;
const POSTING_RECORD_SIZE: usize = 80;

/// A key match returned only after the packed canonical k-mer was checked.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KeyMatch {
    pub key: WitnessKey,
    pub key_index: u64,
    pub block_id: u32,
}

/// Raw bounded data for a selected witness. The reader does not interpret
/// metagenome or position codecs; callers can connect the codec appropriate
/// to their collection profile.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RawPosting {
    pub key: KeyMatch,
    pub header: PostingHeader,
    pub posting: Vec<u8>,
    pub positions: Option<Vec<u8>>,
}

#[derive(Debug, Error)]
pub enum RouterReaderError {
    #[error(transparent)]
    Format(#[from] RouterFormatError),
    #[error("router resource read failed: {0}")]
    Resource(String),
    #[error("router local mmap failed: {0}")]
    Io(String),
    #[error("router object metadata size {actual} does not match declared size {declared}")]
    ObjectSize { actual: u64, declared: u64 },
    #[error("router section {0} is unavailable")]
    SectionUnavailable(u32),
    #[error("router key index {0} is out of bounds")]
    KeyIndexOutOfBounds(u64),
    #[error("router key block {0} failed its checksum")]
    KeyBlockChecksum(u32),
    #[error("router posting {0} failed its checksum")]
    PostingChecksum(u64),
}

trait RouterSource: Send + Sync {
    fn size(&self) -> Result<u64, RouterReaderError>;
    fn read(&self, range: ObjectRange) -> Result<Vec<u8>, RouterReaderError>;
    fn full_bytes(&self) -> Option<&[u8]> {
        None
    }
}

struct RangeSource<R> {
    resource: R,
    size: u64,
}

impl<R: RangeReader> RouterSource for RangeSource<R> {
    fn size(&self) -> Result<u64, RouterReaderError> {
        Ok(self.size)
    }

    fn read(&self, range: ObjectRange) -> Result<Vec<u8>, RouterReaderError> {
        let range = ByteRange::new(range.offset, range.length)
            .map_err(|error| RouterReaderError::Resource(error.to_string()))?;
        self.resource
            .read_range(range)
            .map_err(|error| RouterReaderError::Resource(error.to_string()))
    }
}

struct MmapSource {
    mmap: Mmap,
}

impl RouterSource for MmapSource {
    fn size(&self) -> Result<u64, RouterReaderError> {
        u64::try_from(self.mmap.len())
            .map_err(|_| RouterReaderError::Format(RouterFormatError::OffsetOverflow))
    }

    fn read(&self, range: ObjectRange) -> Result<Vec<u8>, RouterReaderError> {
        let start = checked_usize(range.offset, "mmap range offset")?;
        let length = checked_usize(range.length, "mmap range length")?;
        let end = start
            .checked_add(length)
            .ok_or(RouterFormatError::OffsetOverflow)?;
        self.mmap
            .get(start..end)
            .map(ToOwned::to_owned)
            .ok_or_else(|| {
                RouterReaderError::Format(RouterFormatError::RangeOutOfBounds {
                    what: "mmap range",
                    offset: range.offset,
                    length: range.length,
                    object_size: self.mmap.len() as u64,
                })
            })
    }

    fn full_bytes(&self) -> Option<&[u8]> {
        Some(&self.mmap)
    }
}

/// A parsed immutable router. Directories and membership data are read once;
/// exact key blocks and posting payloads remain lazy range reads.
pub struct RouterReader {
    source: Arc<dyn RouterSource>,
    header: RouterHeader,
    directory: SectionDirectory,
    metagenomes: Vec<MetagenomeEntry>,
    schemes: Vec<crate::router::WitnessScheme>,
    filter: MembershipFilter,
    key_config: KeyConfig,
    key_blocks: Vec<KeyBlockDescriptor>,
    prefix_bits: u8,
    prefixes: Vec<PrefixEntry>,
    posting_headers: Vec<PostingHeader>,
    key_section: SectionDescriptor,
    posting_section: SectionDescriptor,
    position_section: Option<SectionDescriptor>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct KeyConfig {
    width: usize,
    store_jamhash: bool,
    prefix_bits: u8,
    key_count: u64,
}

impl RouterReader {
    /// Open through any checked local or remote range reader. The object
    /// checksum remains a deferred operation for remote resources so opening
    /// does not trigger a full-object download.
    pub fn open<R: RangeReader + 'static>(resource: R) -> Result<Self, RouterReaderError> {
        let metadata = resource
            .metadata()
            .map_err(|error| RouterReaderError::Resource(error.to_string()))?;
        let size = metadata.size;
        let source: Arc<dyn RouterSource> = Arc::new(RangeSource { resource, size });
        Self::open_source(source)
    }

    /// Open an immutable local file with a read-only memory map. The full
    /// identity checksum is checked before returning, because the mmap makes
    /// that validation local and does not add a remote request.
    pub fn open_mmap(path: impl AsRef<Path>) -> Result<Self, RouterReaderError> {
        let file =
            File::open(path.as_ref()).map_err(|error| RouterReaderError::Io(error.to_string()))?;
        let mmap = unsafe { Mmap::map(&file) }
            .map_err(|error| RouterReaderError::Io(error.to_string()))?;
        let source: Arc<dyn RouterSource> = Arc::new(MmapSource { mmap });
        let reader = Self::open_source(source)?;
        reader.verify_identity()?;
        Ok(reader)
    }

    fn open_source(source: Arc<dyn RouterSource>) -> Result<Self, RouterReaderError> {
        let actual_size = source.size()?;
        let header_bytes = source.read(ObjectRange::new(0, SUPERBLOCK_SIZE as u64)?)?;
        let header = decode_header(&header_bytes)?;
        if actual_size != header.object_size {
            return Err(RouterReaderError::ObjectSize {
                actual: actual_size,
                declared: header.object_size,
            });
        }
        let directory_range = range_in_object(
            header.section_directory_offset,
            header.section_directory_length,
            header.object_size,
            "section directory",
        )?;
        let directory_bytes = source.read(directory_range)?;
        if checksum(&directory_bytes) != header.section_directory_checksum {
            return Err(RouterReaderError::Format(
                RouterFormatError::ChecksumMismatch("section directory".to_string()),
            ));
        }
        let directory = decode_section_directory(
            &directory_bytes,
            header.section_count,
            header.object_size,
            header.section_directory_offset,
            header.section_directory_length,
        )?;
        let metagenome_section = directory.required(SECTION_METAGENOMES)?;
        let scheme_section = directory.required(SECTION_SCHEMES)?;
        let filter_section = directory.required(SECTION_FILTER)?;
        let prefix_section = directory.required(SECTION_PREFIX)?;
        let key_section = directory.required(SECTION_KEYS)?;
        let posting_headers_section = directory.required(SECTION_POSTING_HEADERS)?;
        let posting_section = directory.required(SECTION_POSTING_PAYLOAD)?;
        let metagenome_bytes =
            read_checked_section(&*source, metagenome_section, header.object_size)?;
        let scheme_bytes = read_checked_section(&*source, scheme_section, header.object_size)?;
        let filter_bytes = read_checked_section(&*source, filter_section, header.object_size)?;
        let prefix_bytes = read_checked_section(&*source, prefix_section, header.object_size)?;
        let posting_header_bytes =
            read_checked_section(&*source, posting_headers_section, header.object_size)?;
        if checksum(&metagenome_bytes) != header.catalog_checksum {
            return Err(RouterReaderError::Format(
                RouterFormatError::ChecksumMismatch("metagenome identity".to_string()),
            ));
        }
        if checksum(&scheme_bytes) != header.witness_scheme_checksum {
            return Err(RouterReaderError::Format(
                RouterFormatError::ChecksumMismatch("witness scheme identity".to_string()),
            ));
        }
        let metagenomes = decode_metagenomes(&metagenome_bytes)?;
        let schemes = decode_schemes(&scheme_bytes)?;
        if u32::try_from(metagenomes.len()).ok() != Some(header.metagenome_count)
            || u32::try_from(schemes.len()).ok() != Some(header.scheme_count)
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidDirectory(
                    "superblock table counts do not match sections".to_string(),
                ),
            ));
        }
        let filter = decode_filter(&filter_bytes)?;
        let key_bytes = read_key_directory(&*source, key_section, header.object_size)?;
        let (key_config, key_blocks) = decode_key_directory(&key_bytes, key_section.length)?;
        let (prefix_bits, prefixes) = decode_prefix_directory(&prefix_bytes)?;
        if prefix_bits != key_config.prefix_bits {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidPrefixDirectory,
            ));
        }
        validate_prefixes(
            &prefixes,
            &key_blocks,
            key_section.length,
            key_config.width,
            key_config.store_jamhash,
            prefix_bits,
        )?;
        let posting_headers = decode_posting_headers(
            &posting_header_bytes,
            posting_section.length,
            directory
                .get(SECTION_POSITION_PAYLOAD)
                .map_or(0, |entry| entry.length),
        )?;
        if u64::try_from(posting_headers.len()).ok() != Some(key_config.key_count) {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidPostingHeaders,
            ));
        }
        let position_section = directory.get(SECTION_POSITION_PAYLOAD);
        if posting_headers
            .iter()
            .any(|header| header.position_length != 0)
            && position_section.is_none()
        {
            return Err(RouterReaderError::SectionUnavailable(
                SECTION_POSITION_PAYLOAD,
            ));
        }
        Ok(Self {
            source,
            header,
            directory,
            metagenomes,
            schemes,
            filter,
            key_config,
            key_blocks,
            prefix_bits,
            prefixes,
            posting_headers,
            key_section,
            posting_section,
            position_section,
        })
    }

    #[must_use]
    pub fn header(&self) -> &RouterHeader {
        &self.header
    }

    #[must_use]
    pub fn sections(&self) -> &SectionDirectory {
        &self.directory
    }

    #[must_use]
    pub fn metagenomes(&self) -> &[MetagenomeEntry] {
        &self.metagenomes
    }

    #[must_use]
    pub fn schemes(&self) -> &[crate::router::WitnessScheme] {
        &self.schemes
    }

    #[must_use]
    pub fn key_count(&self) -> u64 {
        self.key_config.key_count
    }

    #[must_use]
    pub fn filter_maybe_contains(&self, key: WitnessKey) -> bool {
        self.filter.contains_hash(key.jamhash)
    }

    /// Verify the whole identity binding. This operation is deliberately
    /// explicit for remote readers because it necessarily reads the object.
    pub fn verify_identity(&self) -> Result<(), RouterReaderError> {
        if let Some(bytes) = self.source.full_bytes() {
            verify_object_checksum(bytes, &self.header)?;
            return Ok(());
        }
        let bytes = self
            .source
            .read(ObjectRange::new(0, self.header.object_size)?)?;
        verify_object_checksum(&bytes, &self.header)?;
        Ok(())
    }

    /// Look up one exact key. A filter miss returns before any exact-key block
    /// range is requested; a filter hit still verifies packed equality.
    pub fn lookup(&self, key: WitnessKey) -> Result<Option<KeyMatch>, RouterReaderError> {
        if !self.filter.contains_hash(key.jamhash) {
            return Ok(None);
        }
        let prefix = hash_prefix(key.jamhash, self.prefix_bits)?;
        let mut seen_blocks = BTreeSet::new();
        for entry in self.prefixes.iter().filter(|entry| entry.prefix == prefix) {
            if key.jamhash < entry.first_hash || key.jamhash > entry.last_hash {
                continue;
            }
            if !seen_blocks.insert(entry.block_id) {
                continue;
            }
            let descriptor =
                self.key_blocks
                    .get(entry.block_id as usize)
                    .ok_or(RouterReaderError::Format(
                        RouterFormatError::InvalidPrefixDirectory,
                    ))?;
            let block_range = range_in_object(
                self.key_section
                    .offset
                    .checked_add(descriptor.offset)
                    .ok_or(RouterFormatError::OffsetOverflow)?,
                descriptor.length,
                self.header.object_size,
                "key block",
            )?;
            let bytes = self.source.read(block_range)?;
            if checksum(&bytes) != descriptor.checksum {
                return Err(RouterReaderError::KeyBlockChecksum(descriptor.block_id));
            }
            let record_bytes =
                self.key_config.width + if self.key_config.store_jamhash { 8 } else { 0 };
            if bytes.len() != descriptor.key_count as usize * record_bytes {
                return Err(RouterReaderError::Format(
                    RouterFormatError::InvalidKeyBlock,
                ));
            }
            for local in 0..descriptor.key_count as usize {
                let offset = local * record_bytes;
                let packed = decode_packed(&bytes[offset..offset + self.key_config.width]);
                let jamhash_offset = offset + self.key_config.width;
                let stored_hash = if self.key_config.store_jamhash {
                    get_u64(&bytes, jamhash_offset)?
                } else {
                    crate::jamhash_u64_v1(packed)
                };
                if stored_hash != key.jamhash {
                    continue;
                }
                let checked = match WitnessKey::checked(packed, stored_hash) {
                    Ok(checked) => checked,
                    // Equal hashes are only a routing hint. A different
                    // packed canonical key is a non-match, never evidence.
                    Err(crate::router::RouterContractError::HashMismatch) => continue,
                    Err(_) => {
                        return Err(RouterReaderError::Format(
                            RouterFormatError::InvalidKeyBlock,
                        ));
                    }
                };
                if checked.packed == key.packed {
                    return Ok(Some(KeyMatch {
                        key: checked,
                        key_index: descriptor.first_key_index + local as u64,
                        block_id: descriptor.block_id,
                    }));
                }
            }
        }
        Ok(None)
    }

    pub fn posting_header(&self, key_index: u64) -> Result<PostingHeader, RouterReaderError> {
        self.posting_headers
            .get(checked_usize(key_index, "key index")?)
            .copied()
            .ok_or(RouterReaderError::KeyIndexOutOfBounds(key_index))
    }

    /// Return bounded raw posting and optional position slices, validating the
    /// per-key checksum after both slices are fetched.
    pub fn read_posting(&self, key: KeyMatch) -> Result<RawPosting, RouterReaderError> {
        let header = self.posting_header(key.key_index)?;
        let posting_range = range_in_object(
            self.posting_section
                .offset
                .checked_add(header.posting_offset)
                .ok_or(RouterFormatError::OffsetOverflow)?,
            header.posting_length,
            self.header.object_size,
            "posting payload",
        )?;
        let posting = self.source.read(posting_range)?;
        let positions = if header.position_length != 0 {
            let section = self
                .position_section
                .ok_or(RouterReaderError::SectionUnavailable(
                    SECTION_POSITION_PAYLOAD,
                ))?;
            let range = range_in_object(
                section
                    .offset
                    .checked_add(header.position_offset)
                    .ok_or(RouterFormatError::OffsetOverflow)?,
                header.position_length,
                self.header.object_size,
                "position payload",
            )?;
            Some(self.source.read(range)?)
        } else {
            None
        };
        let mut digest_input =
            Vec::with_capacity(posting.len() + positions.as_ref().map_or(0, Vec::len));
        digest_input.extend_from_slice(&posting);
        if let Some(positions) = &positions {
            digest_input.extend_from_slice(positions);
        }
        if checksum(&digest_input) != header.checksum {
            return Err(RouterReaderError::PostingChecksum(key.key_index));
        }
        Ok(RawPosting {
            key,
            header,
            posting,
            positions,
        })
    }
}

fn read_checked_section(
    source: &dyn RouterSource,
    descriptor: SectionDescriptor,
    object_size: u64,
) -> Result<Vec<u8>, RouterReaderError> {
    let range = range_in_object(descriptor.offset, descriptor.length, object_size, "section")?;
    let bytes = source.read(range)?;
    if bytes.len() != checked_usize(descriptor.length, "section length")? {
        return Err(RouterReaderError::Format(RouterFormatError::Truncated {
            expected: checked_usize(descriptor.length, "section length")?,
            actual: bytes.len(),
        }));
    }
    verify_section_checksum(&bytes, descriptor)?;
    Ok(bytes)
}

fn read_key_directory(
    source: &dyn RouterSource,
    descriptor: SectionDescriptor,
    object_size: u64,
) -> Result<Vec<u8>, RouterReaderError> {
    let header_range = range_in_object(
        descriptor.offset,
        KEY_HEADER_SIZE as u64,
        object_size,
        "key header",
    )?;
    let header = source.read(header_range)?;
    if header.len() < KEY_HEADER_SIZE || get_u32(&header, 0)? != 1 {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidKeyBlock,
        ));
    }
    let block_count = get_u32(&header, 8)?;
    let block_bytes = u64::from(block_count)
        .checked_mul(KEY_BLOCK_RECORD_SIZE as u64)
        .ok_or(RouterFormatError::OffsetOverflow)?;
    let directory_offset = get_u64(&header, 20)?;
    let directory_end = checked_end(directory_offset, block_bytes)?;
    if directory_offset < KEY_HEADER_SIZE as u64 || directory_end > descriptor.length {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidKeyBlock,
        ));
    }
    let range = range_in_object(
        descriptor.offset,
        directory_end,
        object_size,
        "key directory",
    )?;
    source.read(range)
}

fn decode_key_directory(
    bytes: &[u8],
    section_length: u64,
) -> Result<(KeyConfig, Vec<KeyBlockDescriptor>), RouterReaderError> {
    if bytes.len() < KEY_HEADER_SIZE || get_u32(bytes, 0)? != 1 {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidKeyBlock,
        ));
    }
    let width = match bytes[4] {
        6 | 8 => bytes[4] as usize,
        other => {
            return Err(RouterReaderError::Format(
                RouterFormatError::UnsupportedKeyWidth(other),
            ));
        }
    };
    let store_jamhash = match bytes[5] {
        0 => false,
        1 => true,
        _ => {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidKeyBlock,
            ));
        }
    };
    let prefix_bits = bytes[6];
    if bytes[7] != 0 {
        return Err(RouterReaderError::Format(
            RouterFormatError::NonZeroReserved("key header"),
        ));
    }
    let block_count = usize::try_from(get_u32(bytes, 8)?)
        .map_err(|_| RouterReaderError::Format(RouterFormatError::CountTooLarge("key blocks")))?;
    let key_count = get_u64(bytes, 12)?;
    let directory_offset = checked_usize(get_u64(bytes, 20)?, "key directory")?;
    let data_offset = checked_usize(get_u64(bytes, 28)?, "key data")?;
    let expected_data_offset = directory_offset
        .checked_add(
            block_count
                .checked_mul(KEY_BLOCK_RECORD_SIZE)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    if directory_offset != KEY_HEADER_SIZE
        || data_offset != expected_data_offset
        || data_offset > bytes.len()
    {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidKeyBlock,
        ));
    }
    let mut blocks = Vec::with_capacity(block_count);
    let mut previous_hash = 0;
    let mut previous_index = 0;
    let record_bytes = u64::try_from(width + if store_jamhash { 8 } else { 0 })
        .map_err(|_| RouterFormatError::OffsetOverflow)?;
    for index in 0..block_count {
        let offset = directory_offset + index * KEY_BLOCK_RECORD_SIZE;
        let descriptor = KeyBlockDescriptor {
            block_id: get_u32(bytes, offset)?,
            first_key_index: get_u64(bytes, offset + 8)?,
            key_count: get_u32(bytes, offset + 16)?,
            offset: get_u64(bytes, offset + 24)?,
            length: get_u64(bytes, offset + 32)?,
            first_hash: get_u64(bytes, offset + 40)?,
            last_hash: get_u64(bytes, offset + 48)?,
            checksum: array_32(bytes, offset + 56)?,
        };
        if get_u32(bytes, offset + 4)? != 0
            || get_u32(bytes, offset + 20)? != 0
            || descriptor.block_id != index as u32
            || descriptor.key_count == 0
            || descriptor.first_hash > descriptor.last_hash
            || (index != 0
                && (descriptor.first_key_index < previous_index
                    || descriptor.first_hash < previous_hash))
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidKeyBlock,
            ));
        }
        if (index == 0 && descriptor.first_key_index != 0)
            || (index != 0 && descriptor.first_key_index != previous_index)
            || descriptor.offset < data_offset as u64
            || checked_end(descriptor.offset, descriptor.length)? > section_length
            || descriptor.length
                != checked_end(
                    0,
                    u64::from(descriptor.key_count)
                        .checked_mul(record_bytes)
                        .ok_or(RouterFormatError::OffsetOverflow)?,
                )?
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidKeyBlock,
            ));
        }
        previous_hash = descriptor.last_hash;
        previous_index = checked_end(descriptor.first_key_index, u64::from(descriptor.key_count))?;
        blocks.push(descriptor);
    }
    if blocks.last().map_or(key_count != 0, |block| {
        checked_end(block.first_key_index, u64::from(block.key_count)) != Ok(key_count)
    }) {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidKeyBlock,
        ));
    }
    Ok((
        KeyConfig {
            width,
            store_jamhash,
            prefix_bits,
            key_count,
        },
        blocks,
    ))
}

fn decode_prefix_directory(bytes: &[u8]) -> Result<(u8, Vec<PrefixEntry>), RouterReaderError> {
    if bytes.len() < PREFIX_HEADER_SIZE || get_u32(bytes, 0)? != 1 {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidPrefixDirectory,
        ));
    }
    let prefix_bits = bytes[4];
    let entry_count = usize::try_from(get_u32(bytes, 8)?).map_err(|_| {
        RouterReaderError::Format(RouterFormatError::CountTooLarge("prefix entries"))
    })?;
    let records_offset = checked_usize(get_u64(bytes, 12)?, "prefix records")?;
    let expected_length = records_offset
        .checked_add(
            entry_count
                .checked_mul(PREFIX_RECORD_SIZE)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    if records_offset != PREFIX_HEADER_SIZE
        || expected_length != bytes.len()
        || prefix_bits == 0
        || prefix_bits > 32
        || bytes[5..8].iter().any(|byte| *byte != 0)
    {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidPrefixDirectory,
        ));
    }
    let mut entries = Vec::with_capacity(entry_count);
    let mut previous = None;
    for index in 0..entry_count {
        let offset = records_offset + index * PREFIX_RECORD_SIZE;
        let entry = PrefixEntry {
            prefix: get_u32(bytes, offset)?,
            block_id: get_u32(bytes, offset + 4)?,
            first_key_index: get_u64(bytes, offset + 8)?,
            key_count: get_u32(bytes, offset + 16)?,
            offset: get_u64(bytes, offset + 24)?,
            length: get_u64(bytes, offset + 32)?,
            first_hash: get_u64(bytes, offset + 40)?,
            last_hash: get_u64(bytes, offset + 48)?,
        };
        if (prefix_bits < 32 && entry.prefix >= (1u64 << prefix_bits) as u32)
            || entry.key_count == 0
            || entry.first_hash > entry.last_hash
            || previous.is_some_and(|previous| (entry.prefix, entry.block_id) < previous)
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidPrefixDirectory,
            ));
        }
        previous = Some((entry.prefix, entry.block_id));
        entries.push(entry);
    }
    Ok((prefix_bits, entries))
}

fn validate_prefixes(
    entries: &[PrefixEntry],
    blocks: &[KeyBlockDescriptor],
    key_section_length: u64,
    width: usize,
    store_jamhash: bool,
    prefix_bits: u8,
) -> Result<(), RouterReaderError> {
    let record_bytes = u64::try_from(width + if store_jamhash { 8 } else { 0 })
        .map_err(|_| RouterFormatError::OffsetOverflow)?;
    for entry in entries {
        let block = blocks
            .get(entry.block_id as usize)
            .ok_or(RouterReaderError::Format(
                RouterFormatError::InvalidPrefixDirectory,
            ))?;
        let entry_end = checked_end(entry.first_key_index, u64::from(entry.key_count))?;
        let block_end = checked_end(block.first_key_index, u64::from(block.key_count))?;
        let entry_byte_end = checked_end(entry.offset, entry.length)?;
        let block_byte_end = checked_end(block.offset, block.length)?;
        let expected_entry_length = u64::from(entry.key_count)
            .checked_mul(record_bytes)
            .ok_or(RouterFormatError::OffsetOverflow)?;
        if entry.first_key_index < block.first_key_index
            || entry_end > block_end
            || entry.offset < block.offset
            || entry_byte_end > block_byte_end
            || entry.length != expected_entry_length
            || entry_byte_end > key_section_length
            || hash_prefix(entry.first_hash, prefix_bits)? != entry.prefix
            || hash_prefix(entry.last_hash, prefix_bits)? != entry.prefix
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidPrefixDirectory,
            ));
        }
    }
    Ok(())
}

fn decode_posting_headers(
    bytes: &[u8],
    posting_section_length: u64,
    position_section_length: u64,
) -> Result<Vec<PostingHeader>, RouterReaderError> {
    if bytes.len() < POSTING_HEADER_SIZE
        || get_u32(bytes, 0)? != 1
        || get_u32(bytes, 4)? as usize != POSTING_RECORD_SIZE
    {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidPostingHeaders,
        ));
    }
    let count = usize::try_from(get_u64(bytes, 8)?).map_err(|_| {
        RouterReaderError::Format(RouterFormatError::CountTooLarge("posting headers"))
    })?;
    let records_offset = checked_usize(get_u64(bytes, 16)?, "posting records")?;
    let expected_length = records_offset
        .checked_add(
            count
                .checked_mul(POSTING_RECORD_SIZE)
                .ok_or(RouterFormatError::OffsetOverflow)?,
        )
        .ok_or(RouterFormatError::OffsetOverflow)?;
    if records_offset != POSTING_HEADER_SIZE || expected_length != bytes.len() {
        return Err(RouterReaderError::Format(
            RouterFormatError::InvalidPostingHeaders,
        ));
    }
    let mut headers = Vec::with_capacity(count);
    for index in 0..count {
        let offset = records_offset + index * POSTING_RECORD_SIZE;
        let header = PostingHeader {
            document_frequency: get_u32(bytes, offset)?,
            posting_count: get_u32(bytes, offset + 4)?,
            flags: get_u32(bytes, offset + 8)?,
            posting_offset: get_u64(bytes, offset + 16)?,
            posting_length: get_u64(bytes, offset + 24)?,
            position_offset: get_u64(bytes, offset + 32)?,
            position_length: get_u64(bytes, offset + 40)?,
            checksum: array_32(bytes, offset + 48)?,
        };
        if get_u32(bytes, offset + 12)? != 0
            || checked_end(header.posting_offset, header.posting_length)? > posting_section_length
            || checked_end(header.position_offset, header.position_length)?
                > position_section_length
            || (header.position_length != 0
                && header.flags & PostingHeader::FLAG_POSITION_BEARING == 0)
        {
            return Err(RouterReaderError::Format(
                RouterFormatError::InvalidPostingHeaders,
            ));
        }
        headers.push(header);
    }
    Ok(headers)
}

fn decode_packed(bytes: &[u8]) -> u64 {
    let mut raw = [0u8; 8];
    raw[..bytes.len()].copy_from_slice(bytes);
    u64::from_le_bytes(raw)
}
