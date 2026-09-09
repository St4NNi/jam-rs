use crate::owner_postings::MAX_KEYS_PER_BLOCK;
use std::io;
use thiserror::Error;

pub(crate) const OWNER_MAGIC: [u8; 8] = *b"JOWNER\0\0";
pub(crate) const OWNER_VERSION: u16 = 2;
pub(crate) const OWNER_HEADER_SIZE: usize = 512;
pub(crate) const OWNER_PAGE_SIZE: u64 = 4096;
pub(crate) const OWNER_SECTION_COUNT: usize = 8;
pub(crate) const OWNER_DOCUMENT_SIZE: u32 = 96;
pub(crate) const OWNER_CONTIG_SIZE: u32 = 40;
pub(crate) const OWNER_BLOCK_SIZE: u32 = 64;
const MAX_CHECKSUM_LEVELS: usize = 16;
pub(crate) const OWNER_SECTION_TABLE: usize = 136;
pub(crate) const OWNER_SECTION_DESCRIPTOR_SIZE: usize = 24;
pub(crate) const COMPLETE_RANGE: u32 = 1;
pub(crate) const HAS_METADATA: u32 = 2;

pub(crate) fn metadata_digest(parts: [&[u8]; 4]) -> [u8; 32] {
    use sha2::{Digest, Sha256};
    let mut digest = Sha256::new();
    digest.update(b"JOWNER-METADATA-V2\0");
    for part in parts {
        digest.update((part.len() as u64).to_le_bytes());
        digest.update(part);
    }
    digest.finalize().into()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u16)]
pub(crate) enum OwnerSection {
    Strings = 1,
    Documents = 2,
    Contigs = 3,
    Gzi = 4,
    BlockDirectory = 5,
    HotPostings = 6,
    ColdPostings = 7,
    PageChecksums = 8,
}

impl OwnerSection {
    pub(crate) const ALL: [Self; OWNER_SECTION_COUNT] = [
        Self::Strings,
        Self::Documents,
        Self::Contigs,
        Self::Gzi,
        Self::BlockDirectory,
        Self::HotPostings,
        Self::ColdPostings,
        Self::PageChecksums,
    ];

    pub(crate) const fn record_size(self) -> u32 {
        match self {
            Self::Documents => OWNER_DOCUMENT_SIZE,
            Self::Contigs => OWNER_CONTIG_SIZE,
            Self::BlockDirectory => OWNER_BLOCK_SIZE,
            Self::PageChecksums => 32,
            _ => 0,
        }
    }

    fn from_code(value: u16) -> Result<Self, OwnerReaderError> {
        Self::ALL
            .into_iter()
            .find(|section| *section as u16 == value)
            .ok_or(OwnerReaderError::Invalid("section kind"))
    }

    pub(crate) const fn index(self) -> usize {
        self as usize - 1
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OwnerSectionDescriptor {
    pub(crate) kind: OwnerSection,
    pub(crate) record_size: u32,
    pub(crate) offset: u64,
    pub(crate) length: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct ChecksumLevel {
    pub(crate) offset: u64,
    pub(crate) hash_count: u64,
    pub(crate) page_count: u64,
}

pub(crate) fn checksum_layout(data_pages: u64) -> Result<Vec<ChecksumLevel>, OwnerReaderError> {
    let mut levels = Vec::new();
    let mut hash_count = data_pages;
    let mut offset = 0u64;
    loop {
        let page_count = hash_count.max(1).div_ceil(OWNER_PAGE_SIZE / 32);
        levels.push(ChecksumLevel {
            offset,
            hash_count,
            page_count,
        });
        if levels.len() > MAX_CHECKSUM_LEVELS {
            return Err(OwnerReaderError::Invalid("checksum depth"));
        }
        offset = offset
            .checked_add(
                page_count
                    .checked_mul(OWNER_PAGE_SIZE)
                    .ok_or(OwnerReaderError::Invalid("checksum length"))?,
            )
            .ok_or(OwnerReaderError::Invalid("checksum length"))?;
        if page_count == 1 {
            break;
        }
        hash_count = page_count;
    }
    Ok(levels)
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OwnerHeader {
    pub(crate) flags: u32,
    pub(crate) k: u8,
    pub(crate) rescue_k15: bool,
    pub(crate) minimizer_window: u16,
    pub(crate) owner_ordinal: u32,
    pub(crate) owner_count: u32,
    pub(crate) first_key: u64,
    pub(crate) last_key: u64,
    pub(crate) key_count: u64,
    pub(crate) occurrence_count: u64,
    pub(crate) document_count: u32,
    pub(crate) contig_count: u32,
    pub(crate) generation_id: [u8; 32],
    pub(crate) body_sha256: [u8; 32],
    pub(crate) checksum_root_sha256: [u8; 32],
    pub(crate) metadata_sha256: [u8; 32],
    pub(crate) sections: [OwnerSectionDescriptor; OWNER_SECTION_COUNT],
}

impl OwnerHeader {
    pub(crate) fn complete_range(&self) -> bool {
        self.flags & COMPLETE_RANGE != 0
    }

    pub(crate) fn has_metadata(&self) -> bool {
        self.flags & HAS_METADATA != 0
    }

    pub(crate) fn encode(&self) -> Result<[u8; OWNER_HEADER_SIZE], OwnerReaderError> {
        self.validate(self.file_len()?)?;
        let mut out = [0; OWNER_HEADER_SIZE];
        out[..8].copy_from_slice(&OWNER_MAGIC);
        put_u16(&mut out, 8, OWNER_VERSION);
        put_u16(&mut out, 10, OWNER_HEADER_SIZE as u16);
        put_u32(&mut out, 12, self.flags);
        out[16] = self.k;
        out[17] = u8::from(self.rescue_k15);
        put_u16(&mut out, 20, self.minimizer_window);
        put_u16(&mut out, 22, OWNER_SECTION_COUNT as u16);
        put_u32(&mut out, 24, self.owner_ordinal);
        put_u32(&mut out, 28, self.owner_count);
        put_u64(&mut out, 32, self.first_key);
        put_u64(&mut out, 40, self.last_key);
        put_u64(&mut out, 48, self.key_count);
        put_u64(&mut out, 56, self.occurrence_count);
        put_u32(&mut out, 64, self.document_count);
        put_u32(&mut out, 68, self.contig_count);
        out[72..104].copy_from_slice(&self.generation_id);
        out[104..136].copy_from_slice(&self.body_sha256);
        out[328..360].copy_from_slice(&self.checksum_root_sha256);
        out[360..392].copy_from_slice(&self.metadata_sha256);
        for (index, section) in self.sections.iter().enumerate() {
            let start = OWNER_SECTION_TABLE + index * OWNER_SECTION_DESCRIPTOR_SIZE;
            put_u16(&mut out, start, section.kind as u16);
            put_u32(&mut out, start + 4, section.record_size);
            put_u64(&mut out, start + 8, section.offset);
            put_u64(&mut out, start + 16, section.length);
        }
        Ok(out)
    }

    pub(crate) fn decode(bytes: &[u8], file_len: u64) -> Result<Self, OwnerReaderError> {
        if bytes.len() != OWNER_HEADER_SIZE || bytes[..8] != OWNER_MAGIC {
            return Err(OwnerReaderError::Invalid("header"));
        }
        if read_u16(bytes, 8) != OWNER_VERSION {
            return Err(OwnerReaderError::Invalid("version"));
        }
        if read_u16(bytes, 10) != OWNER_HEADER_SIZE as u16
            || read_u16(bytes, 22) != OWNER_SECTION_COUNT as u16
            || bytes[18..20].iter().any(|byte| *byte != 0)
            || bytes[392..].iter().any(|byte| *byte != 0)
        {
            return Err(OwnerReaderError::Invalid("header reservation"));
        }
        let mut sections = Vec::with_capacity(OWNER_SECTION_COUNT);
        for index in 0..OWNER_SECTION_COUNT {
            let start = OWNER_SECTION_TABLE + index * OWNER_SECTION_DESCRIPTOR_SIZE;
            if read_u16(bytes, start + 2) != 0 {
                return Err(OwnerReaderError::Invalid("section flags"));
            }
            sections.push(OwnerSectionDescriptor {
                kind: OwnerSection::from_code(read_u16(bytes, start))?,
                record_size: read_u32(bytes, start + 4),
                offset: read_u64(bytes, start + 8),
                length: read_u64(bytes, start + 16),
            });
        }
        let header = Self {
            flags: read_u32(bytes, 12),
            k: bytes[16],
            rescue_k15: match bytes[17] {
                0 => false,
                1 => true,
                _ => return Err(OwnerReaderError::Invalid("rescue seed flag")),
            },
            minimizer_window: read_u16(bytes, 20),
            owner_ordinal: read_u32(bytes, 24),
            owner_count: read_u32(bytes, 28),
            first_key: read_u64(bytes, 32),
            last_key: read_u64(bytes, 40),
            key_count: read_u64(bytes, 48),
            occurrence_count: read_u64(bytes, 56),
            document_count: read_u32(bytes, 64),
            contig_count: read_u32(bytes, 68),
            generation_id: bytes[72..104].try_into().expect("owner generation ID"),
            body_sha256: bytes[104..136].try_into().expect("owner body digest"),
            checksum_root_sha256: bytes[328..360].try_into().expect("owner checksum root"),
            metadata_sha256: bytes[360..392].try_into().expect("owner metadata digest"),
            sections: sections.try_into().expect("fixed owner section count"),
        };
        header.validate(file_len)?;
        Ok(header)
    }

    pub(crate) fn section(&self, section: OwnerSection) -> OwnerSectionDescriptor {
        self.sections[section.index()]
    }

    fn file_len(&self) -> Result<u64, OwnerReaderError> {
        let last = self.sections[OWNER_SECTION_COUNT - 1];
        last.offset
            .checked_add(last.length)
            .ok_or(OwnerReaderError::Invalid("file length"))
    }

    fn validate(&self, file_len: u64) -> Result<(), OwnerReaderError> {
        if self.flags & !(COMPLETE_RANGE | HAS_METADATA) != 0
            || !(1..=32).contains(&self.k)
            || (self.rescue_k15 && self.k != 21)
            || self.minimizer_window == 0
            || self.owner_count == 0
            || self.owner_ordinal >= self.owner_count
            || self.first_key > self.last_key
            || self.generation_id == [0; 32]
            || self.body_sha256 == [0; 32]
            || self.checksum_root_sha256 == [0; 32]
            || self.metadata_sha256 == [0; 32]
        {
            return Err(OwnerReaderError::Invalid("header values"));
        }
        if self.document_count == 0 || self.contig_count == 0 {
            return Err(OwnerReaderError::Invalid("metadata owner"));
        }
        let mut previous_end = OWNER_HEADER_SIZE as u64;
        for (expected, section) in OwnerSection::ALL.into_iter().zip(self.sections) {
            if section.kind != expected
                || section.record_size != expected.record_size()
                || section.offset < previous_end
                || !section.offset.is_multiple_of(OWNER_PAGE_SIZE)
            {
                return Err(OwnerReaderError::Invalid("section layout"));
            }
            previous_end = section
                .offset
                .checked_add(section.length)
                .ok_or(OwnerReaderError::Invalid("section range"))?;
            if previous_end > file_len {
                return Err(OwnerReaderError::Invalid("section range"));
            }
        }
        let document_bytes = if self.has_metadata() {
            u64::from(self.document_count) * u64::from(OWNER_DOCUMENT_SIZE)
        } else {
            0
        };
        let contig_bytes = if self.has_metadata() {
            u64::from(self.contig_count) * u64::from(OWNER_CONTIG_SIZE)
        } else {
            0
        };
        let block_count =
            self.section(OwnerSection::BlockDirectory).length / u64::from(OWNER_BLOCK_SIZE);
        if previous_end != file_len
            || self.section(OwnerSection::Documents).length != document_bytes
            || self.section(OwnerSection::Contigs).length > contig_bytes
            || !self
                .section(OwnerSection::Contigs)
                .length
                .is_multiple_of(u64::from(OWNER_CONTIG_SIZE))
            || !self
                .section(OwnerSection::BlockDirectory)
                .length
                .is_multiple_of(u64::from(OWNER_BLOCK_SIZE))
            || (self.key_count == 0) != (block_count == 0)
            || self.key_count < block_count
            || self.key_count > block_count.saturating_mul(MAX_KEYS_PER_BLOCK as u64)
        {
            return Err(OwnerReaderError::Invalid("section length"));
        }
        let checksums = self.section(OwnerSection::PageChecksums);
        let checksum_levels = checksum_layout(checksums.offset / OWNER_PAGE_SIZE - 1)?;
        let checksum_length = checksum_levels.iter().try_fold(0u64, |total, level| {
            total
                .checked_add(
                    level
                        .page_count
                        .checked_mul(OWNER_PAGE_SIZE)
                        .ok_or(OwnerReaderError::Invalid("checksum length"))?,
                )
                .ok_or(OwnerReaderError::Invalid("checksum length"))
        })?;
        if checksums.length != checksum_length {
            return Err(OwnerReaderError::Invalid("checksum length"));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OwnerSeed {
    pub(crate) packed_key: u64,
    pub(crate) document_frequency: u64,
    pub(crate) block_ordinal: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct OwnerDocument {
    pub(crate) metagenome_id: u32,
    pub(crate) occurrence_count: u64,
    pub(crate) seed_key: u64,
    pub(crate) owner_ordinal: u32,
    pub(crate) block_ordinal: u64,
    pub(crate) member_ordinal: u64,
}

#[derive(Clone, Copy)]
pub(crate) struct BlockRecord {
    pub(crate) first_key: u64,
    pub(crate) last_key: u64,
    pub(crate) hot_offset: u64,
    pub(crate) hot_length: u64,
    pub(crate) cold_offset: u64,
    pub(crate) cold_length: u64,
    pub(crate) key_count: u32,
}

impl BlockRecord {
    pub(crate) fn decode(bytes: &[u8]) -> Result<Self, OwnerReaderError> {
        if bytes.len() != OWNER_BLOCK_SIZE as usize || bytes[52..].iter().any(|byte| *byte != 0) {
            return Err(OwnerReaderError::Invalid("block record"));
        }
        Ok(Self {
            first_key: read_u64(bytes, 0),
            last_key: read_u64(bytes, 8),
            hot_offset: read_u64(bytes, 16),
            hot_length: read_u64(bytes, 24),
            cold_offset: read_u64(bytes, 32),
            cold_length: read_u64(bytes, 40),
            key_count: read_u32(bytes, 48),
        })
    }
}

#[derive(Clone, Copy)]
pub(crate) struct DocumentRecord {
    pub(crate) name_offset: u32,
    pub(crate) name_length: u32,
    pub(crate) uri_offset: u32,
    pub(crate) uri_length: u32,
    pub(crate) bgzf_bytes: u64,
    pub(crate) bgzf_sha256: [u8; 32],
    pub(crate) contig_start: u32,
    pub(crate) contig_count: u32,
    pub(crate) gzi_offset: u64,
    pub(crate) gzi_length: u64,
    pub(crate) original_contig_start: u32,
    pub(crate) original_contig_count: u32,
    pub(crate) locus_bits: u8,
}

impl DocumentRecord {
    pub(crate) fn decode(bytes: &[u8]) -> Result<Self, OwnerReaderError> {
        if bytes.len() != OWNER_DOCUMENT_SIZE as usize
            || bytes[89..].iter().any(|byte| *byte != 0)
            || !(1..=64).contains(&bytes[88])
        {
            return Err(OwnerReaderError::Invalid("document record"));
        }
        Ok(Self {
            name_offset: read_u32(bytes, 0),
            name_length: read_u32(bytes, 4),
            uri_offset: read_u32(bytes, 8),
            uri_length: read_u32(bytes, 12),
            bgzf_bytes: read_u64(bytes, 16),
            bgzf_sha256: bytes[24..56].try_into().expect("BGZF digest"),
            contig_start: read_u32(bytes, 56),
            contig_count: read_u32(bytes, 60),
            gzi_offset: read_u64(bytes, 64),
            gzi_length: read_u64(bytes, 72),
            original_contig_start: read_u32(bytes, 80),
            original_contig_count: read_u32(bytes, 84),
            locus_bits: bytes[88],
        })
    }
}

pub(crate) fn read_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().expect("owner u16"))
}
pub(crate) fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("owner u32"))
}
pub(crate) fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("owner u64"))
}
pub(crate) fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}
pub(crate) fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}
pub(crate) fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

#[derive(Debug, Error)]
pub enum OwnerReaderError {
    #[error("owner index I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("invalid owner index: {0}")]
    Invalid(&'static str),
    #[error("owner index page checksum mismatch")]
    ChecksumMismatch,
    #[error("owner index key {0} is outside the converted ranges")]
    KeyNotCovered(u64),
    #[error(transparent)]
    Postings(#[from] crate::owner_postings::OwnerPostingsError),
}
