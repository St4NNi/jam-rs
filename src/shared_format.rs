use crate::jidx::{put_u16, put_u32, put_u64, sha256};
use crate::owner_format::checksum_layout;
use std::io;
use thiserror::Error;

pub(crate) use crate::jidx::{read_u32, read_u64};
pub(crate) const PAGE_BYTES: u64 = 4096;
pub(crate) const HEADER_BYTES: usize = 4096;
pub(crate) const CORE_ROW_BYTES: u64 = 24;
pub(crate) const GROUP_ROW_BYTES: u64 = 32;
pub(crate) const MEMBER_ROW_BYTES: u64 = 24;
pub(crate) const OCCURRENCE_ROW_BYTES: u64 = 24;
pub(crate) const MULTIPLE_CORE: u32 = 1 << 31;
pub(crate) const CORE_MASK: u32 = (1 << 30) - 1;
pub(crate) const CORE_PREFIX_BOUNDARIES: usize = 65_537;
const HEADER_HASH_OFFSET: usize = 304;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(usize)]
pub(crate) enum Section {
    Strings,
    Documents,
    Contigs,
    Gzi,
    Cores,
    Groups,
    Members,
    References,
    Occurrences,
    Checksums,
    CorePayloads,
    CorePrefixes,
    CoreFilter,
}

impl Section {
    pub(crate) const ALL: [Self; 10] = [
        Self::Strings,
        Self::Documents,
        Self::Contigs,
        Self::Gzi,
        Self::Cores,
        Self::Groups,
        Self::Members,
        Self::References,
        Self::Occurrences,
        Self::Checksums,
    ];

    const SPLIT: [Self; 12] = [
        Self::Strings,
        Self::Documents,
        Self::Contigs,
        Self::Gzi,
        Self::CorePrefixes,
        Self::Cores,
        Self::CorePayloads,
        Self::Groups,
        Self::Members,
        Self::References,
        Self::Occurrences,
        Self::Checksums,
    ];

    const FILTERED: [Self; 13] = [
        Self::Strings,
        Self::Documents,
        Self::Contigs,
        Self::Gzi,
        Self::CorePrefixes,
        Self::Cores,
        Self::CorePayloads,
        Self::Groups,
        Self::Members,
        Self::References,
        Self::Occurrences,
        Self::CoreFilter,
        Self::Checksums,
    ];

    pub(crate) fn row_bytes(self) -> u64 {
        match self {
            Self::Documents => 80,
            Self::Contigs => 40,
            Self::Cores => CORE_ROW_BYTES,
            Self::Groups => GROUP_ROW_BYTES,
            Self::Members => MEMBER_ROW_BYTES,
            Self::References => 8,
            Self::Occurrences => OCCURRENCE_ROW_BYTES,
            Self::Checksums => 32,
            Self::CorePrefixes => 4,
            _ => 1,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) struct SectionRange {
    pub(crate) offset: u64,
    pub(crate) length: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SharedHeader {
    pub(crate) version: u16,
    pub(crate) core_payload_bytes: u8,
    pub(crate) window: u16,
    pub(crate) core_count: u64,
    pub(crate) occurrence_count: u64,
    pub(crate) document_count: u32,
    pub(crate) contig_count: u32,
    pub(crate) source_bases: u64,
    pub(crate) manifest_sha256: [u8; 32],
    pub(crate) body_sha256: [u8; 32],
    pub(crate) checksum_root_sha256: [u8; 32],
    pub(crate) sections: [SectionRange; 13],
    pub(crate) filter_source_sha256: [u8; 32],
}

impl SharedHeader {
    pub(crate) fn section_order(&self) -> &'static [Section] {
        if self.version == 4 {
            &Section::FILTERED
        } else if self.version >= 3 {
            &Section::SPLIT
        } else {
            &Section::ALL
        }
    }

    fn hash_offset(&self) -> usize {
        HEADER_HASH_OFFSET
            + match self.version {
                4 => 48,
                3 => 32,
                _ => 0,
            }
    }

    pub(crate) fn id_bytes(&self) -> usize {
        if self.document_count <= u32::from(u8::MAX) {
            1
        } else if self.document_count <= u32::from(u16::MAX) {
            2
        } else {
            4
        }
    }

    pub(crate) fn row_bytes(&self, section: Section) -> u64 {
        if self.version >= 3 {
            match section {
                Section::Cores => return 4,
                Section::CorePayloads => return u64::from(self.core_payload_bytes),
                _ => {}
            }
        }
        if matches!(self.version, 2..=4) {
            match section {
                Section::Groups => return 13 + self.id_bytes() as u64,
                Section::Members => return 8 + self.id_bytes() as u64,
                Section::References => return 4,
                _ => {}
            }
        }
        section.row_bytes()
    }

    pub(crate) fn section(&self, section: Section) -> SectionRange {
        self.sections[section as usize]
    }

    pub(crate) fn encode(&self) -> Result<[u8; HEADER_BYTES], SharedError> {
        let end = self.section(Section::Checksums);
        self.validate(
            end.offset
                .checked_add(end.length)
                .ok_or(SharedError::Invalid("file length"))?,
        )?;
        let mut out = [0; HEADER_BYTES];
        out[..8].copy_from_slice(b"JSHARED\0");
        put_u16(&mut out, 8, self.version);
        put_u16(&mut out, 10, HEADER_BYTES as u16);
        put_u16(&mut out, 12, self.window);
        out[14] = self.core_payload_bytes;
        put_u64(&mut out, 16, self.core_count);
        put_u64(&mut out, 24, self.occurrence_count);
        put_u32(&mut out, 32, self.document_count);
        put_u32(&mut out, 36, self.contig_count);
        put_u64(&mut out, 40, self.source_bases);
        out[48..80].copy_from_slice(&self.manifest_sha256);
        out[80..112].copy_from_slice(&self.body_sha256);
        out[112..144].copy_from_slice(&self.checksum_root_sha256);
        for (index, section) in self
            .sections
            .iter()
            .take(self.section_order().len())
            .enumerate()
        {
            put_u64(&mut out, 144 + 16 * index, section.offset);
            put_u64(&mut out, 152 + 16 * index, section.length);
        }
        if self.version == 4 {
            out[384..416].copy_from_slice(&self.filter_source_sha256);
        }
        let digest = sha256(&out);
        out[self.hash_offset()..self.hash_offset() + 32].copy_from_slice(&digest);
        Ok(out)
    }

    pub(crate) fn decode(bytes: &[u8], file_bytes: u64) -> Result<Self, SharedError> {
        let filtered = bytes.get(8) == Some(&4);
        let split = filtered || bytes.get(8) == Some(&3);
        let hash_offset = HEADER_HASH_OFFSET
            + if filtered {
                48
            } else if split {
                32
            } else {
                0
            };
        if bytes.len() != HEADER_BYTES
            || &bytes[..8] != b"JSHARED\0"
            || !matches!(bytes[8], 1..=4)
            || bytes[9..12] != [0, 0, 16]
            || bytes[15] != 0
            || (!split && bytes[14] != 0)
            || bytes[if filtered { 416 } else { hash_offset + 32 }..]
                .iter()
                .any(|&byte| byte != 0)
        {
            return Err(SharedError::Invalid("header"));
        }
        let mut checked = [0; HEADER_BYTES];
        checked.copy_from_slice(bytes);
        checked[hash_offset..hash_offset + 32].fill(0);
        if sha256(&checked) != bytes[hash_offset..hash_offset + 32] {
            return Err(SharedError::ChecksumMismatch);
        }
        let header = Self {
            version: u16::from_le_bytes(bytes[8..10].try_into().unwrap()),
            core_payload_bytes: bytes[14],
            window: u16::from_le_bytes(bytes[12..14].try_into().unwrap()),
            core_count: read_u64(bytes, 16),
            occurrence_count: read_u64(bytes, 24),
            document_count: read_u32(bytes, 32),
            contig_count: read_u32(bytes, 36),
            source_bases: read_u64(bytes, 40),
            manifest_sha256: bytes[48..80].try_into().unwrap(),
            body_sha256: bytes[80..112].try_into().unwrap(),
            checksum_root_sha256: bytes[112..144].try_into().unwrap(),
            filter_source_sha256: if filtered {
                bytes[384..416].try_into().unwrap()
            } else {
                [0; 32]
            },
            sections: std::array::from_fn(|index| {
                if index
                    >= if filtered {
                        13
                    } else if split {
                        12
                    } else {
                        10
                    }
                {
                    return SectionRange::default();
                }
                SectionRange {
                    offset: read_u64(bytes, 144 + 16 * index),
                    length: read_u64(bytes, 152 + 16 * index),
                }
            }),
        };
        header.validate(file_bytes)?;
        Ok(header)
    }

    fn validate(&self, file_bytes: u64) -> Result<(), SharedError> {
        if !matches!(
            (self.version, self.core_payload_bytes),
            (1 | 2, 0) | (3 | 4, 13 | 21)
        ) || self.window == 0
            || self.document_count == 0
            || self.contig_count == 0
            || self.core_count > self.occurrence_count
            || self.version >= 3 && self.core_count > u64::from(u32::MAX)
        {
            return Err(SharedError::Invalid("header counts"));
        }
        let mut previous = HEADER_BYTES as u64;
        for &kind in self.section_order() {
            let section = self.section(kind);
            let expected = previous
                .checked_next_multiple_of(PAGE_BYTES)
                .ok_or(SharedError::Invalid("section padding"))?;
            if section.offset != expected || !section.length.is_multiple_of(self.row_bytes(kind)) {
                return Err(SharedError::Invalid("section layout"));
            }
            previous = section
                .offset
                .checked_add(section.length)
                .filter(|&end| end <= file_bytes)
                .ok_or(SharedError::Invalid("section extent"))?;
        }
        if previous != file_bytes
            || self.section(Section::Documents).length != u64::from(self.document_count) * 80
            || self.section(Section::Contigs).length != u64::from(self.contig_count) * 40
            || self.core_count.checked_mul(self.row_bytes(Section::Cores))
                != Some(self.section(Section::Cores).length)
            || self.version >= 3
                && (self
                    .core_count
                    .checked_mul(u64::from(self.core_payload_bytes))
                    != Some(self.section(Section::CorePayloads).length)
                    || self.section(Section::CorePrefixes).length
                        != CORE_PREFIX_BOUNDARIES as u64 * 4)
        {
            return Err(SharedError::Invalid("section counts"));
        }
        let checksums = self.section(Section::Checksums);
        let levels = checksum_layout(checksums.offset / PAGE_BYTES - 1)
            .map_err(|_| SharedError::Invalid("checksum layout"))?;
        let top = levels
            .last()
            .ok_or(SharedError::Invalid("checksum layout"))?;
        if checksums.length != top.offset + top.page_count * PAGE_BYTES {
            return Err(SharedError::Invalid("checksum extent"));
        }
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum SharedError {
    #[error("shared-anchor I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("invalid shared-anchor {0}")]
    Invalid(&'static str),
    #[error("shared-anchor checksum mismatch")]
    ChecksumMismatch,
    #[error("shared-anchor file changed after opening")]
    SourceChanged,
    #[error("shared-anchor decoded work exceeds the declared byte budget")]
    ResourceLimit,
    #[error(transparent)]
    Index(#[from] crate::jidx_reader::JidxReaderError),
    #[error(transparent)]
    Seed(#[from] crate::jidx_builder::JidxBuildError),
    #[error(transparent)]
    Bgzf(#[from] crate::bgzf::BgzfError),
    #[error(transparent)]
    Record(#[from] crate::jidx::JidxError),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> SharedHeader {
        let lengths = [5, 80, 40, 8, 24, 0, 0, 0, 0, 4096];
        let mut offset = HEADER_BYTES as u64;
        let sections = std::array::from_fn(|index| {
            if index >= lengths.len() {
                return SectionRange::default();
            }
            offset = offset.next_multiple_of(PAGE_BYTES);
            let section = SectionRange {
                offset,
                length: lengths[index],
            };
            offset += section.length;
            section
        });
        SharedHeader {
            version: 1,
            core_payload_bytes: 0,
            window: 64,
            core_count: 1,
            occurrence_count: 1,
            document_count: 1,
            contig_count: 1,
            source_bases: 100,
            manifest_sha256: [1; 32],
            body_sha256: [2; 32],
            checksum_root_sha256: [3; 32],
            sections,
            filter_source_sha256: [0; 32],
        }
    }

    #[test]
    fn header_roundtrip_and_counts_are_checked() {
        let header = fixture();
        let end = header.section(Section::Checksums);
        assert_eq!(
            SharedHeader::decode(&header.encode().unwrap(), end.offset + end.length).unwrap(),
            header
        );
        let mut invalid = header.clone();
        invalid.core_count = u64::MAX;
        invalid.occurrence_count = u64::MAX;
        assert!(invalid.encode().is_err());
        let mut invalid = header.clone();
        invalid.sections[Section::Members as usize].length = u64::MAX;
        assert!(invalid.encode().is_err());
        let mut invalid = header;
        invalid.sections[Section::Contigs as usize].offset -= 4096;
        assert!(invalid.encode().is_err());
    }

    #[test]
    fn header_authentication_includes_reserved_bytes() {
        let header = fixture();
        let end = header.section(Section::Checksums);
        let mut encoded = header.encode().unwrap();
        encoded[48] ^= 1;
        assert!(matches!(
            SharedHeader::decode(&encoded, end.offset + end.length),
            Err(SharedError::ChecksumMismatch)
        ));
        let mut encoded = header.encode().unwrap();
        encoded[HEADER_BYTES - 1] = 1;
        assert!(SharedHeader::decode(&encoded, end.offset + end.length).is_err());
    }

    #[test]
    fn packed_rows_use_checked_document_widths() {
        let mut header = fixture();
        header.version = 2;
        let end = header.section(Section::Checksums);
        assert_eq!(
            SharedHeader::decode(&header.encode().unwrap(), end.offset + end.length).unwrap(),
            header
        );
        for (documents, width) in [(255, 1), (256, 2), (65535, 2), (65536, 4), (u32::MAX, 4)] {
            header.document_count = documents;
            assert_eq!(header.id_bytes(), width);
            assert_eq!(header.row_bytes(Section::Groups), 13 + width as u64);
            assert_eq!(header.row_bytes(Section::Members), 8 + width as u64);
            assert_eq!(header.row_bytes(Section::References), 4);
        }
    }

    #[test]
    fn split_core_header_binds_width_sections_and_prefix_count() {
        for width in [13, 21] {
            let mut header = fixture();
            header.version = 3;
            header.core_payload_bytes = width;
            header.sections[Section::Cores as usize].length = 4;
            header.sections[Section::CorePayloads as usize].length = u64::from(width);
            header.sections[Section::CorePrefixes as usize].length =
                CORE_PREFIX_BOUNDARIES as u64 * 4;
            let mut offset = HEADER_BYTES as u64;
            for &kind in header.section_order() {
                offset = offset.next_multiple_of(PAGE_BYTES);
                header.sections[kind as usize].offset = offset;
                offset += header.section(kind).length;
            }
            let checksums = header.section(Section::Checksums);
            let levels = checksum_layout(checksums.offset / PAGE_BYTES - 1).unwrap();
            let top = levels.last().unwrap();
            header.sections[Section::Checksums as usize].length =
                top.offset + top.page_count * PAGE_BYTES;
            let end = header.section(Section::Checksums);
            let encoded = header.encode().unwrap();
            assert_eq!(
                SharedHeader::decode(&encoded, end.offset + end.length).unwrap(),
                header
            );
            assert_eq!(header.row_bytes(Section::Cores), 4);
            assert_eq!(header.row_bytes(Section::CorePayloads), u64::from(width));
            assert_eq!(header.row_bytes(Section::Groups), 14);
            assert_eq!(header.row_bytes(Section::Members), 9);
            for invalid in [0, 12, 17, 24] {
                let mut bad = header.clone();
                bad.core_payload_bytes = invalid;
                assert!(bad.encode().is_err());
            }
            let mut bad = header.clone();
            bad.sections[Section::CorePrefixes as usize].length -= 4;
            assert!(bad.encode().is_err());
            let mut bad = header.clone();
            bad.sections[Section::CorePayloads as usize].length += 1;
            assert!(bad.encode().is_err());
            let mut bad = encoded;
            bad[8] = 4;
            assert!(SharedHeader::decode(&bad, end.offset + end.length).is_err());
        }
    }
}
