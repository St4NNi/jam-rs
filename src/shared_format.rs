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
    pub(crate) window: u16,
    pub(crate) core_count: u64,
    pub(crate) occurrence_count: u64,
    pub(crate) document_count: u32,
    pub(crate) contig_count: u32,
    pub(crate) source_bases: u64,
    pub(crate) manifest_sha256: [u8; 32],
    pub(crate) body_sha256: [u8; 32],
    pub(crate) checksum_root_sha256: [u8; 32],
    pub(crate) sections: [SectionRange; 10],
}

impl SharedHeader {
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
        put_u16(&mut out, 8, 1);
        put_u16(&mut out, 10, HEADER_BYTES as u16);
        put_u16(&mut out, 12, self.window);
        put_u64(&mut out, 16, self.core_count);
        put_u64(&mut out, 24, self.occurrence_count);
        put_u32(&mut out, 32, self.document_count);
        put_u32(&mut out, 36, self.contig_count);
        put_u64(&mut out, 40, self.source_bases);
        out[48..80].copy_from_slice(&self.manifest_sha256);
        out[80..112].copy_from_slice(&self.body_sha256);
        out[112..144].copy_from_slice(&self.checksum_root_sha256);
        for (index, section) in self.sections.iter().enumerate() {
            put_u64(&mut out, 144 + 16 * index, section.offset);
            put_u64(&mut out, 152 + 16 * index, section.length);
        }
        let digest = sha256(&out);
        out[HEADER_HASH_OFFSET..HEADER_HASH_OFFSET + 32].copy_from_slice(&digest);
        Ok(out)
    }

    pub(crate) fn decode(bytes: &[u8], file_bytes: u64) -> Result<Self, SharedError> {
        if bytes.len() != HEADER_BYTES
            || &bytes[..8] != b"JSHARED\0"
            || bytes[8..12] != [1, 0, 0, 16]
            || bytes[14..16] != [0, 0]
            || bytes[HEADER_HASH_OFFSET + 32..]
                .iter()
                .any(|&byte| byte != 0)
        {
            return Err(SharedError::Invalid("header"));
        }
        let mut checked = [0; HEADER_BYTES];
        checked.copy_from_slice(bytes);
        checked[HEADER_HASH_OFFSET..HEADER_HASH_OFFSET + 32].fill(0);
        if sha256(&checked) != bytes[HEADER_HASH_OFFSET..HEADER_HASH_OFFSET + 32] {
            return Err(SharedError::ChecksumMismatch);
        }
        let header = Self {
            window: u16::from_le_bytes(bytes[12..14].try_into().unwrap()),
            core_count: read_u64(bytes, 16),
            occurrence_count: read_u64(bytes, 24),
            document_count: read_u32(bytes, 32),
            contig_count: read_u32(bytes, 36),
            source_bases: read_u64(bytes, 40),
            manifest_sha256: bytes[48..80].try_into().unwrap(),
            body_sha256: bytes[80..112].try_into().unwrap(),
            checksum_root_sha256: bytes[112..144].try_into().unwrap(),
            sections: std::array::from_fn(|index| SectionRange {
                offset: read_u64(bytes, 144 + 16 * index),
                length: read_u64(bytes, 152 + 16 * index),
            }),
        };
        header.validate(file_bytes)?;
        Ok(header)
    }

    fn validate(&self, file_bytes: u64) -> Result<(), SharedError> {
        if self.window == 0
            || self.document_count == 0
            || self.contig_count == 0
            || self.core_count > self.occurrence_count
        {
            return Err(SharedError::Invalid("header counts"));
        }
        let mut previous = HEADER_BYTES as u64;
        for (kind, section) in Section::ALL.into_iter().zip(self.sections) {
            let expected = previous
                .checked_next_multiple_of(PAGE_BYTES)
                .ok_or(SharedError::Invalid("section padding"))?;
            if section.offset != expected || !section.length.is_multiple_of(kind.row_bytes()) {
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
            || self.core_count.checked_mul(CORE_ROW_BYTES)
                != Some(self.section(Section::Cores).length)
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
            offset = offset.next_multiple_of(PAGE_BYTES);
            let section = SectionRange {
                offset,
                length: lengths[index],
            };
            offset += section.length;
            section
        });
        SharedHeader {
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
}
