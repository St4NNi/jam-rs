use sha2::{Digest, Sha256};
use thiserror::Error;

pub const MAGIC: [u8; 8] = *b"JIDX\0\0\0\0";
pub const VERSION: u16 = 1;
pub const HEADER_SIZE: usize = 320;
pub const SECTION_COUNT: usize = 6;
pub const DOCUMENT_RECORD_SIZE: u32 = 160;
pub const CONTIG_RECORD_SIZE: u32 = 40;
pub const SEED_RECORD_SIZE: u32 = 40;
pub const DOCUMENT_POSTING_SIZE: u32 = 4;
pub const CONTIG_POSTING_SIZE: u32 = 16;

const SECTION_TABLE_OFFSET: usize = 144;
const SECTION_DESCRIPTOR_SIZE: usize = 24;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SeedScheme {
    CanonicalKmer,
}

impl SeedScheme {
    const fn code(self) -> u8 {
        match self {
            Self::CanonicalKmer => 1,
        }
    }

    fn from_code(code: u8) -> Result<Self, JidxError> {
        match code {
            1 => Ok(Self::CanonicalKmer),
            _ => Err(JidxError::Invalid("seed scheme")),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PostingCodec {
    Raw,
}

impl PostingCodec {
    const fn code(self) -> u8 {
        match self {
            Self::Raw => 1,
        }
    }

    fn from_code(code: u8) -> Result<Self, JidxError> {
        match code {
            1 => Ok(Self::Raw),
            _ => Err(JidxError::Invalid("posting codec")),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FilterKind {
    None,
}

impl FilterKind {
    const fn code(self) -> u8 {
        match self {
            Self::None => 0,
        }
    }

    fn from_code(code: u8) -> Result<Self, JidxError> {
        match code {
            0 => Ok(Self::None),
            _ => Err(JidxError::Invalid("filter kind")),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u16)]
pub enum SectionKind {
    Strings = 1,
    Documents = 2,
    Contigs = 3,
    Seeds = 4,
    DocumentPostings = 5,
    ContigPostings = 6,
}

impl SectionKind {
    pub const ALL: [Self; SECTION_COUNT] = [
        Self::Strings,
        Self::Documents,
        Self::Contigs,
        Self::Seeds,
        Self::DocumentPostings,
        Self::ContigPostings,
    ];

    const fn record_size(self) -> u32 {
        match self {
            Self::Strings => 0,
            Self::Documents => DOCUMENT_RECORD_SIZE,
            Self::Contigs => CONTIG_RECORD_SIZE,
            Self::Seeds => SEED_RECORD_SIZE,
            Self::DocumentPostings => DOCUMENT_POSTING_SIZE,
            Self::ContigPostings => CONTIG_POSTING_SIZE,
        }
    }

    fn from_code(code: u16) -> Result<Self, JidxError> {
        match code {
            1 => Ok(Self::Strings),
            2 => Ok(Self::Documents),
            3 => Ok(Self::Contigs),
            4 => Ok(Self::Seeds),
            5 => Ok(Self::DocumentPostings),
            6 => Ok(Self::ContigPostings),
            _ => Err(JidxError::Invalid("section kind")),
        }
    }

    const fn index(self) -> usize {
        self as usize - 1
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SectionDescriptor {
    pub kind: SectionKind,
    pub record_size: u32,
    pub offset: u64,
    pub length: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Header {
    pub k: u8,
    pub seed_scheme: SeedScheme,
    pub posting_codec: PostingCodec,
    pub filter: FilterKind,
    pub document_count: u32,
    pub contig_count: u32,
    pub seed_count: u64,
    pub occurrence_count: u64,
    pub jam_sha256: [u8; 32],
    pub manifest_sha256: [u8; 32],
    pub body_sha256: [u8; 32],
    pub sections: [SectionDescriptor; SECTION_COUNT],
}

impl Header {
    pub fn encode(&self) -> Result<[u8; HEADER_SIZE], JidxError> {
        self.validate_layout(self.file_len()?)?;
        let mut bytes = [0; HEADER_SIZE];
        bytes[..8].copy_from_slice(&MAGIC);
        put_u16(&mut bytes, 8, VERSION);
        put_u16(&mut bytes, 10, HEADER_SIZE as u16);
        bytes[16] = self.k;
        bytes[17] = self.seed_scheme.code();
        bytes[18] = self.posting_codec.code();
        bytes[19] = self.filter.code();
        put_u16(&mut bytes, 20, SECTION_COUNT as u16);
        put_u32(&mut bytes, 24, self.document_count);
        put_u32(&mut bytes, 28, self.contig_count);
        put_u64(&mut bytes, 32, self.seed_count);
        put_u64(&mut bytes, 40, self.occurrence_count);
        bytes[48..80].copy_from_slice(&self.jam_sha256);
        bytes[80..112].copy_from_slice(&self.manifest_sha256);
        bytes[112..144].copy_from_slice(&self.body_sha256);
        for (index, section) in self.sections.iter().enumerate() {
            let start = SECTION_TABLE_OFFSET + index * SECTION_DESCRIPTOR_SIZE;
            put_u16(&mut bytes, start, section.kind as u16);
            put_u32(&mut bytes, start + 4, section.record_size);
            put_u64(&mut bytes, start + 8, section.offset);
            put_u64(&mut bytes, start + 16, section.length);
        }
        Ok(bytes)
    }

    pub fn decode(file: &[u8]) -> Result<Self, JidxError> {
        if file.len() < HEADER_SIZE {
            return Err(JidxError::FileTooSmall {
                expected: HEADER_SIZE,
                actual: file.len(),
            });
        }
        let bytes = &file[..HEADER_SIZE];
        if bytes[..8] != MAGIC {
            return Err(JidxError::Invalid("magic"));
        }
        let version = read_u16(bytes, 8);
        if version != VERSION {
            return Err(JidxError::UnsupportedVersion(version));
        }
        if read_u16(bytes, 10) != HEADER_SIZE as u16 {
            return Err(JidxError::Invalid("header size"));
        }
        if read_u32(bytes, 12) != 0 || read_u16(bytes, 22) != 0 {
            return Err(JidxError::Invalid("header flags"));
        }
        if read_u16(bytes, 20) != SECTION_COUNT as u16 || bytes[288..].iter().any(|byte| *byte != 0)
        {
            return Err(JidxError::Invalid("header reservation"));
        }

        let mut sections = Vec::with_capacity(SECTION_COUNT);
        for index in 0..SECTION_COUNT {
            let start = SECTION_TABLE_OFFSET + index * SECTION_DESCRIPTOR_SIZE;
            if read_u16(bytes, start + 2) != 0 {
                return Err(JidxError::Invalid("section flags"));
            }
            sections.push(SectionDescriptor {
                kind: SectionKind::from_code(read_u16(bytes, start))?,
                record_size: read_u32(bytes, start + 4),
                offset: read_u64(bytes, start + 8),
                length: read_u64(bytes, start + 16),
            });
        }
        let header = Self {
            k: bytes[16],
            seed_scheme: SeedScheme::from_code(bytes[17])?,
            posting_codec: PostingCodec::from_code(bytes[18])?,
            filter: FilterKind::from_code(bytes[19])?,
            document_count: read_u32(bytes, 24),
            contig_count: read_u32(bytes, 28),
            seed_count: read_u64(bytes, 32),
            occurrence_count: read_u64(bytes, 40),
            jam_sha256: bytes[48..80].try_into().expect("JIDX jam digest"),
            manifest_sha256: bytes[80..112].try_into().expect("JIDX manifest digest"),
            body_sha256: bytes[112..144].try_into().expect("JIDX body digest"),
            sections: sections.try_into().expect("fixed JIDX section count"),
        };
        header.validate_layout(
            u64::try_from(file.len()).map_err(|_| JidxError::Invalid("file length"))?,
        )?;
        if sha256(&file[HEADER_SIZE..]) != header.body_sha256 {
            return Err(JidxError::ChecksumMismatch);
        }
        Ok(header)
    }

    pub fn section(&self, kind: SectionKind) -> SectionDescriptor {
        self.sections[kind.index()]
    }

    fn file_len(&self) -> Result<u64, JidxError> {
        let section = self.sections[SECTION_COUNT - 1];
        section
            .offset
            .checked_add(section.length)
            .ok_or(JidxError::Invalid("section range"))
    }

    fn validate_layout(&self, file_len: u64) -> Result<(), JidxError> {
        if !(1..=32).contains(&self.k) {
            return Err(JidxError::Invalid("k-mer size"));
        }
        if self.document_count == 0
            || self.jam_sha256 == [0; 32]
            || self.manifest_sha256 == [0; 32]
            || self.body_sha256 == [0; 32]
        {
            return Err(JidxError::Invalid("required metadata"));
        }
        let mut previous_end = HEADER_SIZE as u64;
        for (index, section) in self.sections.iter().enumerate() {
            let expected = SectionKind::ALL[index];
            if section.kind != expected || section.record_size != expected.record_size() {
                return Err(JidxError::Invalid("section descriptor"));
            }
            if section.offset < previous_end || section.offset % 8 != 0 {
                return Err(JidxError::Invalid("section order"));
            }
            previous_end = section
                .offset
                .checked_add(section.length)
                .ok_or(JidxError::Invalid("section range"))?;
            if previous_end > file_len {
                return Err(JidxError::Invalid("section range"));
            }
        }
        if previous_end != file_len {
            return Err(JidxError::Invalid("trailing bytes"));
        }
        self.expect_length(SectionKind::Documents, u64::from(self.document_count))?;
        self.expect_length(SectionKind::Contigs, u64::from(self.contig_count))?;
        self.expect_length(SectionKind::Seeds, self.seed_count)?;
        self.expect_length(SectionKind::ContigPostings, self.occurrence_count)?;
        let documents = self.section(SectionKind::DocumentPostings);
        if !documents
            .length
            .is_multiple_of(u64::from(DOCUMENT_POSTING_SIZE))
        {
            return Err(JidxError::Invalid("document postings length"));
        }
        Ok(())
    }

    fn expect_length(&self, kind: SectionKind, count: u64) -> Result<(), JidxError> {
        let section = self.section(kind);
        let expected = count
            .checked_mul(u64::from(section.record_size))
            .ok_or(JidxError::Invalid("section length"))?;
        if section.length != expected {
            return Err(JidxError::Invalid("section length"));
        }
        Ok(())
    }
}

pub fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

fn read_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().expect("JIDX u16"))
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("JIDX u32"))
}

fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("JIDX u64"))
}

fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

#[derive(Debug, Error)]
pub enum JidxError {
    #[error("JIDX file is too small: expected {expected} bytes, got {actual}")]
    FileTooSmall { expected: usize, actual: usize },
    #[error("unsupported JIDX version {0}")]
    UnsupportedVersion(u16),
    #[error("invalid JIDX {0}")]
    Invalid(&'static str),
    #[error("JIDX body checksum mismatch")]
    ChecksumMismatch,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (Header, Vec<u8>) {
        let lengths = [8, 160, 40, 40, 4, 16];
        let mut offset = HEADER_SIZE as u64;
        let sections = std::array::from_fn(|index| {
            offset = offset.next_multiple_of(8);
            let section = SectionDescriptor {
                kind: SectionKind::ALL[index],
                record_size: SectionKind::ALL[index].record_size(),
                offset,
                length: lengths[index],
            };
            offset += lengths[index];
            section
        });
        let mut file = vec![0; offset as usize];
        for (index, byte) in file[HEADER_SIZE..].iter_mut().enumerate() {
            *byte = index as u8;
        }
        let header = Header {
            k: 21,
            seed_scheme: SeedScheme::CanonicalKmer,
            posting_codec: PostingCodec::Raw,
            filter: FilterKind::None,
            document_count: 1,
            contig_count: 1,
            seed_count: 1,
            occurrence_count: 1,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
            body_sha256: sha256(&file[HEADER_SIZE..]),
            sections,
        };
        file[..HEADER_SIZE].copy_from_slice(&header.encode().unwrap());
        (header, file)
    }

    #[test]
    fn header_roundtrip_validates_all_sections() {
        let (header, file) = fixture();
        assert_eq!(Header::decode(&file).unwrap(), header);
    }

    #[test]
    fn body_corruption_fails_closed() {
        let (_, mut file) = fixture();
        *file.last_mut().unwrap() ^= 1;
        assert!(matches!(
            Header::decode(&file),
            Err(JidxError::ChecksumMismatch)
        ));
    }

    #[test]
    fn overlapping_sections_are_rejected() {
        let (mut header, _) = fixture();
        header.sections[1].offset = header.sections[0].offset;
        assert!(matches!(header.encode(), Err(JidxError::Invalid(_))));
    }
}
