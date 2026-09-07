use crate::jidx::{
    CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_RECORD_SIZE, DocumentRecord, HEADER_SIZE, Header,
    JidxError, PAGE_SIZE, SectionKind, StringRef, sha256,
};
use memmap2::{Mmap, MmapOptions};
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Mutex;
use thiserror::Error;

pub use crate::jidx_postings::{SeedEntry, SeedOccurrence};

pub type MetagenomeId = u32;
pub type ContigId = u32;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Metagenome<'a> {
    pub id: MetagenomeId,
    pub name: &'a str,
    pub bgzf_uri: &'a str,
    pub bgzf_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub contig_start: ContigId,
    pub contig_count: u32,
    pub gzi: &'a [u8],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Contig<'a> {
    pub id: ContigId,
    pub metagenome_id: MetagenomeId,
    pub name: &'a str,
    pub length: u64,
    pub fasta_offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
}

pub struct JidxReader {
    mmap: Mmap,
    header: Header,
    verified_pages: Mutex<HashSet<u64>>,
}

impl JidxReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, JidxReaderError> {
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        if file_len < HEADER_SIZE as u64 {
            return Err(JidxError::FileTooSmall {
                expected: HEADER_SIZE,
                actual: usize::try_from(file_len).unwrap_or(usize::MAX),
            }
            .into());
        }
        // SAFETY: the mapping is read-only and retained by the reader for all returned borrows.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let header = Header::decode_header(&mmap[..HEADER_SIZE], file_len)?;
        let reader = Self {
            mmap,
            header,
            verified_pages: Mutex::new(HashSet::new()),
        };
        reader.validate_documents()?;
        Ok(reader)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn verify_checksum(&self) -> Result<(), JidxReaderError> {
        self.header.verify_body(&self.mmap)?;
        self.validate_contigs()?;
        crate::jidx_postings::validate_table(self)?;
        for index in 0..self.header.seed_count {
            let seed = crate::jidx_postings::entry(self, index)?;
            self.seed_occurrences(seed)?;
        }
        Ok(())
    }

    pub fn find_seed(&self, packed_key: u64) -> Result<Option<SeedEntry>, JidxReaderError> {
        Ok(crate::jidx_postings::lookup(self, packed_key)?)
    }

    pub fn seed_metagenomes(&self, seed: SeedEntry) -> Result<Vec<MetagenomeId>, JidxReaderError> {
        Ok(crate::jidx_postings::documents(self, seed)?)
    }

    pub fn seed_occurrences(
        &self,
        seed: SeedEntry,
    ) -> Result<Vec<SeedOccurrence>, JidxReaderError> {
        let occurrences = crate::jidx_postings::occurrences(self, seed)?;
        let mut metagenomes = Vec::new();
        for occurrence in &occurrences {
            let contig = self
                .contig(occurrence.contig_id)?
                .ok_or(JidxError::Invalid("missing occurrence contig"))?;
            if occurrence
                .position
                .checked_add(u64::from(self.header.k))
                .is_none_or(|end| end > contig.length)
            {
                return Err(JidxError::Invalid("contig posting position").into());
            }
            if metagenomes.last().copied() != Some(contig.metagenome_id) {
                metagenomes.push(contig.metagenome_id);
            }
        }
        if metagenomes != self.seed_metagenomes(seed)? {
            return Err(JidxError::Invalid("seed document postings").into());
        }
        Ok(occurrences)
    }

    pub(crate) fn metagenome_name(
        &self,
        id: MetagenomeId,
    ) -> Result<Option<&str>, JidxReaderError> {
        if id >= self.header.document_count {
            return Ok(None);
        }
        Ok(Some(self.resolve_string(self.document_record(id)?.name)?))
    }

    pub fn metagenome(&self, id: MetagenomeId) -> Result<Option<Metagenome<'_>>, JidxReaderError> {
        if id >= self.header.document_count {
            return Ok(None);
        }
        let record = self.document_record(id)?;
        let gzi = self.gzi_bytes(record)?;
        Ok(Some(Metagenome {
            id,
            name: self.resolve_string(record.name)?,
            bgzf_uri: self.resolve_string(record.bgzf_uri)?,
            bgzf_bytes: record.bgzf_bytes,
            bgzf_sha256: record.bgzf_sha256,
            contig_start: record.contig_start,
            contig_count: record.contig_count,
            gzi,
        }))
    }

    pub fn contig(&self, id: ContigId) -> Result<Option<Contig<'_>>, JidxReaderError> {
        if id >= self.header.contig_count {
            return Ok(None);
        }
        let record = self.contig_record(id)?;
        let document = self.document_record(record.document_id)?;
        let document_end = document
            .contig_start
            .checked_add(document.contig_count)
            .ok_or(JidxError::Invalid("contig range"))?;
        if id < document.contig_start
            || id >= document_end
            || record.length == 0
            || record.line_bases == 0
            || record.line_width < record.line_bases
            || record.line_width > record.line_bases.saturating_add(2)
        {
            return Err(JidxError::Invalid("contig metadata").into());
        }
        Ok(Some(Contig {
            id,
            metagenome_id: record.document_id,
            name: self.resolve_string(record.name)?,
            length: record.length,
            fasta_offset: record.fasta_offset,
            line_bases: record.line_bases,
            line_width: record.line_width,
        }))
    }

    fn contig_record(&self, id: ContigId) -> Result<ContigRecord, JidxReaderError> {
        Ok(ContigRecord::decode(self.record_bytes(
            SectionKind::Contigs,
            u64::from(id),
            CONTIG_RECORD_SIZE,
        )?)?)
    }

    fn document_record(&self, id: MetagenomeId) -> Result<DocumentRecord, JidxReaderError> {
        if id >= self.header.document_count {
            return Err(JidxError::Invalid("metagenome ID").into());
        }
        Ok(DocumentRecord::decode(self.record_bytes(
            SectionKind::Documents,
            u64::from(id),
            DOCUMENT_RECORD_SIZE,
        )?)?)
    }

    fn resolve_string(&self, reference: StringRef) -> Result<&str, JidxReaderError> {
        let section = self.header.section(SectionKind::Strings);
        let start = section
            .offset
            .checked_add(u64::from(reference.offset))
            .ok_or(JidxError::Invalid("string offset"))?;
        let end = start
            .checked_add(u64::from(reference.length))
            .ok_or(JidxError::Invalid("string range"))?;
        let section_end = section
            .offset
            .checked_add(section.length)
            .ok_or(JidxError::Invalid("string range"))?;
        if end > section_end {
            return Err(JidxError::Invalid("string range").into());
        }
        let bytes = self.checked_bytes(start, end)?;
        if bytes.is_empty() || bytes.iter().any(|byte| matches!(byte, 0 | b'\n' | b'\r')) {
            return Err(JidxError::Invalid("string value").into());
        }
        Ok(std::str::from_utf8(bytes).map_err(|_| JidxError::Invalid("string encoding"))?)
    }

    fn gzi_bytes(&self, record: DocumentRecord) -> Result<&[u8], JidxReaderError> {
        let section = self.header.section(SectionKind::Gzi);
        let start = section
            .offset
            .checked_add(record.gzi_offset)
            .ok_or(JidxError::Invalid("GZI offset"))?;
        let end = start
            .checked_add(record.gzi_length)
            .ok_or(JidxError::Invalid("GZI range"))?;
        let section_end = section
            .offset
            .checked_add(section.length)
            .ok_or(JidxError::Invalid("GZI range"))?;
        if end > section_end {
            return Err(JidxError::Invalid("GZI range").into());
        }
        Ok(self.checked_bytes(start, end)?)
    }

    pub(crate) fn checked_bytes(&self, start: u64, end: u64) -> Result<&[u8], JidxError> {
        let checksums = self.header.section(SectionKind::BlockChecksums);
        if start > end || start < PAGE_SIZE || end > checksums.offset {
            return Err(JidxError::Invalid("checked range"));
        }
        let start_index =
            usize::try_from(start).map_err(|_| JidxError::Invalid("checked range"))?;
        let end_index = usize::try_from(end).map_err(|_| JidxError::Invalid("checked range"))?;
        let bytes = self
            .mmap
            .get(start_index..end_index)
            .ok_or(JidxError::Invalid("checked range"))?;
        if start == end {
            return Ok(bytes);
        }

        let first_page = start / PAGE_SIZE;
        let last_page = (end - 1) / PAGE_SIZE;
        for page in first_page..=last_page {
            let cached = self
                .verified_pages
                .lock()
                .map_err(|_| JidxError::Invalid("page checksum cache"))?
                .contains(&page);
            if cached {
                continue;
            }

            let page_start = page
                .checked_mul(PAGE_SIZE)
                .ok_or(JidxError::Invalid("page range"))?;
            let page_end = page_start
                .checked_add(PAGE_SIZE)
                .ok_or(JidxError::Invalid("page range"))?;
            if page_end > checksums.offset {
                return Err(JidxError::Invalid("page range"));
            }
            let page_start_index =
                usize::try_from(page_start).map_err(|_| JidxError::Invalid("page range"))?;
            let page_end_index =
                usize::try_from(page_end).map_err(|_| JidxError::Invalid("page range"))?;
            let actual = sha256(
                self.mmap
                    .get(page_start_index..page_end_index)
                    .ok_or(JidxError::Invalid("page range"))?,
            );

            let checksum_index = page
                .checked_sub(1)
                .and_then(|page| page.checked_mul(32))
                .ok_or(JidxError::Invalid("page checksum"))?;
            let checksum_start = checksums
                .offset
                .checked_add(checksum_index)
                .ok_or(JidxError::Invalid("page checksum"))?;
            let checksum_end = checksum_start
                .checked_add(32)
                .ok_or(JidxError::Invalid("page checksum"))?;
            let checksum_start =
                usize::try_from(checksum_start).map_err(|_| JidxError::Invalid("page checksum"))?;
            let checksum_end =
                usize::try_from(checksum_end).map_err(|_| JidxError::Invalid("page checksum"))?;
            if self.mmap.get(checksum_start..checksum_end) != Some(actual.as_slice()) {
                return Err(JidxError::ChecksumMismatch);
            }

            let mut cache = self
                .verified_pages
                .lock()
                .map_err(|_| JidxError::Invalid("page checksum cache"))?;
            if cache.len() >= 4096 {
                cache.clear();
            }
            cache.insert(page);
        }
        Ok(bytes)
    }

    fn record_bytes(
        &self,
        kind: SectionKind,
        index: u64,
        size: u32,
    ) -> Result<&[u8], JidxReaderError> {
        let section = self.header.section(kind);
        let relative = index
            .checked_mul(u64::from(size))
            .ok_or(JidxError::Invalid("record offset"))?;
        let start = section
            .offset
            .checked_add(relative)
            .ok_or(JidxError::Invalid("record offset"))?;
        let end = start
            .checked_add(u64::from(size))
            .ok_or(JidxError::Invalid("record range"))?;
        Ok(self.checked_bytes(start, end)?)
    }

    fn validate_documents(&self) -> Result<(), JidxReaderError> {
        let mut expected_contig = 0u32;
        let mut expected_gzi = 0u64;
        let mut metagenome_names = HashSet::new();
        for id in 0..self.header.document_count {
            let record = self.document_record(id)?;
            let name = self.resolve_string(record.name)?;
            self.resolve_string(record.bgzf_uri)?;
            let gzi_end = record
                .gzi_offset
                .checked_add(record.gzi_length)
                .ok_or(JidxError::Invalid("GZI range"))?;
            if record.contig_start != expected_contig
                || record.contig_count == 0
                || record.bgzf_bytes == 0
                || record.bgzf_sha256 == [0; 32]
                || record.gzi_offset != expected_gzi
                || record.gzi_length < 8
                || !metagenome_names.insert(name)
            {
                return Err(JidxError::Invalid("metagenome metadata").into());
            }
            expected_contig = expected_contig
                .checked_add(record.contig_count)
                .ok_or(JidxError::Invalid("contig range"))?;
            if expected_contig > self.header.contig_count {
                return Err(JidxError::Invalid("contig range").into());
            }
            expected_gzi = gzi_end;
        }
        if expected_contig != self.header.contig_count {
            return Err(JidxError::Invalid("contig coverage").into());
        }
        if expected_gzi != self.header.section(SectionKind::Gzi).length {
            return Err(JidxError::Invalid("GZI coverage").into());
        }
        Ok(())
    }

    fn validate_contigs(&self) -> Result<(), JidxReaderError> {
        let mut contig_names = HashSet::new();
        for id in 0..self.header.contig_count {
            let contig = self
                .contig(id)?
                .ok_or(JidxError::Invalid("missing contig"))?;
            if !contig_names.insert((contig.metagenome_id, contig.name)) {
                return Err(JidxError::Invalid("contig metadata").into());
            }
        }
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum JidxReaderError {
    #[error("JIDX I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Format(#[from] JidxError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx::{
        CONTIG_POSTING_SIZE, DOCUMENT_POSTING_SIZE, SECTION_COUNT, SEED_RECORD_SIZE,
        SectionDescriptor, VERSION, sha256,
    };
    use std::io::{Seek, SeekFrom, Write};
    use std::path::PathBuf;

    fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
        bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
    }

    fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
        bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
    }

    fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
        bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
    }

    fn push_string(strings: &mut Vec<u8>, value: &str) -> (u32, u32) {
        let offset = strings.len() as u32;
        strings.extend_from_slice(value.as_bytes());
        (offset, value.len() as u32)
    }

    fn put_string(bytes: &mut [u8], offset: usize, value: (u32, u32)) {
        put_u32(bytes, offset, value.0);
        put_u32(bytes, offset + 4, value.1);
    }

    fn fixture(
        bad_metagenome_id: bool,
        corrupt_padding: bool,
        bad_document_posting: bool,
        bad_occurrence: bool,
    ) -> (tempfile::TempDir, PathBuf) {
        let mut strings = Vec::new();
        let name = push_string(&mut strings, "doc");
        let bgzf = push_string(&mut strings, "seq.bgz");
        let contig_name = push_string(&mut strings, "contig");

        let mut document = vec![0; DOCUMENT_RECORD_SIZE as usize];
        put_string(&mut document, 0, name);
        put_string(&mut document, 8, bgzf);
        put_u64(&mut document, 16, 100);
        put_u32(&mut document, 24, 0);
        put_u32(&mut document, 28, 1);
        document[32..64].fill(3);
        put_u64(&mut document, 64, 0);
        put_u64(&mut document, 72, 8);

        let mut contig = vec![0; CONTIG_RECORD_SIZE as usize];
        put_u32(&mut contig, 0, u32::from(bad_metagenome_id));
        put_string(&mut contig, 4, contig_name);
        put_u64(&mut contig, 16, 100);
        put_u64(&mut contig, 24, 5);
        put_u32(&mut contig, 32, 100);
        put_u32(&mut contig, 36, 101);

        let mut seed = vec![0; SEED_RECORD_SIZE as usize];
        put_u64(&mut seed, 0, 0x1234);
        put_u64(&mut seed, 8, 0);
        put_u32(&mut seed, 16, 1);
        put_u64(&mut seed, 24, 0);
        put_u64(&mut seed, 32, 1);
        let document_posting = u32::from(bad_document_posting).to_le_bytes().to_vec();
        let mut contig_posting = vec![0; CONTIG_POSTING_SIZE as usize];
        contig_posting[4] = 1;
        put_u64(&mut contig_posting, 8, if bad_occurrence { 90 } else { 2 });
        let payloads = [
            strings,
            document,
            contig,
            seed,
            document_posting,
            contig_posting,
            0u64.to_le_bytes().to_vec(),
        ];
        let sizes = [
            0,
            DOCUMENT_RECORD_SIZE,
            CONTIG_RECORD_SIZE,
            SEED_RECORD_SIZE,
            DOCUMENT_POSTING_SIZE,
            CONTIG_POSTING_SIZE,
            0,
            32,
        ];
        let mut offset = PAGE_SIZE;
        let mut sections = Vec::with_capacity(SECTION_COUNT);
        for (index, payload) in payloads.iter().enumerate() {
            let descriptor = SectionDescriptor {
                kind: SectionKind::ALL[index],
                record_size: sizes[index],
                offset,
                length: payload.len() as u64,
            };
            sections.push(descriptor);
            offset = descriptor
                .offset
                .checked_add(descriptor.length)
                .unwrap()
                .next_multiple_of(PAGE_SIZE);
        }
        let checksum_count = offset / PAGE_SIZE - 1;
        sections.push(SectionDescriptor {
            kind: SectionKind::BlockChecksums,
            record_size: 32,
            offset,
            length: checksum_count * 32,
        });
        let sections: [SectionDescriptor; SECTION_COUNT] = sections.try_into().unwrap();
        let mut file = vec![0; (offset + checksum_count * 32) as usize];
        for (section, payload) in sections.iter().zip(&payloads) {
            let start = section.offset as usize;
            file[start..start + payload.len()].copy_from_slice(payload);
        }
        for page in 1..=checksum_count {
            let start = (page * PAGE_SIZE) as usize;
            let digest = sha256(&file[start..start + PAGE_SIZE as usize]);
            let checksum_start = (offset + (page - 1) * 32) as usize;
            file[checksum_start..checksum_start + 32].copy_from_slice(&digest);
        }

        let body_sha256 = sha256(&file[HEADER_SIZE..]);
        let mut header = [0; HEADER_SIZE];
        header[..8].copy_from_slice(&crate::jidx::MAGIC);
        put_u16(&mut header, 8, VERSION);
        put_u16(&mut header, 10, HEADER_SIZE as u16);
        header[16] = 21;
        header[17] = 2;
        header[18] = 1;
        put_u16(&mut header, 20, SECTION_COUNT as u16);
        put_u32(&mut header, 24, 1);
        put_u32(&mut header, 28, 1);
        put_u64(&mut header, 32, 1);
        put_u64(&mut header, 40, 1);
        header[48..80].fill(1);
        header[80..112].fill(2);
        header[112..144].copy_from_slice(&body_sha256);
        for (index, section) in sections.iter().enumerate() {
            let start = 144 + index * 24;
            put_u16(&mut header, start, section.kind as u16);
            put_u32(&mut header, start + 4, section.record_size);
            put_u64(&mut header, start + 8, section.offset);
            put_u64(&mut header, start + 16, section.length);
        }
        put_u16(&mut header, 336, 16);
        file[..HEADER_SIZE].copy_from_slice(&header);
        if corrupt_padding {
            file[sections[2].offset as usize + sections[2].length as usize] ^= 1;
        }
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("fixture.jidx");
        std::fs::write(&path, file).unwrap();
        (directory, path)
    }

    #[test]
    fn reads_compact_metagenome_and_contig_ids() {
        let (_directory, path) = fixture(false, false, false, false);
        let reader = JidxReader::open(path).unwrap();
        reader.verify_checksum().unwrap();
        assert_eq!(reader.header().document_count, 1);
        let metagenome = reader.metagenome(0).unwrap().unwrap();
        assert_eq!((metagenome.id, metagenome.name), (0, "doc"));
        assert_eq!((metagenome.contig_start, metagenome.contig_count), (0, 1));
        let contig = reader.contig(0).unwrap().unwrap();
        assert_eq!(
            (contig.id, contig.metagenome_id, contig.name),
            (0, 0, "contig")
        );
        assert!(reader.metagenome(1).unwrap().is_none());
        assert!(reader.contig(1).unwrap().is_none());
    }

    #[test]
    fn validates_contig_ownership_on_access_and_audit() {
        let (_directory, path) = fixture(true, false, false, false);
        let reader = JidxReader::open(path).unwrap();
        assert!(reader.contig(0).is_err());
        assert!(reader.verify_checksum().is_err());
    }

    #[test]
    fn checksum_verification_is_explicit() {
        let (_directory, path) = fixture(false, true, false, false);
        let reader = JidxReader::open(path).unwrap();
        assert!(reader.verify_checksum().is_err());
    }

    #[test]
    fn accessed_page_checksum_is_verified_lazily() {
        let (_directory, path) = fixture(false, false, false, false);
        let mut file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.seek(SeekFrom::Start(4 * PAGE_SIZE)).unwrap();
        file.write_all(&[0x35]).unwrap();
        drop(file);

        let reader = JidxReader::open(path).unwrap();
        assert!(matches!(
            reader.find_seed(0x1234),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
    }

    #[test]
    fn exact_packed_seed_lookup_rejects_partial_matches() {
        let (_directory, path) = fixture(false, false, false, false);
        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert_eq!((seed.packed_key, seed.document_frequency), (0x1234, 1));
        assert!(reader.find_seed(0x1233).unwrap().is_none());
        assert!(reader.find_seed(0x1235).unwrap().is_none());
    }

    #[test]
    fn document_postings_are_checked_and_decoded() {
        let (_directory, path) = fixture(false, false, false, false);
        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert_eq!(reader.seed_metagenomes(seed).unwrap(), [0]);

        let (_directory, path) = fixture(false, false, true, false);
        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert!(reader.seed_metagenomes(seed).is_err());
    }

    #[test]
    fn contig_postings_include_checked_position_and_orientation() {
        let (_directory, path) = fixture(false, false, false, false);
        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert_eq!(
            reader.seed_occurrences(seed).unwrap(),
            [SeedOccurrence {
                contig_id: 0,
                position: 2,
                canonical_orientation: true,
            }]
        );

        let (_directory, path) = fixture(false, false, false, true);
        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert!(reader.seed_occurrences(seed).is_err());
    }
}
