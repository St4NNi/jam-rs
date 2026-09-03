use crate::jidx::{
    CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_RECORD_SIZE, DocumentRecord, HEADER_SIZE, Header,
    JidxError, SectionKind,
};
use memmap2::{Mmap, MmapOptions};
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::path::Path;
use thiserror::Error;

pub use crate::jidx_postings::{SeedEntry, SeedOccurrence};

pub type MetagenomeId = u32;
pub type ContigId = u32;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Metagenome<'a> {
    pub id: MetagenomeId,
    pub name: &'a str,
    pub bgzf_uri: &'a str,
    pub fai_uri: &'a str,
    pub gzi_uri: &'a str,
    pub bgzf_bytes: u64,
    pub fai_bytes: u64,
    pub gzi_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub fai_sha256: [u8; 32],
    pub gzi_sha256: [u8; 32],
    pub contig_start: ContigId,
    pub contig_count: u32,
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
        let reader = Self { mmap, header };
        reader.validate_metadata()?;
        Ok(reader)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn verify_checksum(&self) -> Result<(), JidxReaderError> {
        self.header.verify_body(&self.mmap)?;
        crate::jidx_postings::validate_table(&self.mmap, &self.header)?;
        for index in 0..self.header.seed_count {
            let seed = crate::jidx_postings::entry(&self.mmap, &self.header, index)?;
            self.seed_occurrences(seed)?;
        }
        Ok(())
    }

    pub fn find_seed(&self, packed_key: u64) -> Result<Option<SeedEntry>, JidxReaderError> {
        Ok(crate::jidx_postings::lookup(
            &self.mmap,
            &self.header,
            packed_key,
        )?)
    }

    pub fn seed_metagenomes(&self, seed: SeedEntry) -> Result<Vec<MetagenomeId>, JidxReaderError> {
        Ok(crate::jidx_postings::documents(
            &self.mmap,
            &self.header,
            seed,
        )?)
    }

    pub fn seed_occurrences(
        &self,
        seed: SeedEntry,
    ) -> Result<Vec<SeedOccurrence>, JidxReaderError> {
        let occurrences = crate::jidx_postings::occurrences(&self.mmap, &self.header, seed)?;
        let mut metagenomes = Vec::new();
        for occurrence in &occurrences {
            let record = self.contig_record(occurrence.contig_id)?;
            if occurrence
                .position
                .checked_add(u64::from(self.header.k))
                .is_none_or(|end| end > record.length)
            {
                return Err(JidxError::Invalid("contig posting position").into());
            }
            if metagenomes.last().copied() != Some(record.document_id) {
                metagenomes.push(record.document_id);
            }
        }
        if metagenomes != self.seed_metagenomes(seed)? {
            return Err(JidxError::Invalid("seed document postings").into());
        }
        Ok(occurrences)
    }

    pub fn metagenome(&self, id: MetagenomeId) -> Result<Option<Metagenome<'_>>, JidxReaderError> {
        if id >= self.header.document_count {
            return Ok(None);
        }
        let strings = self.section_bytes(SectionKind::Strings)?;
        let record = DocumentRecord::decode(self.record_bytes(
            SectionKind::Documents,
            u64::from(id),
            DOCUMENT_RECORD_SIZE,
        )?)?;
        Ok(Some(Metagenome {
            id,
            name: record.name.resolve(strings)?,
            bgzf_uri: record.bgzf_uri.resolve(strings)?,
            fai_uri: record.fai_uri.resolve(strings)?,
            gzi_uri: record.gzi_uri.resolve(strings)?,
            bgzf_bytes: record.bgzf_bytes,
            fai_bytes: record.fai_bytes,
            gzi_bytes: record.gzi_bytes,
            bgzf_sha256: record.bgzf_sha256,
            fai_sha256: record.fai_sha256,
            gzi_sha256: record.gzi_sha256,
            contig_start: record.contig_start,
            contig_count: record.contig_count,
        }))
    }

    pub fn contig(&self, id: ContigId) -> Result<Option<Contig<'_>>, JidxReaderError> {
        if id >= self.header.contig_count {
            return Ok(None);
        }
        let strings = self.section_bytes(SectionKind::Strings)?;
        let record = self.contig_record(id)?;
        Ok(Some(Contig {
            id,
            metagenome_id: record.document_id,
            name: record.name.resolve(strings)?,
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

    fn section_bytes(&self, kind: SectionKind) -> Result<&[u8], JidxReaderError> {
        let section = self.header.section(kind);
        let start =
            usize::try_from(section.offset).map_err(|_| JidxError::Invalid("section offset"))?;
        let length =
            usize::try_from(section.length).map_err(|_| JidxError::Invalid("section length"))?;
        self.mmap
            .get(start..start + length)
            .ok_or_else(|| JidxError::Invalid("section range").into())
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
        let start = usize::try_from(start).map_err(|_| JidxError::Invalid("record offset"))?;
        let end = usize::try_from(end).map_err(|_| JidxError::Invalid("record range"))?;
        self.mmap
            .get(start..end)
            .ok_or_else(|| JidxError::Invalid("record range").into())
    }

    fn validate_metadata(&self) -> Result<(), JidxReaderError> {
        let mut expected_contig = 0u32;
        let mut ranges = Vec::with_capacity(self.header.document_count as usize);
        let mut metagenome_names = HashSet::new();
        for id in 0..self.header.document_count {
            let metagenome = self
                .metagenome(id)?
                .ok_or(JidxError::Invalid("missing metagenome"))?;
            if metagenome.contig_start != expected_contig
                || metagenome.contig_count == 0
                || metagenome.bgzf_bytes == 0
                || metagenome.fai_bytes == 0
                || metagenome.gzi_bytes == 0
                || metagenome.bgzf_sha256 == [0; 32]
                || metagenome.fai_sha256 == [0; 32]
                || metagenome.gzi_sha256 == [0; 32]
                || !metagenome_names.insert(metagenome.name)
            {
                return Err(JidxError::Invalid("metagenome metadata").into());
            }
            expected_contig = expected_contig
                .checked_add(metagenome.contig_count)
                .ok_or(JidxError::Invalid("contig range"))?;
            if expected_contig > self.header.contig_count {
                return Err(JidxError::Invalid("contig range").into());
            }
            ranges.push((metagenome.contig_start, expected_contig));
        }
        if expected_contig != self.header.contig_count {
            return Err(JidxError::Invalid("contig coverage").into());
        }

        let mut contig_names = HashSet::new();
        for id in 0..self.header.contig_count {
            let contig = self
                .contig(id)?
                .ok_or(JidxError::Invalid("missing contig"))?;
            let (start, end) = ranges
                .get(contig.metagenome_id as usize)
                .copied()
                .ok_or(JidxError::Invalid("contig metagenome ID"))?;
            if id < start
                || id >= end
                || contig.length == 0
                || contig.line_bases == 0
                || contig.line_width < contig.line_bases
                || contig.line_width > contig.line_bases.saturating_add(2)
                || !contig_names.insert((contig.metagenome_id, contig.name))
            {
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
        CONTIG_POSTING_SIZE, DOCUMENT_POSTING_SIZE, FilterKind, PostingCodec, SEED_RECORD_SIZE,
        SectionDescriptor, SectionKind, SeedScheme, sha256,
    };
    use std::path::PathBuf;

    fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
        bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
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
        let fai = push_string(&mut strings, "seq.bgz.fai");
        let gzi = push_string(&mut strings, "seq.bgz.gzi");
        let contig_name = push_string(&mut strings, "contig");

        let mut document = vec![0; DOCUMENT_RECORD_SIZE as usize];
        put_string(&mut document, 0, name);
        put_string(&mut document, 8, bgzf);
        put_string(&mut document, 16, fai);
        put_string(&mut document, 24, gzi);
        put_u64(&mut document, 32, 100);
        put_u64(&mut document, 40, 20);
        put_u64(&mut document, 48, 16);
        put_u32(&mut document, 56, 0);
        put_u32(&mut document, 60, 1);
        document[64..96].fill(3);
        document[96..128].fill(4);
        document[128..160].fill(5);

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
        ];
        let sizes = [
            0,
            DOCUMENT_RECORD_SIZE,
            CONTIG_RECORD_SIZE,
            SEED_RECORD_SIZE,
            DOCUMENT_POSTING_SIZE,
            CONTIG_POSTING_SIZE,
        ];
        let mut offset = HEADER_SIZE as u64;
        let sections = std::array::from_fn(|index| {
            offset = offset.next_multiple_of(8);
            let descriptor = SectionDescriptor {
                kind: SectionKind::ALL[index],
                record_size: sizes[index],
                offset,
                length: payloads[index].len() as u64,
            };
            offset += descriptor.length;
            descriptor
        });
        let mut file = vec![0; offset as usize];
        for (section, payload) in sections.iter().zip(&payloads) {
            let start = section.offset as usize;
            file[start..start + payload.len()].copy_from_slice(payload);
        }
        let header = Header {
            k: 21,
            seed_scheme: SeedScheme::WindowMinHash,
            posting_codec: PostingCodec::Raw,
            filter: FilterKind::None,
            document_count: 1,
            contig_count: 1,
            seed_count: 1,
            occurrence_count: 1,
            segment_bases: 256,
            seeds_per_segment: 2,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
            body_sha256: sha256(&file[HEADER_SIZE..]),
            sections,
        };
        file[..HEADER_SIZE].copy_from_slice(&header.encode().unwrap());
        if corrupt_padding {
            file[sections[0].offset as usize + sections[0].length as usize] ^= 1;
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
    fn rejects_invalid_contig_ownership() {
        let (_directory, path) = fixture(true, false, false, false);
        assert!(JidxReader::open(path).is_err());
    }

    #[test]
    fn checksum_verification_is_explicit() {
        let (_directory, path) = fixture(false, true, false, false);
        let reader = JidxReader::open(path).unwrap();
        assert!(reader.verify_checksum().is_err());
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
