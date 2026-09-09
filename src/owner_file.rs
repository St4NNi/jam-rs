use crate::jidx::sha256;
use crate::jidx_reader::{Contig, Metagenome};
use crate::owner_format::{
    BlockRecord, DocumentRecord, OWNER_BLOCK_SIZE, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE,
    OWNER_HEADER_SIZE, OWNER_PAGE_SIZE, OwnerDocument, OwnerHeader, OwnerReaderError, OwnerSection,
    OwnerSeed, checksum_layout, read_u32, read_u64,
};
use crate::owner_postings::{
    MAX_KEYS_PER_BLOCK, OwnerHotKey, OwnerHotMember, decode_member, parse_hot,
};
use memmap2::{Mmap, MmapOptions};
use std::fs::File;
use std::io;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

const MAX_DECODED_HOT_BYTES: usize = 64 * 1024 * 1024;
const MAX_HOT_CACHE_BYTES: usize = 256 * 1024 * 1024;
const MAX_DECODED_OCCURRENCE_BYTES: usize = 64 * 1024 * 1024;

pub(crate) struct OwnerFile {
    pub(crate) _file: File,
    pub(crate) mmap: Mmap,
    pub(crate) file_identity: Option<[u64; 7]>,
    pub(crate) header: OwnerHeader,
    verified_pages: Box<[AtomicU64]>,
    hot_cache: Mutex<HotCache>,
    hot_cache_bytes: Arc<AtomicUsize>,
}

#[derive(Default)]
struct HotCache {
    blocks: std::collections::BTreeMap<u64, Arc<Vec<OwnerHotKey>>>,
}

impl OwnerFile {
    pub(crate) fn open(
        path: PathBuf,
        hot_cache_bytes: Arc<AtomicUsize>,
    ) -> Result<Self, OwnerReaderError> {
        let file = File::open(path)?;
        let file_identity = file_identity(&file)?;
        let file_len = file.metadata()?.len();
        if file_len < OWNER_HEADER_SIZE as u64 {
            return Err(OwnerReaderError::Invalid("file size"));
        }
        // SAFETY: the retained read-only mapping requires immutable owner generation files.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let header = OwnerHeader::decode(&mmap[..OWNER_HEADER_SIZE], file_len)?;
        let page_count = header.section(OwnerSection::PageChecksums).offset / OWNER_PAGE_SIZE - 1;
        let verified_pages = (0..page_count.div_ceil(64))
            .map(|_| AtomicU64::new(0))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        let owner = Self {
            _file: file,
            mmap,
            file_identity,
            header,
            verified_pages,
            hot_cache: Mutex::new(HotCache::default()),
            hot_cache_bytes,
        };
        owner.verify_unchanged()?;
        owner.verify_checksum_root()?;
        owner.verify_unchanged()?;
        Ok(owner)
    }

    pub(crate) fn audit_metadata(&self) -> Result<(), OwnerReaderError> {
        let mut document_names = std::collections::HashSet::new();
        let mut expected_contig = 0u32;
        let mut expected_gzi = 0u64;
        for document_id in 0..self.header.document_count {
            let record = self.document_record(document_id)?;
            if record.bgzf_bytes == 0
                || record.bgzf_sha256 == [0; 32]
                || record.contig_count == 0
                || record.contig_start != expected_contig
                || record.gzi_offset != expected_gzi
                || !document_names.insert(self.string(record.name_offset, record.name_length)?)
            {
                return Err(OwnerReaderError::Invalid("document metadata"));
            }
            self.string(record.uri_offset, record.uri_length)?;
            let gzi =
                self.section_bytes(OwnerSection::Gzi, record.gzi_offset, record.gzi_length)?;
            let mut reader = noodles_bgzf::gzi::io::Reader::new(gzi);
            reader.read_index()?;
            expected_gzi = expected_gzi
                .checked_add(record.gzi_length)
                .ok_or(OwnerReaderError::Invalid("GZI metadata"))?;
            let contig_end = record
                .contig_start
                .checked_add(record.contig_count)
                .ok_or(OwnerReaderError::Invalid("contig metadata"))?;
            let mut contig_names = std::collections::HashSet::new();
            for contig_id in record.contig_start..contig_end {
                let contig = self
                    .contig(contig_id)?
                    .ok_or(OwnerReaderError::Invalid("contig metadata"))?;
                if contig.metagenome_id != document_id
                    || contig.length == 0
                    || contig.line_bases == 0
                    || contig.line_width < contig.line_bases
                    || contig.line_width > contig.line_bases.saturating_add(2)
                    || !contig_names.insert(contig.name)
                {
                    return Err(OwnerReaderError::Invalid("contig metadata"));
                }
            }
            expected_contig = contig_end;
        }
        if expected_contig != self.header.contig_count
            || expected_gzi != self.header.section(OwnerSection::Gzi).length
        {
            return Err(OwnerReaderError::Invalid("metadata coverage"));
        }
        Ok(())
    }

    pub(crate) fn document_record(&self, id: u32) -> Result<DocumentRecord, OwnerReaderError> {
        if id >= self.header.document_count {
            return Err(OwnerReaderError::Invalid("document ID"));
        }
        let section = self.header.section(OwnerSection::Documents);
        let start = section.offset + u64::from(id) * u64::from(OWNER_DOCUMENT_SIZE);
        DocumentRecord::decode(self.checked_bytes(start, start + u64::from(OWNER_DOCUMENT_SIZE))?)
    }

    pub(crate) fn metagenome(&self, id: u32) -> Result<Option<Metagenome<'_>>, OwnerReaderError> {
        if id >= self.header.document_count {
            return Ok(None);
        }
        let record = self.document_record(id)?;
        Ok(Some(Metagenome {
            id,
            name: self.string(record.name_offset, record.name_length)?,
            bgzf_uri: self.string(record.uri_offset, record.uri_length)?,
            bgzf_bytes: record.bgzf_bytes,
            bgzf_sha256: record.bgzf_sha256,
            contig_start: record.contig_start,
            contig_count: record.contig_count,
            gzi: self.section_bytes(OwnerSection::Gzi, record.gzi_offset, record.gzi_length)?,
        }))
    }

    pub(crate) fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, OwnerReaderError> {
        if id >= self.header.contig_count {
            return Ok(None);
        }
        let section = self.header.section(OwnerSection::Contigs);
        let start = section.offset + u64::from(id) * u64::from(OWNER_CONTIG_SIZE);
        let bytes = self.checked_bytes(start, start + u64::from(OWNER_CONTIG_SIZE))?;
        let document = read_u32(bytes, 0);
        let owner = self.document_record(document)?;
        let contig_end = owner
            .contig_start
            .checked_add(owner.contig_count)
            .ok_or(OwnerReaderError::Invalid("contig ownership"))?;
        if id < owner.contig_start || id >= contig_end {
            return Err(OwnerReaderError::Invalid("contig ownership"));
        }
        Ok(Some(Contig {
            id,
            metagenome_id: document,
            name: self.string(read_u32(bytes, 4), read_u32(bytes, 8))?,
            length: read_u64(bytes, 16),
            fasta_offset: read_u64(bytes, 24),
            line_bases: read_u32(bytes, 32),
            line_width: read_u32(bytes, 36),
        }))
    }

    fn string(&self, offset: u32, length: u32) -> Result<&str, OwnerReaderError> {
        std::str::from_utf8(self.section_bytes(
            OwnerSection::Strings,
            u64::from(offset),
            u64::from(length),
        )?)
        .map_err(|_| OwnerReaderError::Invalid("metadata string"))
    }

    fn section_bytes(
        &self,
        kind: OwnerSection,
        offset: u64,
        length: u64,
    ) -> Result<&[u8], OwnerReaderError> {
        let section = self.header.section(kind);
        let start = section
            .offset
            .checked_add(offset)
            .ok_or(OwnerReaderError::Invalid("section range"))?;
        let end = start
            .checked_add(length)
            .ok_or(OwnerReaderError::Invalid("section range"))?;
        if end > section.offset + section.length {
            return Err(OwnerReaderError::Invalid("section range"));
        }
        self.checked_bytes(start, end)
    }

    fn checked_bytes(&self, start: u64, end: u64) -> Result<&[u8], OwnerReaderError> {
        self.verify_unchanged()?;
        if start > end || end > self.mmap.len() as u64 {
            return Err(OwnerReaderError::Invalid("file range"));
        }
        if start < self.header.section(OwnerSection::PageChecksums).offset {
            let first = start.max(OWNER_PAGE_SIZE) / OWNER_PAGE_SIZE;
            let last = end.saturating_sub(1) / OWNER_PAGE_SIZE;
            for page in first..=last {
                self.verify_page(page)?;
            }
        }
        Ok(&self.mmap[start as usize..end as usize])
    }

    fn verify_page(&self, page: u64) -> Result<(), OwnerReaderError> {
        if page == 0 {
            return Ok(());
        }
        let index = page - 1;
        let word = &self.verified_pages[index as usize / 64];
        let mask = 1u64 << (index % 64);
        if self.file_identity.is_some() && word.load(Ordering::Acquire) & mask != 0 {
            return Ok(());
        }
        let start = page * OWNER_PAGE_SIZE;
        let checksum_section = self.header.section(OwnerSection::PageChecksums);
        let expected_start = checksum_section.offset + index * 32;
        let expected = &self.mmap[expected_start as usize..expected_start as usize + 32];
        if sha256(&self.mmap[start as usize..start as usize + OWNER_PAGE_SIZE as usize]) != expected
        {
            return Err(OwnerReaderError::ChecksumMismatch);
        }
        self.verify_checksum_chain(index)?;
        self.verify_unchanged()?;
        if self.file_identity.is_some() {
            word.fetch_or(mask, Ordering::Release);
        }
        Ok(())
    }

    fn verify_checksum_chain(&self, leaf_index: u64) -> Result<(), OwnerReaderError> {
        let checksums = self.header.section(OwnerSection::PageChecksums);
        let levels = checksum_layout(checksums.offset / OWNER_PAGE_SIZE - 1)?;
        if leaf_index >= levels[0].hash_count {
            return Err(OwnerReaderError::Invalid("checksum leaf"));
        }
        let mut hash_index = leaf_index;
        for (index, level) in levels.iter().enumerate() {
            let page_index = hash_index / (OWNER_PAGE_SIZE / 32);
            if page_index >= level.page_count {
                return Err(OwnerReaderError::Invalid("checksum chain"));
            }
            let start = checksums
                .offset
                .checked_add(level.offset)
                .and_then(|offset| {
                    page_index
                        .checked_mul(OWNER_PAGE_SIZE)
                        .and_then(|bytes| offset.checked_add(bytes))
                })
                .ok_or(OwnerReaderError::Invalid("checksum chain"))?;
            let digest =
                sha256(&self.mmap[start as usize..start as usize + OWNER_PAGE_SIZE as usize]);
            if let Some(next) = levels.get(index + 1) {
                let expected = checksums
                    .offset
                    .checked_add(next.offset)
                    .and_then(|offset| {
                        page_index
                            .checked_mul(32)
                            .and_then(|bytes| offset.checked_add(bytes))
                    })
                    .ok_or(OwnerReaderError::Invalid("checksum chain"))?;
                if self.mmap[expected as usize..expected as usize + 32] != digest {
                    return Err(OwnerReaderError::ChecksumMismatch);
                }
                hash_index = page_index;
            } else if self.header.checksum_root_sha256 != digest {
                return Err(OwnerReaderError::ChecksumMismatch);
            }
        }
        Ok(())
    }

    fn verify_checksum_root(&self) -> Result<(), OwnerReaderError> {
        let checksums = self.header.section(OwnerSection::PageChecksums);
        let last = checksum_layout(checksums.offset / OWNER_PAGE_SIZE - 1)?
            .pop()
            .expect("checksum layout is nonempty");
        let start = checksums.offset + last.offset;
        if sha256(&self.mmap[start as usize..start as usize + OWNER_PAGE_SIZE as usize])
            != self.header.checksum_root_sha256
        {
            return Err(OwnerReaderError::ChecksumMismatch);
        }
        Ok(())
    }

    pub(crate) fn verify_unchanged(&self) -> io::Result<()> {
        if let Some(expected) = self.file_identity
            && file_identity(&self._file)? != Some(expected)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "owner index changed after opening",
            ));
        }
        Ok(())
    }

    pub(crate) fn verify_checksum(&self) -> Result<(), OwnerReaderError> {
        self.verify_unchanged()?;
        self.verify_checksum_root()?;
        if sha256(&self.mmap[OWNER_HEADER_SIZE..]) != self.header.body_sha256 {
            return Err(OwnerReaderError::ChecksumMismatch);
        }
        self.verify_unchanged()?;
        Ok(())
    }
}

fn file_identity(file: &File) -> io::Result<Option<[u64; 7]>> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        let metadata = file.metadata()?;
        Ok(Some([
            metadata.dev(),
            metadata.ino(),
            metadata.len(),
            metadata.mtime() as u64,
            metadata.mtime_nsec() as u64,
            metadata.ctime() as u64,
            metadata.ctime_nsec() as u64,
        ]))
    }
    #[cfg(not(unix))]
    {
        let _ = file;
        Ok(None)
    }
}
