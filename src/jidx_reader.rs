use crate::jidx::{
    CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_RECORD_SIZE, DocumentRecord, HEADER_SIZE, Header,
    JidxError, PAGE_SIZE, SEED_RECORD_SIZE, SectionKind, StringRef, read_u32, read_u64,
    seed_length, sha256,
};
#[cfg(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64")))]
use memmap2::Advice;
use memmap2::{Mmap, MmapOptions};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use thiserror::Error;

pub use crate::jidx_postings::{SeedDocument, SeedEntry, SeedOccurrence};

pub type MetagenomeId = u32;
pub type ContigId = u32;

pub(crate) const SEED_LOOKUP_BATCH_KEYS: usize = 32_768;
const EXACT_BLOCK_RECORDS: u64 = 512;
const EXACT_BLOCK_HEADER_SIZE: usize = 4096;
const EXACT_BLOCK_MAX_BYTES: u64 = 64 * 1024 * 1024;

pub(crate) struct ExactBlockFence {
    first_keys: Box<[u64]>,
    source_file_len: u64,
    seeds_offset: u64,
    seeds_length: u64,
    seed_count: u64,
    source_header_sha256: [u8; 32],
    source_body_sha256: [u8; 32],
}

impl ExactBlockFence {
    pub(crate) fn block_for(&self, key: u64) -> Option<u64> {
        self.first_keys
            .partition_point(|first| *first <= key)
            .checked_sub(1)
            .and_then(|block| u64::try_from(block).ok())
    }

    pub(crate) fn first_key(&self, block: u64) -> Option<u64> {
        self.first_keys.get(usize::try_from(block).ok()?).copied()
    }
}

struct ExactBlockFenceCache {
    path: PathBuf,
    sha256: [u8; 32],
    fence: Arc<ExactBlockFence>,
}

static EXACT_BLOCK_FENCE_CACHE: OnceLock<Mutex<Option<ExactBlockFenceCache>>> = OnceLock::new();

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FilterCacheIdentity {
    pub(crate) file: [u64; 7],
    pub(crate) header_sha256: [u8; 32],
    pub(crate) body_sha256: [u8; 32],
    pub(crate) section_offset: u64,
    pub(crate) section_length: u64,
}

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
    file: File,
    mmap: Mmap,
    header: Header,
    verified_pages: Box<[AtomicU64]>,
    filter_directory: Option<crate::jidx_filters::FilterDirectory>,
    exact_block_fence: Option<Arc<ExactBlockFence>>,
    filter_cache_identity: OnceLock<Option<FilterCacheIdentity>>,
    resident_filter_body: OnceLock<Option<Arc<Vec<u8>>>>,
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
        let exact_block_fence = load_exact_block_fence(&mmap[..HEADER_SIZE], file_len, &header)?;
        let verified_pages =
            verified_page_cache(header.section(SectionKind::BlockChecksums).offset)?;
        let reader = Self {
            file,
            mmap,
            header,
            verified_pages,
            filter_directory: None,
            exact_block_fence,
            filter_cache_identity: OnceLock::new(),
            resident_filter_body: OnceLock::new(),
        };
        let filter_directory = crate::jidx_filters::load(&reader)?;
        let reader = Self {
            filter_directory: Some(filter_directory),
            ..reader
        };
        reader.validate_documents()?;
        Ok(reader)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub(crate) fn exact_block_fence(&self) -> Option<&ExactBlockFence> {
        self.exact_block_fence.as_deref()
    }

    pub(crate) fn resident_filter_body(
        &self,
        directory: &crate::jidx_filters::FilterDirectory,
    ) -> Result<Option<&Arc<Vec<u8>>>, JidxError> {
        if self.exact_block_fence.is_none() {
            return Ok(None);
        }
        if self.resident_filter_body.get().is_none() {
            let body = crate::jidx_filters::acquire_resident_body(self, directory)?;
            let _ = self.resident_filter_body.set(body);
        }
        Ok(self.resident_filter_body.get().and_then(Option::as_ref))
    }

    pub(crate) fn filter_cache_identity(&self) -> Option<&FilterCacheIdentity> {
        self.exact_block_fence.as_ref()?;
        self.filter_cache_identity
            .get_or_init(|| {
                let file = self.cache_file_identity().ok().flatten()?;
                let header = self.header.encode().ok()?;
                let section = self.header.section(SectionKind::SeedFilters);
                Some(FilterCacheIdentity {
                    file,
                    header_sha256: sha256(&header),
                    body_sha256: self.header.body_sha256,
                    section_offset: section.offset,
                    section_length: section.length,
                })
            })
            .as_ref()
    }

    pub(crate) fn cache_file_identity(&self) -> std::io::Result<Option<[u64; 7]>> {
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            let metadata = self.file.metadata()?;
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
            let _ = &self.file;
            Ok(None)
        }
    }

    pub fn verify_checksum(&self) -> Result<(), JidxReaderError> {
        self.header.verify_body(&self.mmap)?;
        self.validate_contigs()?;
        crate::jidx_postings::validate_table(self)?;
        crate::jidx_filters::audit(
            self,
            self.filter_directory
                .as_ref()
                .ok_or(JidxError::Invalid("missing seed filter directory"))?,
        )?;
        for index in 0..self.header.seed_count {
            let seed = crate::jidx_postings::entry(self, index)?;
            self.seed_occurrences(seed)?;
        }
        Ok(())
    }

    pub(crate) fn verify_query_filter_pages<'a>(
        &self,
        keys: impl IntoIterator<Item = &'a u64>,
    ) -> Result<(), JidxReaderError> {
        #[cfg(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64")))]
        {
            let seeds = self.header.section(SectionKind::Seeds);
            if let (Ok(offset), Ok(length)) =
                (usize::try_from(seeds.offset), usize::try_from(seeds.length))
            {
                let _ = self.mmap.advise_range(Advice::Random, offset, length);
            }
        }
        Ok(crate::jidx_filters::verify_query_pages(
            self,
            self.filter_directory
                .as_ref()
                .ok_or(JidxError::Invalid("missing seed filter directory"))?,
            keys.into_iter().copied(),
        )?)
    }

    #[cfg(any(
        test,
        all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64"))
    ))]
    pub(crate) fn seed_record_page_span(&self, ordinal: u64) -> Option<(usize, usize)> {
        let seeds = self.header.section(SectionKind::Seeds);
        if ordinal >= self.header.seed_count {
            return None;
        }
        let record_start = ordinal
            .checked_mul(u64::from(crate::jidx::SEED_RECORD_SIZE))?
            .checked_add(seeds.offset)?;
        let record_end = record_start.checked_add(u64::from(crate::jidx::SEED_RECORD_SIZE))?;
        let seeds_end = seeds.offset.checked_add(seeds.length)?;
        if record_end > seeds_end {
            return None;
        }
        let page_start = record_start / PAGE_SIZE * PAGE_SIZE;
        let page_end = record_end
            .checked_add(PAGE_SIZE - 1)?
            .checked_div(PAGE_SIZE)?
            .checked_mul(PAGE_SIZE)?;
        let padding_end = seeds_end
            .checked_add(PAGE_SIZE - 1)?
            .checked_div(PAGE_SIZE)?
            .checked_mul(PAGE_SIZE)?;
        let next_section = self.header.section(SectionKind::DocumentPostings);
        let mmap_len = u64::try_from(self.mmap.len()).ok()?;
        if page_end > padding_end || padding_end > next_section.offset || page_end > mmap_len {
            return None;
        }
        Some((
            usize::try_from(page_start).ok()?,
            usize::try_from(page_end.checked_sub(page_start)?).ok()?,
        ))
    }

    pub(crate) fn advise_seed_record_pages(&self, midpoints: &[(u64, usize)]) {
        #[cfg(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64")))]
        {
            let mut previous = None;
            let mut run = None::<(usize, usize)>;
            for &(ordinal, _) in midpoints {
                if previous == Some(ordinal) {
                    continue;
                }
                previous = Some(ordinal);
                let Some((start, length)) = self.seed_record_page_span(ordinal) else {
                    continue;
                };
                let Some(end) = start.checked_add(length) else {
                    continue;
                };
                if let Some((run_start, run_end)) = run {
                    if start <= run_end {
                        run = Some((run_start, run_end.max(end)));
                        continue;
                    }
                    let _ =
                        self.mmap
                            .advise_range(Advice::WillNeed, run_start, run_end - run_start);
                }
                run = Some((start, end));
            }
            if let Some((start, end)) = run {
                let _ = self.mmap.advise_range(Advice::WillNeed, start, end - start);
            }
        }
        #[cfg(not(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64"))))]
        let _ = midpoints;
    }

    #[cfg(any(
        test,
        all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64"))
    ))]
    pub(crate) fn first_document_row_page_span(&self, seed: SeedEntry) -> Option<(usize, usize)> {
        if seed.document_frequency == 0 {
            return None;
        }
        let documents = self.header.section(SectionKind::DocumentPostings);
        let row_start = documents.offset.checked_add(seed.document_offset)?;
        let row_end = row_start.checked_add(u64::from(crate::jidx::DOCUMENT_POSTING_SIZE))?;
        let documents_end = documents.offset.checked_add(documents.length)?;
        if row_end > documents_end {
            return None;
        }
        let page_start = row_start / PAGE_SIZE * PAGE_SIZE;
        let page_end = row_end
            .checked_add(PAGE_SIZE - 1)?
            .checked_div(PAGE_SIZE)?
            .checked_mul(PAGE_SIZE)?;
        let padding_end = documents_end
            .checked_add(PAGE_SIZE - 1)?
            .checked_div(PAGE_SIZE)?
            .checked_mul(PAGE_SIZE)?;
        let next_section = self.header.section(SectionKind::ContigPostings);
        let mmap_len = u64::try_from(self.mmap.len()).ok()?;
        if page_end > padding_end || padding_end > next_section.offset || page_end > mmap_len {
            return None;
        }
        Some((
            usize::try_from(page_start).ok()?,
            usize::try_from(page_end.checked_sub(page_start)?).ok()?,
        ))
    }

    pub(crate) fn advise_first_document_rows(&self, seeds: &[Option<SeedEntry>]) {
        #[cfg(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64")))]
        {
            let mut run = None::<(usize, usize)>;
            for seed in seeds.iter().flatten() {
                let Some((start, length)) = self.first_document_row_page_span(*seed) else {
                    continue;
                };
                let Some(end) = start.checked_add(length) else {
                    continue;
                };
                if let Some((run_start, run_end)) = run {
                    if start >= run_start && start <= run_end {
                        run = Some((run_start, run_end.max(end)));
                        continue;
                    }
                    let _ =
                        self.mmap
                            .advise_range(Advice::WillNeed, run_start, run_end - run_start);
                }
                run = Some((start, end));
            }
            if let Some((start, end)) = run {
                let _ = self.mmap.advise_range(Advice::WillNeed, start, end - start);
            }
        }
        #[cfg(not(all(target_os = "linux", any(target_arch = "x86", target_arch = "x86_64"))))]
        let _ = seeds;
    }

    pub(crate) fn find_seeds_batch(
        &self,
        packed_keys: &[u64],
    ) -> Result<Vec<Option<SeedEntry>>, JidxReaderError> {
        if packed_keys.len() > SEED_LOOKUP_BATCH_KEYS {
            return Err(JidxError::Invalid("seed lookup batch").into());
        }
        let mut output = Vec::new();
        output
            .try_reserve_exact(packed_keys.len())
            .map_err(|_| JidxError::Invalid("seed lookup batch"))?;
        output.resize(packed_keys.len(), None);
        let mut positive_keys = Vec::new();
        positive_keys
            .try_reserve_exact(packed_keys.len())
            .map_err(|_| JidxError::Invalid("seed lookup batch"))?;
        let mut positive_indexes = Vec::new();
        positive_indexes
            .try_reserve_exact(packed_keys.len())
            .map_err(|_| JidxError::Invalid("seed lookup batch"))?;

        let directory = self
            .filter_directory
            .as_ref()
            .ok_or(JidxError::Invalid("missing seed filter directory"))?;
        for (index, &packed_key) in packed_keys.iter().enumerate() {
            seed_length(self.header.k, self.header.rescue_k15, packed_key)?;
            if crate::jidx_filters::contains(self, directory, packed_key)? {
                positive_keys.push(packed_key);
                positive_indexes.push(index);
            }
        }
        let positive_results = crate::jidx_postings::lookup_batch(self, &positive_keys)?;
        for (index, result) in positive_indexes.into_iter().zip(positive_results) {
            output[index] = result;
        }
        Ok(output)
    }

    pub fn find_seed(&self, packed_key: u64) -> Result<Option<SeedEntry>, JidxReaderError> {
        seed_length(self.header.k, self.header.rescue_k15, packed_key)?;
        if !crate::jidx_filters::contains(
            self,
            self.filter_directory
                .as_ref()
                .ok_or(JidxError::Invalid("missing seed filter directory"))?,
            packed_key,
        )? {
            return Ok(None);
        }
        Ok(crate::jidx_postings::lookup(self, packed_key)?)
    }

    pub fn seed_metagenomes(&self, seed: SeedEntry) -> Result<Vec<MetagenomeId>, JidxReaderError> {
        Ok(self
            .seed_documents(seed)?
            .into_iter()
            .map(|document| document.metagenome_id)
            .collect())
    }

    pub fn seed_documents(&self, seed: SeedEntry) -> Result<Vec<SeedDocument>, JidxReaderError> {
        Ok(crate::jidx_postings::documents(self, seed)?)
    }

    pub fn seed_document_occurrences(
        &self,
        seed: SeedEntry,
        document: SeedDocument,
    ) -> Result<Vec<SeedOccurrence>, JidxReaderError> {
        let occurrences = crate::jidx_postings::document_occurrences(self, seed, document)?;
        let k = seed_length(self.header.k, self.header.rescue_k15, seed.packed_key)?;
        for occurrence in &occurrences {
            let contig = self
                .contig(occurrence.contig_id)?
                .ok_or(JidxError::Invalid("missing occurrence contig"))?;
            if contig.metagenome_id != document.metagenome_id
                || occurrence
                    .position
                    .checked_add(u64::from(k))
                    .is_none_or(|end| end > contig.length)
            {
                return Err(JidxError::Invalid("contig posting position").into());
            }
        }
        Ok(occurrences)
    }

    pub fn seed_occurrences(
        &self,
        seed: SeedEntry,
    ) -> Result<Vec<SeedOccurrence>, JidxReaderError> {
        let documents = self.seed_documents(seed)?;
        let total = documents.iter().try_fold(0usize, |total, document| {
            usize::try_from(document.occurrence_count)
                .ok()
                .and_then(|count| total.checked_add(count))
                .ok_or(JidxError::Invalid("occurrence count"))
        })?;
        let mut occurrences = Vec::new();
        occurrences
            .try_reserve_exact(total)
            .map_err(|_| JidxError::Invalid("occurrence count"))?;
        for document in documents {
            occurrences.extend(self.seed_document_occurrences(seed, document)?);
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

    pub(crate) fn document_record(&self, id: MetagenomeId) -> Result<DocumentRecord, JidxError> {
        if id >= self.header.document_count {
            return Err(JidxError::Invalid("metagenome ID"));
        }
        DocumentRecord::decode(self.record_bytes(
            SectionKind::Documents,
            u64::from(id),
            DOCUMENT_RECORD_SIZE,
        )?)
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
            let page_bit = page
                .checked_sub(1)
                .ok_or(JidxError::Invalid("page checksum"))?;
            let word_index =
                usize::try_from(page_bit / 64).map_err(|_| JidxError::Invalid("page checksum"))?;
            let bit_mask = 1u64 << (page_bit % 64);
            let word = self
                .verified_pages
                .get(word_index)
                .ok_or(JidxError::Invalid("page checksum"))?;
            if word.load(Ordering::Relaxed) & bit_mask != 0 {
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
            // The mmap is immutable; this bit memoizes only the completed hash comparison.
            word.fetch_or(bit_mask, Ordering::Relaxed);
        }
        Ok(bytes)
    }

    #[cfg(unix)]
    pub(crate) fn read_verified_filter_body(
        &self,
        identity: &FilterCacheIdentity,
    ) -> Result<Option<Vec<u8>>, JidxError> {
        use std::os::unix::fs::FileExt;

        const DATA_CHUNK: usize = 1024 * 1024;
        const CHECKSUM_CHUNK: usize = DATA_CHUNK / PAGE_SIZE as usize * 32;
        match self.cache_file_identity() {
            Ok(Some(current)) if current == identity.file => {}
            Ok(Some(_)) => return Err(JidxError::Invalid("resident seed filter identity")),
            _ => return Ok(None),
        }
        let section = self.header.section(SectionKind::SeedFilters);
        if identity.section_offset != section.offset
            || identity.section_length != section.length
            || identity.body_sha256 != self.header.body_sha256
            || identity.header_sha256 != sha256(&self.header.encode()?)
        {
            return Err(JidxError::Invalid("resident seed filter identity"));
        }
        let length = usize::try_from(section.length)
            .map_err(|_| JidxError::Invalid("resident seed filter length"))?;
        let mut body = Vec::new();
        if body.try_reserve_exact(length).is_err() {
            return Ok(None);
        }
        body.resize(length, 0);
        let page_count = section
            .length
            .checked_add(PAGE_SIZE - 1)
            .ok_or(JidxError::Invalid("resident seed filter pages"))?
            / PAGE_SIZE;
        let checksums = self.header.section(SectionKind::BlockChecksums);
        let first_page = section.offset / PAGE_SIZE;
        let mut data = Vec::new();
        let mut expected = Vec::new();
        if data.try_reserve_exact(DATA_CHUNK).is_err()
            || expected.try_reserve_exact(CHECKSUM_CHUNK).is_err()
        {
            return Ok(None);
        }
        data.resize(DATA_CHUNK, 0);
        expected.resize(CHECKSUM_CHUNK, 0);

        let mut page_offset = 0u64;
        let mut copied = 0usize;
        while page_offset < page_count {
            let pages = (page_count - page_offset).min((DATA_CHUNK as u64) / PAGE_SIZE);
            let data_length = usize::try_from(pages * PAGE_SIZE).unwrap();
            let checksum_length = usize::try_from(pages * 32).unwrap();
            let data_offset = section
                .offset
                .checked_add(
                    page_offset
                        .checked_mul(PAGE_SIZE)
                        .ok_or(JidxError::Invalid("resident seed filter range"))?,
                )
                .ok_or(JidxError::Invalid("resident seed filter range"))?;
            let checksum_offset = checksums
                .offset
                .checked_add(
                    first_page
                        .checked_add(page_offset)
                        .and_then(|page| page.checked_sub(1))
                        .and_then(|page| page.checked_mul(32))
                        .ok_or(JidxError::Invalid("resident seed filter checksum"))?,
                )
                .ok_or(JidxError::Invalid("resident seed filter checksum"))?;
            if self
                .file
                .read_exact_at(&mut data[..data_length], data_offset)
                .is_err()
                || self
                    .file
                    .read_exact_at(&mut expected[..checksum_length], checksum_offset)
                    .is_err()
            {
                return Ok(None);
            }
            for page in 0..usize::try_from(pages).unwrap() {
                let data_start = page * PAGE_SIZE as usize;
                let checksum_start = page * 32;
                if sha256(&data[data_start..data_start + PAGE_SIZE as usize])
                    != expected[checksum_start..checksum_start + 32]
                {
                    return Err(JidxError::ChecksumMismatch);
                }
            }
            let copy_length = (length - copied).min(data_length);
            body[copied..copied + copy_length].copy_from_slice(&data[..copy_length]);
            copied += copy_length;
            page_offset += pages;
        }
        if copied != length {
            return Err(JidxError::Invalid("resident seed filter length"));
        }
        match self.cache_file_identity() {
            Ok(Some(current)) if current == identity.file => {}
            Ok(Some(_)) => return Err(JidxError::Invalid("resident seed filter identity")),
            _ => return Ok(None),
        }
        if sha256(&self.header.encode()?) != identity.header_sha256 {
            return Err(JidxError::Invalid("resident seed filter identity"));
        }
        Ok(Some(body))
    }

    #[cfg(not(unix))]
    pub(crate) fn read_verified_filter_body(
        &self,
        _identity: &FilterCacheIdentity,
    ) -> Result<Option<Vec<u8>>, JidxError> {
        Ok(None)
    }

    fn record_bytes(&self, kind: SectionKind, index: u64, size: u32) -> Result<&[u8], JidxError> {
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
        self.checked_bytes(start, end)
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

fn load_exact_block_fence(
    source_header: &[u8],
    source_file_len: u64,
    header: &Header,
) -> Result<Option<Arc<ExactBlockFence>>, JidxReaderError> {
    let path = std::env::var_os("JAM_EXACT_BLOCK_FENCE");
    let expected_sha = std::env::var_os("JAM_EXACT_BLOCK_FENCE_SHA256");
    let (path, expected_sha) = match (path, expected_sha) {
        (None, None) => return Ok(None),
        (Some(path), Some(expected_sha)) => (path, expected_sha),
        _ => return Err(JidxError::Invalid("exact block fence environment").into()),
    };
    let path = PathBuf::from(path);
    if !path.is_absolute() {
        return Err(JidxError::Invalid("exact block fence path").into());
    }
    let expected_sha = parse_sha256(
        expected_sha
            .to_str()
            .ok_or(JidxError::Invalid("exact block fence digest"))?,
    )?;
    let cache = EXACT_BLOCK_FENCE_CACHE.get_or_init(|| Mutex::new(None));
    let mut cache = cache
        .lock()
        .map_err(|_| JidxError::Invalid("exact block fence cache"))?;
    if let Some(cached) = cache.as_ref()
        && cached.path == path
        && cached.sha256 == expected_sha
    {
        validate_exact_block_source(&cached.fence, source_header, source_file_len, header)?;
        return Ok(Some(Arc::clone(&cached.fence)));
    }

    let expected_len = header
        .seed_count
        .checked_add(EXACT_BLOCK_RECORDS - 1)
        .and_then(|count| (count / EXACT_BLOCK_RECORDS).checked_mul(8))
        .and_then(|body| body.checked_add(EXACT_BLOCK_HEADER_SIZE as u64))
        .ok_or(JidxError::Invalid("exact block fence length"))?;
    let metadata = std::fs::symlink_metadata(&path)?;
    if expected_len > EXACT_BLOCK_MAX_BYTES || !metadata.is_file() || metadata.len() != expected_len
    {
        return Err(JidxError::Invalid("exact block fence length").into());
    }
    let file = File::open(&path)?;
    let opened = file.metadata()?;
    if !opened.is_file() || opened.len() != expected_len {
        return Err(JidxError::Invalid("exact block fence length").into());
    }
    let mut bytes = Vec::with_capacity(expected_len as usize);
    file.take(expected_len + 1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 != expected_len {
        return Err(JidxError::Invalid("exact block fence length").into());
    }
    if sha256(&bytes) != expected_sha {
        return Err(JidxError::ChecksumMismatch.into());
    }
    let fence = Arc::new(decode_exact_block_fence(
        &bytes,
        source_header,
        source_file_len,
        header,
    )?);
    *cache = Some(ExactBlockFenceCache {
        path,
        sha256: expected_sha,
        fence: Arc::clone(&fence),
    });
    Ok(Some(fence))
}

fn parse_sha256(value: &str) -> Result<[u8; 32], JidxError> {
    if value.len() != 64 || !value.is_ascii() {
        return Err(JidxError::Invalid("exact block fence digest"));
    }
    let mut digest = [0; 32];
    for (index, byte) in digest.iter_mut().enumerate() {
        *byte = u8::from_str_radix(&value[index * 2..index * 2 + 2], 16)
            .map_err(|_| JidxError::Invalid("exact block fence digest"))?;
    }
    Ok(digest)
}

fn validate_exact_block_source(
    fence: &ExactBlockFence,
    source_header: &[u8],
    source_file_len: u64,
    header: &Header,
) -> Result<(), JidxError> {
    let expected_count = header
        .seed_count
        .checked_add(EXACT_BLOCK_RECORDS - 1)
        .ok_or(JidxError::Invalid("exact block fence count"))?
        / EXACT_BLOCK_RECORDS;
    if u64::try_from(fence.first_keys.len()).ok() != Some(expected_count)
        || source_header.len() != HEADER_SIZE
        || fence.source_file_len != source_file_len
        || fence.seeds_offset != header.section(SectionKind::Seeds).offset
        || fence.seeds_length != header.section(SectionKind::Seeds).length
        || fence.seed_count != header.seed_count
        || fence.source_header_sha256 != sha256(source_header)
        || fence.source_body_sha256 != header.body_sha256
    {
        return Err(JidxError::Invalid("exact block fence source"));
    }
    Ok(())
}

fn decode_exact_block_fence(
    bytes: &[u8],
    source_header: &[u8],
    source_file_len: u64,
    header: &Header,
) -> Result<ExactBlockFence, JidxError> {
    if bytes.len() < EXACT_BLOCK_HEADER_SIZE
        || bytes[..8] != *b"JXFENCE1"
        || u16::from_le_bytes(bytes[8..10].try_into().unwrap()) != 1
        || u16::from_le_bytes(bytes[10..12].try_into().unwrap()) != EXACT_BLOCK_HEADER_SIZE as u16
        || read_u32(bytes, 12) != SEED_RECORD_SIZE
        || read_u32(bytes, 16) != EXACT_BLOCK_RECORDS as u32
        || bytes[20..24].iter().any(|byte| *byte != 0)
        || bytes[168..4064].iter().any(|byte| *byte != 0)
        || sha256(&bytes[..4064]) != bytes[4064..4096]
    {
        return Err(JidxError::Invalid("exact block fence header"));
    }
    let seeds = header.section(SectionKind::Seeds);
    let fence_count = header
        .seed_count
        .checked_add(EXACT_BLOCK_RECORDS - 1)
        .ok_or(JidxError::Invalid("exact block fence count"))?
        / EXACT_BLOCK_RECORDS;
    let body_len = fence_count
        .checked_mul(8)
        .ok_or(JidxError::Invalid("exact block fence length"))?;
    let file_len =
        u64::try_from(bytes.len()).map_err(|_| JidxError::Invalid("exact block fence length"))?;
    if read_u64(bytes, 24) != source_file_len
        || read_u64(bytes, 32) != seeds.offset
        || read_u64(bytes, 40) != seeds.length
        || read_u64(bytes, 48) != header.seed_count
        || read_u64(bytes, 56) != fence_count
        || bytes[64..96] != sha256(source_header)
        || bytes[96..128] != header.body_sha256
        || read_u64(bytes, 128) != body_len
        || file_len != EXACT_BLOCK_HEADER_SIZE as u64 + body_len
        || sha256(&bytes[EXACT_BLOCK_HEADER_SIZE..]) != bytes[136..168]
    {
        return Err(JidxError::Invalid("exact block fence binding"));
    }
    let mut first_keys = Vec::new();
    first_keys
        .try_reserve_exact(
            usize::try_from(fence_count)
                .map_err(|_| JidxError::Invalid("exact block fence count"))?,
        )
        .map_err(|_| JidxError::Invalid("exact block fence count"))?;
    for chunk in bytes[EXACT_BLOCK_HEADER_SIZE..].as_chunks::<8>().0 {
        let key = u64::from_le_bytes(*chunk);
        if first_keys.last().is_some_and(|previous| *previous >= key) {
            return Err(JidxError::Invalid("exact block fence order"));
        }
        first_keys.push(key);
    }
    let fence = ExactBlockFence {
        first_keys: first_keys.into_boxed_slice(),
        source_file_len,
        seeds_offset: seeds.offset,
        seeds_length: seeds.length,
        seed_count: header.seed_count,
        source_header_sha256: sha256(source_header),
        source_body_sha256: header.body_sha256,
    };
    validate_exact_block_source(&fence, source_header, source_file_len, header)?;
    Ok(fence)
}

fn verified_page_cache(checksum_offset: u64) -> Result<Box<[AtomicU64]>, JidxError> {
    let data_pages = checksum_offset
        .checked_div(PAGE_SIZE)
        .and_then(|pages| pages.checked_sub(1))
        .ok_or(JidxError::Invalid("page checksum cache"))?;
    let words = data_pages
        .checked_add(63)
        .ok_or(JidxError::Invalid("page checksum cache"))?
        / 64;
    let words = usize::try_from(words).map_err(|_| JidxError::Invalid("page checksum cache"))?;
    let mut cache = Vec::new();
    cache
        .try_reserve_exact(words)
        .map_err(|_| JidxError::Invalid("page checksum cache"))?;
    cache.resize_with(words, || AtomicU64::new(0));
    Ok(cache.into_boxed_slice())
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
    use crate::jidx::{SECTION_COUNT, SEED_RECORD_SIZE, SectionDescriptor, VERSION, sha256};
    use std::io::{Read, Seek, SeekFrom, Write};
    use std::path::PathBuf;
    use xorf::{BinaryFuse8, DmaSerializable, Filter};

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

    fn seed_filters(keys: &[u64]) -> Vec<u8> {
        if keys.is_empty() {
            return Vec::new();
        }
        let filter = BinaryFuse8::try_from(keys).unwrap();
        let fingerprints = filter.dma_fingerprints();
        let filter_length = 20 + fingerprints.len();
        let mut bytes = vec![0; PAGE_SIZE as usize + filter_length];
        put_u64(&mut bytes, 0, keys[0]);
        put_u64(&mut bytes, 8, *keys.last().unwrap());
        put_u64(&mut bytes, 16, PAGE_SIZE);
        put_u64(&mut bytes, 24, filter_length as u64);
        put_u32(&mut bytes, 32, keys.len() as u32);
        put_u16(&mut bytes, 36, 20);
        filter.dma_copy_descriptor_to(&mut bytes[PAGE_SIZE as usize..PAGE_SIZE as usize + 20]);
        bytes[PAGE_SIZE as usize + 20..].copy_from_slice(fingerprints);
        bytes
    }

    fn section_range(path: &Path, kind: SectionKind) -> (u64, u64) {
        let mut header = [0; HEADER_SIZE];
        File::open(path).unwrap().read_exact(&mut header).unwrap();
        let descriptor = 144 + (kind as usize - 1) * 24;
        (
            u64::from_le_bytes(header[descriptor + 8..descriptor + 16].try_into().unwrap()),
            u64::from_le_bytes(header[descriptor + 16..descriptor + 24].try_into().unwrap()),
        )
    }

    fn section_offset(path: &Path, kind: SectionKind) -> u64 {
        section_range(path, kind).0
    }

    fn corrupt_byte(path: &Path, offset: u64) {
        let mut file = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(path)
            .unwrap();
        file.seek(SeekFrom::Start(offset)).unwrap();
        let mut byte = [0];
        file.read_exact(&mut byte).unwrap();
        byte[0] ^= 0xff;
        file.seek(SeekFrom::Start(offset)).unwrap();
        file.write_all(&byte).unwrap();
    }

    fn replace_bytes(path: &Path, offset: u64, replacement: &[u8]) {
        let mut bytes = std::fs::read(path).unwrap();
        let offset = offset as usize;
        bytes[offset..offset + replacement.len()].copy_from_slice(replacement);
        std::fs::write(path, bytes).unwrap();
    }

    fn rewrite_checksums(path: &Path) {
        let mut bytes = std::fs::read(path).unwrap();
        let (checksum_offset, checksum_length) = section_range(path, SectionKind::BlockChecksums);
        let checksum_count = checksum_length / 32;
        for page in 1..=checksum_count {
            let page_start = (page * PAGE_SIZE) as usize;
            let digest = sha256(&bytes[page_start..page_start + PAGE_SIZE as usize]);
            let checksum_start = (checksum_offset + (page - 1) * 32) as usize;
            bytes[checksum_start..checksum_start + 32].copy_from_slice(&digest);
        }
        let body_sha256 = sha256(&bytes[HEADER_SIZE..]);
        bytes[112..144].copy_from_slice(&body_sha256);
        std::fs::write(path, bytes).unwrap();
    }

    fn exact_fence_bytes(path: &Path) -> (Vec<u8>, Header, Vec<u8>) {
        let source = std::fs::read(path).unwrap();
        let header = Header::decode_header(&source[..HEADER_SIZE], source.len() as u64).unwrap();
        let seeds = header.section(SectionKind::Seeds);
        let mut body = Vec::new();
        for ordinal in (0..header.seed_count).step_by(EXACT_BLOCK_RECORDS as usize) {
            let offset = (seeds.offset + ordinal * u64::from(SEED_RECORD_SIZE)) as usize;
            body.extend_from_slice(&source[offset..offset + 8]);
        }
        let mut bytes = vec![0; EXACT_BLOCK_HEADER_SIZE + body.len()];
        bytes[..8].copy_from_slice(b"JXFENCE1");
        put_u16(&mut bytes, 8, 1);
        put_u16(&mut bytes, 10, EXACT_BLOCK_HEADER_SIZE as u16);
        put_u32(&mut bytes, 12, SEED_RECORD_SIZE);
        put_u32(&mut bytes, 16, EXACT_BLOCK_RECORDS as u32);
        put_u64(&mut bytes, 24, source.len() as u64);
        put_u64(&mut bytes, 32, seeds.offset);
        put_u64(&mut bytes, 40, seeds.length);
        put_u64(&mut bytes, 48, header.seed_count);
        put_u64(&mut bytes, 56, body.len() as u64 / 8);
        bytes[64..96].copy_from_slice(&sha256(&source[..HEADER_SIZE]));
        bytes[96..128].copy_from_slice(&header.body_sha256);
        put_u64(&mut bytes, 128, body.len() as u64);
        bytes[136..168].copy_from_slice(&sha256(&body));
        bytes[EXACT_BLOCK_HEADER_SIZE..].copy_from_slice(&body);
        let header_sha = sha256(&bytes[..4064]);
        bytes[4064..4096].copy_from_slice(&header_sha);
        (source, header, bytes)
    }

    fn exact_fence(path: &Path) -> ExactBlockFence {
        let (source, header, bytes) = exact_fence_bytes(path);
        decode_exact_block_fence(&bytes, &source[..HEADER_SIZE], source.len() as u64, &header)
            .unwrap()
    }

    fn fixture(
        bad_metagenome_id: bool,
        corrupt_padding: bool,
        bad_document_posting: bool,
        bad_occurrence: bool,
    ) -> (tempfile::TempDir, PathBuf) {
        fixture_with_string_pages(
            bad_metagenome_id,
            corrupt_padding,
            bad_document_posting,
            bad_occurrence,
            true,
            None,
        )
    }

    fn fixture_with_string_pages(
        bad_metagenome_id: bool,
        corrupt_padding: bool,
        bad_document_posting: bool,
        bad_occurrence: bool,
        include_seed: bool,
        string_pages: Option<u64>,
    ) -> (tempfile::TempDir, PathBuf) {
        let mut strings = Vec::new();
        let name = push_string(&mut strings, "doc");
        let bgzf = push_string(&mut strings, "seq.bgz");
        let contig_name = push_string(&mut strings, "contig");
        if let Some(pages) = string_pages {
            strings.resize((pages * PAGE_SIZE) as usize, 0);
        }

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

        let mut seed = vec![
            0;
            if include_seed {
                SEED_RECORD_SIZE as usize
            } else {
                0
            }
        ];
        let mut document_posting = vec![0; if include_seed { 16 } else { 0 }];
        if include_seed {
            put_u64(&mut seed, 0, 0x1234);
            put_u64(&mut seed, 8, 0);
            put_u32(&mut seed, 16, 1);
            put_u32(&mut document_posting, 0, u32::from(bad_document_posting));
            put_u32(&mut document_posting, 4, 0);
            put_u64(
                &mut document_posting,
                8,
                ((if bad_occurrence { 90 } else { 2 }) << 1) | 1,
            );
        }
        let contig_posting = Vec::new();
        let payloads = [
            strings,
            document,
            contig,
            seed,
            document_posting,
            contig_posting,
            0u64.to_le_bytes().to_vec(),
            seed_filters(if include_seed { &[0x1234] } else { &[] }),
        ];
        write_fixture(
            payloads,
            1,
            1,
            if include_seed { 1 } else { 0 },
            if include_seed { 1 } else { 0 },
            corrupt_padding,
        )
    }

    fn write_fixture(
        payloads: [Vec<u8>; 8],
        document_count: u32,
        contig_count: u32,
        seed_count: u64,
        occurrence_count: u64,
        corrupt_padding: bool,
    ) -> (tempfile::TempDir, PathBuf) {
        let mut offset = PAGE_SIZE;
        let mut sections = Vec::with_capacity(SECTION_COUNT);
        for (index, payload) in payloads.iter().enumerate() {
            let descriptor = SectionDescriptor {
                kind: SectionKind::ALL[index],
                record_size: SectionKind::ALL[index].record_size(),
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
        header[18] = 2;
        header[19] = 1;
        put_u16(&mut header, 20, SECTION_COUNT as u16);
        put_u32(&mut header, 24, document_count);
        put_u32(&mut header, 28, contig_count);
        put_u64(&mut header, 32, seed_count);
        put_u64(&mut header, 40, occurrence_count);
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
        put_u16(&mut header, 360, 16);
        file[..HEADER_SIZE].copy_from_slice(&header);
        if corrupt_padding {
            file[sections[2].offset as usize + sections[2].length as usize] ^= 1;
        }
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("fixture.jidx");
        std::fs::write(&path, file).unwrap();
        (directory, path)
    }

    fn external_fixture() -> (tempfile::TempDir, PathBuf) {
        external_fixture_with_keys([0x1234, 0x1235])
    }

    fn external_fixture_with_keys(keys: [u64; 2]) -> (tempfile::TempDir, PathBuf) {
        let mut strings = Vec::new();
        let doc0 = push_string(&mut strings, "doc0");
        let uri0 = push_string(&mut strings, "doc0.bgz");
        let doc1 = push_string(&mut strings, "doc1");
        let uri1 = push_string(&mut strings, "doc1.bgz");
        let names = [
            push_string(&mut strings, "contig0"),
            push_string(&mut strings, "contig1"),
            push_string(&mut strings, "contig2"),
        ];

        let mut documents = vec![0; 2 * DOCUMENT_RECORD_SIZE as usize];
        for (index, (name, uri, start, count, gzi_offset)) in
            [(doc0, uri0, 0, 2, 0), (doc1, uri1, 2, 1, 8)]
                .into_iter()
                .enumerate()
        {
            let offset = index * DOCUMENT_RECORD_SIZE as usize;
            put_string(&mut documents, offset, name);
            put_string(&mut documents, offset + 8, uri);
            put_u64(&mut documents, offset + 16, 100);
            put_u32(&mut documents, offset + 24, start);
            put_u32(&mut documents, offset + 28, count);
            documents[offset + 32..offset + 64].fill(3 + index as u8);
            put_u64(&mut documents, offset + 64, gzi_offset);
            put_u64(&mut documents, offset + 72, 8);
        }

        let mut contigs = vec![0; 3 * CONTIG_RECORD_SIZE as usize];
        for (index, (document, length)) in
            [(0, u64::MAX), (0, 100), (1, 100)].into_iter().enumerate()
        {
            let offset = index * CONTIG_RECORD_SIZE as usize;
            put_u32(&mut contigs, offset, document);
            put_string(&mut contigs, offset + 4, names[index]);
            put_u64(&mut contigs, offset + 16, length);
            put_u64(&mut contigs, offset + 24, 5);
            put_u32(&mut contigs, offset + 32, 100);
            put_u32(&mut contigs, offset + 36, 101);
        }

        let mut seeds = vec![0; 2 * SEED_RECORD_SIZE as usize];
        put_u64(&mut seeds, 0, keys[0]);
        put_u64(&mut seeds, 8, 0);
        put_u32(&mut seeds, 16, 2);
        put_u64(&mut seeds, SEED_RECORD_SIZE as usize, keys[1]);
        put_u64(&mut seeds, SEED_RECORD_SIZE as usize + 8, 48);
        put_u32(&mut seeds, SEED_RECORD_SIZE as usize + 16, 1);

        let mut document_postings = vec![0; 80];
        put_u32(&mut document_postings, 0, 0);
        put_u32(&mut document_postings, 4, 1);
        put_u64(&mut document_postings, 8, 5);
        put_u32(&mut document_postings, 16, 1);
        put_u32(&mut document_postings, 20, u32::MAX);
        put_u64(&mut document_postings, 24, 0);
        put_u64(&mut document_postings, 32, 2);
        put_u64(&mut document_postings, 40, 4);
        put_u32(&mut document_postings, 48, 0);
        put_u32(&mut document_postings, 52, u32::MAX);
        put_u64(&mut document_postings, 56, 4);
        put_u64(&mut document_postings, 64, 1);
        put_u64(&mut document_postings, 72, 11);

        let cold = vec![
            0x00, 0x05, 0x01, 0x04, 0x00, 0x85, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
            0x01,
        ];
        write_fixture(
            [
                strings,
                documents,
                contigs,
                seeds,
                document_postings,
                cold,
                vec![0; 16],
                seed_filters(&keys),
            ],
            2,
            3,
            2,
            4,
            false,
        )
    }

    fn multi_page_filter_fixture() -> (tempfile::TempDir, PathBuf) {
        let keys = (0..5_000u64).collect::<Vec<_>>();
        let mut strings = Vec::new();
        let name = push_string(&mut strings, "doc");
        let bgzf = push_string(&mut strings, "seq.bgz");
        let contig_name = push_string(&mut strings, "contig");
        let mut document = vec![0; DOCUMENT_RECORD_SIZE as usize];
        put_string(&mut document, 0, name);
        put_string(&mut document, 8, bgzf);
        put_u64(&mut document, 16, 100);
        put_u32(&mut document, 28, 1);
        document[32..64].fill(3);
        put_u64(&mut document, 72, 8);
        let mut contig = vec![0; CONTIG_RECORD_SIZE as usize];
        put_string(&mut contig, 4, contig_name);
        put_u64(&mut contig, 16, 100);
        put_u64(&mut contig, 24, 5);
        put_u32(&mut contig, 32, 100);
        put_u32(&mut contig, 36, 101);
        let mut seeds = vec![0; keys.len() * SEED_RECORD_SIZE as usize];
        let mut document_postings = vec![0; keys.len() * 16];
        for (index, key) in keys.iter().copied().enumerate() {
            let seed = index * SEED_RECORD_SIZE as usize;
            put_u64(&mut seeds, seed, key);
            put_u64(&mut seeds, seed + 8, (index * 16) as u64);
            put_u32(&mut seeds, seed + 16, 1);
            put_u64(&mut document_postings, index * 16 + 8, 5);
        }
        write_fixture(
            [
                strings,
                document,
                contig,
                seeds,
                document_postings,
                Vec::new(),
                0u64.to_le_bytes().to_vec(),
                seed_filters(&keys),
            ],
            1,
            1,
            keys.len() as u64,
            keys.len() as u64,
            false,
        )
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
        let (_directory, path) = external_fixture();
        corrupt_byte(&path, section_offset(&path, SectionKind::ContigPostings));

        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        let documents = reader.seed_documents(seed).unwrap();
        assert!(matches!(
            reader.seed_document_occurrences(seed, documents[1]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
        assert!(matches!(
            reader.seed_document_occurrences(seed, documents[1]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
    }

    #[test]
    fn verified_page_bitmap_retains_more_than_4096_pages() {
        let pages = 4097u64;
        let (_directory, path) =
            fixture_with_string_pages(false, false, false, false, true, Some(pages));
        let reader = JidxReader::open(path).unwrap();
        let strings = reader.header.section(SectionKind::Strings);
        reader
            .checked_bytes(strings.offset, strings.offset + strings.length)
            .unwrap();

        for page in [1, pages] {
            let page_bit = page - 1;
            let word = &reader.verified_pages[(page_bit / 64) as usize];
            assert_ne!(word.load(Ordering::Relaxed) & (1u64 << (page_bit % 64)), 0);
        }
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
    fn filter_false_positive_still_requires_exact_lookup() {
        let (_directory, path) = external_fixture_with_keys([0, 1_000_000]);
        let reader = JidxReader::open(path).unwrap();
        let directory = reader.filter_directory.as_ref().unwrap();
        let false_positive = (1..1_000_000u64)
            .find(|key| crate::jidx_filters::contains(&reader, directory, *key).unwrap())
            .unwrap();
        assert!(reader.find_seed(false_positive).unwrap().is_none());
    }

    #[test]
    fn page_local_filter_matches_binary_fuse_reference() {
        let (_directory, path) = external_fixture_with_keys([0, 1_000_000]);
        let reader = JidxReader::open(path).unwrap();
        let directory = reader.filter_directory.as_ref().unwrap();
        for key in [0, 1_000_000].into_iter().chain(1..100_000u64) {
            assert_eq!(
                crate::jidx_filters::contains(&reader, directory, key).unwrap(),
                crate::jidx_filters::reference_contains(&reader, directory, key).unwrap(),
                "key {key}"
            );
        }
    }

    #[test]
    fn query_filter_page_verification_touches_only_exact_pages() {
        let key = 1_234u64;
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(&path).unwrap();
        let offsets = crate::jidx_filters::fingerprint_offsets(
            &reader,
            reader.filter_directory.as_ref().unwrap(),
            key,
        )
        .unwrap();
        let mut expected = reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        for page in offsets.map(|offset| offset / PAGE_SIZE) {
            let page_bit = page - 1;
            expected[(page_bit / 64) as usize] |= 1u64 << (page_bit % 64);
        }

        reader.verify_query_filter_pages([&key]).unwrap();

        assert_eq!(
            reader
                .verified_pages
                .iter()
                .map(|word| word.load(Ordering::Relaxed))
                .collect::<Vec<_>>(),
            expected
        );

        drop(reader);
        corrupt_byte(&path, offsets[0]);
        let reader = JidxReader::open(path).unwrap();
        assert!(matches!(
            reader.verify_query_filter_pages([&key]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
    }

    #[test]
    fn seed_record_page_spans_are_bounded_and_advice_stays_lazy() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(path).unwrap();
        let seeds = reader.header.section(SectionKind::Seeds);
        let page = usize::try_from(PAGE_SIZE).unwrap();
        let offset = usize::try_from(seeds.offset).unwrap();
        assert_eq!(reader.seed_record_page_span(0), Some((offset, page)));
        assert_eq!(reader.seed_record_page_span(169), Some((offset, page)));
        assert_eq!(reader.seed_record_page_span(170), Some((offset, 2 * page)));

        let final_span = reader
            .seed_record_page_span(reader.header.seed_count - 1)
            .unwrap();
        let final_end = final_span.0 + final_span.1;
        assert!(final_end as u64 >= seeds.offset + seeds.length);
        assert!(final_end as u64 <= reader.header.section(SectionKind::DocumentPostings).offset);
        assert!(final_end <= reader.mmap.len());
        assert_eq!(reader.seed_record_page_span(reader.header.seed_count), None);

        let before = reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        reader.advise_seed_record_pages(&[
            (0, 0),
            (0, 1),
            (169, 2),
            (170, 3),
            (341, 4),
            (342, 5),
            (reader.header.seed_count - 1, 6),
        ]);
        assert_eq!(
            reader
                .verified_pages
                .iter()
                .map(|word| word.load(Ordering::Relaxed))
                .collect::<Vec<_>>(),
            before
        );
    }

    #[test]
    fn first_document_row_spans_are_bounded_and_advice_stays_lazy() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(path).unwrap();
        let documents = reader.header.section(SectionKind::DocumentPostings);
        let page = usize::try_from(PAGE_SIZE).unwrap();
        let offset = usize::try_from(documents.offset).unwrap();
        let seed = SeedEntry {
            packed_key: 0,
            document_frequency: 1,
            document_offset: 0,
        };
        assert_eq!(
            reader.first_document_row_page_span(seed),
            Some((offset, page))
        );
        let boundary = SeedEntry {
            document_offset: PAGE_SIZE - u64::from(crate::jidx::DOCUMENT_POSTING_SIZE),
            ..seed
        };
        assert_eq!(
            reader.first_document_row_page_span(boundary),
            Some((offset, page))
        );
        let adjacent = SeedEntry {
            document_offset: PAGE_SIZE,
            ..seed
        };
        assert_eq!(
            reader.first_document_row_page_span(adjacent),
            Some((offset + page, page))
        );
        let final_seed = SeedEntry {
            document_offset: documents.length - u64::from(crate::jidx::DOCUMENT_POSTING_SIZE),
            ..seed
        };
        let final_span = reader.first_document_row_page_span(final_seed).unwrap();
        let final_end = final_span.0 + final_span.1;
        assert!(final_end as u64 >= documents.offset + documents.length);
        assert!(final_end as u64 <= reader.header.section(SectionKind::ContigPostings).offset);
        assert!(final_end <= reader.mmap.len());
        assert!(
            reader
                .first_document_row_page_span(SeedEntry {
                    document_frequency: 0,
                    ..seed
                })
                .is_none()
        );
        assert!(
            reader
                .first_document_row_page_span(SeedEntry {
                    document_offset: documents.length,
                    ..seed
                })
                .is_none()
        );

        let before = reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        let same_page = SeedEntry {
            document_offset: u64::from(crate::jidx::DOCUMENT_POSTING_SIZE),
            ..seed
        };
        let gap = SeedEntry {
            document_offset: 3 * PAGE_SIZE,
            ..seed
        };
        let decrease = SeedEntry {
            document_offset: 2 * PAGE_SIZE,
            ..seed
        };
        reader.advise_first_document_rows(&[
            None,
            Some(seed),
            Some(seed),
            Some(same_page),
            Some(adjacent),
            Some(gap),
            Some(decrease),
            Some(final_seed),
        ]);
        assert_eq!(
            reader
                .verified_pages
                .iter()
                .map(|word| word.load(Ordering::Relaxed))
                .collect::<Vec<_>>(),
            before
        );
    }

    #[test]
    fn first_document_row_advice_preserves_public_errors() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(&path).unwrap();
        let seed = reader.find_seed(170).unwrap().unwrap();
        let documents = reader.header.section(SectionKind::DocumentPostings);
        drop(reader);
        corrupt_byte(&path, documents.offset + seed.document_offset);

        let reader = JidxReader::open(&path).unwrap();
        let seed = reader.find_seed(170).unwrap().unwrap();
        reader.advise_first_document_rows(&[Some(seed)]);
        assert!(matches!(
            reader.seed_documents(seed),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));

        let invalid = SeedEntry {
            document_offset: documents.length,
            ..seed
        };
        reader.advise_first_document_rows(&[Some(invalid)]);
        assert!(matches!(
            reader.seed_documents(invalid),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "document posting range"
            )))
        ));
    }

    #[test]
    fn batch_seed_lookup_matches_scalar_across_pages_and_gaps() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(path).unwrap();
        assert_eq!(reader.find_seeds_batch(&[]).unwrap(), Vec::new());

        let mut state = 41u64;
        let mut keys = vec![0, 1, 170, 4_998, 4_999, 5_000, (1u64 << 42) - 1];
        for _ in 0..4_090 {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            keys.push(state % 6_000);
        }
        keys.sort_unstable();
        keys.dedup();
        let scalar = keys
            .iter()
            .map(|&key| reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(reader.find_seeds_batch(&keys).unwrap(), scalar);
        let duplicates = [170, 170, 5_000, 170];
        let scalar = duplicates
            .iter()
            .map(|&key| reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(reader.find_seeds_batch(&duplicates).unwrap(), scalar);
        assert_eq!(
            reader
                .find_seeds_batch(&vec![5_000; SEED_LOOKUP_BATCH_KEYS])
                .unwrap(),
            vec![None; SEED_LOOKUP_BATCH_KEYS]
        );
        assert!(
            reader
                .find_seeds_batch(&vec![0; SEED_LOOKUP_BATCH_KEYS + 1])
                .is_err()
        );
        let invalid = 1u64 << 42;
        assert!(matches!(
            reader.find_seed(invalid),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "packed seed key"
            )))
        ));
        assert!(matches!(
            reader.find_seeds_batch(&[invalid]),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "packed seed key"
            )))
        ));
    }

    #[test]
    fn exact_block_lookup_matches_scalar_with_boundaries_gaps_and_duplicates() {
        let (_directory, path) = multi_page_filter_fixture();
        let mut reader = JidxReader::open(&path).unwrap();
        let keys = [4_999, 512, 511, 513, 0, 5_000, 1_024, 511, 4_998, 6_000];
        let scalar = keys
            .iter()
            .map(|&key| reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        reader.exact_block_fence = Some(Arc::new(exact_fence(&path)));
        assert_eq!(reader.find_seeds_batch(&keys).unwrap(), scalar);
    }

    #[test]
    fn resident_filter_is_lazy_reused_and_matches_reference() {
        let (_directory, path) = multi_page_filter_fixture();
        let control = JidxReader::open(&path).unwrap();
        let keys = [0, 1, 170, 2_500, 4_999, 5_000, 6_000];
        let expected = keys
            .iter()
            .map(|&key| {
                crate::jidx_filters::reference_contains(
                    &control,
                    control.filter_directory.as_ref().unwrap(),
                    key,
                )
            })
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        control.verify_query_filter_pages(keys.iter()).unwrap();
        assert!(control.find_seed(170).unwrap().is_some());
        assert_eq!(control.find_seeds_batch(&[170, 6_000]).unwrap()[1], None);
        assert!(control.filter_cache_identity.get().is_none());
        assert!(control.resident_filter_body.get().is_none());

        let fence = Arc::new(exact_fence(&path));
        let mut first = JidxReader::open(&path).unwrap();
        first.exact_block_fence = Some(Arc::clone(&fence));
        assert!(first.resident_filter_body.get().is_none());
        let verified_before = first
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        first.verify_query_filter_pages(keys.iter()).unwrap();
        assert_eq!(
            first
                .verified_pages
                .iter()
                .map(|word| word.load(Ordering::Relaxed))
                .collect::<Vec<_>>(),
            verified_before
        );
        let first_body = first
            .resident_filter_body
            .get()
            .and_then(Option::as_ref)
            .unwrap();
        let actual = keys
            .iter()
            .map(|&key| {
                crate::jidx_filters::contains(&first, first.filter_directory.as_ref().unwrap(), key)
            })
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(actual, expected);

        let mut second = JidxReader::open(&path).unwrap();
        second.exact_block_fence = Some(fence);
        second.verify_query_filter_pages(keys.iter()).unwrap();
        let second_body = second
            .resident_filter_body
            .get()
            .and_then(Option::as_ref)
            .unwrap();
        assert!(Arc::ptr_eq(first_body, second_body));

        let mut scalar = JidxReader::open(&path).unwrap();
        scalar.exact_block_fence = Some(Arc::new(exact_fence(&path)));
        assert!(scalar.resident_filter_body.get().is_none());
        assert!(scalar.find_seed(6_000).unwrap().is_none());
        assert!(Arc::ptr_eq(
            first_body,
            scalar
                .resident_filter_body
                .get()
                .and_then(Option::as_ref)
                .unwrap()
        ));
    }

    #[test]
    fn resident_filter_corruption_never_falls_back() {
        for final_page_padding in [false, true] {
            let (_directory, path) = multi_page_filter_fixture();
            let fence = Arc::new(exact_fence(&path));
            let mut reader = JidxReader::open(&path).unwrap();
            reader.exact_block_fence = Some(fence);
            let filter = reader.header.section(SectionKind::SeedFilters);
            let corruption = if final_page_padding {
                let end = filter.offset + filter.length;
                assert!(end < reader.header.section(SectionKind::BlockChecksums).offset);
                end
            } else {
                filter.offset + PAGE_SIZE + 20
            };
            corrupt_byte(&path, corruption);
            assert!(matches!(
                reader.verify_query_filter_pages([&0]),
                Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
            ));
            assert!(reader.resident_filter_body.get().is_none());
        }
    }

    #[test]
    fn resident_filter_concurrent_first_use_reuses_one_body() {
        let (_directory, path) = multi_page_filter_fixture();
        let fence = Arc::new(exact_fence(&path));
        let barrier = Arc::new(std::sync::Barrier::new(3));
        let mut readers = Vec::new();
        for _ in 0..2 {
            let mut reader = JidxReader::open(&path).unwrap();
            reader.exact_block_fence = Some(Arc::clone(&fence));
            readers.push(reader);
        }
        let handles = readers
            .into_iter()
            .map(|reader| {
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    barrier.wait();
                    reader.verify_query_filter_pages([&170]).unwrap();
                    Arc::clone(
                        reader
                            .resident_filter_body
                            .get()
                            .and_then(Option::as_ref)
                            .unwrap(),
                    )
                })
            })
            .collect::<Vec<_>>();
        barrier.wait();
        let bodies = handles
            .into_iter()
            .map(|handle| handle.join().unwrap())
            .collect::<Vec<_>>();
        assert!(Arc::ptr_eq(&bodies[0], &bodies[1]));
    }

    #[test]
    fn resident_filter_identity_change_is_an_error() {
        let (_directory, path) = multi_page_filter_fixture();
        let fence = Arc::new(exact_fence(&path));
        let mut reader = JidxReader::open(&path).unwrap();
        reader.exact_block_fence = Some(fence);
        assert!(reader.filter_cache_identity().is_some());
        let filter = reader.header.section(SectionKind::SeedFilters);
        corrupt_byte(&path, filter.offset + PAGE_SIZE + 20);
        assert!(matches!(
            reader.verify_query_filter_pages([&0]),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "resident seed filter identity"
            )))
        ));
    }

    #[test]
    fn exact_block_fence_rejects_header_binding_body_and_order_changes() {
        let (_directory, path) = multi_page_filter_fixture();
        let (source, header, valid) = exact_fence_bytes(&path);
        for offset in [12, 16, 48, 64, 96, 128, 136, 168, 4064, 4096 + 8] {
            let mut bytes = valid.clone();
            bytes[offset] ^= 1;
            assert!(
                decode_exact_block_fence(
                    &bytes,
                    &source[..HEADER_SIZE],
                    source.len() as u64,
                    &header
                )
                .is_err()
            );
        }
        assert!(parse_sha256(&format!("a{}b", "é".repeat(31))).is_err());
    }

    #[test]
    fn exact_block_lookup_preserves_mixed_seed_families() {
        let (_directory, path) =
            external_fixture_with_keys([0x1234, crate::jidx::RESCUE_K15_TAG | 0x123]);
        replace_bytes(&path, 362, &[1]);
        let mut reader = JidxReader::open(&path).unwrap();
        let keys = [
            crate::jidx::RESCUE_K15_TAG | 0x124,
            0x1234,
            crate::jidx::RESCUE_K15_TAG | 0x123,
            0x1233,
        ];
        let scalar = keys
            .iter()
            .map(|&key| reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        reader.exact_block_fence = Some(Arc::new(exact_fence(&path)));
        assert_eq!(reader.find_seeds_batch(&keys).unwrap(), scalar);
    }

    #[test]
    fn exact_block_lookup_checks_each_page_before_decoding() {
        for page in 0..3 {
            let (_directory, path) = multi_page_filter_fixture();
            let fence = Arc::new(exact_fence(&path));
            let seeds = section_offset(&path, SectionKind::Seeds);
            corrupt_byte(&path, seeds + page * PAGE_SIZE + 2_000);
            let mut reader = JidxReader::open(&path).unwrap();
            reader.exact_block_fence = Some(Arc::clone(&fence));
            assert!(matches!(
                reader.find_seeds_batch(&[0]),
                Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
            ));
        }
    }

    #[test]
    #[ignore = "requires JAM_EXACT_BLOCK_SMOKE_SOURCE and the two fence environment variables"]
    fn exact_block_external_builder_reader_cross_smoke() {
        let source = std::env::var_os("JAM_EXACT_BLOCK_SMOKE_SOURCE").unwrap();
        let reader = JidxReader::open(source).unwrap();
        assert!(reader.exact_block_fence().is_some());
        let mut keys = (0..reader.header().seed_count)
            .map(|ordinal| {
                crate::jidx_postings::entry(&reader, ordinal)
                    .unwrap()
                    .packed_key
            })
            .collect::<Vec<_>>();
        let present = keys.clone();
        for window in present.windows(2) {
            if let Some(candidate) = window[0].checked_add(1)
                && candidate < window[1]
                && seed_length(reader.header().k, reader.header().rescue_k15, candidate).is_ok()
            {
                keys.push(candidate);
            }
        }
        if let Some(first) = present.first().copied().filter(|first| *first > 0) {
            let candidate = first - 1;
            if seed_length(reader.header().k, reader.header().rescue_k15, candidate).is_ok() {
                keys.push(candidate);
            }
        }
        if let Some(candidate) = present.last().and_then(|last| last.checked_add(1))
            && seed_length(reader.header().k, reader.header().rescue_k15, candidate).is_ok()
        {
            keys.push(candidate);
        }
        keys.extend(present.first().copied());
        keys.extend(present.last().copied());
        keys.reverse();
        for chunk in keys.chunks(SEED_LOOKUP_BATCH_KEYS) {
            let scalar = chunk
                .iter()
                .map(|&key| reader.find_seed(key))
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            assert_eq!(reader.find_seeds_batch(chunk).unwrap(), scalar);
        }
        assert!(
            reader
                .resident_filter_body
                .get()
                .and_then(Option::as_ref)
                .is_some()
        );
    }

    #[test]
    fn batch_seed_lookup_matches_scalar_for_mixed_seed_families() {
        let (_directory, path) =
            external_fixture_with_keys([0x1234, crate::jidx::RESCUE_K15_TAG | 0x123]);
        replace_bytes(&path, 362, &[1]);
        let reader = JidxReader::open(path).unwrap();
        let mut keys = vec![
            0,
            0x1233,
            0x1234,
            0x1235,
            crate::jidx::RESCUE_K15_TAG | 0x122,
            crate::jidx::RESCUE_K15_TAG | 0x123,
            crate::jidx::RESCUE_K15_TAG | 0x124,
            crate::jidx::RESCUE_K15_TAG | ((1u64 << 30) - 1),
        ];
        keys.sort_unstable();
        let scalar = keys
            .iter()
            .map(|&key| reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(reader.find_seeds_batch(&keys).unwrap(), scalar);
    }

    #[test]
    fn batch_filter_negatives_do_not_touch_seed_table() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(&path).unwrap();
        let seeds = reader.header.section(SectionKind::Seeds);
        drop(reader);
        corrupt_byte(&path, seeds.offset);

        let reader = JidxReader::open(path).unwrap();
        assert_eq!(reader.find_seeds_batch(&[5_000]).unwrap(), vec![None]);
    }

    #[test]
    fn batch_seed_lookup_matches_scalar_pages_and_error_chain() {
        let (_directory, path) = multi_page_filter_fixture();
        let keys = [0, 170, 1_234, 2_500, 4_999, 5_000];
        let scalar_reader = JidxReader::open(&path).unwrap();
        let scalar = keys
            .iter()
            .map(|&key| scalar_reader.find_seed(key))
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        let scalar_pages = scalar_reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();

        let batch_reader = JidxReader::open(&path).unwrap();
        assert_eq!(batch_reader.find_seeds_batch(&keys).unwrap(), scalar);
        let batch_pages = batch_reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        assert_eq!(batch_pages, scalar_pages);

        let seeds = batch_reader.header.section(SectionKind::Seeds);
        let first_page = seeds.offset / PAGE_SIZE;
        let last_page = (seeds.offset + seeds.length - 1) / PAGE_SIZE;
        assert!((first_page..=last_page).any(|page| {
            let page_bit = page - 1;
            batch_pages[(page_bit / 64) as usize] & (1u64 << (page_bit % 64)) == 0
        }));
        let probe_reader = JidxReader::open(&path).unwrap();
        let before = probe_reader
            .verified_pages
            .iter()
            .map(|word| word.load(Ordering::Relaxed))
            .collect::<Vec<_>>();
        probe_reader.find_seed(keys[1]).unwrap();
        let accessed_page = (first_page..=last_page)
            .find(|&page| {
                let page_bit = page - 1;
                let mask = 1u64 << (page_bit % 64);
                before[(page_bit / 64) as usize] & mask == 0
                    && probe_reader.verified_pages[(page_bit / 64) as usize].load(Ordering::Relaxed)
                        & mask
                        != 0
            })
            .unwrap();
        drop(probe_reader);
        drop(batch_reader);
        drop(scalar_reader);
        corrupt_byte(&path, accessed_page * PAGE_SIZE);

        let scalar_reader = JidxReader::open(&path).unwrap();
        assert!(matches!(
            scalar_reader.find_seed(keys[1]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
        let batch_reader = JidxReader::open(path).unwrap();
        assert!(matches!(
            batch_reader.find_seeds_batch(&[keys[1]]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
    }

    #[test]
    fn batch_seed_lookup_preserves_cross_page_reservation_error() {
        let (_directory, path) = multi_page_filter_fixture();
        let reader = JidxReader::open(&path).unwrap();
        let seeds = reader.header.section(SectionKind::Seeds);
        drop(reader);
        let record_start = seeds.offset + 170 * u64::from(SEED_RECORD_SIZE);
        assert_eq!(record_start % PAGE_SIZE, PAGE_SIZE - 16);
        replace_bytes(&path, record_start + 20, &1u32.to_le_bytes());
        rewrite_checksums(&path);

        let scalar_reader = JidxReader::open(&path).unwrap();
        assert!(matches!(
            scalar_reader.find_seed(170),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "seed reservation"
            )))
        ));
        let batch_reader = JidxReader::open(path).unwrap();
        assert!(matches!(
            batch_reader.find_seeds_batch(&[170]),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "seed reservation"
            )))
        ));
    }

    #[test]
    fn corrupt_filter_directory_and_payload_remain_errors() {
        let (_directory, directory_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&directory_path, SectionKind::SeedFilters);
        corrupt_byte(&directory_path, filter_offset);
        assert!(JidxReader::open(&directory_path).is_err());
        assert!(JidxReader::open(&directory_path).is_err());

        let (_directory, payload_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&payload_path, SectionKind::SeedFilters);
        corrupt_byte(&payload_path, filter_offset + PAGE_SIZE);
        let reader = JidxReader::open(payload_path).unwrap();
        assert!(reader.find_seed(0x1234).is_err());
        assert!(reader.find_seed(0x1234).is_err());
    }

    #[test]
    fn checksummed_filter_structure_and_audit_fail_closed() {
        let (_directory, range_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&range_path, SectionKind::SeedFilters);
        replace_bytes(
            &range_path,
            filter_offset + 16,
            &(PAGE_SIZE + 1).to_le_bytes(),
        );
        rewrite_checksums(&range_path);
        assert!(JidxReader::open(range_path).is_err());

        let (_directory, key_range_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&key_range_path, SectionKind::SeedFilters);
        replace_bytes(&key_range_path, filter_offset, &0x1235u64.to_le_bytes());
        replace_bytes(&key_range_path, filter_offset + 8, &0x1235u64.to_le_bytes());
        rewrite_checksums(&key_range_path);
        let reader = JidxReader::open(key_range_path).unwrap();
        assert!(reader.verify_checksum().is_err());

        let (_directory, descriptor_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&descriptor_path, SectionKind::SeedFilters);
        replace_bytes(
            &descriptor_path,
            filter_offset + PAGE_SIZE + 8,
            &0u32.to_le_bytes(),
        );
        rewrite_checksums(&descriptor_path);
        let reader = JidxReader::open(descriptor_path).unwrap();
        assert!(reader.find_seed(0x1234).is_err());
        assert!(reader.find_seed(0x1234).is_err());

        let (_directory, padding_path) = fixture(false, false, false, false);
        let filter_offset = section_offset(&padding_path, SectionKind::SeedFilters);
        replace_bytes(&padding_path, filter_offset + 40, &[1]);
        rewrite_checksums(&padding_path);
        let reader = JidxReader::open(padding_path).unwrap();
        assert!(reader.verify_checksum().is_err());

        let (_directory, mismatch_path) = fixture(false, false, false, false);
        let wrong_key = (0..1_000_000u64)
            .find(|key| {
                *key != 0x1234
                    && !BinaryFuse8::try_from(&[*key][..])
                        .unwrap()
                        .contains(&0x1234)
            })
            .unwrap();
        let wrong_filter = seed_filters(&[wrong_key]);
        let (filter_offset, filter_length) =
            section_range(&mismatch_path, SectionKind::SeedFilters);
        assert_eq!(wrong_filter.len() as u64, filter_length);
        replace_bytes(
            &mismatch_path,
            filter_offset + PAGE_SIZE,
            &wrong_filter[PAGE_SIZE as usize..],
        );
        rewrite_checksums(&mismatch_path);
        let reader = JidxReader::open(mismatch_path).unwrap();
        assert!(reader.verify_checksum().is_err());
    }

    #[test]
    fn multi_page_filter_checks_accessed_pages_and_audits_untouched_pages() {
        let key = 1_234u64;
        let (_directory, accessed_path) = multi_page_filter_fixture();
        let (_, filter_length) = section_range(&accessed_path, SectionKind::SeedFilters);
        assert!(filter_length > 2 * PAGE_SIZE);
        let reader = JidxReader::open(&accessed_path).unwrap();
        let offsets = crate::jidx_filters::fingerprint_offsets(
            &reader,
            reader.filter_directory.as_ref().unwrap(),
            key,
        )
        .unwrap();
        drop(reader);
        corrupt_byte(&accessed_path, offsets[0]);
        let reader = JidxReader::open(accessed_path).unwrap();
        assert!(reader.find_seed(key).is_err());
        assert!(reader.find_seed(key).is_err());

        let (_directory, untouched_path) = multi_page_filter_fixture();
        let reader = JidxReader::open(&untouched_path).unwrap();
        let directory = reader.filter_directory.as_ref().unwrap();
        let accessed_pages = crate::jidx_filters::fingerprint_offsets(&reader, directory, key)
            .unwrap()
            .map(|offset| offset / PAGE_SIZE);
        let untouched = (0..5_000u64)
            .filter(|candidate| *candidate != key)
            .find_map(|candidate| {
                crate::jidx_filters::fingerprint_offsets(&reader, directory, candidate)
                    .unwrap()
                    .into_iter()
                    .find(|offset| !accessed_pages.contains(&(offset / PAGE_SIZE)))
            })
            .unwrap();
        drop(reader);
        corrupt_byte(&untouched_path, untouched);
        rewrite_checksums(&untouched_path);
        let reader = JidxReader::open(untouched_path).unwrap();
        assert!(reader.find_seed(key).unwrap().is_some());
        assert!(reader.verify_checksum().is_err());
    }

    #[test]
    fn empty_seed_table_has_an_empty_filter_directory() {
        let (_directory, path) = fixture_with_string_pages(false, false, false, false, false, None);
        let reader = JidxReader::open(path).unwrap();
        assert_eq!(reader.header().seed_count, 0);
        assert!(reader.find_seed(0).unwrap().is_none());
        reader.verify_checksum().unwrap();
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

    #[test]
    fn golden_inline_and_external_groups_decode_independently() {
        let (_directory, path) = external_fixture();
        let reader = JidxReader::open(path).unwrap();
        reader.verify_checksum().unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        let documents = reader.seed_documents(seed).unwrap();
        assert_eq!(
            documents
                .iter()
                .map(|document| (document.metagenome_id, document.occurrence_count))
                .collect::<Vec<_>>(),
            [(0, 1), (1, 2)]
        );
        assert_eq!(
            reader
                .seed_document_occurrences(seed, documents[0])
                .unwrap(),
            [SeedOccurrence {
                contig_id: 1,
                position: 2,
                canonical_orientation: true,
            }]
        );
        assert_eq!(
            reader
                .seed_document_occurrences(seed, documents[1])
                .unwrap(),
            [
                SeedOccurrence {
                    contig_id: 2,
                    position: 5,
                    canonical_orientation: false,
                },
                SeedOccurrence {
                    contig_id: 2,
                    position: 9,
                    canonical_orientation: true,
                },
            ]
        );

        let large = reader.find_seed(0x1235).unwrap().unwrap();
        assert_eq!(
            reader.seed_occurrences(large).unwrap()[0].position,
            (1u64 << 63) + 5
        );
        assert!(
            reader
                .seed_document_occurrences(large, documents[1])
                .is_err()
        );
    }

    #[test]
    fn document_metadata_does_not_decode_an_untouched_cold_group() {
        let (_directory, path) = external_fixture();
        let mut file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.seek(SeekFrom::Start(6 * PAGE_SIZE)).unwrap();
        file.write_all(&[0xff]).unwrap();
        drop(file);

        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        let documents = reader.seed_documents(seed).unwrap();
        assert_eq!(documents[1].occurrence_count, 2);
        assert!(matches!(
            reader.seed_document_occurrences(seed, documents[1]),
            Err(JidxReaderError::Format(JidxError::ChecksumMismatch))
        ));
    }

    #[test]
    fn document_metadata_rejects_per_seed_count_above_header_total() {
        let (_directory, path) = external_fixture();
        let mut file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.seek(SeekFrom::Start(40)).unwrap();
        file.write_all(&2u64.to_le_bytes()).unwrap();
        file.seek(SeekFrom::Start(6 * PAGE_SIZE)).unwrap();
        file.write_all(&[0xff]).unwrap();
        drop(file);

        let reader = JidxReader::open(path).unwrap();
        let seed = reader.find_seed(0x1234).unwrap().unwrap();
        assert!(matches!(
            reader.seed_documents(seed),
            Err(JidxReaderError::Format(JidxError::Invalid(
                "occurrence count"
            )))
        ));
    }
}
