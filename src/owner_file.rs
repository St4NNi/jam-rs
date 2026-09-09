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

    pub(crate) fn find_seeds_batch(
        &self,
        keys: &[u64],
    ) -> Result<Vec<Option<OwnerSeed>>, OwnerReaderError> {
        let mut out = vec![None; keys.len()];
        let mut grouped = std::collections::BTreeMap::<u64, Vec<(usize, u64)>>::new();
        for (position, &key) in keys.iter().enumerate() {
            if key < self.header.first_key || key > self.header.last_key {
                return Err(OwnerReaderError::KeyNotCovered(key));
            }
            let block = self.find_block(key)?;
            let Some((block_ordinal, _)) = block else {
                if !self.header.complete_range() {
                    return Err(OwnerReaderError::KeyNotCovered(key));
                }
                continue;
            };
            grouped
                .entry(block_ordinal)
                .or_default()
                .push((position, key));
        }
        for (block_ordinal, requests) in grouped {
            let record = self.validated_block_record(block_ordinal)?;
            let decoded = self.hot_block(block_ordinal, record)?;
            for (position, key) in requests {
                out[position] = decoded
                    .binary_search_by_key(&key, |entry| entry.key)
                    .ok()
                    .map(|index| OwnerSeed {
                        packed_key: key,
                        document_frequency: decoded[index].document_frequency,
                        block_ordinal,
                    });
                if out[position].is_none() && !self.header.complete_range() {
                    return Err(OwnerReaderError::KeyNotCovered(key));
                }
            }
        }
        Ok(out)
    }

    pub(crate) fn seed_documents(
        &self,
        seed: OwnerSeed,
    ) -> Result<Vec<OwnerDocument>, OwnerReaderError> {
        let record = self.validated_block_record(seed.block_ordinal)?;
        let decoded = self.hot_block(seed.block_ordinal, record)?;
        let hot = decoded
            .binary_search_by_key(&seed.packed_key, |entry| entry.key)
            .ok()
            .map(|index| &decoded[index])
            .ok_or(OwnerReaderError::Invalid("missing cached seed"))?;
        if hot.document_frequency != seed.document_frequency {
            return Err(OwnerReaderError::Invalid("seed frequency"));
        }
        hot.members
            .iter()
            .enumerate()
            .map(|(member_ordinal, member)| {
                Ok(OwnerDocument {
                    metagenome_id: member.document_id,
                    occurrence_count: member.occurrence_count,
                    seed_key: seed.packed_key,
                    owner_ordinal: self.header.owner_ordinal,
                    block_ordinal: seed.block_ordinal,
                    member_ordinal: u64::try_from(member_ordinal)
                        .map_err(|_| OwnerReaderError::Invalid("member count"))?,
                })
            })
            .collect()
    }

    pub(crate) fn document_occurrences(
        &self,
        document: OwnerDocument,
    ) -> Result<Vec<crate::owner_postings::OwnerOccurrence>, OwnerReaderError> {
        if document.owner_ordinal != self.header.owner_ordinal
            || document.block_ordinal >= self.block_count()?
        {
            return Err(OwnerReaderError::Invalid("seed document"));
        }
        let record = self.validated_block_record(document.block_ordinal)?;
        let decoded = self.hot_block(document.block_ordinal, record)?;
        let hot = decoded
            .binary_search_by_key(&document.seed_key, |entry| entry.key)
            .ok()
            .map(|index| &decoded[index])
            .ok_or(OwnerReaderError::Invalid("missing seed document"))?;
        let member: &OwnerHotMember = hot
            .members
            .get(
                usize::try_from(document.member_ordinal)
                    .map_err(|_| OwnerReaderError::Invalid("member count"))?,
            )
            .filter(|member| {
                member.document_id == document.metagenome_id
                    && member.occurrence_count == document.occurrence_count
            })
            .ok_or(OwnerReaderError::Invalid("seed document"))?;
        let decoded_bytes = usize::try_from(member.occurrence_count)
            .ok()
            .and_then(|count| {
                count.checked_mul(std::mem::size_of::<crate::owner_postings::OwnerOccurrence>())
            })
            .ok_or(OwnerReaderError::Invalid("decoded occurrence bytes"))?;
        if decoded_bytes > MAX_DECODED_OCCURRENCE_BYTES {
            return Err(OwnerReaderError::Invalid("decoded occurrence bytes"));
        }
        let member_end = member
            .cold_offset
            .checked_add(member.cold_length)
            .filter(|end| *end <= record.cold_length)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let cold = self.header.section(OwnerSection::ColdPostings);
        let block_start = cold
            .offset
            .checked_add(record.cold_offset)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let start = block_start
            .checked_add(member.cold_offset)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let end = block_start
            .checked_add(member_end)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let occurrences = decode_member(
            self.checked_bytes(start, end)?,
            OwnerHotMember {
                cold_offset: 0,
                ..*member
            },
        )?;
        self.verify_unchanged()?;
        Ok(occurrences)
    }

    fn find_block(&self, key: u64) -> Result<Option<(u64, BlockRecord)>, OwnerReaderError> {
        let count = self.block_count()?;
        let mut low = 0;
        let mut high = count;
        while low < high {
            let middle = low + (high - low) / 2;
            let record = self.validated_block_record(middle)?;
            if record.last_key < key {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        if low == count {
            return Ok(None);
        }
        let record = self.validated_block_record(low)?;
        Ok((record.first_key <= key).then_some((low, record)))
    }

    pub(crate) fn block_count(&self) -> Result<u64, OwnerReaderError> {
        Ok(self.header.section(OwnerSection::BlockDirectory).length / u64::from(OWNER_BLOCK_SIZE))
    }

    pub(crate) fn block_record(&self, ordinal: u64) -> Result<BlockRecord, OwnerReaderError> {
        if ordinal >= self.block_count()? {
            return Err(OwnerReaderError::Invalid("block ordinal"));
        }
        let section = self.header.section(OwnerSection::BlockDirectory);
        let start = section.offset + ordinal * u64::from(OWNER_BLOCK_SIZE);
        let bytes = self.checked_bytes(start, start + u64::from(OWNER_BLOCK_SIZE))?;
        BlockRecord::decode(bytes)
    }

    fn validated_block_record(&self, ordinal: u64) -> Result<BlockRecord, OwnerReaderError> {
        let record = self.block_record(ordinal)?;
        let hot = self.header.section(OwnerSection::HotPostings);
        let cold = self.header.section(OwnerSection::ColdPostings);
        if record.first_key > record.last_key
            || record.first_key < self.header.first_key
            || record.last_key > self.header.last_key
            || record.key_count == 0
            || record.key_count as usize > MAX_KEYS_PER_BLOCK
            || record
                .hot_offset
                .checked_add(record.hot_length)
                .is_none_or(|end| end > hot.length)
            || record
                .cold_offset
                .checked_add(record.cold_length)
                .is_none_or(|end| end > cold.length)
        {
            return Err(OwnerReaderError::Invalid("block directory"));
        }
        if ordinal > 0 && self.block_record(ordinal - 1)?.last_key >= record.first_key {
            return Err(OwnerReaderError::Invalid("block directory"));
        }
        if ordinal + 1 < self.block_count()?
            && record.last_key >= self.block_record(ordinal + 1)?.first_key
        {
            return Err(OwnerReaderError::Invalid("block directory"));
        }
        Ok(record)
    }

    fn block_bytes(&self, block: BlockRecord, hot: bool) -> Result<&[u8], OwnerReaderError> {
        let (section, offset, length) = if hot {
            (
                OwnerSection::HotPostings,
                block.hot_offset,
                block.hot_length,
            )
        } else {
            (
                OwnerSection::ColdPostings,
                block.cold_offset,
                block.cold_length,
            )
        };
        let section = self.header.section(section);
        let start = section
            .offset
            .checked_add(offset)
            .ok_or(OwnerReaderError::Invalid("block range"))?;
        let end = start
            .checked_add(length)
            .ok_or(OwnerReaderError::Invalid("block range"))?;
        if end > section.offset + section.length {
            return Err(OwnerReaderError::Invalid("block range"));
        }
        self.checked_bytes(start, end)
    }

    pub(crate) fn validate_directory(&self) -> Result<(), OwnerReaderError> {
        let mut previous = None;
        let mut keys = 0u64;
        let mut hot_offset = 0u64;
        let mut cold_offset = 0u64;
        for ordinal in 0..self.block_count()? {
            let block = self.block_record(ordinal)?;
            if block.first_key > block.last_key
                || block.first_key < self.header.first_key
                || block.last_key > self.header.last_key
                || previous.is_some_and(|key| key >= block.first_key)
                || block.hot_offset != hot_offset
                || block.cold_offset != cold_offset
                || block.key_count == 0
                || block.key_count as usize > MAX_KEYS_PER_BLOCK
            {
                return Err(OwnerReaderError::Invalid("block directory"));
            }
            keys += u64::from(block.key_count);
            hot_offset = hot_offset
                .checked_add(block.hot_length)
                .ok_or(OwnerReaderError::Invalid("block range"))?;
            cold_offset = cold_offset
                .checked_add(block.cold_length)
                .ok_or(OwnerReaderError::Invalid("block range"))?;
            previous = Some(block.last_key);
        }
        if keys != self.header.key_count
            || hot_offset != self.header.section(OwnerSection::HotPostings).length
            || cold_offset != self.header.section(OwnerSection::ColdPostings).length
        {
            return Err(OwnerReaderError::Invalid("block coverage"));
        }
        Ok(())
    }

    fn parse_hot_bounded(&self, hot: &[u8]) -> Result<Vec<OwnerHotKey>, OwnerReaderError> {
        let decoded_bound = decoded_hot_bound(hot)?;
        if decoded_bound > MAX_DECODED_HOT_BYTES {
            return Err(OwnerReaderError::Invalid("decoded hot postings"));
        }
        Ok(parse_hot(hot)?)
    }

    pub(crate) fn hot_block(
        &self,
        ordinal: u64,
        record: BlockRecord,
    ) -> Result<Arc<Vec<OwnerHotKey>>, OwnerReaderError> {
        self.verify_unchanged()?;
        if let Some(decoded) = self
            .hot_cache
            .lock()
            .map_err(|_| OwnerReaderError::Invalid("hot cache"))?
            .blocks
            .get(&ordinal)
            .cloned()
        {
            return Ok(decoded);
        }
        let hot = self.block_bytes(record, true)?;
        let charge = decoded_hot_bound(hot)?;
        let decoded = Arc::new(self.parse_hot_bounded(hot)?);
        self.verify_unchanged()?;
        if decoded.len() != record.key_count as usize
            || decoded.first().map(|entry| entry.key) != Some(record.first_key)
            || decoded.last().map(|entry| entry.key) != Some(record.last_key)
        {
            return Err(OwnerReaderError::Invalid("block contents"));
        }
        let mut cache = self
            .hot_cache
            .lock()
            .map_err(|_| OwnerReaderError::Invalid("hot cache"))?;
        if let Some(existing) = cache.blocks.get(&ordinal) {
            return Ok(Arc::clone(existing));
        }
        if self.file_identity.is_some()
            && self
                .hot_cache_bytes
                .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |bytes| {
                    bytes
                        .checked_add(charge)
                        .filter(|total| *total <= MAX_HOT_CACHE_BYTES)
                })
                .is_ok()
        {
            cache.blocks.insert(ordinal, Arc::clone(&decoded));
        }
        Ok(decoded)
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

fn decoded_hot_bound(hot: &[u8]) -> Result<usize, OwnerReaderError> {
    hot.len()
        .checked_div(3)
        .and_then(|members| members.checked_mul(std::mem::size_of::<OwnerHotMember>()))
        .and_then(|bytes| {
            MAX_KEYS_PER_BLOCK
                .checked_mul(std::mem::size_of::<OwnerHotKey>())
                .and_then(|keys| bytes.checked_add(keys))
        })
        .ok_or(OwnerReaderError::Invalid("decoded hot postings"))
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
