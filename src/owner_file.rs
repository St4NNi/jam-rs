use crate::jidx::sha256;
use crate::jidx_reader::{Contig, Metagenome};
use crate::owner_format::{
    BlockRecord, DocumentRecord, OWNER_BLOCK_SIZE, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE,
    OWNER_HEADER_SIZE, OWNER_PAGE_SIZE, OwnerDocument, OwnerHeader, OwnerReaderError, OwnerSection,
    OwnerSeed, checksum_layout, read_u32, read_u64,
};
use crate::owner_observer::{self, OwnerReadObserver};
use crate::owner_postings::{
    MAX_KEYS_PER_BLOCK, OwnerAnchor, OwnerHotBlock, OwnerHotKey, OwnerHotMember,
    decode_member_window, find_key, key_members, locate_member, parse_hot,
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
    observer: Arc<OwnerReadObserver>,
}

#[derive(Default)]
struct HotCache {
    blocks: std::collections::BTreeMap<u64, Arc<OwnerHotBlock>>,
}

impl OwnerFile {
    pub(crate) fn open(
        path: PathBuf,
        hot_cache_bytes: Arc<AtomicUsize>,
        observer: Arc<OwnerReadObserver>,
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
            observer,
        };
        owner
            .observer
            .record_open(0, 1, OWNER_HEADER_SIZE as u64, 0, 1);
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
                out[position] = find_key(&decoded, key).map(|entry| OwnerSeed {
                    packed_key: key,
                    document_frequency: entry.document_frequency,
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
        let hot = find_key(&decoded, seed.packed_key)
            .ok_or(OwnerReaderError::Invalid("missing cached seed"))?;
        if hot.document_frequency != seed.document_frequency {
            return Err(OwnerReaderError::Invalid("seed frequency"));
        }
        key_members(&decoded, hot)?
            .iter()
            .map(|member| {
                Ok(OwnerDocument {
                    metagenome_id: member.document_id,
                    occurrence_count: member.occurrence_count,
                    seed_key: seed.packed_key,
                    owner_ordinal: self.header.owner_ordinal,
                    block_ordinal: seed.block_ordinal,
                    member_ordinal: member.member_ordinal,
                })
            })
            .collect()
    }

    pub(crate) fn document_occurrences(
        &self,
        document: OwnerDocument,
        document_widths: &[u8],
    ) -> Result<Vec<(u64, bool)>, OwnerReaderError> {
        if document.owner_ordinal != self.header.owner_ordinal
            || document.block_ordinal >= self.block_count()?
        {
            return Err(OwnerReaderError::Invalid("seed document"));
        }
        let record = self.validated_block_record(document.block_ordinal)?;
        let decoded = self.hot_block(document.block_ordinal, record)?;
        let hot = find_key(&decoded, document.seed_key)
            .ok_or(OwnerReaderError::Invalid("missing seed document"))?;
        let member: &OwnerHotMember = key_members(&decoded, hot)?
            .get(
                usize::try_from(
                    document
                        .member_ordinal
                        .checked_sub(hot.member_start)
                        .ok_or(OwnerReaderError::Invalid("seed document"))?,
                )
                .map_err(|_| OwnerReaderError::Invalid("member count"))?,
            )
            .filter(|member| {
                member.document_id == document.metagenome_id
                    && member.occurrence_count == document.occurrence_count
            })
            .ok_or(OwnerReaderError::Invalid("seed document"))?;
        if member
            .occurrence_count
            .checked_mul(std::mem::size_of::<crate::owner_postings::OwnerOccurrence>() as u64)
            .is_none_or(|bytes| bytes > MAX_DECODED_OCCURRENCE_BYTES as u64)
        {
            return Err(OwnerReaderError::Invalid("decoded occurrence bytes"));
        }
        let window = locate_member(&decoded, member.member_ordinal)?;
        let window_start = window.start_bit / 8;
        let window_end = window.end_bit.div_ceil(8);
        if window_end > record.cold_length {
            return Err(OwnerReaderError::Invalid("cold member range"));
        }
        let cold = self.header.section(OwnerSection::ColdPostings);
        let block_start = cold
            .offset
            .checked_add(record.cold_offset)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let start = block_start
            .checked_add(window_start)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let end = block_start
            .checked_add(window_end)
            .ok_or(OwnerReaderError::Invalid("cold member range"))?;
        let occurrences = decode_member_window(
            self.checked_bytes(start, end)?,
            window_start,
            window,
            &decoded,
            document_widths,
            MAX_DECODED_OCCURRENCE_BYTES,
        )?;
        self.observer.record_cold(
            end - start,
            occurrences.loci.len() as u64 + occurrences.skipped_occurrences,
        );
        self.verify_unchanged()?;
        Ok(occurrences.loci)
    }

    fn find_block(&self, key: u64) -> Result<Option<(u64, BlockRecord)>, OwnerReaderError> {
        let count = self.block_count()?;
        let mut low = 0;
        let mut high = count;
        while low < high {
            self.observer.record_directory(1, 0, 0);
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
        self.observer
            .record_directory(0, 1, u64::from(OWNER_BLOCK_SIZE));
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

    fn parse_hot_bounded(&self, hot: &[u8]) -> Result<OwnerHotBlock, OwnerReaderError> {
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
    ) -> Result<Arc<OwnerHotBlock>, OwnerReaderError> {
        self.verify_unchanged()?;
        if let Some(decoded) = self
            .hot_cache
            .lock()
            .map_err(|_| OwnerReaderError::Invalid("hot cache"))?
            .blocks
            .get(&ordinal)
            .cloned()
        {
            self.observer
                .record_hot_request(self.header.owner_ordinal, ordinal, 0, true);
            return Ok(decoded);
        }
        self.observer.record_hot_request(
            self.header.owner_ordinal,
            ordinal,
            record.hot_length,
            false,
        );
        let hot = self.block_bytes(record, true)?;
        let charge = decoded_hot_bound(hot)?;
        let decoded = Arc::new(self.parse_hot_bounded(hot)?);
        self.observer
            .record_hot_decode(charge as u64, decoded.members.len() as u64);
        self.verify_unchanged()?;
        if decoded.keys.len() != record.key_count as usize
            || decoded.keys.first().map(|entry| entry.key) != Some(record.first_key)
            || decoded.keys.last().map(|entry| entry.key) != Some(record.last_key)
            || decoded.document_count != self.header.document_count
            || decoded.cold_bits.div_ceil(8) != record.cold_length
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
        let metadata_sections = [
            OwnerSection::Strings,
            OwnerSection::Documents,
            OwnerSection::Contigs,
            OwnerSection::Gzi,
        ]
        .map(|kind| self.section_bytes(kind, 0, self.header.section(kind).length));
        let [strings, documents, contigs, gzi] = metadata_sections;
        if crate::owner_format::metadata_digest([strings?, documents?, contigs?, gzi?])
            != self.header.metadata_sha256
        {
            return Err(OwnerReaderError::Invalid("metadata digest"));
        }
        let mut document_names = std::collections::HashSet::new();
        let mut expected_contig = 0u32;
        let mut expected_original = 0u32;
        let mut expected_gzi = 0u64;
        for document_id in 0..self.header.document_count {
            let record = self.document_record(document_id)?;
            if record.bgzf_bytes == 0
                || record.bgzf_sha256 == [0; 32]
                || record.original_contig_count == 0
                || record.contig_start != expected_contig
                || record.original_contig_start != expected_original
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
            let mut previous_id = None;
            let mut previous_end = 0;
            for row in record.contig_start..contig_end {
                let contig = self.contig_row(row)?;
                if contig.metagenome_id != document_id
                    || contig.length == 0
                    || contig.line_bases == 0
                    || contig.line_width < contig.line_bases
                    || contig.line_width > contig.line_bases.saturating_add(2)
                    || !contig_names.insert(contig.name)
                    || previous_id.is_some_and(|id| id >= contig.id)
                    || contig.fasta_offset < previous_end
                {
                    return Err(OwnerReaderError::Invalid("contig metadata"));
                }
                previous_id = Some(contig.id);
                let last = contig.length - 1;
                previous_end = (last / u64::from(contig.line_bases))
                    .checked_mul(u64::from(contig.line_width))
                    .and_then(|offset| offset.checked_add(last % u64::from(contig.line_bases)))
                    .and_then(|offset| offset.checked_add(contig.fasta_offset))
                    .and_then(|offset| offset.checked_add(1))
                    .ok_or(OwnerReaderError::Invalid("contig locus"))?;
                if 64 - (previous_end - 1).leading_zeros() > u32::from(record.locus_bits) {
                    return Err(OwnerReaderError::Invalid("document locus width"));
                }
            }
            expected_contig = contig_end;
            expected_original = expected_original
                .checked_add(record.original_contig_count)
                .ok_or(OwnerReaderError::Invalid("contig metadata"))?;
        }
        if u64::from(expected_contig) * u64::from(OWNER_CONTIG_SIZE)
            != self.header.section(OwnerSection::Contigs).length
            || expected_original != self.header.contig_count
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
        self.observer.record_metadata(1, 0, 0, 0);
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
            contig_start: record.original_contig_start,
            contig_count: record.original_contig_count,
            gzi: self.section_bytes(OwnerSection::Gzi, record.gzi_offset, record.gzi_length)?,
        }))
    }

    pub(crate) fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, OwnerReaderError> {
        if id >= self.header.contig_count {
            return Ok(None);
        }
        if self.header.section(OwnerSection::Contigs).length
            == u64::from(self.header.contig_count) * u64::from(OWNER_CONTIG_SIZE)
        {
            return self.contig_row(id).map(Some);
        }
        let mut low = 0u32;
        let mut high = u32::try_from(
            self.header.section(OwnerSection::Contigs).length / u64::from(OWNER_CONTIG_SIZE),
        )
        .map_err(|_| OwnerReaderError::Invalid("contig count"))?;
        while low < high {
            let middle = low + (high - low) / 2;
            if read_u32(self.contig_bytes(middle)?, 12) < id {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        if u64::from(low) * u64::from(OWNER_CONTIG_SIZE)
            == self.header.section(OwnerSection::Contigs).length
        {
            return Ok(None);
        }
        let contig = self.contig_row(low)?;
        Ok((contig.id == id).then_some(contig))
    }

    fn contig_bytes(&self, row: u32) -> Result<&[u8], OwnerReaderError> {
        self.observer.record_metadata(0, 1, 0, 0);
        self.section_bytes(
            OwnerSection::Contigs,
            u64::from(row) * u64::from(OWNER_CONTIG_SIZE),
            u64::from(OWNER_CONTIG_SIZE),
        )
    }

    fn contig_row(&self, row: u32) -> Result<Contig<'_>, OwnerReaderError> {
        let bytes = self.contig_bytes(row)?;
        let id = read_u32(bytes, 12);
        let document = read_u32(bytes, 0);
        let owner = self.document_record(document)?;
        let contig_end = owner
            .original_contig_start
            .checked_add(owner.original_contig_count)
            .ok_or(OwnerReaderError::Invalid("contig ownership"))?;
        if id < owner.original_contig_start
            || id >= contig_end
            || row < owner.contig_start
            || row - owner.contig_start >= owner.contig_count
        {
            return Err(OwnerReaderError::Invalid("contig ownership"));
        }
        Ok(Contig {
            id,
            metagenome_id: document,
            name: self.string(read_u32(bytes, 4), read_u32(bytes, 8))?,
            length: read_u64(bytes, 16),
            fasta_offset: read_u64(bytes, 24),
            line_bases: read_u32(bytes, 32),
            line_width: read_u32(bytes, 36),
        })
    }

    pub(crate) fn locus_occurrence(
        &self,
        document: u32,
        locus: u64,
        orientation: bool,
        k: u8,
    ) -> Result<crate::jidx_reader::SeedOccurrence, OwnerReaderError> {
        let record = self.document_record(document)?;
        let mut low = record.contig_start;
        let mut high = low
            .checked_add(record.contig_count)
            .ok_or(OwnerReaderError::Invalid("contig range"))?;
        while low < high {
            let middle = low + (high - low) / 2;
            if read_u64(self.contig_bytes(middle)?, 24) <= locus {
                low = middle + 1;
            } else {
                high = middle;
            }
        }
        if low == record.contig_start {
            return Err(OwnerReaderError::Invalid("locus contig"));
        }
        let bytes = self.contig_bytes(low - 1)?;
        let contig_id = read_u32(bytes, 12);
        let offset = locus
            .checked_sub(read_u64(bytes, 24))
            .ok_or(OwnerReaderError::Invalid("locus position"))?;
        let line_bases = u64::from(read_u32(bytes, 32));
        let line_width = u64::from(read_u32(bytes, 36));
        if read_u32(bytes, 0) != document
            || line_bases == 0
            || line_width < line_bases
            || offset % line_width >= line_bases
            || contig_id < record.original_contig_start
            || contig_id - record.original_contig_start >= record.original_contig_count
        {
            return Err(OwnerReaderError::Invalid("locus metadata"));
        }
        let position = (offset / line_width)
            .checked_mul(line_bases)
            .and_then(|position| position.checked_add(offset % line_width))
            .ok_or(OwnerReaderError::Invalid("locus position"))?;
        if position
            .checked_add(u64::from(k))
            .is_none_or(|end| end > read_u64(bytes, 16))
        {
            return Err(OwnerReaderError::Invalid("occurrence position"));
        }
        Ok(crate::jidx_reader::SeedOccurrence {
            contig_id,
            position,
            canonical_orientation: orientation,
        })
    }

    fn string(&self, offset: u32, length: u32) -> Result<&str, OwnerReaderError> {
        self.observer.record_metadata(0, 0, u64::from(length), 0);
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
        if kind == OwnerSection::Gzi {
            self.observer.record_metadata(0, 0, 0, length);
        }
        self.checked_bytes(start, end)
    }

    fn checked_bytes(&self, start: u64, end: u64) -> Result<&[u8], OwnerReaderError> {
        self.verify_unchanged()?;
        if start > end || end > self.mmap.len() as u64 {
            return Err(OwnerReaderError::Invalid("file range"));
        }
        self.observer
            .add(owner_observer::REQUESTED_BYTES, end - start);
        self.observer.add(owner_observer::READ_CALLS, 1);
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
        self.observer.record_integrity(1, 0, 0, 0, 0);
        if self.file_identity.is_some() && word.load(Ordering::Acquire) & mask != 0 {
            return Ok(());
        }
        let start = page * OWNER_PAGE_SIZE;
        let checksum_section = self.header.section(OwnerSection::PageChecksums);
        let expected_start = checksum_section.offset + index * 32;
        let expected = &self.mmap[expected_start as usize..expected_start as usize + 32];
        self.observer.record_integrity(0, 1, 0, OWNER_PAGE_SIZE, 0);
        if sha256(&self.mmap[start as usize..start as usize + OWNER_PAGE_SIZE as usize]) != expected
        {
            return Err(OwnerReaderError::ChecksumMismatch);
        }
        self.verify_checksum_chain(index)?;
        self.verify_unchanged()?;
        if self.file_identity.is_some() {
            if word.fetch_or(mask, Ordering::Release) & mask == 0 {
                self.observer.record_integrity(0, 0, 1, 0, 0);
            }
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
            self.observer.record_integrity(0, 0, 0, 0, 1);
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
        self.observer
            .add(owner_observer::CHECKSUM_ROOT_PAGES_HASHED, 1);
        if sha256(&self.mmap[start as usize..start as usize + OWNER_PAGE_SIZE as usize])
            != self.header.checksum_root_sha256
        {
            return Err(OwnerReaderError::ChecksumMismatch);
        }
        Ok(())
    }

    pub(crate) fn verify_unchanged(&self) -> io::Result<()> {
        if self.file_identity.is_some() {
            self.observer.add(owner_observer::FILE_IDENTITY_CHECKS, 1);
        }
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
    let invalid = || OwnerReaderError::Invalid("decoded hot postings");
    if hot.len() < 40 {
        return Err(invalid());
    }
    let members = usize::try_from(read_u64(hot, 16)).map_err(|_| invalid())?;
    let anchors = read_u32(hot, 32) as usize;
    let keys = u16::from_le_bytes(hot[10..12].try_into().expect("key count")) as usize;
    members
        .checked_mul(std::mem::size_of::<OwnerHotMember>())
        .and_then(|bytes| {
            anchors
                .checked_mul(std::mem::size_of::<OwnerAnchor>())
                .and_then(|anchors| bytes.checked_add(anchors))
        })
        .and_then(|bytes| {
            keys.checked_mul(std::mem::size_of::<OwnerHotKey>())
                .and_then(|keys| bytes.checked_add(keys))
        })
        .and_then(|bytes| bytes.checked_add(std::mem::size_of::<OwnerHotBlock>()))
        .ok_or_else(invalid)
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
