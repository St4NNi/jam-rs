//! Checked selective reader for the self-contained JMA format 1 object.

use crate::jma::contigs::decode_contigs;
use crate::jma::format::{
    HEADER_SIZE, SECTION_CONTIGS, SECTION_FLAG_COMPRESSED, SECTION_METADATA, SECTION_SEEDS,
    SECTION_SEQUENCE_DIRECTORY, SECTION_SEQUENCE_PAYLOAD, SectionDescriptor, checked_end,
    checked_usize, checksum, decode_metadata, get_u32, get_u64,
};
use crate::jma::header::{
    parse_header, parse_section_directory, unique_section, validate_known_sections,
};
use crate::jma::index::{
    SEED_DIRECTORY_HEADER_SIZE_U64, SEED_PAGE_RECORD_SIZE_U64, SEED_SCHEME_RECORD_SIZE_U64,
    SeedIndexDirectory, SeedSchemeIndex, decode_seed_index_directory_with_length, decode_seed_page,
    decode_shared_sequence_directory,
};
use crate::jma::{
    ArchiveHeader, ArchiveReader, ContigId, ContigMetadata, JmaError, JmaResult, SeedLevel,
    SeedOccurrence, SeedQuery, SequenceRange,
};
use crate::resource::{ByteRange, RangeBytes, RangeReader, ResourceMetrics};
use crate::sequence::{
    DEFAULT_MAX_DECODED_BLOCK_BYTES, SequenceBlockRecord, decode_block_range_bounded,
};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_SECTION_COUNT: u32 = 1 << 20;
const MAX_INDEXED_RANGE_BYTES: u64 = 256 * 1024 * 1024;

/// A validated JMA object. Opening reads the superblock, directories, metadata,
/// and contig table only. Seed pages and sequence payloads are fetched on use.
pub struct JmaReader<R> {
    resource: R,
    header: ArchiveHeader,
    metadata_fields: crate::jma::format::ArchiveMetadataFields,
    sections: Vec<SectionDescriptor>,
    contigs: Vec<ContigMetadata>,
    sequence_index: Vec<SequenceBlockRecord>,
    seed_index: SeedIndexDirectory,
    metrics: ReaderMetrics,
}

#[derive(Default)]
struct ReaderMetrics {
    metadata_bytes_read: AtomicU64,
    seed_pages_read: AtomicU64,
    sequence_blocks_read: AtomicU64,
    seed_bytes_read: AtomicU64,
    sequence_bytes_read: AtomicU64,
    decoded_sequence_bases: AtomicU64,
    coalesced_range_requests: AtomicU64,
    decoded_bytes: AtomicU64,
}

struct SequenceBlockBytes<'a> {
    bytes: RangeBytes<'a>,
    payload: std::ops::Range<usize>,
    ambiguity: std::ops::Range<usize>,
}

impl SequenceBlockBytes<'_> {
    fn parts(&self) -> JmaResult<(&[u8], &[u8])> {
        let bytes = self.bytes.as_ref();
        let payload = bytes.get(self.payload.clone()).ok_or_else(|| {
            JmaError::CorruptSection("coalesced sequence payload is truncated".to_string())
        })?;
        let ambiguity = bytes.get(self.ambiguity.clone()).ok_or_else(|| {
            JmaError::CorruptSection("coalesced ambiguity payload is truncated".to_string())
        })?;
        Ok((payload, ambiguity))
    }

    fn len(&self) -> usize {
        self.bytes.as_ref().len()
    }
}

impl ReaderMetrics {
    fn with_metadata_bytes(bytes: u64) -> Self {
        Self {
            metadata_bytes_read: AtomicU64::new(bytes),
            ..Self::default()
        }
    }

    fn add_seed_page(&self, bytes: usize) {
        self.seed_pages_read.fetch_add(1, Ordering::Relaxed);
        self.seed_bytes_read
            .fetch_add(u64::try_from(bytes).unwrap_or(u64::MAX), Ordering::Relaxed);
        self.coalesced_range_requests
            .fetch_add(1, Ordering::Relaxed);
        self.decoded_bytes
            .fetch_add(u64::try_from(bytes).unwrap_or(u64::MAX), Ordering::Relaxed);
    }

    fn add_sequence_block(&self, bytes: usize, decoded_bases: usize) {
        self.sequence_blocks_read.fetch_add(1, Ordering::Relaxed);
        self.sequence_bytes_read
            .fetch_add(u64::try_from(bytes).unwrap_or(u64::MAX), Ordering::Relaxed);
        self.decoded_sequence_bases.fetch_add(
            u64::try_from(decoded_bases).unwrap_or(u64::MAX),
            Ordering::Relaxed,
        );
        self.coalesced_range_requests
            .fetch_add(1, Ordering::Relaxed);
        self.decoded_bytes
            .fetch_add(u64::try_from(bytes).unwrap_or(u64::MAX), Ordering::Relaxed);
    }

    fn snapshot(&self) -> ResourceMetrics {
        ResourceMetrics {
            decoded_bytes: self.decoded_bytes.load(Ordering::Relaxed),
            seed_buckets_read: self.seed_pages_read.load(Ordering::Relaxed),
            sequence_blocks_read: self.sequence_blocks_read.load(Ordering::Relaxed),
            ..ResourceMetrics::default()
        }
    }

    fn archive_snapshot(&self, resource: ResourceMetrics) -> crate::archive::ArchiveMetrics {
        let seed_bytes_read = self.seed_bytes_read.load(Ordering::Relaxed);
        let sequence_bytes_read = self.sequence_bytes_read.load(Ordering::Relaxed);
        crate::archive::ArchiveMetrics {
            resource,
            mapped_bytes: resource.mapped_bytes,
            resident_bytes: resource.resident_bytes,
            metadata_bytes_read: self.metadata_bytes_read.load(Ordering::Relaxed),
            seed_bytes_read,
            sequence_bytes_read,
            decoded_sequence_bases: self.decoded_sequence_bases.load(Ordering::Relaxed),
            coalesced_range_requests: self.coalesced_range_requests.load(Ordering::Relaxed),
        }
    }
}

impl<R: RangeReader> JmaReader<R> {
    /// Opens a local or remote range resource without reading complete seed or
    /// sequence sections.
    pub fn from_resource(resource: R) -> JmaResult<Self> {
        let resource_metadata = resource.metadata()?;
        let header_size = u64::try_from(HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
        if resource_metadata.size < header_size {
            return Err(JmaError::CorruptSection(
                "JMA object is shorter than its superblock".to_string(),
            ));
        }
        let header_bytes = read_exact(&resource, ByteRange::new(0, header_size)?)?;
        let parsed = parse_header(header_bytes.as_ref())?;
        if parsed.object_size != resource_metadata.size {
            return Err(JmaError::CorruptSection(format!(
                "declared object size {} does not match resource size {}",
                parsed.object_size, resource_metadata.size
            )));
        }
        if parsed.section_count > MAX_SECTION_COUNT {
            return Err(JmaError::CorruptSection(format!(
                "section count {} exceeds limit {MAX_SECTION_COUNT}",
                parsed.section_count
            )));
        }
        let directory = read_exact(
            &resource,
            ByteRange::new(
                parsed.section_directory_offset,
                parsed.section_directory_length,
            )?,
        )?;
        if checksum(directory.as_ref()) != parsed.section_directory_checksum {
            return Err(JmaError::ChecksumMismatch("section directory".to_string()));
        }
        let sections = parse_section_directory(
            directory.as_ref(),
            parsed.section_count,
            parsed.section_directory_offset,
            resource_metadata.size,
        )?;
        validate_known_sections(&sections)?;
        let metadata_descriptor = unique_section(&sections, SECTION_METADATA)?;
        let contig_descriptor = unique_section(&sections, SECTION_CONTIGS)?;
        let seed_descriptor = unique_section(&sections, SECTION_SEEDS)?;
        let sequence_directory_descriptor = unique_section(&sections, SECTION_SEQUENCE_DIRECTORY)?;
        let sequence_payload_descriptor = unique_section(&sections, SECTION_SEQUENCE_PAYLOAD)?;

        let metadata_payload = read_section(&resource, metadata_descriptor)?;
        let mut metadata_fields = decode_metadata(metadata_payload.as_ref())?;
        if let Some(metadata_checksum) = metadata_fields.archive_sha256
            && metadata_checksum != parsed.archive_sha256
        {
            return Err(JmaError::CorruptSection(
                "metadata archive checksum does not match superblock".to_string(),
            ));
        }
        // The fixed superblock owns the object checksum so computing it does
        // not create a recursive metadata/section checksum dependency.
        metadata_fields.archive_sha256 =
            (parsed.archive_sha256 != [0; 32]).then_some(parsed.archive_sha256);
        if metadata_fields.source_assembly_sha256 != parsed.archive.source_sha256 {
            return Err(JmaError::CorruptSection(
                "metadata source checksum does not match superblock".to_string(),
            ));
        }
        if metadata_fields.format_version != parsed.archive.format_version {
            return Err(JmaError::UnsupportedVersion(metadata_fields.format_version));
        }
        if metadata_fields.format_identifier != "JMA"
            || metadata_fields.hash_algorithm != "jamhash_u64_v1"
        {
            return Err(JmaError::CorruptSection(
                "metadata format or hash algorithm identifier is unsupported".to_string(),
            ));
        }
        let contig_payload = read_section(&resource, contig_descriptor)?;
        let contigs = decode_contigs(contig_payload.as_ref(), parsed.archive.contig_count)?;
        let total_bases = contigs.iter().try_fold(0u64, |total, contig| {
            total
                .checked_add(contig.length)
                .ok_or(JmaError::OffsetOverflow)
        })?;
        if total_bases != parsed.archive.total_bases {
            return Err(JmaError::CorruptSection(format!(
                "contig bases {total_bases} do not match archive total {}",
                parsed.archive.total_bases
            )));
        }
        let seed_directory_prefix = read_seed_directory_prefix(&resource, seed_descriptor)?;
        let seed_index = decode_seed_index_directory_with_length(
            seed_directory_prefix.as_ref(),
            seed_descriptor.length,
        )?;
        let sequence_directory_payload = read_section(&resource, sequence_directory_descriptor)?;
        let sequence_index = decode_shared_sequence_directory(
            sequence_directory_payload.as_ref(),
            sequence_payload_descriptor.length,
            sequence_payload_descriptor.offset,
            &contigs,
            resource_metadata.size,
        )?;
        validate_sequence_coverage(&sequence_index, &contigs)?;
        let metadata_bytes_read = [
            header_bytes.as_ref().len(),
            directory.as_ref().len(),
            metadata_payload.as_ref().len(),
            contig_payload.as_ref().len(),
            seed_directory_prefix.as_ref().len(),
            sequence_directory_payload.as_ref().len(),
        ]
        .into_iter()
        .try_fold(0u64, |total, length| {
            let length = u64::try_from(length).map_err(|_| JmaError::OffsetOverflow)?;
            total.checked_add(length).ok_or(JmaError::OffsetOverflow)
        })?;

        let mut header = parsed.archive;
        header.min_entropy = metadata_fields.min_entropy;
        header.seed_levels = seed_index
            .schemes
            .iter()
            .map(|scheme| SeedLevel {
                k: u8::try_from(scheme.descriptor.span).unwrap_or(u8::MAX),
                scale: u64::from(scheme.descriptor.density_parameter),
            })
            .collect();
        Ok(Self {
            resource,
            header,
            metadata_fields,
            sections,
            contigs,
            sequence_index,
            seed_index,
            metrics: ReaderMetrics::with_metadata_bytes(metadata_bytes_read),
        })
    }

    pub fn open(resource: R) -> JmaResult<Self> {
        Self::from_resource(resource)
    }

    #[must_use]
    pub fn sections(&self) -> &[SectionDescriptor] {
        &self.sections
    }

    #[must_use]
    pub fn metadata_fields(&self) -> &crate::jma::format::ArchiveMetadataFields {
        &self.metadata_fields
    }

    pub fn seed_schemes(&self) -> impl Iterator<Item = &crate::archive::SeedSchemeDescriptor> {
        self.seed_index
            .schemes
            .iter()
            .map(|scheme| &scheme.descriptor)
    }

    #[must_use]
    pub fn seed_index(&self) -> &SeedIndexDirectory {
        &self.seed_index
    }

    #[must_use]
    pub fn sequence_index(&self) -> &[SequenceBlockRecord] {
        &self.sequence_index
    }

    #[must_use]
    pub fn archive_metrics(&self) -> crate::archive::ArchiveMetrics {
        self.metrics.archive_snapshot(
            self.resource
                .metrics()
                .saturating_add(self.metrics.snapshot()),
        )
    }

    pub fn seed_occurrences_for_scheme(
        &self,
        scheme_id: crate::archive::SeedSchemeId,
        query: SeedQuery,
    ) -> JmaResult<Vec<SeedOccurrence>> {
        let scheme = self
            .seed_index
            .schemes
            .iter()
            .find(|scheme| scheme.descriptor.scheme_id == scheme_id.0)
            .ok_or_else(|| {
                JmaError::CorruptSection(format!("seed scheme {} is unavailable", scheme_id.0))
            })?;
        if scheme.descriptor.span != u16::from(query.k) {
            return Err(JmaError::CorruptSection(
                "seed query span does not match its scheme".to_string(),
            ));
        }
        self.lookup_scheme_query(scheme, query)
    }

    /// Looks up exact keys in one scheme while reading every selected page at
    /// most once. Results retain caller order, including duplicate keys.
    pub fn seed_occurrences_for_scheme_batch(
        &self,
        scheme_id: crate::archive::SeedSchemeId,
        queries: &[SeedQuery],
    ) -> JmaResult<Vec<Vec<SeedOccurrence>>> {
        let scheme = self
            .seed_index
            .schemes
            .iter()
            .find(|scheme| scheme.descriptor.scheme_id == scheme_id.0)
            .ok_or_else(|| {
                JmaError::CorruptSection(format!("seed scheme {} is unavailable", scheme_id.0))
            })?;
        let k = u8::try_from(scheme.descriptor.span).map_err(|_| JmaError::OffsetOverflow)?;
        let mut page_queries = BTreeMap::<usize, Vec<usize>>::new();
        for (query_index, query) in queries.iter().copied().enumerate() {
            if query.k != k {
                return Err(JmaError::CorruptSection(
                    "seed query span does not match its scheme".to_string(),
                ));
            }
            if query.hash == 0 {
                continue;
            }
            let prefix =
                crate::jma::format::hash_prefix(query.hash, scheme.descriptor.bucket_bits)?;
            for (page_index, _page) in scheme.pages.iter().enumerate().filter(|(_, page)| {
                page.hash_prefix == prefix
                    && page.first_hash <= query.hash
                    && page.last_hash >= query.hash
            }) {
                page_queries
                    .entry(page_index)
                    .or_default()
                    .push(query_index);
            }
        }

        let seed_descriptor = unique_section(&self.sections, SECTION_SEEDS)?;
        let mut matches = vec![Vec::new(); queries.len()];
        for (page_index, query_indices) in page_queries {
            let page = scheme.pages.get(page_index).ok_or_else(|| {
                JmaError::CorruptSection("selected seed page is unavailable".to_string())
            })?;
            let records = self.read_seed_page(seed_descriptor, scheme, page)?;
            let mut indices_by_key = BTreeMap::<(u64, u64), Vec<usize>>::new();
            for query_index in query_indices {
                let query = queries[query_index];
                indices_by_key
                    .entry((query.hash, query.canonical_kmer))
                    .or_default()
                    .push(query_index);
            }
            for record in records {
                let Some(query_indices) =
                    indices_by_key.get(&(record.query.hash, record.query.canonical_kmer))
                else {
                    continue;
                };
                self.validate_seed_occurrence(record.occurrence, k)?;
                for &query_index in query_indices {
                    matches[query_index].push(record.occurrence);
                }
            }
        }
        for occurrences in &mut matches {
            occurrences.sort_by_key(|occurrence| {
                (
                    occurrence.contig_id,
                    occurrence.position,
                    occurrence.reverse,
                )
            });
        }
        Ok(matches)
    }

    /// Looks up a seed at a density level. Page and occurrence ranges are
    /// read only after selecting by hash prefix; canonical k-mer equality is
    /// checked after the digest match.
    pub fn seed_occurrences_at_scale(
        &self,
        query: SeedQuery,
        scale: u64,
    ) -> JmaResult<Vec<SeedOccurrence>> {
        if query.hash == 0 {
            return Ok(Vec::new());
        }
        let scheme = self
            .seed_index
            .schemes
            .iter()
            .find(|scheme| {
                scheme.descriptor.span == u16::from(query.k)
                    && u64::from(scheme.descriptor.density_parameter) == scale
            })
            .ok_or_else(|| {
                JmaError::CorruptSection(format!(
                    "seed scale {scale} is unavailable for k={}",
                    query.k
                ))
            })?;
        self.lookup_scheme_query(scheme, query)
    }

    /// Reads and reverse-complements only the requested contig range. The
    /// range is expressed in forward contig coordinates.
    pub fn read_sequence_reverse_complement(
        &self,
        contig_id: ContigId,
        range: SequenceRange,
    ) -> JmaResult<Vec<u8>> {
        let mut bases = self.read_sequence(contig_id, range)?;
        bases.reverse();
        for base in &mut bases {
            *base = complement_base(*base)?;
        }
        Ok(bases)
    }

    fn lookup_scheme_query(
        &self,
        scheme: &SeedSchemeIndex,
        query: SeedQuery,
    ) -> JmaResult<Vec<SeedOccurrence>> {
        let prefix = crate::jma::format::hash_prefix(query.hash, scheme.descriptor.bucket_bits)?;
        let seed_descriptor = unique_section(&self.sections, SECTION_SEEDS)?;
        let mut matches = Vec::new();
        for page in scheme.pages.iter().filter(|page| {
            page.hash_prefix == prefix
                && page.first_hash <= query.hash
                && page.last_hash >= query.hash
        }) {
            let records = self.read_seed_page(seed_descriptor, scheme, page)?;
            for record in records.into_iter().filter(|record| {
                record.query.hash == query.hash
                    && record.query.canonical_kmer == query.canonical_kmer
            }) {
                self.validate_seed_occurrence(record.occurrence, query.k)?;
                matches.push(record.occurrence);
            }
        }
        matches.sort_by_key(|occurrence| {
            (
                occurrence.contig_id,
                occurrence.position,
                occurrence.reverse,
            )
        });
        Ok(matches)
    }

    fn read_seed_page(
        &self,
        seed_descriptor: &SectionDescriptor,
        scheme: &SeedSchemeIndex,
        page: &crate::jma::index::SeedPageIndex,
    ) -> JmaResult<Vec<crate::jma::writer::SeedRecord>> {
        let key_offset = seed_descriptor
            .offset
            .checked_add(page.key_offset)
            .ok_or(JmaError::OffsetOverflow)?;
        let occurrence_offset = seed_descriptor
            .offset
            .checked_add(page.occurrence_offset)
            .ok_or(JmaError::OffsetOverflow)?;
        let range_start = key_offset.min(occurrence_offset);
        let range_end = checked_end(key_offset, page.key_length)?
            .max(checked_end(occurrence_offset, page.occurrence_length)?);
        let combined = self.read_indexed_range(
            range_start,
            range_end
                .checked_sub(range_start)
                .ok_or(JmaError::OffsetOverflow)?,
        )?;
        let combined = combined.as_ref();
        let key_start = usize::try_from(
            key_offset
                .checked_sub(range_start)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .map_err(|_| JmaError::OffsetOverflow)?;
        let key_end = key_start
            .checked_add(usize::try_from(page.key_length).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        let occurrence_start = usize::try_from(
            occurrence_offset
                .checked_sub(range_start)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .map_err(|_| JmaError::OffsetOverflow)?;
        let occurrence_end = occurrence_start
            .checked_add(
                usize::try_from(page.occurrence_length).map_err(|_| JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        let key_bytes = combined.get(key_start..key_end).ok_or_else(|| {
            JmaError::CorruptSection("coalesced seed key range is truncated".to_string())
        })?;
        let occurrence_bytes = combined
            .get(occurrence_start..occurrence_end)
            .ok_or_else(|| {
                JmaError::CorruptSection("coalesced seed occurrence range is truncated".to_string())
            })?;
        let records = decode_seed_page(
            key_bytes,
            occurrence_bytes,
            page,
            u8::try_from(scheme.descriptor.span).map_err(|_| JmaError::OffsetOverflow)?,
        )?;
        self.metrics
            .add_seed_page(key_bytes.len().saturating_add(occurrence_bytes.len()));
        Ok(records)
    }

    fn validate_seed_occurrence(&self, occurrence: SeedOccurrence, k: u8) -> JmaResult<()> {
        let contig = self
            .contigs
            .iter()
            .find(|contig| contig.id == occurrence.contig_id)
            .ok_or(JmaError::UnknownContig(occurrence.contig_id))?;
        let end = occurrence
            .position
            .checked_add(u64::from(k))
            .ok_or(JmaError::OffsetOverflow)?;
        if end > contig.length {
            return Err(JmaError::CorruptSection(
                "seed occurrence exceeds contig length".to_string(),
            ));
        }
        Ok(())
    }

    fn read_sequence_block_parts<'a>(
        &'a self,
        record: &SequenceBlockRecord,
    ) -> JmaResult<SequenceBlockBytes<'a>> {
        let payload_end = record.payload_end().ok_or(JmaError::OffsetOverflow)?;
        let ambiguity_end = record.ambiguity_end().ok_or(JmaError::OffsetOverflow)?;
        let range_start = record.payload_offset.min(record.ambiguity_offset);
        let range_end = payload_end.max(ambiguity_end);
        let range_length = range_end
            .checked_sub(range_start)
            .ok_or(JmaError::OffsetOverflow)?;
        let bytes = self.read_indexed_range(range_start, range_length)?;
        let payload_start = usize::try_from(
            record
                .payload_offset
                .checked_sub(range_start)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .map_err(|_| JmaError::OffsetOverflow)?;
        let payload_end = payload_start
            .checked_add(
                usize::try_from(record.stored_length).map_err(|_| JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        let ambiguity_start = usize::try_from(
            record
                .ambiguity_offset
                .checked_sub(range_start)
                .ok_or(JmaError::OffsetOverflow)?,
        )
        .map_err(|_| JmaError::OffsetOverflow)?;
        let ambiguity_end = ambiguity_start
            .checked_add(
                usize::try_from(record.ambiguity_length).map_err(|_| JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        Ok(SequenceBlockBytes {
            bytes,
            payload: payload_start..payload_end,
            ambiguity: ambiguity_start..ambiguity_end,
        })
    }

    fn read_indexed_range(&self, offset: u64, length: u64) -> JmaResult<RangeBytes<'_>> {
        if length > MAX_INDEXED_RANGE_BYTES {
            return Err(JmaError::CorruptSection(
                "indexed JMA range exceeds the read limit".to_string(),
            ));
        }
        read_exact(&self.resource, ByteRange::new(offset, length)?)
    }
}

impl<R: RangeReader> ArchiveReader for JmaReader<R> {
    fn resource(&self) -> &dyn RangeReader {
        &self.resource
    }

    fn header(&self) -> &ArchiveHeader {
        &self.header
    }

    fn contigs(&self) -> &[ContigMetadata] {
        &self.contigs
    }

    fn read_sequence(&self, contig_id: ContigId, range: SequenceRange) -> JmaResult<Vec<u8>> {
        let contig = self
            .contigs
            .iter()
            .find(|contig| contig.id == contig_id)
            .ok_or(JmaError::UnknownContig(contig_id))?;
        if range.start > range.end || range.end > contig.length {
            return Err(JmaError::InvalidSequenceRange {
                start: range.start,
                end: range.end,
            });
        }
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let mut output =
            Vec::with_capacity(usize::try_from(range.len()).map_err(|_| JmaError::OffsetOverflow)?);
        for entry in self.sequence_index.iter().filter(|entry| {
            entry.contig_id == contig_id
                && entry.base_start < range.end
                && entry.base_end().is_some_and(|end| end > range.start)
        }) {
            let request_start = range.start.max(entry.base_start);
            let request_end = range
                .end
                .min(entry.base_end().ok_or(JmaError::OffsetOverflow)?);
            let block = self.read_sequence_block_parts(entry)?;
            let (payload, ambiguity_payload) = block.parts()?;
            let decoded = decode_block_range_bounded(
                entry,
                payload,
                ambiguity_payload,
                request_start..request_end,
                DEFAULT_MAX_DECODED_BLOCK_BYTES,
            )
            .map_err(|error| JmaError::CorruptSection(error.to_string()))?;
            self.metrics.add_sequence_block(block.len(), decoded.len());
            output.extend_from_slice(&decoded);
        }
        if output.len() != usize::try_from(range.len()).map_err(|_| JmaError::OffsetOverflow)? {
            return Err(JmaError::CorruptSection(
                "sequence range is not covered by indexed blocks".to_string(),
            ));
        }
        Ok(output)
    }

    fn seed_occurrences(&self, query: SeedQuery) -> JmaResult<Vec<SeedOccurrence>> {
        if query.hash == 0 {
            return Ok(Vec::new());
        }
        let scheme = self
            .seed_index
            .schemes
            .iter()
            .find(|scheme| scheme.descriptor.span == u16::from(query.k))
            .ok_or_else(|| {
                JmaError::CorruptSection(format!("seed level k={} is unavailable", query.k))
            })?;
        self.lookup_scheme_query(scheme, query)
    }

    fn metrics(&self) -> ResourceMetrics {
        self.resource
            .metrics()
            .saturating_add(self.metrics.snapshot())
    }
}

fn read_exact<'a, R: RangeReader>(resource: &'a R, range: ByteRange) -> JmaResult<RangeBytes<'a>> {
    let expected = checked_usize(range.length, "resource read length")?;
    let bytes = resource.read_range_bytes(range)?;
    if bytes.as_ref().len() != expected {
        return Err(JmaError::CorruptSection(format!(
            "resource returned {} bytes for requested {expected}",
            bytes.as_ref().len()
        )));
    }
    Ok(bytes)
}

fn read_seed_directory_prefix<'a, R: RangeReader>(
    resource: &'a R,
    descriptor: &SectionDescriptor,
) -> JmaResult<RangeBytes<'a>> {
    if descriptor.length < SEED_DIRECTORY_HEADER_SIZE_U64 {
        return Err(JmaError::CorruptSection(
            "seed directory section is shorter than its header".to_string(),
        ));
    }
    let header = read_exact(
        resource,
        ByteRange::new(descriptor.offset, SEED_DIRECTORY_HEADER_SIZE_U64)?,
    )?;
    let header = header.as_ref();
    let scheme_count = get_u32(header, 4)?;
    let page_count = get_u32(header, 8)?;
    let scheme_offset = get_u64(header, 16)?;
    let page_offset = get_u64(header, 24)?;
    let data_offset = get_u64(header, 32)?;
    let scheme_end = checked_end(
        scheme_offset,
        u64::from(scheme_count)
            .checked_mul(SEED_SCHEME_RECORD_SIZE_U64)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    let page_end = checked_end(
        page_offset,
        u64::from(page_count)
            .checked_mul(SEED_PAGE_RECORD_SIZE_U64)
            .ok_or(JmaError::OffsetOverflow)?,
    )?;
    if scheme_offset < SEED_DIRECTORY_HEADER_SIZE_U64
        || page_offset < scheme_end
        || data_offset < page_end
        || data_offset > descriptor.length
        || data_offset > MAX_INDEXED_RANGE_BYTES
    {
        return Err(JmaError::CorruptSection(
            "seed directory prefix offsets are invalid".to_string(),
        ));
    }
    let prefix = read_exact(resource, ByteRange::new(descriptor.offset, data_offset)?)?;
    if checksum(prefix.as_ref()) != descriptor.checksum {
        return Err(JmaError::ChecksumMismatch(
            "seed directory prefix".to_string(),
        ));
    }
    Ok(prefix)
}

fn read_section<'a, R: RangeReader>(
    resource: &'a R,
    descriptor: &SectionDescriptor,
) -> JmaResult<RangeBytes<'a>> {
    checked_end(descriptor.offset, descriptor.length)?;
    let bytes = read_exact(
        resource,
        ByteRange::new(descriptor.offset, descriptor.length)?,
    )?;
    if checksum(bytes.as_ref()) != descriptor.checksum {
        return Err(JmaError::ChecksumMismatch(format!(
            "section kind {}",
            descriptor.kind
        )));
    }
    if descriptor.flags & SECTION_FLAG_COMPRESSED != 0 {
        return Err(JmaError::CorruptSection(format!(
            "compressed section kind {} is not available in this reader",
            descriptor.kind
        )));
    }
    if descriptor.uncompressed_length != descriptor.length {
        return Err(JmaError::CorruptSection(format!(
            "uncompressed section kind {} has an invalid decoded length",
            descriptor.kind
        )));
    }
    Ok(bytes)
}

fn validate_sequence_coverage(
    blocks: &[SequenceBlockRecord],
    contigs: &[ContigMetadata],
) -> JmaResult<()> {
    let mut cursor = 0usize;
    for contig in contigs {
        let mut expected = 0u64;
        while cursor < blocks.len() && blocks[cursor].contig_id == contig.id {
            if blocks[cursor].base_start != expected {
                return Err(JmaError::CorruptSection(format!(
                    "sequence blocks leave a gap or overlap in contig {} at {}",
                    contig.id, expected
                )));
            }
            expected = blocks[cursor].base_end().ok_or(JmaError::OffsetOverflow)?;
            cursor += 1;
        }
        if expected != contig.length {
            return Err(JmaError::CorruptSection(format!(
                "sequence blocks cover {expected} of contig {}, expected {}",
                contig.id, contig.length
            )));
        }
    }
    if cursor != blocks.len() {
        return Err(JmaError::CorruptSection(
            "sequence directory contains an unknown contig or unsorted block".to_string(),
        ));
    }
    Ok(())
}

fn complement_base(base: u8) -> JmaResult<u8> {
    let value = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' | b'-' => base.to_ascii_uppercase(),
        other => {
            return Err(JmaError::CorruptSection(format!(
                "unsupported ambiguity base 0x{other:02x}"
            )));
        }
    };
    Ok(value)
}
