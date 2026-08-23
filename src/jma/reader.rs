//! Checked random-access reader for JMA v1 resources.

use crate::jma::contigs::decode_contigs;
use crate::jma::format::{
    HEADER_SIZE, SECTION_CONTIGS, SECTION_SEEDS_K21, SECTION_SEEDS_K31, SECTION_SEQUENCES,
    SectionDescriptor, checked_end, checked_usize, checksum,
};
use crate::jma::header::{parse_header, parse_section_directory};
use crate::jma::index::{
    JmaSidecarIndex, SeedLevelIndex, SequenceBlockIndex, bucket_for, decode_index, level_for,
    validate_against_archive,
};
use crate::jma::sequence::{
    SequenceBlock, decode_range, decode_sequence_block, decode_sequence_blocks, validate_blocks,
};
use crate::jma::writer::{SeedRecord, SeedSection, decode_seed_records, decode_seed_section};
use crate::jma::{
    ArchiveHeader, ArchiveReader, ContigId, ContigMetadata, JmaError, JmaResult, SeedOccurrence,
    SeedQuery, SequenceRange,
};
use crate::resource::{ByteRange, RangeReader, ResourceMetrics};
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_SECTION_COUNT: u32 = 1 << 20;

/// A fully validated JMA archive backed by a range-readable resource.
pub struct JmaReader<R> {
    resource: R,
    index_resource: Option<R>,
    header: ArchiveHeader,
    sections: Vec<SectionDescriptor>,
    contigs: Vec<ContigMetadata>,
    sequence_blocks: Vec<SequenceBlock>,
    seed_sections: Vec<SeedSection>,
    sidecar: Option<JmaSidecarIndex>,
    sequence_index: Vec<SequenceBlockIndex>,
    seed_levels_index: Vec<SeedLevelIndex>,
    metrics: ReaderMetrics,
}

#[derive(Default)]
struct ReaderMetrics {
    seed_buckets_read: AtomicU64,
    sequence_blocks_read: AtomicU64,
    decoded_bytes: AtomicU64,
}

impl ReaderMetrics {
    fn add_seed_bucket(&self) {
        self.seed_buckets_read.fetch_add(1, Ordering::Relaxed);
    }

    fn add_sequence_block(&self) {
        self.sequence_blocks_read.fetch_add(1, Ordering::Relaxed);
    }

    fn add_decoded_bytes(&self, bytes: usize) {
        self.decoded_bytes
            .fetch_add(u64::try_from(bytes).unwrap_or(u64::MAX), Ordering::Relaxed);
    }

    fn snapshot(&self) -> ResourceMetrics {
        ResourceMetrics {
            decoded_bytes: self.decoded_bytes.load(Ordering::Relaxed),
            seed_buckets_read: self.seed_buckets_read.load(Ordering::Relaxed),
            sequence_blocks_read: self.sequence_blocks_read.load(Ordering::Relaxed),
            ..ResourceMetrics::default()
        }
    }
}

impl<R: RangeReader> JmaReader<R> {
    /// Opens and validates the header, directory, required sections, and
    /// checksums. Sequence and seed payloads are decoded once and then reused
    /// for subsequent range and lookup operations.
    pub fn from_resource(resource: R) -> JmaResult<Self> {
        let metadata = resource.metadata()?;
        let header_bytes = read_exact(&resource, ByteRange::new(0, HEADER_SIZE as u64)?)?;
        let parsed = parse_header(&header_bytes)?;
        if parsed.section_count > MAX_SECTION_COUNT {
            return Err(JmaError::CorruptSection(format!(
                "section count {} exceeds the JMA v1 limit {MAX_SECTION_COUNT}",
                parsed.section_count
            )));
        }
        let directory_end = checked_end(
            parsed.section_directory_offset,
            parsed.section_directory_length,
        )?;
        if directory_end > metadata.size {
            return Err(JmaError::CorruptSection(format!(
                "section directory ends at {directory_end}, resource has {} bytes",
                metadata.size
            )));
        }
        let directory_length =
            checked_usize(parsed.section_directory_length, "section directory length")?;
        let directory = read_exact(
            &resource,
            ByteRange::new(
                parsed.section_directory_offset,
                parsed.section_directory_length,
            )?,
        )?;
        if directory.len() != directory_length {
            return Err(JmaError::CorruptSection(
                "section directory read has an unexpected length".to_string(),
            ));
        }
        let sections = parse_section_directory(
            &directory,
            parsed.section_count,
            parsed.section_directory_offset,
            metadata.size,
        )?;
        validate_section_kinds(&sections)?;
        let contig_descriptor = find_unique_section(&sections, SECTION_CONTIGS)?;
        let sequence_descriptor = find_unique_section(&sections, SECTION_SEQUENCES)?;
        let contig_payload = read_section(&resource, contig_descriptor)?;
        let contigs = decode_contigs(&contig_payload, parsed.archive.contig_count)?;
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
        let sequence_payload = read_section(&resource, sequence_descriptor)?;
        let sequence_blocks = decode_sequence_blocks(&sequence_payload)?;
        validate_blocks(&sequence_blocks, &contigs)?;

        let mut seed_sections = Vec::new();
        for kind in [SECTION_SEEDS_K21, SECTION_SEEDS_K31] {
            if let Some(descriptor) = sections.iter().find(|entry| entry.kind == kind) {
                let payload = read_section(&resource, descriptor)?;
                let section = decode_seed_section(&payload)?;
                let expected_k = if kind == SECTION_SEEDS_K21 { 21 } else { 31 };
                if section.k != expected_k {
                    return Err(JmaError::CorruptSection(format!(
                        "seed section kind {kind} contains k={}",
                        section.k
                    )));
                }
                seed_sections.push(section);
            }
        }
        let actual_levels = seed_sections
            .iter()
            .flat_map(|section| section.levels.iter().map(|level| level.level))
            .collect::<Vec<_>>();
        if actual_levels != parsed.archive.seed_levels {
            return Err(JmaError::CorruptSection(
                "header seed levels do not match seed sections".to_string(),
            ));
        }

        Ok(Self {
            resource,
            index_resource: None,
            header: parsed.archive,
            sections,
            contigs,
            sequence_blocks,
            seed_sections,
            sidecar: None,
            sequence_index: Vec::new(),
            seed_levels_index: Vec::new(),
            metrics: ReaderMetrics::default(),
        })
    }

    /// Opens a checksum-bound sidecar without reading the sequence or seed
    /// payloads. The archive and sidecar use the same concrete range-reader
    /// type so local, HTTP, and S3 resources retain identical semantics.
    pub fn open_indexed(archive: R, index: R) -> JmaResult<Self> {
        let metadata = archive.metadata()?;
        let header_bytes = read_exact(&archive, ByteRange::new(0, HEADER_SIZE as u64)?)?;
        let parsed = parse_header(&header_bytes)?;
        if parsed.section_count > MAX_SECTION_COUNT {
            return Err(JmaError::CorruptSection(format!(
                "section count {} exceeds the JMA v1 limit {MAX_SECTION_COUNT}",
                parsed.section_count
            )));
        }
        let directory_end = checked_end(
            parsed.section_directory_offset,
            parsed.section_directory_length,
        )?;
        if directory_end > metadata.size {
            return Err(JmaError::CorruptSection(format!(
                "section directory ends at {directory_end}, resource has {} bytes",
                metadata.size
            )));
        }
        let directory = read_exact(
            &archive,
            ByteRange::new(
                parsed.section_directory_offset,
                parsed.section_directory_length,
            )?,
        )?;
        let sections = parse_section_directory(
            &directory,
            parsed.section_count,
            parsed.section_directory_offset,
            metadata.size,
        )?;
        validate_section_kinds(&sections)?;
        let contig_descriptor = find_unique_section(&sections, SECTION_CONTIGS)?;
        let contig_payload = read_section(&archive, contig_descriptor)?;
        let contigs = decode_contigs(&contig_payload, parsed.archive.contig_count)?;
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

        let index_metadata = index.metadata()?;
        if index_metadata.size > crate::jma::index::MAX_INDEX_BYTES {
            return Err(JmaError::CorruptSection(format!(
                "JMA sidecar is too large: {} bytes",
                index_metadata.size
            )));
        }
        let index_bytes = read_exact(&index, ByteRange::new(0, index_metadata.size)?)?;
        let sidecar = decode_index(&index_bytes)?;
        validate_against_archive(
            &sidecar,
            metadata.size,
            &header_bytes,
            &parsed.archive,
            &sections,
            &contigs,
        )?;

        Ok(Self {
            resource: archive,
            index_resource: Some(index),
            header: parsed.archive,
            sections,
            contigs,
            sequence_blocks: Vec::new(),
            seed_sections: Vec::new(),
            sequence_index: sidecar.sequence_blocks.clone(),
            seed_levels_index: sidecar.seed_levels.clone(),
            sidecar: Some(sidecar),
            metrics: ReaderMetrics::default(),
        })
    }

    /// Alias with an explicit name for callers that construct readers from a
    /// local or remote `RangeReader`.
    pub fn open(resource: R) -> JmaResult<Self> {
        Self::from_resource(resource)
    }

    #[must_use]
    pub fn sections(&self) -> &[SectionDescriptor] {
        &self.sections
    }

    /// Looks up a seed at an explicit density level. The exact hash and
    /// canonical packed k-mer are both compared; hash equality alone is not a
    /// valid occurrence match.
    pub fn seed_occurrences_at_scale(
        &self,
        query: SeedQuery,
        scale: u64,
    ) -> JmaResult<Vec<SeedOccurrence>> {
        if query.hash == 0 {
            return Ok(Vec::new());
        }
        if self.sidecar.is_some() {
            return self.indexed_seed_occurrences(query, scale);
        }
        let section = self
            .seed_sections
            .iter()
            .find(|section| section.k == query.k)
            .ok_or_else(|| {
                JmaError::CorruptSection(format!("seed level k={} is unavailable", query.k))
            })?;
        let level = section
            .levels
            .iter()
            .find(|level| level.level.scale == scale)
            .ok_or_else(|| {
                JmaError::CorruptSection(format!(
                    "seed scale {scale} is unavailable for k={}",
                    query.k
                ))
            })?;
        Ok(matching_occurrences(&level.records, query))
    }

    fn indexed_seed_occurrences(
        &self,
        query: SeedQuery,
        scale: u64,
    ) -> JmaResult<Vec<SeedOccurrence>> {
        let sidecar = self.sidecar.as_ref().ok_or_else(|| {
            JmaError::CorruptSection("indexed seed lookup has no sidecar".to_string())
        })?;
        let level = level_for(sidecar, query.k, scale).ok_or_else(|| {
            JmaError::CorruptSection(format!(
                "seed scale {scale} is unavailable for k={}",
                query.k
            ))
        })?;
        let Some(bucket) = bucket_for(sidecar, query.k, scale, query.hash) else {
            return Ok(Vec::new());
        };
        let bytes = self.read_indexed_range(bucket.offset, bucket.length)?;
        if checksum(&bytes) != bucket.checksum {
            return Err(JmaError::ChecksumMismatch(format!(
                "seed bucket k={} scale={} prefix={}",
                bucket.k, bucket.scale, bucket.hash_prefix
            )));
        }
        if bucket.length != bucket.record_count.saturating_mul(32)
            || bytes.len() != usize::try_from(bucket.length).unwrap_or(usize::MAX)
        {
            return Err(JmaError::CorruptSection(
                "seed bucket byte length does not match record count".to_string(),
            ));
        }
        let records = decode_seed_records(&bytes, query.k)?;
        if records.len() != usize::try_from(bucket.record_count).unwrap_or(usize::MAX) {
            return Err(JmaError::CorruptSection(
                "seed bucket record count does not match payload".to_string(),
            ));
        }
        let prefix = bucket.hash_prefix;
        if records.iter().any(|record| {
            record.query.k != query.k
                || (record.query.hash >> (64 - u32::from(crate::jma::index::HASH_PREFIX_BITS)))
                    as u16
                    != prefix
        }) {
            return Err(JmaError::CorruptSection(
                "seed bucket contains a record outside its hash prefix".to_string(),
            ));
        }
        self.metrics.add_seed_bucket();
        self.metrics.add_decoded_bytes(bytes.len());
        let _ = level;
        Ok(matching_occurrences(&records, query))
    }

    fn read_indexed_sequence(
        &self,
        contig_id: ContigId,
        range: SequenceRange,
    ) -> JmaResult<Vec<u8>> {
        let metadata = self
            .contigs
            .iter()
            .find(|contig| contig.id == contig_id)
            .ok_or(JmaError::UnknownContig(contig_id))?;
        if range.end > metadata.length {
            return Err(JmaError::InvalidSequenceRange {
                start: range.start,
                end: range.end,
            });
        }
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let mut blocks = Vec::new();
        for entry in self
            .sequence_index
            .iter()
            .filter(|entry| entry.contig_id == contig_id)
        {
            let block_end = entry
                .start
                .checked_add(entry.base_length)
                .ok_or(JmaError::OffsetOverflow)?;
            if entry.start >= range.end || block_end <= range.start {
                continue;
            }
            let bytes = self.read_indexed_range(entry.offset, entry.length)?;
            if checksum(&bytes) != entry.checksum {
                return Err(JmaError::ChecksumMismatch(format!(
                    "sequence block contig={} start={}",
                    entry.contig_id, entry.start
                )));
            }
            let block = decode_sequence_block(&bytes)?;
            if block.contig_id != entry.contig_id
                || block.start != entry.start
                || block.base_length != entry.base_length
            {
                return Err(JmaError::CorruptSection(
                    "sequence sidecar metadata does not match fetched block".to_string(),
                ));
            }
            self.metrics.add_sequence_block();
            self.metrics.add_decoded_bytes(bytes.len());
            blocks.push(block);
        }
        decode_range(&blocks, &self.contigs, contig_id, range)
    }

    fn read_indexed_range(&self, offset: u64, length: u64) -> JmaResult<Vec<u8>> {
        if length > crate::jma::index::MAX_INDEXED_RANGE_BYTES {
            return Err(JmaError::CorruptSection(
                "indexed JMA range exceeds the read limit".to_string(),
            ));
        }
        let range = ByteRange::new(offset, length)?;
        read_exact(&self.resource, range)
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
        if self.sidecar.is_some() {
            return self.read_indexed_sequence(contig_id, range);
        }
        decode_range(&self.sequence_blocks, &self.contigs, contig_id, range)
    }

    fn seed_occurrences(&self, query: SeedQuery) -> JmaResult<Vec<SeedOccurrence>> {
        if query.hash == 0 {
            return Ok(Vec::new());
        }
        if self.sidecar.is_some() {
            let level = self
                .seed_levels_index
                .iter()
                .find(|level| level.k == query.k)
                .ok_or_else(|| {
                    JmaError::CorruptSection(format!("seed level k={} is unavailable", query.k))
                })?;
            return self.seed_occurrences_at_scale(query, level.scale);
        }
        let section = self
            .seed_sections
            .iter()
            .find(|section| section.k == query.k)
            .ok_or_else(|| {
                JmaError::CorruptSection(format!("seed level k={} is unavailable", query.k))
            })?;
        let level = section.levels.first().ok_or_else(|| {
            JmaError::CorruptSection(format!("seed level k={} has no density data", query.k))
        })?;
        Ok(matching_occurrences(&level.records, query))
    }

    fn metrics(&self) -> ResourceMetrics {
        let mut metrics = self.resource.metrics();
        if let Some(index_resource) = self.index_resource.as_ref() {
            metrics = metrics.saturating_add(index_resource.metrics());
        }
        metrics.saturating_add(self.metrics.snapshot())
    }
}

fn read_exact<R: RangeReader>(resource: &R, range: ByteRange) -> JmaResult<Vec<u8>> {
    let expected = checked_usize(range.length, "resource read length")?;
    let bytes = resource.read_range(range)?;
    if bytes.len() != expected {
        return Err(JmaError::CorruptSection(format!(
            "resource returned {} bytes for requested {expected}",
            bytes.len()
        )));
    }
    Ok(bytes)
}

fn read_section<R: RangeReader>(
    resource: &R,
    descriptor: &SectionDescriptor,
) -> JmaResult<Vec<u8>> {
    checked_end(descriptor.offset, descriptor.length)?;
    let bytes = read_exact(
        resource,
        ByteRange::new(descriptor.offset, descriptor.length)?,
    )?;
    if checksum(&bytes) != descriptor.checksum {
        return Err(JmaError::ChecksumMismatch(format!(
            "section kind {}",
            descriptor.kind
        )));
    }
    if descriptor.uncompressed_length != descriptor.length {
        return Err(JmaError::CorruptSection(format!(
            "compressed section kind {} is unsupported in JMA v1",
            descriptor.kind
        )));
    }
    Ok(bytes)
}

fn find_unique_section(sections: &[SectionDescriptor], kind: u32) -> JmaResult<&SectionDescriptor> {
    let mut matches = sections.iter().filter(|entry| entry.kind == kind);
    let first = matches.next().ok_or_else(|| {
        JmaError::CorruptSection(format!("required section kind {kind} is missing"))
    })?;
    if matches.next().is_some() {
        return Err(JmaError::CorruptSection(format!(
            "section kind {kind} occurs more than once"
        )));
    }
    Ok(first)
}

fn validate_section_kinds(sections: &[SectionDescriptor]) -> JmaResult<()> {
    for section in sections {
        if !matches!(
            section.kind,
            SECTION_CONTIGS | SECTION_SEQUENCES | SECTION_SEEDS_K21 | SECTION_SEEDS_K31
        ) {
            return Err(JmaError::CorruptSection(format!(
                "unsupported JMA section kind {}",
                section.kind
            )));
        }
    }
    if sections
        .iter()
        .filter(|section| section.kind == SECTION_SEEDS_K21)
        .count()
        > 1
        || sections
            .iter()
            .filter(|section| section.kind == SECTION_SEEDS_K31)
            .count()
            > 1
    {
        return Err(JmaError::CorruptSection(
            "seed section kind occurs more than once".to_string(),
        ));
    }
    Ok(())
}

fn matching_occurrences(records: &[SeedRecord], query: SeedQuery) -> Vec<SeedOccurrence> {
    let mut matches = records
        .iter()
        .filter(|record| {
            record.query.hash == query.hash && record.query.canonical_kmer == query.canonical_kmer
        })
        .map(|record| record.occurrence)
        .collect::<Vec<_>>();
    matches.sort_by_key(|occurrence| {
        (
            occurrence.contig_id,
            occurrence.position,
            occurrence.reverse,
        )
    });
    matches
}
