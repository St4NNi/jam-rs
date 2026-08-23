//! Checked random-access reader for JMA v1 resources.

use crate::jma::contigs::decode_contigs;
use crate::jma::format::{
    HEADER_SIZE, SECTION_CONTIGS, SECTION_SEEDS_K21, SECTION_SEEDS_K31, SECTION_SEQUENCES,
    SectionDescriptor, checked_end, checked_usize, checksum,
};
use crate::jma::header::{parse_header, parse_section_directory};
use crate::jma::sequence::{SequenceBlock, decode_range, decode_sequence_blocks, validate_blocks};
use crate::jma::writer::{SeedRecord, SeedSection, decode_seed_section};
use crate::jma::{
    ArchiveHeader, ArchiveReader, ContigId, ContigMetadata, JmaError, JmaResult, SeedOccurrence,
    SeedQuery, SequenceRange,
};
use crate::resource::{ByteRange, RangeReader};

const MAX_SECTION_COUNT: u32 = 1 << 20;

/// A fully validated JMA archive backed by a range-readable resource.
pub struct JmaReader<R> {
    resource: R,
    header: ArchiveHeader,
    sections: Vec<SectionDescriptor>,
    contigs: Vec<ContigMetadata>,
    sequence_blocks: Vec<SequenceBlock>,
    seed_sections: Vec<SeedSection>,
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
            header: parsed.archive,
            sections,
            contigs,
            sequence_blocks,
            seed_sections,
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
        decode_range(&self.sequence_blocks, &self.contigs, contig_id, range)
    }

    fn seed_occurrences(&self, query: SeedQuery) -> JmaResult<Vec<SeedOccurrence>> {
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
