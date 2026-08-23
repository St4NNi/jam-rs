//! Deterministic JMA v1 archive and seed-section writer.

use crate::jma::contigs::encode_contigs;
use crate::jma::format::{
    HEADER_SIZE, SECTION_CONTIGS, SECTION_SEEDS_K21, SECTION_SEEDS_K31, SECTION_SEQUENCES,
    SECTION_VERSION, SectionDescriptor, checksum, encode_header, encode_section_directory,
};
use crate::jma::sequence::{SequenceBlock, encode_sequence_blocks, validate_blocks};
use crate::jma::{
    ArchiveHeader, ContigMetadata, JMA_FORMAT_VERSION, JmaError, JmaResult, SeedLevel,
    SeedOccurrence, SeedQuery,
};
use std::io::Write;

const SEED_RECORD_SIZE: usize = 32;

/// One exact seed occurrence in a JMA seed section.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedRecord {
    pub query: SeedQuery,
    pub occurrence: SeedOccurrence,
}

/// Records retained for one density level. Levels must be nested by the
/// caller; the writer preserves their explicit scale identities.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedLevelData {
    pub level: SeedLevel,
    pub records: Vec<SeedRecord>,
}

/// One seed section, normally k=31 primary or k=21 rescue.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SeedSection {
    pub k: u8,
    pub levels: Vec<SeedLevelData>,
}

/// Input to the deterministic archive writer. Sequence blocks must cover the
/// declared contigs; no sequence is silently discarded.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArchiveParts {
    pub flags: u32,
    pub source_sha256: [u8; 32],
    pub contigs: Vec<ContigMetadata>,
    pub sequence_blocks: Vec<SequenceBlock>,
    pub seed_sections: Vec<SeedSection>,
}

/// Compatibility-friendly alias for callers that call the input an archive.
pub type JmaArchive = ArchiveParts;

/// Encodes one seed section with exact canonical k-mer values in each record.
pub fn encode_seed_section(section: &SeedSection) -> JmaResult<Vec<u8>> {
    validate_seed_section(section)?;
    let level_count = u32::try_from(section.levels.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.push(section.k);
    bytes.extend_from_slice(&[0u8; 3]);
    bytes.extend_from_slice(&level_count.to_le_bytes());
    for level in &section.levels {
        bytes.extend_from_slice(&level.level.scale.to_le_bytes());
        let mut records = level.records.clone();
        records.sort_by_key(|record| {
            (
                record.query.hash,
                record.query.canonical_kmer,
                record.occurrence.contig_id,
                record.occurrence.position,
                record.occurrence.reverse,
            )
        });
        let count = u64::try_from(records.len()).map_err(|_| JmaError::OffsetOverflow)?;
        bytes.extend_from_slice(&count.to_le_bytes());
        for record in &records {
            bytes.extend_from_slice(&record.query.hash.to_le_bytes());
            bytes.extend_from_slice(&record.query.canonical_kmer.to_le_bytes());
            bytes.extend_from_slice(&record.occurrence.contig_id.to_le_bytes());
            bytes.push(u8::from(record.occurrence.reverse));
            bytes.extend_from_slice(&[0u8; 3]);
            bytes.extend_from_slice(&record.occurrence.position.to_le_bytes());
        }
    }
    Ok(bytes)
}

/// Decodes one complete seed section.
pub fn decode_seed_section(bytes: &[u8]) -> JmaResult<SeedSection> {
    if bytes.len() < 12 {
        return Err(JmaError::CorruptSection(
            "seed section header is truncated".to_string(),
        ));
    }
    let version = u32::from_le_bytes(bytes[0..4].try_into().expect("checked length"));
    if version != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported seed section version".to_string(),
        ));
    }
    let k = bytes[4];
    if !(1..=32).contains(&k) {
        return Err(JmaError::CorruptSection(format!("invalid seed k={k}")));
    }
    if bytes[5] != 0 || bytes[6] != 0 || bytes[7] != 0 {
        return Err(JmaError::CorruptSection(
            "non-zero seed section reserved bytes".to_string(),
        ));
    }
    let level_count = u32::from_le_bytes(bytes[8..12].try_into().expect("checked length"));
    let minimum_level_bytes = (level_count as usize)
        .checked_mul(16)
        .ok_or(JmaError::OffsetOverflow)?;
    let minimum_length = 12usize
        .checked_add(minimum_level_bytes)
        .ok_or(JmaError::OffsetOverflow)?;
    if bytes.len() < minimum_length {
        return Err(JmaError::CorruptSection(format!(
            "seed section is too short for {level_count} levels"
        )));
    }
    let mut cursor = 12usize;
    let mut levels = Vec::with_capacity(level_count as usize);
    for _ in 0..level_count {
        let header_end = cursor.checked_add(16).ok_or(JmaError::OffsetOverflow)?;
        if header_end > bytes.len() {
            return Err(JmaError::CorruptSection(
                "seed level header is truncated".to_string(),
            ));
        }
        let scale = u64::from_le_bytes(
            bytes[cursor..cursor + 8]
                .try_into()
                .expect("checked length"),
        );
        let record_count = u64::from_le_bytes(
            bytes[cursor + 8..cursor + 16]
                .try_into()
                .expect("checked length"),
        );
        if scale == 0 {
            return Err(JmaError::CorruptSection(
                "seed level scale must be non-zero".to_string(),
            ));
        }
        cursor = header_end;
        let record_count_usize = usize::try_from(record_count)
            .map_err(|_| JmaError::CorruptSection("seed record count is too large".to_string()))?;
        let record_bytes = record_count_usize
            .checked_mul(SEED_RECORD_SIZE)
            .ok_or(JmaError::OffsetOverflow)?;
        let end = cursor
            .checked_add(record_bytes)
            .ok_or(JmaError::OffsetOverflow)?;
        if end > bytes.len() {
            return Err(JmaError::CorruptSection(
                "seed records are truncated".to_string(),
            ));
        }
        let mut records = Vec::with_capacity(record_count_usize);
        for _ in 0..record_count_usize {
            let hash = u64::from_le_bytes(
                bytes[cursor..cursor + 8]
                    .try_into()
                    .expect("checked length"),
            );
            let canonical_kmer = u64::from_le_bytes(
                bytes[cursor + 8..cursor + 16]
                    .try_into()
                    .expect("checked length"),
            );
            let contig_id = u32::from_le_bytes(
                bytes[cursor + 16..cursor + 20]
                    .try_into()
                    .expect("checked length"),
            );
            let reverse = match bytes[cursor + 20] {
                0 => false,
                1 => true,
                value => {
                    return Err(JmaError::CorruptSection(format!(
                        "invalid seed orientation flag {value}"
                    )));
                }
            };
            if bytes[cursor + 21] != 0 || bytes[cursor + 22] != 0 || bytes[cursor + 23] != 0 {
                return Err(JmaError::CorruptSection(
                    "non-zero seed record reserved bytes".to_string(),
                ));
            }
            let position = u64::from_le_bytes(
                bytes[cursor + 24..cursor + 32]
                    .try_into()
                    .expect("checked length"),
            );
            records.push(SeedRecord {
                query: SeedQuery {
                    k,
                    hash,
                    canonical_kmer,
                },
                occurrence: SeedOccurrence {
                    contig_id,
                    position,
                    reverse,
                },
            });
            cursor += SEED_RECORD_SIZE;
        }
        levels.push(SeedLevelData {
            level: SeedLevel { k, scale },
            records,
        });
    }
    if cursor != bytes.len() {
        return Err(JmaError::CorruptSection(format!(
            "seed section has {} trailing bytes",
            bytes.len() - cursor
        )));
    }
    let section = SeedSection { k, levels };
    validate_seed_section(&section)?;
    Ok(section)
}

/// Decodes a byte range containing only fixed-size seed records.  This is
/// intentionally separate from [`decode_seed_section`], which also expects a
/// complete section header and all density levels.
pub(crate) fn decode_seed_records(bytes: &[u8], k: u8) -> JmaResult<Vec<SeedRecord>> {
    if !(1..=32).contains(&k) {
        return Err(JmaError::CorruptSection(format!("invalid seed k={k}")));
    }
    if !bytes.len().is_multiple_of(SEED_RECORD_SIZE) {
        return Err(JmaError::CorruptSection(
            "seed record range is not a whole number of records".to_string(),
        ));
    }
    let mut records = Vec::with_capacity(bytes.len() / SEED_RECORD_SIZE);
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        let hash = u64::from_le_bytes(
            bytes[cursor..cursor + 8]
                .try_into()
                .expect("record range checked"),
        );
        let canonical_kmer = u64::from_le_bytes(
            bytes[cursor + 8..cursor + 16]
                .try_into()
                .expect("record range checked"),
        );
        let contig_id = u32::from_le_bytes(
            bytes[cursor + 16..cursor + 20]
                .try_into()
                .expect("record range checked"),
        );
        let reverse = match bytes[cursor + 20] {
            0 => false,
            1 => true,
            value => {
                return Err(JmaError::CorruptSection(format!(
                    "invalid seed orientation flag {value}"
                )));
            }
        };
        if bytes[cursor + 21] != 0 || bytes[cursor + 22] != 0 || bytes[cursor + 23] != 0 {
            return Err(JmaError::CorruptSection(
                "non-zero seed record reserved bytes".to_string(),
            ));
        }
        let position = u64::from_le_bytes(
            bytes[cursor + 24..cursor + 32]
                .try_into()
                .expect("record range checked"),
        );
        if !valid_canonical_kmer(k, canonical_kmer) {
            return Err(JmaError::CorruptSection(
                "seed record contains an invalid packed k-mer".to_string(),
            ));
        }
        records.push(SeedRecord {
            query: SeedQuery {
                k,
                hash,
                canonical_kmer,
            },
            occurrence: SeedOccurrence {
                contig_id,
                position,
                reverse,
            },
        });
        cursor += SEED_RECORD_SIZE;
    }
    Ok(records)
}

/// Encodes an entire JMA v1 archive in memory. The output order and all
/// record ordering are deterministic for a deterministic `ArchiveParts` input.
pub fn encode_archive(parts: &ArchiveParts) -> JmaResult<Vec<u8>> {
    encode_archive_with_min_entropy(parts, None)
}

/// Encodes an archive and records the entropy threshold used while building
/// its positional seed sections.  The public `encode_archive` entry point
/// remains compatible for callers that do not have a threshold to record.
pub fn encode_archive_with_min_entropy(
    parts: &ArchiveParts,
    min_entropy: Option<f64>,
) -> JmaResult<Vec<u8>> {
    let contig_count = u32::try_from(parts.contigs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let total_bases = parts.contigs.iter().try_fold(0u64, |total, contig| {
        total
            .checked_add(contig.length)
            .ok_or(JmaError::OffsetOverflow)
    })?;
    validate_contig_table(&parts.contigs)?;
    validate_blocks(&parts.sequence_blocks, &parts.contigs)?;
    validate_complete_coverage(&parts.sequence_blocks, &parts.contigs)?;

    let contig_payload = encode_contigs(&parts.contigs)?;
    let sequence_payload = encode_sequence_blocks(&parts.sequence_blocks)?;
    let mut payloads = vec![
        (SECTION_CONTIGS, contig_payload),
        (SECTION_SEQUENCES, sequence_payload),
    ];
    let mut seed_levels = Vec::new();
    let mut sections = parts.seed_sections.clone();
    sections.sort_by_key(|section| section.k);
    for section in &sections {
        let kind = match section.k {
            21 => SECTION_SEEDS_K21,
            31 => SECTION_SEEDS_K31,
            other => {
                return Err(JmaError::CorruptSection(format!(
                    "JMA v1 seed section k={other} is not k=21 or k=31"
                )));
            }
        };
        let payload = encode_seed_section(section)?;
        payloads.push((kind, payload));
        seed_levels.extend(section.levels.iter().map(|level| level.level));
    }
    if seed_levels.len() > 2 {
        return Err(JmaError::CorruptSection(
            "JMA v1 header cannot describe more than two seed levels".to_string(),
        ));
    }
    if payloads.iter().enumerate().any(|(index, (kind, _))| {
        payloads
            .iter()
            .take(index)
            .any(|(previous_kind, _)| previous_kind == kind)
    }) {
        return Err(JmaError::CorruptSection(
            "duplicate JMA section kind".to_string(),
        ));
    }

    let directory_offset = HEADER_SIZE as u64;
    let directory_length = u64::try_from(
        payloads
            .len()
            .checked_mul(crate::jma::format::SECTION_ENTRY_SIZE)
            .ok_or(JmaError::OffsetOverflow)?,
    )
    .map_err(|_| JmaError::OffsetOverflow)?;
    let mut next_offset = directory_offset
        .checked_add(directory_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let mut descriptors = Vec::with_capacity(payloads.len());
    for (kind, payload) in &payloads {
        let length = u64::try_from(payload.len()).map_err(|_| JmaError::OffsetOverflow)?;
        descriptors.push(SectionDescriptor {
            kind: *kind,
            flags: 0,
            offset: next_offset,
            length,
            uncompressed_length: length,
            checksum: checksum(payload),
        });
        next_offset = next_offset
            .checked_add(length)
            .ok_or(JmaError::OffsetOverflow)?;
    }
    let archive_header = ArchiveHeader {
        format_version: JMA_FORMAT_VERSION,
        flags: parts.flags,
        contig_count,
        total_bases,
        source_sha256: parts.source_sha256,
        seed_levels,
        algorithm_id: Some("jam-seed-chain-align-v1".to_string()),
        algorithm_version: Some(1),
        min_entropy,
    };
    let header = encode_header(
        &archive_header,
        u32::try_from(descriptors.len()).map_err(|_| JmaError::OffsetOverflow)?,
        directory_offset,
        directory_length,
    )?;
    let directory = encode_section_directory(&descriptors)?;
    let capacity = usize::try_from(next_offset).map_err(|_| JmaError::OffsetOverflow)?;
    let mut output = Vec::with_capacity(capacity);
    output.extend_from_slice(&header);
    output.extend_from_slice(&directory);
    for (_, payload) in payloads {
        output.extend_from_slice(&payload);
    }
    Ok(output)
}

/// Streams an already encoded archive to a writer. The JMA error type keeps
/// I/O details redacted and leaves the caller's writer ownership untouched.
pub fn write_archive<W: Write>(writer: &mut W, parts: &ArchiveParts) -> JmaResult<()> {
    write_archive_with_min_entropy(writer, parts, None)
}

/// Streams an archive while recording the threshold used for its seed
/// sections.  This is the builder-facing form of [`write_archive`].
pub fn write_archive_with_min_entropy<W: Write>(
    writer: &mut W,
    parts: &ArchiveParts,
    min_entropy: Option<f64>,
) -> JmaResult<()> {
    let bytes = encode_archive_with_min_entropy(parts, min_entropy)?;
    writer
        .write_all(&bytes)
        .map_err(|error| JmaError::CorruptSection(format!("archive write failed: {error}")))
}

fn validate_contig_table(contigs: &[ContigMetadata]) -> JmaResult<()> {
    let mut previous = None;
    for contig in contigs {
        if previous.is_some_and(|previous_id| contig.id <= previous_id) {
            return Err(JmaError::CorruptSection(
                "contig identifiers must be strictly increasing".to_string(),
            ));
        }
        previous = Some(contig.id);
    }
    Ok(())
}

fn validate_complete_coverage(
    blocks: &[SequenceBlock],
    contigs: &[ContigMetadata],
) -> JmaResult<()> {
    for contig in contigs {
        let mut next = 0u64;
        for block in blocks.iter().filter(|block| block.contig_id == contig.id) {
            if block.start != next {
                return Err(JmaError::CorruptSection(format!(
                    "sequence blocks leave a gap or overlap in contig {} at {}",
                    contig.id, next
                )));
            }
            next = block.end().ok_or(JmaError::OffsetOverflow)?;
        }
        if next != contig.length {
            return Err(JmaError::CorruptSection(format!(
                "sequence blocks cover {next} of contig {}, expected {}",
                contig.id, contig.length
            )));
        }
    }
    Ok(())
}

fn validate_seed_section(section: &SeedSection) -> JmaResult<()> {
    if !(1..=32).contains(&section.k) {
        return Err(JmaError::CorruptSection(format!(
            "invalid seed section k={}",
            section.k
        )));
    }
    let mut previous_scale = None;
    for level in &section.levels {
        if level.level.k != section.k || level.level.scale == 0 {
            return Err(JmaError::CorruptSection(
                "seed level identity does not match its section".to_string(),
            ));
        }
        if previous_scale.is_some_and(|previous| level.level.scale >= previous) {
            return Err(JmaError::CorruptSection(
                "seed levels must be ordered from densest to sparsest".to_string(),
            ));
        }
        previous_scale = Some(level.level.scale);
        for record in &level.records {
            if record.query.k != section.k
                || !valid_canonical_kmer(section.k, record.query.canonical_kmer)
            {
                return Err(JmaError::CorruptSection(
                    "seed record query does not match its section".to_string(),
                ));
            }
        }
    }
    Ok(())
}

fn valid_canonical_kmer(k: u8, value: u64) -> bool {
    k == 32 || value < (1u64 << (2 * u32::from(k)))
}
