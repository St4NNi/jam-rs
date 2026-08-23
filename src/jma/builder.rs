//! Deterministic construction of JMA v1 archives from FASTA/FASTQ input.
//!
//! The builder parses each contig once, retaining the normalized sequence long
//! enough to construct packed blocks and positional seed records.  The output
//! is an `ArchiveParts` value consumed by the checked JMA writer.  Source
//! provenance is the SHA-256 of the input file bytes (compressed bytes when a
//! compressed FASTA/FASTQ path is supplied); it is not the checksum of the JMA
//! container.

use crate::jma::seed_builder::{
    EmbeddedSampleSketch, SeedBuildConfig, SeedInput, build_seed_sections,
};
use crate::jma::sequence_builder::{SequenceBuildConfig, build_sequence_blocks};
use crate::jma::writer::{ArchiveParts, write_archive_with_min_entropy};
use crate::jma::{ContigMetadata, JmaError, JmaResult};
use needletail::{Sequence, parse_fastx_file};
use sha2::{Digest, Sha256};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

/// Build flags reserved for future JMA metadata.  No flags are set by the
/// current complete local builder.
pub const DEFAULT_FLAGS: u32 = 0;

/// Configuration for a complete JMA archive build.
#[derive(Clone, Copy, Debug)]
pub struct ArchiveBuildConfig {
    pub block_bases: usize,
    pub k31_scale: u64,
    pub k21_scale: Option<u64>,
    pub min_entropy: Option<f64>,
    pub flags: u32,
}

impl Default for ArchiveBuildConfig {
    fn default() -> Self {
        Self {
            block_bases: SequenceBuildConfig::default().block_bases,
            k31_scale: SeedBuildConfig::default().k31_scale,
            k21_scale: SeedBuildConfig::default().k21_scale,
            min_entropy: SeedBuildConfig::default().min_entropy,
            flags: DEFAULT_FLAGS,
        }
    }
}

/// Basic counts and embedded sketch evidence from an archive build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArchiveBuildStats {
    pub source_sha256: [u8; 32],
    pub contig_count: u32,
    pub total_bases: u64,
    pub k31_seed_count: u64,
    pub k21_seed_count: u64,
    pub sample_sketches: Vec<EmbeddedSampleSketch>,
}

/// Output of `build_archive_from_fasta` before it is serialized.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ArchiveBuildOutput {
    pub parts: ArchiveParts,
    pub stats: ArchiveBuildStats,
}

#[derive(Clone, Debug)]
struct InputContig {
    name: String,
    sequence: Vec<u8>,
}

/// Reads FASTA/FASTQ (including needletail-supported compressed paths) and
/// constructs deterministic JMA writer parts.
pub fn build_archive_from_fasta<P: AsRef<Path>>(
    input: P,
    config: ArchiveBuildConfig,
) -> JmaResult<ArchiveBuildOutput> {
    validate_config(config)?;
    let input_path = input.as_ref();
    let source_sha256 = checksum_file(input_path)?;
    let records = read_contigs(input_path)?;

    let mut contigs = Vec::with_capacity(records.len());
    let mut sequence_blocks = Vec::new();
    let mut seed_inputs = Vec::with_capacity(records.len());
    let mut total_bases = 0u64;

    for (index, record) in records.iter().enumerate() {
        let id = u32::try_from(index).map_err(|_| JmaError::OffsetOverflow)?;
        let length = u64::try_from(record.sequence.len()).map_err(|_| JmaError::OffsetOverflow)?;
        total_bases = total_bases
            .checked_add(length)
            .ok_or(JmaError::OffsetOverflow)?;
        contigs.push(ContigMetadata {
            id,
            name: record.name.clone(),
            length,
        });
        sequence_blocks.extend(build_sequence_blocks(
            id,
            &record.sequence,
            SequenceBuildConfig {
                block_bases: config.block_bases,
            },
        )?);
        seed_inputs.push(SeedInput {
            contig_id: id,
            sequence: &record.sequence,
        });
    }

    let seed_result = build_seed_sections(
        &seed_inputs,
        SeedBuildConfig {
            k31_scale: config.k31_scale,
            k21_scale: config.k21_scale,
            min_entropy: config.min_entropy,
        },
    )?;
    let k31_seed_count = seed_result
        .sections
        .iter()
        .find(|section| section.k == 31)
        .map_or(0, |section| {
            section
                .levels
                .iter()
                .map(|level| level.records.len() as u64)
                .sum()
        });
    let k21_seed_count = seed_result
        .sections
        .iter()
        .find(|section| section.k == 21)
        .map_or(0, |section| {
            section
                .levels
                .iter()
                .map(|level| level.records.len() as u64)
                .sum()
        });
    let contig_count = u32::try_from(contigs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let stats = ArchiveBuildStats {
        source_sha256,
        contig_count,
        total_bases,
        k31_seed_count,
        k21_seed_count,
        sample_sketches: seed_result.sample_sketches,
    };
    Ok(ArchiveBuildOutput {
        parts: ArchiveParts {
            flags: config.flags,
            source_sha256,
            contigs,
            sequence_blocks,
            seed_sections: seed_result.sections,
        },
        stats,
    })
}

/// Builds and writes a complete archive.  Serialization is delegated to the
/// checked deterministic JMA writer.
pub fn write_archive_from_fasta<P: AsRef<Path>, Q: AsRef<Path>>(
    input: P,
    output: Q,
    config: ArchiveBuildConfig,
) -> JmaResult<ArchiveBuildStats> {
    let built = build_archive_from_fasta(input, config)?;
    let mut file = File::create(output.as_ref())
        .map_err(|error| JmaError::CorruptSection(format!("cannot create JMA output: {error}")))?;
    write_archive_with_min_entropy(&mut file, &built.parts, config.min_entropy)?;
    file.flush()
        .map_err(|error| JmaError::CorruptSection(format!("cannot flush JMA output: {error}")))?;
    Ok(built.stats)
}

/// Alias for callers that use the shorter archive-builder name.
pub fn build<P: AsRef<Path>>(
    input: P,
    config: ArchiveBuildConfig,
) -> JmaResult<ArchiveBuildOutput> {
    build_archive_from_fasta(input, config)
}

fn validate_config(config: ArchiveBuildConfig) -> JmaResult<()> {
    if config.block_bases == 0 {
        return Err(JmaError::CorruptSection(
            "sequence block size must be greater than zero".to_string(),
        ));
    }
    if config.k31_scale == 0 || config.k21_scale == Some(0) {
        return Err(JmaError::CorruptSection(
            "seed scales must be greater than zero".to_string(),
        ));
    }
    Ok(())
}

fn checksum_file(path: &Path) -> JmaResult<[u8; 32]> {
    let mut file = File::open(path).map_err(|error| {
        JmaError::CorruptSection(format!("cannot open input {}: {error}", path.display()))
    })?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 1024 * 1024];
    loop {
        let count = file.read(&mut buffer).map_err(|error| {
            JmaError::CorruptSection(format!("cannot read input {}: {error}", path.display()))
        })?;
        if count == 0 {
            break;
        }
        digest.update(&buffer[..count]);
    }
    let mut checksum = [0u8; 32];
    checksum.copy_from_slice(&digest.finalize());
    Ok(checksum)
}

fn read_contigs(path: &Path) -> JmaResult<Vec<InputContig>> {
    let mut reader = match parse_fastx_file(path) {
        Ok(reader) => reader,
        Err(error) if error.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            return Ok(Vec::new());
        }
        Err(error) => {
            return Err(JmaError::CorruptSection(format!(
                "cannot parse input {}: {error}",
                path.display()
            )));
        }
    };

    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| {
            JmaError::CorruptSection(format!("cannot parse input {}: {error}", path.display()))
        })?;
        let name = std::str::from_utf8(record.id())
            .map_err(|error| {
                JmaError::CorruptSection(format!("input contig name is not UTF-8: {error}"))
            })?
            .trim()
            .to_string();
        if name.is_empty() {
            return Err(JmaError::CorruptSection(
                "input contains a contig with an empty name".to_string(),
            ));
        }
        let sequence = record.normalize(false).into_owned();
        records.push(InputContig { name, sequence });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jma::ArchiveReader;
    use crate::jma::reader::JmaReader;
    use crate::jma::writer::encode_archive;
    use crate::resource::{
        ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
    };
    use std::io::Cursor;
    use tempfile::NamedTempFile;

    struct MemoryResource {
        locator: ResourceLocator,
        bytes: Vec<u8>,
    }

    impl RangeReader for MemoryResource {
        fn locator(&self) -> &ResourceLocator {
            &self.locator
        }

        fn metadata(&self) -> Result<ResourceMetadata, ResourceError> {
            Ok(ResourceMetadata {
                size: self.bytes.len() as u64,
                etag: None,
                last_modified: None,
                accepts_ranges: true,
            })
        }

        fn read_range(&self, range: ByteRange) -> Result<Vec<u8>, ResourceError> {
            let end = range.end()?;
            if end > self.bytes.len() as u64 {
                return Err(ResourceError::RangeOutOfBounds {
                    offset: range.offset,
                    length: range.length,
                    size: self.bytes.len() as u64,
                });
            }
            Ok(self.bytes[range.offset as usize..end as usize].to_vec())
        }

        fn stream(&self) -> Result<Box<dyn Read + Send>, ResourceError> {
            Ok(Box::new(Cursor::new(self.bytes.clone())))
        }

        fn metrics(&self) -> ResourceMetrics {
            ResourceMetrics::default()
        }
    }

    fn fasta(contents: &[u8]) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(contents).unwrap();
        file
    }

    #[test]
    fn empty_input_builds_a_valid_empty_archive() {
        let file = fasta(b"");
        let output = build_archive_from_fasta(file.path(), ArchiveBuildConfig::default()).unwrap();
        assert_eq!(output.stats.contig_count, 0);
        assert_eq!(output.stats.total_bases, 0);
        let bytes = encode_archive(&output.parts).unwrap();
        assert_eq!(&bytes[..4], b"JMA\0");
    }

    #[test]
    fn archive_is_deterministic_and_roundtrips() {
        let file = fasta(
            b">first\nACGTNacgt\n>second\nACGTTGCATGTCAGTAGGCATCAGTACCGATGCTAGCTAGGCTAACGTTACGATCGATGCA\n",
        );
        let config = ArchiveBuildConfig {
            block_bases: 4,
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        };
        let first = build_archive_from_fasta(file.path(), config).unwrap();
        let second = build_archive_from_fasta(file.path(), config).unwrap();
        let bytes_one = encode_archive(&first.parts).unwrap();
        let bytes_two = encode_archive(&second.parts).unwrap();
        assert_eq!(bytes_one, bytes_two);
        assert_eq!(first.stats, second.stats);

        let resource = MemoryResource {
            locator: ResourceLocator::parse("file://jma-builder-test").unwrap(),
            bytes: bytes_one,
        };
        let reader = JmaReader::from_resource(resource).unwrap();
        assert_eq!(reader.contigs().len(), 2);
        assert_eq!(reader.contigs()[0].name, "first");
        assert_eq!(
            reader
                .read_sequence(0, crate::jma::SequenceRange::new(0, 9).unwrap())
                .unwrap(),
            b"ACGTNACGT"
        );
    }

    #[test]
    fn malformed_input_is_not_reported_as_empty() {
        let file = fasta(b"not-a-fasta-record\n");
        let error =
            build_archive_from_fasta(file.path(), ArchiveBuildConfig::default()).unwrap_err();
        assert!(error.to_string().contains("parse"));
    }
}
