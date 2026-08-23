use jam_rs::jma::JmaError;
use jam_rs::jma::contigs::decode_contigs;
use jam_rs::jma::format::{
    HEADER_SIZE, SECTION_ENTRY_SIZE, checksum, decode_header, decode_section_directory,
};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::sequence::{decode_sequence_blocks, unpack_bases};
use jam_rs::jma::writer::decode_seed_section;
use jam_rs::jma::writer::{ArchiveParts, encode_archive};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use std::io::{Cursor, Read};

struct BytesResource {
    locator: ResourceLocator,
    bytes: Vec<u8>,
}

impl BytesResource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("file://jma-corrupt").unwrap(),
            bytes,
        }
    }
}

impl RangeReader for BytesResource {
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

fn fixture() -> ArchiveParts {
    ArchiveParts {
        flags: 0,
        source_sha256: [0; 32],
        contigs: Vec::new(),
        sequence_blocks: Vec::new(),
        seed_sections: Vec::new(),
    }
}

#[test]
fn header_corruption_is_rejected_by_checksum() {
    let mut bytes = encode_archive(&fixture()).unwrap();
    bytes[8] ^= 1;
    assert!(matches!(
        decode_header(&bytes[..HEADER_SIZE]),
        Err(JmaError::ChecksumMismatch(section)) if section == "header"
    ));
}

#[test]
fn section_payload_corruption_is_rejected_by_checksum() {
    let mut bytes = encode_archive(&fixture()).unwrap();
    let directory_offset = u64::from_le_bytes(bytes[64..72].try_into().unwrap()) as usize;
    let first_section_offset = u64::from_le_bytes(
        bytes[directory_offset + 8..directory_offset + 16]
            .try_into()
            .unwrap(),
    ) as usize;
    bytes[first_section_offset] ^= 1;
    assert!(matches!(
        JmaReader::from_resource(BytesResource::new(bytes)),
        Err(JmaError::ChecksumMismatch(section)) if section.contains("section kind")
    ));
}

#[test]
fn malformed_directory_length_and_truncated_entry_are_rejected() {
    assert!(decode_section_directory(&[], 1).is_err());
    assert!(decode_section_directory(&[0; SECTION_ENTRY_SIZE - 1], 1).is_err());
}

fn refresh_header_checksum(bytes: &mut [u8]) {
    bytes[80..112].fill(0);
    let digest = checksum(&bytes[..HEADER_SIZE]);
    bytes[80..112].copy_from_slice(&digest);
}

#[test]
fn unsupported_version_and_truncated_header_are_rejected() {
    let bytes = encode_archive(&fixture()).unwrap();
    assert!(matches!(
        decode_header(&bytes[..HEADER_SIZE - 1]),
        Err(JmaError::CorruptSection(message)) if message.contains("truncated")
    ));

    let mut unsupported = bytes;
    unsupported[4..6].copy_from_slice(&2u16.to_le_bytes());
    refresh_header_checksum(&mut unsupported);
    assert!(matches!(
        decode_header(&unsupported[..HEADER_SIZE]),
        Err(JmaError::UnsupportedVersion(2))
    ));
}

#[test]
fn unknown_algorithm_tag_and_invalid_entropy_are_rejected() {
    let mut unknown_tag = encode_archive(&fixture()).unwrap();
    unknown_tag[144] = b'X';
    refresh_header_checksum(&mut unknown_tag);
    assert!(matches!(
        decode_header(&unknown_tag[..HEADER_SIZE]),
        Err(JmaError::CorruptSection(message)) if message.contains("unknown JMA algorithm tag")
    ));

    let mut invalid_entropy = encode_archive(&fixture()).unwrap();
    invalid_entropy[152..160].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    refresh_header_checksum(&mut invalid_entropy);
    assert!(matches!(
        decode_header(&invalid_entropy[..HEADER_SIZE]),
        Err(JmaError::CorruptSection(message)) if message.contains("invalid JMA minimum entropy")
    ));

    let mut entropy_without_tag = encode_archive(&fixture()).unwrap();
    entropy_without_tag[144..152].fill(0);
    entropy_without_tag[152..160].copy_from_slice(&1.0f64.to_le_bytes());
    refresh_header_checksum(&mut entropy_without_tag);
    assert!(matches!(
        decode_header(&entropy_without_tag[..HEADER_SIZE]),
        Err(JmaError::CorruptSection(message)) if message.contains("without a JMA algorithm tag")
    ));
}

#[test]
fn section_offsets_outside_resource_are_rejected() {
    let mut bytes = encode_archive(&fixture()).unwrap();
    let directory_offset = u64::from_le_bytes(bytes[64..72].try_into().unwrap()) as usize;
    let invalid_offset = bytes.len() as u64 + 1;
    bytes[directory_offset + 8..directory_offset + 16]
        .copy_from_slice(&invalid_offset.to_le_bytes());
    refresh_header_checksum(&mut bytes);
    assert!(matches!(
        JmaReader::from_resource(BytesResource::new(bytes)),
        Err(JmaError::CorruptSection(message)) if message.contains("ends at")
    ));
}

#[test]
fn oversized_declared_counts_and_lengths_do_not_allocate_or_panic() {
    let mut contig_header = Vec::new();
    contig_header.extend_from_slice(&1u32.to_le_bytes());
    contig_header.extend_from_slice(&u32::MAX.to_le_bytes());
    assert!(decode_contigs(&contig_header, u32::MAX).is_err());

    let mut sequence_header = Vec::new();
    sequence_header.extend_from_slice(&1u32.to_le_bytes());
    sequence_header.extend_from_slice(&u32::MAX.to_le_bytes());
    assert!(decode_sequence_blocks(&sequence_header).is_err());

    let mut seed_header = Vec::new();
    seed_header.extend_from_slice(&1u32.to_le_bytes());
    seed_header.push(31);
    seed_header.extend_from_slice(&[0; 3]);
    seed_header.extend_from_slice(&u32::MAX.to_le_bytes());
    assert!(decode_seed_section(&seed_header).is_err());

    assert!(unpack_bases(&[], &[], u64::MAX).is_err());

    let mut archive = encode_archive(&fixture()).unwrap();
    let excessive_count = (1u32 << 20) + 1;
    archive[24..28].copy_from_slice(&excessive_count.to_le_bytes());
    archive[72..80]
        .copy_from_slice(&(u64::from(excessive_count) * SECTION_ENTRY_SIZE as u64).to_le_bytes());
    refresh_header_checksum(&mut archive);
    assert!(matches!(
        JmaReader::from_resource(BytesResource::new(archive)),
        Err(JmaError::CorruptSection(message)) if message.contains("section count")
    ));
}

#[test]
fn single_byte_mutations_never_panic_in_archive_reader() {
    let bytes = encode_archive(&fixture()).unwrap();
    for index in 0..bytes.len() {
        let mut mutated = bytes.clone();
        mutated[index] ^= 0x01;
        let result = std::panic::catch_unwind(|| {
            let _ = JmaReader::from_resource(BytesResource::new(mutated));
        });
        assert!(result.is_ok(), "mutation at byte {index} panicked");
    }
}
