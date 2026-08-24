use jam_rs::jma::ContigMetadata;
use jam_rs::jma::JmaError;
use jam_rs::jma::format::{
    HEADER_SIZE, SECTION_ENTRY_SIZE, SECTION_FLAG_COMPRESSED, SECTION_FLAG_REQUIRED,
    SECTION_SEQUENCE_PAYLOAD, SectionDescriptor, checksum, decode_header, decode_section_directory,
    encode_section_directory,
};
use jam_rs::jma::header::{parse_section_directory, validate_known_sections};
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::writer::{ArchiveParts, encode_archive};
use jam_rs::jma::{ArchiveReader, SequenceRange};
use jam_rs::resource::{
    ByteRange, RangeReader, ResourceError, ResourceLocator, ResourceMetadata, ResourceMetrics,
};
use jam_rs::sequence::{BlockCodec, encode_sequence_block};
use std::io::{Cursor, Read};

struct BytesResource {
    locator: ResourceLocator,
    bytes: Vec<u8>,
}
impl BytesResource {
    fn new(bytes: Vec<u8>) -> Self {
        Self {
            locator: ResourceLocator::parse("file:///jma-corrupt").unwrap(),
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
    let sequence = b"ACGTACGTACGTACGT";
    let sequence_block = encode_sequence_block(0, 0, sequence, BlockCodec::Raw2Bit).unwrap();
    ArchiveParts {
        flags: 0,
        source_sha256: [1; 32],
        contigs: vec![ContigMetadata {
            id: 0,
            name: "c".into(),
            length: sequence.len() as u64,
        }],
        sequence_blocks: vec![sequence_block],
        seed_sections: Vec::new(),
    }
}

#[test]
fn superblock_and_section_checksum_corruption_fails_closed() {
    let bytes = encode_archive(&fixture()).unwrap();
    let mut bad_header = bytes.clone();
    bad_header[16] ^= 1;
    assert!(matches!(
        decode_header(&bad_header[..HEADER_SIZE]),
        Err(JmaError::ChecksumMismatch(_))
    ));

    let mut bad_section = bytes;
    let directory_offset = u64::from_le_bytes(bad_section[24..32].try_into().unwrap()) as usize;
    let section_offset = u64::from_le_bytes(
        bad_section[directory_offset + 8..directory_offset + 16]
            .try_into()
            .unwrap(),
    ) as usize;
    bad_section[section_offset] ^= 1;
    assert!(matches!(
        JmaReader::open(BytesResource::new(bad_section)),
        Err(JmaError::ChecksumMismatch(_))
    ));

    let mut bad_directory = encode_archive(&fixture()).unwrap();
    let directory_offset = u64::from_le_bytes(bad_directory[24..32].try_into().unwrap()) as usize;
    bad_directory[directory_offset] ^= 1;
    assert!(matches!(
        JmaReader::open(BytesResource::new(bad_directory)),
        Err(JmaError::ChecksumMismatch(section)) if section == "section directory"
    ));
}

#[test]
fn section_directory_rejects_bounds_overlap_and_unknown_required() {
    assert!(decode_section_directory(&[0; SECTION_ENTRY_SIZE - 1], 1).is_err());
    let entry = SectionDescriptor {
        kind: 99,
        flags: SECTION_FLAG_REQUIRED,
        offset: 400,
        length: 4,
        uncompressed_length: 4,
        checksum: checksum(b"test"),
    };
    let mut encoded = Vec::new();
    encoded.extend_from_slice(&entry.kind.to_le_bytes());
    encoded.extend_from_slice(&entry.flags.to_le_bytes());
    encoded.extend_from_slice(&entry.offset.to_le_bytes());
    encoded.extend_from_slice(&entry.length.to_le_bytes());
    encoded.extend_from_slice(&entry.uncompressed_length.to_le_bytes());
    encoded.extend_from_slice(&entry.checksum);
    let parsed = parse_section_directory(&encoded, 1, 256, 512).unwrap();
    assert!(validate_known_sections(&parsed).is_err());
    assert!(parse_section_directory(&encoded, 1, 256, 302).is_err());
}

#[test]
fn sequence_block_checksum_rejects_corrupt_raw_or_zstd_payload() {
    for codec in [BlockCodec::Raw2Bit, BlockCodec::Zstd2Bit] {
        let block = encode_sequence_block(0, 0, b"ACGTRYKMBVDHN", codec).unwrap();
        let mut bytes = encode_archive(&ArchiveParts {
            flags: 0,
            source_sha256: [2; 32],
            contigs: vec![ContigMetadata {
                id: 0,
                name: "c".into(),
                length: 13,
            }],
            sequence_blocks: vec![block],
            seed_sections: Vec::new(),
        })
        .unwrap();
        let header = decode_header(&bytes[..HEADER_SIZE]).unwrap();
        let directory_start = usize::try_from(header.section_directory_offset).unwrap();
        let directory_end =
            directory_start + usize::try_from(header.section_directory_length).unwrap();
        let sections =
            decode_section_directory(&bytes[directory_start..directory_end], header.section_count)
                .unwrap();
        let payload = sections
            .iter()
            .find(|section| section.kind == SECTION_SEQUENCE_PAYLOAD)
            .unwrap();
        bytes[usize::try_from(payload.offset).unwrap()] ^= 1;
        let reader = JmaReader::open(BytesResource::new(bytes)).unwrap();
        assert!(
            reader
                .read_sequence(0, SequenceRange::new(0, 13).unwrap())
                .is_err()
        );
    }
}

#[test]
fn required_section_flags_and_compression_lengths_are_checked_conditionally() {
    let compressed = SectionDescriptor {
        kind: 99,
        flags: SECTION_FLAG_COMPRESSED,
        offset: 400,
        length: 4,
        uncompressed_length: 128,
        checksum: checksum(b"test"),
    };
    let encoded = encode_section_directory(&[compressed]).unwrap();
    assert!(decode_section_directory(&encoded, 1).is_ok());

    let mut uncompressed = compressed;
    uncompressed.flags = 0;
    assert!(encode_section_directory(&[uncompressed]).is_err());

    let missing_required = [
        SectionDescriptor {
            kind: 1,
            flags: 0,
            ..compressed
        },
        SectionDescriptor {
            kind: 2,
            flags: 0,
            ..compressed
        },
        SectionDescriptor {
            kind: 3,
            flags: 0,
            ..compressed
        },
        SectionDescriptor {
            kind: 4,
            flags: 0,
            ..compressed
        },
        SectionDescriptor {
            kind: 5,
            flags: 0,
            ..compressed
        },
        SectionDescriptor {
            kind: 6,
            flags: 0,
            ..compressed
        },
    ];
    assert!(validate_known_sections(&missing_required).is_err());
    let duplicate = [
        SectionDescriptor {
            kind: 1,
            flags: SECTION_FLAG_REQUIRED,
            ..compressed
        },
        SectionDescriptor {
            kind: 1,
            flags: 0,
            ..compressed
        },
    ];
    assert!(validate_known_sections(&duplicate).is_err());
    let unsupported_flags = [SectionDescriptor {
        kind: 1,
        flags: 8,
        ..compressed
    }];
    assert!(validate_known_sections(&unsupported_flags).is_err());
}
