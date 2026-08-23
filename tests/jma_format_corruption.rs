use jam_rs::jma::JmaError;
use jam_rs::jma::format::{
    HEADER_SIZE, SECTION_ENTRY_SIZE, decode_header, decode_section_directory,
};
use jam_rs::jma::reader::JmaReader;
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
