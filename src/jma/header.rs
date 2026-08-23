//! JMA v1 fixed-header and section-directory validation.

use crate::jma::format::{
    HEADER_SIZE, ParsedHeader, SECTION_ENTRY_SIZE, SectionDescriptor, checked_end, decode_header,
    decode_section_directory,
};
use crate::jma::{JmaError, JmaResult};

/// Decodes a fixed header and ensures the header's directory arithmetic is
/// representable. Resource-size checks are performed by the reader once
/// resource metadata is available.
pub fn parse_header(bytes: &[u8]) -> JmaResult<ParsedHeader> {
    let parsed = decode_header(bytes)?;
    if parsed.section_directory_offset < HEADER_SIZE as u64 {
        return Err(JmaError::CorruptSection(format!(
            "section directory starts before the {HEADER_SIZE}-byte header"
        )));
    }
    let expected_length = (parsed.section_count as u64)
        .checked_mul(SECTION_ENTRY_SIZE as u64)
        .ok_or(JmaError::OffsetOverflow)?;
    if parsed.section_directory_length != expected_length {
        return Err(JmaError::CorruptSection(format!(
            "section directory length {} does not match count {}",
            parsed.section_directory_length, parsed.section_count
        )));
    }
    checked_end(
        parsed.section_directory_offset,
        parsed.section_directory_length,
    )?;
    Ok(parsed)
}

/// Decodes and validates a section directory against the complete resource
/// size. Entries may be listed in any order, but their byte ranges may not
/// overlap the header or directory and may not overlap one another.
pub fn parse_section_directory(
    bytes: &[u8],
    section_count: u32,
    directory_offset: u64,
    resource_size: u64,
) -> JmaResult<Vec<SectionDescriptor>> {
    let entries = decode_section_directory(bytes, section_count)?;
    let directory_end = checked_end(directory_offset, bytes.len() as u64)?;
    if directory_end > resource_size {
        return Err(JmaError::CorruptSection(format!(
            "section directory ends at {directory_end}, resource has {resource_size} bytes"
        )));
    }

    let mut ordered = entries.clone();
    ordered.sort_by_key(|entry| (entry.offset, entry.length, entry.kind));
    let mut previous_end = directory_end;
    for entry in &ordered {
        let end = checked_end(entry.offset, entry.length)?;
        if entry.offset < directory_end {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} overlaps the section directory",
                entry.kind
            )));
        }
        if end > resource_size {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} ends at {end}, resource has {resource_size} bytes",
                entry.kind
            )));
        }
        if entry.offset < previous_end {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} overlaps another section",
                entry.kind
            )));
        }
        previous_end = end;
    }
    Ok(entries)
}
