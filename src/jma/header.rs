//! Checked superblock and section-directory handling for JMA format 1.

use crate::jma::format::{
    HEADER_SIZE, SECTION_CONTIGS, SECTION_FLAG_COMPRESSED, SECTION_FLAG_REQUIRED,
    SECTION_GEAR_DIRECTORY, SECTION_METADATA, SECTION_SEEDS, SECTION_SEQUENCE_DIRECTORY,
    SECTION_SEQUENCE_PAYLOAD, SECTION_SKETCH, SectionDescriptor, checked_end, decode_header,
    decode_section_directory,
};
use crate::jma::{JmaError, JmaResult};

/// Decodes a fixed superblock and validates directory arithmetic. Resource
/// size is checked by [`parse_section_directory`] once metadata is available.
pub fn parse_header(bytes: &[u8]) -> JmaResult<crate::jma::format::ParsedHeader> {
    let parsed = decode_header(bytes)?;
    let header_size = u64::try_from(HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    if parsed.section_directory_offset < header_size {
        return Err(JmaError::CorruptSection(format!(
            "section directory starts before the {HEADER_SIZE}-byte superblock"
        )));
    }
    let expected_length = u64::from(parsed.section_count)
        .checked_mul(
            u64::try_from(crate::jma::format::SECTION_ENTRY_SIZE)
                .map_err(|_| JmaError::OffsetOverflow)?,
        )
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
/// size. Entries may be listed in any order, but may not overlap the
/// superblock, directory, or one another. Unknown optional entries are kept so
/// callers can skip them; unknown required entries are rejected by
/// [`validate_known_sections`].
pub fn parse_section_directory(
    bytes: &[u8],
    section_count: u32,
    directory_offset: u64,
    resource_size: u64,
) -> JmaResult<Vec<SectionDescriptor>> {
    let entries = decode_section_directory(bytes, section_count)?;
    let directory_length = u64::try_from(bytes.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let directory_end = checked_end(directory_offset, directory_length)?;
    let header_size = u64::try_from(HEADER_SIZE).map_err(|_| JmaError::OffsetOverflow)?;
    if directory_offset < header_size || directory_end > resource_size {
        return Err(JmaError::CorruptSection(format!(
            "section directory [{directory_offset}, {directory_end}) is outside resource size {resource_size}"
        )));
    }

    let mut ordered = entries.clone();
    ordered.sort_by_key(|entry| (entry.offset, entry.length, entry.kind));
    let mut previous_end = directory_end;
    for entry in &ordered {
        let end = checked_end(entry.offset, entry.length)?;
        if entry.offset < header_size {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} overlaps the superblock",
                entry.kind
            )));
        }
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

/// Rejects duplicate required sections and unknown required kinds while
/// leaving optional extensions available for a caller that understands them.
pub fn validate_known_sections(sections: &[SectionDescriptor]) -> JmaResult<()> {
    let mut required = std::collections::BTreeSet::new();
    let intrinsic_required = [
        SECTION_METADATA,
        SECTION_CONTIGS,
        SECTION_SKETCH,
        SECTION_SEEDS,
        SECTION_SEQUENCE_DIRECTORY,
        SECTION_SEQUENCE_PAYLOAD,
    ];
    for section in sections {
        if section.flags & !(SECTION_FLAG_REQUIRED | SECTION_FLAG_COMPRESSED) != 0 {
            return Err(JmaError::CorruptSection(format!(
                "section kind {} has unsupported flags 0x{:08x}",
                section.kind, section.flags
            )));
        }
        let intrinsic = intrinsic_required.contains(&section.kind);
        if intrinsic && section.flags & SECTION_FLAG_REQUIRED == 0 {
            return Err(JmaError::CorruptSection(format!(
                "required section kind {} is not marked required",
                section.kind
            )));
        }
        if intrinsic && !required.insert(section.kind) {
            return Err(JmaError::CorruptSection(format!(
                "duplicate required section kind {}",
                section.kind
            )));
        }
        let known = matches!(
            section.kind,
            SECTION_METADATA
                | SECTION_CONTIGS
                | SECTION_SKETCH
                | SECTION_SEEDS
                | SECTION_SEQUENCE_DIRECTORY
                | SECTION_SEQUENCE_PAYLOAD
                | SECTION_GEAR_DIRECTORY
        );
        if !known && section.flags & SECTION_FLAG_REQUIRED != 0 {
            return Err(JmaError::CorruptSection(format!(
                "unknown required section kind {}",
                section.kind
            )));
        }
    }
    for kind in intrinsic_required {
        if !sections.iter().any(|entry| entry.kind == kind) {
            return Err(JmaError::CorruptSection(format!(
                "required section kind {kind} is missing"
            )));
        }
    }
    Ok(())
}

/// Returns the unique descriptor for a required section.
pub fn unique_section(sections: &[SectionDescriptor], kind: u32) -> JmaResult<&SectionDescriptor> {
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
