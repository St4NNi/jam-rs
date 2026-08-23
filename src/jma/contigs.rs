//! Contig-table encoding for JMA v1.

use crate::jma::format::{SECTION_VERSION, checked_slice, get_u32, get_u64};
use crate::jma::{ContigId, ContigMetadata, JmaError, JmaResult};

/// Encodes contig identifiers, lengths, and UTF-8 names with explicit
/// little-endian fields.
pub fn encode_contigs(contigs: &[ContigMetadata]) -> JmaResult<Vec<u8>> {
    let count = u32::try_from(contigs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.extend_from_slice(&count.to_le_bytes());
    let mut previous_id = None;
    for contig in contigs {
        if previous_id.is_some_and(|previous| contig.id <= previous) {
            return Err(JmaError::CorruptSection(
                "contig identifiers must be strictly increasing".to_string(),
            ));
        }
        previous_id = Some(contig.id);
        let name = contig.name.as_bytes();
        let name_len = u32::try_from(name.len()).map_err(|_| JmaError::OffsetOverflow)?;
        bytes.extend_from_slice(&contig.id.to_le_bytes());
        bytes.extend_from_slice(&contig.length.to_le_bytes());
        bytes.extend_from_slice(&name_len.to_le_bytes());
        bytes.extend_from_slice(name);
    }
    Ok(bytes)
}

/// Decodes a complete contig section and verifies the declared count.
pub fn decode_contigs(bytes: &[u8], expected_count: u32) -> JmaResult<Vec<ContigMetadata>> {
    if bytes.len() < 8 {
        return Err(JmaError::CorruptSection(
            "contig section header is truncated".to_string(),
        ));
    }
    if get_u32(bytes, 0)? != SECTION_VERSION {
        return Err(JmaError::CorruptSection(
            "unsupported contig section version".to_string(),
        ));
    }
    let count = get_u32(bytes, 4)?;
    if count != expected_count {
        return Err(JmaError::CorruptSection(format!(
            "contig section count {count} does not match archive count {expected_count}"
        )));
    }
    let minimum_record_bytes = (count as usize)
        .checked_mul(16)
        .ok_or(JmaError::OffsetOverflow)?;
    let minimum_length = 8usize
        .checked_add(minimum_record_bytes)
        .ok_or(JmaError::OffsetOverflow)?;
    if bytes.len() < minimum_length {
        return Err(JmaError::CorruptSection(format!(
            "contig section is too short for {count} records"
        )));
    }
    let mut cursor = 8usize;
    let mut contigs = Vec::with_capacity(count as usize);
    let mut previous_id = None;
    for _ in 0..count {
        let id = get_u32(bytes, cursor)?;
        cursor = cursor.checked_add(4).ok_or(JmaError::OffsetOverflow)?;
        let length = get_u64(bytes, cursor)?;
        cursor = cursor.checked_add(8).ok_or(JmaError::OffsetOverflow)?;
        let name_len = get_u32(bytes, cursor)? as u64;
        cursor = cursor.checked_add(4).ok_or(JmaError::OffsetOverflow)?;
        if previous_id.is_some_and(|previous| id <= previous) {
            return Err(JmaError::CorruptSection(
                "contig identifiers are not strictly increasing".to_string(),
            ));
        }
        previous_id = Some(id);
        let name = checked_slice(bytes, cursor, name_len, "contig name")?;
        let name = std::str::from_utf8(name)
            .map_err(|_| JmaError::CorruptSection("contig name is not UTF-8".to_string()))?
            .to_string();
        cursor = cursor
            .checked_add(usize::try_from(name_len).map_err(|_| JmaError::OffsetOverflow)?)
            .ok_or(JmaError::OffsetOverflow)?;
        contigs.push(ContigMetadata { id, name, length });
    }
    if cursor != bytes.len() {
        return Err(JmaError::CorruptSection(format!(
            "contig section has {} trailing bytes",
            bytes.len() - cursor
        )));
    }
    Ok(contigs)
}

/// Finds metadata by identifier without permitting a caller to index a table
/// using an unchecked on-disk value.
pub fn find_contig(contigs: &[ContigMetadata], id: ContigId) -> JmaResult<&ContigMetadata> {
    contigs
        .iter()
        .find(|contig| contig.id == id)
        .ok_or(JmaError::UnknownContig(id))
}
