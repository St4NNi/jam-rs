//! Contig table and name-string-table encoding for JMA format 1.

use crate::jma::format::{SECTION_VERSION, checked_slice, get_u32, get_u64};
use crate::jma::{ContigId, ContigMetadata, JmaError, JmaResult};

const HEADER_BYTES: usize = 20;
const RECORD_BYTES: usize = 32;

/// Encodes a fixed-width contig table followed by one deterministic UTF-8
/// name string table. Names are never parsed on a hot seed lookup path.
pub fn encode_contigs(contigs: &[ContigMetadata]) -> JmaResult<Vec<u8>> {
    let count = u32::try_from(contigs.len()).map_err(|_| JmaError::OffsetOverflow)?;
    let mut names = Vec::new();
    let mut records = Vec::with_capacity(contigs.len());
    let mut previous_id = None;
    for contig in contigs {
        if previous_id.is_some_and(|previous| contig.id <= previous) {
            return Err(JmaError::CorruptSection(
                "contig identifiers must be strictly increasing".to_string(),
            ));
        }
        previous_id = Some(contig.id);
        let offset = u64::try_from(names.len()).map_err(|_| JmaError::OffsetOverflow)?;
        let name = contig.name.as_bytes();
        let length = u32::try_from(name.len()).map_err(|_| JmaError::OffsetOverflow)?;
        names.extend_from_slice(name);
        records.push((contig.id, contig.length, offset, length));
    }
    let mut bytes = Vec::with_capacity(
        HEADER_BYTES
            .checked_add(
                RECORD_BYTES
                    .checked_mul(records.len())
                    .ok_or(JmaError::OffsetOverflow)?,
            )
            .and_then(|size| size.checked_add(names.len()))
            .ok_or(JmaError::OffsetOverflow)?,
    );
    bytes.extend_from_slice(&SECTION_VERSION.to_le_bytes());
    bytes.extend_from_slice(&count.to_le_bytes());
    bytes.extend_from_slice(
        &u64::try_from(names.len())
            .map_err(|_| JmaError::OffsetOverflow)?
            .to_le_bytes(),
    );
    bytes.extend_from_slice(&0u32.to_le_bytes());
    for (id, length, name_offset, name_length) in records {
        bytes.extend_from_slice(&id.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&length.to_le_bytes());
        bytes.extend_from_slice(&name_offset.to_le_bytes());
        bytes.extend_from_slice(&name_length.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
    }
    bytes.extend_from_slice(&names);
    Ok(bytes)
}

/// Decodes the complete contig and name tables and verifies all offsets.
pub fn decode_contigs(bytes: &[u8], expected_count: u32) -> JmaResult<Vec<ContigMetadata>> {
    if bytes.len() < HEADER_BYTES {
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
    let names_length = get_u64(bytes, 8)?;
    let records_length = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(RECORD_BYTES))
        .ok_or(JmaError::OffsetOverflow)?;
    let names_start = HEADER_BYTES
        .checked_add(records_length)
        .ok_or(JmaError::OffsetOverflow)?;
    let names = checked_slice(bytes, names_start, names_length, "contig name table")?;
    let names_end = names_start
        .checked_add(usize::try_from(names_length).map_err(|_| JmaError::OffsetOverflow)?)
        .ok_or(JmaError::OffsetOverflow)?;
    if names_end != bytes.len() {
        return Err(JmaError::CorruptSection(format!(
            "contig name table leaves {} trailing bytes",
            bytes.len().saturating_sub(names_end)
        )));
    }
    let mut contigs = Vec::with_capacity(usize::try_from(count).unwrap_or(0));
    let mut previous_id = None;
    for index in 0..usize::try_from(count).map_err(|_| JmaError::OffsetOverflow)? {
        let cursor = HEADER_BYTES
            .checked_add(
                index
                    .checked_mul(RECORD_BYTES)
                    .ok_or(JmaError::OffsetOverflow)?,
            )
            .ok_or(JmaError::OffsetOverflow)?;
        let id = get_u32(bytes, cursor)?;
        let length = get_u64(bytes, cursor + 8)?;
        let name_offset = get_u64(bytes, cursor + 16)?;
        let name_length = u64::from(get_u32(bytes, cursor + 24)?);
        if previous_id.is_some_and(|previous| id <= previous) {
            return Err(JmaError::CorruptSection(
                "contig identifiers are not strictly increasing".to_string(),
            ));
        }
        previous_id = Some(id);
        let name_bytes = checked_slice(
            names,
            usize::try_from(name_offset).map_err(|_| JmaError::OffsetOverflow)?,
            name_length,
            "contig name",
        )?;
        let name = std::str::from_utf8(name_bytes)
            .map_err(|_| JmaError::CorruptSection("contig name is not UTF-8".to_string()))?
            .to_string();
        contigs.push(ContigMetadata { id, name, length });
    }
    Ok(contigs)
}

pub fn find_contig(contigs: &[ContigMetadata], id: ContigId) -> JmaResult<&ContigMetadata> {
    contigs
        .iter()
        .find(|contig| contig.id == id)
        .ok_or(JmaError::UnknownContig(id))
}
