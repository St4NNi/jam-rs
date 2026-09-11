use crate::jidx::{read_u32, read_u64};
use crate::shared_file::SharedFile;
use crate::shared_format::{MULTIPLE_CORE, Section, SharedError};
use crate::shared_reader::SharedReader;
use crate::shared_writer::{SharedBuildStats, publish};
use serde::Serialize;
use std::collections::BTreeSet;
use std::path::Path;

#[derive(Debug, Serialize)]
pub struct SharedPackStats {
    pub build: SharedBuildStats,
    pub source_body_sha256: [u8; 32],
    pub inline_member_groups: u64,
    pub direct_placement_members: u64,
    pub logical_memberships: u64,
    pub logical_references: u64,
}

pub fn repack_shared_index(
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
) -> Result<SharedPackStats, SharedError> {
    let source = SharedFile::open(input.as_ref(), false)?;
    source.verify_checksum()?;
    if source.header.version != 1 {
        return Err(SharedError::Invalid("packing source version"));
    }
    let mut header = source.header.clone();
    header.version = 2;
    let width = header.id_bytes();
    let mut sections: [Vec<u8>; 10] = std::array::from_fn(|_| Vec::new());
    for kind in [
        Section::Strings,
        Section::Documents,
        Section::Contigs,
        Section::Gzi,
        Section::Cores,
        Section::Occurrences,
    ] {
        let length = header.section(kind).length;
        if length > 1024 * 1024 * 1024 {
            return Err(SharedError::ResourceLimit);
        }
        sections[kind as usize] = source.section(kind, 0, length)?.to_vec();
    }
    let groups = source.section(Section::Groups, 0, header.section(Section::Groups).length)?;
    let members = source.section(Section::Members, 0, header.section(Section::Members).length)?;
    let references = source.section(
        Section::References,
        0,
        header.section(Section::References).length,
    )?;
    let occurrence_rows = header.section(Section::Occurrences).length / 24;
    let mut inline_member_groups = 0;
    let mut direct_placement_members = 0;
    let mut logical_memberships = 0u64;
    let mut logical_references = 0u64;
    let mut expected_member = 0u64;
    let mut expected_reference = 0u64;
    for group in groups.chunks_exact(32) {
        let code = read_u64(group, 0);
        let first = read_u64(group, 8);
        let count = read_u64(group, 16);
        let member_count = read_u32(group, 24);
        let tag = (code >> 62) as u8;
        let payload = code & ((1 << 62) - 1);
        if !matches!(
            (tag, payload),
            (0, 0) | (1, 0..=4095) | (2, 0..=0xffff_ffff)
        ) || member_count == 0
            || count == 0
            || count > header.occurrence_count
            || read_u32(group, 28) != 0
            || first != expected_member
        {
            return Err(SharedError::Invalid("packing group"));
        }
        let end = first
            .checked_add(u64::from(member_count))
            .filter(|end| *end <= members.len() as u64 / 24)
            .ok_or(SharedError::Invalid("packing member extent"))?;
        let inline = member_count == 1;
        let mut locator = narrow(
            sections[Section::Members as usize].len() as u64 / header.row_bytes(Section::Members),
        )?;
        let mut value = member_count;
        let mut previous = None;
        let mut total = 0u64;
        for member in members[first as usize * 24..end as usize * 24].chunks_exact(24) {
            let id = read_u32(member, 0);
            let first = read_u64(member, 8);
            let count = read_u64(member, 16);
            if id >= header.document_count
                || previous.is_some_and(|before| before >= id)
                || read_u32(member, 4) != 0
                || count == 0
                || first != expected_reference
            {
                return Err(SharedError::Invalid("packing member"));
            }
            let end = first
                .checked_add(count)
                .filter(|end| *end <= references.len() as u64 / 8)
                .ok_or(SharedError::Invalid("packing reference extent"))?;
            let stored_count = narrow(count)?;
            let mut stored_first = narrow(sections[Section::References as usize].len() as u64 / 4)?;
            for raw in references[first as usize * 8..end as usize * 8].chunks_exact(8) {
                let ordinal = read_u64(raw, 0);
                if ordinal >= occurrence_rows {
                    return Err(SharedError::Invalid("packing occurrence reference"));
                }
                let ordinal = narrow(ordinal)?;
                if count == 1 {
                    stored_first = ordinal;
                    direct_placement_members += 1;
                } else {
                    sections[Section::References as usize]
                        .extend_from_slice(&ordinal.to_le_bytes());
                }
            }
            if inline {
                locator = stored_first;
                value = id;
                inline_member_groups += 1;
            } else {
                let target = &mut sections[Section::Members as usize];
                target.extend_from_slice(&id.to_le_bytes()[..width]);
                target.extend_from_slice(&stored_first.to_le_bytes());
                target.extend_from_slice(&stored_count.to_le_bytes());
            }
            total = total.checked_add(count).ok_or(SharedError::ResourceLimit)?;
            previous = Some(id);
            expected_reference = end;
        }
        if total != count {
            return Err(SharedError::Invalid("packing group count"));
        }
        let target = &mut sections[Section::Groups as usize];
        target.extend_from_slice(&(payload as u32).to_le_bytes());
        target.push(tag | if inline { 4 } else { 0 });
        target.extend_from_slice(&value.to_le_bytes()[..width]);
        target.extend_from_slice(&locator.to_le_bytes());
        target.extend_from_slice(&narrow(count)?.to_le_bytes());
        logical_memberships += u64::from(member_count);
        logical_references = logical_references
            .checked_add(count)
            .ok_or(SharedError::ResourceLimit)?;
        expected_member = end;
        if sections.iter().map(Vec::capacity).sum::<usize>() > 1024 * 1024 * 1024 {
            return Err(SharedError::ResourceLimit);
        }
    }
    if expected_member != members.len() as u64 / 24
        || expected_reference != references.len() as u64 / 8
    {
        return Err(SharedError::Invalid("packing unreferenced payload"));
    }
    let singletons = sections[Section::Cores as usize]
        .chunks_exact(24)
        .filter(|row| read_u32(row, 0) & MULTIPLE_CORE == 0)
        .count() as u64;
    let metadata = SharedReader::open(input)?;
    let mut seen = BTreeSet::new();
    let mut bgzf_bytes = 0u64;
    for id in 0..header.document_count {
        let document = metadata
            .metagenome(id)?
            .ok_or(SharedError::Invalid("packing document"))?;
        if seen.insert(document.bgzf_uri) {
            bgzf_bytes = bgzf_bytes
                .checked_add(document.bgzf_bytes)
                .ok_or(SharedError::ResourceLimit)?;
        }
    }
    source.verify_unchanged()?;
    Ok(SharedPackStats {
        build: publish(header, output.as_ref(), sections, bgzf_bytes, singletons)?,
        source_body_sha256: source.header.body_sha256,
        inline_member_groups,
        direct_placement_members,
        logical_memberships,
        logical_references,
    })
}

fn narrow(value: u64) -> Result<u32, SharedError> {
    u32::try_from(value).map_err(|_| SharedError::ResourceLimit)
}

#[cfg(test)]
mod tests {
    #[test]
    fn packed_locators_reject_overflow() {
        assert_eq!(super::narrow(u64::from(u32::MAX)).unwrap(), u32::MAX);
        assert!(super::narrow(u64::from(u32::MAX) + 1).is_err());
        assert!(super::narrow(u64::MAX).is_err());
    }
}
