use crate::jidx::{read_u32, read_u64};
use crate::shared_file::SharedFile;
use crate::shared_format::{
    CORE_PREFIX_BOUNDARIES, CORE_ROW_BYTES, MULTIPLE_CORE, Section, SharedError,
};
use crate::shared_reader::{CoreKind, CoreRow, SharedReader};
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

#[derive(Debug, Serialize)]
pub struct SharedCorePackStats {
    pub build: SharedBuildStats,
    pub source_version: u16,
    pub source_body_sha256: [u8; 32],
    pub source_identity: Option<[u64; 7]>,
    pub core_payload_bytes: u8,
    pub hot_core_bytes: u64,
    pub cold_core_bytes: u64,
    pub core_prefix_bytes: u64,
}

pub fn repack_shared_index(
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
) -> Result<SharedPackStats, SharedError> {
    let input = input.as_ref();
    let source = SharedFile::open(input, false)?;
    source.verify_checksum()?;
    if source.header.version != 1 {
        return Err(SharedError::Invalid("packing source version"));
    }
    let mut header = source.header.clone();
    header.version = 2;
    let width = header.id_bytes();
    let mut sections: [Vec<u8>; 12] = std::array::from_fn(|_| Vec::new());
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
    for group in groups.as_chunks::<32>().0 {
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
        for member in members[first as usize * 24..end as usize * 24]
            .as_chunks::<24>()
            .0
        {
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
            for raw in references[first as usize * 8..end as usize * 8]
                .as_chunks::<8>()
                .0
            {
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
        .as_chunks::<24>()
        .0
        .iter()
        .filter(|row| read_u32(*row, 0) & MULTIPLE_CORE == 0)
        .count() as u64;
    let bgzf_bytes = bgzf_bytes_once(input, header.document_count)?;
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

pub fn repack_shared_cores(
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
) -> Result<SharedCorePackStats, SharedError> {
    let input = input.as_ref();
    let source = SharedFile::open(input, false)?;
    source.verify_checksum()?;
    if source.header.version != 2 {
        return Err(SharedError::Invalid("core packing source version"));
    }
    let source_identity = source.identity();
    let mut header = source.header.clone();
    let group_count = header.section(Section::Groups).length / header.row_bytes(Section::Groups);
    let cores = source.section(Section::Cores, 0, header.section(Section::Cores).length)?;
    let mut previous = None;
    let mut wide = false;
    let mut singletons = 0u64;
    for bytes in cores.as_chunks::<{ CORE_ROW_BYTES as usize }>().0 {
        let row = CoreRow::decode(bytes, header.contig_count, group_count)?;
        if previous.is_some_and(|core| core >= row.core) {
            return Err(SharedError::Invalid("core packing order"));
        }
        previous = Some(row.core);
        match row.kind {
            CoreKind::Singleton { .. } => {
                singletons += 1;
                wide |= core_payload_needs_wide(row.kind);
            }
            CoreKind::Repeated { .. } => wide |= core_payload_needs_wide(row.kind),
        }
    }
    let core_payload_bytes = if wide { 21 } else { 13 };
    let core_count = usize::try_from(header.core_count).map_err(|_| SharedError::ResourceLimit)?;
    let hot_bytes = core_count
        .checked_mul(4)
        .ok_or(SharedError::ResourceLimit)?;
    let cold_bytes = core_count
        .checked_mul(core_payload_bytes as usize)
        .ok_or(SharedError::ResourceLimit)?;
    let prefix_bytes = CORE_PREFIX_BOUNDARIES
        .checked_mul(4)
        .ok_or(SharedError::ResourceLimit)?;
    let mut sections: [Vec<u8>; 12] = std::array::from_fn(|_| Vec::new());
    for kind in [
        Section::Strings,
        Section::Documents,
        Section::Contigs,
        Section::Gzi,
        Section::Groups,
        Section::Members,
        Section::References,
        Section::Occurrences,
    ] {
        let length = header.section(kind).length;
        if length > 1024 * 1024 * 1024 {
            return Err(SharedError::ResourceLimit);
        }
        sections[kind as usize] = source.section(kind, 0, length)?.to_vec();
        check_capacity(&sections)?;
    }
    let added_capacity = hot_bytes
        .checked_add(cold_bytes)
        .and_then(|bytes| bytes.checked_add(prefix_bytes))
        .ok_or(SharedError::ResourceLimit)?;
    check_additional_capacity(&sections, added_capacity)?;
    reserve(&mut sections[Section::Cores as usize], hot_bytes)?;
    reserve(&mut sections[Section::CorePayloads as usize], cold_bytes)?;
    reserve(&mut sections[Section::CorePrefixes as usize], prefix_bytes)?;
    check_capacity(&sections)?;
    let count = u32::try_from(header.core_count).map_err(|_| SharedError::ResourceLimit)?;
    let mut next_prefix = 0usize;
    for (ordinal, bytes) in cores
        .as_chunks::<{ CORE_ROW_BYTES as usize }>()
        .0
        .iter()
        .enumerate()
    {
        let row = CoreRow::decode(bytes, header.contig_count, group_count)?;
        let prefix = (row.core >> 14) as usize;
        append_prefix_boundaries(
            &mut sections[Section::CorePrefixes as usize],
            &mut next_prefix,
            prefix,
            ordinal as u32,
        );
        let tagged = row.core
            | if matches!(row.kind, CoreKind::Repeated { .. }) {
                MULTIPLE_CORE
            } else {
                0
            };
        sections[Section::Cores as usize].extend_from_slice(&tagged.to_le_bytes());
        encode_core_payload(
            &mut sections[Section::CorePayloads as usize],
            row.kind,
            wide,
        )?;
    }
    while next_prefix < CORE_PREFIX_BOUNDARIES {
        sections[Section::CorePrefixes as usize].extend_from_slice(&count.to_le_bytes());
        next_prefix += 1;
    }
    if sections[Section::Cores as usize].len() != hot_bytes
        || sections[Section::CorePayloads as usize].len() != cold_bytes
        || sections[Section::CorePrefixes as usize].len() != prefix_bytes
    {
        return Err(SharedError::Invalid("core packing lengths"));
    }
    check_capacity(&sections)?;
    let bgzf_bytes = bgzf_bytes_once(input, header.document_count)?;
    source.verify_unchanged()?;
    header.version = 3;
    header.core_payload_bytes = core_payload_bytes;
    Ok(SharedCorePackStats {
        build: publish(header, output.as_ref(), sections, bgzf_bytes, singletons)?,
        source_version: source.header.version,
        source_body_sha256: source.header.body_sha256,
        source_identity,
        core_payload_bytes,
        hot_core_bytes: hot_bytes as u64,
        cold_core_bytes: cold_bytes as u64,
        core_prefix_bytes: prefix_bytes as u64,
    })
}

fn core_payload_needs_wide(kind: CoreKind) -> bool {
    match kind {
        CoreKind::Singleton { position, .. } => u32::try_from(position).is_err(),
        CoreKind::Repeated {
            first_group,
            occurrence_count,
            ..
        } => u32::try_from(first_group).is_err() || u32::try_from(occurrence_count).is_err(),
    }
}

fn append_prefix_boundaries(
    target: &mut Vec<u8>,
    next_prefix: &mut usize,
    through: usize,
    ordinal: u32,
) {
    while *next_prefix <= through {
        target.extend_from_slice(&ordinal.to_le_bytes());
        *next_prefix += 1;
    }
}

fn encode_core_payload(
    target: &mut Vec<u8>,
    kind: CoreKind,
    wide: bool,
) -> Result<(), SharedError> {
    match (wide, kind) {
        (
            false,
            CoreKind::Singleton {
                context,
                contig_id,
                flags,
                position,
            },
        ) => {
            target.extend_from_slice(&context.to_le_bytes());
            target.extend_from_slice(&contig_id.to_le_bytes());
            target.extend_from_slice(&narrow(position)?.to_le_bytes());
            target.push(flags as u8);
        }
        (
            false,
            CoreKind::Repeated {
                first_group,
                group_count,
                occurrence_count,
            },
        ) => {
            target.extend_from_slice(&narrow(first_group)?.to_le_bytes());
            target.extend_from_slice(&group_count.to_le_bytes());
            target.extend_from_slice(&narrow(occurrence_count)?.to_le_bytes());
            target.push(0);
        }
        (
            true,
            CoreKind::Singleton {
                context,
                contig_id,
                flags,
                position,
            },
        ) => {
            target.extend_from_slice(&context.to_le_bytes());
            target.extend_from_slice(&contig_id.to_le_bytes());
            target.extend_from_slice(&position.to_le_bytes());
            target.extend_from_slice(&0u32.to_le_bytes());
            target.push(flags as u8);
        }
        (
            true,
            CoreKind::Repeated {
                first_group,
                group_count,
                occurrence_count,
            },
        ) => {
            target.extend_from_slice(&first_group.to_le_bytes());
            target.extend_from_slice(&group_count.to_le_bytes());
            target.extend_from_slice(&occurrence_count.to_le_bytes());
            target.push(0);
        }
    }
    Ok(())
}

fn bgzf_bytes_once(input: &Path, document_count: u32) -> Result<u64, SharedError> {
    let metadata = SharedReader::open(input)?;
    let mut seen = BTreeSet::new();
    let mut bgzf_bytes = 0u64;
    for id in 0..document_count {
        let document = metadata
            .metagenome(id)?
            .ok_or(SharedError::Invalid("packing document"))?;
        if seen.insert(document.bgzf_uri) {
            bgzf_bytes = bgzf_bytes
                .checked_add(document.bgzf_bytes)
                .ok_or(SharedError::ResourceLimit)?;
        }
    }
    Ok(bgzf_bytes)
}

fn reserve(target: &mut Vec<u8>, bytes: usize) -> Result<(), SharedError> {
    target
        .try_reserve_exact(bytes)
        .map_err(|_| SharedError::ResourceLimit)
}

fn check_capacity(sections: &[Vec<u8>; 12]) -> Result<(), SharedError> {
    if retained_capacity(sections)? > 1024 * 1024 * 1024 {
        return Err(SharedError::ResourceLimit);
    }
    Ok(())
}

fn check_additional_capacity(
    sections: &[Vec<u8>; 12],
    additional: usize,
) -> Result<(), SharedError> {
    if retained_capacity(sections)?
        .checked_add(additional)
        .is_none_or(|bytes| bytes > 1024 * 1024 * 1024)
    {
        return Err(SharedError::ResourceLimit);
    }
    Ok(())
}

fn retained_capacity(sections: &[Vec<u8>; 12]) -> Result<usize, SharedError> {
    sections.iter().try_fold(0usize, |total, section| {
        total
            .checked_add(section.capacity())
            .ok_or(SharedError::ResourceLimit)
    })
}

fn narrow(value: u64) -> Result<u32, SharedError> {
    u32::try_from(value).map_err(|_| SharedError::ResourceLimit)
}

#[cfg(test)]
mod tests {
    use crate::shared_reader::CoreKind;

    #[test]
    fn packed_locators_reject_overflow() {
        assert_eq!(super::narrow(u64::from(u32::MAX)).unwrap(), u32::MAX);
        assert!(super::narrow(u64::from(u32::MAX) + 1).is_err());
        assert!(super::narrow(u64::MAX).is_err());
    }

    #[test]
    fn core_payloads_use_exact_narrow_and_wide_layouts() {
        let singleton = CoreKind::Singleton {
            context: 0x1122_3344,
            contig_id: 0x5566_7788,
            flags: 3,
            position: 0x99aa_bbcc,
        };
        let repeated = CoreKind::Repeated {
            first_group: 0x1122_3344,
            group_count: 0x5566_7788,
            occurrence_count: 0x99aa_bbcc,
        };
        let mut bytes = Vec::new();
        super::encode_core_payload(&mut bytes, singleton, false).unwrap();
        super::encode_core_payload(&mut bytes, repeated, false).unwrap();
        assert_eq!(
            bytes,
            [
                0x44, 0x33, 0x22, 0x11, 0x88, 0x77, 0x66, 0x55, 0xcc, 0xbb, 0xaa, 0x99, 3, 0x44,
                0x33, 0x22, 0x11, 0x88, 0x77, 0x66, 0x55, 0xcc, 0xbb, 0xaa, 0x99, 0,
            ]
        );
        bytes.clear();
        super::encode_core_payload(&mut bytes, singleton, true).unwrap();
        super::encode_core_payload(&mut bytes, repeated, true).unwrap();
        assert_eq!(bytes.len(), 42);
        assert_eq!(&bytes[8..16], &0x99aa_bbccu64.to_le_bytes());
        assert_eq!(&bytes[16..21], &[0, 0, 0, 0, 3]);
        assert_eq!(&bytes[21..29], &0x1122_3344u64.to_le_bytes());
        assert_eq!(&bytes[33..41], &0x99aa_bbccu64.to_le_bytes());
        assert_eq!(bytes[41], 0);
    }

    #[test]
    fn narrow_core_payloads_reject_wide_values() {
        let mut bytes = Vec::new();
        assert!(
            super::encode_core_payload(
                &mut bytes,
                CoreKind::Singleton {
                    context: 0,
                    contig_id: 0,
                    flags: 0,
                    position: u64::from(u32::MAX) + 1,
                },
                false,
            )
            .is_err()
        );
        assert!(
            super::encode_core_payload(
                &mut bytes,
                CoreKind::Repeated {
                    first_group: u64::from(u32::MAX) + 1,
                    group_count: 1,
                    occurrence_count: 2,
                },
                false,
            )
            .is_err()
        );
        assert!(
            super::encode_core_payload(
                &mut bytes,
                CoreKind::Repeated {
                    first_group: 0,
                    group_count: 1,
                    occurrence_count: u64::from(u32::MAX) + 1,
                },
                false,
            )
            .is_err()
        );
    }

    #[test]
    fn core_payload_mode_checks_synthetic_overflow_boundaries() {
        assert!(!super::core_payload_needs_wide(CoreKind::Singleton {
            context: 0,
            contig_id: 0,
            flags: 0,
            position: u64::from(u32::MAX),
        }));
        assert!(super::core_payload_needs_wide(CoreKind::Singleton {
            context: 0,
            contig_id: 0,
            flags: 0,
            position: u64::from(u32::MAX) + 1,
        }));
        assert!(super::core_payload_needs_wide(CoreKind::Repeated {
            first_group: u64::from(u32::MAX) + 1,
            group_count: 1,
            occurrence_count: 2,
        }));
        assert!(super::core_payload_needs_wide(CoreKind::Repeated {
            first_group: 0,
            group_count: 1,
            occurrence_count: u64::from(u32::MAX) + 1,
        }));
    }

    #[test]
    fn prefix_boundaries_cover_empty_and_edge_ranges() {
        let mut bytes = Vec::new();
        let mut next = 0;
        super::append_prefix_boundaries(&mut bytes, &mut next, 0, 0);
        super::append_prefix_boundaries(&mut bytes, &mut next, 2, 1);
        super::append_prefix_boundaries(&mut bytes, &mut next, 65_535, 2);
        assert_eq!(next, 65_536);
        assert_eq!(bytes.len(), 65_536 * 4);
        assert_eq!(u32::from_le_bytes(bytes[0..4].try_into().unwrap()), 0);
        assert_eq!(u32::from_le_bytes(bytes[4..8].try_into().unwrap()), 1);
        assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 1);
        assert_eq!(
            u32::from_le_bytes(bytes[65_535 * 4..65_536 * 4].try_into().unwrap()),
            2
        );
        super::append_prefix_boundaries(&mut bytes, &mut next, 65_536, 3);
        assert_eq!(bytes.len(), super::CORE_PREFIX_BOUNDARIES * 4);
        assert_eq!(
            u32::from_le_bytes(bytes[65_536 * 4..].try_into().unwrap()),
            3
        );
    }
}
