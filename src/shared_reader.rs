use crate::jidx::sha256;
pub use crate::shared_file::FileReadStats;
use crate::shared_file::SharedFile;
use crate::shared_format::{
    CORE_MASK, CORE_ROW_BYTES, GROUP_ROW_BYTES, MEMBER_ROW_BYTES, MULTIPLE_CORE, Section,
    SharedError, read_u32, read_u64,
};
use crate::shared_seed::SharedKey;
use serde::Serialize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_DECODED_RESULT_BYTES: usize = 64 * 1024 * 1024;

pub struct SharedReader {
    file: SharedFile,
    core_inspections: AtomicU64,
    group_inspections: AtomicU64,
    member_inspections: AtomicU64,
    references_decoded: AtomicU64,
    positions_decoded: AtomicU64,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct SharedReadStats {
    pub file: FileReadStats,
    pub core_descriptor_inspections: u64,
    pub group_descriptor_inspections: u64,
    pub member_descriptor_inspections: u64,
    pub references_decoded: u64,
    pub physical_positions_decoded: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SharedGroup {
    identity: HandleIdentity,
    key: SharedKey,
    core_ordinal: u64,
    location: GroupLocation,
    member_count: u32,
    occurrence_count: u64,
}

impl SharedGroup {
    pub fn key(self) -> SharedKey {
        self.key
    }

    pub fn member_count(self) -> u32 {
        self.member_count
    }

    pub fn occurrence_count(self) -> u64 {
        self.occurrence_count
    }

    pub fn core_ordinal(self) -> u64 {
        self.core_ordinal
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct HandleIdentity {
    file: Option<[u64; 7]>,
    body_sha256: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum GroupLocation {
    Singleton {
        core_ordinal: u64,
    },
    Repeated {
        group_ordinal: u64,
        first_member: u64,
    },
}

impl SharedReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, SharedError> {
        Ok(Self {
            file: SharedFile::open(path)?,
            core_inspections: AtomicU64::new(0),
            group_inspections: AtomicU64::new(0),
            member_inspections: AtomicU64::new(0),
            references_decoded: AtomicU64::new(0),
            positions_decoded: AtomicU64::new(0),
        })
    }

    pub fn window(&self) -> u16 {
        self.file.header.window
    }

    pub fn core_count(&self) -> u64 {
        self.file.header.core_count
    }

    pub fn occurrence_count(&self) -> u64 {
        self.file.header.occurrence_count
    }

    pub fn document_count(&self) -> u32 {
        self.file.header.document_count
    }

    pub fn contig_count(&self) -> u32 {
        self.file.header.contig_count
    }

    pub fn source_bases(&self) -> u64 {
        self.file.header.source_bases
    }

    pub fn manifest_sha256(&self) -> [u8; 32] {
        self.file.header.manifest_sha256
    }

    pub fn body_sha256(&self) -> [u8; 32] {
        self.file.header.body_sha256
    }

    pub fn header_sha256(&self) -> Result<[u8; 32], SharedError> {
        Ok(sha256(&self.file.header.encode()?))
    }

    pub fn cache_identity(&self) -> Result<Option<[u64; 7]>, SharedError> {
        self.file.verify_unchanged()?;
        Ok(self.file.identity())
    }

    pub fn verify_checksum(&self) -> Result<(), SharedError> {
        self.file.verify_checksum()
    }

    pub fn identity(&self) -> Option<[u64; 7]> {
        self.file.identity()
    }

    pub fn stats(&self) -> SharedReadStats {
        SharedReadStats {
            file: self.file.stats(),
            core_descriptor_inspections: self.core_inspections.load(Ordering::Relaxed),
            group_descriptor_inspections: self.group_inspections.load(Ordering::Relaxed),
            member_descriptor_inspections: self.member_inspections.load(Ordering::Relaxed),
            references_decoded: self.references_decoded.load(Ordering::Relaxed),
            physical_positions_decoded: self.positions_decoded.load(Ordering::Relaxed),
        }
    }

    pub fn find(&self, key: SharedKey) -> Result<Option<SharedGroup>, SharedError> {
        self.begin_operation()?;
        let group = self.find_unchecked(key)?;
        self.file.verify_unchanged()?;
        Ok(group)
    }

    pub fn find_many(&self, keys: &[SharedKey]) -> Result<Vec<Option<SharedGroup>>, SharedError> {
        self.begin_operation()?;
        admit_result(keys.len(), size_of::<Option<SharedGroup>>())?;
        let mut groups = Vec::new();
        groups
            .try_reserve_exact(keys.len())
            .map_err(|_| SharedError::ResourceLimit)?;
        for &key in keys {
            groups.push(self.find_unchecked(key)?);
        }
        self.file.verify_unchanged()?;
        Ok(groups)
    }

    pub fn core_group(
        &self,
        core_ordinal: u64,
        expected_core: u32,
    ) -> Result<SharedGroup, SharedError> {
        self.begin_operation()?;
        if expected_core > CORE_MASK || core_ordinal >= self.file.header.core_count {
            return Err(SharedError::Invalid("core handle"));
        }
        let row = self.core_row(core_ordinal)?;
        if row.core != expected_core {
            return Err(SharedError::Invalid("core handle"));
        }
        let group = self
            .group_from_core(core_ordinal, row, SharedKey::core(expected_core))?
            .ok_or(SharedError::Invalid("core group"))?;
        self.file.verify_unchanged()?;
        Ok(group)
    }

    fn begin_operation(&self) -> Result<(), SharedError> {
        self.file.verify_unchanged()
    }

    fn handle_identity(&self) -> HandleIdentity {
        HandleIdentity {
            file: self.file.identity(),
            body_sha256: self.file.header.body_sha256,
        }
    }

    fn validate_group(&self, group: SharedGroup) -> Result<(), SharedError> {
        self.begin_operation()?;
        if group.identity != self.handle_identity() {
            return Err(SharedError::Invalid("group handle"));
        }
        Ok(())
    }

    fn find_unchecked(&self, key: SharedKey) -> Result<Option<SharedGroup>, SharedError> {
        key.context_code()
            .ok_or(SharedError::Invalid("shared key"))?;
        let mut low = 0;
        let mut high = self.file.header.core_count;
        while low < high {
            let middle = low + (high - low) / 2;
            let row = self.core_row(middle)?;
            match row.core.cmp(&key.core) {
                std::cmp::Ordering::Less => low = middle + 1,
                std::cmp::Ordering::Greater => high = middle,
                std::cmp::Ordering::Equal => return self.group_from_core(middle, row, key),
            }
        }
        Ok(None)
    }

    fn core_row(&self, ordinal: u64) -> Result<CoreRow, SharedError> {
        let bytes = self.file.record(Section::Cores, ordinal, CORE_ROW_BYTES)?;
        self.core_inspections.fetch_add(1, Ordering::Relaxed);
        CoreRow::decode(
            bytes,
            self.file.header.contig_count,
            self.file.header.section(Section::Groups).length / GROUP_ROW_BYTES,
        )
    }

    fn group_from_core(
        &self,
        core_ordinal: u64,
        row: CoreRow,
        key: SharedKey,
    ) -> Result<Option<SharedGroup>, SharedError> {
        if let CoreKind::Singleton { context, flags, .. } = row.kind {
            let matches = match key.length {
                15 => true,
                21 => flags & 2 != 0 && key.context == context >> 20,
                31 => flags & 4 != 0 && key.context == context,
                _ => false,
            };
            return Ok(matches.then(|| SharedGroup {
                identity: self.handle_identity(),
                key,
                core_ordinal,
                location: GroupLocation::Singleton { core_ordinal },
                member_count: 1,
                occurrence_count: 1,
            }));
        }
        let CoreKind::Repeated {
            first_group,
            group_count,
            occurrence_count,
        } = row.kind
        else {
            unreachable!();
        };
        let wanted = key.context_code().unwrap();
        let mut low = 0u64;
        let mut high = u64::from(group_count);
        while low < high {
            let middle = low + (high - low) / 2;
            let ordinal = first_group
                .checked_add(middle)
                .ok_or(SharedError::Invalid("group ordinal"))?;
            let group = self.group_row(ordinal)?;
            match group.context_code.cmp(&wanted) {
                std::cmp::Ordering::Less => low = middle + 1,
                std::cmp::Ordering::Greater => high = middle,
                std::cmp::Ordering::Equal => {
                    if wanted == 0 && group.occurrence_count != occurrence_count {
                        return Err(SharedError::Invalid("core occurrence count"));
                    }
                    return Ok(Some(SharedGroup {
                        identity: self.handle_identity(),
                        key,
                        core_ordinal,
                        location: GroupLocation::Repeated {
                            group_ordinal: ordinal,
                            first_member: group.first_member,
                        },
                        member_count: group.member_count,
                        occurrence_count: group.occurrence_count,
                    }));
                }
            }
        }
        Ok(None)
    }

    fn group_row(&self, ordinal: u64) -> Result<GroupRow, SharedError> {
        let bytes = self
            .file
            .record(Section::Groups, ordinal, GROUP_ROW_BYTES)?;
        self.group_inspections.fetch_add(1, Ordering::Relaxed);
        let context_code = read_u64(bytes, 0);
        let first_member = read_u64(bytes, 8);
        let occurrence_count = read_u64(bytes, 16);
        let member_count = read_u32(bytes, 24);
        if !valid_context_code(context_code)
            || member_count == 0
            || occurrence_count == 0
            || read_u32(bytes, 28) != 0
            || first_member
                .checked_add(u64::from(member_count))
                .is_none_or(|end| {
                    end > self.file.header.section(Section::Members).length / MEMBER_ROW_BYTES
                })
        {
            return Err(SharedError::Invalid("group row"));
        }
        Ok(GroupRow {
            context_code,
            first_member,
            occurrence_count,
            member_count,
        })
    }
}

#[derive(Clone, Copy)]
struct CoreRow {
    core: u32,
    kind: CoreKind,
}

#[derive(Clone, Copy)]
enum CoreKind {
    Singleton {
        context: u32,
        contig_id: u32,
        flags: u32,
        position: u64,
    },
    Repeated {
        first_group: u64,
        group_count: u32,
        occurrence_count: u64,
    },
}

impl CoreRow {
    fn decode(bytes: &[u8], contig_count: u32, total_groups: u64) -> Result<Self, SharedError> {
        let word = read_u32(bytes, 0);
        if word & !(MULTIPLE_CORE | CORE_MASK) != 0 {
            return Err(SharedError::Invalid("core row"));
        }
        let core = word & CORE_MASK;
        let kind = if word & MULTIPLE_CORE == 0 {
            let context = read_u32(bytes, 4);
            let flags = read_u32(bytes, 12);
            if flags & !7 != 0
                || flags & 4 != 0 && flags & 2 == 0
                || flags & 2 == 0 && context != 0
                || flags & 4 == 0 && context & ((1 << 20) - 1) != 0
                || read_u32(bytes, 8) >= contig_count
            {
                return Err(SharedError::Invalid("singleton core"));
            }
            CoreKind::Singleton {
                context,
                contig_id: read_u32(bytes, 8),
                flags,
                position: read_u64(bytes, 16),
            }
        } else {
            let occurrence_count =
                u64::from(read_u32(bytes, 8)) | (u64::from(read_u32(bytes, 12)) << 32);
            let group_count = read_u32(bytes, 4);
            let first_group = read_u64(bytes, 16);
            if group_count == 0
                || occurrence_count < 2
                || first_group
                    .checked_add(u64::from(group_count))
                    .is_none_or(|end| end > total_groups)
            {
                return Err(SharedError::Invalid("repeated core"));
            }
            CoreKind::Repeated {
                first_group,
                group_count,
                occurrence_count,
            }
        };
        Ok(Self { core, kind })
    }
}

struct GroupRow {
    context_code: u64,
    first_member: u64,
    occurrence_count: u64,
    member_count: u32,
}

fn valid_context_code(code: u64) -> bool {
    match code >> 62 {
        0 => code == 0,
        1 => code & ((1 << 62) - 1) < 1 << 12,
        2 => code & ((1 << 62) - 1) <= u64::from(u32::MAX),
        _ => false,
    }
}

fn admit_result(count: usize, row_bytes: usize) -> Result<(), SharedError> {
    if count
        .checked_mul(row_bytes)
        .is_none_or(|bytes| bytes > MAX_DECODED_RESULT_BYTES)
    {
        return Err(SharedError::ResourceLimit);
    }
    Ok(())
}
