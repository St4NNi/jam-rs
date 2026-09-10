use crate::jidx::{
    CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_RECORD_SIZE, DocumentRecord, StringRef, sha256,
};
use crate::jidx_reader::{Contig, Metagenome, SeedOccurrence};
pub use crate::shared_file::FileReadStats;
use crate::shared_file::SharedFile;
use crate::shared_format::{
    CORE_MASK, CORE_ROW_BYTES, GROUP_ROW_BYTES, MEMBER_ROW_BYTES, MULTIPLE_CORE,
    OCCURRENCE_ROW_BYTES, Section, SharedError, read_u32, read_u64,
};
use crate::shared_seed::{SharedKey, SharedSeed};
use serde::Serialize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

const MAX_DECODED_RESULT_BYTES: usize = 64 * 1024 * 1024;

pub struct SharedReader {
    file: SharedFile,
    observed: bool,
    core_inspections: AtomicU64,
    group_inspections: AtomicU64,
    member_inspections: AtomicU64,
    references_decoded: AtomicU64,
    positions_decoded: AtomicU64,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct SharedReadStats {
    pub observed: bool,
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

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct SharedOccurrenceStorage {
    body_sha256: [u8; 32],
    kind: u8,
    first: u64,
    count: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SharedMember {
    identity: HandleIdentity,
    group: GroupLocation,
    pub metagenome_id: u32,
    first_reference: u64,
    occurrence_count: u64,
}

impl SharedMember {
    pub fn occurrence_count(self) -> u64 {
        self.occurrence_count
    }

    pub fn occurrence_storage_identity(self) -> SharedOccurrenceStorage {
        let (kind, first) = match self.group {
            GroupLocation::Singleton { core_ordinal } => (2, core_ordinal),
            GroupLocation::Repeated { .. } => (3, self.first_reference),
        };
        SharedOccurrenceStorage {
            body_sha256: self.identity.body_sha256,
            kind,
            first,
            count: self.occurrence_count,
        }
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
        Self::open_inner(path, false)
    }

    pub fn open_observed(path: impl AsRef<Path>) -> Result<Self, SharedError> {
        Self::open_inner(path, true)
    }

    fn open_inner(path: impl AsRef<Path>, observed: bool) -> Result<Self, SharedError> {
        let file = SharedFile::open(path, observed)?;
        Ok(Self {
            file,
            observed,
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
            observed: self.observed,
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
        admit_result(groups.capacity(), size_of::<Option<SharedGroup>>())?;
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
        self.group_at(core_ordinal, SharedKey::core(expected_core))?
            .ok_or(SharedError::Invalid("core group"))
    }

    pub fn group_at(
        &self,
        core_ordinal: u64,
        key: SharedKey,
    ) -> Result<Option<SharedGroup>, SharedError> {
        self.begin_operation()?;
        key.context_code()
            .ok_or(SharedError::Invalid("shared key"))?;
        if core_ordinal >= self.file.header.core_count {
            return Err(SharedError::Invalid("core handle"));
        }
        let row = self.core_row(core_ordinal)?;
        if row.core != key.core {
            return Err(SharedError::Invalid("core handle"));
        }
        let group = self.group_from_core(core_ordinal, row, key)?;
        self.file.verify_unchanged()?;
        Ok(group)
    }

    pub fn members(&self, group: SharedGroup) -> Result<Vec<SharedMember>, SharedError> {
        self.validate_group(group)?;
        admit_result(group.member_count as usize, size_of::<SharedMember>())?;
        let members = match group.location {
            GroupLocation::Singleton { .. } => vec![self.singleton_member(group)?],
            GroupLocation::Repeated { first_member, .. } => {
                let mut members = Vec::new();
                members
                    .try_reserve_exact(group.member_count as usize)
                    .map_err(|_| SharedError::ResourceLimit)?;
                admit_result(members.capacity(), size_of::<SharedMember>())?;
                let mut previous = None;
                let mut occurrences = 0u64;
                for offset in 0..u64::from(group.member_count) {
                    let member = self.member_row(group, first_member + offset)?;
                    if previous.is_some_and(|id| id >= member.metagenome_id) {
                        return Err(SharedError::Invalid("member order"));
                    }
                    previous = Some(member.metagenome_id);
                    occurrences = occurrences
                        .checked_add(member.occurrence_count)
                        .ok_or(SharedError::Invalid("member occurrence count"))?;
                    members.push(member);
                }
                if occurrences != group.occurrence_count {
                    return Err(SharedError::Invalid("member occurrence count"));
                }
                members
            }
        };
        self.file.verify_unchanged()?;
        Ok(members)
    }

    pub fn member(
        &self,
        group: SharedGroup,
        metagenome_id: u32,
    ) -> Result<Option<SharedMember>, SharedError> {
        self.validate_group(group)?;
        if metagenome_id >= self.file.header.document_count {
            return Ok(None);
        }
        let GroupLocation::Repeated { first_member, .. } = group.location else {
            let member = self.singleton_member(group)?;
            let result = (member.metagenome_id == metagenome_id).then_some(member);
            self.file.verify_unchanged()?;
            return Ok(result);
        };
        let mut low = 0u64;
        let mut high = u64::from(group.member_count);
        while low < high {
            let middle = low + (high - low) / 2;
            let member = self.member_row(group, first_member + middle)?;
            match member.metagenome_id.cmp(&metagenome_id) {
                std::cmp::Ordering::Less => low = middle + 1,
                std::cmp::Ordering::Greater => high = middle,
                std::cmp::Ordering::Equal => {
                    self.file.verify_unchanged()?;
                    return Ok(Some(member));
                }
            }
        }
        self.file.verify_unchanged()?;
        Ok(None)
    }

    pub fn member_occurrences(
        &self,
        group: SharedGroup,
        member: SharedMember,
    ) -> Result<Vec<SeedOccurrence>, SharedError> {
        self.validate_member(group, member)?;
        let count =
            usize::try_from(member.occurrence_count).map_err(|_| SharedError::ResourceLimit)?;
        admit_result(count, size_of::<SeedOccurrence>())?;
        let mut occurrences = Vec::new();
        occurrences
            .try_reserve_exact(count)
            .map_err(|_| SharedError::ResourceLimit)?;
        admit_result(occurrences.capacity(), size_of::<SeedOccurrence>())?;
        let mut start = 0;
        while start < member.occurrence_count {
            let block = self.occurrence_block_unchecked(group, member, start, 4096)?;
            if block.is_empty() {
                return Err(SharedError::Invalid("occurrence progress"));
            }
            start += block.len() as u64;
            occurrences.extend(block);
        }
        self.file.verify_unchanged()?;
        Ok(occurrences)
    }

    pub fn occurrence_block(
        &self,
        group: SharedGroup,
        member: SharedMember,
        start: u64,
        limit: usize,
    ) -> Result<Vec<SeedOccurrence>, SharedError> {
        self.validate_member(group, member)?;
        let block = self.occurrence_block_unchecked(group, member, start, limit)?;
        self.file.verify_unchanged()?;
        Ok(block)
    }

    pub fn metagenome(&self, id: u32) -> Result<Option<Metagenome<'_>>, SharedError> {
        self.begin_operation()?;
        if id >= self.file.header.document_count {
            return Ok(None);
        }
        let record = self.document_record(id)?;
        let result = Metagenome {
            id,
            name: self.resolve_string(record.name)?,
            bgzf_uri: self.resolve_string(record.bgzf_uri)?,
            bgzf_bytes: record.bgzf_bytes,
            bgzf_sha256: record.bgzf_sha256,
            contig_start: record.contig_start,
            contig_count: record.contig_count,
            gzi: self
                .file
                .section(Section::Gzi, record.gzi_offset, record.gzi_length)?,
        };
        self.file.verify_unchanged()?;
        Ok(Some(result))
    }

    pub fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, SharedError> {
        self.begin_operation()?;
        if id >= self.file.header.contig_count {
            return Ok(None);
        }
        let record = self.contig_record(id)?;
        let document = self.document_record(record.document_id)?;
        let document_end = document
            .contig_start
            .checked_add(document.contig_count)
            .ok_or(SharedError::Invalid("contig range"))?;
        if id < document.contig_start
            || id >= document_end
            || record.length == 0
            || record.line_bases == 0
            || record.line_width < record.line_bases
            || record.line_width > record.line_bases.saturating_add(2)
        {
            return Err(SharedError::Invalid("contig metadata"));
        }
        let result = Contig {
            id,
            metagenome_id: record.document_id,
            name: self.resolve_string(record.name)?,
            length: record.length,
            fasta_offset: record.fasta_offset,
            line_bases: record.line_bases,
            line_width: record.line_width,
        };
        self.file.verify_unchanged()?;
        Ok(Some(result))
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

    fn validate_member(&self, group: SharedGroup, member: SharedMember) -> Result<(), SharedError> {
        self.validate_group(group)?;
        if member.identity != group.identity || member.group != group.location {
            return Err(SharedError::Invalid("member handle"));
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
        self.observe(&self.core_inspections, 1);
        CoreRow::decode(
            bytes,
            self.file.header.contig_count,
            self.file.header.section(Section::Groups).length / GROUP_ROW_BYTES,
        )
    }

    fn singleton_member(&self, group: SharedGroup) -> Result<SharedMember, SharedError> {
        let GroupLocation::Singleton { core_ordinal } = group.location else {
            return Err(SharedError::Invalid("singleton group"));
        };
        let row = self.core_row(core_ordinal)?;
        let CoreKind::Singleton { contig_id, .. } = row.kind else {
            return Err(SharedError::Invalid("singleton group"));
        };
        let metagenome_id = self.contig_record(contig_id)?.document_id;
        Ok(SharedMember {
            identity: group.identity,
            group: group.location,
            metagenome_id,
            first_reference: core_ordinal,
            occurrence_count: 1,
        })
    }

    fn member_row(&self, group: SharedGroup, ordinal: u64) -> Result<SharedMember, SharedError> {
        let GroupLocation::Repeated {
            first_member,
            group_ordinal: _,
        } = group.location
        else {
            return Err(SharedError::Invalid("member group"));
        };
        let member_end = first_member
            .checked_add(u64::from(group.member_count))
            .ok_or(SharedError::Invalid("member ordinal"))?;
        if ordinal < first_member || ordinal >= member_end {
            return Err(SharedError::Invalid("member ordinal"));
        }
        let bytes = self
            .file
            .record(Section::Members, ordinal, MEMBER_ROW_BYTES)?;
        self.observe(&self.member_inspections, 1);
        let metagenome_id = read_u32(bytes, 0);
        let first_reference = read_u64(bytes, 8);
        let occurrence_count = read_u64(bytes, 16);
        let reference_count = self.file.header.section(Section::References).length / 8;
        if metagenome_id >= self.file.header.document_count
            || read_u32(bytes, 4) != 0
            || occurrence_count == 0
            || occurrence_count > group.occurrence_count
            || first_reference
                .checked_add(occurrence_count)
                .is_none_or(|end| end > reference_count)
        {
            return Err(SharedError::Invalid("member row"));
        }
        Ok(SharedMember {
            identity: group.identity,
            group: group.location,
            metagenome_id,
            first_reference,
            occurrence_count,
        })
    }

    fn occurrence_block_unchecked(
        &self,
        group: SharedGroup,
        member: SharedMember,
        start: u64,
        limit: usize,
    ) -> Result<Vec<SeedOccurrence>, SharedError> {
        if limit == 0 || limit > 4096 || start > member.occurrence_count {
            return Err(SharedError::Invalid("occurrence block"));
        }
        if start == member.occurrence_count {
            return Ok(Vec::new());
        }
        let count = usize::try_from((member.occurrence_count - start).min(limit as u64))
            .map_err(|_| SharedError::ResourceLimit)?;
        let flank = u64::from((group.key.length - 15) / 2);
        admit_result(count, size_of::<SeedOccurrence>())?;
        let mut output = Vec::new();
        output
            .try_reserve_exact(count)
            .map_err(|_| SharedError::ResourceLimit)?;
        admit_result(output.capacity(), size_of::<SeedOccurrence>())?;
        if let GroupLocation::Singleton { core_ordinal } = group.location {
            if start != 0 {
                return Err(SharedError::Invalid("singleton occurrence"));
            }
            let row = self.core_row(core_ordinal)?;
            let CoreKind::Singleton {
                contig_id,
                flags,
                position,
                ..
            } = row.kind
            else {
                return Err(SharedError::Invalid("singleton occurrence"));
            };
            self.validate_position(member.metagenome_id, contig_id, position, flank)?;
            self.observe(&self.positions_decoded, 1);
            output.push(SeedOccurrence {
                contig_id,
                position,
                canonical_orientation: flags & 1 != 0,
            });
            return Ok(output);
        }
        let first = member
            .first_reference
            .checked_add(start)
            .ok_or(SharedError::Invalid("reference range"))?;
        let references = self.file.section(
            Section::References,
            first
                .checked_mul(8)
                .ok_or(SharedError::Invalid("reference range"))?,
            count as u64 * 8,
        )?;
        self.observe(&self.references_decoded, count as u64);
        let occurrence_rows =
            self.file.header.section(Section::Occurrences).length / OCCURRENCE_ROW_BYTES;
        let (references, remainder) = references.as_chunks::<8>();
        if !remainder.is_empty() {
            return Err(SharedError::Invalid("reference records"));
        }
        for &raw in references {
            let ordinal = u64::from_le_bytes(raw);
            if ordinal >= occurrence_rows {
                return Err(SharedError::Invalid("occurrence reference"));
            }
            let bytes = self
                .file
                .record(Section::Occurrences, ordinal, OCCURRENCE_ROW_BYTES)?;
            let context = read_u32(bytes, 0);
            let contig_id = read_u32(bytes, 4);
            let flags = read_u32(bytes, 8);
            let position = read_u64(bytes, 16);
            if flags & !7 != 0
                || flags & 4 != 0 && flags & 2 == 0
                || flags & 2 == 0 && context != 0
                || flags & 4 == 0 && context & ((1 << 20) - 1) != 0
                || read_u32(bytes, 12) != 0
            {
                return Err(SharedError::Invalid("occurrence row"));
            }
            let seed = SharedSeed {
                core: group.key.core,
                context,
                flags: flags as u8,
                position,
            };
            if seed.key(group.key.length) != Some(group.key) {
                return Err(SharedError::Invalid("occurrence context"));
            }
            self.validate_position(member.metagenome_id, contig_id, position, flank)?;
            self.observe(&self.positions_decoded, 1);
            output.push(SeedOccurrence {
                contig_id,
                position,
                canonical_orientation: flags & 1 != 0,
            });
        }
        Ok(output)
    }

    fn validate_position(
        &self,
        metagenome_id: u32,
        contig_id: u32,
        position: u64,
        flank: u64,
    ) -> Result<(), SharedError> {
        if contig_id >= self.file.header.contig_count {
            return Err(SharedError::Invalid("occurrence contig"));
        }
        let contig = self.contig_record(contig_id)?;
        if contig.document_id != metagenome_id
            || position < flank
            || position
                .checked_add(15 + flank)
                .is_none_or(|end| end > contig.length)
        {
            return Err(SharedError::Invalid("occurrence position"));
        }
        Ok(())
    }

    fn document_record(&self, id: u32) -> Result<DocumentRecord, SharedError> {
        if id >= self.file.header.document_count {
            return Err(SharedError::Invalid("metagenome ID"));
        }
        Ok(DocumentRecord::decode(self.file.record(
            Section::Documents,
            u64::from(id),
            u64::from(DOCUMENT_RECORD_SIZE),
        )?)?)
    }

    fn contig_record(&self, id: u32) -> Result<ContigRecord, SharedError> {
        if id >= self.file.header.contig_count {
            return Err(SharedError::Invalid("contig ID"));
        }
        Ok(ContigRecord::decode(self.file.record(
            Section::Contigs,
            u64::from(id),
            u64::from(CONTIG_RECORD_SIZE),
        )?)?)
    }

    fn resolve_string(&self, reference: StringRef) -> Result<&str, SharedError> {
        let bytes = self.file.section(
            Section::Strings,
            u64::from(reference.offset),
            u64::from(reference.length),
        )?;
        std::str::from_utf8(bytes).map_err(|_| SharedError::Invalid("metadata string"))
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
        self.observe(&self.group_inspections, 1);
        let context_code = read_u64(bytes, 0);
        let first_member = read_u64(bytes, 8);
        let occurrence_count = read_u64(bytes, 16);
        let member_count = read_u32(bytes, 24);
        if !valid_context_code(context_code)
            || member_count == 0
            || occurrence_count == 0
            || occurrence_count > self.file.header.occurrence_count
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

    fn observe(&self, counter: &AtomicU64, count: u64) {
        if self.observed {
            counter.fetch_add(count, Ordering::Relaxed);
        }
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
