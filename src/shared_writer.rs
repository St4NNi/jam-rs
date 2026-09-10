use crate::bgzf::BgzfReader;
use crate::jidx::{
    ContigRecord, DocumentRecord, StringRef, put_u32, put_u64, sha256, sha256_reader,
};
use crate::jidx_reader::JidxReader;
use crate::owner_format::checksum_layout;
use crate::shared_format::{
    HEADER_BYTES, MULTIPLE_CORE, PAGE_BYTES, Section, SectionRange, SharedError, SharedHeader,
};
use crate::shared_seed::{SharedSeed, select_shared_seeds};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

const BUILD_BYTES: usize = 1024 * 1024 * 1024;
const SOURCE_CHUNK_BYTES: u64 = 1024 * 1024;

#[derive(Clone, Copy, Debug)]
pub(crate) struct IndexedSeed {
    pub(crate) member: u32,
    pub(crate) contig: u32,
    pub(crate) seed: SharedSeed,
}

#[derive(Debug, Serialize)]
pub struct SharedBuildStats {
    pub window: u16,
    pub source_bases: u64,
    pub core_count: u64,
    pub occurrence_count: u64,
    pub singleton_cores: u64,
    pub context_groups: u64,
    pub member_descriptors: u64,
    pub occurrence_references: u64,
    pub repeated_core_positions: u64,
    pub index_bytes: u64,
    pub bgzf_bytes_once: u64,
    pub complete_query_ready_bytes: u64,
    pub section_bytes: BTreeMap<String, u64>,
}

pub fn build_shared_index(
    reference: impl AsRef<Path>,
    output: impl AsRef<Path>,
    window: u16,
) -> Result<SharedBuildStats, SharedError> {
    if window == 0 {
        return Err(SharedError::Invalid("minimizer window"));
    }
    let reference = JidxReader::open(reference)?;
    reference.verify_checksum()?;
    let mut seeds = Vec::new();
    let halo = u64::from(window) + 21;
    for member in 0..reference.header().document_count {
        let source = reference
            .metagenome(member)?
            .ok_or(SharedError::Invalid("source metadata"))?;
        let mut reader = BgzfReader::open(source, None, true)?;
        let end = source
            .contig_start
            .checked_add(source.contig_count)
            .ok_or(SharedError::Invalid("contig IDs"))?;
        for id in source.contig_start..end {
            let contig = reference
                .contig(id)?
                .ok_or(SharedError::Invalid("contig metadata"))?;
            let mut core_start = 0;
            while core_start < contig.length {
                let core_end = core_start
                    .saturating_add(SOURCE_CHUNK_BYTES)
                    .min(contig.length);
                let start = core_start.saturating_sub(halo);
                let end = core_end.saturating_add(halo).min(contig.length);
                let sequence = reader.read_contig_range(contig, start, end)?;
                for mut seed in select_shared_seeds(&sequence, window)? {
                    seed.position = seed
                        .position
                        .checked_add(start)
                        .ok_or(SharedError::Invalid("source position"))?;
                    if (core_start..core_end).contains(&seed.position) {
                        if seeds.len() >= BUILD_BYTES / std::mem::size_of::<IndexedSeed>() / 4 {
                            return Err(SharedError::ResourceLimit);
                        }
                        seeds.push(IndexedSeed {
                            member,
                            contig: id,
                            seed,
                        });
                    }
                }
                core_start = core_end;
            }
        }
    }
    write_shared_index(&reference, output.as_ref(), window, &mut seeds)
}

pub(crate) fn write_shared_index(
    reference: &JidxReader,
    output: &Path,
    window: u16,
    seeds: &mut [IndexedSeed],
) -> Result<SharedBuildStats, SharedError> {
    let mut sections: [Vec<u8>; 10] = std::array::from_fn(|_| Vec::new());
    let (source_bases, bgzf_bytes_once) = metadata(reference, &mut sections)?;
    seeds.sort_unstable_by_key(|row| (row.seed.core, row.member, row.contig, row.seed.position));
    for rows in seeds.windows(2) {
        if rows[0].seed.core == rows[1].seed.core
            && rows[0].contig == rows[1].contig
            && rows[0].seed.position == rows[1].seed.position
        {
            return Err(SharedError::Invalid("duplicate physical anchor"));
        }
    }
    let mut singletons = 0;
    for rows in seeds.chunk_by(|left, right| left.seed.core == right.seed.core) {
        for row in rows {
            let contig = reference
                .contig(row.contig)?
                .ok_or(SharedError::Invalid("anchor contig"))?;
            if row.member != contig.metagenome_id
                || row.seed.core >= 1 << 30
                || row.seed.flags & !7 != 0
                || row
                    .seed
                    .position
                    .checked_add(15)
                    .is_none_or(|end| end > contig.length)
                || row.seed.flags & 4 != 0 && row.seed.flags & 2 == 0
            {
                return Err(SharedError::Invalid("anchor location"));
            }
        }
        let mut core = [0u8; 24];
        put_u32(&mut core, 0, rows[0].seed.core);
        if rows.len() == 1 {
            let row = rows[0];
            put_u32(&mut core, 4, row.seed.context);
            put_u32(&mut core, 8, row.contig);
            put_u32(&mut core, 12, u32::from(row.seed.flags));
            put_u64(&mut core, 16, row.seed.position);
            singletons += 1;
        } else {
            let first_group = sections[Section::Groups as usize].len() as u64 / 32;
            let first_position = sections[Section::Occurrences as usize].len() as u64 / 24;
            let mut groups = BTreeMap::<u64, BTreeMap<u32, Vec<u64>>>::new();
            for (offset, row) in rows.iter().enumerate() {
                let mut occurrence = [0; 24];
                put_u32(&mut occurrence, 0, row.seed.context);
                put_u32(&mut occurrence, 4, row.contig);
                put_u32(&mut occurrence, 8, u32::from(row.seed.flags));
                put_u64(&mut occurrence, 16, row.seed.position);
                sections[Section::Occurrences as usize].extend_from_slice(&occurrence);
                for length in [15, 21, 31] {
                    if let Some(key) = row.seed.key(length) {
                        groups
                            .entry(
                                key.context_code()
                                    .ok_or(SharedError::Invalid("context key"))?,
                            )
                            .or_default()
                            .entry(row.member)
                            .or_default()
                            .push(first_position + offset as u64);
                    }
                }
            }
            put_u32(&mut core, 0, rows[0].seed.core | MULTIPLE_CORE);
            put_u32(
                &mut core,
                4,
                u32::try_from(groups.len()).map_err(|_| SharedError::ResourceLimit)?,
            );
            put_u64(&mut core, 8, rows.len() as u64);
            put_u64(&mut core, 16, first_group);
            for (code, members) in groups {
                let mut group = [0; 32];
                put_u64(&mut group, 0, code);
                put_u64(
                    &mut group,
                    8,
                    sections[Section::Members as usize].len() as u64 / 24,
                );
                let count = members
                    .values()
                    .try_fold(0u64, |sum, values| sum.checked_add(values.len() as u64))
                    .ok_or(SharedError::ResourceLimit)?;
                put_u64(&mut group, 16, count);
                put_u32(
                    &mut group,
                    24,
                    u32::try_from(members.len()).map_err(|_| SharedError::ResourceLimit)?,
                );
                for (id, references) in members {
                    let mut member = [0; 24];
                    put_u32(&mut member, 0, id);
                    put_u64(
                        &mut member,
                        8,
                        sections[Section::References as usize].len() as u64 / 8,
                    );
                    put_u64(&mut member, 16, references.len() as u64);
                    sections[Section::Members as usize].extend_from_slice(&member);
                    for ordinal in references {
                        sections[Section::References as usize]
                            .extend_from_slice(&ordinal.to_le_bytes());
                    }
                }
                sections[Section::Groups as usize].extend_from_slice(&group);
            }
        }
        sections[Section::Cores as usize].extend_from_slice(&core);
        if sections
            .iter()
            .map(|section| section.capacity())
            .sum::<usize>()
            > BUILD_BYTES
        {
            return Err(SharedError::ResourceLimit);
        }
    }
    publish(
        reference,
        output,
        window,
        sections,
        source_bases,
        bgzf_bytes_once,
        seeds.len() as u64,
        singletons,
    )
}

fn metadata(
    reference: &JidxReader,
    sections: &mut [Vec<u8>; 10],
) -> Result<(u64, u64), SharedError> {
    let mut bases = 0u64;
    let mut objects = BTreeSet::new();
    let mut sequence_bytes = 0u64;
    for id in 0..reference.header().document_count {
        let source = reference
            .metagenome(id)?
            .ok_or(SharedError::Invalid("source metadata"))?;
        let record = DocumentRecord {
            name: string(&mut sections[Section::Strings as usize], source.name)?,
            bgzf_uri: string(&mut sections[Section::Strings as usize], source.bgzf_uri)?,
            bgzf_bytes: source.bgzf_bytes,
            contig_start: source.contig_start,
            contig_count: source.contig_count,
            bgzf_sha256: source.bgzf_sha256,
            gzi_offset: sections[Section::Gzi as usize].len() as u64,
            gzi_length: source.gzi.len() as u64,
        };
        if objects.insert(source.bgzf_uri.to_owned()) {
            sequence_bytes = sequence_bytes
                .checked_add(source.bgzf_bytes)
                .ok_or(SharedError::ResourceLimit)?;
        }
        sections[Section::Documents as usize].extend_from_slice(&record.encode());
        sections[Section::Gzi as usize].extend_from_slice(source.gzi);
    }
    for id in 0..reference.header().contig_count {
        let contig = reference
            .contig(id)?
            .ok_or(SharedError::Invalid("contig metadata"))?;
        let record = ContigRecord {
            document_id: contig.metagenome_id,
            name: string(&mut sections[Section::Strings as usize], contig.name)?,
            length: contig.length,
            fasta_offset: contig.fasta_offset,
            line_bases: contig.line_bases,
            line_width: contig.line_width,
        };
        sections[Section::Contigs as usize].extend_from_slice(&record.encode());
        bases = bases
            .checked_add(contig.length)
            .ok_or(SharedError::ResourceLimit)?;
    }
    Ok((bases, sequence_bytes))
}

fn string(strings: &mut Vec<u8>, value: &str) -> Result<StringRef, SharedError> {
    let record = StringRef {
        offset: u32::try_from(strings.len()).map_err(|_| SharedError::ResourceLimit)?,
        length: u32::try_from(value.len()).map_err(|_| SharedError::ResourceLimit)?,
    };
    strings.extend_from_slice(value.as_bytes());
    Ok(record)
}

fn publish(
    reference: &JidxReader,
    output: &Path,
    window: u16,
    sections: [Vec<u8>; 10],
    source_bases: u64,
    bgzf_bytes_once: u64,
    occurrences: u64,
    singletons: u64,
) -> Result<SharedBuildStats, SharedError> {
    let mut offset = HEADER_BYTES as u64;
    let mut ranges = [SectionRange::default(); 10];
    for kind in Section::ALL {
        offset = offset
            .checked_next_multiple_of(PAGE_BYTES)
            .ok_or(SharedError::ResourceLimit)?;
        let length = sections[kind as usize].len() as u64;
        ranges[kind as usize] = SectionRange { offset, length };
        offset = offset
            .checked_add(length)
            .ok_or(SharedError::ResourceLimit)?;
    }
    let checksum_start = ranges[Section::Checksums as usize].offset;
    let data_pages = checksum_start / PAGE_BYTES - 1;
    let levels =
        checksum_layout(data_pages).map_err(|_| SharedError::Invalid("checksum layout"))?;
    let top = levels.last().ok_or(SharedError::Invalid("checksum root"))?;
    let checksum_bytes = top.offset + top.page_count * PAGE_BYTES;
    ranges[Section::Checksums as usize].length = checksum_bytes;
    let file_bytes = checksum_start
        .checked_add(checksum_bytes)
        .ok_or(SharedError::ResourceLimit)?;
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut temporary = tempfile::Builder::new()
        .prefix(".jam-shared-")
        .tempfile_in(parent)?;
    let file = temporary.as_file_mut();
    file.set_len(file_bytes)?;
    for kind in Section::ALL {
        file.seek(SeekFrom::Start(ranges[kind as usize].offset))?;
        file.write_all(&sections[kind as usize])?;
    }
    let mut hashes = Vec::with_capacity(data_pages as usize);
    let mut page = [0; PAGE_BYTES as usize];
    file.seek(SeekFrom::Start(PAGE_BYTES))?;
    for _ in 0..data_pages {
        file.read_exact(&mut page)?;
        hashes.push(sha256(&page));
    }
    for level in &levels {
        file.seek(SeekFrom::Start(checksum_start + level.offset))?;
        let mut parents = Vec::new();
        for chunk in hashes.chunks(128) {
            page.fill(0);
            for (index, hash) in chunk.iter().enumerate() {
                page[32 * index..32 * index + 32].copy_from_slice(hash);
            }
            file.write_all(&page)?;
            parents.push(sha256(&page));
        }
        hashes = parents;
    }
    let root = *hashes
        .first()
        .ok_or(SharedError::Invalid("empty checksum root"))?;
    file.seek(SeekFrom::Start(HEADER_BYTES as u64))?;
    let body_sha256 = sha256_reader(&mut *file)?;
    let header = SharedHeader {
        window,
        core_count: ranges[Section::Cores as usize].length / 24,
        occurrence_count: occurrences,
        document_count: reference.header().document_count,
        contig_count: reference.header().contig_count,
        source_bases,
        manifest_sha256: reference.header().manifest_sha256,
        body_sha256,
        checksum_root_sha256: root,
        sections: ranges,
    };
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&header.encode()?)?;
    file.sync_all()?;
    temporary
        .persist_noclobber(output)
        .map_err(|error| SharedError::Io(error.error))?;
    std::fs::File::open(parent)?.sync_all()?;
    Ok(SharedBuildStats {
        window,
        source_bases,
        core_count: header.core_count,
        occurrence_count: occurrences,
        singleton_cores: singletons,
        context_groups: ranges[Section::Groups as usize].length / 32,
        member_descriptors: ranges[Section::Members as usize].length / 24,
        occurrence_references: ranges[Section::References as usize].length / 8,
        repeated_core_positions: ranges[Section::Occurrences as usize].length / 24,
        index_bytes: file_bytes,
        bgzf_bytes_once,
        complete_query_ready_bytes: file_bytes
            .checked_add(bgzf_bytes_once)
            .ok_or(SharedError::ResourceLimit)?,
        section_bytes: Section::ALL
            .into_iter()
            .map(|kind| (format!("{kind:?}"), ranges[kind as usize].length))
            .collect(),
    })
}
