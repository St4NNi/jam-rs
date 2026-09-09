use crate::jidx::{sha256, sha256_reader};
use crate::jidx_writer::sync_directory;
use crate::owner_format::{
    COMPLETE_RANGE, HAS_METADATA, OWNER_BLOCK_SIZE, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE,
    OWNER_HEADER_SIZE, OWNER_PAGE_SIZE, OWNER_SECTION_COUNT, OwnerHeader, OwnerSection,
    OwnerSectionDescriptor, checksum_layout, put_u32, put_u64,
};
use crate::owner_postings::{MAX_KEYS_PER_BLOCK, OwnerKey, encode_block};
use serde::Serialize;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

const MAX_PROTOTYPE_KEYS: usize = 1_000_000;
pub(crate) const MAX_PROTOTYPE_ENCODED_BYTES: u64 = 512 * 1024 * 1024;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerKeyRange {
    pub first: u64,
    pub last: u64,
    pub complete: bool,
}

#[derive(Clone, Debug)]
pub struct OwnerMetadata {
    pub metagenomes: Vec<OwnerMetagenomeInput>,
}

#[derive(Clone, Debug)]
pub struct OwnerMetagenomeInput {
    pub name: String,
    pub bgzf_uri: String,
    pub bgzf_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub gzi: Vec<u8>,
    pub contigs: Vec<OwnerContigInput>,
}

#[derive(Clone, Debug)]
pub struct OwnerContigInput {
    pub name: String,
    pub length: u64,
    pub fasta_offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
}

pub struct OwnerWriteInput<'a> {
    pub owner_ordinal: u32,
    pub owner_count: u32,
    pub range: OwnerKeyRange,
    pub k: u8,
    pub rescue_k15: bool,
    pub minimizer_window: u16,
    pub generation_id: [u8; 32],
    pub document_count: u32,
    pub contig_count: u32,
    pub keys: &'a [OwnerKey],
    pub metadata: Option<&'a OwnerMetadata>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerWriteStats {
    pub keys: u64,
    pub occurrences: u64,
    pub hot_bytes: u64,
    pub cold_bytes: u64,
    pub file_bytes: u64,
    pub header_sha256: [u8; 32],
    pub file_sha256: [u8; 32],
}

fn validate_input(input: &OwnerWriteInput<'_>) -> Result<(), OwnerWriteError> {
    if input.owner_count == 0
        || input.owner_ordinal >= input.owner_count
        || input.range.first > input.range.last
        || !(1..=32).contains(&input.k)
        || (input.rescue_k15 && input.k != 21)
        || input.minimizer_window == 0
        || input.generation_id == [0; 32]
        || input.document_count == 0
        || input.contig_count == 0
        || input.keys.len() > MAX_PROTOTYPE_KEYS
    {
        return Err(OwnerWriteError::Invalid("writer input"));
    }
    if input.keys.windows(2).any(|pair| pair[0].key >= pair[1].key)
        || input
            .keys
            .iter()
            .any(|key| key.key < input.range.first || key.key > input.range.last)
    {
        return Err(OwnerWriteError::Invalid("key range or order"));
    }
    if let Some(metadata) = input.metadata {
        let mut bytes = 0;
        for metagenome in &metadata.metagenomes {
            bytes = validate_prototype_payload_bytes(&[
                bytes,
                u64::from(OWNER_DOCUMENT_SIZE),
                metagenome.name.len() as u64,
                metagenome.bgzf_uri.len() as u64,
                metagenome.gzi.len() as u64,
            ])?;
            for contig in &metagenome.contigs {
                bytes = validate_prototype_payload_bytes(&[
                    bytes,
                    u64::from(OWNER_CONTIG_SIZE),
                    contig.name.len() as u64,
                ])?;
            }
        }
        let metadata_contigs = metadata
            .metagenomes
            .iter()
            .try_fold(0usize, |total, entry| {
                total
                    .checked_add(entry.contigs.len())
                    .ok_or(OwnerWriteError::Invalid("metadata counts"))
            })?;
        if metadata.metagenomes.len() != input.document_count as usize
            || metadata_contigs != input.contig_count as usize
        {
            return Err(OwnerWriteError::Invalid("metadata counts"));
        }
    }
    for key in input.keys {
        crate::jidx::seed_length(input.k, input.rescue_k15, key.key)
            .map_err(|_| OwnerWriteError::Invalid("seed key"))?;
        for member in &key.members {
            if member.document_id >= input.document_count {
                return Err(OwnerWriteError::Invalid("document ID"));
            }
        }
    }
    Ok(())
}

fn encode_metadata(
    metadata: &OwnerMetadata,
    strings: &mut Vec<u8>,
    documents: &mut Vec<u8>,
    contigs: &mut Vec<u8>,
    gzi: &mut Vec<u8>,
) -> Result<(), OwnerWriteError> {
    let mut names = HashSet::new();
    let mut contig_start = 0u32;
    for (document_id, metagenome) in metadata.metagenomes.iter().enumerate() {
        if !names.insert(metagenome.name.as_str())
            || metagenome.bgzf_bytes == 0
            || metagenome.bgzf_sha256 == [0; 32]
            || metagenome.contigs.is_empty()
        {
            return Err(OwnerWriteError::Invalid("metagenome metadata"));
        }
        let mut gzi_reader = noodles_bgzf::gzi::io::Reader::new(metagenome.gzi.as_slice());
        gzi_reader.read_index()?;
        let name = push_string(strings, &metagenome.name)?;
        let uri = push_string(strings, &metagenome.bgzf_uri)?;
        let mut record = [0; OWNER_DOCUMENT_SIZE as usize];
        put_u32(&mut record, 0, name.0);
        put_u32(&mut record, 4, name.1);
        put_u32(&mut record, 8, uri.0);
        put_u32(&mut record, 12, uri.1);
        put_u64(&mut record, 16, metagenome.bgzf_bytes);
        record[24..56].copy_from_slice(&metagenome.bgzf_sha256);
        put_u32(&mut record, 56, contig_start);
        put_u32(
            &mut record,
            60,
            u32::try_from(metagenome.contigs.len())
                .map_err(|_| OwnerWriteError::Invalid("contig count"))?,
        );
        put_u64(&mut record, 64, gzi.len() as u64);
        put_u64(&mut record, 72, metagenome.gzi.len() as u64);
        documents.extend_from_slice(&record);
        let mut contig_names = HashSet::new();
        for contig in &metagenome.contigs {
            if !contig_names.insert(contig.name.as_str())
                || contig.length == 0
                || contig.line_bases == 0
                || contig.line_width < contig.line_bases
                || contig.line_width > contig.line_bases.saturating_add(2)
            {
                return Err(OwnerWriteError::Invalid("contig metadata"));
            }
            let name = push_string(strings, &contig.name)?;
            let mut record = [0; OWNER_CONTIG_SIZE as usize];
            put_u32(&mut record, 0, document_id as u32);
            put_u32(&mut record, 4, name.0);
            put_u32(&mut record, 8, name.1);
            put_u64(&mut record, 16, contig.length);
            put_u64(&mut record, 24, contig.fasta_offset);
            put_u32(&mut record, 32, contig.line_bases);
            put_u32(&mut record, 36, contig.line_width);
            contigs.extend_from_slice(&record);
        }
        contig_start = contig_start
            .checked_add(
                u32::try_from(metagenome.contigs.len())
                    .map_err(|_| OwnerWriteError::Invalid("contig count"))?,
            )
            .ok_or(OwnerWriteError::Invalid("contig count"))?;
        gzi.extend_from_slice(&metagenome.gzi);
    }
    Ok(())
}

fn push_string(output: &mut Vec<u8>, value: &str) -> Result<(u32, u32), OwnerWriteError> {
    if value.is_empty() || value.bytes().any(|byte| matches!(byte, 0 | b'\n' | b'\r')) {
        return Err(OwnerWriteError::Invalid("metadata string"));
    }
    let offset = u32::try_from(output.len()).map_err(|_| OwnerWriteError::Invalid("strings"))?;
    let length = u32::try_from(value.len()).map_err(|_| OwnerWriteError::Invalid("strings"))?;
    output.extend_from_slice(value.as_bytes());
    Ok((offset, length))
}

pub(crate) fn validate_prototype_payload_bytes(lengths: &[u64]) -> Result<u64, OwnerWriteError> {
    let total = lengths.iter().try_fold(0u64, |total, length| {
        total
            .checked_add(*length)
            .ok_or(OwnerWriteError::Invalid("prototype encoded bytes"))
    })?;
    if total > MAX_PROTOTYPE_ENCODED_BYTES {
        return Err(OwnerWriteError::Invalid("prototype encoded bytes"));
    }
    Ok(total)
}

fn build_checksum_tree(
    input: &mut BufReader<File>,
    data_pages: u64,
) -> Result<(Vec<u8>, [u8; 32]), OwnerWriteError> {
    let levels = checksum_layout(data_pages)?;
    let capacity = levels.iter().try_fold(0usize, |total, level| {
        level
            .page_count
            .checked_mul(OWNER_PAGE_SIZE)
            .and_then(|bytes| usize::try_from(bytes).ok())
            .and_then(|bytes| total.checked_add(bytes))
            .ok_or(OwnerWriteError::Invalid("checksum length"))
    })?;
    let mut tree = Vec::with_capacity(capacity);
    let mut level = Vec::new();
    let mut page = [0; OWNER_PAGE_SIZE as usize];
    for _ in 0..data_pages {
        input.read_exact(&mut page)?;
        level.extend_from_slice(&sha256(&page));
    }
    pad_checksum_level(&mut level)?;
    loop {
        tree.extend_from_slice(&level);
        if level.len() == OWNER_PAGE_SIZE as usize {
            return Ok((tree, sha256(&level)));
        }
        let mut next = Vec::with_capacity(level.len() / OWNER_PAGE_SIZE as usize * 32);
        for page in level.as_chunks::<{ OWNER_PAGE_SIZE as usize }>().0 {
            next.extend_from_slice(&sha256(page));
        }
        pad_checksum_level(&mut next)?;
        level = next;
    }
}

fn pad_checksum_level(level: &mut Vec<u8>) -> Result<(), OwnerWriteError> {
    let length = level
        .len()
        .max(1)
        .checked_add(OWNER_PAGE_SIZE as usize - 1)
        .ok_or(OwnerWriteError::Invalid("checksum length"))?
        / OWNER_PAGE_SIZE as usize
        * OWNER_PAGE_SIZE as usize;
    level.resize(length, 0);
    Ok(())
}

fn section_layout(
    lengths: [u64; OWNER_SECTION_COUNT],
) -> Result<([OwnerSectionDescriptor; OWNER_SECTION_COUNT], u64), OwnerWriteError> {
    let mut offset = OWNER_HEADER_SIZE as u64;
    let mut sections = Vec::with_capacity(OWNER_SECTION_COUNT);
    for (kind, mut length) in OwnerSection::ALL.into_iter().zip(lengths) {
        offset = offset
            .checked_add(OWNER_PAGE_SIZE - 1)
            .ok_or(OwnerWriteError::Invalid("file length"))?
            & !(OWNER_PAGE_SIZE - 1);
        if kind == OwnerSection::PageChecksums {
            length = checksum_layout(offset / OWNER_PAGE_SIZE - 1)?
                .iter()
                .try_fold(0u64, |total, level| {
                    total
                        .checked_add(
                            level
                                .page_count
                                .checked_mul(OWNER_PAGE_SIZE)
                                .ok_or(OwnerWriteError::Invalid("checksum length"))?,
                        )
                        .ok_or(OwnerWriteError::Invalid("checksum length"))
                })?;
        }
        sections.push(OwnerSectionDescriptor {
            kind,
            record_size: kind.record_size(),
            offset,
            length,
        });
        offset = offset
            .checked_add(length)
            .ok_or(OwnerWriteError::Invalid("file length"))?;
    }
    Ok((sections.try_into().expect("fixed section count"), offset))
}

fn write_padding(output: &mut BufWriter<&mut File>, target: u64) -> io::Result<()> {
    let count = target
        .checked_sub(output.stream_position()?)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "owner section overlap"))?;
    output.write_all(&vec![0; count])
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[derive(Debug, Error)]
pub enum OwnerWriteError {
    #[error("owner index I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("invalid owner index: {0}")]
    Invalid(&'static str),
    #[error(transparent)]
    Postings(#[from] crate::owner_postings::OwnerPostingsError),
    #[error(transparent)]
    Header(#[from] crate::owner_format::OwnerReaderError),
    #[error("owner manifest JSON failed: {0}")]
    Json(#[from] serde_json::Error),
}
