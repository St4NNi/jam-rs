use crate::jidx::{sha256, sha256_reader};
use crate::jidx_writer::sync_directory;
use crate::owner_format::{
    COMPLETE_RANGE, HAS_METADATA, OWNER_BLOCK_SIZE, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE,
    OWNER_HEADER_SIZE, OWNER_PAGE_SIZE, OWNER_SECTION_COUNT, OwnerHeader, OwnerSection,
    OwnerSectionDescriptor, checksum_layout, put_u32, put_u64,
};
use crate::owner_postings::{
    MAX_KEYS_PER_BLOCK, OwnerKey, OwnerOccurrence, OwnerPostingsError, encode_block,
};
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerMetadata {
    pub metagenomes: Vec<OwnerMetagenomeInput>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerMetagenomeInput {
    pub name: String,
    pub bgzf_uri: String,
    pub bgzf_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub gzi: Vec<u8>,
    pub original_contig_start: u32,
    pub original_contig_count: u32,
    pub locus_bits: u8,
    pub contigs: Vec<OwnerContigInput>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerContigInput {
    pub local_contig: u32,
    pub name: String,
    pub length: u64,
    pub fasta_offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
}

impl OwnerMetadata {
    pub fn locus(
        &self,
        document: u32,
        occurrence: OwnerOccurrence,
    ) -> Result<u64, OwnerPostingsError> {
        let invalid = || OwnerPostingsError("locus metadata");
        let metagenome = self
            .metagenomes
            .get(document as usize)
            .ok_or_else(invalid)?;
        let index = metagenome
            .contigs
            .binary_search_by_key(&occurrence.local_contig, |contig| contig.local_contig)
            .map_err(|_| invalid())?;
        let contig = &metagenome.contigs[index];
        if occurrence.position >= contig.length
            || contig.line_bases == 0
            || contig.line_width < contig.line_bases
        {
            return Err(invalid());
        }
        let line_bases = u64::from(contig.line_bases);
        (occurrence.position / line_bases)
            .checked_mul(u64::from(contig.line_width))
            .and_then(|offset| offset.checked_add(occurrence.position % line_bases))
            .and_then(|offset| contig.fasta_offset.checked_add(offset))
            .ok_or_else(invalid)
    }
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
    pub loci: &'a OwnerMetadata,
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

#[derive(Serialize)]
struct PublishedManifest<'a> {
    version: u32,
    complete: bool,
    metadata_owner: u32,
    owners: &'a [PublishedOwner],
}

#[derive(Serialize)]
struct PublishedOwner {
    path: PathBuf,
    header_sha256: String,
}

pub fn publish_owner_manifest(
    path: impl AsRef<Path>,
    owners: &[PathBuf],
    metadata_owner: u32,
    complete: bool,
) -> Result<[u8; 32], OwnerWriteError> {
    let path = path.as_ref();
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut entries = Vec::with_capacity(owners.len());
    for owner in owners {
        let final_path = if owner.is_absolute() {
            owner.clone()
        } else {
            parent.join(owner)
        };
        let relative = final_path
            .strip_prefix(parent)
            .map_err(|_| OwnerWriteError::Invalid("owner manifest path"))?;
        if relative.as_os_str().is_empty()
            || relative
                .components()
                .any(|part| !matches!(part, std::path::Component::Normal(_)))
        {
            return Err(OwnerWriteError::Invalid("owner manifest path"));
        }
        let mut file = File::open(&final_path)?;
        let mut header = [0; OWNER_HEADER_SIZE];
        file.read_exact(&mut header)?;
        OwnerHeader::decode(&header, file.metadata()?.len())?;
        entries.push(PublishedOwner {
            path: relative.to_path_buf(),
            header_sha256: hex(&sha256(&header)),
        });
    }
    let bytes = serde_json::to_vec(&PublishedManifest {
        version: 1,
        complete,
        metadata_owner,
        owners: &entries,
    })?;
    let mut temporary = tempfile::Builder::new()
        .prefix(".owner-manifest-")
        .tempfile_in(parent)?;
    temporary.write_all(&bytes)?;
    temporary.flush()?;
    temporary.as_file().sync_all()?;
    crate::owner_reader::OwnerReader::open(temporary.path())?.verify_checksum()?;
    temporary
        .persist_noclobber(path)
        .map_err(|error| error.error)?;
    sync_directory(parent)?;
    Ok(sha256(&bytes))
}

pub fn write_owner(
    path: impl AsRef<Path>,
    input: OwnerWriteInput<'_>,
) -> Result<OwnerWriteStats, OwnerWriteError> {
    validate_input(&input)?;
    let mut strings = Vec::new();
    let mut documents = Vec::new();
    let mut contigs = Vec::new();
    let mut gzi = Vec::new();
    encode_metadata(
        input.loci,
        &mut strings,
        &mut documents,
        &mut contigs,
        &mut gzi,
    )?;
    let metadata_sha256 =
        crate::owner_format::metadata_digest([&strings, &documents, &contigs, &gzi]);
    if input.metadata.is_none() {
        strings = Vec::new();
        documents = Vec::new();
        contigs = Vec::new();
        gzi = Vec::new();
    }
    let mut directory = Vec::new();
    let mut hot = Vec::new();
    let mut cold = Vec::new();
    let mut occurrence_count = 0u64;
    let widths = input
        .loci
        .metagenomes
        .iter()
        .map(|document| document.locus_bits)
        .collect::<Vec<_>>();
    for block in input.keys.chunks(MAX_KEYS_PER_BLOCK) {
        let encoded = encode_block(block, &widths, |document, occurrence| {
            input.loci.locus(document, occurrence)
        })?;
        let first = block.first().expect("nonempty block").key;
        let last = block.last().expect("nonempty block").key;
        let mut record = [0; OWNER_BLOCK_SIZE as usize];
        put_u64(&mut record, 0, first);
        put_u64(&mut record, 8, last);
        put_u64(&mut record, 16, hot.len() as u64);
        put_u64(&mut record, 24, encoded.hot.len() as u64);
        put_u64(&mut record, 32, cold.len() as u64);
        put_u64(&mut record, 40, encoded.cold.len() as u64);
        put_u32(&mut record, 48, block.len() as u32);
        directory.extend_from_slice(&record);
        let posting_bytes = hot
            .len()
            .checked_add(cold.len())
            .and_then(|length| length.checked_add(encoded.hot.len()))
            .and_then(|length| length.checked_add(encoded.cold.len()))
            .ok_or(OwnerWriteError::Invalid("prototype posting bytes"))?;
        validate_prototype_payload_bytes(&[
            strings.len() as u64,
            documents.len() as u64,
            contigs.len() as u64,
            gzi.len() as u64,
            directory.len() as u64,
            posting_bytes as u64,
        ])?;
        hot.extend_from_slice(&encoded.hot);
        cold.extend_from_slice(&encoded.cold);
        for key in block {
            for member in &key.members {
                occurrence_count = occurrence_count
                    .checked_add(member.occurrences.len() as u64)
                    .ok_or(OwnerWriteError::Invalid("occurrence count"))?;
            }
        }
    }
    let lengths = [
        strings.len() as u64,
        documents.len() as u64,
        contigs.len() as u64,
        gzi.len() as u64,
        directory.len() as u64,
        hot.len() as u64,
        cold.len() as u64,
        0,
    ];
    validate_prototype_payload_bytes(&lengths[..7])?;
    let (sections, expected_file_bytes) = section_layout(lengths)?;
    let path = path.as_ref();
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let mut temporary = tempfile::Builder::new()
        .prefix(".owner-")
        .tempfile_in(parent)?;
    let checksum_input = temporary.reopen()?;
    let mut output = BufWriter::with_capacity(1024 * 1024, temporary.as_file_mut());
    output.write_all(&[0; OWNER_HEADER_SIZE])?;
    for (section, payload) in sections[..7].iter().zip([
        &strings, &documents, &contigs, &gzi, &directory, &hot, &cold,
    ]) {
        write_padding(&mut output, section.offset)?;
        output.write_all(payload)?;
    }
    write_padding(&mut output, sections[7].offset)?;
    output.flush()?;
    let mut checksum_input = BufReader::new(checksum_input);
    checksum_input.seek(SeekFrom::Start(OWNER_PAGE_SIZE))?;
    let (checksum_tree, checksum_root_sha256) = build_checksum_tree(
        &mut checksum_input,
        sections[7].offset / OWNER_PAGE_SIZE - 1,
    )?;
    if checksum_tree.len() as u64 != sections[7].length {
        return Err(OwnerWriteError::Invalid("checksum length"));
    }
    output.write_all(&checksum_tree)?;
    if output.stream_position()? != expected_file_bytes {
        return Err(OwnerWriteError::Invalid("written length"));
    }
    output.flush()?;
    output.seek(SeekFrom::Start(OWNER_HEADER_SIZE as u64))?;
    let body_sha256 = sha256_reader(&mut **output.get_mut())?;
    let header = OwnerHeader {
        flags: (if input.range.complete {
            COMPLETE_RANGE
        } else {
            0
        }) | (if input.metadata.is_some() {
            HAS_METADATA
        } else {
            0
        }),
        k: input.k,
        rescue_k15: input.rescue_k15,
        minimizer_window: input.minimizer_window,
        owner_ordinal: input.owner_ordinal,
        owner_count: input.owner_count,
        first_key: input.range.first,
        last_key: input.range.last,
        key_count: input.keys.len() as u64,
        occurrence_count,
        document_count: input.document_count,
        contig_count: input.contig_count,
        generation_id: input.generation_id,
        body_sha256,
        checksum_root_sha256,
        metadata_sha256,
        sections,
    };
    let header_bytes = header.encode()?;
    output.seek(SeekFrom::Start(0))?;
    output.write_all(&header_bytes)?;
    output.flush()?;
    output.get_ref().sync_all()?;
    output.seek(SeekFrom::Start(0))?;
    let file_sha256 = sha256_reader(&mut **output.get_mut())?;
    drop(output);
    temporary
        .persist_noclobber(path)
        .map_err(|error| error.error)?;
    sync_directory(parent)?;
    Ok(OwnerWriteStats {
        keys: input.keys.len() as u64,
        occurrences: occurrence_count,
        hot_bytes: hot.len() as u64,
        cold_bytes: cold.len() as u64,
        file_bytes: expected_file_bytes,
        header_sha256: sha256(&header_bytes),
        file_sha256,
    })
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
        || input.loci.metagenomes.len() != input.document_count as usize
        || input
            .metadata
            .is_some_and(|metadata| !std::ptr::eq(metadata, input.loci) && metadata != input.loci)
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
    {
        let metadata = input.loci;
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
                    .checked_add(entry.original_contig_count as usize)
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
    let mut original_start = 0u32;
    for (document_id, metagenome) in metadata.metagenomes.iter().enumerate() {
        if !names.insert(metagenome.name.as_str())
            || metagenome.bgzf_bytes == 0
            || metagenome.bgzf_sha256 == [0; 32]
            || metagenome.original_contig_count == 0
            || metagenome.original_contig_start != original_start
            || !(1..=64).contains(&metagenome.locus_bits)
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
        put_u32(&mut record, 80, metagenome.original_contig_start);
        put_u32(&mut record, 84, metagenome.original_contig_count);
        record[88] = metagenome.locus_bits;
        documents.extend_from_slice(&record);
        let mut contig_names = HashSet::new();
        let mut previous_contig = None;
        let mut previous_end = 0;
        for contig in &metagenome.contigs {
            if !contig_names.insert(contig.name.as_str())
                || contig.length == 0
                || contig.line_bases == 0
                || contig.line_width < contig.line_bases
                || contig.line_width > contig.line_bases.saturating_add(2)
                || contig.local_contig >= metagenome.original_contig_count
                || previous_contig.is_some_and(|previous| previous >= contig.local_contig)
                || contig.fasta_offset < previous_end
            {
                return Err(OwnerWriteError::Invalid("contig metadata"));
            }
            let name = push_string(strings, &contig.name)?;
            let mut record = [0; OWNER_CONTIG_SIZE as usize];
            put_u32(&mut record, 0, document_id as u32);
            put_u32(&mut record, 4, name.0);
            put_u32(&mut record, 8, name.1);
            put_u32(
                &mut record,
                12,
                metagenome
                    .original_contig_start
                    .checked_add(contig.local_contig)
                    .ok_or(OwnerWriteError::Invalid("contig identity"))?,
            );
            put_u64(&mut record, 16, contig.length);
            put_u64(&mut record, 24, contig.fasta_offset);
            put_u32(&mut record, 32, contig.line_bases);
            put_u32(&mut record, 36, contig.line_width);
            contigs.extend_from_slice(&record);
            previous_contig = Some(contig.local_contig);
            let last_position = contig.length - 1;
            previous_end = (last_position / u64::from(contig.line_bases))
                .checked_mul(u64::from(contig.line_width))
                .and_then(|offset| offset.checked_add(last_position % u64::from(contig.line_bases)))
                .and_then(|offset| offset.checked_add(contig.fasta_offset))
                .and_then(|last| last.checked_add(1))
                .ok_or(OwnerWriteError::Invalid("contig locus"))?;
            if 64 - (previous_end - 1).leading_zeros() > u32::from(metagenome.locus_bits) {
                return Err(OwnerWriteError::Invalid("document locus width"));
            }
        }
        contig_start = contig_start
            .checked_add(
                u32::try_from(metagenome.contigs.len())
                    .map_err(|_| OwnerWriteError::Invalid("contig count"))?,
            )
            .ok_or(OwnerWriteError::Invalid("contig count"))?;
        gzi.extend_from_slice(&metagenome.gzi);
        original_start = original_start
            .checked_add(metagenome.original_contig_count)
            .ok_or(OwnerWriteError::Invalid("contig count"))?;
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
