//! Independent local Jam Index part data.

use super::manifest::ScreenSelectionPolicy;
use super::signature::{MetagenomeSignatureBuilder, SignatureSelectionError};
use crate::format::bucket_id;
use crate::jam_index::external::{ExternalError, ExternalSource, SequenceAccess, read_fai};
use crate::sequence::{
    EncodedContig, SequenceError, decode_ambiguity_payload, decode_range, encode_ambiguity_payload,
    encode_contig,
};
use memmap2::Mmap;
use needletail::{Sequence, parse_fastx_file};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{File, OpenOptions};
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

const MAGIC: [u8; 8] = *b"JAMIDX1P";
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 512;
#[allow(dead_code)]
const SOURCE_RECORD_SIZE: usize = 192;
const METAGENOME_RECORD_SIZE: usize = 96;
const CONTIG_RECORD_SIZE: usize = 128;
#[allow(dead_code)]
const POSTING_HEADER_SIZE: usize = 16;
const SIGNATURE_RECORD_SIZE: usize = 16;
const HEADER_CHECKSUM_OFFSET: usize = 280;

const SIGNATURE_CONTIG: u32 = 1;
const SIGNATURE_WHOLE_METAGENOME: u32 = 1 << 1;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeSource {
    pub metagenome_id: String,
    pub sequence_path: PathBuf,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PartScreenSample {
    pub metagenome_id: String,
    pub hashes: Vec<u64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PartWriteResult {
    pub screen_samples: Vec<PartScreenSample>,
    pub metagenome_count: u32,
    pub contig_count: u64,
    pub total_bases: u64,
    pub estimated_signature_count: u64,
    pub posting_count: u64,
    pub signature_record_count: u64,
    pub contig_signature_bytes: u64,
    pub source_reference_bytes: u64,
    pub packed_sequence_bytes: u64,
    pub data_file_bytes: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PartAccess {
    PlainFai,
    Bgzf,
    Sequential,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PartMetagenome {
    pub metagenome_id: String,
    pub first_contig: u32,
    pub contig_count: u32,
    pub total_bases: u64,
    pub screen_hash_count: u32,
    pub source_path: PathBuf,
    pub source_size: u64,
    pub source_sha256: [u8; 32],
    pub access: PartAccess,
    pub fai_path: Option<PathBuf>,
    pub fai_sha256: Option<[u8; 32]>,
    pub gzi_path: Option<PathBuf>,
    pub gzi_sha256: Option<[u8; 32]>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PartContig {
    pub contig_id: u32,
    pub metagenome_id: u32,
    pub name: String,
    pub base_count: u64,
    pub source_ordinal: u32,
    pub fai_offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
    pub sequence_sha256: [u8; 32],
    pub signature_count: u32,
    packed_offset: u64,
    packed_length: u64,
    ambiguity_offset: u64,
    ambiguity_length: u64,
    sequence_checksum: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SignatureHit {
    pub contig_id: u32,
    pub contig_selected: bool,
    pub whole_metagenome_selected: bool,
}

#[derive(Clone, Copy, Debug)]
struct Header {
    metagenome_count: u32,
    contig_count: u64,
    total_bases: u64,
    signature_count: u64,
    sequence_offset: u64,
    sequence_length: u64,
    signature_offset: u64,
    signature_length: u64,
    metagenome_offset: u64,
    metagenome_length: u64,
    contig_offset: u64,
    contig_length: u64,
    string_offset: u64,
    string_length: u64,
    sequence_checksum: [u8; 32],
    signature_checksum: [u8; 32],
    metagenome_checksum: [u8; 32],
    contig_checksum: [u8; 32],
    string_checksum: [u8; 32],
}

#[derive(Clone, Debug)]
struct MetagenomeBuildRecord {
    metagenome_id: u32,
    first_contig: u32,
    contig_count: u32,
    screen_hash_count: u32,
    id_offset: u64,
    id_length: u32,
    total_bases: u64,
    source_sha256: [u8; 32],
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
struct SourceBuildRecord {
    path_offset: u64,
    path_length: u32,
    fai_offset: u64,
    fai_length: u32,
    gzi_offset: u64,
    gzi_length: u32,
    source_size: u64,
    access: PartAccess,
    source_sha256: [u8; 32],
    fai_sha256: [u8; 32],
    gzi_sha256: [u8; 32],
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
struct ContigBuildRecord {
    contig_id: u32,
    metagenome_id: u32,
    name_offset: u64,
    name_length: u32,
    base_count: u64,
    source_ordinal: u32,
    fai_offset: u64,
    line_bases: u32,
    line_width: u32,
    sequence_sha256: [u8; 32],
    signature_count: u32,
    packed_offset: u64,
    packed_length: u64,
    ambiguity_offset: u64,
    ambiguity_length: u64,
    sequence_checksum: [u8; 32],
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct SignatureRecord {
    hash: u64,
    contig_id: u32,
    flags: u32,
}

#[allow(dead_code)]
#[derive(Clone, Copy, Debug)]
struct PostingHeader {
    offset: u64,
    count: u32,
    length: u32,
}

pub fn write_part(
    output: impl AsRef<Path>,
    sources: &[MetagenomeSource],
    policy: &ScreenSelectionPolicy,
) -> Result<PartWriteResult, JamIndexPartError> {
    policy
        .validate()
        .map_err(|error| JamIndexPartError::InvalidInput(error.to_string()))?;
    if sources.is_empty() {
        return Err(JamIndexPartError::InvalidInput(
            "a part requires at least one metagenome".to_string(),
        ));
    }
    let mut seen_ids = BTreeSet::new();
    if sources.iter().any(|source| {
        source.metagenome_id.trim().is_empty()
            || !source.sequence_path.is_file()
            || !seen_ids.insert(source.metagenome_id.clone())
    }) {
        return Err(JamIndexPartError::InvalidInput(
            "metagenome IDs must be unique and sequence paths must be local files".to_string(),
        ));
    }

    let output = output.as_ref();
    let mut file = OpenOptions::new()
        .create_new(true)
        .read(true)
        .write(true)
        .open(output)?;
    file.write_all(&vec![0u8; HEADER_SIZE])?;
    let sequence_offset = HEADER_SIZE as u64;
    let mut sequence_hasher = Sha256::new();
    let mut strings = Vec::new();
    let mut metagenomes = Vec::with_capacity(sources.len());
    let mut contigs = Vec::new();
    let mut signatures = Vec::new();
    let mut screen_samples = Vec::with_capacity(sources.len());
    let mut total_bases = 0u64;
    let mut packed_sequence_bytes = 0u64;
    let mut estimated_signature_count = 0u64;

    for (metagenome_index, source) in sources.iter().enumerate() {
        let metagenome_id =
            u32::try_from(metagenome_index).map_err(|_| JamIndexPartError::Overflow)?;
        let first_contig = u32::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let id_offset = append_string(&mut strings, &source.metagenome_id)?;
        let id_length =
            u32::try_from(source.metagenome_id.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let source_sha256 = sha256_file(&source.sequence_path)?;
        let mut selector = MetagenomeSignatureBuilder::new(policy.clone())?;
        let mut reader =
            parse_fastx_file(&source.sequence_path).map_err(|error| JamIndexPartError::Parse {
                path: source.sequence_path.clone(),
                message: error.to_string(),
            })?;
        let mut metagenome_bases = 0u64;
        let mut metagenome_contigs = 0u32;
        let mut metagenome_contig_ids = Vec::new();

        while let Some(record) = reader.next() {
            let record = record.map_err(|error| JamIndexPartError::Parse {
                path: source.sequence_path.clone(),
                message: error.to_string(),
            })?;
            let sequence = record.seq();
            let contig_id =
                u32::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let local_contig_id = metagenome_contigs;
            let name = std::str::from_utf8(record.id()).map_err(|_| {
                JamIndexPartError::InvalidInput(format!(
                    "contig name in {} is not UTF-8",
                    source.sequence_path.display()
                ))
            })?;
            let name_offset = append_string(&mut strings, name)?;
            let name_length = u32::try_from(name.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let selected = selector.add_contig(&sequence)?;
            estimated_signature_count =
                estimated_signature_count.saturating_add(u64::from(selected.requested_budget));
            for hash in selected.hashes {
                signatures.push(SignatureRecord {
                    hash,
                    contig_id,
                    flags: SIGNATURE_CONTIG,
                });
            }

            let encoded = encode_contig(&sequence)?;
            let ambiguity = encode_ambiguity_payload(&encoded.ambiguities)?;
            let packed_offset = file.stream_position()?;
            file.write_all(&encoded.two_bit)?;
            sequence_hasher.update(&encoded.two_bit);
            let packed_length =
                u64::try_from(encoded.two_bit.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let ambiguity_offset = file.stream_position()?;
            file.write_all(&ambiguity)?;
            sequence_hasher.update(&ambiguity);
            let ambiguity_length =
                u64::try_from(ambiguity.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let sequence_checksum = checksum_pair(&encoded.two_bit, &ambiguity);
            packed_sequence_bytes = packed_sequence_bytes.saturating_add(packed_length);
            metagenome_bases = metagenome_bases.saturating_add(encoded.base_count);
            total_bases = total_bases.saturating_add(encoded.base_count);
            metagenome_contigs = metagenome_contigs
                .checked_add(1)
                .ok_or(JamIndexPartError::Overflow)?;
            metagenome_contig_ids.push(contig_id);
            contigs.push(ContigBuildRecord {
                contig_id,
                metagenome_id,
                name_offset,
                name_length,
                base_count: encoded.base_count,
                source_ordinal: local_contig_id,
                fai_offset: 0,
                line_bases: 0,
                line_width: 0,
                sequence_sha256: sha256(&record.normalize(true)),
                signature_count: 0,
                packed_offset,
                packed_length,
                ambiguity_offset,
                ambiguity_length,
                sequence_checksum,
            });
            debug_assert_eq!(
                usize::try_from(local_contig_id).ok(),
                Some(metagenome_contig_ids.len() - 1)
            );
        }
        if metagenome_contigs == 0 {
            return Err(JamIndexPartError::InvalidInput(format!(
                "metagenome {} contains no contigs",
                source.metagenome_id
            )));
        }
        estimated_signature_count =
            estimated_signature_count.saturating_add(u64::from(policy.whole_metagenome_budget));
        let selected = selector.finish();
        for (hash, local_contig_id) in &selected.whole_hash_contigs {
            let contig_id = *metagenome_contig_ids
                .get(usize::try_from(*local_contig_id).map_err(|_| JamIndexPartError::Overflow)?)
                .ok_or(JamIndexPartError::Overflow)?;
            signatures.push(SignatureRecord {
                hash: *hash,
                contig_id,
                flags: SIGNATURE_WHOLE_METAGENOME,
            });
        }
        let screen_hash_count =
            u32::try_from(selected.union_hashes.len()).map_err(|_| JamIndexPartError::Overflow)?;
        screen_samples.push(PartScreenSample {
            metagenome_id: source.metagenome_id.clone(),
            hashes: selected.union_hashes,
        });
        metagenomes.push(MetagenomeBuildRecord {
            metagenome_id,
            first_contig,
            contig_count: metagenome_contigs,
            screen_hash_count,
            id_offset,
            id_length,
            total_bases: metagenome_bases,
            source_sha256,
        });
    }

    signatures.sort_unstable();
    let mut merged = Vec::<SignatureRecord>::with_capacity(signatures.len());
    for record in signatures {
        if let Some(previous) = merged.last_mut()
            && previous.hash == record.hash
            && previous.contig_id == record.contig_id
        {
            previous.flags |= record.flags;
        } else {
            merged.push(record);
        }
    }
    let mut signature_counts = vec![0u32; contigs.len()];
    for record in &merged {
        let count = signature_counts
            .get_mut(usize::try_from(record.contig_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::Overflow)?;
        *count = count.checked_add(1).ok_or(JamIndexPartError::Overflow)?;
    }
    for (contig, count) in contigs.iter_mut().zip(signature_counts) {
        contig.signature_count = count;
    }

    let sequence_end = align_file(&mut file, &mut sequence_hasher, 8)?;
    let sequence_length = sequence_end
        .checked_sub(sequence_offset)
        .ok_or(JamIndexPartError::Overflow)?;
    let signature_offset = sequence_end;
    let signature_bytes = encode_signatures(&merged);
    file.write_all(&signature_bytes)?;
    let metagenome_offset = align_file_plain(&mut file, 8)?;
    let metagenome_bytes = encode_metagenomes(&metagenomes);
    file.write_all(&metagenome_bytes)?;
    let contig_offset = align_file_plain(&mut file, 8)?;
    let contig_bytes = encode_contigs(&contigs);
    file.write_all(&contig_bytes)?;
    let string_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&strings)?;
    let object_size = file.stream_position()?;

    let header = Header {
        metagenome_count: u32::try_from(metagenomes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_count: u64::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?,
        total_bases,
        signature_count: u64::try_from(merged.len()).map_err(|_| JamIndexPartError::Overflow)?,
        sequence_offset,
        sequence_length,
        signature_offset,
        signature_length: u64::try_from(signature_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        metagenome_offset,
        metagenome_length: u64::try_from(metagenome_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_offset,
        contig_length: u64::try_from(contig_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        string_offset,
        string_length: u64::try_from(strings.len()).map_err(|_| JamIndexPartError::Overflow)?,
        sequence_checksum: sequence_hasher.finalize().into(),
        signature_checksum: sha256(&signature_bytes),
        metagenome_checksum: sha256(&metagenome_bytes),
        contig_checksum: sha256(&contig_bytes),
        string_checksum: sha256(&strings),
    };
    let header_bytes = encode_header(header);
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&header_bytes)?;
    file.sync_all()?;
    Ok(PartWriteResult {
        screen_samples,
        metagenome_count: header.metagenome_count,
        contig_count: header.contig_count,
        total_bases,
        estimated_signature_count,
        posting_count: header.signature_count,
        signature_record_count: header.signature_count,
        contig_signature_bytes: header.signature_length,
        source_reference_bytes: 0,
        packed_sequence_bytes,
        data_file_bytes: object_size,
    })
}

pub fn write_external_part(
    output: impl AsRef<Path>,
    sources: &[MetagenomeSource],
    policy: &ScreenSelectionPolicy,
) -> Result<PartWriteResult, JamIndexPartError> {
    policy
        .validate()
        .map_err(|error| JamIndexPartError::InvalidInput(error.to_string()))?;
    if sources.is_empty() {
        return Err(JamIndexPartError::InvalidInput(
            "a part requires at least one metagenome".to_string(),
        ));
    }
    let mut seen_ids = BTreeSet::new();
    if sources.iter().any(|source| {
        source.metagenome_id.trim().is_empty()
            || !source.sequence_path.is_file()
            || !seen_ids.insert(source.metagenome_id.clone())
    }) {
        return Err(JamIndexPartError::InvalidInput(
            "metagenome IDs must be unique and sequence paths must be local files".to_string(),
        ));
    }

    let mut strings = Vec::new();
    let mut source_records = Vec::with_capacity(sources.len());
    let mut metagenomes = Vec::with_capacity(sources.len());
    let mut contigs = Vec::new();
    let mut postings = BTreeMap::<(u64, u32), BTreeSet<u32>>::new();
    let mut screen_samples = Vec::with_capacity(sources.len());
    let mut total_bases = 0u64;
    let mut estimated_signature_count = 0u64;

    for (metagenome_index, source) in sources.iter().enumerate() {
        let metagenome_id =
            u32::try_from(metagenome_index).map_err(|_| JamIndexPartError::Overflow)?;
        let first_contig = u32::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let id_offset = append_string(&mut strings, &source.metagenome_id)?;
        let id_length =
            u32::try_from(source.metagenome_id.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let external = ExternalSource::detect(&source.sequence_path);
        let fai = external.fai_path.as_ref().map(read_fai).transpose()?;
        let (path_offset, path_length) = append_path(&mut strings, Some(&external.path))?;
        let (fai_offset, fai_length) = append_path(&mut strings, external.fai_path.as_deref())?;
        let (gzi_offset, gzi_length) = append_path(&mut strings, external.gzi_path.as_deref())?;
        let source_sha256 = sha256_file(&source.sequence_path)?;
        source_records.push(SourceBuildRecord {
            path_offset,
            path_length,
            fai_offset,
            fai_length,
            gzi_offset,
            gzi_length,
            source_size: std::fs::metadata(&source.sequence_path)?.len(),
            access: part_access(external.access),
            source_sha256,
            fai_sha256: external
                .fai_path
                .as_deref()
                .map(sha256_file)
                .transpose()?
                .unwrap_or([0; 32]),
            gzi_sha256: external
                .gzi_path
                .as_deref()
                .map(sha256_file)
                .transpose()?
                .unwrap_or([0; 32]),
        });
        let mut selector = MetagenomeSignatureBuilder::new(policy.clone())?;
        let mut reader =
            parse_fastx_file(&source.sequence_path).map_err(|error| JamIndexPartError::Parse {
                path: source.sequence_path.clone(),
                message: error.to_string(),
            })?;
        let mut metagenome_bases = 0u64;
        let mut metagenome_contigs = 0u32;
        let mut metagenome_contig_ids = Vec::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|error| JamIndexPartError::Parse {
                path: source.sequence_path.clone(),
                message: error.to_string(),
            })?;
            let sequence = record.normalize(true);
            let contig_id =
                u32::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let source_ordinal = metagenome_contigs;
            let name = std::str::from_utf8(record.id()).map_err(|_| {
                JamIndexPartError::InvalidInput(format!(
                    "contig name in {} is not UTF-8",
                    source.sequence_path.display()
                ))
            })?;
            let name_offset = append_string(&mut strings, name)?;
            let name_length = u32::try_from(name.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let length = u64::try_from(sequence.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let indexed = fai.as_ref().and_then(|records| {
                records.get(usize::try_from(source_ordinal).unwrap_or(usize::MAX))
            });
            if indexed.is_some_and(|indexed| indexed.name != name || indexed.length != length) {
                return Err(JamIndexPartError::InvalidInput(format!(
                    "FAI does not match {} contig {}",
                    source.sequence_path.display(),
                    name
                )));
            }
            let selected = selector.add_contig(&sequence)?;
            estimated_signature_count =
                estimated_signature_count.saturating_add(u64::from(selected.requested_budget));
            for hash in selected.hashes {
                postings
                    .entry((hash, metagenome_id))
                    .or_default()
                    .insert(contig_id);
            }
            let (fai_offset, line_bases, line_width) = indexed.map_or((0, 0, 0), |indexed| {
                (indexed.offset, indexed.line_bases, indexed.line_width)
            });
            metagenome_bases = metagenome_bases.saturating_add(length);
            total_bases = total_bases.saturating_add(length);
            metagenome_contigs = metagenome_contigs
                .checked_add(1)
                .ok_or(JamIndexPartError::Overflow)?;
            metagenome_contig_ids.push(contig_id);
            contigs.push(ContigBuildRecord {
                contig_id,
                metagenome_id,
                name_offset,
                name_length,
                base_count: length,
                source_ordinal,
                fai_offset,
                line_bases,
                line_width,
                sequence_sha256: sha256(&sequence),
                signature_count: 0,
                packed_offset: 0,
                packed_length: 0,
                ambiguity_offset: 0,
                ambiguity_length: 0,
                sequence_checksum: [0; 32],
            });
        }
        if metagenome_contigs == 0 {
            return Err(JamIndexPartError::InvalidInput(format!(
                "metagenome {} contains no contigs",
                source.metagenome_id
            )));
        }
        if fai.as_ref().is_some_and(|records| {
            records.len() != usize::try_from(metagenome_contigs).unwrap_or(usize::MAX)
        }) {
            return Err(JamIndexPartError::InvalidInput(format!(
                "FAI contig count does not match {}",
                source.sequence_path.display()
            )));
        }
        estimated_signature_count =
            estimated_signature_count.saturating_add(u64::from(policy.whole_metagenome_budget));
        let selected = selector.finish();
        for (hash, source_ordinal) in &selected.whole_hash_contigs {
            let contig_id = *metagenome_contig_ids
                .get(usize::try_from(*source_ordinal).map_err(|_| JamIndexPartError::Overflow)?)
                .ok_or(JamIndexPartError::Overflow)?;
            postings
                .entry((*hash, metagenome_id))
                .or_default()
                .insert(contig_id);
        }
        let screen_hash_count =
            u32::try_from(selected.union_hashes.len()).map_err(|_| JamIndexPartError::Overflow)?;
        screen_samples.push(PartScreenSample {
            metagenome_id: source.metagenome_id.clone(),
            hashes: selected.union_hashes,
        });
        metagenomes.push(MetagenomeBuildRecord {
            metagenome_id,
            first_contig,
            contig_count: metagenome_contigs,
            screen_hash_count,
            id_offset,
            id_length,
            total_bases: metagenome_bases,
            source_sha256,
        });
    }

    let mut entries = screen_samples
        .iter()
        .enumerate()
        .flat_map(|(sample_id, sample)| {
            sample.hashes.iter().map(move |hash| {
                (
                    bucket_id(*hash),
                    *hash,
                    u32::try_from(sample_id).unwrap_or(u32::MAX),
                )
            })
        })
        .collect::<Vec<_>>();
    entries.sort_unstable();
    let mut headers = Vec::with_capacity(entries.len());
    let mut payload = Vec::new();
    for (_, hash, sample_id) in entries {
        let contig_ids = postings
            .get(&(hash, sample_id))
            .ok_or(JamIndexPartError::Corrupt("screen posting binding"))?;
        let offset = u64::try_from(payload.len()).map_err(|_| JamIndexPartError::Overflow)?;
        encode_deltas(contig_ids, &mut payload)?;
        let length = u32::try_from(
            u64::try_from(payload.len())
                .map_err(|_| JamIndexPartError::Overflow)?
                .saturating_sub(offset),
        )
        .map_err(|_| JamIndexPartError::Overflow)?;
        headers.push(PostingHeader {
            offset,
            count: u32::try_from(contig_ids.len()).map_err(|_| JamIndexPartError::Overflow)?,
            length,
        });
    }
    let source_bytes = encode_sources(&source_records);
    let posting_bytes = encode_postings(&headers, &payload);
    let metagenome_bytes = encode_metagenomes(&metagenomes);
    let contig_bytes = encode_external_contigs(&contigs);
    let output = output.as_ref();
    let mut file = OpenOptions::new()
        .create_new(true)
        .read(true)
        .write(true)
        .open(output)?;
    file.write_all(&vec![0u8; HEADER_SIZE])?;
    let source_offset = file.stream_position()?;
    file.write_all(&source_bytes)?;
    let posting_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&posting_bytes)?;
    let metagenome_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&metagenome_bytes)?;
    let contig_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&contig_bytes)?;
    let string_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&strings)?;
    let object_size = file.stream_position()?;
    let header = Header {
        metagenome_count: u32::try_from(metagenomes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_count: u64::try_from(contigs.len()).map_err(|_| JamIndexPartError::Overflow)?,
        total_bases,
        signature_count: u64::try_from(headers.len()).map_err(|_| JamIndexPartError::Overflow)?,
        sequence_offset: source_offset,
        sequence_length: u64::try_from(source_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        signature_offset: posting_offset,
        signature_length: u64::try_from(posting_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        metagenome_offset,
        metagenome_length: u64::try_from(metagenome_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_offset,
        contig_length: u64::try_from(contig_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        string_offset,
        string_length: u64::try_from(strings.len()).map_err(|_| JamIndexPartError::Overflow)?,
        sequence_checksum: sha256(&source_bytes),
        signature_checksum: sha256(&posting_bytes),
        metagenome_checksum: sha256(&metagenome_bytes),
        contig_checksum: sha256(&contig_bytes),
        string_checksum: sha256(&strings),
    };
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&encode_header(header))?;
    file.sync_all()?;
    Ok(PartWriteResult {
        screen_samples,
        metagenome_count: header.metagenome_count,
        contig_count: header.contig_count,
        total_bases,
        estimated_signature_count,
        posting_count: header.signature_count,
        signature_record_count: 0,
        contig_signature_bytes: header.signature_length,
        source_reference_bytes: header.sequence_length,
        packed_sequence_bytes: 0,
        data_file_bytes: object_size,
    })
}

pub struct JamIndexPartReader {
    mmap: Mmap,
    header: Header,
    metagenomes: Vec<PartMetagenome>,
    contigs: Vec<PartContig>,
}

impl JamIndexPartReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, JamIndexPartError> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        if mmap.len() < HEADER_SIZE {
            return Err(JamIndexPartError::Corrupt("truncated header"));
        }
        let header = decode_header(&mmap[..HEADER_SIZE])?;
        validate_range(&mmap, header.sequence_offset, header.sequence_length)?;
        let signature = section(&mmap, header.signature_offset, header.signature_length)?;
        let metagenome = section(&mmap, header.metagenome_offset, header.metagenome_length)?;
        let contig = section(&mmap, header.contig_offset, header.contig_length)?;
        let strings = section(&mmap, header.string_offset, header.string_length)?;
        if sha256(signature) != header.signature_checksum
            || sha256(metagenome) != header.metagenome_checksum
            || sha256(contig) != header.contig_checksum
            || sha256(strings) != header.string_checksum
        {
            return Err(JamIndexPartError::Corrupt("section checksum mismatch"));
        }
        validate_section_lengths(header)?;
        let metagenomes = decode_metagenomes(metagenome, strings, header.metagenome_count)?;
        let contigs = decode_contigs(contig, strings, header.contig_count, mmap.len() as u64)?;
        validate_directories(&metagenomes, &contigs, header.total_bases)?;
        validate_signatures(signature, header.signature_count, header.contig_count)?;
        Ok(Self {
            mmap,
            header,
            metagenomes,
            contigs,
        })
    }

    #[must_use]
    pub fn metagenomes(&self) -> &[PartMetagenome] {
        &self.metagenomes
    }

    #[must_use]
    pub fn contigs(&self) -> &[PartContig] {
        &self.contigs
    }

    #[must_use]
    pub const fn total_bases(&self) -> u64 {
        self.header.total_bases
    }

    #[must_use]
    pub const fn signature_record_count(&self) -> u64 {
        self.header.signature_count
    }

    pub fn signature_hits(&self, hash: u64) -> Result<Vec<SignatureHit>, JamIndexPartError> {
        let mut hits = Vec::new();
        self.visit_signature_hits(hash, &mut |hit| hits.push(hit))?;
        Ok(hits)
    }

    pub fn visit_signature_hits(
        &self,
        hash: u64,
        visitor: &mut impl FnMut(SignatureHit),
    ) -> Result<(), JamIndexPartError> {
        let bytes = section(
            &self.mmap,
            self.header.signature_offset,
            self.header.signature_length,
        )?;
        let count = usize::try_from(self.header.signature_count)
            .map_err(|_| JamIndexPartError::Overflow)?;
        let mut left = 0usize;
        let mut right = count;
        while left < right {
            let middle = left + (right - left) / 2;
            if signature_hash(bytes, middle) < hash {
                left = middle + 1;
            } else {
                right = middle;
            }
        }
        while left < count && signature_hash(bytes, left) == hash {
            let record = decode_signature(bytes, left)?;
            visitor(SignatureHit {
                contig_id: record.contig_id,
                contig_selected: record.flags & SIGNATURE_CONTIG != 0,
                whole_metagenome_selected: record.flags & SIGNATURE_WHOLE_METAGENOME != 0,
            });
            left += 1;
        }
        Ok(())
    }

    pub fn metagenome_contigs(
        &self,
        metagenome_id: u32,
    ) -> Result<std::ops::Range<u32>, JamIndexPartError> {
        let metagenome = self
            .metagenomes
            .get(usize::try_from(metagenome_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::UnknownMetagenome(metagenome_id))?;
        Ok(metagenome.first_contig
            ..metagenome
                .first_contig
                .checked_add(metagenome.contig_count)
                .ok_or(JamIndexPartError::Overflow)?)
    }

    pub fn read_contig(&self, contig_id: u32) -> Result<Vec<u8>, JamIndexPartError> {
        let contig = self
            .contigs
            .get(usize::try_from(contig_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::UnknownContig(contig_id))?;
        let packed = section(&self.mmap, contig.packed_offset, contig.packed_length)?;
        let ambiguity = section(&self.mmap, contig.ambiguity_offset, contig.ambiguity_length)?;
        if checksum_pair(packed, ambiguity) != contig.sequence_checksum {
            return Err(JamIndexPartError::Corrupt("contig checksum mismatch"));
        }
        let encoded = EncodedContig {
            base_count: contig.base_count,
            two_bit: packed.to_vec(),
            ambiguities: decode_ambiguity_payload(ambiguity)?,
        };
        decode_range(&encoded, 0..contig.base_count).map_err(Into::into)
    }

    pub fn verify_sequence_section(&self) -> Result<(), JamIndexPartError> {
        let sequence = section(
            &self.mmap,
            self.header.sequence_offset,
            self.header.sequence_length,
        )?;
        if sha256(sequence) != self.header.sequence_checksum {
            return Err(JamIndexPartError::Corrupt(
                "sequence section checksum mismatch",
            ));
        }
        Ok(())
    }
}

fn append_string(strings: &mut Vec<u8>, value: &str) -> Result<u64, JamIndexPartError> {
    let offset = u64::try_from(strings.len()).map_err(|_| JamIndexPartError::Overflow)?;
    strings.extend_from_slice(value.as_bytes());
    Ok(offset)
}

fn append_path(
    strings: &mut Vec<u8>,
    path: Option<&Path>,
) -> Result<(u64, u32), JamIndexPartError> {
    let Some(path) = path else {
        return Ok((0, 0));
    };
    let value = path
        .to_str()
        .ok_or_else(|| JamIndexPartError::InvalidInput("source path is not UTF-8".to_string()))?;
    Ok((
        append_string(strings, value)?,
        u32::try_from(value.len()).map_err(|_| JamIndexPartError::Overflow)?,
    ))
}

fn part_access(access: SequenceAccess) -> PartAccess {
    match access {
        SequenceAccess::PlainFai => PartAccess::PlainFai,
        SequenceAccess::Bgzf => PartAccess::Bgzf,
        SequenceAccess::Sequential => PartAccess::Sequential,
    }
}

fn encode_sources(records: &[SourceBuildRecord]) -> Vec<u8> {
    let mut bytes = vec![0u8; records.len() * SOURCE_RECORD_SIZE];
    for (index, record) in records.iter().enumerate() {
        let offset = index * SOURCE_RECORD_SIZE;
        put_u64(&mut bytes, offset, record.path_offset);
        put_u32(&mut bytes, offset + 8, record.path_length);
        put_u32(
            &mut bytes,
            offset + 12,
            match record.access {
                PartAccess::PlainFai => 1,
                PartAccess::Bgzf => 2,
                PartAccess::Sequential => 3,
            },
        );
        put_u64(&mut bytes, offset + 16, record.fai_offset);
        put_u32(&mut bytes, offset + 24, record.fai_length);
        put_u64(&mut bytes, offset + 32, record.gzi_offset);
        put_u32(&mut bytes, offset + 40, record.gzi_length);
        put_u64(&mut bytes, offset + 48, record.source_size);
        bytes[offset + 56..offset + 88].copy_from_slice(&record.source_sha256);
        bytes[offset + 88..offset + 120].copy_from_slice(&record.fai_sha256);
        bytes[offset + 120..offset + 152].copy_from_slice(&record.gzi_sha256);
    }
    bytes
}

fn encode_external_contigs(records: &[ContigBuildRecord]) -> Vec<u8> {
    let mut bytes = vec![0u8; records.len() * CONTIG_RECORD_SIZE];
    for (index, record) in records.iter().enumerate() {
        let offset = index * CONTIG_RECORD_SIZE;
        put_u32(&mut bytes, offset, record.contig_id);
        put_u32(&mut bytes, offset + 4, record.metagenome_id);
        put_u64(&mut bytes, offset + 8, record.name_offset);
        put_u32(&mut bytes, offset + 16, record.name_length);
        put_u32(&mut bytes, offset + 20, record.source_ordinal);
        put_u64(&mut bytes, offset + 24, record.base_count);
        put_u64(&mut bytes, offset + 32, record.fai_offset);
        put_u32(&mut bytes, offset + 40, record.line_bases);
        put_u32(&mut bytes, offset + 44, record.line_width);
        bytes[offset + 48..offset + 80].copy_from_slice(&record.sequence_sha256);
    }
    bytes
}

fn encode_postings(headers: &[PostingHeader], payload: &[u8]) -> Vec<u8> {
    let header_bytes = headers.len() * POSTING_HEADER_SIZE;
    let mut bytes = vec![0u8; header_bytes + payload.len()];
    for (index, header) in headers.iter().enumerate() {
        let offset = index * POSTING_HEADER_SIZE;
        put_u64(&mut bytes, offset, header.offset);
        put_u32(&mut bytes, offset + 8, header.count);
        put_u32(&mut bytes, offset + 12, header.length);
    }
    bytes[header_bytes..].copy_from_slice(payload);
    bytes
}

fn encode_deltas(
    contig_ids: &BTreeSet<u32>,
    payload: &mut Vec<u8>,
) -> Result<(), JamIndexPartError> {
    let mut previous = 0u32;
    for (index, contig_id) in contig_ids.iter().copied().enumerate() {
        let mut delta = if index == 0 {
            contig_id
        } else {
            contig_id
                .checked_sub(previous)
                .ok_or(JamIndexPartError::Overflow)?
        };
        loop {
            let mut byte = u8::try_from(delta & 0x7f).expect("seven-bit varint chunk");
            delta >>= 7;
            if delta != 0 {
                byte |= 0x80;
            }
            payload.push(byte);
            if delta == 0 {
                break;
            }
        }
        previous = contig_id;
    }
    Ok(())
}

fn align_file(
    file: &mut File,
    hasher: &mut Sha256,
    alignment: u64,
) -> Result<u64, JamIndexPartError> {
    let offset = file.stream_position()?;
    let padding = (alignment - offset % alignment) % alignment;
    let zeros = vec![0u8; usize::try_from(padding).map_err(|_| JamIndexPartError::Overflow)?];
    file.write_all(&zeros)?;
    hasher.update(&zeros);
    file.stream_position().map_err(Into::into)
}

fn align_file_plain(file: &mut File, alignment: u64) -> Result<u64, JamIndexPartError> {
    let offset = file.stream_position()?;
    let padding = (alignment - offset % alignment) % alignment;
    file.write_all(&vec![
        0u8;
        usize::try_from(padding)
            .map_err(|_| JamIndexPartError::Overflow)?
    ])?;
    file.stream_position().map_err(Into::into)
}

fn encode_header(header: Header) -> [u8; HEADER_SIZE] {
    let mut bytes = [0u8; HEADER_SIZE];
    bytes[..8].copy_from_slice(&MAGIC);
    put_u16(&mut bytes, 8, VERSION);
    put_u16(&mut bytes, 10, HEADER_SIZE as u16);
    put_u32(&mut bytes, 16, header.metagenome_count);
    put_u64(&mut bytes, 24, header.contig_count);
    put_u64(&mut bytes, 32, header.total_bases);
    put_u64(&mut bytes, 40, header.signature_count);
    put_u64(&mut bytes, 48, header.sequence_offset);
    put_u64(&mut bytes, 56, header.sequence_length);
    put_u64(&mut bytes, 64, header.signature_offset);
    put_u64(&mut bytes, 72, header.signature_length);
    put_u64(&mut bytes, 80, header.metagenome_offset);
    put_u64(&mut bytes, 88, header.metagenome_length);
    put_u64(&mut bytes, 96, header.contig_offset);
    put_u64(&mut bytes, 104, header.contig_length);
    put_u64(&mut bytes, 112, header.string_offset);
    put_u64(&mut bytes, 120, header.string_length);
    bytes[128..160].copy_from_slice(&header.sequence_checksum);
    bytes[160..192].copy_from_slice(&header.signature_checksum);
    bytes[192..224].copy_from_slice(&header.metagenome_checksum);
    bytes[224..256].copy_from_slice(&header.contig_checksum);
    bytes[256..280].copy_from_slice(&header.string_checksum[..24]);
    bytes[312..320].copy_from_slice(&header.string_checksum[24..]);
    let checksum = sha256(&bytes);
    bytes[HEADER_CHECKSUM_OFFSET..HEADER_CHECKSUM_OFFSET + 32].copy_from_slice(&checksum);
    bytes
}

fn decode_header(bytes: &[u8]) -> Result<Header, JamIndexPartError> {
    if bytes.len() != HEADER_SIZE
        || bytes[..8] != MAGIC
        || get_u16(bytes, 8) != VERSION
        || usize::from(get_u16(bytes, 10)) != HEADER_SIZE
        || bytes[12..16].iter().any(|byte| *byte != 0)
        || bytes[320..].iter().any(|byte| *byte != 0)
    {
        return Err(JamIndexPartError::Corrupt("invalid part header"));
    }
    let mut checked = [0u8; HEADER_SIZE];
    checked.copy_from_slice(bytes);
    let expected = array_32(bytes, HEADER_CHECKSUM_OFFSET);
    checked[HEADER_CHECKSUM_OFFSET..HEADER_CHECKSUM_OFFSET + 32].fill(0);
    if sha256(&checked) != expected {
        return Err(JamIndexPartError::Corrupt("header checksum mismatch"));
    }
    let mut string_checksum = [0u8; 32];
    string_checksum[..24].copy_from_slice(&bytes[256..280]);
    string_checksum[24..].copy_from_slice(&bytes[312..320]);
    Ok(Header {
        metagenome_count: get_u32(bytes, 16),
        contig_count: get_u64(bytes, 24),
        total_bases: get_u64(bytes, 32),
        signature_count: get_u64(bytes, 40),
        sequence_offset: get_u64(bytes, 48),
        sequence_length: get_u64(bytes, 56),
        signature_offset: get_u64(bytes, 64),
        signature_length: get_u64(bytes, 72),
        metagenome_offset: get_u64(bytes, 80),
        metagenome_length: get_u64(bytes, 88),
        contig_offset: get_u64(bytes, 96),
        contig_length: get_u64(bytes, 104),
        string_offset: get_u64(bytes, 112),
        string_length: get_u64(bytes, 120),
        sequence_checksum: array_32(bytes, 128),
        signature_checksum: array_32(bytes, 160),
        metagenome_checksum: array_32(bytes, 192),
        contig_checksum: array_32(bytes, 224),
        string_checksum,
    })
}

fn encode_metagenomes(records: &[MetagenomeBuildRecord]) -> Vec<u8> {
    let mut bytes = vec![0u8; records.len() * METAGENOME_RECORD_SIZE];
    for (index, record) in records.iter().enumerate() {
        let offset = index * METAGENOME_RECORD_SIZE;
        put_u32(&mut bytes, offset, record.metagenome_id);
        put_u32(&mut bytes, offset + 4, record.first_contig);
        put_u32(&mut bytes, offset + 8, record.contig_count);
        put_u32(&mut bytes, offset + 12, record.screen_hash_count);
        put_u64(&mut bytes, offset + 16, record.id_offset);
        put_u32(&mut bytes, offset + 24, record.id_length);
        put_u64(&mut bytes, offset + 32, record.total_bases);
        bytes[offset + 40..offset + 72].copy_from_slice(&record.source_sha256);
    }
    bytes
}

fn encode_contigs(records: &[ContigBuildRecord]) -> Vec<u8> {
    let mut bytes = vec![0u8; records.len() * CONTIG_RECORD_SIZE];
    for (index, record) in records.iter().enumerate() {
        let offset = index * CONTIG_RECORD_SIZE;
        put_u32(&mut bytes, offset, record.contig_id);
        put_u32(&mut bytes, offset + 4, record.metagenome_id);
        put_u64(&mut bytes, offset + 8, record.name_offset);
        put_u32(&mut bytes, offset + 16, record.name_length);
        put_u32(&mut bytes, offset + 20, record.signature_count);
        put_u64(&mut bytes, offset + 24, record.base_count);
        put_u64(&mut bytes, offset + 32, record.packed_offset);
        put_u64(&mut bytes, offset + 40, record.packed_length);
        put_u64(&mut bytes, offset + 48, record.ambiguity_offset);
        put_u64(&mut bytes, offset + 56, record.ambiguity_length);
        bytes[offset + 64..offset + 96].copy_from_slice(&record.sequence_checksum);
    }
    bytes
}

fn encode_signatures(records: &[SignatureRecord]) -> Vec<u8> {
    let mut bytes = vec![0u8; records.len() * SIGNATURE_RECORD_SIZE];
    for (index, record) in records.iter().enumerate() {
        let offset = index * SIGNATURE_RECORD_SIZE;
        put_u64(&mut bytes, offset, record.hash);
        put_u32(&mut bytes, offset + 8, record.contig_id);
        put_u32(&mut bytes, offset + 12, record.flags);
    }
    bytes
}

fn decode_metagenomes(
    bytes: &[u8],
    strings: &[u8],
    count: u32,
) -> Result<Vec<PartMetagenome>, JamIndexPartError> {
    if bytes.len()
        != usize::try_from(count).map_err(|_| JamIndexPartError::Overflow)? * METAGENOME_RECORD_SIZE
    {
        return Err(JamIndexPartError::Corrupt("metagenome directory length"));
    }
    (0..count)
        .map(|id| {
            let offset = usize::try_from(id).map_err(|_| JamIndexPartError::Overflow)?
                * METAGENOME_RECORD_SIZE;
            if get_u32(bytes, offset) != id
                || bytes[offset + 28..offset + 32]
                    .iter()
                    .any(|byte| *byte != 0)
                || bytes[offset + 72..offset + 96]
                    .iter()
                    .any(|byte| *byte != 0)
            {
                return Err(JamIndexPartError::Corrupt("metagenome directory record"));
            }
            let mut source_sha256 = [0u8; 32];
            source_sha256.copy_from_slice(&bytes[offset + 40..offset + 72]);
            Ok(PartMetagenome {
                metagenome_id: decode_string(
                    strings,
                    get_u64(bytes, offset + 16),
                    get_u32(bytes, offset + 24),
                )?,
                first_contig: get_u32(bytes, offset + 4),
                contig_count: get_u32(bytes, offset + 8),
                screen_hash_count: get_u32(bytes, offset + 12),
                total_bases: get_u64(bytes, offset + 32),
                source_path: PathBuf::new(),
                source_size: 0,
                source_sha256,
                access: PartAccess::Sequential,
                fai_path: None,
                fai_sha256: None,
                gzi_path: None,
                gzi_sha256: None,
            })
        })
        .collect()
}

fn decode_contigs(
    bytes: &[u8],
    strings: &[u8],
    count: u64,
    object_size: u64,
) -> Result<Vec<PartContig>, JamIndexPartError> {
    let count = usize::try_from(count).map_err(|_| JamIndexPartError::Overflow)?;
    if bytes.len() != count * CONTIG_RECORD_SIZE {
        return Err(JamIndexPartError::Corrupt("contig directory length"));
    }
    (0..count)
        .map(|id| {
            let offset = id * CONTIG_RECORD_SIZE;
            let contig_id = get_u32(bytes, offset);
            let packed_offset = get_u64(bytes, offset + 32);
            let packed_length = get_u64(bytes, offset + 40);
            let ambiguity_offset = get_u64(bytes, offset + 48);
            let ambiguity_length = get_u64(bytes, offset + 56);
            if usize::try_from(contig_id).ok() != Some(id)
                || bytes[offset + 96..offset + 128]
                    .iter()
                    .any(|byte| *byte != 0)
                || packed_offset
                    .checked_add(packed_length)
                    .is_none_or(|end| end > object_size)
                || ambiguity_offset
                    .checked_add(ambiguity_length)
                    .is_none_or(|end| end > object_size)
            {
                return Err(JamIndexPartError::Corrupt("contig directory record"));
            }
            Ok(PartContig {
                contig_id,
                metagenome_id: get_u32(bytes, offset + 4),
                name: decode_string(
                    strings,
                    get_u64(bytes, offset + 8),
                    get_u32(bytes, offset + 16),
                )?,
                base_count: get_u64(bytes, offset + 24),
                source_ordinal: 0,
                fai_offset: 0,
                line_bases: 0,
                line_width: 0,
                sequence_sha256: array_32(bytes, offset + 64),
                signature_count: get_u32(bytes, offset + 20),
                packed_offset,
                packed_length,
                ambiguity_offset,
                ambiguity_length,
                sequence_checksum: array_32(bytes, offset + 64),
            })
        })
        .collect()
}

fn validate_directories(
    metagenomes: &[PartMetagenome],
    contigs: &[PartContig],
    total_bases: u64,
) -> Result<(), JamIndexPartError> {
    let observed_bases = contigs.iter().map(|contig| contig.base_count).sum::<u64>();
    if observed_bases != total_bases {
        return Err(JamIndexPartError::Corrupt("total base count mismatch"));
    }
    for (metagenome_id, metagenome) in metagenomes.iter().enumerate() {
        let end = metagenome
            .first_contig
            .checked_add(metagenome.contig_count)
            .ok_or(JamIndexPartError::Overflow)?;
        let range = usize::try_from(metagenome.first_contig)
            .map_err(|_| JamIndexPartError::Overflow)?
            ..usize::try_from(end).map_err(|_| JamIndexPartError::Overflow)?;
        let selected = contigs
            .get(range)
            .ok_or(JamIndexPartError::Corrupt("metagenome contig range"))?;
        if selected
            .iter()
            .any(|contig| usize::try_from(contig.metagenome_id).ok() != Some(metagenome_id))
            || selected.iter().map(|contig| contig.base_count).sum::<u64>()
                != metagenome.total_bases
        {
            return Err(JamIndexPartError::Corrupt("metagenome contig binding"));
        }
    }
    Ok(())
}

fn validate_signatures(
    bytes: &[u8],
    count: u64,
    contig_count: u64,
) -> Result<(), JamIndexPartError> {
    let count = usize::try_from(count).map_err(|_| JamIndexPartError::Overflow)?;
    if bytes.len() != count * SIGNATURE_RECORD_SIZE {
        return Err(JamIndexPartError::Corrupt("signature table length"));
    }
    let mut previous = None;
    for index in 0..count {
        let record = decode_signature(bytes, index)?;
        if record.hash == 0
            || u64::from(record.contig_id) >= contig_count
            || record.flags == 0
            || record.flags & !(SIGNATURE_CONTIG | SIGNATURE_WHOLE_METAGENOME) != 0
            || previous.is_some_and(|old| old >= (record.hash, record.contig_id))
        {
            return Err(JamIndexPartError::Corrupt("signature table ordering"));
        }
        previous = Some((record.hash, record.contig_id));
    }
    Ok(())
}

fn validate_section_lengths(header: Header) -> Result<(), JamIndexPartError> {
    if header.metagenome_length
        != u64::from(header.metagenome_count) * METAGENOME_RECORD_SIZE as u64
        || header.contig_length != header.contig_count * CONTIG_RECORD_SIZE as u64
        || header.signature_length != header.signature_count * SIGNATURE_RECORD_SIZE as u64
    {
        return Err(JamIndexPartError::Corrupt("section length mismatch"));
    }
    Ok(())
}

fn decode_signature(bytes: &[u8], index: usize) -> Result<SignatureRecord, JamIndexPartError> {
    let offset = index
        .checked_mul(SIGNATURE_RECORD_SIZE)
        .ok_or(JamIndexPartError::Overflow)?;
    let record = bytes
        .get(offset..offset + SIGNATURE_RECORD_SIZE)
        .ok_or(JamIndexPartError::Corrupt("signature record"))?;
    Ok(SignatureRecord {
        hash: get_u64(record, 0),
        contig_id: get_u32(record, 8),
        flags: get_u32(record, 12),
    })
}

fn signature_hash(bytes: &[u8], index: usize) -> u64 {
    let offset = index * SIGNATURE_RECORD_SIZE;
    get_u64(bytes, offset)
}

fn decode_string(strings: &[u8], offset: u64, length: u32) -> Result<String, JamIndexPartError> {
    let start = usize::try_from(offset).map_err(|_| JamIndexPartError::Overflow)?;
    let end = start
        .checked_add(length as usize)
        .ok_or(JamIndexPartError::Overflow)?;
    std::str::from_utf8(
        strings
            .get(start..end)
            .ok_or(JamIndexPartError::Corrupt("string range"))?,
    )
    .map(str::to_string)
    .map_err(|_| JamIndexPartError::Corrupt("string UTF-8"))
}

fn section(bytes: &[u8], offset: u64, length: u64) -> Result<&[u8], JamIndexPartError> {
    let start = usize::try_from(offset).map_err(|_| JamIndexPartError::Overflow)?;
    let end = start
        .checked_add(usize::try_from(length).map_err(|_| JamIndexPartError::Overflow)?)
        .ok_or(JamIndexPartError::Overflow)?;
    bytes
        .get(start..end)
        .ok_or(JamIndexPartError::Corrupt("section range"))
}

fn validate_range(bytes: &[u8], offset: u64, length: u64) -> Result<(), JamIndexPartError> {
    section(bytes, offset, length).map(|_| ())
}

fn sha256_file(path: &Path) -> Result<[u8; 32], JamIndexPartError> {
    let mut reader = BufReader::with_capacity(1024 * 1024, File::open(path)?);
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        digest.update(&buffer[..read]);
    }
    Ok(digest.finalize().into())
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    Sha256::digest(bytes).into()
}

fn checksum_pair(left: &[u8], right: &[u8]) -> [u8; 32] {
    let mut digest = Sha256::new();
    digest.update(left);
    digest.update(right);
    digest.finalize().into()
}

fn put_u16(bytes: &mut [u8], offset: usize, value: u16) {
    bytes[offset..offset + 2].copy_from_slice(&value.to_le_bytes());
}

fn put_u32(bytes: &mut [u8], offset: usize, value: u32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

fn put_u64(bytes: &mut [u8], offset: usize, value: u64) {
    bytes[offset..offset + 8].copy_from_slice(&value.to_le_bytes());
}

fn get_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().expect("u16 field"))
}

fn get_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("u32 field"))
}

fn get_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("u64 field"))
}

fn array_32(bytes: &[u8], offset: usize) -> [u8; 32] {
    bytes[offset..offset + 32]
        .try_into()
        .expect("checksum field")
}

#[derive(Debug, Error)]
pub enum JamIndexPartError {
    #[error("Jam Index part I/O failed: {0}")]
    Io(#[from] std::io::Error),
    #[error("Jam Index part sequence failed: {0}")]
    Sequence(#[from] SequenceError),
    #[error("Jam Index external sequence failed: {0}")]
    External(String),
    #[error(transparent)]
    Signature(#[from] SignatureSelectionError),
    #[error("Jam Index part parse failed for {path}: {message}")]
    Parse { path: PathBuf, message: String },
    #[error("invalid Jam Index part input: {0}")]
    InvalidInput(String),
    #[error("corrupt Jam Index part: {0}")]
    Corrupt(&'static str),
    #[error("unknown Jam Index metagenome {0}")]
    UnknownMetagenome(u32),
    #[error("unknown Jam Index contig {0}")]
    UnknownContig(u32),
    #[error("Jam Index part coordinate overflow")]
    Overflow,
}

impl From<ExternalError> for JamIndexPartError {
    fn from(error: ExternalError) -> Self {
        Self::External(error.to_string())
    }
}
