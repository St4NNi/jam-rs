//! Independent local Jam Index part data.

use super::manifest::ScreenSelectionPolicy;
use super::signature::{MetagenomeSignatureBuilder, SignatureSelectionError};
use crate::format::{BUCKET_COUNT, bucket_id};
use crate::jam_index::external::{ContigRequest, read_selected as read_external};
use crate::jam_index::external::{ExternalError, ExternalSource, SequenceAccess, read_fai};
use memmap2::Mmap;
use needletail::{Sequence, parse_fastx_file, parse_fastx_reader};
use sha2::{Digest, Sha256};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use thiserror::Error;

const MAGIC: [u8; 8] = *b"JAMIDX2P";
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 512;
const SOURCE_RECORD_SIZE: usize = 160;
const METAGENOME_RECORD_SIZE: usize = 48;
const CONTIG_LENGTH_SIZE: usize = 4;
const EXCEPTIONAL_LENGTH_RECORD_SIZE: usize = 16;
const MAPPING_RECORD_SIZE: usize = 8;
const MAPPING_OVERFLOW: u32 = 1 << 31;
const MAPPING_WHOLE_SAMPLE: u32 = 1 << 30;
const MAPPING_SPATIAL: u32 = 1 << 29;
const MAPPING_COUNT_MASK: u32 = (1 << 29) - 1;
const SIGNATURE_RUN_RECORD_SIZE: usize = 24;
const SIGNATURE_RUN_RECORD_LIMIT: usize = 1_000_000;
const HEADER_CHECKSUM_OFFSET: usize = 280;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeSource {
    pub metagenome_id: String,
    pub sequence_path: PathBuf,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PublishedMetagenomeSource {
    pub metagenome_id: String,
    pub sequence_path: PathBuf,
    pub source_size: u64,
    pub source_sha256: [u8; 32],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StagedMetagenomeSource {
    pub metagenome_id: String,
    pub staged_sequence_path: PathBuf,
    pub published_sequence_path: PathBuf,
    pub expected_source_size: u64,
    pub expected_source_sha256: Option<[u8; 32]>,
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
    pub contig_posting_bytes: u64,
    pub source_reference_bytes: u64,
    pub metagenome_directory_bytes: u64,
    pub contig_length_bytes: u64,
    pub exceptional_length_bytes: u64,
    pub string_table_bytes: u64,
    pub data_file_bytes: u64,
    pub contig_signature_histogram: BTreeMap<u32, u64>,
    pub single_contig_mappings: u64,
    pub overflow_mappings: u64,
    pub overflow_contigs: u64,
    pub maximum_overflow_count: u32,
    pub signature_run_count: u32,
    pub signature_run_record_limit: u64,
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
    pub first_contig: u64,
    pub contig_count: u32,
    pub total_bases: u64,
    pub screen_hash_count: u32,
    pub exceptional_length_start: u32,
    pub exceptional_length_count: u32,
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
pub struct LoadedPartContig {
    pub name: String,
    pub bases: Arc<[u8]>,
}

#[derive(Clone, Copy, Debug)]
struct Header {
    format_version: u16,
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
    first_contig: u64,
    contig_count: u32,
    screen_hash_count: u32,
    id_offset: u64,
    id_length: u32,
    exceptional_start: u32,
    exceptional_count: u32,
    total_bases: u64,
}

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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ExceptionalLengthRecord {
    ordinal: u64,
    length: u64,
}

#[derive(Clone, Debug, Default)]
struct PostingBuildRecord {
    contig_ids: BTreeSet<u32>,
    spatial: bool,
    whole_sample: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PostingKind {
    pub spatial: bool,
    pub whole_sample: bool,
}

struct DecodedMapping {
    contig_ids: Vec<u32>,
    overflow_range: Option<std::ops::Range<usize>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SignatureRunRecord {
    hash: u64,
    sample_id: u32,
    contig_id: u32,
    spatial: bool,
    whole_sample: bool,
}

impl SignatureRunRecord {
    fn key(self) -> (usize, u64, u32, u32) {
        (
            bucket_id(self.hash),
            self.hash,
            self.sample_id,
            self.contig_id,
        )
    }
}

struct SignatureRunBuilder {
    directory: tempfile::TempDir,
    records: Vec<SignatureRunRecord>,
    runs: Vec<PathBuf>,
}

struct StreamingMappings {
    _directory: tempfile::TempDir,
    records_path: PathBuf,
    payload_path: PathBuf,
    signature_count: u64,
    stats: MappingStats,
    run_count: u32,
}

impl SignatureRunBuilder {
    fn new(parent: &Path) -> Result<Self, JamIndexPartError> {
        Ok(Self {
            directory: tempfile::Builder::new()
                .prefix(".jam-index-signature-runs-")
                .tempdir_in(parent)?,
            records: Vec::with_capacity(SIGNATURE_RUN_RECORD_LIMIT),
            runs: Vec::new(),
        })
    }

    fn push(&mut self, record: SignatureRunRecord) -> Result<(), JamIndexPartError> {
        self.records.push(record);
        if self.records.len() >= SIGNATURE_RUN_RECORD_LIMIT {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), JamIndexPartError> {
        if self.records.is_empty() {
            return Ok(());
        }
        self.records.sort_unstable_by_key(|record| record.key());
        let path = self
            .directory
            .path()
            .join(format!("run-{:06}.bin", self.runs.len()));
        let mut writer = BufWriter::new(
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(&path)?,
        );
        for record in self.records.drain(..) {
            write_signature_run_record(&mut writer, record)?;
        }
        writer.flush()?;
        self.runs.push(path);
        Ok(())
    }

    fn finish(
        mut self,
        screen_samples: &mut [PartScreenSample],
    ) -> Result<StreamingMappings, JamIndexPartError> {
        self.flush()?;
        let run_count = u32::try_from(self.runs.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let records_path = self.directory.path().join("mapping-records.bin");
        let payload_path = self.directory.path().join("mapping-payload.bin");
        let mut mapping_records = BufWriter::new(
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(&records_path)?,
        );
        let mut mapping_payload = BufWriter::new(
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(&payload_path)?,
        );
        let mut readers = self
            .runs
            .iter()
            .map(|path| File::open(path).map(BufReader::new))
            .collect::<Result<Vec<_>, _>>()?;
        let mut heap = BinaryHeap::new();
        for (run_index, reader) in readers.iter_mut().enumerate() {
            if let Some(record) = read_signature_run_record(reader)? {
                heap.push(Reverse((
                    record.key(),
                    record.spatial,
                    record.whole_sample,
                    run_index,
                )));
            }
        }
        let mut current_key = None;
        let mut current = PostingBuildRecord::default();
        let mut signature_count = 0u64;
        let mut payload_bytes = 0u64;
        let mut stats = MappingStats::default();
        while let Some(Reverse((key, spatial, whole_sample, run_index))) = heap.pop() {
            let (_, hash, sample_id, contig_id) = key;
            if current_key != Some((hash, sample_id)) {
                if let Some((previous_hash, previous_sample)) = current_key {
                    write_stream_mapping(
                        &mut mapping_records,
                        &mut mapping_payload,
                        &mut payload_bytes,
                        &mut stats,
                        screen_samples,
                        previous_hash,
                        previous_sample,
                        &current,
                    )?;
                    signature_count = signature_count.saturating_add(1);
                    current = PostingBuildRecord::default();
                }
                current_key = Some((hash, sample_id));
            }
            current.spatial |= spatial;
            current.whole_sample |= whole_sample;
            current.contig_ids.insert(contig_id);
            if let Some(next) = read_signature_run_record(&mut readers[run_index])? {
                heap.push(Reverse((
                    next.key(),
                    next.spatial,
                    next.whole_sample,
                    run_index,
                )));
            }
        }
        if let Some((hash, sample_id)) = current_key {
            write_stream_mapping(
                &mut mapping_records,
                &mut mapping_payload,
                &mut payload_bytes,
                &mut stats,
                screen_samples,
                hash,
                sample_id,
                &current,
            )?;
            signature_count = signature_count.saturating_add(1);
        }
        mapping_records.flush()?;
        mapping_payload.flush()?;
        Ok(StreamingMappings {
            _directory: self.directory,
            records_path,
            payload_path,
            signature_count,
            stats,
            run_count,
        })
    }
}

fn write_signature_run_record(
    writer: &mut impl Write,
    record: SignatureRunRecord,
) -> Result<(), JamIndexPartError> {
    let mut bytes = [0u8; SIGNATURE_RUN_RECORD_SIZE];
    put_u64(&mut bytes, 0, record.hash);
    put_u32(&mut bytes, 8, record.sample_id);
    put_u32(&mut bytes, 12, record.contig_id);
    put_u32(
        &mut bytes,
        16,
        u32::from(record.spatial) | (u32::from(record.whole_sample) << 1),
    );
    writer.write_all(&bytes)?;
    Ok(())
}

fn read_signature_run_record(
    reader: &mut impl Read,
) -> Result<Option<SignatureRunRecord>, JamIndexPartError> {
    let mut bytes = [0u8; SIGNATURE_RUN_RECORD_SIZE];
    let mut read = 0usize;
    while read < bytes.len() {
        let count = reader.read(&mut bytes[read..])?;
        if count == 0 {
            if read == 0 {
                return Ok(None);
            }
            return Err(JamIndexPartError::Corrupt("truncated signature run"));
        }
        read += count;
    }
    if bytes[20..].iter().any(|byte| *byte != 0) || get_u32(&bytes, 16) & !3 != 0 {
        return Err(JamIndexPartError::Corrupt("signature run record"));
    }
    let flags = get_u32(&bytes, 16);
    Ok(Some(SignatureRunRecord {
        hash: get_u64(&bytes, 0),
        sample_id: get_u32(&bytes, 8),
        contig_id: get_u32(&bytes, 12),
        spatial: flags & 1 != 0,
        whole_sample: flags & 2 != 0,
    }))
}

#[allow(clippy::too_many_arguments)]
fn write_stream_mapping(
    records: &mut impl Write,
    payload: &mut impl Write,
    payload_bytes: &mut u64,
    stats: &mut MappingStats,
    screen_samples: &mut [PartScreenSample],
    hash: u64,
    sample_id: u32,
    posting: &PostingBuildRecord,
) -> Result<(), JamIndexPartError> {
    let sample = screen_samples
        .get_mut(usize::try_from(sample_id).map_err(|_| JamIndexPartError::Overflow)?)
        .ok_or(JamIndexPartError::Overflow)?;
    sample.hashes.push(hash);
    let count = u32::try_from(posting.contig_ids.len()).map_err(|_| JamIndexPartError::Overflow)?;
    if count == 0 || count > MAPPING_COUNT_MASK || (!posting.spatial && !posting.whole_sample) {
        return Err(JamIndexPartError::Overflow);
    }
    let mut metadata = count;
    if posting.spatial {
        metadata |= MAPPING_SPATIAL;
    }
    if posting.whole_sample {
        metadata |= MAPPING_WHOLE_SAMPLE;
    }
    let value = if count == 1 {
        stats.single = stats.single.saturating_add(1);
        *posting
            .contig_ids
            .first()
            .ok_or(JamIndexPartError::Corrupt("empty inline mapping"))?
    } else {
        stats.overflow = stats.overflow.saturating_add(1);
        stats.overflow_contigs = stats.overflow_contigs.saturating_add(u64::from(count));
        stats.maximum_overflow_count = stats.maximum_overflow_count.max(count);
        metadata |= MAPPING_OVERFLOW;
        let offset = u32::try_from(*payload_bytes).map_err(|_| JamIndexPartError::Overflow)?;
        let mut encoded = Vec::new();
        encode_deltas(&posting.contig_ids, &mut encoded)?;
        payload.write_all(&encoded)?;
        *payload_bytes =
            (*payload_bytes).saturating_add(u64::try_from(encoded.len()).unwrap_or(u64::MAX));
        offset
    };
    let mut bytes = [0u8; MAPPING_RECORD_SIZE];
    put_u32(&mut bytes, 0, value);
    put_u32(&mut bytes, 4, metadata);
    records.write_all(&bytes)?;
    Ok(())
}

pub fn write_part(
    output: impl AsRef<Path>,
    sources: &[MetagenomeSource],
    policy: &ScreenSelectionPolicy,
) -> Result<PartWriteResult, JamIndexPartError> {
    let staged = sources
        .iter()
        .map(|source| {
            Ok(StagedMetagenomeSource {
                metagenome_id: source.metagenome_id.clone(),
                staged_sequence_path: source.sequence_path.clone(),
                published_sequence_path: source.sequence_path.clone(),
                expected_source_size: std::fs::metadata(&source.sequence_path)?.len(),
                expected_source_sha256: None,
            })
        })
        .collect::<Result<Vec<_>, JamIndexPartError>>()?;
    write_part_staged(output, &staged, policy)
}

pub fn write_part_staged(
    output: impl AsRef<Path>,
    sources: &[StagedMetagenomeSource],
    policy: &ScreenSelectionPolicy,
) -> Result<PartWriteResult, JamIndexPartError> {
    let output = output.as_ref();
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
            || !source.staged_sequence_path.is_file()
            || !source.published_sequence_path.is_file()
            || !seen_ids.insert(source.metagenome_id.clone())
    }) {
        return Err(JamIndexPartError::InvalidInput(
            "metagenome IDs must be unique and sequence paths must be local files".to_string(),
        ));
    }

    let mut strings = Vec::new();
    let mut source_records = Vec::with_capacity(sources.len());
    let mut metagenomes = Vec::with_capacity(sources.len());
    let mut contig_lengths = Vec::new();
    let mut exceptional_lengths = Vec::new();
    let mut contig_signature_histogram = BTreeMap::<u32, u64>::new();
    let mut signature_runs =
        SignatureRunBuilder::new(output.parent().unwrap_or_else(|| Path::new(".")))?;
    let mut screen_samples = Vec::with_capacity(sources.len());
    let mut total_bases = 0u64;
    let mut estimated_signature_count = 0u64;

    for (metagenome_index, source) in sources.iter().enumerate() {
        let metagenome_id =
            u32::try_from(metagenome_index).map_err(|_| JamIndexPartError::Overflow)?;
        let first_contig =
            u64::try_from(contig_lengths.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let exceptional_start =
            u32::try_from(exceptional_lengths.len()).map_err(|_| JamIndexPartError::Overflow)?;
        let id_offset = append_string(&mut strings, &source.metagenome_id)?;
        let id_length =
            u32::try_from(source.metagenome_id.len()).map_err(|_| JamIndexPartError::Overflow)?;
        screen_samples.push(PartScreenSample {
            metagenome_id: source.metagenome_id.clone(),
            hashes: Vec::new(),
        });
        let staged = ExternalSource::detect(&source.staged_sequence_path);
        let published = ExternalSource::detect(&source.published_sequence_path);
        if part_access(staged.access) != part_access(published.access)
            || std::fs::metadata(&staged.path)?.len() != source.expected_source_size
            || std::fs::metadata(&published.path)?.len() != source.expected_source_size
        {
            return Err(JamIndexPartError::SourceIdentity(metagenome_id));
        }
        let fai = staged.fai_path.as_ref().map(read_fai).transpose()?;
        let (path_offset, path_length) = append_path(&mut strings, Some(&published.path))?;
        let (fai_offset, fai_length) = append_path(&mut strings, published.fai_path.as_deref())?;
        let (gzi_offset, gzi_length) = append_path(&mut strings, published.gzi_path.as_deref())?;
        let source_size = source.expected_source_size;
        source_records.push(SourceBuildRecord {
            path_offset,
            path_length,
            fai_offset,
            fai_length,
            gzi_offset,
            gzi_length,
            source_size,
            access: part_access(published.access),
            source_sha256: [0; 32],
            fai_sha256: published
                .fai_path
                .as_deref()
                .map(sha256_file)
                .transpose()?
                .unwrap_or([0; 32]),
            gzi_sha256: published
                .gzi_path
                .as_deref()
                .map(sha256_file)
                .transpose()?
                .unwrap_or([0; 32]),
        });
        let mut selector = MetagenomeSignatureBuilder::new_streaming(policy.clone())?;
        let raw_digest = Arc::new(Mutex::new(Sha256::new()));
        let raw_bytes = Arc::new(AtomicU64::new(0));
        let raw = HashingReader {
            inner: File::open(&source.staged_sequence_path)?,
            digest: Arc::clone(&raw_digest),
            bytes: Arc::clone(&raw_bytes),
        };
        let mut reader = parse_fastx_reader(raw).map_err(|error| JamIndexPartError::Parse {
            path: source.staged_sequence_path.clone(),
            message: error.to_string(),
        })?;
        let mut metagenome_bases = 0u64;
        let mut metagenome_contigs = 0u32;
        while let Some(record) = reader.next() {
            let record = record.map_err(|error| JamIndexPartError::Parse {
                path: source.staged_sequence_path.clone(),
                message: error.to_string(),
            })?;
            let sequence = record.normalize(true);
            let source_ordinal = metagenome_contigs;
            let name = std::str::from_utf8(record.id()).map_err(|_| {
                JamIndexPartError::InvalidInput(format!(
                    "contig name in {} is not UTF-8",
                    source.staged_sequence_path.display()
                ))
            })?;
            let length = u64::try_from(sequence.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let indexed = fai.as_ref().and_then(|records| {
                records.get(usize::try_from(source_ordinal).unwrap_or(usize::MAX))
            });
            if indexed.is_some_and(|indexed| indexed.name != name || indexed.length != length) {
                return Err(JamIndexPartError::InvalidInput(format!(
                    "FAI does not match {} contig {}",
                    source.staged_sequence_path.display(),
                    name
                )));
            }
            let selected = selector.add_contig(&sequence)?;
            let screen_signature_count =
                u32::try_from(selected.hashes.len()).map_err(|_| JamIndexPartError::Overflow)?;
            *contig_signature_histogram
                .entry(screen_signature_count)
                .or_default() += 1;
            estimated_signature_count =
                estimated_signature_count.saturating_add(u64::from(screen_signature_count));
            for hash in selected.hashes {
                signature_runs.push(SignatureRunRecord {
                    hash,
                    sample_id: metagenome_id,
                    contig_id: source_ordinal,
                    spatial: true,
                    whole_sample: false,
                })?;
            }
            metagenome_bases = metagenome_bases.saturating_add(length);
            total_bases = total_bases.saturating_add(length);
            metagenome_contigs = metagenome_contigs
                .checked_add(1)
                .ok_or(JamIndexPartError::Overflow)?;
            push_contig_length(&mut contig_lengths, &mut exceptional_lengths, length)?;
        }
        drop(reader);
        if raw_bytes.load(Ordering::Relaxed) != source_size {
            return Err(JamIndexPartError::SourceIdentity(metagenome_id));
        }
        let source_sha256: [u8; 32] = raw_digest
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone()
            .finalize()
            .into();
        if source
            .expected_source_sha256
            .is_some_and(|expected| expected != source_sha256)
        {
            return Err(JamIndexPartError::SourceIdentity(metagenome_id));
        }
        source_records
            .last_mut()
            .ok_or(JamIndexPartError::Corrupt("source build record"))?
            .source_sha256 = source_sha256;
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
                source.staged_sequence_path.display()
            )));
        }
        estimated_signature_count =
            estimated_signature_count.saturating_add(u64::from(policy.whole_metagenome_budget));
        let selected = selector.finish();
        for (hash, source_ordinal) in &selected.whole_hash_contigs {
            if *source_ordinal >= metagenome_contigs {
                return Err(JamIndexPartError::Overflow);
            }
            signature_runs.push(SignatureRunRecord {
                hash: *hash,
                sample_id: metagenome_id,
                contig_id: *source_ordinal,
                spatial: false,
                whole_sample: true,
            })?;
        }
        metagenomes.push(MetagenomeBuildRecord {
            first_contig,
            contig_count: metagenome_contigs,
            screen_hash_count: 0,
            id_offset,
            id_length,
            exceptional_start,
            exceptional_count: u32::try_from(exceptional_lengths.len())
                .map_err(|_| JamIndexPartError::Overflow)?
                .checked_sub(exceptional_start)
                .ok_or(JamIndexPartError::Overflow)?,
            total_bases: metagenome_bases,
        });
    }
    let mappings = signature_runs.finish(&mut screen_samples)?;
    for (metagenome, sample) in metagenomes.iter_mut().zip(&screen_samples) {
        metagenome.screen_hash_count =
            u32::try_from(sample.hashes.len()).map_err(|_| JamIndexPartError::Overflow)?;
    }
    write_encoded_part_streaming(
        output,
        source_records,
        metagenomes,
        contig_lengths,
        exceptional_lengths,
        contig_signature_histogram,
        mappings,
        strings,
        screen_samples,
        total_bases,
        estimated_signature_count,
    )
}

pub fn merge_part_fragments(
    output: impl AsRef<Path>,
    fragments: &[(PathBuf, PathBuf)],
    published_sources: &BTreeMap<String, PublishedMetagenomeSource>,
    contig_signature_histogram: BTreeMap<u32, u64>,
) -> Result<PartWriteResult, JamIndexPartError> {
    if fragments.is_empty() {
        return Err(JamIndexPartError::InvalidInput(
            "part merge requires at least one fragment".to_string(),
        ));
    }
    let mut strings = Vec::new();
    let mut source_records = Vec::new();
    let mut metagenomes = Vec::new();
    let mut contig_lengths = Vec::new();
    let mut exceptional_lengths = Vec::new();
    let mut screen_samples = Vec::<PartScreenSample>::new();
    let mut seen_metagenomes = BTreeSet::new();
    let mut total_bases = 0u64;
    let mut estimated_signature_count = 0u64;
    let mut opened = Vec::with_capacity(fragments.len());
    let mut signature_runs =
        SignatureRunBuilder::new(output.as_ref().parent().unwrap_or_else(|| Path::new(".")))?;

    for (screen_path, data_path) in fragments {
        let screen = crate::reader::JamReader::open(screen_path)?;
        let data = JamIndexPartReader::open(data_path)?;
        if screen.sample_names().len() != data.metagenomes().len() {
            return Err(JamIndexPartError::Corrupt(
                "fragment screen/data sample count",
            ));
        }
        let sample_offset =
            u32::try_from(metagenomes.len()).map_err(|_| JamIndexPartError::Overflow)?;
        for (local_sample, metagenome) in data.metagenomes().iter().enumerate() {
            if screen
                .sample_name(u32::try_from(local_sample).map_err(|_| JamIndexPartError::Overflow)?)
                != Some(metagenome.metagenome_id.as_str())
                || !seen_metagenomes.insert(metagenome.metagenome_id.clone())
            {
                return Err(JamIndexPartError::Corrupt(
                    "fragment screen/data sample binding",
                ));
            }
            let published = published_sources
                .get(&metagenome.metagenome_id)
                .ok_or_else(|| {
                    JamIndexPartError::InvalidInput(format!(
                        "published source is missing for {}",
                        metagenome.metagenome_id
                    ))
                })?;
            if published.metagenome_id != metagenome.metagenome_id {
                return Err(JamIndexPartError::Corrupt("published source binding"));
            }
            let id_offset = append_string(&mut strings, &metagenome.metagenome_id)?;
            let id_length = u32::try_from(metagenome.metagenome_id.len())
                .map_err(|_| JamIndexPartError::Overflow)?;
            let external = ExternalSource::detect(&published.sequence_path);
            if std::fs::metadata(&external.path)?.len() != published.source_size {
                return Err(JamIndexPartError::SourceIdentity(
                    u32::try_from(local_sample).unwrap_or(u32::MAX),
                ));
            }
            let (path_offset, path_length) = append_path(&mut strings, Some(&external.path))?;
            let (fai_offset, fai_length) = append_path(&mut strings, external.fai_path.as_deref())?;
            let (gzi_offset, gzi_length) = append_path(&mut strings, external.gzi_path.as_deref())?;
            source_records.push(SourceBuildRecord {
                path_offset,
                path_length,
                fai_offset,
                fai_length,
                gzi_offset,
                gzi_length,
                source_size: published.source_size,
                access: part_access(external.access),
                source_sha256: published.source_sha256,
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
            let first_contig =
                u64::try_from(contig_lengths.len()).map_err(|_| JamIndexPartError::Overflow)?;
            let exceptional_start = u32::try_from(exceptional_lengths.len())
                .map_err(|_| JamIndexPartError::Overflow)?;
            for contig_id in 0..metagenome.contig_count {
                push_contig_length(
                    &mut contig_lengths,
                    &mut exceptional_lengths,
                    data.contig_length(
                        u32::try_from(local_sample).map_err(|_| JamIndexPartError::Overflow)?,
                        contig_id,
                    )?,
                )?;
            }
            metagenomes.push(MetagenomeBuildRecord {
                first_contig,
                contig_count: metagenome.contig_count,
                screen_hash_count: metagenome.screen_hash_count,
                id_offset,
                id_length,
                exceptional_start,
                exceptional_count: u32::try_from(exceptional_lengths.len())
                    .map_err(|_| JamIndexPartError::Overflow)?
                    .checked_sub(exceptional_start)
                    .ok_or(JamIndexPartError::Overflow)?,
                total_bases: metagenome.total_bases,
            });
            screen_samples.push(PartScreenSample {
                metagenome_id: metagenome.metagenome_id.clone(),
                hashes: Vec::new(),
            });
            total_bases = total_bases.saturating_add(metagenome.total_bases);
        }
        estimated_signature_count = estimated_signature_count.saturating_add(data.posting_count());
        opened.push((screen, data, sample_offset));
    }

    for bucket in 0..BUCKET_COUNT {
        for (screen, data, sample_offset) in &opened {
            let ordinal_start = screen.bucket_entry_ordinal_start(bucket);
            for (entry_index, entry) in screen.bucket_entries(bucket).iter().enumerate() {
                let ordinal = ordinal_start
                    .checked_add(
                        u64::try_from(entry_index).map_err(|_| JamIndexPartError::Overflow)?,
                    )
                    .ok_or(JamIndexPartError::Overflow)?;
                let kind = data.posting_kind(ordinal)?;
                let sample_id = sample_offset
                    .checked_add(entry.sample_id)
                    .ok_or(JamIndexPartError::Overflow)?;
                for contig_id in data.posting(ordinal)? {
                    signature_runs.push(SignatureRunRecord {
                        hash: entry.hash,
                        sample_id,
                        contig_id,
                        spatial: kind.spatial,
                        whole_sample: kind.whole_sample,
                    })?;
                }
            }
        }
    }
    if published_sources.len() != screen_samples.len() {
        return Err(JamIndexPartError::InvalidInput(
            "published source set does not match merged fragments".to_string(),
        ));
    }
    let mappings = signature_runs.finish(&mut screen_samples)?;
    write_encoded_part_streaming(
        output.as_ref(),
        source_records,
        metagenomes,
        contig_lengths,
        exceptional_lengths,
        contig_signature_histogram,
        mappings,
        strings,
        screen_samples,
        total_bases,
        estimated_signature_count,
    )
}

#[allow(clippy::too_many_arguments)]
fn write_encoded_part_streaming(
    output: &Path,
    source_records: Vec<SourceBuildRecord>,
    metagenomes: Vec<MetagenomeBuildRecord>,
    contig_lengths: Vec<u32>,
    exceptional_lengths: Vec<ExceptionalLengthRecord>,
    contig_signature_histogram: BTreeMap<u32, u64>,
    mappings: StreamingMappings,
    strings: Vec<u8>,
    screen_samples: Vec<PartScreenSample>,
    total_bases: u64,
    estimated_signature_count: u64,
) -> Result<PartWriteResult, JamIndexPartError> {
    if mappings.signature_count > estimated_signature_count {
        return Err(JamIndexPartError::Corrupt("signature estimate underflow"));
    }
    let source_bytes = encode_sources(&source_records);
    let metagenome_bytes = encode_metagenomes(&metagenomes);
    let contig_bytes = encode_contig_lengths(&contig_lengths, &exceptional_lengths);
    let mapping_records_bytes = std::fs::metadata(&mappings.records_path)?.len();
    let mapping_payload_bytes = std::fs::metadata(&mappings.payload_path)?.len();
    let signature_length = mapping_records_bytes.saturating_add(mapping_payload_bytes);
    if mapping_records_bytes
        != mappings
            .signature_count
            .saturating_mul(MAPPING_RECORD_SIZE as u64)
    {
        return Err(JamIndexPartError::Corrupt("mapping record count"));
    }
    let mut file = OpenOptions::new()
        .create_new(true)
        .read(true)
        .write(true)
        .open(output)?;
    file.write_all(&vec![0u8; HEADER_SIZE])?;
    let source_offset = file.stream_position()?;
    file.write_all(&source_bytes)?;
    let posting_offset = align_file_plain(&mut file, 8)?;
    std::io::copy(&mut File::open(&mappings.records_path)?, &mut file)?;
    std::io::copy(&mut File::open(&mappings.payload_path)?, &mut file)?;
    let metagenome_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&metagenome_bytes)?;
    let contig_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&contig_bytes)?;
    let string_offset = align_file_plain(&mut file, 8)?;
    file.write_all(&strings)?;
    let object_size = file.stream_position()?;
    file.flush()?;
    let header = Header {
        format_version: VERSION,
        metagenome_count: u32::try_from(metagenomes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_count: u64::try_from(contig_lengths.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        total_bases,
        signature_count: mappings.signature_count,
        sequence_offset: source_offset,
        sequence_length: u64::try_from(source_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        signature_offset: posting_offset,
        signature_length,
        metagenome_offset,
        metagenome_length: u64::try_from(metagenome_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        contig_offset,
        contig_length: u64::try_from(contig_bytes.len())
            .map_err(|_| JamIndexPartError::Overflow)?,
        string_offset,
        string_length: u64::try_from(strings.len()).map_err(|_| JamIndexPartError::Overflow)?,
        sequence_checksum: sha256(&source_bytes),
        signature_checksum: sha256_file_range(output, posting_offset, signature_length)?,
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
        estimated_signature_count: header.signature_count,
        posting_count: header.signature_count,
        contig_posting_bytes: header.signature_length,
        source_reference_bytes: header.sequence_length,
        metagenome_directory_bytes: header.metagenome_length,
        contig_length_bytes: header
            .contig_count
            .saturating_mul(CONTIG_LENGTH_SIZE as u64),
        exceptional_length_bytes: u64::try_from(exceptional_lengths.len())
            .unwrap_or(u64::MAX)
            .saturating_mul(EXCEPTIONAL_LENGTH_RECORD_SIZE as u64),
        string_table_bytes: header.string_length,
        data_file_bytes: object_size,
        contig_signature_histogram,
        single_contig_mappings: mappings.stats.single,
        overflow_mappings: mappings.stats.overflow,
        overflow_contigs: mappings.stats.overflow_contigs,
        maximum_overflow_count: mappings.stats.maximum_overflow_count,
        signature_run_count: mappings.run_count,
        signature_run_record_limit: SIGNATURE_RUN_RECORD_LIMIT as u64,
    })
}

pub struct JamIndexPartReader {
    mmap: Mmap,
    header: Header,
    metagenomes: Vec<PartMetagenome>,
    exceptional_lengths: BTreeMap<u64, u64>,
}

#[derive(Clone, Debug)]
pub struct PartReadResult {
    pub contigs: BTreeMap<u32, LoadedPartContig>,
    pub source_bytes: u64,
}

impl JamIndexPartReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, JamIndexPartError> {
        crate::profiling::add_counter("part_open_count", 1);
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        crate::profiling::add_counter(
            "part_mapped_bytes",
            u64::try_from(mmap.len()).unwrap_or(u64::MAX),
        );
        crate::profiling::add_counter("part_bytes_remapped", 0);
        if mmap.len() < HEADER_SIZE {
            return Err(JamIndexPartError::Corrupt("truncated header"));
        }
        let header = decode_header(&mmap[..HEADER_SIZE])?;
        let source = section(&mmap, header.sequence_offset, header.sequence_length)?;
        let posting = section(&mmap, header.signature_offset, header.signature_length)?;
        let metagenome = section(&mmap, header.metagenome_offset, header.metagenome_length)?;
        let contig = section(&mmap, header.contig_offset, header.contig_length)?;
        let strings = section(&mmap, header.string_offset, header.string_length)?;
        if sha256(source) != header.sequence_checksum
            || sha256(posting) != header.signature_checksum
            || sha256(metagenome) != header.metagenome_checksum
            || sha256(contig) != header.contig_checksum
            || sha256(strings) != header.string_checksum
        {
            return Err(JamIndexPartError::Corrupt("section checksum mismatch"));
        }
        validate_external_lengths(header)?;
        let sources = decode_sources(source, strings, header.metagenome_count)?;
        let mut metagenomes = decode_metagenomes(metagenome, strings, header.metagenome_count)?;
        for (metagenome, source) in metagenomes.iter_mut().zip(sources) {
            metagenome.source_path = source.path;
            metagenome.source_size = source.source_size;
            metagenome.source_sha256 = source.source_sha256;
            metagenome.access = source.access;
            metagenome.fai_path = source.fai_path;
            metagenome.fai_sha256 = source.fai_sha256;
            metagenome.gzi_path = source.gzi_path;
            metagenome.gzi_sha256 = source.gzi_sha256;
        }
        let exceptional_lengths = decode_contig_lengths(contig, header.contig_count)?;
        validate_directories(
            &metagenomes,
            contig,
            &exceptional_lengths,
            header.contig_count,
            header.total_bases,
        )?;
        validate_postings(posting, header.signature_count, &metagenomes)?;
        crate::profiling::add_counter("part_validation_count", 1);
        Ok(Self {
            mmap,
            header,
            metagenomes,
            exceptional_lengths,
        })
    }

    pub fn remap_source(
        &mut self,
        metagenome_id: &str,
        source_path: impl AsRef<Path>,
        source_sha256: [u8; 32],
    ) -> Result<(), JamIndexPartError> {
        let metagenome_index = self
            .metagenomes
            .iter()
            .position(|metagenome| metagenome.metagenome_id == metagenome_id)
            .ok_or_else(|| {
                JamIndexPartError::InvalidInput(format!(
                    "source override references unknown metagenome {metagenome_id}"
                ))
            })?;
        let external = ExternalSource::detect(source_path.as_ref());
        let metagenome = &self.metagenomes[metagenome_index];
        if source_sha256 != metagenome.source_sha256
            || std::fs::metadata(&external.path)?.len() != metagenome.source_size
            || part_access(external.access) != metagenome.access
        {
            return Err(JamIndexPartError::SourceIdentity(
                u32::try_from(metagenome_index).unwrap_or(u32::MAX),
            ));
        }
        let local_id = u32::try_from(metagenome_index).map_err(|_| JamIndexPartError::Overflow)?;
        check_sidecar(
            external.fai_path.as_deref(),
            metagenome.fai_sha256,
            local_id,
        )?;
        check_sidecar(
            external.gzi_path.as_deref(),
            metagenome.gzi_sha256,
            local_id,
        )?;
        let metagenome = &mut self.metagenomes[metagenome_index];
        metagenome.source_path = external.path;
        metagenome.fai_path = external.fai_path;
        metagenome.gzi_path = external.gzi_path;
        Ok(())
    }

    pub fn metagenomes(&self) -> &[PartMetagenome] {
        &self.metagenomes
    }

    pub const fn total_bases(&self) -> u64 {
        self.header.total_bases
    }

    pub const fn posting_count(&self) -> u64 {
        self.header.signature_count
    }

    pub fn object_size(&self) -> u64 {
        u64::try_from(self.mmap.len()).unwrap_or(u64::MAX)
    }

    pub fn posting(&self, ordinal: u64) -> Result<Vec<u32>, JamIndexPartError> {
        crate::profiling::add_counter("posting_records_read", 1);
        let bytes = section(
            &self.mmap,
            self.header.signature_offset,
            self.header.signature_length,
        )?;
        decode_mapping(bytes, self.header.signature_count, ordinal)
            .map(|decoded| decoded.contig_ids)
    }

    pub fn posting_kind(&self, ordinal: u64) -> Result<PostingKind, JamIndexPartError> {
        if ordinal >= self.header.signature_count {
            return Err(JamIndexPartError::UnknownPosting(ordinal));
        }
        let bytes = section(
            &self.mmap,
            self.header.signature_offset,
            self.header.signature_length,
        )?;
        let record = usize::try_from(ordinal)
            .map_err(|_| JamIndexPartError::Overflow)?
            .checked_mul(MAPPING_RECORD_SIZE)
            .ok_or(JamIndexPartError::Overflow)?;
        let metadata = get_u32(bytes, record + 4);
        Ok(PostingKind {
            spatial: metadata & MAPPING_SPATIAL != 0,
            whole_sample: metadata & MAPPING_WHOLE_SAMPLE != 0,
        })
    }

    pub fn metagenome_contigs(
        &self,
        metagenome_id: u32,
    ) -> Result<std::ops::Range<u32>, JamIndexPartError> {
        let metagenome = self
            .metagenomes
            .get(usize::try_from(metagenome_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::UnknownMetagenome(metagenome_id))?;
        Ok(0..metagenome.contig_count)
    }

    pub fn contig_length(
        &self,
        metagenome_id: u32,
        contig_id: u32,
    ) -> Result<u64, JamIndexPartError> {
        let metagenome = self
            .metagenomes
            .get(usize::try_from(metagenome_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::UnknownMetagenome(metagenome_id))?;
        if contig_id >= metagenome.contig_count {
            return Err(JamIndexPartError::UnknownContig(contig_id));
        }
        let ordinal = metagenome
            .first_contig
            .checked_add(u64::from(contig_id))
            .ok_or(JamIndexPartError::Overflow)?;
        decode_contig_length(
            section(
                &self.mmap,
                self.header.contig_offset,
                self.header.contig_length,
            )?,
            self.header.contig_count,
            &self.exceptional_lengths,
            ordinal,
        )
    }

    pub fn read_contigs(
        &self,
        metagenome_id: u32,
        contig_ids: &[u32],
    ) -> Result<PartReadResult, JamIndexPartError> {
        let (_, external) = self.checked_external_source(metagenome_id)?;
        let allowed = self.metagenome_contigs(metagenome_id)?;
        let indexed = external.fai_path.as_ref().map(read_fai).transpose()?;
        if indexed
            .as_ref()
            .is_some_and(|records| records.len() != allowed.len())
        {
            return Err(JamIndexPartError::SourceIdentity(metagenome_id));
        }
        let requests = contig_ids
            .iter()
            .copied()
            .map(|contig_id| {
                if !allowed.contains(&contig_id) {
                    return Err(JamIndexPartError::UnknownContig(contig_id));
                }
                let length = self.contig_length(metagenome_id, contig_id)?;
                let record = indexed.as_ref().and_then(|records| {
                    records.get(usize::try_from(contig_id).unwrap_or(usize::MAX))
                });
                if record.is_some_and(|record| record.length != length) {
                    return Err(JamIndexPartError::ContigChecksum(contig_id));
                }
                Ok(ContigRequest {
                    contig_id,
                    source_ordinal: contig_id,
                    name: record.map(|record| record.name.clone()),
                    length,
                    offset: record.map_or(0, |record| record.offset),
                    line_bases: record.map_or(0, |record| record.line_bases),
                    line_width: record.map_or(0, |record| record.line_width),
                })
            })
            .collect::<Result<Vec<_>, JamIndexPartError>>()?;
        let loaded = read_external(&external, &requests)?;
        let mut contigs = BTreeMap::new();
        for loaded in loaded.contigs {
            contigs.insert(
                loaded.contig_id,
                LoadedPartContig {
                    name: loaded.name,
                    bases: Arc::from(loaded.bases),
                },
            );
        }
        Ok(PartReadResult {
            contigs,
            source_bytes: loaded.source_bytes,
        })
    }

    pub fn stream_contigs<E, F>(
        &self,
        metagenome_id: u32,
        batch_size: usize,
        mut emit: F,
    ) -> Result<u64, E>
    where
        E: From<JamIndexPartError>,
        F: FnMut(BTreeMap<u32, LoadedPartContig>) -> Result<(), E>,
    {
        if batch_size == 0 {
            return Err(E::from(JamIndexPartError::InvalidInput(
                "stream batch size must be positive".to_string(),
            )));
        }
        let (metagenome, external) = self
            .checked_external_source(metagenome_id)
            .map_err(E::from)?;
        let allowed = self.metagenome_contigs(metagenome_id).map_err(E::from)?;
        let mut reader = parse_fastx_file(&external.path).map_err(|error| {
            E::from(JamIndexPartError::Parse {
                path: external.path.clone(),
                message: error.to_string(),
            })
        })?;
        let mut ordinal = 0u32;
        let mut batch = BTreeMap::new();
        while let Some(record) = reader.next() {
            let record = record.map_err(|error| {
                E::from(JamIndexPartError::Parse {
                    path: external.path.clone(),
                    message: error.to_string(),
                })
            })?;
            let contig_id = ordinal;
            if contig_id >= allowed.end {
                return Err(E::from(JamIndexPartError::SourceIdentity(metagenome_id)));
            }
            let name = std::str::from_utf8(record.id()).map_err(|_| {
                E::from(JamIndexPartError::InvalidInput(format!(
                    "contig name in {} is not UTF-8",
                    external.path.display()
                )))
            })?;
            let sequence = record.normalize(true).into_owned();
            if self
                .contig_length(metagenome_id, contig_id)
                .map_err(E::from)?
                != u64::try_from(sequence.len()).unwrap_or(u64::MAX)
            {
                return Err(E::from(JamIndexPartError::ContigChecksum(contig_id)));
            }
            batch.insert(
                contig_id,
                LoadedPartContig {
                    name: name.to_string(),
                    bases: Arc::from(sequence),
                },
            );
            ordinal = ordinal
                .checked_add(1)
                .ok_or_else(|| E::from(JamIndexPartError::Overflow))?;
            if batch.len() == batch_size {
                emit(std::mem::take(&mut batch))?;
            }
        }
        if ordinal != metagenome.contig_count {
            return Err(E::from(JamIndexPartError::SourceIdentity(metagenome_id)));
        }
        if !batch.is_empty() {
            emit(batch)?;
        }
        Ok(metagenome.source_size)
    }

    fn checked_external_source(
        &self,
        metagenome_id: u32,
    ) -> Result<(&PartMetagenome, ExternalSource), JamIndexPartError> {
        let metagenome = self
            .metagenomes
            .get(usize::try_from(metagenome_id).map_err(|_| JamIndexPartError::Overflow)?)
            .ok_or(JamIndexPartError::UnknownMetagenome(metagenome_id))?;
        let external = ExternalSource::detect(&metagenome.source_path);
        if std::fs::metadata(&external.path)?.len() != metagenome.source_size
            || part_access(external.access) != metagenome.access
            || external.fai_path != metagenome.fai_path
            || external.gzi_path != metagenome.gzi_path
        {
            return Err(JamIndexPartError::SourceIdentity(metagenome_id));
        }
        check_sidecar(
            metagenome.fai_path.as_deref(),
            metagenome.fai_sha256,
            metagenome_id,
        )?;
        check_sidecar(
            metagenome.gzi_path.as_deref(),
            metagenome.gzi_sha256,
            metagenome_id,
        )?;
        Ok((metagenome, external))
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

fn push_contig_length(
    lengths: &mut Vec<u32>,
    exceptional: &mut Vec<ExceptionalLengthRecord>,
    length: u64,
) -> Result<(), JamIndexPartError> {
    let ordinal = u64::try_from(lengths.len()).map_err(|_| JamIndexPartError::Overflow)?;
    if length < u64::from(u32::MAX) {
        lengths.push(u32::try_from(length).map_err(|_| JamIndexPartError::Overflow)?);
    } else {
        lengths.push(u32::MAX);
        exceptional.push(ExceptionalLengthRecord { ordinal, length });
    }
    Ok(())
}

fn encode_contig_lengths(lengths: &[u32], exceptional: &[ExceptionalLengthRecord]) -> Vec<u8> {
    let length_bytes = lengths.len().saturating_mul(CONTIG_LENGTH_SIZE);
    let mut bytes = vec![
        0u8;
        length_bytes.saturating_add(
            exceptional
                .len()
                .saturating_mul(EXCEPTIONAL_LENGTH_RECORD_SIZE),
        )
    ];
    for (index, length) in lengths.iter().copied().enumerate() {
        put_u32(&mut bytes, index * CONTIG_LENGTH_SIZE, length);
    }
    for (index, record) in exceptional.iter().enumerate() {
        let offset = length_bytes + index * EXCEPTIONAL_LENGTH_RECORD_SIZE;
        put_u64(&mut bytes, offset, record.ordinal);
        put_u64(&mut bytes, offset + 8, record.length);
    }
    bytes
}

#[derive(Clone, Copy, Debug, Default)]
struct MappingStats {
    single: u64,
    overflow: u64,
    overflow_contigs: u64,
    maximum_overflow_count: u32,
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

struct DecodedSource {
    path: PathBuf,
    source_size: u64,
    access: PartAccess,
    source_sha256: [u8; 32],
    fai_path: Option<PathBuf>,
    fai_sha256: Option<[u8; 32]>,
    gzi_path: Option<PathBuf>,
    gzi_sha256: Option<[u8; 32]>,
}

fn decode_sources(
    bytes: &[u8],
    strings: &[u8],
    count: u32,
) -> Result<Vec<DecodedSource>, JamIndexPartError> {
    let count = usize::try_from(count).map_err(|_| JamIndexPartError::Overflow)?;
    if bytes.len() != count * SOURCE_RECORD_SIZE {
        return Err(JamIndexPartError::Corrupt("source directory length"));
    }
    (0..count)
        .map(|index| {
            let offset = index * SOURCE_RECORD_SIZE;
            if bytes[offset + 152..offset + SOURCE_RECORD_SIZE]
                .iter()
                .any(|byte| *byte != 0)
            {
                return Err(JamIndexPartError::Corrupt("source directory record"));
            }
            let access = match get_u32(bytes, offset + 12) {
                1 => PartAccess::PlainFai,
                2 => PartAccess::Bgzf,
                3 => PartAccess::Sequential,
                _ => return Err(JamIndexPartError::Corrupt("source access mode")),
            };
            let fai_path = decode_path(
                strings,
                get_u64(bytes, offset + 16),
                get_u32(bytes, offset + 24),
            )?;
            let gzi_path = decode_path(
                strings,
                get_u64(bytes, offset + 32),
                get_u32(bytes, offset + 40),
            )?;
            let fai_checksum = array_32(bytes, offset + 88);
            let gzi_checksum = array_32(bytes, offset + 120);
            if fai_path.is_none() != (fai_checksum == [0; 32])
                || gzi_path.is_none() != (gzi_checksum == [0; 32])
                || matches!(access, PartAccess::PlainFai) && fai_path.is_none()
                || matches!(access, PartAccess::Bgzf) && (fai_path.is_none() || gzi_path.is_none())
            {
                return Err(JamIndexPartError::Corrupt("source sidecar binding"));
            }
            Ok(DecodedSource {
                path: PathBuf::from(decode_string(
                    strings,
                    get_u64(bytes, offset),
                    get_u32(bytes, offset + 8),
                )?),
                source_size: get_u64(bytes, offset + 48),
                access,
                source_sha256: array_32(bytes, offset + 56),
                fai_path,
                fai_sha256: (fai_checksum != [0; 32]).then_some(fai_checksum),
                gzi_path,
                gzi_sha256: (gzi_checksum != [0; 32]).then_some(gzi_checksum),
            })
        })
        .collect()
}

fn decode_contig_lengths(
    bytes: &[u8],
    count: u64,
) -> Result<BTreeMap<u64, u64>, JamIndexPartError> {
    let length_bytes = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(CONTIG_LENGTH_SIZE))
        .ok_or(JamIndexPartError::Overflow)?;
    let exceptional_bytes = bytes
        .len()
        .checked_sub(length_bytes)
        .ok_or(JamIndexPartError::Corrupt("contig length range"))?;
    if exceptional_bytes % EXCEPTIONAL_LENGTH_RECORD_SIZE != 0 {
        return Err(JamIndexPartError::Corrupt("exceptional length table"));
    }
    let mut exceptional = BTreeMap::new();
    for index in 0..exceptional_bytes / EXCEPTIONAL_LENGTH_RECORD_SIZE {
        let offset = length_bytes + index * EXCEPTIONAL_LENGTH_RECORD_SIZE;
        let ordinal = get_u64(bytes, offset);
        let length = get_u64(bytes, offset + 8);
        let word_offset = usize::try_from(ordinal)
            .ok()
            .and_then(|ordinal| ordinal.checked_mul(CONTIG_LENGTH_SIZE))
            .ok_or(JamIndexPartError::Overflow)?;
        if ordinal >= count
            || length < u64::from(u32::MAX)
            || get_u32(bytes, word_offset) != u32::MAX
            || exceptional.insert(ordinal, length).is_some()
        {
            return Err(JamIndexPartError::Corrupt("exceptional length record"));
        }
    }
    let sentinels = (0..usize::try_from(count).map_err(|_| JamIndexPartError::Overflow)?)
        .filter(|index| get_u32(bytes, index * CONTIG_LENGTH_SIZE) == u32::MAX)
        .count();
    if sentinels != exceptional.len() {
        return Err(JamIndexPartError::Corrupt("exceptional length binding"));
    }
    Ok(exceptional)
}

fn decode_contig_length(
    bytes: &[u8],
    count: u64,
    exceptional: &BTreeMap<u64, u64>,
    ordinal: u64,
) -> Result<u64, JamIndexPartError> {
    if ordinal >= count {
        return Err(JamIndexPartError::Overflow);
    }
    let offset = usize::try_from(ordinal)
        .ok()
        .and_then(|ordinal| ordinal.checked_mul(CONTIG_LENGTH_SIZE))
        .ok_or(JamIndexPartError::Overflow)?;
    let value = get_u32(bytes, offset);
    if value == u32::MAX {
        exceptional
            .get(&ordinal)
            .copied()
            .ok_or(JamIndexPartError::Corrupt("missing exceptional length"))
    } else {
        Ok(u64::from(value))
    }
}

fn validate_postings(
    bytes: &[u8],
    count: u64,
    metagenomes: &[PartMetagenome],
) -> Result<(), JamIndexPartError> {
    let maximum_contigs = metagenomes
        .iter()
        .map(|metagenome| metagenome.contig_count)
        .max()
        .ok_or(JamIndexPartError::Corrupt("empty metagenome directory"))?;
    let mut expected_overflow_start = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(MAPPING_RECORD_SIZE))
        .ok_or(JamIndexPartError::Overflow)?;
    for ordinal in 0..count {
        let decoded = decode_mapping(bytes, count, ordinal)?;
        let posting = decoded.contig_ids;
        let overflow_range = decoded.overflow_range.map(|range| (range.start, range.end));
        if posting.is_empty()
            || posting
                .iter()
                .any(|contig_id| *contig_id >= maximum_contigs)
        {
            return Err(JamIndexPartError::Corrupt("contig posting range"));
        }
        if let Some((start, end)) = overflow_range {
            if start != expected_overflow_start {
                return Err(JamIndexPartError::Corrupt("mapping overflow order"));
            }
            expected_overflow_start = end;
        }
    }
    if expected_overflow_start != bytes.len() {
        return Err(JamIndexPartError::Corrupt("mapping overflow length"));
    }
    Ok(())
}

fn decode_mapping(
    bytes: &[u8],
    count: u64,
    ordinal: u64,
) -> Result<DecodedMapping, JamIndexPartError> {
    if ordinal >= count {
        return Err(JamIndexPartError::UnknownPosting(ordinal));
    }
    let record_bytes = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(MAPPING_RECORD_SIZE))
        .ok_or(JamIndexPartError::Overflow)?;
    if record_bytes > bytes.len() {
        return Err(JamIndexPartError::Corrupt("mapping record range"));
    }
    let record = usize::try_from(ordinal)
        .map_err(|_| JamIndexPartError::Overflow)?
        .checked_mul(MAPPING_RECORD_SIZE)
        .ok_or(JamIndexPartError::Overflow)?;
    let value = get_u32(bytes, record);
    let metadata = get_u32(bytes, record + 4);
    let expected = metadata & MAPPING_COUNT_MASK;
    if expected == 0 {
        return Err(JamIndexPartError::Corrupt("empty mapping record"));
    }
    if metadata & MAPPING_OVERFLOW == 0 {
        if expected != 1 {
            return Err(JamIndexPartError::Corrupt("inline mapping count"));
        }
        return Ok(DecodedMapping {
            contig_ids: vec![value],
            overflow_range: None,
        });
    }
    let start = record_bytes
        .checked_add(usize::try_from(value).map_err(|_| JamIndexPartError::Overflow)?)
        .ok_or(JamIndexPartError::Overflow)?;
    let mut cursor = start;
    let mut values =
        Vec::with_capacity(usize::try_from(expected).map_err(|_| JamIndexPartError::Overflow)?);
    let mut previous = 0u32;
    while values.len() < usize::try_from(expected).unwrap_or(usize::MAX) {
        let mut delta = 0u32;
        let mut shift = 0u32;
        loop {
            let byte = *bytes
                .get(cursor)
                .ok_or(JamIndexPartError::Corrupt("truncated mapping overflow"))?;
            cursor += 1;
            delta |= u32::from(byte & 0x7f)
                .checked_shl(shift)
                .ok_or(JamIndexPartError::Overflow)?;
            if byte & 0x80 == 0 {
                break;
            }
            shift = shift.checked_add(7).ok_or(JamIndexPartError::Overflow)?;
            if shift >= 32 {
                return Err(JamIndexPartError::Corrupt("mapping varint overflow"));
            }
        }
        let contig_id = if values.is_empty() {
            delta
        } else {
            previous
                .checked_add(delta)
                .ok_or(JamIndexPartError::Overflow)?
        };
        if !values.is_empty() && contig_id <= previous {
            return Err(JamIndexPartError::Corrupt("mapping order"));
        }
        values.push(contig_id);
        previous = contig_id;
    }
    Ok(DecodedMapping {
        contig_ids: values,
        overflow_range: Some(start..cursor),
    })
}

fn validate_external_lengths(header: Header) -> Result<(), JamIndexPartError> {
    let header_bytes = header
        .signature_count
        .checked_mul(MAPPING_RECORD_SIZE as u64)
        .ok_or(JamIndexPartError::Overflow)?;
    let normal_lengths = header
        .contig_count
        .checked_mul(CONTIG_LENGTH_SIZE as u64)
        .ok_or(JamIndexPartError::Overflow)?;
    if header.sequence_length != u64::from(header.metagenome_count) * SOURCE_RECORD_SIZE as u64
        || header.metagenome_length
            != u64::from(header.metagenome_count) * METAGENOME_RECORD_SIZE as u64
        || header.contig_length < normal_lengths
        || !(header.contig_length - normal_lengths)
            .is_multiple_of(EXCEPTIONAL_LENGTH_RECORD_SIZE as u64)
        || header.signature_length < header_bytes
    {
        return Err(JamIndexPartError::Corrupt("section length mismatch"));
    }
    Ok(())
}

fn decode_path(
    strings: &[u8],
    offset: u64,
    length: u32,
) -> Result<Option<PathBuf>, JamIndexPartError> {
    if length == 0 {
        return Ok(None);
    }
    decode_string(strings, offset, length)
        .map(PathBuf::from)
        .map(Some)
}

fn check_sidecar(
    path: Option<&Path>,
    expected: Option<[u8; 32]>,
    metagenome_id: u32,
) -> Result<(), JamIndexPartError> {
    match (path, expected) {
        (Some(path), Some(expected)) if sha256_sidecar(path)? == expected => Ok(()),
        (None, None) => Ok(()),
        _ => Err(JamIndexPartError::SourceIdentity(metagenome_id)),
    }
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
    put_u16(&mut bytes, 8, header.format_version);
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
    if bytes.len() != HEADER_SIZE {
        return Err(JamIndexPartError::Corrupt("invalid part header"));
    }
    let format_version = get_u16(bytes, 8);
    if bytes[..8] != MAGIC
        || format_version != VERSION
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
        format_version,
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
        put_u64(&mut bytes, offset, record.first_contig);
        put_u32(&mut bytes, offset + 8, record.contig_count);
        put_u32(&mut bytes, offset + 12, record.screen_hash_count);
        put_u64(&mut bytes, offset + 16, record.id_offset);
        put_u32(&mut bytes, offset + 24, record.id_length);
        put_u32(&mut bytes, offset + 28, record.exceptional_start);
        put_u32(&mut bytes, offset + 32, record.exceptional_count);
        put_u64(&mut bytes, offset + 40, record.total_bases);
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
            if bytes[offset + 36..offset + 40]
                .iter()
                .any(|byte| *byte != 0)
            {
                return Err(JamIndexPartError::Corrupt("metagenome directory record"));
            }
            Ok(PartMetagenome {
                metagenome_id: decode_string(
                    strings,
                    get_u64(bytes, offset + 16),
                    get_u32(bytes, offset + 24),
                )?,
                first_contig: get_u64(bytes, offset),
                contig_count: get_u32(bytes, offset + 8),
                screen_hash_count: get_u32(bytes, offset + 12),
                exceptional_length_start: get_u32(bytes, offset + 28),
                exceptional_length_count: get_u32(bytes, offset + 32),
                total_bases: get_u64(bytes, offset + 40),
                source_path: PathBuf::new(),
                source_size: 0,
                source_sha256: [0; 32],
                access: PartAccess::Sequential,
                fai_path: None,
                fai_sha256: None,
                gzi_path: None,
                gzi_sha256: None,
            })
        })
        .collect()
}

fn validate_directories(
    metagenomes: &[PartMetagenome],
    contigs: &[u8],
    exceptional: &BTreeMap<u64, u64>,
    contig_count: u64,
    total_bases: u64,
) -> Result<(), JamIndexPartError> {
    let mut expected_contig = 0u64;
    let mut expected_exceptional = 0u32;
    let mut observed_bases = 0u64;
    for metagenome in metagenomes {
        if metagenome.first_contig != expected_contig
            || metagenome.exceptional_length_start != expected_exceptional
            || metagenome.contig_count == 0
        {
            return Err(JamIndexPartError::Corrupt("metagenome contig binding"));
        }
        let end = metagenome
            .first_contig
            .checked_add(u64::from(metagenome.contig_count))
            .ok_or(JamIndexPartError::Overflow)?;
        let mut metagenome_bases = 0u64;
        let mut exceptional_count = 0u32;
        for ordinal in metagenome.first_contig..end {
            let length = decode_contig_length(contigs, contig_count, exceptional, ordinal)?;
            metagenome_bases = metagenome_bases.saturating_add(length);
            exceptional_count =
                exceptional_count.saturating_add(u32::from(length >= u64::from(u32::MAX)));
        }
        if metagenome_bases != metagenome.total_bases
            || exceptional_count != metagenome.exceptional_length_count
        {
            return Err(JamIndexPartError::Corrupt("metagenome contig binding"));
        }
        expected_contig = end;
        expected_exceptional = expected_exceptional
            .checked_add(exceptional_count)
            .ok_or(JamIndexPartError::Overflow)?;
        observed_bases = observed_bases.saturating_add(metagenome_bases);
    }
    if expected_contig != contig_count
        || usize::try_from(expected_exceptional).ok() != Some(exceptional.len())
        || observed_bases != total_bases
    {
        return Err(JamIndexPartError::Corrupt("total base count mismatch"));
    }
    Ok(())
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

fn sha256_file_range(path: &Path, offset: u64, length: u64) -> Result<[u8; 32], JamIndexPartError> {
    let mut file = File::open(path)?;
    file.seek(SeekFrom::Start(offset))?;
    let mut reader = BufReader::with_capacity(1024 * 1024, file.take(length));
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    let mut observed = 0u64;
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        observed = observed.saturating_add(u64::try_from(read).unwrap_or(u64::MAX));
        digest.update(&buffer[..read]);
    }
    if observed != length {
        return Err(JamIndexPartError::Corrupt("truncated checksum range"));
    }
    Ok(digest.finalize().into())
}

struct HashingReader<R> {
    inner: R,
    digest: Arc<Mutex<Sha256>>,
    bytes: Arc<AtomicU64>,
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        if read != 0 {
            self.digest
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .update(&buffer[..read]);
            self.bytes
                .fetch_add(u64::try_from(read).unwrap_or(u64::MAX), Ordering::Relaxed);
        }
        Ok(read)
    }
}

fn sha256_sidecar(path: &Path) -> Result<[u8; 32], JamIndexPartError> {
    let mut reader = File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = [0u8; 16 * 1024];
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
    #[error("Jam Index fragment screen failed: {0}")]
    Screen(#[from] crate::reader::ReaderError),
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
    #[error("unknown Jam Index posting {0}")]
    UnknownPosting(u64),
    #[error("Jam Index source identity changed for metagenome {0}")]
    SourceIdentity(u32),
    #[error("Jam Index external contig {0} checksum mismatch")]
    ContigChecksum(u32),
    #[error("Jam Index part coordinate overflow")]
    Overflow,
}

impl From<ExternalError> for JamIndexPartError {
    fn from(error: ExternalError) -> Self {
        Self::External(error.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exceptional_contig_lengths_use_bounded_overflow_records() {
        let mut lengths = Vec::new();
        let mut exceptional = Vec::new();
        push_contig_length(&mut lengths, &mut exceptional, 123).unwrap();
        push_contig_length(&mut lengths, &mut exceptional, u64::from(u32::MAX) + 17).unwrap();
        assert_eq!(lengths, vec![123, u32::MAX]);
        assert_eq!(exceptional.len(), 1);
        let encoded = encode_contig_lengths(&lengths, &exceptional);
        let decoded = decode_contig_lengths(&encoded, 2).unwrap();
        assert_eq!(decode_contig_length(&encoded, 2, &decoded, 0).unwrap(), 123);
        assert_eq!(
            decode_contig_length(&encoded, 2, &decoded, 1).unwrap(),
            u64::from(u32::MAX) + 17
        );
    }
}
