//! Candidate-only reads from external FASTA resources.

use flate2::read::MultiGzDecoder;
use needletail::{Sequence, parse_fastx_file};
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SequenceAccess {
    PlainFai,
    Bgzf,
    Sequential,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct FaiRecord {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ExternalSource {
    pub path: PathBuf,
    pub access: SequenceAccess,
    pub fai_path: Option<PathBuf>,
    pub gzi_path: Option<PathBuf>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ContigRequest {
    pub contig_id: u32,
    pub source_ordinal: u32,
    pub name: Option<String>,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct LoadedSequence {
    pub contig_id: u32,
    pub name: String,
    pub bases: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ReadResult {
    pub contigs: Vec<LoadedSequence>,
    pub source_bytes: u64,
}

impl ExternalSource {
    pub(crate) fn detect(path: impl AsRef<Path>) -> Self {
        let path = path.as_ref().to_path_buf();
        let fai_path = sidecar_path(&path, ".fai");
        let gzi_path = sidecar_path(&path, ".gzi");
        let compressed = path
            .extension()
            .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"));
        let access = if fai_path.is_file() && gzi_path.is_file() {
            SequenceAccess::Bgzf
        } else if fai_path.is_file() && !compressed {
            SequenceAccess::PlainFai
        } else {
            SequenceAccess::Sequential
        };
        Self {
            path,
            access,
            fai_path: fai_path.is_file().then_some(fai_path),
            gzi_path: gzi_path.is_file().then_some(gzi_path),
        }
    }
}

pub(crate) fn read_fai(path: impl AsRef<Path>) -> Result<Vec<FaiRecord>, ExternalError> {
    let text = std::fs::read_to_string(path)?;
    let mut records = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 5 {
            return Err(ExternalError::InvalidFai(line_index + 1));
        }
        let record = FaiRecord {
            name: fields[0].to_string(),
            length: fields[1]
                .parse()
                .map_err(|_| ExternalError::InvalidFai(line_index + 1))?,
            offset: fields[2]
                .parse()
                .map_err(|_| ExternalError::InvalidFai(line_index + 1))?,
            line_bases: fields[3]
                .parse()
                .map_err(|_| ExternalError::InvalidFai(line_index + 1))?,
            line_width: fields[4]
                .parse()
                .map_err(|_| ExternalError::InvalidFai(line_index + 1))?,
        };
        if record.name.is_empty()
            || record.length == 0
            || record.line_bases == 0
            || record.line_width < record.line_bases
        {
            return Err(ExternalError::InvalidFai(line_index + 1));
        }
        records.push(record);
    }
    if records.is_empty() {
        return Err(ExternalError::InvalidFai(0));
    }
    Ok(records)
}

pub(crate) fn read_selected(
    source: &ExternalSource,
    requests: &[ContigRequest],
) -> Result<ReadResult, ExternalError> {
    if requests.is_empty() {
        return Ok(ReadResult {
            contigs: Vec::new(),
            source_bytes: 0,
        });
    }
    match source.access {
        SequenceAccess::PlainFai => read_plain(source, requests),
        SequenceAccess::Bgzf => read_bgzf(source, requests),
        SequenceAccess::Sequential => read_stream(source, requests),
    }
}

fn read_plain(
    source: &ExternalSource,
    requests: &[ContigRequest],
) -> Result<ReadResult, ExternalError> {
    let mut file = File::open(&source.path)?;
    let mut contigs = Vec::with_capacity(requests.len());
    let mut source_bytes = 0u64;
    for request in requests {
        let span = raw_span(request)?;
        file.seek(SeekFrom::Start(request.offset))?;
        let mut raw = vec![0u8; usize::try_from(span).map_err(|_| ExternalError::Overflow)?];
        file.read_exact(&mut raw)?;
        source_bytes = source_bytes.saturating_add(span);
        contigs.push(LoadedSequence {
            contig_id: request.contig_id,
            name: request
                .name
                .clone()
                .ok_or(ExternalError::InvalidCoordinates)?,
            bases: strip_lines(raw, request.length)?,
        });
    }
    Ok(ReadResult {
        contigs,
        source_bytes,
    })
}

fn read_bgzf(
    source: &ExternalSource,
    requests: &[ContigRequest],
) -> Result<ReadResult, ExternalError> {
    let gzi_path = source.gzi_path.as_ref().ok_or(ExternalError::MissingGzi)?;
    let entries = read_gzi(gzi_path)?;
    let mut contigs = Vec::with_capacity(requests.len());
    let mut source_bytes = 0u64;
    for request in requests {
        let &(compressed, uncompressed) = entries
            .iter()
            .rev()
            .find(|(_, uncompressed)| *uncompressed <= request.offset)
            .ok_or(ExternalError::InvalidGzi)?;
        let mut file = {
            let _profile = crate::profiling::scope("gzip_open");
            File::open(&source.path)?
        };
        file.seek(SeekFrom::Start(compressed))?;
        let counter = Arc::new(AtomicU64::new(0));
        let counted = CountReader {
            inner: file,
            count: Arc::clone(&counter),
        };
        let mut decoder = MultiGzDecoder::new(counted);
        let (span, raw) = {
            let _profile = crate::profiling::scope("gzip_decompression");
            io::copy(
                &mut decoder
                    .by_ref()
                    .take(request.offset.saturating_sub(uncompressed)),
                &mut io::sink(),
            )?;
            let span = raw_span(request)?;
            let mut raw = Vec::with_capacity(usize::try_from(span).unwrap_or(0));
            decoder.by_ref().take(span).read_to_end(&mut raw)?;
            (span, raw)
        };
        if u64::try_from(raw.len()).unwrap_or(u64::MAX) != span {
            return Err(ExternalError::TruncatedSequence);
        }
        source_bytes = source_bytes.saturating_add(counter.load(Ordering::Relaxed));
        let bases = {
            let _profile = crate::profiling::scope("selected_contig_extraction");
            strip_lines(raw, request.length)?
        };
        contigs.push(LoadedSequence {
            contig_id: request.contig_id,
            name: request
                .name
                .clone()
                .ok_or(ExternalError::InvalidCoordinates)?,
            bases,
        });
    }
    Ok(ReadResult {
        contigs,
        source_bytes,
    })
}

fn read_stream(
    source: &ExternalSource,
    requests: &[ContigRequest],
) -> Result<ReadResult, ExternalError> {
    let requested = requests
        .iter()
        .map(|request| (request.source_ordinal, request))
        .collect::<BTreeMap<_, _>>();
    let mut found = BTreeMap::new();
    let mut reader = {
        let _profile = crate::profiling::scope("gzip_open");
        parse_fastx_file(&source.path).map_err(|error| ExternalError::Parse {
            path: source.path.clone(),
            message: error.to_string(),
        })?
    };
    let mut ordinal = 0u32;
    let _decompression_profile = crate::profiling::scope("gzip_decompression");
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| ExternalError::Parse {
            path: source.path.clone(),
            message: error.to_string(),
        })?;
        if let Some(request) = requested.get(&ordinal) {
            let _profile = crate::profiling::scope("selected_contig_extraction");
            let name = std::str::from_utf8(record.id()).map_err(|_| ExternalError::InvalidName)?;
            let bases = record.normalize(true).into_owned();
            if request
                .name
                .as_deref()
                .is_some_and(|expected| name != expected)
                || u64::try_from(bases.len()).unwrap_or(u64::MAX) != request.length
            {
                return Err(ExternalError::ContigMismatch(request.contig_id));
            }
            found.insert(
                ordinal,
                LoadedSequence {
                    contig_id: request.contig_id,
                    name: name.to_string(),
                    bases,
                },
            );
        }
        ordinal = ordinal.checked_add(1).ok_or(ExternalError::Overflow)?;
    }
    drop(_decompression_profile);
    if found.len() != requested.len() {
        return Err(ExternalError::MissingContig);
    }
    Ok(ReadResult {
        contigs: requested
            .keys()
            .map(|ordinal| found.remove(ordinal).expect("selected contig was checked"))
            .collect(),
        source_bytes: std::fs::metadata(&source.path)?.len(),
    })
}

fn read_gzi(path: impl AsRef<Path>) -> Result<Vec<(u64, u64)>, ExternalError> {
    let bytes = std::fs::read(path)?;
    if bytes.len() < 8 {
        return Err(ExternalError::InvalidGzi);
    }
    let count = u64::from_le_bytes(bytes[..8].try_into().expect("eight-byte prefix"));
    let expected = usize::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(16))
        .and_then(|body| body.checked_add(8))
        .ok_or(ExternalError::Overflow)?;
    if bytes.len() != expected {
        return Err(ExternalError::InvalidGzi);
    }
    let mut entries = vec![(0, 0)];
    let (chunks, remainder) = bytes[8..].as_chunks::<16>();
    debug_assert!(remainder.is_empty());
    for chunk in chunks {
        let compressed = u64::from_le_bytes(chunk[..8].try_into().expect("eight-byte field"));
        let uncompressed = u64::from_le_bytes(chunk[8..].try_into().expect("eight-byte field"));
        if entries
            .last()
            .is_some_and(|previous| compressed <= previous.0 || uncompressed <= previous.1)
        {
            return Err(ExternalError::InvalidGzi);
        }
        entries.push((compressed, uncompressed));
    }
    Ok(entries)
}

fn raw_span(request: &ContigRequest) -> Result<u64, ExternalError> {
    if request.length == 0 || request.line_bases == 0 || request.line_width < request.line_bases {
        return Err(ExternalError::InvalidCoordinates);
    }
    let line_bases = u64::from(request.line_bases);
    let line_width = u64::from(request.line_width);
    let complete = request.length.saturating_sub(1) / line_bases;
    complete
        .checked_mul(line_width)
        .and_then(|bytes| bytes.checked_add((request.length - 1) % line_bases + 1))
        .ok_or(ExternalError::Overflow)
}

fn strip_lines(raw: Vec<u8>, expected: u64) -> Result<Vec<u8>, ExternalError> {
    let bases = raw
        .into_iter()
        .filter(|base| *base != b'\n' && *base != b'\r')
        .map(|base| match base.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect::<Vec<_>>();
    if u64::try_from(bases.len()).unwrap_or(u64::MAX) != expected {
        return Err(ExternalError::ContigLength);
    }
    Ok(bases)
}

fn sidecar_path(path: &Path, suffix: &str) -> PathBuf {
    let mut value = OsString::from(path.as_os_str());
    value.push(suffix);
    PathBuf::from(value)
}

struct CountReader<R> {
    inner: R,
    count: Arc<AtomicU64>,
}

impl<R: Read> Read for CountReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.count
            .fetch_add(u64::try_from(read).unwrap_or(u64::MAX), Ordering::Relaxed);
        Ok(read)
    }
}

#[derive(Debug, Error)]
pub(crate) enum ExternalError {
    #[error("external sequence I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("invalid FAI record at line {0}")]
    InvalidFai(usize),
    #[error("invalid GZI index")]
    InvalidGzi,
    #[error("BGZF source has no GZI index")]
    MissingGzi,
    #[error("external sequence coordinate overflow")]
    Overflow,
    #[error("invalid external sequence coordinates")]
    InvalidCoordinates,
    #[error("external sequence payload is truncated")]
    TruncatedSequence,
    #[error("external contig length does not match metadata")]
    ContigLength,
    #[error("external contig name is not UTF-8")]
    InvalidName,
    #[error("external contig {0} does not match metadata")]
    ContigMismatch(u32),
    #[error("selected external contig is missing")]
    MissingContig,
    #[error("external sequence parse failed for {path}: {message}")]
    Parse { path: PathBuf, message: String },
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    const FASTA: &[u8] = b">one\nACGTAC\nGT\n>two\nTTTT\nCCCC\n";

    fn fai_text() -> &'static str {
        "one\t8\t5\t6\t7\ntwo\t8\t20\t4\t5\n"
    }

    fn requests() -> Vec<ContigRequest> {
        vec![
            ContigRequest {
                contig_id: 7,
                source_ordinal: 0,
                name: Some("one".to_string()),
                length: 8,
                offset: 5,
                line_bases: 6,
                line_width: 7,
            },
            ContigRequest {
                contig_id: 9,
                source_ordinal: 1,
                name: Some("two".to_string()),
                length: 8,
                offset: 20,
                line_bases: 4,
                line_width: 5,
            },
        ]
    }

    fn gzip_block(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn plain_fai_reads() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("source.fa");
        std::fs::write(&path, FASTA).unwrap();
        let fai_path = sidecar_path(&path, ".fai");
        std::fs::write(&fai_path, fai_text()).unwrap();
        assert_eq!(read_fai(fai_path).unwrap().len(), 2);
        let source = ExternalSource::detect(&path);
        assert_eq!(source.access, SequenceAccess::PlainFai);
        let loaded = read_selected(&source, &requests()).unwrap();
        assert_eq!(loaded.source_bytes, 18);
        assert_eq!(loaded.contigs[0].bases, b"ACGTACGT");
        assert_eq!(loaded.contigs[1].bases, b"TTTTCCCC");
    }

    #[test]
    fn bgzf_fai_reads() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("source.fa.gz");
        let first = gzip_block(&FASTA[..15]);
        let second = gzip_block(&FASTA[15..]);
        let mut compressed = first.clone();
        compressed.extend_from_slice(&second);
        std::fs::write(&path, compressed).unwrap();
        std::fs::write(sidecar_path(&path, ".fai"), fai_text()).unwrap();
        let mut gzi = 1u64.to_le_bytes().to_vec();
        gzi.extend_from_slice(&u64::try_from(first.len()).unwrap().to_le_bytes());
        gzi.extend_from_slice(&15u64.to_le_bytes());
        std::fs::write(sidecar_path(&path, ".gzi"), gzi).unwrap();
        let source = ExternalSource::detect(&path);
        assert_eq!(source.access, SequenceAccess::Bgzf);
        let loaded = read_selected(&source, &requests()).unwrap();
        assert_eq!(loaded.contigs[0].bases, b"ACGTACGT");
        assert_eq!(loaded.contigs[1].bases, b"TTTTCCCC");
        assert!(loaded.source_bytes < std::fs::metadata(path).unwrap().len() * 2);
    }

    #[test]
    fn gzip_stream_reads() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("source.fa.gz");
        std::fs::write(&path, gzip_block(FASTA)).unwrap();
        let source = ExternalSource::detect(&path);
        assert_eq!(source.access, SequenceAccess::Sequential);
        let loaded = read_selected(&source, &requests()[1..]).unwrap();
        assert_eq!(loaded.contigs.len(), 1);
        assert_eq!(loaded.contigs[0].bases, b"TTTTCCCC");
        assert_eq!(loaded.source_bytes, std::fs::metadata(path).unwrap().len());
    }
}
