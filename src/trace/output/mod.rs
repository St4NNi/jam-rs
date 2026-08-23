//! Incremental JSONL output for the trace workflow.
//!
//! The writer deliberately owns only one serialized record at a time. This
//! keeps a large one-plasmid-to-many-metagenome run bounded by the size of its
//! largest result rather than by the number of metagenomes in the run.

use crate::trace::TRACE_JSON_SCHEMA_VERSION;
use crate::trace::model::{TraceMetagenomeResult, TraceRecord, TraceRunFooter, TraceRunHeader};
use serde::Serialize;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;
use thiserror::Error;

/// The JSONL schema emitted by this module.
pub const SCHEMA_VERSION: &str = TRACE_JSON_SCHEMA_VERSION;

/// Errors raised while validating or writing a trace stream.
#[derive(Debug, Error)]
pub enum TraceOutputError {
    #[error("trace output I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("trace record serialization failed: {0}")]
    Json(#[from] serde_json::Error),
    #[error("trace output must begin with a run_header record")]
    HeaderRequired,
    #[error("trace output already contains a run_header record")]
    HeaderAlreadyWritten,
    #[error("trace output already contains a run_footer record")]
    FooterAlreadyWritten,
    #[error("trace output cannot be finished before a run_footer record")]
    FooterRequired,
    #[error("trace record has schema version {found:?}; expected {expected:?}")]
    SchemaVersionMismatch {
        expected: &'static str,
        found: String,
    },
    #[error("trace record has run ID {found:?}; expected {expected:?}")]
    RunIdMismatch { expected: String, found: String },
    #[error("trace output path has no file name: {0:?}")]
    MissingFileName(PathBuf),
}

pub type TraceOutputResult<T> = Result<T, TraceOutputError>;

#[derive(Clone, Debug, Default)]
struct StreamState {
    run_id: Option<String>,
    header_written: bool,
    footer_written: bool,
    records_written: u64,
}

/// An incremental JSONL writer over any `Write` implementation.
///
/// `W` can be a regular file, an in-memory buffer used by tests, or a
/// `zstd::stream::write::Encoder`. The writer emits one newline-terminated
/// JSON object per call and flushes after every record so a long-running
/// search leaves a usable prefix if a later metagenome fails.
pub struct TraceJsonlWriter<W: Write> {
    writer: BufWriter<W>,
    state: StreamState,
}

impl<W: Write> TraceJsonlWriter<W> {
    /// Construct a writer with a bounded output buffer.
    #[must_use]
    pub fn new(writer: W) -> Self {
        Self {
            writer: BufWriter::with_capacity(1024 * 1024, writer),
            state: StreamState::default(),
        }
    }

    /// Write one validated record and flush it to the underlying writer.
    pub fn write_record(&mut self, record: &TraceRecord) -> TraceOutputResult<()> {
        self.validate_record(record)?;

        // Serialize before touching the stream. This makes a serialization
        // error unable to leave a truncated JSON line behind.
        let mut encoded = serde_json::to_vec(record)?;
        encoded.push(b'\n');
        self.writer.write_all(&encoded)?;
        self.writer.flush()?;

        match record {
            TraceRecord::RunHeader(header) => {
                self.state.run_id = Some(header.run_id.clone());
                self.state.header_written = true;
            }
            TraceRecord::RunFooter(_) => self.state.footer_written = true,
            TraceRecord::MetagenomeResult(_) => {}
        }
        self.state.records_written += 1;
        Ok(())
    }

    /// Write the required first record.
    pub fn write_header(&mut self, header: &TraceRunHeader) -> TraceOutputResult<()> {
        self.write_record(&TraceRecord::RunHeader(header.clone()))
    }

    /// Write one per-metagenome result.
    pub fn write_metagenome_result(
        &mut self,
        result: &TraceMetagenomeResult,
    ) -> TraceOutputResult<()> {
        self.write_record(&TraceRecord::MetagenomeResult(result.clone()))
    }

    /// Write the required final record.
    pub fn write_footer(&mut self, footer: &TraceRunFooter) -> TraceOutputResult<()> {
        self.write_record(&TraceRecord::RunFooter(footer.clone()))
    }

    /// Flush buffered bytes without changing stream state.
    pub fn flush(&mut self) -> TraceOutputResult<()> {
        self.writer.flush().map_err(TraceOutputError::from)
    }

    /// Finish the JSONL stream and recover the underlying writer.
    ///
    /// A normal completed stream must contain both a header and a footer. A
    /// caller handling an interrupted run can intentionally drop the writer
    /// after flushing its last successful result instead.
    pub fn finish(mut self) -> TraceOutputResult<W> {
        if !self.state.header_written {
            return Err(TraceOutputError::HeaderRequired);
        }
        if !self.state.footer_written {
            return Err(TraceOutputError::FooterRequired);
        }
        self.writer.flush()?;
        self.writer
            .into_inner()
            .map_err(|error| TraceOutputError::Io(error.into_error()))
    }

    /// Number of records successfully written to this stream.
    #[must_use]
    pub fn records_written(&self) -> u64 {
        self.state.records_written
    }

    /// Whether a complete header/footer stream has been written.
    #[must_use]
    pub fn is_finished(&self) -> bool {
        self.state.header_written && self.state.footer_written
    }

    fn validate_record(&self, record: &TraceRecord) -> TraceOutputResult<()> {
        let (schema_version, run_id) = match record {
            TraceRecord::RunHeader(header) => {
                if self.state.header_written {
                    return Err(TraceOutputError::HeaderAlreadyWritten);
                }
                (&header.schema_version, &header.run_id)
            }
            TraceRecord::MetagenomeResult(result) => {
                if !self.state.header_written {
                    return Err(TraceOutputError::HeaderRequired);
                }
                if self.state.footer_written {
                    return Err(TraceOutputError::FooterAlreadyWritten);
                }
                (&result.schema_version, &result.run_id)
            }
            TraceRecord::RunFooter(footer) => {
                if !self.state.header_written {
                    return Err(TraceOutputError::HeaderRequired);
                }
                if self.state.footer_written {
                    return Err(TraceOutputError::FooterAlreadyWritten);
                }
                (&footer.schema_version, &footer.run_id)
            }
        };

        if schema_version != SCHEMA_VERSION {
            return Err(TraceOutputError::SchemaVersionMismatch {
                expected: SCHEMA_VERSION,
                found: schema_version.clone(),
            });
        }
        if let Some(expected) = &self.state.run_id
            && expected != run_id
        {
            return Err(TraceOutputError::RunIdMismatch {
                expected: expected.clone(),
                found: run_id.clone(),
            });
        }
        Ok(())
    }
}

/// A file-backed trace writer. A `.zst` or `.zstd` suffix selects a
/// Zstandard stream; every other suffix writes plain UTF-8 JSONL.
pub struct TraceFileWriter {
    inner: TraceFileWriterInner,
}

enum TraceFileWriterInner {
    Plain {
        writer: TraceJsonlWriter<File>,
        output: AtomicOutput,
    },
    Zstd {
        writer: TraceJsonlWriter<zstd::stream::write::Encoder<'static, File>>,
        output: AtomicOutput,
    },
}

/// A temporary output in the destination directory that is published only
/// after a complete JSONL stream has been finalized successfully.
struct AtomicOutput {
    destination: PathBuf,
    temporary: Option<NamedTempFile>,
    committed: bool,
}

impl AtomicOutput {
    fn create(destination: &Path) -> TraceOutputResult<(Self, File)> {
        let parent = destination
            .parent()
            .filter(|parent| !parent.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."));
        let temporary = tempfile::Builder::new()
            .prefix(".jam-trace-output-")
            .tempfile_in(parent)?;
        let file = temporary.reopen()?;
        Ok((
            Self {
                destination: destination.to_path_buf(),
                temporary: Some(temporary),
                committed: false,
            },
            file,
        ))
    }

    fn commit(mut self, file: File) -> TraceOutputResult<()> {
        file.sync_all()?;
        drop(file);
        let temporary = self.temporary.take().ok_or_else(|| {
            TraceOutputError::Io(io::Error::other("atomic output temporary file is missing"))
        })?;
        match temporary.persist(&self.destination) {
            Ok(_) => {
                self.committed = true;
                Ok(())
            }
            Err(error) => {
                self.temporary = Some(error.file);
                Err(error.error.into())
            }
        }
    }
}

impl Drop for AtomicOutput {
    fn drop(&mut self) {
        if self.committed {
            return;
        }
        // Dropping the retained NamedTempFile removes the unpublished file.
        let _ = self.temporary.take();
    }
}

impl TraceFileWriter {
    /// Open a plain or Zstandard-compressed JSONL output file.
    pub fn create(path: impl AsRef<Path>) -> TraceOutputResult<Self> {
        let path = path.as_ref();
        let (output, file) = AtomicOutput::create(path)?;
        let is_zstd = path
            .extension()
            .and_then(|extension| extension.to_str())
            .is_some_and(|extension| {
                extension.eq_ignore_ascii_case("zst") || extension.eq_ignore_ascii_case("zstd")
            });
        if is_zstd {
            let encoder = zstd::stream::write::Encoder::new(file, 0)?;
            Ok(Self {
                inner: TraceFileWriterInner::Zstd {
                    writer: TraceJsonlWriter::new(encoder),
                    output,
                },
            })
        } else {
            Ok(Self {
                inner: TraceFileWriterInner::Plain {
                    writer: TraceJsonlWriter::new(file),
                    output,
                },
            })
        }
    }

    pub fn write_record(&mut self, record: &TraceRecord) -> TraceOutputResult<()> {
        match &mut self.inner {
            TraceFileWriterInner::Plain { writer, .. } => writer.write_record(record),
            TraceFileWriterInner::Zstd { writer, .. } => writer.write_record(record),
        }
    }

    pub fn write_header(&mut self, header: &TraceRunHeader) -> TraceOutputResult<()> {
        match &mut self.inner {
            TraceFileWriterInner::Plain { writer, .. } => writer.write_header(header),
            TraceFileWriterInner::Zstd { writer, .. } => writer.write_header(header),
        }
    }

    pub fn write_metagenome_result(
        &mut self,
        result: &TraceMetagenomeResult,
    ) -> TraceOutputResult<()> {
        match &mut self.inner {
            TraceFileWriterInner::Plain { writer, .. } => writer.write_metagenome_result(result),
            TraceFileWriterInner::Zstd { writer, .. } => writer.write_metagenome_result(result),
        }
    }

    pub fn write_footer(&mut self, footer: &TraceRunFooter) -> TraceOutputResult<()> {
        match &mut self.inner {
            TraceFileWriterInner::Plain { writer, .. } => writer.write_footer(footer),
            TraceFileWriterInner::Zstd { writer, .. } => writer.write_footer(footer),
        }
    }

    pub fn flush(&mut self) -> TraceOutputResult<()> {
        match &mut self.inner {
            TraceFileWriterInner::Plain { writer, .. } => writer.flush(),
            TraceFileWriterInner::Zstd { writer, .. } => writer.flush(),
        }
    }

    /// Finish the JSONL record stream and, for compressed output, the zstd
    /// frame as well.
    pub fn finish(self) -> TraceOutputResult<()> {
        match self.inner {
            TraceFileWriterInner::Plain { writer, output } => {
                let file = writer.finish()?;
                output.commit(file)
            }
            TraceFileWriterInner::Zstd { writer, output } => {
                let encoder = writer.finish()?;
                let file = encoder.finish()?;
                output.commit(file)
            }
        }
    }
}

/// Convenience constructor for path-based callers.
pub fn create(path: impl AsRef<Path>) -> TraceOutputResult<TraceFileWriter> {
    let path = path.as_ref();
    if path.file_name().is_none() {
        return Err(TraceOutputError::MissingFileName(path.to_path_buf()));
    }
    TraceFileWriter::create(path)
}

/// Serialize one record to a string for logging or a caller-owned stream.
/// The returned string is newline terminated like the file writer output.
pub fn record_json<T: Serialize>(record: &T) -> TraceOutputResult<String> {
    let mut line = serde_json::to_string(record)?;
    line.push('\n');
    Ok(line)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resource::ResourceMetrics;
    use crate::trace::config::SensitivityConfig;
    use crate::trace::model::{
        CandidatePerformanceCounters, CoordinateModel, InputResource, QueryKind, TopologyEvidence,
        TopologyRequested, TraceStatus,
    };
    use std::io::Cursor;

    fn header() -> TraceRunHeader {
        TraceRunHeader {
            schema_version: SCHEMA_VERSION.to_string(),
            run_id: "run-1".to_string(),
            jam_rs_version: "0.10.0".to_string(),
            source_commit: Some("abc".to_string()),
            started_at_utc: "2026-01-01T00:00:00Z".to_string(),
            command: vec!["jam".to_string(), "trace".to_string()],
            plasmid_id: "plasmid-1".to_string(),
            plasmid_length: 12,
            query_kind: QueryKind::Plasmid,
            topology_requested: TopologyRequested::Auto,
            threads: 1,
            io_concurrency: 1,
            sensitivity: SensitivityConfig::default(),
            algorithms: crate::trace::algorithm_identifiers(),
            algorithm: crate::trace::TraceAlgorithmMetadata::for_sensitivity(
                SensitivityConfig::default(),
            ),
            inputs: vec![InputResource {
                role: "plasmid".to_string(),
                redacted_locator: "file:///plasmid.fa".to_string(),
                sha256: None,
            }],
        }
    }

    fn result() -> TraceMetagenomeResult {
        TraceMetagenomeResult {
            schema_version: SCHEMA_VERSION.to_string(),
            run_id: "run-1".to_string(),
            plasmid_id: "plasmid-1".to_string(),
            metagenome_id: "sample-1".to_string(),
            query_kind: QueryKind::Plasmid,
            topology_requested: TopologyRequested::Auto,
            coordinate_model: CoordinateModel::Undetermined,
            topology_evidence: TopologyEvidence::Insufficient,
            algorithms: crate::trace::algorithm_identifiers(),
            algorithm: crate::trace::TraceAlgorithmMetadata::for_sensitivity(
                SensitivityConfig::default(),
            ),
            status: TraceStatus::NoCandidate,
            candidate: None,
            alignments: Vec::new(),
            primary_fragment_mosaic: None,
            topology: None,
            rescue_rounds: Vec::new(),
            performance_counters: CandidatePerformanceCounters::default(),
            coverage: None,
            warnings: Vec::new(),
            failures: Vec::new(),
            resource_metrics: ResourceMetrics::default(),
        }
    }

    fn footer() -> TraceRunFooter {
        TraceRunFooter {
            schema_version: SCHEMA_VERSION.to_string(),
            run_id: "run-1".to_string(),
            completed_at_utc: "2026-01-01T00:00:01Z".to_string(),
            metagenomes_total: 1,
            metagenomes_with_candidates: 0,
            metagenomes_aligned: 0,
            metagenomes_failed: 0,
            alignments_total: 0,
            resource_metrics: ResourceMetrics::default(),
        }
    }

    #[test]
    fn writes_ordered_newline_delimited_records_incrementally() {
        let mut writer = TraceJsonlWriter::new(Cursor::new(Vec::new()));
        writer.write_header(&header()).unwrap();
        writer.write_metagenome_result(&result()).unwrap();
        writer.write_footer(&footer()).unwrap();
        let output = writer.finish().unwrap().into_inner();
        let lines: Vec<_> = String::from_utf8(output)
            .unwrap()
            .lines()
            .map(serde_json::from_str::<serde_json::Value>)
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0]["record_type"], "run_header");
        assert_eq!(lines[1]["record_type"], "metagenome_result");
        assert_eq!(lines[2]["record_type"], "run_footer");
    }

    #[test]
    fn enforces_header_footer_and_run_identity() {
        let mut writer = TraceJsonlWriter::new(Cursor::new(Vec::new()));
        assert_eq!(
            writer
                .write_metagenome_result(&result())
                .unwrap_err()
                .to_string(),
            "trace output must begin with a run_header record"
        );
        writer.write_header(&header()).unwrap();
        let mut wrong = result();
        wrong.run_id = "other".to_string();
        assert!(matches!(
            writer.write_metagenome_result(&wrong),
            Err(TraceOutputError::RunIdMismatch { .. })
        ));
        assert!(matches!(
            writer.finish(),
            Err(TraceOutputError::FooterRequired)
        ));
    }

    #[test]
    fn serializes_schema_version_and_record_type() {
        let json = record_json(&TraceRecord::RunHeader(header())).unwrap();
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(value["record_type"], "run_header");
        assert_eq!(value["schema_version"], SCHEMA_VERSION);
    }
}
