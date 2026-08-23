use jam_rs::resource::ResourceMetrics;
use jam_rs::trace::config::SensitivityConfig;
use jam_rs::trace::model::{
    InputResource, QueryKind, TopologyRequested, TraceRunFooter, TraceRunHeader,
};
use jam_rs::trace::output::{SCHEMA_VERSION, TraceFileWriter, TraceJsonlWriter, TraceOutputError};
use std::fs;
use std::io::{self, Read, Write};

fn header() -> TraceRunHeader {
    TraceRunHeader {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "atomic-run".to_string(),
        jam_rs_version: "0.10.0".to_string(),
        source_commit: Some("test".to_string()),
        started_at_utc: "2026-01-01T00:00:00Z".to_string(),
        command: vec!["jam".to_string(), "trace".to_string()],
        plasmid_id: "plasmid-1".to_string(),
        plasmid_length: 12,
        query_kind: QueryKind::Plasmid,
        topology_requested: TopologyRequested::Auto,
        threads: 1,
        io_concurrency: 1,
        sensitivity: SensitivityConfig::default(),
        algorithms: jam_rs::trace::algorithm_identifiers(),
        algorithm: jam_rs::trace::TraceAlgorithmMetadata::for_sensitivity(
            SensitivityConfig::default(),
        ),
        inputs: vec![InputResource {
            role: "plasmid".to_string(),
            redacted_locator: "file:///plasmid.fa".to_string(),
            sha256: None,
        }],
    }
}

fn footer() -> TraceRunFooter {
    TraceRunFooter {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "atomic-run".to_string(),
        completed_at_utc: "2026-01-01T00:00:01Z".to_string(),
        metagenomes_total: 0,
        metagenomes_with_candidates: 0,
        metagenomes_aligned: 0,
        metagenomes_failed: 0,
        alignments_total: 0,
        resource_metrics: ResourceMetrics::default(),
    }
}

fn temp_entry_count(path: &std::path::Path) -> usize {
    fs::read_dir(path).unwrap().count()
}

#[test]
fn plain_file_is_published_only_after_successful_finish() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("trace.jsonl");
    let mut writer = TraceFileWriter::create(&path).unwrap();

    assert!(!path.exists());
    assert_eq!(temp_entry_count(directory.path()), 1);
    writer.write_header(&header()).unwrap();
    writer.write_footer(&footer()).unwrap();
    writer.flush().unwrap();
    assert!(!path.exists());

    writer.finish().unwrap();
    assert!(path.is_file());
    assert_eq!(temp_entry_count(directory.path()), 1);
    let text = fs::read_to_string(path).unwrap();
    assert_eq!(text.lines().count(), 2);
}

#[test]
fn successful_finish_replaces_existing_output_only_at_commit() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("trace.jsonl");
    let original = b"previous complete output\n";
    fs::write(&path, original).unwrap();

    let mut writer = TraceFileWriter::create(&path).unwrap();
    writer.write_header(&header()).unwrap();
    writer.write_footer(&footer()).unwrap();
    assert_eq!(fs::read(&path).unwrap(), original);

    writer.finish().unwrap();
    assert_ne!(fs::read(&path).unwrap(), original);
    assert_eq!(temp_entry_count(directory.path()), 1);
}

#[test]
fn zstd_file_is_published_only_after_encoder_finishes() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("trace.jsonl.zst");
    let mut writer = TraceFileWriter::create(&path).unwrap();
    writer.write_header(&header()).unwrap();
    writer.write_footer(&footer()).unwrap();
    assert!(!path.exists());

    writer.finish().unwrap();
    assert_eq!(temp_entry_count(directory.path()), 1);
    let file = fs::File::open(path).unwrap();
    let mut decoder = zstd::stream::read::Decoder::new(file).unwrap();
    let mut text = String::new();
    decoder.read_to_string(&mut text).unwrap();
    assert_eq!(text.lines().count(), 2);
}

#[test]
fn failed_stream_preserves_existing_output_and_removes_temp_file() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("trace.jsonl");
    let original = b"previous complete output\n";
    fs::write(&path, original).unwrap();

    let mut writer = TraceFileWriter::create(&path).unwrap();
    writer.write_header(&header()).unwrap();
    assert_eq!(fs::read(&path).unwrap(), original);
    assert_eq!(temp_entry_count(directory.path()), 2);
    assert!(matches!(
        writer.finish(),
        Err(TraceOutputError::FooterRequired)
    ));

    assert_eq!(fs::read(&path).unwrap(), original);
    assert_eq!(temp_entry_count(directory.path()), 1);
}

#[test]
fn failed_final_rename_preserves_existing_directory_and_removes_temp_file() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("trace-output");
    fs::create_dir(&path).unwrap();

    let mut writer = TraceFileWriter::create(&path).unwrap();
    writer.write_header(&header()).unwrap();
    writer.write_footer(&footer()).unwrap();
    let error = writer.finish().err().unwrap();
    assert!(matches!(error, TraceOutputError::Io(_)));
    assert!(path.is_dir());
    assert_eq!(temp_entry_count(directory.path()), 1);
}

struct FailingWriter {
    writes_before_failure: usize,
}

impl Write for FailingWriter {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.writes_before_failure == 0 {
            return Err(io::Error::other("injected output failure"));
        }
        self.writes_before_failure -= 1;
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn injected_writer_failure_is_returned_without_advancing_stream_state() {
    let mut writer = TraceJsonlWriter::new(FailingWriter {
        writes_before_failure: 0,
    });
    let error = writer.write_header(&header()).unwrap_err();
    assert!(matches!(error, TraceOutputError::Io(_)));
    assert_eq!(writer.records_written(), 0);
    assert!(!writer.is_finished());
}
