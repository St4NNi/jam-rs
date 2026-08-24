use jam_rs::resource::ResourceMetrics;
use jam_rs::trace::config::SensitivityConfig;
use jam_rs::trace::model::{
    CandidatePerformanceCounters, CoordinateModel, InputResource, QueryKind, TopologyEvidence,
    TopologyRequested, TraceMetagenomeResult, TraceRecord, TraceRunFooter, TraceRunHeader,
    TraceStatus,
};
use jam_rs::trace::output::{SCHEMA_VERSION, TraceFileWriter, TraceJsonlWriter, TraceOutputError};
use serde_json::Value;
use std::io::{Cursor, Read};

fn header() -> TraceRunHeader {
    TraceRunHeader {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "fixture-run".to_string(),
        jam_rs_version: "0.10.0".to_string(),
        source_commit: Some("fixture-commit".to_string()),
        started_at_utc: "2026-01-01T00:00:00Z".to_string(),
        command: vec!["jam".to_string(), "trace".to_string()],
        plasmid_id: "plasmid-fixture".to_string(),
        plasmid_length: 32,
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
            redacted_locator: "https://example.org/catalog/plasmid.fa".to_string(),
            sha256: Some("abc".to_string()),
        }],
    }
}

fn result() -> TraceMetagenomeResult {
    TraceMetagenomeResult {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "fixture-run".to_string(),
        plasmid_id: "plasmid-fixture".to_string(),
        metagenome_id: "metagenome-fixture".to_string(),
        query_kind: QueryKind::Plasmid,
        topology_requested: TopologyRequested::Auto,
        coordinate_model: CoordinateModel::Undetermined,
        topology_evidence: TopologyEvidence::Insufficient,
        algorithms: jam_rs::trace::algorithm_identifiers(),
        algorithm: jam_rs::trace::TraceAlgorithmMetadata::for_sensitivity(
            SensitivityConfig::default(),
        ),
        status: TraceStatus::NoCandidate,
        candidate: None,
        alignments: Vec::new(),
        primary_fragment_mosaic: None,
        topology: None,
        rescue_rounds: Vec::new(),
        stages: Vec::new(),
        alignment_retries: Vec::new(),
        performance_counters: CandidatePerformanceCounters::default(),
        coverage: None,
        warnings: Vec::new(),
        failures: Vec::new(),
        archive_metrics: None,
        resource_metrics: ResourceMetrics::default(),
    }
}

fn footer() -> TraceRunFooter {
    TraceRunFooter {
        schema_version: SCHEMA_VERSION.to_string(),
        run_id: "fixture-run".to_string(),
        completed_at_utc: "2026-01-01T00:00:01Z".to_string(),
        metagenomes_total: 1,
        metagenomes_with_candidates: 0,
        metagenomes_aligned: 0,
        metagenomes_failed: 0,
        alignments_total: 0,
        resource_metrics: ResourceMetrics::default(),
    }
}

fn record_values() -> Vec<Value> {
    [
        TraceRecord::RunHeader(header()),
        TraceRecord::MetagenomeResult(result()),
        TraceRecord::RunFooter(footer()),
    ]
    .into_iter()
    .map(|record| serde_json::to_value(record).unwrap())
    .collect()
}

#[test]
fn schema_fixture_declares_all_record_kinds_and_shared_fields() {
    let schema: Value =
        serde_json::from_str(include_str!("../schemas/jam-trace-v2.schema.json")).unwrap();
    assert_eq!(
        schema["$schema"],
        "https://json-schema.org/draft/2020-12/schema"
    );
    assert_eq!(
        schema["$defs"]["alignment"]["properties"]["cigar"]["type"],
        "string"
    );
    assert_eq!(
        schema["$defs"]["coverage"]["properties"]["gaps"]["type"],
        "array"
    );
    assert_eq!(
        schema["$defs"]["algorithm_metadata"]["properties"]["id"]["const"],
        jam_rs::trace::TRACE_ALGORITHM_ID
    );
    assert_eq!(
        schema["$defs"]["algorithm_identifiers"]["properties"]["trace_workflow"]["const"],
        jam_rs::trace::config::TRACE_WORKFLOW_ID
    );
    assert_eq!(
        schema["$defs"]["fragment_mosaic"]["properties"]["mosaic_algorithm"]["const"],
        jam_rs::trace::config::MOSAIC_ALGORITHM_ID
    );
    let kinds = schema["properties"]["record_type"]["enum"]
        .as_array()
        .unwrap();
    assert_eq!(
        kinds,
        &[
            Value::String("run_header".to_string()),
            Value::String("metagenome_result".to_string()),
            Value::String("run_footer".to_string()),
        ]
    );
}

#[test]
fn serialized_fixture_records_have_schema_and_required_payloads() {
    let values = record_values();
    assert_eq!(values.len(), 3);
    assert_eq!(values[0]["record_type"], "run_header");
    assert_eq!(values[1]["record_type"], "metagenome_result");
    assert_eq!(values[2]["record_type"], "run_footer");
    for value in &values {
        assert_eq!(value["schema_version"], SCHEMA_VERSION);
        assert_eq!(value["run_id"], "fixture-run");
    }
    assert_eq!(
        values[0]["inputs"][0]["redacted_locator"],
        "https://example.org/catalog/plasmid.fa"
    );
    assert_eq!(
        values[0]["algorithm"]["id"],
        jam_rs::trace::TRACE_ALGORITHM_ID
    );
    assert_eq!(values[0]["algorithm"], values[1]["algorithm"]);
    assert_eq!(values[0]["query_id"], "plasmid-fixture");
    assert!(values[0].get("plasmid_id").is_none());
    assert_eq!(values[0]["query_kind"], "plasmid");
    assert_eq!(values[0]["topology_requested"], "auto");
    assert_eq!(values[1]["coordinate_model"], "undetermined");
    assert_eq!(values[0]["algorithms"], values[1]["algorithms"]);
    assert_eq!(
        values[0]["algorithm"]["parameters"]["alignment_mode"],
        "local_affine_chain_corridor_with_semiglobal_refinement"
    );
    assert!(values[1].get("alignments").unwrap().is_array());
    assert!(values[2].get("resource_metrics").unwrap().is_object());
}

#[test]
fn plain_jsonl_writer_is_incremental_and_stable() {
    let mut writer = TraceJsonlWriter::new(Cursor::new(Vec::<u8>::new()));
    writer.write_header(&header()).unwrap();
    writer.write_metagenome_result(&result()).unwrap();
    writer.write_footer(&footer()).unwrap();
    assert_eq!(writer.records_written(), 3);
    assert!(writer.is_finished());
    let bytes = writer.finish().unwrap().into_inner();
    let lines: Vec<_> = String::from_utf8(bytes)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .collect();
    assert_eq!(
        lines
            .iter()
            .map(|value| value["record_type"].as_str().unwrap())
            .collect::<Vec<_>>(),
        vec!["run_header", "metagenome_result", "run_footer"]
    );
}

#[test]
fn malformed_stream_metadata_is_rejected_before_writing() {
    let mut writer = TraceJsonlWriter::new(Cursor::new(Vec::<u8>::new()));
    let mut wrong = header();
    wrong.schema_version = "0.9.0".to_string();
    assert!(matches!(
        writer.write_header(&wrong),
        Err(TraceOutputError::SchemaVersionMismatch { .. })
    ));
    assert_eq!(writer.records_written(), 0);
}

#[test]
fn zstd_file_writer_round_trips_the_same_record_sequence() {
    let temporary = tempfile::tempdir().unwrap();
    let path = temporary.path().join("trace.jsonl.zst");
    let mut writer = TraceFileWriter::create(&path).unwrap();
    writer.write_header(&header()).unwrap();
    writer.write_metagenome_result(&result()).unwrap();
    writer.write_footer(&footer()).unwrap();
    writer.finish().unwrap();

    let file = std::fs::File::open(path).unwrap();
    let mut decoder = zstd::stream::read::Decoder::new(file).unwrap();
    let mut text = String::new();
    decoder.read_to_string(&mut text).unwrap();
    assert_eq!(text.lines().count(), 3);
    assert_eq!(
        text.lines()
            .map(|line| serde_json::from_str::<Value>(line).unwrap()["record_type"].clone())
            .collect::<Vec<_>>(),
        [
            Value::String("run_header".to_string()),
            Value::String("metagenome_result".to_string()),
            Value::String("run_footer".to_string())
        ]
    );
}
