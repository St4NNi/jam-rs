use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexTraceConfig, MetagenomeSource, build_jam_index, trace_index,
};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use std::fs;
use tempfile::Builder;

#[test]
fn exact_trace_works() {
    let sources = Builder::new()
        .prefix("index-trace-source-")
        .tempdir_in("target")
        .unwrap();
    let outputs = Builder::new()
        .prefix("index-trace-output-")
        .tempdir_in("target")
        .unwrap();
    let query = b"ACGTTGCATGTCAGTACGATCGTACGTTAGCTAGCTGACTGATCGTAGCTAGTCGATCGTACGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTAGCTGATCGTAGCTAGTCGATGCTA";
    let source = sources.path().join("target.fasta");
    fs::write(
        &source,
        format!(">target\n{}\n", String::from_utf8_lossy(query)),
    )
    .unwrap();
    let index = outputs.path().join("index");
    build_jam_index(
        &index,
        &[MetagenomeSource {
            metagenome_id: "target".into(),
            sequence_path: source,
        }],
        &JamIndexBuildConfig {
            source_manifest_sha256: "a".repeat(64),
            ..JamIndexBuildConfig::default()
        },
    )
    .unwrap();
    let mut config = JamIndexTraceConfig::default();
    config.runner.sensitivity = SensitivityConfig::for_profile(SensitivityProfile::Sensitive);
    config.runner.threads = 1;
    config.runner.io_concurrency = 1;
    config.runner.memory_budget_bytes = 512 * 1024 * 1024;
    let output = trace_index(&index, "query", query, &config).unwrap();
    assert_eq!(output.metagenomes.len(), 1);
    let result = &output.metagenomes[0];
    assert_eq!(result.metagenome_id, "target");
    assert_eq!(
        result.coverage.as_ref().unwrap().supported_bases,
        query.len() as u64
    );
    assert!(
        result
            .alignments
            .iter()
            .any(|alignment| alignment.cigar == format!("{}=", query.len()))
    );
    assert_eq!(result.archive_metrics.as_ref().unwrap().seed_bytes_read, 0);
    assert_eq!(output.trace_metrics.selected_contigs, 1);
}
