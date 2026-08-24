use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexTraceConfig, MetagenomeSource, build_jam_index, trace_index,
};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use jam_rs::trace::model::{Strand, TopologyRequested};
use needletail::Sequence;
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
    let reverse_source = sources.path().join("reverse.fasta");
    let rotated_source = sources.path().join("rotated.fasta");
    let reverse = query.reverse_complement();
    let split = 53;
    let mut rotated = query[split..].to_vec();
    rotated.extend_from_slice(&query[..split]);
    fs::write(
        &source,
        format!(">target\n{}\n", String::from_utf8_lossy(query)),
    )
    .unwrap();
    fs::write(
        &reverse_source,
        format!(">reverse\n{}\n", String::from_utf8_lossy(&reverse)),
    )
    .unwrap();
    fs::write(
        &rotated_source,
        format!(">rotated\n{}\n", String::from_utf8_lossy(&rotated)),
    )
    .unwrap();
    let index = outputs.path().join("index");
    build_jam_index(
        &index,
        &[
            MetagenomeSource {
                metagenome_id: "target".into(),
                sequence_path: source,
            },
            MetagenomeSource {
                metagenome_id: "reverse".into(),
                sequence_path: reverse_source,
            },
            MetagenomeSource {
                metagenome_id: "rotated".into(),
                sequence_path: rotated_source,
            },
        ],
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
    config.runner.topology_requested = TopologyRequested::Circular;
    let output = trace_index(&index, "query", query, &config).unwrap();
    assert_eq!(output.metagenomes.len(), 3);
    for result in &output.metagenomes {
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
    }
    let reverse = output
        .metagenomes
        .iter()
        .find(|result| result.metagenome_id == "reverse")
        .unwrap();
    assert!(
        reverse
            .alignments
            .iter()
            .any(|alignment| alignment.strand == Strand::Reverse)
    );
    let rotated = output
        .metagenomes
        .iter()
        .find(|result| result.metagenome_id == "rotated")
        .unwrap();
    assert!(
        rotated
            .alignments
            .iter()
            .any(|alignment| alignment.origin_crossing && alignment.query_segments.len() == 2)
    );
    assert_eq!(output.trace_metrics.selected_contigs, 3);
}
