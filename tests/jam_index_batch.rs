use flate2::Compression;
use flate2::write::GzEncoder;
use jam_rs::jam_index::{
    JamIndexBatchConfig, JamIndexBatchQuery, JamIndexBatchStatusKind, JamIndexBuildConfig,
    JamIndexTraceConfig, MetagenomeSource, build_jam_index, execute_batch, plan_batch, trace_index,
};
use jam_rs::trace::model::{QueryKind, TopologyRequested, TraceMetagenomeResult};
use sha2::{Digest, Sha256};
use std::io::Write;

fn sequence(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            b"ACGT"[(state as usize) & 3]
        })
        .collect()
}

fn write_gzip(path: &std::path::Path, records: &[(&str, &[u8])]) {
    let file = std::fs::File::create(path).unwrap();
    let mut writer = GzEncoder::new(file, Compression::fast());
    for (name, sequence) in records {
        writeln!(writer, ">{name}").unwrap();
        writer.write_all(sequence).unwrap();
        writer.write_all(b"\n").unwrap();
    }
    writer.finish().unwrap();
}

fn checksum(sequence: &[u8]) -> String {
    format!("{:x}", Sha256::digest(sequence))
}

fn query(id: &str, sequence: Vec<u8>) -> JamIndexBatchQuery {
    JamIndexBatchQuery {
        query_id: id.to_string(),
        original_header: format!("{id} retained description"),
        sequence_sha256: checksum(&sequence),
        sequence,
        query_kind: QueryKind::Plasmid,
        topology_requested: TopologyRequested::Linear,
    }
}

fn biological(result: &TraceMetagenomeResult) -> serde_json::Value {
    serde_json::json!({
        "query_id": result.plasmid_id,
        "metagenome_id": result.metagenome_id,
        "query_kind": result.query_kind,
        "topology_requested": result.topology_requested,
        "coordinate_model": result.coordinate_model,
        "topology_evidence": result.topology_evidence,
        "status": result.status,
        "candidate": result.candidate,
        "alignments": result.alignments,
        "mosaic": result.primary_fragment_mosaic,
        "topology": result.topology,
        "coverage": result.coverage,
        "warnings": result.warnings,
        "failures": result.failures,
    })
}

#[test]
fn batch_reuses_gzip_group_and_matches_single_query_results() {
    let sources = tempfile::Builder::new()
        .prefix("jam-index-batch-sources-")
        .tempdir_in("target")
        .unwrap();
    let outputs = tempfile::Builder::new()
        .prefix("jam-index-batch-output-")
        .tempdir_in("target")
        .unwrap();
    let first = sequence(640, 0x1234_5678_9abc_def0);
    let second = sequence(720, 0xfedc_ba98_7654_3210);
    let unrelated = sequence(560, 0x0ddc_0ffe_e123_4567);
    let background = sequence(2_000, 0xa5a5_5a5a_1234_4321);
    let positive_path = sources.path().join("positive.fa.gz");
    let negative_path = sources.path().join("negative.fa.gz");
    write_gzip(
        &positive_path,
        &[
            ("noise", &background),
            ("first", &first),
            ("second", &second),
        ],
    );
    write_gzip(&negative_path, &[("negative", &unrelated)]);
    let index = outputs.path().join("index");
    build_jam_index(
        &index,
        &[
            MetagenomeSource {
                metagenome_id: "positive".to_string(),
                sequence_path: positive_path,
            },
            MetagenomeSource {
                metagenome_id: "negative".to_string(),
                sequence_path: negative_path,
            },
        ],
        &JamIndexBuildConfig {
            target_parts: 2,
            parallel_parts: 2,
            source_manifest_sha256: "a".repeat(64),
            ..JamIndexBuildConfig::default()
        },
    )
    .unwrap();
    let queries = vec![
        query("query-first", first.clone()),
        query("query-second", second.clone()),
        query("query-absent", sequence(512, 0x9999_8888_7777_6666)),
    ];
    let config = JamIndexBatchConfig {
        trace: JamIndexTraceConfig {
            parallel_candidates: 2,
            ..JamIndexTraceConfig::default()
        },
        max_group_contig_bases: 16 * 1024 * 1024,
        fallback_contigs_per_chunk: 2,
        ..JamIndexBatchConfig::default()
    };
    let plan = plan_batch(&index, &queries, &config).unwrap();
    assert_eq!(plan.screen_parts_opened, 2);
    assert_eq!(
        plan.groups.len(),
        1,
        "noncandidate source must not be grouped"
    );
    assert!(
        plan.groups
            .iter()
            .any(|group| { group.metagenome_id == "positive" && group.work.len() == 2 })
    );
    let mut emitted = Vec::new();
    let execution = execute_batch(&plan, &config, "batch-run", |result| {
        emitted.push(result);
        Ok(())
    })
    .unwrap();
    assert_eq!(execution.metrics.source_open_count, 1);
    assert_eq!(
        execution
            .statuses
            .iter()
            .find(|status| status.query_id == "query-absent")
            .unwrap()
            .status,
        JamIndexBatchStatusKind::NoCandidate
    );
    for input in &queries[..2] {
        let single = trace_index(&index, &input.query_id, &input.sequence, &config.trace)
            .unwrap()
            .metagenomes
            .into_iter()
            .find(|result| result.metagenome_id == "positive")
            .unwrap();
        let batch = emitted
            .iter()
            .find(|result| {
                result.plasmid_id == input.query_id && result.metagenome_id == "positive"
            })
            .unwrap();
        assert_eq!(biological(batch), biological(&single));
    }

    let many = (0..100)
        .map(|index| query(&format!("query-{index:03}"), first.clone()))
        .collect::<Vec<_>>();
    let many_plan = plan_batch(&index, &many, &config).unwrap();
    assert_eq!(many_plan.screen_parts_opened, 2);
    assert_eq!(many_plan.groups.len(), 1);
    let mut emitted = 0u64;
    let many_execution = execute_batch(&many_plan, &config, "many-run", |_| {
        emitted += 1;
        Ok(())
    })
    .unwrap();
    assert_eq!(emitted, 100);
    assert_eq!(many_execution.statuses.len(), 100);
    assert_eq!(many_execution.metrics.source_open_count, 1);
}
