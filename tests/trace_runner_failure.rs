use jam_rs::resource::ResourceOpenOptions;
use jam_rs::trace::catalog::{CatalogEntry, TraceCatalog};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use jam_rs::trace::runner::{TraceQuery, TraceRunner, TraceRunnerConfig};
use jam_rs::trace::screen::CandidateSearchConfig;
use std::fs;
use std::process::Command;

#[path = "support/remote/mod.rs"]
mod remote;

fn dna(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            b"ACGT"[((state >> 62) & 3) as usize]
        })
        .collect()
}

fn runner_config(threads: usize) -> TraceRunnerConfig {
    TraceRunnerConfig {
        sensitivity: SensitivityConfig::for_profile(SensitivityProfile::Sensitive),
        candidates: CandidateSearchConfig {
            min_shared_hashes: 1,
            min_plasmid_containment: 0.0,
            min_metagenome_containment: 0.0,
            top_candidates: 8,
        },
        resources: ResourceOpenOptions {
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
        threads,
        max_alignments_per_candidate: 8,
    }
}

fn build_database(directory: &std::path::Path, sequence: &[u8]) -> std::path::PathBuf {
    let input = directory.join("candidate.fa");
    fs::write(
        &input,
        format!(
            ">candidate\n{}\n",
            String::from_utf8(sequence.to_vec()).unwrap()
        ),
    )
    .unwrap();
    let database = directory.join("candidates.jam");
    let sketch = Command::new(env!("CARGO_BIN_EXE_jam"))
        .arg("--silent")
        .arg("sketch")
        .arg(&input)
        .arg("--singleton")
        .arg("--output")
        .arg(&database)
        .arg("--kmer-size")
        .arg("21")
        .arg("--fscale")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        sketch.status.success(),
        "sketch failed: {}",
        String::from_utf8_lossy(&sketch.stderr)
    );
    database
}

#[test]
fn one_failed_candidate_does_not_hide_a_success_and_order_is_stable() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-failures-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(19, 2_000);
    let success = directory.path().join("success.fa");
    fs::write(
        &success,
        format!(
            ">success\n{}\n",
            String::from_utf8(sequence.clone()).unwrap()
        ),
    )
    .unwrap();
    let failing = directory.path().join("failing.fa");
    fs::write(
        &failing,
        format!(
            ">failing\n{}\n",
            String::from_utf8(sequence.clone()).unwrap()
        ),
    )
    .unwrap();

    let database = directory.path().join("candidates.jam");
    let sketch = Command::new(env!("CARGO_BIN_EXE_jam"))
        .arg("--silent")
        .arg("sketch")
        .arg(&success)
        .arg(&failing)
        .arg("--singleton")
        .arg("--output")
        .arg(&database)
        .arg("--kmer-size")
        .arg("21")
        .arg("--fscale")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        sketch.status.success(),
        "sketch failed: {}",
        String::from_utf8_lossy(&sketch.stderr)
    );

    let catalog = TraceCatalog::from_entries(vec![
        CatalogEntry {
            metagenome_id: "success".to_string(),
            jma: None,
            jma_index: None,
            raw: Some(success.to_string_lossy().into_owned()),
        },
        CatalogEntry {
            metagenome_id: "failing".to_string(),
            jma: None,
            jma_index: None,
            raw: Some(
                directory
                    .path()
                    .join("missing.fa")
                    .to_string_lossy()
                    .into_owned(),
            ),
        },
    ])
    .unwrap();
    let query = TraceQuery {
        plasmid_id: "trace-plasmid".to_string(),
        plasmid_sequence: sequence,
        database: database.to_string_lossy().into_owned(),
        catalog,
    };

    let single = TraceRunner::new(runner_config(1))
        .unwrap()
        .run(&query)
        .unwrap();
    let parallel = TraceRunner::new(runner_config(4))
        .unwrap()
        .run(&query)
        .unwrap();
    assert_eq!(single.metagenomes.len(), 2);
    assert_eq!(
        single
            .metagenomes
            .iter()
            .map(|result| result.metagenome_id.as_str())
            .collect::<Vec<_>>(),
        vec!["failing", "success"]
    );
    let failed = &single.metagenomes[0];
    assert_eq!(failed.status, jam_rs::trace::model::TraceStatus::Failed);
    assert_eq!(failed.failures[0].stage, "raw");
    assert_eq!(failed.failures[0].code, "resource_io_error");
    assert!(!failed.failures[0].retryable);
    assert_eq!(failed.resource_metrics.metadata_requests, 0);

    let succeeded = &single.metagenomes[1];
    assert_eq!(
        succeeded.status,
        jam_rs::trace::model::TraceStatus::Complete
    );
    assert!(!succeeded.alignments.is_empty());
    assert_eq!(single.counters.failures, 1);
    assert_eq!(parallel.counters.failures, 1);
    for (left, right) in single.metagenomes.iter().zip(&parallel.metagenomes) {
        assert_eq!(left.metagenome_id, right.metagenome_id);
        assert_eq!(left.status, right.status);
        assert_eq!(left.failures, right.failures);
        assert_eq!(left.alignments, right.alignments);
        assert_eq!(left.coverage, right.coverage);
    }
}

#[test]
fn remote_http_failures_keep_codes_retryability_and_metrics() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-http-failures-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(27, 2_000);
    let database = build_database(directory.path(), &sequence);
    for (status, code, retryable) in [
        (403, "access_denied", false),
        (404, "missing_object", false),
        (503, "retryable_server_error", true),
    ] {
        let fixture =
            remote::HttpFixture::new(b">not-a-valid-assembly".to_vec()).with_status(status);
        let catalog = TraceCatalog::from_entries(vec![CatalogEntry {
            metagenome_id: "candidate".to_string(),
            jma: None,
            jma_index: None,
            raw: Some(fixture.url("/assembly.fa?token=redacted")),
        }])
        .unwrap();
        let result = TraceRunner::new(runner_config(1))
            .unwrap()
            .run(&TraceQuery {
                plasmid_id: "trace-plasmid".to_string(),
                plasmid_sequence: sequence.clone(),
                database: database.to_string_lossy().into_owned(),
                catalog,
            })
            .unwrap();
        let result = &result.metagenomes[0];
        assert_eq!(result.status, jam_rs::trace::model::TraceStatus::Failed);
        assert_eq!(result.failures[0].code, code);
        assert_eq!(result.failures[0].retryable, retryable);
        assert!(result.resource_metrics.metadata_requests > 0);
        assert!(result.resource_metrics.head_requests > 0);
        assert!(
            !result.failures[0]
                .resource
                .as_ref()
                .unwrap()
                .contains("token")
        );
    }
}

#[test]
fn jma_failure_falls_back_to_raw_and_merges_available_metrics() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-jma-fallback-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(31, 2_000);
    let database = build_database(directory.path(), &sequence);
    let raw = directory.path().join("assembly.fa");
    fs::write(
        &raw,
        format!(
            ">candidate\n{}\n",
            String::from_utf8(sequence.clone()).unwrap()
        ),
    )
    .unwrap();
    let invalid_jma = directory.path().join("invalid.jma");
    fs::write(&invalid_jma, vec![0_u8; 160]).unwrap();
    let catalog = TraceCatalog::from_entries(vec![CatalogEntry {
        metagenome_id: "candidate".to_string(),
        jma: Some(invalid_jma.to_string_lossy().into_owned()),
        jma_index: None,
        raw: Some(raw.to_string_lossy().into_owned()),
    }])
    .unwrap();
    let result = TraceRunner::new(runner_config(1))
        .unwrap()
        .run(&TraceQuery {
            plasmid_id: "trace-plasmid".to_string(),
            plasmid_sequence: sequence,
            database: database.to_string_lossy().into_owned(),
            catalog,
        })
        .unwrap();
    let result = &result.metagenomes[0];
    assert_eq!(result.status, jam_rs::trace::model::TraceStatus::Partial);
    assert_eq!(result.failures[0].stage, "jma");
    assert_eq!(result.failures[0].code, "jma_invalid_magic");
    assert!(!result.alignments.is_empty());
    assert!(result.resource_metrics.stream_requests > 0);
}
