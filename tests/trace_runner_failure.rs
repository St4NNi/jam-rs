use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::trace::catalog::{CatalogEntry, TraceCatalog};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use jam_rs::trace::runner::{TraceQuery, TraceRunner, TraceRunnerConfig};
use jam_rs::trace::screen::CandidateSearchConfig;
use sha2::{Digest, Sha256};
use std::fs;
use std::path::Path;
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
        io_concurrency: threads.max(1),
        max_alignments_per_candidate: 8,
        query_kind: jam_rs::trace::model::QueryKind::Unknown,
        topology_requested: jam_rs::trace::model::TopologyRequested::Auto,
        topology_margin_bases: 100,
    }
}

fn build_database(directory: &Path, inputs: &[&Path]) -> std::path::PathBuf {
    let database = directory.join("candidates.jam");
    let mut command = Command::new(env!("CARGO_BIN_EXE_jam"));
    command.arg("--silent").arg("sketch");
    for input in inputs {
        command.arg(input);
    }
    let sketch = command
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

fn entry(metagenome_id: &str, resource_uri: String, bytes: &[u8]) -> CatalogEntry {
    CatalogEntry {
        metagenome_id: metagenome_id.to_string(),
        resource_uri,
        sha256: format!("{:x}", Sha256::digest(bytes)),
        etag: None,
        object_version: None,
        label: None,
        source_assembly_sha256: None,
    }
}

fn write_fasta(path: &Path, name: &str, sequence: &[u8]) {
    fs::write(
        path,
        format!(
            ">{name}\n{}\n",
            String::from_utf8(sequence.to_vec()).unwrap()
        ),
    )
    .unwrap();
}

#[test]
fn one_failed_jma_candidate_does_not_hide_a_success_and_order_is_stable() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-failures-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(19, 2_000);
    let success_input = directory.path().join("success.fa");
    let failing_input = directory.path().join("failing.fa");
    write_fasta(&success_input, "success-contig", &sequence);
    write_fasta(&failing_input, "failing-contig", &sequence);
    let database = build_database(directory.path(), &[&success_input, &failing_input]);
    let success_archive = directory.path().join("success.jma");
    write_archive_from_fasta(
        &success_input,
        &success_archive,
        ArchiveBuildConfig {
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();
    let success_bytes = fs::read(&success_archive).unwrap();
    let catalog = TraceCatalog::from_entries(vec![
        entry(
            "success-contig",
            success_archive.to_string_lossy().into_owned(),
            &success_bytes,
        ),
        entry(
            "failing-contig",
            directory
                .path()
                .join("missing.jma")
                .to_string_lossy()
                .into_owned(),
            &[0],
        ),
    ])
    .unwrap();
    let query = TraceQuery {
        plasmid_id: "trace-query".to_string(),
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
        vec!["failing-contig", "success-contig"]
    );
    assert_eq!(
        single.metagenomes[0].status,
        jam_rs::trace::model::TraceStatus::Failed
    );
    assert_eq!(single.metagenomes[0].failures[0].stage, "jma");
    assert_eq!(single.metagenomes[0].failures[0].code, "resource_io_error");
    assert!(single.metagenomes[0].alignments.is_empty());
    assert_eq!(
        single.metagenomes[1].status,
        jam_rs::trace::model::TraceStatus::Complete
    );
    assert!(!single.metagenomes[1].alignments.is_empty());
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
fn remote_jma_failures_keep_codes_retryability_and_metrics() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-http-failures-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(27, 2_000);
    let input = directory.path().join("candidate.fa");
    write_fasta(&input, "candidate-contig", &sequence);
    let database = build_database(directory.path(), &[&input]);
    for (status, code, retryable) in [
        (403, "access_denied", false),
        (404, "missing_object", false),
        (503, "retryable_server_error", true),
    ] {
        let fixture = remote::HttpFixture::new(Vec::new()).with_status(status);
        let catalog = TraceCatalog::from_entries(vec![entry(
            "candidate-contig",
            fixture.url("/candidate.jma?token=redacted"),
            &[],
        )])
        .unwrap();
        let output = TraceRunner::new(runner_config(1))
            .unwrap()
            .run(&TraceQuery {
                plasmid_id: "trace-query".to_string(),
                plasmid_sequence: sequence.clone(),
                database: database.to_string_lossy().into_owned(),
                catalog,
            })
            .unwrap();
        let result = &output.metagenomes[0];
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
fn unavailable_seed_scheme_fails_without_an_alternate_sequence_backend() {
    let directory = tempfile::Builder::new()
        .prefix("trace-runner-seed-mismatch-")
        .tempdir_in("target")
        .unwrap();
    let sequence = dna(37, 2_000);
    let input = directory.path().join("candidate.fa");
    write_fasta(&input, "candidate-contig", &sequence);
    let database = build_database(directory.path(), &[&input]);
    let archive = directory.path().join("candidate.jma");
    write_archive_from_fasta(
        &input,
        &archive,
        ArchiveBuildConfig {
            k31_scale: 200,
            k21_scale: Some(200),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();
    let bytes = fs::read(&archive).unwrap();
    let catalog = TraceCatalog::from_entries(vec![entry(
        "candidate-contig",
        archive.to_string_lossy().into_owned(),
        &bytes,
    )])
    .unwrap();
    let output = TraceRunner::new(runner_config(1))
        .unwrap()
        .run(&TraceQuery {
            plasmid_id: "trace-query".to_string(),
            plasmid_sequence: sequence,
            database: database.to_string_lossy().into_owned(),
            catalog,
        })
        .unwrap();
    let result = &output.metagenomes[0];
    assert_eq!(result.status, jam_rs::trace::model::TraceStatus::Failed);
    assert_eq!(result.failures.len(), 1);
    assert_eq!(result.failures[0].code, "seed_level_mismatch");
    assert!(result.alignments.is_empty());
    assert!(result.coverage.is_none());
}
