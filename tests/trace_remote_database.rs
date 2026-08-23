#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::trace::catalog::{CatalogEntry, TraceCatalog};
use jam_rs::trace::config::{SensitivityConfig, SensitivityProfile};
use jam_rs::trace::runner::{RunnerError, TraceQuery, TraceRunner, TraceRunnerConfig};
use jam_rs::trace::screen::CandidateSearchConfig;
use remote::HttpFixture;
use std::fs;
use std::process::Command;

fn deterministic_dna(seed: u64, length: usize) -> String {
    let mut state = seed;
    let bases = *b"ACGT";
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            bases[(state >> 62) as usize] as char
        })
        .collect()
}

#[test]
fn remote_candidate_database_is_identity_cached_before_mmap() {
    let directory = tempfile::Builder::new()
        .prefix("trace-remote-database-")
        .tempdir_in("target")
        .unwrap();
    let sequence = deterministic_dna(41, 2_000);
    let assembly = directory.path().join("assembly.fa");
    fs::write(&assembly, format!(">trace-contig\n{sequence}\n")).unwrap();
    let database = directory.path().join("metagenomes.jam");
    let output = Command::new(env!("CARGO_BIN_EXE_jam"))
        .arg("--silent")
        .arg("sketch")
        .arg(&assembly)
        .arg("--output")
        .arg(&database)
        .arg("--kmer-size")
        .arg("21")
        .arg("--fscale")
        .arg("1")
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "sketch failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let archive = directory.path().join("assembly.jma");
    write_archive_from_fasta(
        &assembly,
        &archive,
        ArchiveBuildConfig {
            block_bases: 512,
            k31_scale: 1,
            k21_scale: Some(1),
            ..ArchiveBuildConfig::default()
        },
    )
    .unwrap();
    let catalog = TraceCatalog::from_entries(vec![CatalogEntry {
        metagenome_id: "assembly.fa".to_string(),
        jma: Some(archive.to_string_lossy().into_owned()),
        raw: None,
    }])
    .unwrap();
    let fixture = HttpFixture::new(fs::read(&database).unwrap());
    let database_url = fixture.url("/metagenomes.jam?token=do-not-record");
    let cache_dir = directory.path().join("cache");
    let resources = ResourceOpenOptions {
        cache_dir: Some(cache_dir.clone()),
        cache_block_bytes: 1024,
        max_cache_bytes: 1024 * 1024,
        max_retries: 0,
        ..ResourceOpenOptions::default()
    };
    let runner = TraceRunner::new(TraceRunnerConfig {
        sensitivity: SensitivityConfig::for_profile(SensitivityProfile::Sensitive),
        candidates: CandidateSearchConfig {
            min_shared_hashes: 1,
            min_plasmid_containment: 0.0,
            min_metagenome_containment: 0.0,
            top_candidates: 1,
        },
        resources: resources.clone(),
        threads: 1,
        max_alignments_per_candidate: 8,
    })
    .unwrap();
    let query = TraceQuery {
        plasmid_id: "trace-plasmid".to_string(),
        plasmid_sequence: sequence.into_bytes(),
        database: database_url.clone(),
        catalog,
    };

    let first = runner.run(&query).unwrap();
    assert_eq!(first.search.candidates.len(), 1);
    assert!(!first.metagenomes[0].alignments.is_empty());
    assert!(first.database_metrics.remote_bytes > 0);
    assert_eq!(
        first.database_input.redacted_locator,
        fixture.url("/metagenomes.jam")
    );
    assert!(first.database_input.sha256.is_some());
    let get_requests = fixture
        .requests()
        .iter()
        .filter(|request| request.method == "GET")
        .count();

    let second = runner.run(&query).unwrap();
    assert_eq!(second.search.candidates.len(), 1);
    assert_eq!(second.database_metrics.remote_bytes, 0);
    assert_eq!(
        fixture
            .requests()
            .iter()
            .filter(|request| request.method == "GET")
            .count(),
        get_requests
    );
    assert!(
        cache_dir
            .read_dir()
            .unwrap()
            .filter_map(Result::ok)
            .any(|entry| entry
                .path()
                .extension()
                .is_some_and(|extension| extension == "jam"))
    );

    let no_cache = TraceRunner::new(TraceRunnerConfig {
        resources: ResourceOpenOptions {
            cache_dir: None,
            ..resources
        },
        ..runner.config().clone()
    })
    .unwrap();
    assert!(matches!(
        no_cache.run(&query),
        Err(RunnerError::DatabaseCache { .. })
    ));
}
