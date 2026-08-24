use jam_rs::archive::SeedSchemeDescriptor;
use jam_rs::jma::builder::{ArchiveBuildConfig, write_archive_from_fasta};
use jam_rs::jma::reader::JmaReader;
use jam_rs::resource::local::LocalResource;
use jam_rs::resource::{ResourceMetrics, ResourceOpenOptions};
use jam_rs::trace::catalog::{CatalogEntry, TraceCatalog};
use jam_rs::trace::config::{SeedSensitivity, SensitivityConfig, SensitivityProfile};
use jam_rs::trace::model::{QueryKind, TopologyRequested};
use jam_rs::trace::runner::{TraceQuery, TraceRunner, TraceRunnerConfig};
use jam_rs::trace::screen::CandidateSearchConfig;
use jam_rs::trace::seeds::extract_seed_level;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::fs;
use std::path::Path;
use std::process::Command;

fn deterministic_dna(mut state: u64, length: usize) -> Vec<u8> {
    let bases = *b"ACGT";
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            bases[(state >> 62) as usize]
        })
        .collect()
}

fn fasta_record(name: &str, sequence: &[u8]) -> String {
    format!(">{}\n{}\n", name, String::from_utf8_lossy(sequence))
}

fn build_database(assembly: &Path, database: &Path) {
    let output = Command::new(env!("CARGO_BIN_EXE_jam"))
        .arg("--silent")
        .arg("sketch")
        .arg(assembly)
        .arg("--output")
        .arg(database)
        .arg("--kmer-size")
        .arg("21")
        .arg("--fscale")
        .arg("1")
        .output()
        .expect("run jam sketch");
    assert!(
        output.status.success(),
        "sketch failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn build_archive(
    assembly: &Path,
    archive: &Path,
    block_bases: usize,
    k31_scale: u64,
    k21_scale: u64,
) -> (Vec<u8>, Vec<SeedSchemeDescriptor>) {
    write_archive_from_fasta(
        assembly,
        archive,
        ArchiveBuildConfig {
            sequence_policy: jam_rs::sequence::SequenceBlockPolicy::Fixed { block_bases },
            k31_scale,
            k21_scale: Some(k21_scale),
            ..ArchiveBuildConfig::default()
        },
    )
    .expect("build indexed metagenome archive");
    let archive_bytes = fs::read(archive).expect("read JMA archive");
    let reader = JmaReader::from_resource(
        LocalResource::from_path(archive, ResourceOpenOptions::default())
            .expect("open local JMA resource"),
    )
    .expect("open embedded JMA index");
    (archive_bytes, reader.seed_schemes().copied().collect())
}

fn runner(sensitivity: SensitivityConfig) -> TraceRunner {
    TraceRunner::new(TraceRunnerConfig {
        sensitivity,
        candidates: CandidateSearchConfig {
            min_shared_hashes: 1,
            min_plasmid_containment: 0.0,
            min_metagenome_containment: 0.0,
            top_candidates: 1,
        },
        resources: ResourceOpenOptions {
            // Local block reads retain bounded accounting without spawning a
            // process per byte as the HTTP transport does.
            cache_block_bytes: 64,
            max_cache_bytes: 16 * 1024 * 1024,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
        threads: 1,
        io_concurrency: 1,
        max_alignments_per_candidate: 64,
        query_kind: QueryKind::Phage,
        topology_requested: TopologyRequested::Linear,
        topology_margin_bases: 10,
        memory_budget_bytes: 1024 * 1024 * 1024,
    })
    .expect("construct trace runner")
}

fn local_catalog(archive: &Path, archive_bytes: &[u8], missing_root: &Path) -> TraceCatalog {
    TraceCatalog::from_entries(vec![
        CatalogEntry {
            metagenome_id: "target.fa".to_string(),
            resource_uri: archive.to_string_lossy().into_owned(),
            sha256: format!("{:x}", Sha256::digest(archive_bytes)),
            etag: None,
            object_version: None,
            label: None,
            source_assembly_sha256: None,
        },
        CatalogEntry {
            // This resource is intentionally absent. It must not be opened
            // because it is not present in the candidate .jam result.
            metagenome_id: "not-a-candidate.fa".to_string(),
            resource_uri: missing_root
                .join("never-opened.jma")
                .to_string_lossy()
                .into_owned(),
            sha256: "00".repeat(32),
            etag: None,
            object_version: None,
            label: None,
            source_assembly_sha256: None,
        },
    ])
    .expect("construct local resource catalog")
}

fn expected_seed_bucket_reads(
    schemes: &[SeedSchemeDescriptor],
    query: &[u8],
    levels: &[SeedSensitivity],
) -> u64 {
    let mut seen = HashSet::new();
    let mut expected = 0_u64;
    for level in levels {
        let extracted = extract_seed_level(query, *level).expect("extract query seeds");
        for seed in extracted.seeds {
            let key = (
                level.k,
                level.scale,
                seed.hash,
                seed.canonical_kmer,
                seed.position,
            );
            if seen.insert(key)
                && schemes.iter().any(|scheme| {
                    scheme.span == u16::from(level.k)
                        && u64::from(scheme.density_parameter) == level.scale
                })
            {
                expected = expected.saturating_add(1);
            }
        }
    }
    expected
}

fn rounded_bytes(length: usize, block_bytes: u64) -> u64 {
    let length = length as u64;
    length.div_ceil(block_bytes) * block_bytes
}

fn assert_selective_metrics(
    metrics: ResourceMetrics,
    archive_len: usize,
    cache_block_bytes: u64,
    expected_seed_reads: u64,
    expected_sequence_reads: u64,
) {
    assert!(
        metrics.range_requests > 0,
        "indexed path did not use ranges"
    );
    assert!(
        metrics.range_requests < 512,
        "selective-I/O fixture made {} local range requests",
        metrics.range_requests
    );
    assert_eq!(
        metrics.stream_requests, 0,
        "indexed path streamed an object"
    );
    assert_eq!(
        metrics.full_object_fallbacks, 0,
        "indexed path fell back to a full object"
    );
    assert!(metrics.seed_buckets_read > 0);
    assert!(metrics.seed_buckets_read <= expected_seed_reads);
    assert_eq!(metrics.sequence_blocks_read, expected_sequence_reads);

    let archive_cost = rounded_bytes(archive_len, cache_block_bytes);
    assert!(
        metrics.requested_bytes < archive_cost,
        "indexed archive requested {} bytes of {archive_cost}",
        metrics.requested_bytes
    );
}

#[test]
fn generic_query_uses_local_indexed_jma_and_skips_missing_noncandidate() {
    let directory = tempfile::Builder::new()
        .prefix("trace-mge-selective-local-")
        .tempdir_in("target")
        .expect("temporary selective-I/O fixture directory");
    let query = deterministic_dna(17, 2_048);
    let unrelated = deterministic_dna(91, 2_048);
    let assembly = directory.path().join("target.fa");
    fs::write(
        &assembly,
        format!(
            "{}{}",
            fasta_record("target-contig", &query),
            fasta_record("unrelated-contig", &unrelated)
        ),
    )
    .unwrap();
    let database = directory.path().join("metagenomes.jam");
    build_database(&assembly, &database);
    let archive_path = directory.path().join("target.jma");
    let (archive, index) = build_archive(&assembly, &archive_path, 4_096, 100, 200);
    let catalog = local_catalog(&archive_path, &archive, &directory.path().join("missing"));
    let sensitivity = SensitivityConfig::for_profile(SensitivityProfile::Sensitive);
    let expected_seed_reads = expected_seed_bucket_reads(
        &index,
        &query,
        &[sensitivity.primary, sensitivity.rescue.unwrap()],
    );
    let output = runner(sensitivity)
        .run(&TraceQuery {
            plasmid_id: "query-phage".to_string(),
            plasmid_sequence: query.clone(),
            database: database.to_string_lossy().into_owned(),
            catalog,
        })
        .expect("run generic query trace");
    assert_eq!(output.search.candidates.len(), 1);
    let result = &output.metagenomes[0];
    assert_eq!(result.query_kind, QueryKind::Phage);
    assert!(!result.alignments.is_empty());
    assert_eq!(result.alignments[0].contig_id, "target-contig");
    assert!(result.coverage.as_ref().unwrap().supported_bases > 0);
    assert_eq!(result.failures.len(), 0);
    assert_selective_metrics(
        result.resource_metrics,
        archive.len(),
        64,
        expected_seed_reads,
        1,
    );
}

#[test]
fn gap_rescue_uses_only_local_fragment_blocks() {
    let directory = tempfile::Builder::new()
        .prefix("trace-mge-gap-rescue-local-")
        .tempdir_in("target")
        .expect("temporary gap-rescue fixture directory");
    let query = deterministic_dna(321, 90);
    let prefix = query[..33].to_vec();
    let rescue_fragment = query[60..82].to_vec();
    let unrelated = deterministic_dna(987, 64);
    let assembly = directory.path().join("target.fa");
    fs::write(
        &assembly,
        format!(
            "{}{}{}",
            fasta_record("prefix-contig", &prefix),
            fasta_record("rescue-contig", &rescue_fragment),
            fasta_record("unrelated-contig", &unrelated)
        ),
    )
    .unwrap();
    let database = directory.path().join("metagenomes.jam");
    build_database(&assembly, &database);
    let archive_path = directory.path().join("target.jma");
    let (archive, index) = build_archive(&assembly, &archive_path, 128, 1, 1);
    let catalog = local_catalog(&archive_path, &archive, &directory.path().join("missing"));

    let mut sensitivity = SensitivityConfig::for_profile(SensitivityProfile::Sensitive);
    sensitivity.primary = SeedSensitivity {
        k: 31,
        scale: 1,
        max_occurrences: 256,
    };
    sensitivity.rescue = Some(SeedSensitivity {
        k: 21,
        scale: 1,
        max_occurrences: 256,
    });
    sensitivity.gap_rescue.dense_primary = None;
    sensitivity.gap_rescue.max_rounds = 2;
    sensitivity.gap_rescue.min_gap_bases = 20;
    sensitivity.gap_rescue.flank_bases = 16;
    let expected_seed_reads = expected_seed_bucket_reads(
        &index,
        &query,
        &[sensitivity.primary, sensitivity.rescue.unwrap()],
    );
    let output = runner(sensitivity)
        .run(&TraceQuery {
            plasmid_id: "query-linear-phage".to_string(),
            plasmid_sequence: query,
            database: database.to_string_lossy().into_owned(),
            catalog,
        })
        .expect("run gap-rescue trace");
    let result = &output.metagenomes[0];
    assert_eq!(result.query_kind, QueryKind::Phage);
    assert!(
        result
            .rescue_rounds
            .iter()
            .any(|round| round.round > 1 && round.new_query_bases_supported > 0),
        "the k=21 rescue round did not add the synthetic fragment"
    );
    assert!(
        result
            .alignments
            .iter()
            .all(|alignment| alignment.contig_id == "prefix-contig"
                || alignment.contig_id == "rescue-contig")
    );
    assert_selective_metrics(
        result.resource_metrics,
        archive.len(),
        64,
        expected_seed_reads,
        2,
    );
}
