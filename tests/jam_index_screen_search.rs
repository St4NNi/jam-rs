use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexScreenConfig, MetagenomeSource, ScreenSelectionPolicy,
    append_jam_index, build_jam_index, prepare_screen_query, search_jam_index,
};
use std::fs;
use tempfile::Builder;

fn directory(prefix: &str) -> tempfile::TempDir {
    Builder::new().prefix(prefix).tempdir_in("target").unwrap()
}

fn sequence(length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            b"ACGT"[(state as usize) & 3]
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn write_source(root: &std::path::Path, id: &str, contigs: &[Vec<u8>]) -> MetagenomeSource {
    let path = root.join(format!("{id}.fasta"));
    let mut fasta = String::new();
    for (index, sequence) in contigs.iter().enumerate() {
        fasta.push_str(&format!(
            ">{id}-{index}\n{}\n",
            String::from_utf8_lossy(sequence)
        ));
    }
    fs::write(&path, fasta).unwrap();
    MetagenomeSource {
        metagenome_id: id.to_string(),
        sequence_path: path,
    }
}

fn build_fixture() -> (tempfile::TempDir, tempfile::TempDir, Vec<u8>) {
    let source_root = directory("jam-index-search-source-");
    let output_root = directory("jam-index-search-output-");
    let query = sequence(160, 0x1234_5678_9abc_def0);
    let sources = vec![
        write_source(source_root.path(), "target", std::slice::from_ref(&query)),
        write_source(source_root.path(), "distractor-a", &[sequence(2_000, 17)]),
        write_source(source_root.path(), "distractor-b", &[sequence(2_000, 29)]),
    ];
    let output = output_root.path().join("index");
    build_jam_index(
        &output,
        &sources,
        &JamIndexBuildConfig {
            max_part_bases: 2_200,
            max_part_signatures: 10_000,
            parallel_parts: 2,
            source_manifest_sha256: "a".repeat(64),
            ..JamIndexBuildConfig::default()
        },
    )
    .unwrap();
    (source_root, output_root, query)
}

#[test]
fn exact_160_base_standalone_contig_is_selected_across_parallel_parts() {
    let (_sources, output_root, query) = build_fixture();
    let output = output_root.path().join("index");
    let prepared =
        prepare_screen_query(&query, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let result = search_jam_index(
        &output,
        &prepared,
        JamIndexScreenConfig {
            top_candidates: 2,
            accumulator_capacity: 8,
            parallel_parts: 2,
            ..JamIndexScreenConfig::default()
        },
    )
    .unwrap();
    assert_eq!(result.candidates[0].metagenome_id, "target");
    assert!(result.candidates[0].shared_hash_count >= 4);
    assert!(!result.candidates[0].shared_hashes.is_empty());
    assert_eq!(result.metrics.parts_searched, 2);
    assert!(result.metrics.maximum_accumulator_entries <= 8);
}

#[test]
fn reverse_query_and_append_only_part_produce_the_same_target() {
    let (source_root, output_root, query) = build_fixture();
    let output = output_root.path().join("index");
    let appended = write_source(
        source_root.path(),
        "appended-target",
        std::slice::from_ref(&query),
    );
    append_jam_index(
        &output,
        &[appended],
        &JamIndexBuildConfig {
            max_part_bases: 10_000,
            max_part_signatures: 10_000,
            parallel_parts: 1,
            source_manifest_sha256: "b".repeat(64),
            ..JamIndexBuildConfig::default()
        },
    )
    .unwrap();
    let reverse = reverse_complement(&query);
    let prepared =
        prepare_screen_query(&reverse, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let result = search_jam_index(
        &output,
        &prepared,
        JamIndexScreenConfig {
            top_candidates: 4,
            accumulator_capacity: 8,
            parallel_parts: 3,
            ..JamIndexScreenConfig::default()
        },
    )
    .unwrap();
    let names = result
        .candidates
        .iter()
        .map(|candidate| candidate.metagenome_id.as_str())
        .collect::<Vec<_>>();
    assert!(names.contains(&"target"));
    assert!(names.contains(&"appended-target"));
    assert!(
        result
            .candidates
            .iter()
            .all(|candidate| !candidate.shared_hashes.is_empty())
    );
}

#[test]
fn globally_unselected_candidates_do_not_retain_shared_hash_vectors() {
    let (_sources, output_root, query) = build_fixture();
    let output = output_root.path().join("index");
    let prepared =
        prepare_screen_query(&query, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let result = search_jam_index(
        &output,
        &prepared,
        JamIndexScreenConfig {
            top_candidates: 1,
            accumulator_capacity: 4,
            parallel_parts: 2,
            ..JamIndexScreenConfig::default()
        },
    )
    .unwrap();
    assert_eq!(result.candidates.len(), 1);
    assert_eq!(result.metrics.selected_candidates, 1);
    assert!(result.metrics.maximum_accumulator_entries <= 4);
}
