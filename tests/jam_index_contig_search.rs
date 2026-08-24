use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexCandidate, JamIndexContigSearchConfig, MetagenomeSource,
    QueryHashOccurrence, ScreenSelectionPolicy, SharedScreenHash, build_jam_index,
    prepare_screen_query, search_jam_index, select_candidate_contigs,
};
use std::fs;
use tempfile::Builder;

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

fn fixture() -> (tempfile::TempDir, tempfile::TempDir, Vec<u8>) {
    let source_root = Builder::new()
        .prefix("jam-index-contig-source-")
        .tempdir_in("target")
        .unwrap();
    let output_root = Builder::new()
        .prefix("jam-index-contig-output-")
        .tempdir_in("target")
        .unwrap();
    let query = sequence(160, 0x1234_5678_9abc_def0);
    let target_path = source_root.path().join("target.fasta");
    fs::write(
        &target_path,
        format!(
            ">noise\n{}\n>exact-target\n{}\n>weak\n{}\n",
            String::from_utf8_lossy(&sequence(2_000, 3)),
            String::from_utf8_lossy(&query),
            String::from_utf8_lossy(&sequence(400, 7))
        ),
    )
    .unwrap();
    let distractor_path = source_root.path().join("distractor.fasta");
    fs::write(
        &distractor_path,
        format!(
            ">distractor\n{}\n",
            String::from_utf8_lossy(&sequence(2_000, 11))
        ),
    )
    .unwrap();
    let output = output_root.path().join("index");
    build_jam_index(
        &output,
        &[
            MetagenomeSource {
                metagenome_id: "target".to_string(),
                sequence_path: target_path,
            },
            MetagenomeSource {
                metagenome_id: "distractor".to_string(),
                sequence_path: distractor_path,
            },
        ],
        &JamIndexBuildConfig {
            max_part_bases: 3_000,
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
fn selected_hashes_map_to_the_exact_contig_without_positions() {
    let (_source, output, query) = fixture();
    let root = output.path().join("index");
    let prepared =
        prepare_screen_query(&query, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let screen = search_jam_index(&root, &prepared, Default::default()).unwrap();
    let contigs = select_candidate_contigs(
        &root,
        &prepared,
        &screen.candidates,
        JamIndexContigSearchConfig {
            initial_contigs_per_candidate: 1,
            max_contigs_per_candidate: 4,
            accumulator_capacity: 8,
            parallel_parts: 2,
            ..JamIndexContigSearchConfig::default()
        },
    )
    .unwrap();
    let target = contigs
        .plans
        .iter()
        .find(|plan| plan.metagenome_id == "target")
        .unwrap();
    assert_eq!(target.initial_contigs()[0].contig_name, "exact-target");
    assert!(target.sequential_fallback_range.is_none());
    assert!(target.ranked_contigs.len() <= 4);
    assert!(contigs.metrics.maximum_accumulator_entries <= 8);
}

#[test]
fn weaker_contigs_are_retained_only_as_a_bounded_expansion_suffix() {
    let (_source, output, query) = fixture();
    let root = output.path().join("index");
    let prepared =
        prepare_screen_query(&query, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let screen = search_jam_index(&root, &prepared, Default::default()).unwrap();
    let result = select_candidate_contigs(
        &root,
        &prepared,
        &screen.candidates,
        JamIndexContigSearchConfig {
            initial_contigs_per_candidate: 1,
            max_contigs_per_candidate: 2,
            accumulator_capacity: 4,
            parallel_parts: 1,
            ..JamIndexContigSearchConfig::default()
        },
    )
    .unwrap();
    assert!(
        result
            .plans
            .iter()
            .all(|plan| plan.initial_contigs().len() <= 1)
    );
    assert!(
        result
            .plans
            .iter()
            .all(|plan| plan.weaker_contigs().len() <= 1)
    );
}

#[test]
fn only_a_strong_no_hit_candidate_receives_a_sequential_range() {
    let (_source, output, query) = fixture();
    let root = output.path().join("index");
    let prepared =
        prepare_screen_query(&query, &ScreenSelectionPolicy::default_signatures()).unwrap();
    let missing = SharedScreenHash {
        hash: 1,
        document_frequency: 1,
        occurrences: vec![QueryHashOccurrence {
            packed_kmer: 1,
            query_position: 0,
            query_reverse: false,
            query_window_id: 0,
        }],
    };
    let candidate = |shared_hash_count| JamIndexCandidate {
        part_id: 0,
        metagenome_local_id: 0,
        metagenome_id: "target".to_string(),
        rank: 1,
        shared_hash_count,
        supported_query_windows: 1,
        longest_supported_window_run: 1,
        weighted_hash_sum: 1.0,
        shared_hashes: vec![missing.clone()],
    };
    let weak = select_candidate_contigs(
        &root,
        &prepared,
        &[candidate(1)],
        JamIndexContigSearchConfig {
            strong_candidate_shared_hashes: 4,
            ..JamIndexContigSearchConfig::default()
        },
    )
    .unwrap();
    assert!(weak.plans[0].sequential_fallback_range.is_none());
    let strong = select_candidate_contigs(
        &root,
        &prepared,
        &[candidate(4)],
        JamIndexContigSearchConfig {
            strong_candidate_shared_hashes: 4,
            ..JamIndexContigSearchConfig::default()
        },
    )
    .unwrap();
    assert!(strong.plans[0].sequential_fallback_range.is_some());
}
