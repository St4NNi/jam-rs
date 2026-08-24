use jam_rs::jam_index::{
    JamIndexBuildConfig, JamIndexPartReader, MetagenomeSource, ScreenSelectionPolicy,
    append_jam_index, build_jam_index, load_manifest,
};
use jam_rs::reader::JamReader;
use std::fs;
use tempfile::Builder;

fn directory(prefix: &str) -> tempfile::TempDir {
    Builder::new().prefix(prefix).tempdir_in("target").unwrap()
}

fn source(root: &std::path::Path, id: &str, bases: usize) -> MetagenomeSource {
    let path = root.join(format!("{id}.fasta"));
    let sequence = (0..bases)
        .map(|index| b"ACGT"[(index * 13 + index / 7) & 3] as char)
        .collect::<String>();
    fs::write(&path, format!(">{id}-contig\n{sequence}\n")).unwrap();
    MetagenomeSource {
        metagenome_id: id.to_string(),
        sequence_path: path,
    }
}

fn config(max_bases: u64, max_signatures: u64) -> JamIndexBuildConfig {
    JamIndexBuildConfig {
        max_part_bases: max_bases,
        max_part_signatures: max_signatures,
        parallel_parts: 2,
        source_manifest_sha256: "a".repeat(64),
        ..JamIndexBuildConfig::default()
    }
}

#[test]
fn build_splits_by_bases_and_every_part_searches_independently() {
    let source_dir = directory("jam-index-build-source-");
    let output_parent = directory("jam-index-build-output-");
    let output = output_parent.path().join("index");
    let sources = vec![
        source(source_dir.path(), "mg-a", 300),
        source(source_dir.path(), "mg-b", 300),
        source(source_dir.path(), "mg-c", 300),
    ];
    let stats = build_jam_index(&output, &sources, &config(500, 10_000)).unwrap();
    assert_eq!(stats.new_parts, 3);
    assert_eq!(stats.total_metagenomes, 3);
    let manifest = load_manifest(&output).unwrap();
    assert_eq!(manifest.parts.len(), 3);
    for part in &manifest.parts {
        let root = output.join(&part.directory);
        let screen = JamReader::open(root.join(&part.screen_file)).unwrap();
        let data = JamIndexPartReader::open(root.join(&part.data_file)).unwrap();
        assert_eq!(screen.sample_names().len(), 1);
        assert_eq!(data.metagenomes().len(), 1);
        assert_eq!(data.read_contig(0).unwrap().len(), 300);
    }
}

#[test]
fn append_publishes_only_new_parts_and_preserves_old_bytes() {
    let source_dir = directory("jam-index-append-source-");
    let output_parent = directory("jam-index-append-output-");
    let output = output_parent.path().join("index");
    let first = source(source_dir.path(), "mg-a", 400);
    build_jam_index(
        &output,
        std::slice::from_ref(&first),
        &config(1_000, 10_000),
    )
    .unwrap();
    let before = load_manifest(&output).unwrap();
    let old_screen = fs::read(
        output
            .join(&before.parts[0].directory)
            .join(&before.parts[0].screen_file),
    )
    .unwrap();
    let second = source(source_dir.path(), "mg-b", 500);
    let mut append_config = config(1_000, 10_000);
    append_config.source_manifest_sha256 = "b".repeat(64);
    let stats = append_jam_index(&output, &[second], &append_config).unwrap();
    assert_eq!(stats.new_parts, 1);
    let after = load_manifest(&output).unwrap();
    assert_eq!(after.parts.len(), 2);
    assert_eq!(before.parts[0], after.parts[0]);
    assert_eq!(
        old_screen,
        fs::read(
            output
                .join(&after.parts[0].directory)
                .join(&after.parts[0].screen_file)
        )
        .unwrap()
    );
    assert!(append_jam_index(&output, &[first], &append_config).is_err());
}

#[test]
fn append_rejects_policy_change_and_signature_bound_splits() {
    let source_dir = directory("jam-index-policy-source-");
    let output_parent = directory("jam-index-policy-output-");
    let output = output_parent.path().join("index");
    let sources = vec![
        source(source_dir.path(), "mg-a", 10_000),
        source(source_dir.path(), "mg-b", 10_000),
    ];
    let signature_limited = config(1_000_000, 600);
    build_jam_index(&output, &sources, &signature_limited).unwrap();
    assert_eq!(load_manifest(&output).unwrap().parts.len(), 2);

    let mut changed = signature_limited;
    changed.selection_policy = ScreenSelectionPolicy::smaller_signatures();
    changed.source_manifest_sha256 = "c".repeat(64);
    let third = source(source_dir.path(), "mg-c", 100);
    assert!(append_jam_index(&output, &[third], &changed).is_err());
}
