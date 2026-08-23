#[path = "trace_fixture_support_mge/mod.rs"]
mod fixtures;

use sha2::{Digest, Sha256};
use std::collections::BTreeSet;

#[test]
fn construction_truth_covers_required_query_classes_and_topologies() {
    let manifest = fixtures::load_manifest();
    let kinds = manifest
        .queries
        .iter()
        .map(|query| query.kind.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        kinds,
        BTreeSet::from(["other", "phage", "plasmid", "unknown"])
    );
    let topologies = manifest
        .queries
        .iter()
        .map(|query| query.topology_requested.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        topologies,
        BTreeSet::from(["auto", "circular", "linear", "unknown"])
    );
    for id in [
        "circular_origin_crossing",
        "circular_rotated_reference",
        "linear_exact",
        "phage_direct_repeat_one_copy",
        "phage_inverted_repeat_core",
        "synthetic_transposon_only",
        "synthetic_resistance_only",
        "phage_integrated_fragment",
        "circular_duplicate",
        "circular_unrelated",
    ] {
        let _ = fixtures::case(&manifest, id);
    }
}

#[test]
fn truth_intervals_are_bounded_and_do_not_claim_topology_from_unknown_inputs() {
    let manifest = fixtures::load_manifest();
    for case in &manifest.cases {
        let query = fixtures::query(&manifest, &case.query);
        assert!(case.expected.unique_query_coverage_bases <= query.length);
        for interval in &case.query_intervals {
            assert!(interval.start <= interval.end);
            assert!(interval.end <= query.length);
        }
        if case.topology_requested == "linear" {
            assert_ne!(case.expected.coordinate_model, "wrap");
        }
        if case.topology_requested == "unknown" {
            assert_eq!(case.expected.topology_evidence, "undetermined");
        }
    }
}

#[test]
fn checked_files_match_the_committed_fixture_manifest() {
    let checksums = fixtures::read_fixture("CHECKSUMS.sha256");
    for line in checksums
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
    {
        let (expected, file) = line.split_once("  ").expect("checksum line");
        let actual = format!("{:x}", Sha256::digest(fixtures::read_bytes(file)));
        assert_eq!(actual, expected, "fixture checksum changed for {file}");
    }
}

#[test]
fn fasta_edge_cases_remain_explicit() {
    assert!(
        fixtures::read_fasta("assembly_empty.fa")
            .unwrap()
            .is_empty()
    );
    assert!(fixtures::read_fasta("assembly_malformed.fa").is_err());
}
