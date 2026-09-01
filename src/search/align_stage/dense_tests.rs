//! End-to-end tests of stage 3 on repeats, low-complexity runs and short conserved blocks.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::tests::{align, align_located, build, fragments, rule, single_contig_entries};
use super::*;
use crate::align::{Interval, Strand, format_cigar};
use crate::index::Shard;
use crate::search::ContigWindows;
use crate::seeds::reverse_complement;
use crate::test_support::{random_bases, rng};

/// Every contig of an entry as one window.
fn whole_contigs(shard: &Shard, entry: u32, query_length: u64) -> LocatedEntry {
    let packed = shard.packed(entry).unwrap();
    let windows = (0..packed.contig_count())
        .map(|contig| {
            let range = packed.contig_range(contig).unwrap();
            let whole = 0..range.end - range.start;
            ContigWindows {
                contig: contig as u32,
                ranges: vec![whole],
            }
        })
        .collect();
    LocatedEntry {
        entry,
        coverage_bound: query_length,
        containment_estimate: 1.0,
        windows,
    }
}

fn cigars(aligned: &AlignedEntry) -> Vec<String> {
    fragments(aligned)
        .iter()
        .map(|fragment| format_cigar(&fragment.alignment.runs))
        .collect()
}

#[test]
fn costly_gap_split() {
    // The entry holds query[..700], 1,000 unrelated bases and query[700..]. Bridging the gap
    // costs more than the second block earns, so the first alignment is the first block alone;
    // the second block is aligned from the anchors it leaves out.
    let directory = tempfile::tempdir().unwrap();
    let query = random_bases(&mut rng(71), 1200);
    let contig = [
        &query[..700],
        &random_bases(&mut rng(72), 1000)[..],
        &query[700..],
    ]
    .concat();
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("split", &contig)]),
    );
    let aligned = align(&shard, 0, &query, false, &rule("1"));
    assert_eq!(aligned.decision.mosaic.covered_bases, 1200);
    let mut spans: Vec<(u64, u64)> = fragments(&aligned)
        .iter()
        .map(|f| (f.query_segments[0].start, f.query_segments[0].end))
        .collect();
    spans.sort();
    // The second block may take a base or two of the gap that match by chance.
    assert_eq!(spans.len(), 2, "{spans:?}");
    assert_eq!(spans[0], (0, 700));
    assert!(
        (697..=700).contains(&spans[1].0) && spans[1].1 == 1200,
        "{spans:?}"
    );
}

#[test]
fn ambiguity_codes_align() {
    let directory = tempfile::tempdir().unwrap();
    let query = [
        random_bases(&mut rng(51), 1000),
        vec![b'R'; 1000],
        random_bases(&mut rng(52), 1000),
    ]
    .concat();
    let reverse = reverse_complement(&query);
    let entries = single_contig_entries(&[("forward", &query), ("reverse", &reverse)]);
    let shard = build(directory.path(), &entries);
    for circular in [false, true] {
        for entry in [0, 1] {
            let lower = query.to_ascii_lowercase();
            let aligned = align(&shard, entry, &lower, circular, &rule("1"));
            assert!(aligned.decision.qualifies, "circular {circular}");
            assert_eq!(aligned.decision.mosaic.covered_bases, 3000);
            assert_eq!(aligned.decision.counts.matches, 3000);
        }
    }
}

#[test]
fn blocks_join_gap() {
    let directory = tempfile::tempdir().unwrap();
    let query = [vec![b'A'; 600], vec![b'C'; 600]].concat();
    let contig = [vec![b'A'; 600], vec![b'G'; 300], vec![b'C'; 600]].concat();
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("gap", &contig)]),
    );
    let aligned = align_located(
        &shard,
        &query,
        false,
        &whole_contigs(&shard, 0, 1200),
        &rule("1"),
    );
    assert!(fragments(&aligned).iter().any(|f| {
        format_cigar(&f.alignment.runs) == "600=300D600="
            && f.alignment.contig_interval == Interval::new(0, 1500).unwrap()
    }));
}

#[test]
fn conserved_blocks_cigar() {
    let directory = tempfile::tempdir().unwrap();
    let query = [vec![b'A'; 200], vec![b'C'; 200]].concat();
    let contig = [vec![b'A'; 200], vec![b'G'; 20], vec![b'C'; 200]].concat();
    let reverse = reverse_complement(&contig);
    let entries = single_contig_entries(&[("forward", &contig), ("reverse", &reverse)]);
    let shard = build(directory.path(), &entries);
    for circular in [false, true] {
        for entry in [0, 1] {
            let located = whole_contigs(&shard, entry, 400);
            let aligned = align_located(&shard, &query, circular, &located, &rule("1"));
            assert!(
                fragments(&aligned).iter().any(|f| {
                    format_cigar(&f.alignment.runs) == "200=20D200="
                        && f.alignment.contig_interval == Interval::new(0, 420).unwrap()
                        && f.query_segments.iter().map(|part| part.len()).sum::<u64>() == 400
                }),
                "entry {entry} circular {circular}: {:?}",
                cigars(&aligned)
            );
        }
    }
}

#[test]
fn periodic_mappings_kept() {
    // 4 kb of ACGT against itself: every diagonal a multiple of four apart is an exact
    // mapping, on both strands.
    let directory = tempfile::tempdir().unwrap();
    let query = b"ACGT".repeat(1000);
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("periodic", &query)]),
    );
    let aligned = align_located(
        &shard,
        &query,
        false,
        &whole_contigs(&shard, 0, 4_000),
        &rule("0.5"),
    );
    // Chains need 40 bases (the least chain score); fragments need 100 columns.
    let per_strand = |minimum: u64| 2 * (4_000 - minimum) / 4 + 1;
    assert_eq!(aligned.work.chains, 2 * per_strand(40));
    assert_eq!(aligned.work.redundant_chains, 0);
    let decision = &aligned.decision;
    let mappings = decision.mosaic.primary.len() + decision.mosaic.alternatives.len();
    assert_eq!(mappings as u64, 2 * per_strand(100));
    assert_eq!(
        decision.dropped_fragments as u64,
        2 * (per_strand(40) - per_strand(100))
    );
    assert!(fragments(&aligned).iter().all(|f| {
        let counts = f.alignment.counts();
        format_cigar(&f.alignment.runs) == format!("{}=", counts.matches)
    }));
}

#[test]
fn circular_repeat_locations() {
    for query in [vec![b'A'; 1200], b"ACGT".repeat(300)] {
        let directory = tempfile::tempdir().unwrap();
        let contig = [query.clone(), query.clone()].concat();
        let shard = build(
            directory.path(),
            &single_contig_entries(&[("copies", &contig)]),
        );
        let located = whole_contigs(&shard, 0, 1200);
        let aligned = align_located(&shard, &query, true, &located, &rule("1"));
        let locations: std::collections::BTreeSet<_> = fragments(&aligned)
            .iter()
            .filter(|f| {
                f.alignment.strand == Strand::Forward && f.alignment.counts().matches == 1200
            })
            .map(|f| {
                (
                    f.alignment.contig_interval.start,
                    f.alignment.contig_interval.end,
                )
            })
            .collect();
        let expected = (0..=1200).map(|start| (start, start + 1200)).collect();
        assert_eq!(locations, expected);
    }
}

#[test]
fn imperfect_repeat_copies() {
    // 4 kb of ACGT with one N against the same contig: the anchors would exceed the chaining
    // budget, so the copies are walked along their diagonals instead.
    let directory = tempfile::tempdir().unwrap();
    let mut query = b"ACGT".repeat(1000);
    query[2_000] = b'N';
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("dense", &query)]),
    );
    for circular in [false, true] {
        let located = whole_contigs(&shard, 0, 4_000);
        let aligned = align_located(&shard, &query, circular, &located, &rule("0.5"));
        assert_eq!(aligned.work.anchors, 0);
        assert_eq!(aligned.decision.mosaic.covered_bases, 4_000);
        assert!(aligned.decision.identity() > 0.999);
        let found = fragments(&aligned);
        assert!(
            found
                .iter()
                .all(|f| f.alignment.counts().substitutions <= 2)
        );
        // Every diagonal a multiple of four apart, on both strands: a linear query has one per
        // offset with at least 100 columns, a circular one each offset modulo its length.
        let expected = if circular {
            2 * 4_000 / 4
        } else {
            2 * (2 * (4_000 - 100) / 4 + 1)
        };
        assert!(
            found.len() <= expected && found.len() + 8 >= expected,
            "circular {circular}: {} of {expected}",
            found.len()
        );
    }
}

#[test]
fn low_complexity_aligned() {
    let directory = tempfile::tempdir().unwrap();
    // A homopolymer query and entry: no ordinary anchor exists.
    let poly = vec![b'A'; 400];
    // A mostly low-complexity query whose entry holds only that part, with clustered
    // substitutions in one block, inside unrelated sequence. The repeat stays short, as every
    // walk of it is aligned.
    let repeat: Vec<u8> = b"AAAT".repeat(60);
    let query = [repeat.clone(), random_bases(&mut rng(71), 50)].concat();
    let mut copy = repeat.clone();
    for at in [81, 121, 161] {
        copy[at] = if copy[at] == b'A' { b'C' } else { b'G' };
    }
    let source = [
        random_bases(&mut rng(72), 1000),
        copy,
        random_bases(&mut rng(73), 1000),
    ]
    .concat();
    let entries = single_contig_entries(&[("poly", &poly), ("repeat", &source)]);
    let shard = build(directory.path(), &entries);
    for circular in [false, true] {
        let found = align(&shard, 0, &poly, circular, &rule("0.5"));
        assert_eq!(found.decision.mosaic.covered_bases, 400);
        let found = align(&shard, 1, &query, circular, &rule("0.5"));
        let covered = found.decision.mosaic.covered_bases;
        assert!(covered >= 237, "{covered}");
        // A circular query may also align around its origin, across its unrelated tail.
        assert!(found.decision.qualifies);
        // Three substitutions in 240 bases leave at most 0.9875.
        assert!(found.decision.identity() >= if circular { 0.8 } else { 0.98 });
    }
}

#[test]
fn clustered_change_kept() {
    let directory = tempfile::tempdir().unwrap();
    let left = [random_bases(&mut rng(42), 998), b"GC".to_vec()].concat();
    let right = [b"GC".to_vec(), random_bases(&mut rng(43), 998)].concat();
    let original = b"AAAT".repeat(45);
    let mut changed = original.clone();
    changed[90] = b'C';
    let query = [left.clone(), changed, right.clone()].concat();
    let entries = vec![(
        "mosaic".to_string(),
        vec![
            ("left".to_string(), left),
            ("low".to_string(), original),
            ("right".to_string(), right),
        ],
    )];
    let shard = build(directory.path(), &entries);
    let any_fragment = DecisionRule {
        min_fragment_identity: 0.0,
        ..rule("1")
    };
    for circular in [false, true] {
        let located = whole_contigs(&shard, 0, query.len() as u64);
        let found = align_located(&shard, &query, circular, &located, &any_fragment);
        assert_eq!(found.decision.mosaic.covered_bases, 2180);
        assert!(fragments(&found).iter().any(|f| f.contig_index == 1
            && f.alignment.contig_interval == Interval::new(0, 180).unwrap()
            && f.alignment.counts().substitutions == 1));
    }
}

/// Searches random flanks around `middle_query` in entries "forward" and "reverse" whose three
/// contigs are the flanks and `middle_contig`, on the given and the opposite strand, and checks
/// that both give the whole query with the full gapped middle fragment.
fn gapped_middle_found(
    middle_query: &[u8],
    middle_contig: &[u8],
    middle: impl Fn(&crate::align::Alignment) -> bool,
) {
    let directory = tempfile::tempdir().unwrap();
    let (left, right) = (
        random_bases(&mut rng(61), 300),
        random_bases(&mut rng(62), 300),
    );
    let query = [left.as_slice(), middle_query, &right].concat();
    let contigs = |sequences: [Vec<u8>; 3]| {
        let names = ["c0", "c1", "c2"].map(String::from);
        names.into_iter().zip(sequences).collect::<Vec<_>>()
    };
    let entries = vec![
        (
            "forward".to_string(),
            contigs([left.clone(), middle_contig.to_vec(), right.clone()]),
        ),
        (
            "reverse".to_string(),
            contigs([
                reverse_complement(&right),
                reverse_complement(middle_contig),
                reverse_complement(&left),
            ]),
        ),
    ];
    let shard = build(directory.path(), &entries);
    for circular in [false, true] {
        for entry in [0, 1] {
            let context = format!("entry {entry} circular {circular}");
            let located = whole_contigs(&shard, entry, query.len() as u64);
            let found = align_located(&shard, &query, circular, &located, &rule("1"));
            let covered = found.decision.mosaic.covered_bases;
            assert_eq!(covered, query.len() as u64, "{context}");
            let contig = Interval::new(0, middle_contig.len() as u64).unwrap();
            assert!(
                fragments(&found).iter().any(|f| f.contig_index == 1
                    && f.alignment.contig_interval == contig
                    && middle(&f.alignment)),
                "{context}: {:?}",
                cigars(&found)
            );
        }
    }
}

#[test]
fn short_blocks_join() {
    // Conserved blocks of 30 bases, each shorter than an anchor chain, with 20 contig-only bases
    // between them.
    gapped_middle_found(
        &[[b'A'; 30], [b'C'; 30], [b'G'; 30]].concat(),
        &[
            &[b'A'; 30][..],
            &[b'T'; 20],
            &[b'C'; 30],
            &[b'T'; 20],
            &[b'G'; 30],
        ]
        .concat(),
        |alignment| format_cigar(&alignment.runs) == "30=20D30=20D30=",
    );
}

#[test]
fn substitution_beside_gap() {
    // The base after the first block differs, so no exact seed sits next to its end: 80
    // matches, one substitution and 39 contig-only bases.
    gapped_middle_found(
        &[&[b'A'; 40][..], b"T", &[b'C'; 40]].concat(),
        &[[b'A'; 40], [b'G'; 40], [b'C'; 40]].concat(),
        |alignment| {
            let counts = alignment.counts();
            let found = (
                counts.matches,
                counts.substitutions,
                counts.insertions,
                counts.deletions,
            );
            found == (80, 1, 0, 39)
        },
    );
}

#[test]
fn weak_join_columns() {
    // Blocks of 25 bases around 50 contig-only bases: one block has fewer than 100 columns,
    // the join has 100.
    gapped_middle_found(
        &[[b'A'; 25], [b'C'; 25]].concat(),
        &[&[b'A'; 25][..], &[b'G'; 50], &[b'C'; 25]].concat(),
        |alignment| format_cigar(&alignment.runs) == "25=50D25=",
    );
}

#[test]
fn weak_join_global() {
    // Three query bases against 50 contig bases between the blocks: three substitutions and 47
    // contig-only bases score above three query-only and 50 contig-only bases.
    gapped_middle_found(
        &[&[b'A'; 25][..], &[b'T'; 3], &[b'C'; 25]].concat(),
        &[&[b'A'; 25][..], &[b'G'; 50], &[b'C'; 25]].concat(),
        |alignment| {
            let counts = alignment.counts();
            let found = (
                counts.matches,
                counts.substitutions,
                counts.insertions,
                counts.deletions,
                alignment.score,
            );
            found == (50, 3, 0, 47, 39)
        },
    );
}

#[test]
fn short_joins_count() {
    // Blocks of 17 and 18 bases: 35 matches score 70; 65 contig-only bases score -70 and 100
    // score -105.
    gapped_middle_found(
        &[&[b'A'; 17][..], &[b'C'; 18]].concat(),
        &[&[b'A'; 17][..], &[b'G'; 65], &[b'C'; 18]].concat(),
        |alignment| format_cigar(&alignment.runs) == "17=65D18=" && alignment.score == 0,
    );
    gapped_middle_found(
        &[&[b'A'; 17][..], &[b'C'; 18]].concat(),
        &[&[b'A'; 17][..], &[b'G'; 100], &[b'C'; 18]].concat(),
        |alignment| format_cigar(&alignment.runs) == "17=100D18=" && alignment.score == -35,
    );
}
