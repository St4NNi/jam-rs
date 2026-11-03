//! Tests of the reach of prepared contigs and of entries dropped once their fragments cannot
//! reach the coverage.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::tests::{EntryContigs, build, located};
use super::*;
use crate::align::Interval;
use crate::test_support::{no_stop, random_bases, rng, with_substitutions};

/// Contigs that hold pieces of `query`: copies, substituted copies on both strands, a run of
/// one base and a piece inside unrelated bases, each between random flanks. A piece starts
/// every `step` query bases; flanks are `flank` bases or up to 1,800 more.
fn piece_contigs(query: &[u8], seed: u64, step: usize, flank: usize) -> Vec<(String, Vec<u8>)> {
    let mut random = rng(seed);
    let mut contigs = Vec::new();
    for (index, start) in (0..query.len() - 1_500).step_by(step).enumerate() {
        let length = 400 + 300 * (index % 4);
        let piece = &query[start..start + length];
        let piece = match index % 4 {
            0 => piece.to_vec(),
            1 => with_substitutions(&mut random, piece, 0.9),
            2 => crate::seeds::reverse_complement(&with_substitutions(&mut random, piece, 0.85)),
            _ => [piece, &[b'A'; 120][..], piece].concat(),
        };
        let flank = flank + 900 * (index % 3);
        let bases = [
            random_bases(&mut random, flank),
            piece,
            random_bases(&mut random, flank),
        ]
        .concat();
        contigs.push((format!("c{index}"), bases));
    }
    contigs
}

/// The query bases covered by `intervals`.
fn covered(intervals: Vec<Interval>) -> u64 {
    Interval::total_length(&Interval::union(intervals))
}

#[test]
fn reach_holds_fragments() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(5);
    let query = [
        random_bases(&mut random, 9_000),
        vec![b'A'; 150],
        random_bases(&mut random, 3_000),
    ]
    .concat();
    let entries: Vec<EntryContigs> = (0..3)
        .map(|seed| {
            (
                format!("e{seed}"),
                piece_contigs(&query, 40 + seed, 1_700, 200),
            )
        })
        .collect();
    let shard = build(directory.path(), &entries);
    let config = AlignStageConfig::default();
    let rule = DecisionRule::default();
    let prepared = PreparedQuery::new(&query, false, &config).unwrap();
    let workspace = AlignStageWorkspace::default();
    let mut checked = 0;
    for entry in 0..entries.len() as u32 {
        // Small flanks leave cut window edges in the longer contigs.
        for flank in [100, 3_000] {
            let located = located(&shard, entry, &query, flank);
            let packed = shard.packed(entry).unwrap();
            for windows in &located.alignment_windows {
                let job = contig::ContigJob::new(
                    &packed,
                    windows.contig,
                    &prepared,
                    &config,
                    &rule,
                    no_stop(),
                )
                .unwrap();
                let mut work = AlignWork::default();
                let round = job.prepare(&windows.ranges, &workspace, &mut work).unwrap();
                let Some(reach) = job.reach(round.prepared.as_deref().unwrap()) else {
                    continue;
                };
                let mut fragments = place::Collected::new(config.max_fragment_bytes);
                job.finish(round, &workspace, &mut fragments, &mut work)
                    .unwrap();
                let segments: Vec<Interval> = (fragments.into_list().unwrap().iter())
                    .flat_map(|fragment| fragment.query_segments.clone())
                    .collect();
                let reached = covered(reach.clone());
                assert_eq!(covered([reach, segments.clone()].concat()), reached);
                checked += usize::from(!segments.is_empty());
            }
        }
    }
    assert!(checked >= 10, "{checked}");
}

#[test]
fn short_entries_stop() {
    let directory = tempfile::tempdir().unwrap();
    let mut random = rng(9);
    let query = random_bases(&mut random, 80_000);
    let whole = with_substitutions(&mut random, &query, 0.97);
    // Short flanks keep the reach of each piece small: past the first batch it falls short.
    let mut pieces = piece_contigs(&query, 3, 2_500, 100);
    pieces
        .iter_mut()
        .for_each(|(_, bases)| bases.truncate(bases.len().min(1_900)));
    let entries: Vec<EntryContigs> = vec![
        ("pieces".to_string(), pieces),
        ("whole".to_string(), vec![("c".to_string(), whole)]),
    ];
    let shard = build(directory.path(), &entries);
    let config = AlignStageConfig::default();
    let prepared = PreparedQuery::new(&query, false, &config).unwrap();
    let workspace = AlignStageWorkspace::default();
    for rule in [
        DecisionRule::default(),
        DecisionRule {
            min_identity: 0.8,
            min_fragment_identity: 0.8,
            ..DecisionRule::default()
        },
    ] {
        for entry in 0..2 {
            let located = located(&shard, entry, &query, 3_000);
            let packed = shard.packed(entry).unwrap();
            let aligned = align_packed_entry(
                &packed,
                &prepared,
                &located,
                &config,
                &rule,
                &workspace,
                no_stop(),
            )
            .unwrap();
            let (fragments, all_work) = entry_fragments(
                &packed,
                &prepared,
                &located,
                &config,
                &rule,
                &workspace,
                no_stop(),
            )
            .unwrap();
            let low = &prepared.low_complexity;
            let full = rule
                .decide(prepared.length(), fragments, low, || false, |_| 0)
                .unwrap();
            // The pieces cover less than half the query: aligning stops before the last ones.
            assert_eq!(aligned.decision.qualifies, entry == 1);
            assert_eq!(aligned.decision.qualifies, full.qualifies);
            if entry == 1 {
                assert_eq!(aligned.decision, full);
                assert_eq!(aligned.work, all_work);
            } else {
                assert!(aligned.work.alignments < all_work.alignments, "{aligned:?}");
            }
        }
    }
}
