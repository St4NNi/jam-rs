//! Stage 3 with results computed ahead on several threads gives what one loop gives alone.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::dense_tests::whole_contigs;
use super::tests::{align_located, build, single_contig_entries};
use super::walks::{AHEAD, COMPUTED_AHEAD};
use crate::decide::DecisionRule;
use crate::test_support::{random_bases, rng, with_substitutions};
use rand::RngExt;
use std::sync::atomic::Ordering;

/// A query of `length` bases or a little more: short unique stretches between copies of two
/// repeat elements at 90% identity and tandem repeats at 93%. The contig is a copy of the query
/// at 97% identity with a short insertion or deletion every few hundred bases.
fn repeat_rich(seed: u64, length: usize) -> (Vec<u8>, Vec<u8>) {
    let mut random = rng(seed);
    let elements: Vec<Vec<u8>> = (0..2).map(|_| random_bases(&mut random, 300)).collect();
    let mut query = Vec::new();
    while query.len() < length {
        query.extend(random_bases(&mut random, 30));
        let element = &elements[random.random_range(0..elements.len())];
        query.extend(with_substitutions(&mut random, element, 0.9));
        let period = random.random_range(2..6);
        let motif = random_bases(&mut random, period);
        let tandem: Vec<u8> = motif.iter().copied().cycle().take(400).collect();
        query.extend(with_substitutions(&mut random, &tandem, 0.93));
    }
    let mut contig = Vec::new();
    for piece in with_substitutions(&mut random, &query, 0.97).chunks(500) {
        contig.extend_from_slice(piece);
        match random.random_range(0..3) {
            0 => contig.extend(random_bases(&mut random, 3)),
            1 => contig.truncate(contig.len() - 2),
            _ => {}
        }
    }
    (query, contig)
}

/// A homopolymer contig, as in the screen of SPIRE assemblies, and a query with two runs of C:
/// many kernel boxes, many with the same bases at other positions.
fn homopolymer(seed: u64) -> (Vec<u8>, Vec<u8>) {
    let mut query = random_bases(&mut rng(seed), 600);
    query[200..260].fill(b'C');
    query[400..430].fill(b'C');
    (query, vec![b'C'; 4_000])
}

#[test]
fn computed_ahead_same() {
    let rule = DecisionRule {
        min_query_coverage: "0.05".parse().unwrap(),
        min_identity: 0.8,
        min_fragment_columns: 40,
        ..DecisionRule::default()
    };
    for (query, contig) in [repeat_rich(101, 10_000), homopolymer(91)] {
        let directory = tempfile::tempdir().unwrap();
        let shard = build(
            directory.path(),
            &single_contig_entries(&[("repeats", &contig)]),
        );
        let located = whole_contigs(&shard, 0, query.len() as u64);
        for circular in [false, true] {
            let run = |mode: u8, threads: usize| {
                AHEAD.store(mode, Ordering::Relaxed);
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                pool.install(|| align_located(&shard, &query, circular, &located, &rule))
            };
            let alone = run(2, 1);
            let before = COMPUTED_AHEAD.load(Ordering::Relaxed);
            for threads in [1, 4] {
                let ahead = run(1, threads);
                assert_eq!(ahead, alone, "circular {circular}, {threads} threads");
            }
            AHEAD.store(0, Ordering::Relaxed);
            // The window has walks and alignments, and results were computed ahead for it.
            let work = alone.work;
            assert!(work.alignments > 20 && work.chains > 10, "{work:?}");
            assert!(COMPUTED_AHEAD.load(Ordering::Relaxed) > before);
        }
    }
}
