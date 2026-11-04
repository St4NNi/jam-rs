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
/// repeat elements at 90% identity and short tandem repeats at 93%. The contig is a copy of the
/// query at 97% identity with a short insertion or deletion every few thousand bases.
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
        let tandem: Vec<u8> = motif.iter().copied().cycle().take(150).collect();
        query.extend(with_substitutions(&mut random, &tandem, 0.93));
    }
    let mut contig = Vec::new();
    for piece in with_substitutions(&mut random, &query, 0.97).chunks(3_000) {
        contig.extend_from_slice(piece);
        match random.random_range(0..3) {
            0 => contig.extend(random_bases(&mut random, 3)),
            1 => contig.truncate(contig.len() - 2),
            _ => {}
        }
    }
    (query, contig)
}

#[test]
fn computed_ahead_same() {
    let (query, contig) = repeat_rich(101, 10_000);
    let directory = tempfile::tempdir().unwrap();
    let shard = build(
        directory.path(),
        &single_contig_entries(&[("repeats", &contig)]),
    );
    let located = whole_contigs(&shard, 0, query.len() as u64);
    let rule = DecisionRule {
        min_query_coverage: "0.5".parse().unwrap(),
        min_identity: 0.8,
        min_fragment_columns: 40,
        ..DecisionRule::default()
    };
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
            assert_eq!(
                run(1, threads),
                alone,
                "circular {circular}, {threads} threads"
            );
        }
        AHEAD.store(0, Ordering::Relaxed);
        // The window has walks, joins and gapped alignments, and their results were computed
        // ahead.
        let work = alone.work;
        eprintln!("{work:?} {} fragments", alone.decision.mosaic.primary.len());
        assert!(work.alignments > 100 && work.chains > 10, "{work:?}");
        assert!(COMPUTED_AHEAD.load(Ordering::Relaxed) > before);
    }
}
