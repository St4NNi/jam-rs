//! Tests of walks: periodic and circular queries, budgets, gapped alignments and joins.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::align::{
    AlignmentConfig, AlignmentWorkspace, AnchorConfig, BandConfig, ChainConfig, DenseResults,
    DenseWalk, JoinContext, JoinLimits, JoinScratch, QueryTable, Strand, Walk, align_dense_walk,
    dense_walks, format_cigar, holds_walk, walk_alignment,
};
use crate::test_support::{allocation_peak, random_bases, rng};

fn walks_of(query: &[u8], window: &[u8], circular: bool) -> (Vec<Walk>, Vec<u8>, QueryTable) {
    let table = QueryTable::new(
        query,
        circular,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let dense = dense_walks(&table, window, 0, BandConfig::default(), 64 << 20).unwrap();
    let (walks, reverse) = (dense.walks, dense.reverse);
    (walks, reverse, table)
}

fn linear(query: &[u8]) -> QueryTable {
    QueryTable::new(
        query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap()
}

#[test]
fn periodic_diagonals_once() {
    // 10 kb of ACGT with one N, against itself: every diagonal a multiple of four apart is one
    // ungapped alignment across both N positions, on both strands. N never matches, not even N.
    let mut query = b"ACGT".repeat(2500);
    query[5000] = b'N';
    let (walks, reverse, table) = walks_of(&query, &query, false);
    let mut lines = std::collections::BTreeSet::new();
    for walk in &walks {
        let diagonal = i64::from(walk.contig) - i64::from(walk.query);
        assert!(
            lines.insert((walk.strand.is_reverse(), diagonal)),
            "one walk per diagonal"
        );
        assert_eq!(diagonal.rem_euclid(4), 0);
        // The whole overlap, but a local alignment leaves out an N in its first or last column.
        let overlap = 10_000 - diagonal.unsigned_abs();
        assert!(overlap - u64::from(walk.length) <= 2, "{walk:?}");
        let oriented = if walk.strand.is_reverse() {
            &reverse
        } else {
            &query
        };
        let alignment = walk_alignment(walk, oriented, 0, &table, BandConfig::default()).unwrap();
        alignment.validate().unwrap();
        assert!(alignment.counts().substitutions <= 2);
        let cigar = format_cigar(&alignment.runs);
        assert!(cigar.ends_with('=') && !cigar.starts_with("1X"), "{cigar}");
        if (walk.strand, diagonal) == (Strand::Forward, 0) {
            assert_eq!(cigar, "5000=1X4999=");
        }
    }
    let chained = walks.iter().filter(|walk| walk.length >= 40).count();
    assert_eq!(chained, 2 * (2 * (10_000 - 40) / 4 + 1));
}

#[test]
fn circular_walks_bounded() {
    let mut query = b"ACGT".repeat(500);
    query[1000] = b'N';
    let window = [query.clone(), query.clone()].concat();
    let (walks, reverse, table) = walks_of(&query, &window, true);
    assert!(!walks.is_empty());
    assert!(walks.iter().all(|walk| walk.length <= 2000));
    assert!(walks.iter().any(|walk| walk.length == 2000));
    for walk in &walks {
        let oriented = if walk.strand.is_reverse() {
            &reverse
        } else {
            &window
        };
        let alignment = walk_alignment(walk, oriented, 0, &table, BandConfig::default()).unwrap();
        alignment.validate().unwrap();
        assert_eq!(alignment.query_interval.len(), u64::from(walk.length));
        assert!(holds_walk(&alignment, walk, window.len(), 0, 2000));
    }
}

#[test]
fn walks_end_copy() {
    // A repeat copy inside unrelated sequence: the walk stops at the copy's end, not after the
    // X-drop tail, and starts inside it give no second walk.
    let element = b"ACG".repeat(400);
    let query = [b"GATTACATTGCCAGTTCAGG".to_vec(), element.clone()].concat();
    let window = [
        random_bases(&mut rng(7), 3000),
        element,
        random_bases(&mut rng(8), 3000),
    ]
    .concat();
    let (walks, _, _) = walks_of(&query, &window, false);
    let longest = walks
        .iter()
        .filter(|walk| !walk.strand.is_reverse())
        .max_by_key(|walk| walk.length)
        .unwrap();
    assert!((1200..1210).contains(&longest.length), "{longest:?}");
    let diagonal = i64::from(longest.contig) - i64::from(longest.query);
    let same = walks.iter().filter(|walk| {
        !walk.strand.is_reverse() && i64::from(walk.contig) - i64::from(walk.query) == diagonal
    });
    assert_eq!(same.count(), 1);
}

#[test]
fn starts_budget_fails() {
    let query = b"ACGT".repeat(2500);
    let table = linear(&query);
    assert!(dense_walks(&table, &query, 0, BandConfig::default(), 1 << 10).is_err());
}

#[test]
fn optimum_reuse_matches() {
    let query = [vec![b'A'; 100], vec![b'C'; 100]].concat();
    let contig = [vec![b'A'; 100], vec![b'G'; 20], vec![b'C'; 100]].concat();
    let table = linear(&query);
    let mut workspace = AlignmentWorkspace::default();
    let mut saved = DenseResults::default();
    let mut align = |walk: &Walk, config, saved: &mut DenseResults| {
        let input = DenseWalk {
            walk,
            later: &[],
            window: &contig,
            offset: 0,
            table: &table,
        };
        align_dense_walk(input, config, &mut workspace, saved).unwrap()
    };
    let config = BandConfig::default();
    let walk = Walk {
        strand: Strand::Forward,
        contig: 0,
        query: 0,
        length: 200,
    };
    let first = align(&walk, config, &mut saved).unwrap();
    assert_eq!(format_cigar(&first.runs), "100=20D100=");
    let other = Walk { contig: 4, ..walk };
    assert!(align(&other, config, &mut saved).is_none());
    assert_eq!(
        Some(first),
        align(&other, config, &mut DenseResults::default())
    );
    let reverse = Walk {
        strand: Strand::Reverse,
        ..walk
    };
    assert!(align(&reverse, config, &mut saved).is_some());
    let narrow = BandConfig {
        max_drift: 1,
        alignment: AlignmentConfig {
            band_width: 1,
            ..config.alignment
        },
    };
    let shifted = Walk {
        contig: 100,
        ..walk
    };
    assert!(align(&shifted, narrow, &mut saved).is_some());
}

#[test]
fn shared_kernels_match() {
    // Tandem repeats with changes give many walks on near diagonals; aligning them in one box
    // must give every walk the alignment of its own box, on both query kinds.
    let mut generator = rng(11);
    let repeat = |unit: &[u8], count: usize| unit.repeat(count);
    let mut query = random_bases(&mut generator, 3_000);
    query.extend(repeat(b"AT", 60));
    query.extend(repeat(b"GGC", 30));
    query.extend(random_bases(&mut generator, 3_000));
    let mut contig = random_bases(&mut generator, 1_000);
    contig.extend_from_slice(&query[2_500..3_000]);
    contig.extend(repeat(b"AT", 40));
    contig.push(b'C');
    contig.extend(repeat(b"AT", 30));
    contig.extend(repeat(b"GGC", 12));
    contig.extend_from_slice(b"TT");
    contig.extend(repeat(b"GGC", 15));
    contig.extend_from_slice(&query[3_210..3_700]);
    contig.extend(random_bases(&mut generator, 1_000));
    for circular in [false, true] {
        let (_, table) = QueryTable::split(
            &query,
            circular,
            AnchorConfig::default(),
            ChainConfig::default(),
        )
        .unwrap();
        let walks = dense_walks(&table, &contig, 0, BandConfig::default(), 64 << 20)
            .unwrap()
            .walks;
        assert!(walks.len() > 100, "{}", walks.len());
        let mut workspace = AlignmentWorkspace::default();
        let config = BandConfig {
            max_drift: 64,
            ..BandConfig::default()
        };
        let mut align_all = |shared: bool| {
            let mut results = DenseResults::default();
            let walks = &walks;
            (0..walks.len())
                .map(|index| {
                    let input = DenseWalk {
                        walk: &walks[index],
                        later: if shared { &walks[index + 1..] } else { &[] },
                        window: &contig,
                        offset: 7,
                        table: &table,
                    };
                    align_dense_walk(input, config, &mut workspace, &mut results).unwrap()
                })
                .collect::<Vec<_>>()
        };
        let alone = align_all(false);
        assert!(alone.iter().filter(|found| found.is_some()).count() > 100);
        assert_eq!(align_all(true), alone);
    }
}

#[test]
fn gapped_walk_widens() {
    let query = random_bases(&mut rng(7), 10_000);
    let mut contig = random_bases(&mut rng(8), 300);
    contig.extend_from_slice(&query[3_000..5_000]);
    contig.extend(random_bases(&mut rng(9), 20));
    contig.extend_from_slice(&query[5_000..7_000]);
    contig.extend(random_bases(&mut rng(10), 300));
    let walk = Walk {
        strand: Strand::Forward,
        contig: 1_300,
        query: 4_000,
        length: 300,
    };
    // The first window, 300 bases and twice the margin on each side, holds only a part of it.
    let config = BandConfig {
        max_drift: 128,
        ..BandConfig::default()
    };
    let table = linear(&query);
    let mut workspace = AlignmentWorkspace::default();
    let input = DenseWalk {
        walk: &walk,
        later: &[],
        window: &contig,
        offset: 0,
        table: &table,
    };
    let found = align_dense_walk(input, config, &mut workspace, &mut DenseResults::default())
        .unwrap()
        .unwrap();
    assert!(
        found.counts().deletions >= 20,
        "{}",
        format_cigar(&found.runs)
    );
    let aligned = found.query_interval;
    assert!(
        aligned.start <= 3_008 && aligned.end >= 6_992,
        "{aligned:?}"
    );
}

#[test]
fn short_blocks_join() {
    // Two 30-base blocks of the query with a 6-base piece between them that differs: each block
    // is shorter than the 50 columns a join must reach, but together they reach it.
    let query = random_bases(&mut rng(31), 400);
    let mut contig = query[100..166].to_vec();
    for base in &mut contig[30..36] {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    let table = linear(&query);
    let config = BandConfig::default();
    let walks = [
        Walk {
            strand: Strand::Forward,
            contig: 0,
            query: 100,
            length: 30,
        },
        Walk {
            strand: Strand::Forward,
            contig: 36,
            query: 136,
            length: 30,
        },
    ];
    let context = JoinContext {
        table: &table,
        band: config,
        limits: JoinLimits {
            min_columns: 50,
            min_identity: 0.8,
        },
    };
    let mut scratch = JoinScratch::default();
    scratch.prepare(&walks, &context);
    let joined = scratch
        .join_through(&walks, 0, &contig, 1_000, &context)
        .unwrap()
        .unwrap();
    joined.validate().unwrap();
    assert_eq!(joined.query_interval.start, 100);
    assert_eq!(joined.query_interval.end, 166);
    assert_eq!(joined.contig_interval.start, 1_000);
    assert!(JoinScratch::bytes(2, &context) > 0);
    // Past the walks a join can link, each more walk only adds its two orders.
    let (many, more) = (1 << 20, 1 << 21);
    let added = JoinScratch::bytes(more, &context) - JoinScratch::bytes(many, &context);
    assert_eq!(added, (more - many) * 2 * size_of::<u32>());
}

#[test]
fn widened_result_kept() {
    // The first range [0, 3176) ends 176 bases after the path 1000=1D2000=, so the walk widens
    // to [0, 9704) with the same query start and contig columns. A proof over the first range
    // must not answer the second.
    let query = random_bases(&mut rng(9), 10_000);
    let window = [
        &query[..1000],
        b"T",
        &query[1000..3000],
        &random_bases(&mut rng(10), 999),
    ]
    .concat();
    let table = linear(&query);
    let walk = Walk {
        strand: Strand::Forward,
        contig: 0,
        query: 0,
        length: 1000,
    };
    let input = DenseWalk {
        walk: &walk,
        later: &[],
        window: &window,
        offset: 0,
        table: &table,
    };
    let (mut workspace, mut results) = (AlignmentWorkspace::default(), DenseResults::default());
    let config = BandConfig::default();
    let found = align_dense_walk(input, config, &mut workspace, &mut results).unwrap();
    assert_eq!(format_cigar(&found.unwrap().runs), "1000=1D2000=");
}

/// Window k-mers that reach the entropy limit are not looked up in the table of the k-mers
/// below it; the walks stay the same.
#[test]
fn low_walks_unchanged() {
    let mut generator = rng(31);
    let repeat = |unit: &[u8], times: usize| unit.repeat(times);
    let mut query = random_bases(&mut generator, 2_000);
    for unit in [&b"A"[..], b"AT", b"GGC", b"AAAAAAAAAAAAG", b"CCCCTCCC"] {
        query.extend(repeat(unit, 40));
        query.extend(random_bases(&mut generator, 500));
    }
    let mut contig = random_bases(&mut generator, 800);
    contig.extend_from_slice(&query[1_500..4_000]);
    contig.extend(repeat(b"TTTTTTTTTTTTC", 9));
    contig.extend(repeat(b"GGGGAGGG", 12));
    contig.extend(random_bases(&mut generator, 800));
    for circular in [false, true] {
        let split = || {
            QueryTable::split(
                &query,
                circular,
                AnchorConfig::default(),
                ChainConfig::default(),
            )
            .unwrap()
            .1
        };
        let walks = |table: &QueryTable| {
            let dense = dense_walks(table, &contig, 0, BandConfig::default(), 64 << 20).unwrap();
            (dense.walks, dense.reverse)
        };
        let limited = walks(&split());
        assert!(limited.0.len() > 20, "{}", limited.0.len());
        assert_eq!(limited, walks(&split().without_limit()));
    }
}

#[test]
fn join_scratch_within_price() {
    // 187 walks of 15 bases, one base apart on one diagonal: one join links all of them. A later
    // window of 10 walks holds no more than its own price.
    let query = random_bases(&mut rng(33), 3_000);
    let table = linear(&query);
    let walks: Vec<Walk> = (0..187u32)
        .map(|at| Walk {
            strand: Strand::Forward,
            contig: 16 * at,
            query: 16 * at,
            length: 15,
        })
        .collect();
    let context = JoinContext {
        table: &table,
        band: BandConfig::default(),
        limits: JoinLimits {
            min_columns: 50,
            min_identity: 0.5,
        },
    };
    let mut scratch = JoinScratch::default();
    // Preparing the next window holds no more at any time than the change it makes: the last
    // window's buffers are freed before the new ones are taken.
    let prepare = |scratch: &mut JoinScratch, walks: &[Walk]| {
        let before = scratch.held_bytes() as isize;
        let peak = allocation_peak(|| scratch.prepare(walks, &context));
        assert!(
            peak <= (scratch.held_bytes() as isize - before).max(0),
            "{peak}"
        );
    };
    prepare(&mut scratch, &walks);
    prepare(&mut scratch, &walks);
    prepare(&mut scratch, &walks[..100]);
    prepare(&mut scratch, &walks);
    let joined = scratch
        .join_through(&walks, 93, &query, 0, &context)
        .unwrap()
        .unwrap();
    assert_eq!(joined.query_interval.start, 0);
    // The price of the walks, without the piece alignment.
    let price =
        |walks: usize| JoinScratch::bytes(walks, &context) - JoinScratch::bytes(0, &context);
    assert!(scratch.walk_bytes() <= price(walks.len()));
    scratch.prepare(&walks[..10], &context);
    scratch
        .join_through(&walks[..10], 5, &query, 0, &context)
        .unwrap();
    assert!(scratch.walk_bytes() <= price(10));
}

#[test]
fn join_matrix_within_price() {
    // Two-walk windows whose joins align pieces of 800 and then 1,039 bases, then a third
    // window: all buffers together, the piece matrix too, stay within the price of two walks,
    // and no window change holds more at any time than the change it makes.
    let query = random_bases(&mut rng(34), 3_000);
    let table = linear(&query);
    let context = JoinContext {
        table: &table,
        band: BandConfig::default(),
        limits: JoinLimits {
            min_columns: 50,
            min_identity: 0.0,
        },
    };
    let pair = |gap: u32| {
        [0, 30 + gap].map(|start| Walk {
            strand: Strand::Forward,
            contig: start,
            query: start,
            length: 30,
        })
    };
    let price = JoinScratch::bytes(2, &context);
    let mut scratch = JoinScratch::default();
    for gap in [800, 1_039, 800] {
        let walks = pair(gap);
        let before = scratch.held_bytes() as isize;
        let peak = allocation_peak(|| scratch.prepare(&walks, &context));
        assert!(
            peak <= (scratch.held_bytes() as isize - before).max(0),
            "{peak}"
        );
        let joined = scratch
            .join_through(&walks, 0, &query, 0, &context)
            .unwrap();
        assert!(joined.is_some());
        assert!(scratch.held_bytes() <= price, "{}", scratch.held_bytes());
    }
}

#[test]
fn join_transitions_within_price() {
    // Windows of 2 and 100,000 walks with reaches of 1,039 and 527 bases, joined or not: at
    // every change of window all buffers together, the piece matrix too, stay within the larger
    // price of the two windows, and after it within the price of the new one.
    let query = random_bases(&mut rng(35), 3_000);
    let table = linear(&query);
    let context = |max_drift: u32| JoinContext {
        table: &table,
        band: BandConfig {
            max_drift,
            ..BandConfig::default()
        },
        limits: JoinLimits {
            min_columns: 50,
            min_identity: 0.0,
        },
    };
    let (wide, narrow) = (context(1_024), context(512));
    // Two walks 800 bases apart, beyond the reach of 527 bases, and 500 bases apart.
    let pair_of = |second: u32| -> Vec<Walk> {
        [0, second]
            .map(|start| Walk {
                strand: Strand::Forward,
                contig: start,
                query: start,
                length: 30,
            })
            .to_vec()
    };
    let (pair, near_pair) = (pair_of(830), pair_of(530));
    let many: Vec<Walk> = (0..100_000u32)
        .map(|at| Walk {
            strand: Strand::Forward,
            contig: 7 * at,
            query: (13 * at) % 2_900,
            length: 20,
        })
        .collect();
    // (walks, context, whether a join is tried and then whether it is made)
    let windows = [
        (&pair, &wide, None),
        (&many, &narrow, None),
        (&pair, &wide, Some(true)),
        (&many, &narrow, None),
        (&near_pair, &narrow, Some(true)),
        (&pair, &wide, Some(true)),
        (&many, &wide, None),
        (&pair, &narrow, Some(false)),
    ];
    let mut scratch = JoinScratch::default();
    let mut last_price = 0;
    for (walks, context, join) in windows {
        let price = JoinScratch::bytes(walks.len(), context);
        let before = scratch.held_bytes() as isize;
        let peak = allocation_peak(|| scratch.prepare(walks, context));
        assert!(
            before + peak <= price.max(last_price) as isize,
            "{before} {peak}"
        );
        if let Some(made) = join {
            let joined = scratch.join_through(walks, 0, &query, 0, context).unwrap();
            assert_eq!(joined.is_some(), made);
        }
        assert!(scratch.held_bytes() <= price, "{}", scratch.held_bytes());
        last_price = price;
    }
}

#[test]
fn parallel_walks_same() {
    // A query of diverged repeat copies and tandem runs, and a window of its copies on both
    // strands: the starts and walks found on four threads equal those of one loop, also for a
    // circular query whose walks pass its origin, and the budget fails the same way.
    use super::dense::PARALLEL_FROM;
    use crate::seeds::reverse_complement;
    use crate::test_support::with_substitutions;
    use std::sync::atomic::Ordering;
    let mut random = rng(41);
    let element = random_bases(&mut random, 300);
    let mut query = Vec::new();
    while query.len() < 20_000 {
        query.extend(random_bases(&mut random, 40));
        query.extend(with_substitutions(&mut random, &element, 0.92));
        query.extend(with_substitutions(&mut random, &b"AGT".repeat(40), 0.95));
    }
    let mut window = Vec::new();
    for copy in 0..6 {
        let changed = with_substitutions(&mut random, &query, 0.97);
        window.extend(match copy % 2 {
            0 => changed,
            _ => reverse_complement(&changed),
        });
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    for circular in [false, true] {
        let table = QueryTable::new(
            &query,
            circular,
            AnchorConfig::default(),
            ChainConfig::default(),
        )
        .unwrap();
        let run = |from: usize, budget: usize| {
            PARALLEL_FROM.store(from, Ordering::Relaxed);
            let config = BandConfig::default();
            pool.install(|| dense_walks(&table, &window, 0, config, budget))
        };
        let alone = run(usize::MAX, 256 << 20).unwrap();
        let spread = run(0, 256 << 20).unwrap();
        assert_eq!(spread.walks, alone.walks, "circular {circular}");
        assert_eq!(spread.reverse, alone.reverse);
        assert!(alone.walks.len() > 50_000, "{}", alone.walks.len());
        let small = (
            run(usize::MAX, 1 << 20).unwrap_err(),
            run(0, 1 << 20).unwrap_err(),
        );
        assert_eq!(small.0.to_string(), small.1.to_string());
    }
    PARALLEL_FROM.store(1 << 20, Ordering::Relaxed);
}
