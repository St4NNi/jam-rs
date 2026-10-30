//! Tests of chaining: copies, circular origins, repeats, periodic queries and drift.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::anchors::{AnchorConfig, find_anchors};
use crate::seeds::reverse_complement;
use crate::test_support::{random_bases, rng};
use rand::RngExt;

fn chains(query: &[u8], circular: bool, window: &[u8]) -> Vec<Chain> {
    let table = QueryTable::new(
        query,
        circular,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    window_chains(&table, window, &mut Vec::new(), &|| false)
        .unwrap()
        .chains
}

fn position(anchor: Anchor) -> (u32, u32) {
    (anchor.contig, anchor.query)
}

fn diagonals(chain: &Chain) -> Vec<i64> {
    let mut diagonals: Vec<i64> = chain
        .anchors
        .iter()
        .zip(&chain.offsets)
        .map(|(anchor, &offset)| i64::from(anchor.contig) - i64::from(offset))
        .collect();
    diagonals.dedup();
    diagonals
}

#[test]
fn cyclic_whole_periods() {
    let mut query = b"ACGT".repeat(300);
    let window = query.clone();
    let whole = QueryTable::new(
        &query,
        true,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    assert!(whole.periodic_pair(&window).is_some());
    query.push(b'A');
    let linear = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    assert_eq!(
        linear
            .periodic_pair(&window)
            .map(|periodic| periodic.period),
        Some(4)
    );
    let circular = QueryTable::new(
        &query,
        true,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    assert!(circular.periodic_pair(&window).is_none());
    let random = QueryTable::new(
        &random_bases(&mut rng(91), 20_000),
        true,
        AnchorConfig::default(),
        ChainConfig::default(),
    );
    assert!(random.unwrap().periodic_pair(&window).is_none());
}

#[test]
fn periodic_copy_chains() {
    let query = b"ACGT".repeat(300);
    let found = chains(&query, false, &query);
    assert!(found.iter().any(|chain| {
        !chain.strand.is_reverse()
            && position(chain.anchors[0]) == (0, 0)
            && chain.anchors.last().unwrap().query >= 1180
    }));
}

#[test]
fn dense_product_walked() {
    let mut query = b"ACGT".repeat(25_000);
    query[50_000] = b'N';
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut scratch = Vec::new();
    // Every k-mer occurs about 25,000 times per strand: none gives anchors, all are left to
    // dense walks, and nothing of the 1.25 billion anchor product is held.
    let found = window_chains(&table, &query, &mut scratch, &|| false).unwrap();
    let (chains, above) = (found.chains, found.dense_above);
    assert!(chains.is_empty());
    assert_eq!(above, Some(0));
    assert_eq!(scratch.capacity(), 0);
}

#[test]
fn periodic_chains_enumerated() {
    let query = b"ACGT".repeat(25_000);
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut scratch = Vec::new();
    let found = window_chains(&table, &query, &mut scratch, &|| false).unwrap();
    assert!(found.chains.len() > 99_000);
    assert_eq!(found.dense_above, None);
    assert_eq!(scratch.capacity(), 0);
    assert!(found.chains.iter().all(|chain| chain.anchors.len() == 3));
}

#[test]
fn chains_both_strands() {
    let query = random_bases(&mut rng(1), 3000);
    let window = [
        random_bases(&mut rng(2), 500),
        query[1000..2500].to_vec(),
        random_bases(&mut rng(3), 500),
    ]
    .concat();
    let found = chains(&query, false, &window);
    assert_eq!(found.len(), 1);
    let chain = &found[0];
    assert!(!chain.strand.is_reverse());
    assert_eq!(position(chain.anchors[0]), (500, 1000));
    assert!(chain.anchors.last().unwrap().query >= 2480);
    assert_eq!(diagonals(chain), [500]);

    let window = [
        random_bases(&mut rng(4), 300),
        reverse_complement(&query[200..1800]),
    ]
    .concat();
    let found = chains(&query, false, &window);
    assert_eq!(found.len(), 1);
    assert!(found[0].strand.is_reverse());
    // In the reverse complement of the window the copy starts at 0 and holds query 200..1800.
    let first = found[0].anchors[0];
    assert_eq!(i64::from(first.contig) - i64::from(first.query), -200);
    assert!(found[0].anchors.last().unwrap().query >= 1780);
}

#[test]
fn origin_chain_splits() {
    let query = random_bases(&mut rng(5), 2000);
    let window = [query[1500..].to_vec(), query[..600].to_vec()].concat();
    let found = chains(&query, true, &window);
    assert_eq!(found.len(), 1);
    assert_eq!(position(found[0].anchors[0]), (0, 1500));
    assert!(*found[0].offsets.last().unwrap() >= 1080);
    // One diagonal across the origin, counted from the first anchor at window position 0.
    assert_eq!(diagonals(&found[0]), [0]);
    let linear = chains(&query, false, &window);
    assert_eq!(linear.len(), 2);
    let mut parts: Vec<(u32, u32)> = linear
        .iter()
        .map(|chain| (chain.anchors[0].query, *chain.offsets.last().unwrap()))
        .collect();
    parts.sort_unstable();
    assert_eq!((parts[0].0, parts[1].0), (0, 1500));
    assert!(parts[0].1 <= 600 - 15 && parts[1].1 <= 500 - 15);
}

#[test]
fn repeats_chain_apart() {
    let repeat = random_bases(&mut rng(6), 800);
    let query = [
        repeat.clone(),
        random_bases(&mut rng(7), 1000),
        repeat.clone(),
    ]
    .concat();
    let found = chains(&query, false, &repeat);
    let mut starts: Vec<u32> = found.iter().map(|chain| chain.anchors[0].query).collect();
    starts.sort_unstable();
    assert_eq!(starts, [0, 1800]);
    // Two copies of a circular query in a row: no chain covers more than one query length.
    let plasmid = random_bases(&mut rng(8), 1500);
    let window = [plasmid.clone(), plasmid.clone()].concat();
    let found = chains(&plasmid, true, &window);
    assert!(found.len() >= 2);
    assert!(
        found
            .iter()
            .all(|chain| chain.offsets.last().unwrap() + 15 <= 1500)
    );
}

#[test]
fn mutated_copy_chains() {
    let query = random_bases(&mut rng(9), 4000);
    let mut copy = query.clone();
    for position in (25..copy.len()).step_by(60) {
        copy[position] = if copy[position] == b'A' { b'C' } else { b'A' };
    }
    let found = chains(&query, false, &copy);
    assert_eq!(found.len(), 1);
    assert!(found[0].anchors.len() > 100);
    let poly = [vec![b'A'; 300], random_bases(&mut rng(10), 100)].concat();
    let (table, low) = QueryTable::split(
        &poly,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut anchors = Vec::new();
    find_anchors(&table, &[b'A'; 500], &mut anchors, None, 0);
    assert!(anchors.is_empty());
    find_anchors(&low, &[b'A'; 500], &mut anchors, None, 0);
    assert!(!anchors.is_empty());
    assert!(
        AnchorConfig {
            stride: 3,
            ..AnchorConfig::default()
        }
        .validate()
        .is_err()
    );
}

#[test]
fn drift_splits_chain() {
    let chain = Chain {
        strand: Strand::Forward,
        anchors: [
            (0, 0),
            (100, 100),
            (200, 200),
            (900, 300),
            (1000, 400),
            (1100, 500),
        ]
        .map(|(contig, query)| Anchor {
            strand: Strand::Forward,
            contig,
            query,
        })
        .to_vec(),
        offsets: vec![0, 100, 200, 300, 400, 500],
        score: 90,
    };
    let parts = split_by_drift(&chain, 64, 3);
    assert_eq!(parts.len(), 2);
    assert_eq!(position(parts[1].anchors[0]), (900, 300));
    assert_eq!(parts[1].offsets, [0, 100, 200]);
    assert_eq!(split_by_drift(&chain, 1_000, 3).len(), 1);
}

#[test]
fn split_by_entropy() {
    let mut generator = rng(12);
    let query = [
        random_bases(&mut generator, 2000),
        b"AT".repeat(50),
        vec![b'A'; 40],
        b"AACAAG".repeat(10),
        random_bases(&mut generator, 2000),
    ]
    .concat();
    let config = AnchorConfig::default();
    let (table, low) = QueryTable::split(&query, false, config, ChainConfig::default()).unwrap();
    let mut low_kmers = 0;
    crate::seeds::kmers(&query, config.a, |kmer| {
        let mut counts = [0u32; 4];
        for at in 0..config.a {
            counts[((kmer.forward >> (2 * at)) & 3) as usize] += 1;
        }
        let entropy: f64 = counts
            .iter()
            .filter(|&&count| count > 0)
            .map(|&count| {
                let share = f64::from(count) / f64::from(config.a);
                -share * share.log2()
            })
            .sum();
        let passes = entropy >= config.min_entropy;
        low_kmers += usize::from(!passes);
        assert_eq!(table.multiplicity(kmer.forward) > 0, passes);
        assert_eq!(low.multiplicity(kmer.forward) > 0, !passes);
    });
    assert!(low_kmers > 100);
}

#[test]
fn copy_across_pieces() {
    // The query is scanned in pieces of 65,536 bases; a copy across the first border chains
    // as one, with anchors on both sides and in the last a - 1 bases before it.
    let query = random_bases(&mut rng(11), 150_000);
    let found = chains(&query, false, &query[60_000..72_000]);
    assert_eq!(found.len(), 1);
    let queries: Vec<u32> = found[0].anchors.iter().map(|anchor| anchor.query).collect();
    assert!(queries.iter().any(|&at| at < 65_000));
    assert!(queries.iter().any(|&at| (65_522..65_536).contains(&at)));
    assert!(queries.iter().any(|&at| at > 66_000));
    assert!(queries.len() > 2_000);
}

#[test]
fn retained_capacity_ignored() {
    // A buffer kept from an earlier, larger window must not change which anchors this window
    // takes: a capacity that alone fills the budget used to send every k-mer to dense walks.
    let query = random_bases(&mut rng(17), 2000);
    let window = [&random_bases(&mut rng(18), 500)[..], &query[200..1800]].concat();
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let fresh = window_chains(&table, &window, &mut Vec::new(), &|| false).unwrap();
    let mut retained = Vec::with_capacity(CHAINING_BYTES / std::mem::size_of::<Anchor>());
    let reused = window_chains(&table, &window, &mut retained, &|| false).unwrap();
    assert_eq!(fresh.dense_above, None);
    assert_eq!(reused.dense_above, fresh.dense_above);
    assert_eq!(reused.chains, fresh.chains);
    assert!(retained.capacity() < CHAINING_BYTES / std::mem::size_of::<Anchor>());
}

#[test]
fn chaining_work_bounded() {
    // One 15-mer at 40 query positions 6,000 apart and 50 window positions 20 apart: no
    // predecessor is valid, as the query steps are 0 or above the largest gap, yet each anchor
    // examines copies before it. The work ends in an error, never in fewer chains.
    let query = random_bases(&mut rng(10), 240_000);
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let copies: Vec<Anchor> = (0..50u32)
        .flat_map(|j| {
            (0..40u32).map(move |i| Anchor {
                strand: Strand::Forward,
                contig: 20 * j,
                query: 6000 * i,
            })
        })
        .collect();
    let run = |anchors: &[Anchor], limit: u64, stop: bool| {
        let stopped = move || stop;
        let mut examined = Examined::new(limit, &stopped);
        let mut chains = Vec::new();
        let result = chain_anchors(&table, anchors, &mut chains, &mut examined).map(|()| chains);
        (result, examined.count)
    };
    let (result, work) = run(&copies, u64::MAX, false);
    assert_eq!(result, Ok(Vec::new()));
    // Fewer than the 40^2 x 50 x 49 / 2 = 1,960,000 predecessors before them, as once the work
    // passes the lookback of every anchor only those near the diagonal are examined.
    assert!(work < 1_960_000);
    assert_eq!(run(&copies, work, false).0, Ok(Vec::new()));
    assert_eq!(
        run(&copies, work - 1, false).0,
        Err(AlignmentError::ChainingBudget(Budgeted::ChainingWork))
    );
    // Copies at one query position are near the diagonal of each other: millions are examined.
    let line: Vec<Anchor> = (0..3_000u32)
        .map(|contig| Anchor {
            strand: Strand::Forward,
            contig,
            query: 0,
        })
        .collect();
    assert_eq!(run(&line, u64::MAX, true).0, Err(AlignmentError::Stopped));
}

#[test]
fn single_anchor_chains() {
    let query = random_bases(&mut rng(31), 1_600);
    let mut random = rng(32);
    // Eighteen query bases at a looked-up window position give exactly one anchor.
    let window = [
        random_bases(&mut random, 300),
        query[100..118].to_vec(),
        random_bases(&mut random, 300),
    ]
    .concat();
    let short = chains(&query[..1_000], false, &window);
    assert_eq!(short.len(), 1);
    assert_eq!(position(short[0].anchors[0]), (300, 100));
    assert!(chains(&query, false, &window).is_empty());
}

/// The best link of every anchor as found by examining every predecessor one by one.
fn links_one_by_one(table: &QueryTable, anchors: &[Anchor]) -> Vec<(i64, usize, i64)> {
    let a = i64::from(table.anchors().a);
    let max_gap = i64::from(table.chains().max_gap);
    let length = i64::from(table.query_length());
    let mut links: Vec<(i64, usize, i64)> = Vec::new();
    for (current, anchor) in anchors.iter().enumerate() {
        let (mut best, mut from, mut best_span) = (a, usize::MAX, 0);
        let end = anchors[..current].partition_point(|earlier| earlier.contig < anchor.contig);
        let mut eligible = 0;
        for earlier in (0..end).rev() {
            let contig_step = i64::from(anchor.contig) - i64::from(anchors[earlier].contig);
            if contig_step > max_gap {
                break;
            }
            let mut query_step = i64::from(anchor.query) - i64::from(anchors[earlier].query);
            if table.circular() && query_step <= 0 {
                query_step += length;
            }
            if query_step <= 0 || query_step > max_gap || links[earlier].2 + query_step + a > length
            {
                continue;
            }
            let drift = (contig_step - query_step).abs();
            if drift > DRIFT_BASES + contig_step.min(query_step) / DRIFT_SHARE {
                continue;
            }
            eligible += 1;
            let gap_cost = if drift == 0 {
                0
            } else {
                drift / DRIFT_COST + 1 + i64::from(64 - (drift as u64).leading_zeros()) / 2
            };
            let gained = links[earlier].0 + contig_step.min(query_step).min(a) - gap_cost;
            if gained > best {
                (best, from, best_span) = (gained, earlier, links[earlier].2 + query_step);
            }
            if eligible == LOOKBACK {
                break;
            }
        }
        links.push((best, from, best_span));
    }
    links
}

#[test]
fn runs_link_as_one_by_one() {
    // Runs of up to 80 anchors at one window position, some near a few diagonals and the rest
    // anywhere on the query, so runs are searched and many steps wrap around a circular origin.
    let query = random_bases(&mut rng(40), 9_000);
    let mut random = rng(41);
    for circular in [false, true] {
        let table = QueryTable::new(
            &query,
            circular,
            AnchorConfig::default(),
            ChainConfig::default(),
        )
        .unwrap();
        let mut anchors = Vec::new();
        for contig in (0..12_000u32).step_by(7) {
            for _ in 0..random.random_range(0..80u32) {
                let query = if random.random_range(0..3) == 0 {
                    let diagonal = [0, 3_000, 8_500][random.random_range(0..3)];
                    (contig + diagonal + random.random_range(0..120)) % 9_000
                } else {
                    random.random_range(0..9_000)
                };
                anchors.push(Anchor {
                    strand: Strand::Forward,
                    contig,
                    query,
                });
            }
        }
        anchors.sort_unstable();
        anchors.dedup();
        let mut examined = Examined::new(u64::MAX, &|| false);
        let links = chain_links(&table, &anchors, &mut examined).unwrap();
        let found: Vec<(i64, usize, i64)> = links
            .iter()
            .map(|link| (link.score, link.previous, link.span))
            .collect();
        assert_eq!(found, links_one_by_one(&table, &anchors));
        assert!(found.iter().filter(|link| link.1 != usize::MAX).count() > anchors.len() / 4);
    }
}

#[test]
#[ignore = "timing: run in release with --nocapture"]
fn dense_runs_timed() {
    // As the slowest window of SAMEA6809557: 22,293 bases with anchors at every third position,
    // about 44 per position, against a query of 83,586 bases.
    let query = random_bases(&mut rng(50), 83_586);
    let mut random = rng(51);
    for circular in [false, true] {
        let table = QueryTable::new(
            &query,
            circular,
            AnchorConfig::default(),
            ChainConfig::default(),
        )
        .unwrap();
        let mut anchors = Vec::new();
        for contig in (0..22_293u32).step_by(3) {
            for _ in 0..random.random_range(0..88u32) {
                anchors.push(Anchor {
                    strand: Strand::Forward,
                    contig,
                    query: random.random_range(0..83_586),
                });
            }
        }
        anchors.sort_unstable();
        anchors.dedup();
        let clock = std::time::Instant::now();
        let expected = links_one_by_one(&table, &anchors);
        let one_by_one = clock.elapsed();
        let clock = std::time::Instant::now();
        let mut examined = Examined::new(u64::MAX, &|| false);
        let links = chain_links(&table, &anchors, &mut examined).unwrap();
        let searched = clock.elapsed();
        let found: Vec<(i64, usize, i64)> = links
            .iter()
            .map(|link| (link.score, link.previous, link.span))
            .collect();
        assert_eq!(found, expected);
        println!(
            "circular {circular} anchors {} one by one {one_by_one:?} searched {searched:?} examined {}",
            anchors.len(),
            examined.count
        );
    }
}
