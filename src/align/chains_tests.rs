//! Tests of chaining: copies, circular origins, repeats, periodic queries and drift.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::anchors::{AnchorConfig, find_anchors};
use crate::seeds::reverse_complement;
use crate::test_support::{random_bases, rng};

fn chains(query: &[u8], circular: bool, window: &[u8]) -> Vec<Chain> {
    let table = QueryTable::new(
        query,
        circular,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    window_chains(&table, window, &mut Vec::new())
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
    let found = window_chains(&table, &query, &mut scratch).unwrap();
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
    let found = window_chains(&table, &query, &mut scratch).unwrap();
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
    let table = QueryTable::new(
        &poly,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut anchors = Vec::new();
    find_anchors(&table, &[b'A'; 500], &mut anchors, None, 0);
    assert!(anchors.is_empty());
    let low = QueryTable::low_complexity(
        &poly,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
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
