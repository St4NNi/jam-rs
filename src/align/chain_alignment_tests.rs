//! Tests of chain alignment: copies, divergent stretches, drift, indels and circular queries.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::align::chains::tests::{random_sequence, reverse_complement};
use crate::align::{
    Anchor, AnchorConfig, Chain, ChainConfig, EditOperation, Strand, format_cigar, split_by_drift,
    window_chains,
};

/// Substitutes every base whose roll is below `per_thousand` of 1,000.
fn substitute(sequence: &mut [u8], per_thousand: u64, mut state: u64) {
    for base in sequence.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        if state % 1000 < per_thousand {
            *base = match *base {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
    }
}

/// The alignments of every chain part of `window`, best chains first. Windows start at contig
/// position 0.
fn alignments(query: &[u8], window: &[u8], circular: bool) -> Vec<Alignment> {
    let config = BandConfig::default();
    let table = QueryTable::new(
        query,
        circular,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut workspace = AlignmentWorkspace::default();
    let chains = window_chains(&table, window, &mut Vec::new())
        .unwrap()
        .chains;
    let mut found = Vec::new();
    for chain in chains {
        for part in split_by_drift(&chain, config.max_drift, 3) {
            let reach = query.len() as u64;
            let aligned = align_chain(&part, window, 0, &table, config, &mut workspace, reach);
            if let Some(alignment) = aligned.unwrap() {
                alignment.validate().unwrap();
                found.push(alignment);
            }
        }
    }
    found
}

#[test]
fn widens_without_anchors() {
    // The copy differs from the query at every 12th base outside the middle 600 bases, so only
    // the middle holds anchors. The alignment still reaches both ends: the first window around
    // the chain is much narrower than the 3,000 bases that pair.
    let query = random_sequence(3000, 91);
    let mut window = query.clone();
    for at in (0..3000).filter(|at| !(1200..1800).contains(at) && at % 12 == 6) {
        window[at] = match window[at] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let found = alignments(&query, &window, false);
    assert_eq!(found.len(), 1);
    assert!(found[0].query_interval.len() >= 2900, "{:?}", found[0]);
}

#[test]
fn divergent_tail_kept() {
    let query = random_sequence(10_000, 42);
    let mut contig = query.clone();
    for at in (6004..10_000).step_by(5) {
        contig[at] = match contig[at] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
    }
    let found = alignments(&query, &contig, false);
    assert!(
        found
            .iter()
            .any(|alignment| alignment.query_interval == Interval::new(0, 10_000).unwrap())
    );
}

#[test]
fn copies_both_strands() {
    let query = random_sequence(3000, 1);
    let window = [
        random_sequence(700, 2),
        query[500..2500].to_vec(),
        random_sequence(300, 3),
    ]
    .concat();
    let found = alignments(&query, &window, false);
    assert_eq!(found.len(), 1);
    // A random flank base can match by chance and extend a local alignment by a few bases.
    let aligned = found[0].query_interval;
    assert!((495..=500).contains(&aligned.start) && (2500..=2505).contains(&aligned.end));
    assert_eq!(found[0].contig_interval.start, aligned.start + 200);
    assert_eq!(found[0].strand, Strand::Forward);
    assert_eq!(found[0].counts().matches, aligned.len());

    let reverse = reverse_complement(&query[100..1900]);
    let found = alignments(&query, &reverse, false);
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].strand, Strand::Reverse);
    let aligned = found[0].query_interval;
    assert!((95..=100).contains(&aligned.start) && aligned.end == 1900);
    assert_eq!(found[0].contig_interval.start, 0);
}

#[test]
fn diverged_block_kept() {
    let query = random_sequence(3000, 4);
    let mut copy = query.clone();
    substitute(&mut copy[1000..2000], 250, 5);
    let found = alignments(&query, &copy, false);
    assert_eq!(found[0].query_interval, Interval::new(0, 3000).unwrap());
    let identity = found[0].identity();
    assert!((0.88..0.95).contains(&identity), "{identity}");
}

#[test]
fn spans_circular_origin() {
    let query = random_sequence(2000, 6);
    let window = [query[1500..].to_vec(), query[..1500].to_vec()].concat();
    let found = alignments(&query, &window, true);
    assert_eq!(found.len(), 1);
    // Query positions run on over the origin.
    assert_eq!(found[0].query_interval, Interval::new(1500, 3500).unwrap());
    assert_eq!(format_cigar(&found[0].runs), "2000=");
}

#[test]
fn whole_linear_copy() {
    // The first anchor of the reverse copy is at query 1: the one base before it must still
    // align, as the window has no base after the copy.
    let query = random_sequence(3000, 9);
    let found = alignments(&query, &reverse_complement(&query), true);
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].query_interval, Interval::new(0, 3000).unwrap());
    assert_eq!(format_cigar(&found[0].runs), "3000=");
}

#[test]
fn drift_long_deletion() {
    // A band wider than the drift margin and a chain drifting 1,000 bases: the alignment pairs
    // the whole query with 1,000 contig bases between.
    let query = random_sequence(3000, 21);
    let window = [
        random_sequence(2000, 22),
        query[..1500].to_vec(),
        random_sequence(1000, 23),
        query[1500..].to_vec(),
        random_sequence(2000, 24),
    ]
    .concat();
    let config = BandConfig {
        alignment: AlignmentConfig {
            band_width: 1_100,
            ..AlignmentConfig::default()
        },
        ..BandConfig::default()
    };
    let anchors = [0u32, 500, 1000, 1500, 2000, 2500].map(|query| Anchor {
        strand: Strand::Forward,
        contig: if query < 1500 {
            2000 + query
        } else {
            3000 + query
        },
        query,
    });
    let chain = Chain {
        strand: Strand::Forward,
        anchors: anchors.to_vec(),
        offsets: anchors.iter().map(|anchor| anchor.query).collect(),
        score: 100,
    };
    let table = QueryTable::new(
        &query,
        false,
        AnchorConfig::default(),
        ChainConfig::default(),
    )
    .unwrap();
    let mut workspace = AlignmentWorkspace::default();
    let alignment = align_chain(&chain, &window, 0, &table, config, &mut workspace, 3000)
        .unwrap()
        .unwrap();
    assert_eq!(alignment.query_interval, Interval::new(0, 3000).unwrap());
    assert_eq!(alignment.counts().deletions, 1000);
}

#[test]
fn long_insert_splits() {
    let query = random_sequence(4000, 7);
    let copy = [
        &query[..1500],
        &query[1530..2500],
        b"ACGTTGCAAGCTAGGCTTAA",
        &query[2500..],
    ]
    .concat();
    let found = alignments(&query, &copy, false);
    assert_eq!(found.len(), 1);
    let runs = &found[0].runs;
    assert!(
        runs.iter()
            .any(|run| run.operation == EditOperation::Insertion)
    );
    assert!(
        runs.iter()
            .any(|run| run.operation == EditOperation::Deletion)
    );
    assert_eq!(found[0].query_interval, Interval::new(0, 4000).unwrap());
    let inserted = [
        &query[..2000],
        &random_sequence(2000, 8)[..],
        &query[2000..],
    ]
    .concat();
    let found = alignments(&query, &inserted, false);
    let mut intervals: Vec<_> = found
        .iter()
        .map(|alignment| alignment.query_interval)
        .collect();
    intervals.sort_unstable();
    assert_eq!(intervals.len(), 2);
    assert!(intervals[0].end >= 1990 && intervals[1].start <= 2010);
}

#[test]
fn periodic_exact_diagonals() {
    let query = b"ACGT".repeat(300);
    let found = alignments(&query, &query, false);
    let mappings: std::collections::BTreeSet<_> = found
        .iter()
        .filter(|alignment| alignment.counts().matches >= 100)
        .map(|alignment| {
            let contig = if alignment.strand == Strand::Reverse {
                1200 - alignment.contig_interval.end
            } else {
                alignment.contig_interval.start
            };
            (
                alignment.strand,
                contig as i64 - alignment.query_interval.start as i64,
            )
        })
        .collect();
    let expected = [Strand::Forward, Strand::Reverse]
        .into_iter()
        .flat_map(|strand| (-1100..=1100).step_by(4).map(move |d| (strand, d)))
        .collect();
    assert_eq!(mappings, expected);
}
