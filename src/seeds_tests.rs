//! Tests of seed selection, query k-mers and the entropy filter against direct computations.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::test_support::{random_bases, rng};

fn scheme(k: u8, w: u16, min_entropy: f64) -> SeedScheme {
    SeedScheme { k, w, min_entropy }
}

fn collect(sequence: &[u8], scheme: SeedScheme) -> Vec<Kmer> {
    let mut seeds = Vec::new();
    Seeder::new(scheme)
        .unwrap()
        .seeds(sequence, |seed| seeds.push(seed));
    seeds
}

fn encode(kmer: &[u8]) -> (u64, u64) {
    let forward = kmer.iter().fold(0u64, |value, byte| {
        (value << 2) | u64::from(base_code(*byte).unwrap())
    });
    let reverse = kmer.iter().rev().fold(0u64, |value, byte| {
        (value << 2) | (3 - u64::from(base_code(*byte).unwrap()))
    });
    (forward, reverse)
}

fn entropy_of(kmer: &[u8]) -> f64 {
    let mut counts = [0usize; 4];
    for byte in kmer {
        counts[usize::from(base_code(*byte).unwrap())] += 1;
    }
    shannon_entropy(counts, kmer.len())
}

/// The leftmost minimum of every full minimizer window of `hashes` (one clean run starting at `offset`),
/// or of the whole run when it is shorter than one.
fn sliding_minima(hashes: &[u64], offset: usize, w: usize, selected: &mut Vec<u64>) {
    if hashes.is_empty() {
        return;
    }
    let first = selected.len();
    for end in w..=hashes.len() {
        let slice = &hashes[end - w..end];
        let best = (0..w).min_by_key(|i| (slice[*i], *i)).unwrap();
        let position = (offset + end - w + best) as u64;
        if selected[first..].last() != Some(&position) {
            selected.push(position);
        }
    }
    if hashes.len() < w {
        let best = (0..hashes.len()).min_by_key(|i| (hashes[*i], *i)).unwrap();
        selected.push((offset + best) as u64);
    }
}

/// Clean runs of `sequence` with the offset of their first base.
fn clean_runs(sequence: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    sequence
        .split(|byte| base_code(*byte).is_none())
        .map(move |run| (run.as_ptr() as usize - sequence.as_ptr() as usize, run))
}

fn reference(sequence: &[u8], k: usize, w: usize) -> Vec<u64> {
    let mut selected = Vec::new();
    for (offset, run) in clean_runs(sequence).filter(|(_, run)| run.len() >= k) {
        let hashes: Vec<u64> = run
            .windows(k)
            .map(|kmer| {
                let (forward, reverse) = encode(kmer);
                jamhash::jamhash_u64(forward.min(reverse))
            })
            .collect();
        sliding_minima(&hashes, offset, w, &mut selected);
    }
    selected
}

#[test]
fn seeds_match_minima() {
    for (round, (k, w)) in [(5, 1), (7, 4), (15, 10), (21, 50), (25, 100), (32, 7)]
        .into_iter()
        .enumerate()
    {
        let mut sequence = random_bases(&mut rng(0x1234_5678 + round as u64), 5000);
        for at in [700, 701, 2500, 4990] {
            sequence[at] = b'N';
        }
        sequence[3000..3010].copy_from_slice(b"acgtRYacgt");
        let positions: Vec<u64> = collect(&sequence, scheme(k, w, 0.0))
            .iter()
            .map(|seed| seed.position)
            .collect();
        assert_eq!(
            positions,
            reference(&sequence, k as usize, w as usize),
            "k {k}"
        );
    }
}

#[test]
fn single_kmer_windows() {
    let mut sequence = random_bases(&mut rng(0x9e37_79b9), 3000);
    sequence[1000..1060].fill(b'A');
    sequence[1500] = b'N';
    sequence[2000..2004].copy_from_slice(b"acgt");
    let mut expected = Vec::new();
    for start in 0..=sequence.len() - 25 {
        let kmer = &sequence[start..start + 25];
        if kmer.iter().any(|byte| base_code(*byte).is_none()) || entropy_of(kmer) < 1.0 {
            continue;
        }
        let (forward, reverse) = encode(kmer);
        expected.push(Kmer {
            hash: jamhash::jamhash_u64(forward.min(reverse)),
            position: start as u64,
            reverse: reverse < forward,
        });
    }
    assert!(expected.len() > 2500 && expected.len() < 2976);
    let ordinary: Vec<Kmer> = collect(&sequence, scheme(25, 1, 1.0))
        .into_iter()
        .filter(|seed| entropy_of(&sequence[seed.position as usize..][..25]) >= 1.0)
        .collect();
    assert_eq!(ordinary, expected);
}

#[test]
fn strands_agree() {
    let sequence = random_bases(&mut rng(99), 20_000);
    let reverse = reverse_complement(&sequence);
    let scheme = scheme(25, 100, 0.0);
    let forward = collect(&sequence, scheme);
    let mut mirrored: Vec<Kmer> = collect(&reverse, scheme)
        .into_iter()
        .map(|seed| Kmer {
            position: sequence.len() as u64 - 25 - seed.position,
            reverse: !seed.reverse,
            ..seed
        })
        .collect();
    mirrored.reverse();
    assert_eq!(forward, mirrored);
    let density = forward.len() as f64 / sequence.len() as f64;
    assert!((0.015..0.025).contains(&density), "density {density}");
    let gaps = forward
        .windows(2)
        .map(|pair| pair[1].position - pair[0].position);
    assert!(gaps.max().unwrap() <= 100);
}

#[test]
fn shared_windows_agree() {
    let left = random_bases(&mut rng(5), 3000);
    let mut right = random_bases(&mut rng(6), 2000);
    right.extend_from_slice(&left[1000..2000]);
    let scheme = scheme(21, 50, 1.0);
    let hashes = |seeds: Vec<Kmer>, from: u64, to: u64| -> Vec<u64> {
        seeds
            .into_iter()
            .filter(|seed| (from..to).contains(&seed.position))
            .map(|seed| seed.hash)
            .collect()
    };
    // Seeds whose whole minimizer window lies in the shared stretch are selected in both sequences.
    let inner = hashes(collect(&left, scheme), 1070, 1930);
    let other = hashes(collect(&right, scheme), 2000, 3000);
    assert!(!inner.is_empty());
    assert!(inner.iter().all(|hash| other.contains(hash)));
}

#[test]
fn low_entropy_skipped() {
    let low = scheme(21, 10, 1.0);
    let ordinary = |sequence: &[u8], scheme| {
        collect(sequence, scheme)
            .into_iter()
            .filter(|seed| {
                let kmer = &sequence[seed.position as usize..seed.position as usize + 21];
                let (forward, reverse) = encode(kmer);
                seed.hash == jamhash::jamhash_u64(forward.min(reverse))
            })
            .count()
    };
    assert_eq!(ordinary(&[b'A'; 500], low), 0);
    // An odd k-mer of a dinucleotide repeat holds 11 and 10 of two bases, below 1 bit.
    let repeat = b"AT".repeat(100);
    assert_eq!(ordinary(&repeat, low), 0);
    assert!(ordinary(&repeat, scheme(21, 10, 0.9)) > 0);
    assert!(ordinary(&[b'A'; 500], scheme(21, 10, 0.0)) > 0);
}

#[test]
fn scheme_ranges_checked() {
    assert!(scheme(0, 10, 1.0).validate().is_err());
    assert!(scheme(33, 10, 1.0).validate().is_err());
    assert!(scheme(21, 0, 1.0).validate().is_err());
    assert!(scheme(21, 10, 2.5).validate().is_err());
    assert!(scheme(21, 10, f64::NAN).validate().is_err());
    assert!(scheme(21, 10, 1.0).validate().is_ok());
}

#[test]
fn salted_run_seeds() {
    let scheme = scheme(21, 100, 1.0);
    let flank = random_bases(&mut rng(3), 500);
    let with_run = |length: usize| {
        let mut sequence = flank.clone();
        sequence.extend(std::iter::repeat_n(b'A', length));
        sequence.extend(random_bases(&mut rng(4), 500));
        sequence
    };
    let salted = |sequence: &[u8]| -> Vec<Kmer> {
        collect(sequence, scheme)
            .into_iter()
            .filter(|seed| seed.hash == jamhash::jamhash_u64(0) ^ LOW_COMPLEXITY_SALT)
            .collect()
    };
    // A run shorter than a minimizer window still gets one salted seed.
    assert!(!salted(&with_run(50)).is_empty());
    let sequence = with_run(1200);
    let seeds = salted(&sequence);
    assert!(seeds.len() >= 12);
    assert!(
        seeds
            .iter()
            .all(|seed| (480..500 + 1200).contains(&(seed.position as usize)))
    );
    // The ordinary seeds of the same sequence never include that k-mer.
    assert!(
        collect(&sequence, scheme)
            .iter()
            .all(|seed| seed.hash != jamhash::jamhash_u64(0))
    );
    // The query side salts every k-mer below the limit, so it finds the run's seed.
    let mut query = Vec::new();
    Seeder::new(scheme)
        .unwrap()
        .query_kmers(&[b'A'; 1200], |kmer| query.push(kmer.hash));
    assert_eq!(query.len(), 1200 - 20);
    assert!(query.iter().all(|hash| *hash == seeds[0].hash));
}

#[test]
fn salted_seed_minima() {
    // Random sequence, homopolymers, skewed two-base repeats and breaks, against every k-mer's
    // entropy computed alone.
    let mut sequence = random_bases(&mut rng(7), 3000);
    sequence.extend(std::iter::repeat_n(b'T', 400));
    sequence.extend(b"AAAT".repeat(150));
    sequence.extend(random_bases(&mut rng(8), 50));
    sequence.push(b'N');
    sequence.extend(b"CCCG".repeat(60));
    sequence.extend(random_bases(&mut rng(9), 2000));
    for (k, min_entropy) in [(21u8, 1.0), (25, 1.0), (15, 1.2)] {
        let k_bases = usize::from(k);
        let mut expected = Vec::new();
        for (offset, run) in clean_runs(&sequence).filter(|(_, run)| run.len() >= k_bases) {
            let kmers: Vec<(bool, u64)> = run
                .windows(k_bases)
                .map(|kmer| {
                    let (forward, reverse) = encode(kmer);
                    (
                        entropy_of(kmer) < min_entropy,
                        jamhash::jamhash_u64(forward.min(reverse)),
                    )
                })
                .collect();
            for group in kmers.chunk_by(|a, b| a.0 == b.0).filter(|group| group[0].0) {
                let start = (group.as_ptr() as usize - kmers.as_ptr() as usize)
                    / std::mem::size_of::<(bool, u64)>();
                let hashes: Vec<u64> = group.iter().map(|(_, hash)| *hash).collect();
                sliding_minima(&hashes, offset + start, 100, &mut expected);
            }
        }
        let seeder = Seeder::new(scheme(k, 100, min_entropy)).unwrap();
        let mut found = Vec::new();
        seeder.seeds(&sequence, |seed| {
            let kmer = &sequence[seed.position as usize..seed.position as usize + k_bases];
            if !seeder.passes_entropy(encode(kmer).0) {
                found.push(seed.position);
            }
        });
        assert_eq!(found, expected, "k {k}");
        assert!(!found.is_empty());
    }
}

#[test]
fn query_kmers_complete() {
    let mut query = random_bases(&mut rng(11), 400);
    query[200] = b'N';
    let seeder = Seeder::new(scheme(15, 10, 1.0)).unwrap();
    let mut kmers = Vec::new();
    seeder.query_kmers(&query, |kmer| kmers.push(kmer));
    let expected: Vec<u64> = (0..=400 - 15)
        .filter(|start| !(*start..start + 15).contains(&200))
        .map(|start| start as u64)
        .collect();
    assert_eq!(
        kmers.iter().map(|kmer| kmer.position).collect::<Vec<_>>(),
        expected
    );
    for kmer in &kmers {
        let start = kmer.position as usize;
        let (forward, reverse) = encode(&query[start..start + 15]);
        assert_eq!(kmer.reverse, reverse < forward);
        assert_eq!(kmer.hash, jamhash::jamhash_u64(forward.min(reverse)));
    }
}

#[test]
fn neighbours_cover_substitutions() {
    let query = random_bases(&mut rng(12), 60);
    let seeder = Seeder::new(scheme(15, 10, 1.0)).unwrap();
    let mut neighbours = std::collections::HashSet::new();
    let mut count = 0;
    seeder.query_neighbours(&query, |kmer| {
        neighbours.insert((kmer.position, kmer.hash, kmer.reverse));
        count += 1;
    });
    assert_eq!(count, (60 - 15 + 1) * 15 * 3);
    for at in 0..60 {
        for base in *b"ACGT" {
            if base == query[at] {
                continue;
            }
            let mut changed = query.clone();
            changed[at] = base;
            seeder.query_kmers(&changed, |kmer| {
                let start = kmer.position as usize;
                if (start..start + 15).contains(&at) {
                    assert!(neighbours.contains(&(kmer.position, kmer.hash, kmer.reverse)));
                }
            });
        }
    }
}

#[test]
fn two_neighbours_cover() {
    let query = random_bases(&mut rng(14), 40);
    let seeder = Seeder::new(scheme(15, 10, 1.0)).unwrap();
    let mut neighbours = std::collections::HashSet::new();
    let mut count = 0;
    seeder.query_neighbours_two(&query, |kmer| {
        neighbours.insert((kmer.position, kmer.hash, kmer.reverse));
        count += 1;
    });
    // 9 * C(15, 2) per k-mer, all different.
    assert_eq!(count, (40 - 15 + 1) * 945);
    assert_eq!(neighbours.len(), count);
    for (first, second) in [(0, 1), (3, 17), (20, 39), (5, 6), (10, 24)] {
        let mut changed = query.clone();
        for at in [first, second] {
            changed[at] = if changed[at] == b'A' { b'C' } else { b'A' };
        }
        seeder.query_kmers(&changed, |kmer| {
            let start = kmer.position as usize;
            let range = start..start + 15;
            if range.contains(&first) && range.contains(&second) {
                assert!(neighbours.contains(&(kmer.position, kmer.hash, kmer.reverse)));
            }
        });
    }
}

#[test]
fn entropy_counts_match() {
    let sequence = random_bases(&mut rng(13), 2000);
    for k in [1u8, 7, 25, 32] {
        let seeder = Seeder::new(scheme(k, 10, 1.3)).unwrap();
        for kmer in sequence.windows(usize::from(k)).step_by(7) {
            assert_eq!(
                seeder.passes_entropy(encode(kmer).0),
                entropy_of(kmer) >= 1.3,
                "k {k}"
            );
        }
    }
}

#[test]
fn base_helpers_agree() {
    assert_eq!(base_code(b'g'), Some(2));
    assert_eq!(base_code(b'N'), None);
    assert_eq!(complement(b'a'), b'T');
    assert_eq!(reverse_complement(b"AACGTNr"), b"YNACGTT");
    let mut found = Vec::new();
    kmers(b"ACGTNacgta", 3, |code| found.push(code));
    let starts: Vec<usize> = found.iter().map(|code| code.position).collect();
    assert_eq!(starts, [0, 1, 5, 6, 7]);
    assert_eq!(
        (found[0].forward, found[0].reverse),
        (0b00_01_10, 0b01_10_11)
    );
    assert_eq!(found[0].canonical(), 0b00_01_10);
}

#[test]
fn repeats_masked() {
    let mut random = rng(0xd05);
    for period in [&b"A"[..], b"AT", b"ACG", b"AAGT", b"ACGTTG"] {
        let mut repeat: Vec<u8> = period.iter().copied().cycle().take(2_000).collect();
        // Scattered substitutions, about one in 60 bases, leave the repeat masked.
        for at in (30..2_000).step_by(60) {
            repeat[at] = if repeat[at] == b'C' { b'G' } else { b'C' };
        }
        assert_eq!(
            low_complexity_mask(&repeat),
            vec![0..2_000],
            "period {period:?}"
        );
        let flanked = [
            random_bases(&mut random, 500),
            repeat[..300].to_vec(),
            random_bases(&mut random, 500),
        ]
        .concat();
        let masked = low_complexity_mask(&flanked);
        assert_eq!(masked.len(), 1, "period {period:?}");
        assert!(masked[0].start <= 520 && masked[0].end >= 780, "{masked:?}");
        assert!(masked[0].start >= 440 && masked[0].end <= 860, "{masked:?}");
    }
}

#[test]
fn ordinary_unmasked() {
    let mut random = rng(0xd06);
    let bases = random_bases(&mut random, 100_000);
    assert!(low_complexity_mask(&bases).is_empty());
    // A short run inside random bases and unknown bases do not make a window simple.
    let run = [
        &bases[..500],
        &[b'A'; 12][..],
        &[b'N'; 40][..],
        &bases[500..1_000],
    ]
    .concat();
    assert!(low_complexity_mask(&run).is_empty());
    assert!(low_complexity_mask(b"AA").is_empty());
    assert_eq!(low_complexity_mask(&[b'T'; 30]), vec![0..30]);
}

#[test]
fn indel_neighbours_cover() {
    let query = random_bases(&mut rng(15), 50);
    let seeder = Seeder::new(scheme(15, 10, 1.0)).unwrap();
    let collect = |substituted: bool| {
        let mut found = std::collections::HashSet::new();
        let mut count = 0;
        seeder.query_neighbours_indel(&query, substituted, |kmer| {
            found.insert((kmer.position, kmer.hash, kmer.reverse));
            count += 1;
        });
        (found, count)
    };
    let (plain, count) = collect(false);
    // 13 places times 4 bases for an added base, 14 places for a removed one.
    assert_eq!(count, (50 - 14 + 1) * 52 + (50 - 16 + 1) * 14);
    let (changed, more) = collect(true);
    assert_eq!(more, count * 46);
    // A seed holding an added or a removed query base inside it, and one substitution.
    for at in [1, 9, 20, 33, 48] {
        let added = [&query[..at], b"G", &query[at..]].concat();
        let removed = [&query[..at], &query[at + 1..]].concat();
        for (copy, end) in [(added, 14), (removed, 15)] {
            let mut substituted = copy.clone();
            let other = (at + 5) % copy.len();
            substituted[other] = if copy[other] == b'A' { b'C' } else { b'A' };
            seeder.query_kmers(&copy, |kmer| {
                let start = kmer.position as usize;
                if start < at && at < start + end {
                    let key = (kmer.position, kmer.hash, kmer.reverse);
                    assert!(plain.contains(&key), "at {at}");
                }
            });
            seeder.query_kmers(&substituted, |kmer| {
                let start = kmer.position as usize;
                if start < at && at < start + end && (start..start + 15).contains(&other) {
                    assert!(changed.contains(&(kmer.position, kmer.hash, kmer.reverse)));
                }
            });
        }
    }
}

/// A k-mer of a clean run and whether it reaches the entropy limit, for `queue_seeds`.
#[derive(Clone, Copy)]
struct Candidate {
    kmer: Kmer,
    passes: bool,
}

/// The sliding minimum as a queue of rising hashes, for `queue_seeds`.
struct QueueMinimum {
    w: usize,
    queue: std::collections::VecDeque<Candidate>,
    kmers: usize,
    last: Option<u64>,
}

impl QueueMinimum {
    fn new(w: usize) -> Self {
        let queue = std::collections::VecDeque::new();
        let (kmers, last) = (0, None);
        Self {
            w,
            queue,
            kmers,
            last,
        }
    }

    fn push(&mut self, candidate: Candidate) -> Option<Candidate> {
        while (self.queue.back()).is_some_and(|back| back.kmer.hash > candidate.kmer.hash) {
            self.queue.pop_back();
        }
        self.queue.push_back(candidate);
        self.kmers += 1;
        if self.kmers < self.w {
            return None;
        }
        let first_start = candidate.kmer.position + 1 - self.w as u64;
        while (self.queue.front()).is_some_and(|front| front.kmer.position < first_start) {
            self.queue.pop_front();
        }
        let best = *self.queue.front()?;
        if self.last == Some(best.kmer.position) {
            return None;
        }
        self.last = Some(best.kmer.position);
        Some(best)
    }

    fn finish(&mut self) -> Option<Candidate> {
        let short = (self.kmers < self.w)
            .then(|| self.queue.front().copied())
            .flatten();
        self.queue.clear();
        (self.kmers, self.last) = (0, None);
        short
    }
}

/// The seeds of a contig as the first seed selection chose them, k-mer by k-mer with two
/// queues: the reference that `Seeder::seeds` must match exactly, in order.
fn queue_seeds(seeder: &Seeder, contig: &[u8]) -> Vec<Kmer> {
    let mut seeds = Vec::new();
    let w = usize::from(seeder.scheme.w);
    let (mut all, mut low) = (QueueMinimum::new(w), QueueMinimum::new(w));
    let mut low_seeds = Vec::new();
    let mut emit_all = |selected: Option<Candidate>| {
        if let Some(selected) = selected.filter(|selected| selected.passes) {
            seeds.push(selected.kmer);
        }
    };
    let mut keep_low = |selected: Option<Candidate>| {
        if let Some(selected) = selected {
            let hash = selected.kmer.hash ^ LOW_COMPLEXITY_SALT;
            low_seeds.push(Kmer {
                hash,
                ..selected.kmer
            });
        }
    };
    scan(contig, seeder.scheme.k, |step| match step {
        Step::Break => {
            emit_all(all.finish());
            keep_low(low.finish());
        }
        Step::Kmer(code) => {
            let candidate = Candidate {
                kmer: Kmer {
                    hash: jamhash::jamhash_u64(code.canonical()),
                    position: code.position as u64,
                    reverse: code.reverse < code.forward,
                },
                passes: seeder.passes_entropy(code.canonical()),
            };
            emit_all(all.push(candidate));
            if candidate.passes {
                keep_low(low.finish());
            } else {
                keep_low(low.push(candidate));
            }
        }
    });
    emit_all(all.finish());
    keep_low(low.finish());
    seeds.extend(low_seeds);
    seeds
}

/// Random contigs with breaks, lower case, IUPAC codes, homopolymers, short tandem repeats and
/// two-base stretches, for comparing seed selections.
fn mixed_contig(draws: &mut rand::rngs::StdRng, length: usize) -> Vec<u8> {
    let mut contig = Vec::with_capacity(length + 600);
    while contig.len() < length {
        match rand::RngExt::random_range(draws, 0..12) {
            0 => contig.extend(std::iter::repeat_n(
                b"ACGT"[contig.len() % 4],
                1 + contig.len() % 300,
            )),
            1 => {
                let unit_length = rand::RngExt::random_range(draws, 1..7);
                let unit = random_bases(draws, unit_length);
                let times = rand::RngExt::random_range(draws, 1..80);
                contig.extend(unit.repeat(times));
            }
            2 => {
                let length = rand::RngExt::random_range(draws, 1..400);
                for _ in 0..length {
                    contig.push(b"AT"[rand::RngExt::random_range(draws, 0..2)]);
                }
            }
            3 => contig.push(b"NRYnx-"[rand::RngExt::random_range(draws, 0..6)]),
            4 => contig.extend(random_bases(draws, 40).to_ascii_lowercase()),
            _ => {
                let length = rand::RngExt::random_range(draws, 1..600);
                contig.extend(random_bases(draws, length));
            }
        }
    }
    contig
}

/// `Seeder::seeds` gives exactly the seeds of the queue selection, in the same order, for many
/// schemes and contigs, including contigs shorter than one k-mer or one window.
#[test]
fn seeds_match_queue_selection() {
    let mut draws = rng(0xc0ffee);
    let schemes = [
        (1, 1, 0.0),
        (3, 2, 0.5),
        (15, 10, 1.2),
        (21, 100, 1.0),
        (25, 100, 1.0),
        (25, 100, 1.9),
        (25, 1, 1.0),
        (31, 1000, 1.5),
        (32, 7, 2.0),
        (16, 50, 1.0),
    ];
    let (mut seeds, mut salted) = (0, 0);
    for (k, w, min_entropy) in schemes {
        let seeder = Seeder::new(scheme(k, w, min_entropy)).unwrap();
        for length in [0, 1, 5, 24, 25, 26, 130, 3000, 40_000] {
            let contig = mixed_contig(&mut draws, length);
            let mut found = Vec::new();
            seeder.seeds(&contig, |seed| found.push(seed));
            let expected = queue_seeds(&seeder, &contig);
            assert_eq!(found, expected, "k {k} w {w} length {length}");
            seeds += found.len();
            let low = |seed: &&Kmer| {
                let at = seed.position as usize;
                !seeder.passes_entropy(encode(&contig[at..at + usize::from(k)]).0)
            };
            salted += found.iter().filter(low).count();
        }
    }
    assert!(
        seeds > 20_000 && salted > 1_000,
        "{seeds} seeds, {salted} salted"
    );
}
