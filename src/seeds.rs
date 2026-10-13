//! Canonical k-mers, their jamhash values, and the seeds of a contig: the smallest k-mer of each
//! minimizer window, with the low-complexity filter.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

//! A k-mer below the entropy limit is never an ordinary seed; its runs get salted seeds of their
//! own, and the query side salts such query k-mers too, so both sides meet.

use serde::{Deserialize, Serialize};
use std::ops::Range;
use thiserror::Error;

/// Errors of the seeds module.
#[derive(Clone, Debug, Error, PartialEq)]
pub enum SeedError {
    /// The seed scheme has a value outside its range.
    #[error("invalid seed scheme: {0}")]
    Scheme(&'static str),
}

/// How seeds are chosen: k-mer size, minimizer window and entropy limit.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SeedScheme {
    /// Bases per k-mer, 1 to 32.
    pub k: u8,
    /// K-mers per minimizer window (w), at least 1.
    pub w: u16,
    /// The lowest Shannon entropy (bits, 0 to 2) of an ordinary seed.
    pub min_entropy: f64,
}

impl Default for SeedScheme {
    fn default() -> Self {
        Self {
            k: 25,
            w: 100,
            min_entropy: 1.0,
        }
    }
}

impl SeedScheme {
    /// Checks that every value lies in its range.
    pub fn validate(&self) -> Result<(), SeedError> {
        if !(1..=32).contains(&self.k) {
            return Err(SeedError::Scheme("k must be 1 to 32"));
        }
        if self.w == 0 {
            return Err(SeedError::Scheme(
                "w (k-mers per minimizer window) must be at least 1",
            ));
        }
        if !(0.0..=2.0).contains(&self.min_entropy) {
            return Err(SeedError::Scheme("the entropy limit must be 0 to 2"));
        }
        Ok(())
    }
}

/// A hashed canonical k-mer at a position of a sequence.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Kmer {
    /// jamhash of the canonical k-mer, salted for k-mers below the entropy limit.
    pub hash: u64,
    /// Start of the k-mer in the sequence.
    pub position: u64,
    /// The canonical k-mer is the reverse complement of the k-mer read on the sequence.
    pub reverse: bool,
}

/// Chooses seeds and query k-mers with one seed scheme.
pub struct Seeder {
    scheme: SeedScheme,
    /// Whether a k-mer passes the entropy limit, by its count of A, C and G (T is the rest).
    passes: Vec<bool>,
}

impl Seeder {
    /// Prepares the entropy table of a valid scheme.
    pub fn new(scheme: SeedScheme) -> Result<Self, SeedError> {
        scheme.validate()?;
        let k = usize::from(scheme.k);
        let side = k + 1;
        let mut passes = vec![false; side * side * side];
        for a in 0..=k {
            for c in 0..=k - a {
                for g in 0..=k - a - c {
                    let entropy = shannon_entropy([a, c, g, k - a - c - g], k);
                    passes[(a * side + c) * side + g] = entropy >= scheme.min_entropy;
                }
            }
        }
        Ok(Self { scheme, passes })
    }

    /// The scheme this seeder uses.
    pub fn scheme(&self) -> SeedScheme {
        self.scheme
    }

    /// Whether a 2-bit k-mer reaches the entropy limit.
    fn passes_entropy(&self, kmer: u64) -> bool {
        let k = u32::from(self.scheme.k);
        let mask = kmer_mask(k) & 0x5555_5555_5555_5555;
        let (high, low) = ((kmer >> 1) & mask, kmer & mask);
        let t = (high & low).count_ones();
        let g = (high & !low).count_ones();
        let c = (!high & low & mask).count_ones();
        let a = k - t - g - c;
        let side = k as usize + 1;
        self.passes[(a as usize * side + c as usize) * side + g as usize]
    }

    /// Calls `emit` for every seed of a contig: the smallest hash (leftmost on ties) of each
    /// minimizer window of a clean run, or of a whole clean run shorter than one. Selected k-mers
    /// below the entropy limit give no seed; their runs select salted seeds of their own.
    pub fn seeds(&self, contig: &[u8], mut emit: impl FnMut(Kmer)) {
        let (k, w) = (usize::from(self.scheme.k), usize::from(self.scheme.w));
        let mut run = RunKmers::default();
        let mut low_seeds = Vec::new();
        let mut start = 0;
        let breaks = (0..contig.len()).filter(|&at| CODE[usize::from(contig[at])] == OTHER);
        for end in breaks.chain([contig.len()]) {
            if end - start >= k {
                self.run_kmers(&contig[start..end], start, &mut run);
                minimizers(&run.hashes, w, |at| {
                    if run.passes(at) {
                        emit(run.kmer(at));
                    }
                });
                let mut at = 0;
                while at < run.hashes.len() {
                    let first = at;
                    while at < run.hashes.len() && !run.passes(at) {
                        at += 1;
                    }
                    minimizers(&run.hashes[first..at], w, |low| {
                        let kmer = run.kmer(first + low);
                        let hash = kmer.hash ^ LOW_COMPLEXITY_SALT;
                        low_seeds.push(Kmer { hash, ..kmer });
                    });
                    at += 1;
                }
            }
            start = end + 1;
        }
        low_seeds.into_iter().for_each(emit);
    }

    /// Fills `kmers` with the hash, strand and entropy test of every k-mer of a clean run that
    /// starts at `offset` of its contig; the base counts of each k-mer follow it along the run.
    fn run_kmers(&self, run: &[u8], offset: usize, kmers: &mut RunKmers) {
        let k = usize::from(self.scheme.k);
        let mask = kmer_mask(k as u32);
        let side = k + 1;
        let (mut forward, mut reverse, mut counts) = (0u64, 0u64, [0usize; 4]);
        kmers.hashes.clear();
        kmers.flags.clear();
        kmers.offset = offset;
        kmers.k = k;
        for (index, &byte) in run.iter().enumerate() {
            let code = CODE[usize::from(byte)];
            counts[usize::from(code)] += 1;
            if index >= k {
                counts[usize::from(CODE[usize::from(run[index - k])])] -= 1;
            }
            let code = u64::from(code);
            forward = ((forward << 2) | code) & mask;
            reverse = (reverse >> 2) | ((3 - code) << (2 * (k - 1)));
            if index + 1 < k {
                continue;
            }
            // The canonical k-mer is the reverse complement when it is the smaller: its A, C and
            // G are the T, G and C of the k-mer read on the run.
            let backward = reverse < forward;
            let [a, c, g, t] = counts;
            let counted = if backward { (t, g, c) } else { (a, c, g) };
            let passes = self.passes[(counted.0 * side + counted.1) * side + counted.2];
            kmers
                .hashes
                .push(jamhash::jamhash_u64(forward.min(reverse)));
            kmers
                .flags
                .push(u8::from(backward) | (u8::from(passes) << 1));
        }
    }

    /// Calls `emit` for every k-mer of a query, in order. A k-mer below the entropy limit gets the
    /// salted hash, as the seeds of low-complexity runs do.
    pub fn query_kmers(&self, query: &[u8], mut emit: impl FnMut(Kmer)) {
        kmers(query, self.scheme.k, |code| emit(self.query_kmer(code)));
    }

    /// Calls `emit` for every k-mer with one substitution against a k-mer of the query (3 k per
    /// query k-mer), at the start of the query k-mer, so a seed with one mismatch is found too.
    pub fn query_neighbours(&self, query: &[u8], mut emit: impl FnMut(Kmer)) {
        let k = usize::from(self.scheme.k);
        kmers(query, self.scheme.k, |code| {
            for base in 0..k {
                for change in 1..4u64 {
                    emit(self.query_kmer(KmerCode {
                        forward: code.forward ^ (change << (2 * (k - 1 - base))),
                        reverse: code.reverse ^ (change << (2 * base)),
                        ..code
                    }));
                }
            }
        });
    }

    /// Calls `emit` for every k-mer with two substitutions against a k-mer of the query (9 k (k - 1)
    /// / 2 per query k-mer), at the start of the query k-mer, so a seed with two mismatches is found.
    pub fn query_neighbours_two(&self, query: &[u8], mut emit: impl FnMut(Kmer)) {
        let k = usize::from(self.scheme.k);
        kmers(query, self.scheme.k, |code| {
            for first in 0..k {
                for second in first + 1..k {
                    for a in 1..4u64 {
                        for b in 1..4u64 {
                            let forward =
                                (a << (2 * (k - 1 - first))) ^ (b << (2 * (k - 1 - second)));
                            let reverse = (a << (2 * first)) ^ (b << (2 * second));
                            emit(self.query_kmer(KmerCode {
                                forward: code.forward ^ forward,
                                reverse: code.reverse ^ reverse,
                                ..code
                            }));
                        }
                    }
                }
            }
        });
    }

    /// Calls `emit` for every k-mer one indel away from the query at a start of the query: a
    /// base inserted into, or removed from, the inside of k bases. With `substituted`, also every
    /// k-mer one substitution away from those.
    pub fn query_neighbours_indel(&self, query: &[u8], substituted: bool, emit: impl FnMut(Kmer)) {
        self.indel_neighbours(query, 0..query.len(), substituted, emit);
    }

    /// As `query_neighbours_indel`, for the first start of `query` only.
    pub fn query_neighbours_indel_at(
        &self,
        query: &[u8],
        substituted: bool,
        emit: impl FnMut(Kmer),
    ) {
        self.indel_neighbours(query, 0..1, substituted, emit);
    }

    /// The k-mers one indel away (and with `substituted` one substitution more) at `starts`.
    fn indel_neighbours(
        &self,
        query: &[u8],
        starts: Range<usize>,
        substituted: bool,
        mut emit: impl FnMut(Kmer),
    ) {
        let k = usize::from(self.scheme.k);
        let codes: Vec<u8> = query.iter().map(|&byte| CODE[usize::from(byte)]).collect();
        let clean = |range: Range<usize>| {
            range.end <= codes.len() && codes[range].iter().all(|&code| code != OTHER)
        };
        let mut push = |position: usize, bases: &mut dyn Iterator<Item = u8>| {
            let forward = bases.fold(0u64, |word, code| (word << 2) | u64::from(code));
            let reverse = reverse_code(forward, k);
            let code = KmerCode {
                position,
                forward,
                reverse,
            };
            emit(self.query_kmer(code));
            if substituted {
                for base in 0..k {
                    for change in 1..4u64 {
                        emit(self.query_kmer(KmerCode {
                            forward: forward ^ (change << (2 * (k - 1 - base))),
                            reverse: reverse ^ (change << (2 * base)),
                            position,
                        }));
                    }
                }
            }
        };
        for start in starts.start..starts.end.min(codes.len()) {
            if clean(start..start + k - 1) {
                let kept = &codes[start..start + k - 1];
                for at in 1..k - 1 {
                    for added in 0..4u8 {
                        let mut bases = kept[..at]
                            .iter()
                            .copied()
                            .chain([added])
                            .chain(kept[at..].iter().copied());
                        push(start, &mut bases);
                    }
                }
            }
            if clean(start..start + k + 1) {
                let longer = &codes[start..start + k + 1];
                for at in 1..k {
                    let mut bases = longer[..at].iter().chain(&longer[at + 1..]).copied();
                    push(start, &mut bases);
                }
            }
        }
    }

    /// The hash and strand of one query k-mer, salted below the entropy limit.
    fn query_kmer(&self, code: KmerCode) -> Kmer {
        let mut hash = jamhash::jamhash_u64(code.canonical());
        if !self.passes_entropy(code.canonical()) {
            hash ^= LOW_COMPLEXITY_SALT;
        }
        Kmer {
            hash,
            position: code.position as u64,
            reverse: code.reverse < code.forward,
        }
    }
}

/// A k-mer of a sequence in 2-bit code (A 0, C 1, G 2, T 3), read on both strands.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct KmerCode {
    /// Start of the k-mer in the sequence.
    pub position: usize,
    /// The k-mer as read on the sequence, its first base in the highest bits.
    pub forward: u64,
    /// The reverse complement of the k-mer in the same code.
    pub reverse: u64,
}

impl KmerCode {
    /// The canonical k-mer: the smaller of the k-mer and its reverse complement.
    pub fn canonical(self) -> u64 {
        self.forward.min(self.reverse)
    }
}

/// Calls `visit` for every k-mer of `sequence` that holds only A, C, G and T (either case).
/// `k` is 1 to 32.
pub fn kmers(sequence: &[u8], k: u8, mut visit: impl FnMut(KmerCode)) {
    scan(sequence, k, |step| {
        if let Step::Kmer(code) = step {
            visit(code);
        }
    });
}

/// The 2-bit code of a base (A 0, C 1, G 2, T 3, either case), `None` for any other byte.
pub fn base_code(byte: u8) -> Option<u8> {
    let code = CODE[usize::from(byte)];
    (code != OTHER).then_some(code)
}

/// The complement of an IUPAC base, in upper case; other bytes stay as they are.
pub fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    }
}

/// The reverse complement of `bases`.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&base| complement(base)).collect()
}

/// The low-complexity stretches of `sequence`, ascending and disjoint: the union of the DUST
/// windows whose triplet score passes `DUST_LIMIT`. A window holds `DUST_WINDOW` bases, or the
/// whole sequence when it is shorter; triplets with other bytes than A, C, G and T add no pairs.
pub fn low_complexity_mask(sequence: &[u8]) -> Vec<Range<usize>> {
    let window = DUST_WINDOW.min(sequence.len());
    let mut masked: Vec<Range<usize>> = Vec::new();
    if window < 3 {
        return masked;
    }
    let triplet = |at: usize| -> Option<usize> {
        let mut codes = sequence[at..at + 3].iter().map(|byte| base_code(*byte));
        codes.try_fold(0, |value, code| Some(value * 4 + usize::from(code?)))
    };
    let per_window = window - 2;
    let pairs_limit = DUST_LIMIT * (per_window as u64 - 1);
    let (mut counts, mut score) = ([0u64; 64], 0u64);
    for end in 0..sequence.len() - 2 {
        if let Some(code) = triplet(end) {
            score += counts[code];
            counts[code] += 1;
        }
        if end >= per_window
            && let Some(code) = triplet(end - per_window)
        {
            counts[code] -= 1;
            score -= counts[code];
        }
        // The score pairs equal triplets: sum of c (c - 1) / 2 over the triplet counts c.
        if end + 1 >= per_window && score * 10 > pairs_limit {
            let start = end + 1 - per_window;
            match masked.last_mut() {
                Some(last) if last.end >= start => last.end = start + window,
                _ => masked.push(start..start + window),
            }
        }
    }
    masked
}

/// Bases in one DUST window, as in BLAST and minimap2 (sdust).
const DUST_WINDOW: usize = 64;
/// The DUST score limit, ten times the mean pairs of equal triplets per triplet (sdust `-t`).
const DUST_LIMIT: u64 = 20;

/// Mixed into the hash of a seed of a low-complexity run, so it differs from the hash of the same
/// k-mer as an ordinary seed.
const LOW_COMPLEXITY_SALT: u64 = 0x6c6f_7720_636f_6d70;

/// The code of every byte that is not A, C, G or T.
const OTHER: u8 = 4;

/// The 2-bit code of every byte, `OTHER` for bytes other than A, C, G and T.
const CODE: [u8; 256] = {
    let mut code = [OTHER; 256];
    let mut index = 0;
    while index < 4 {
        code[b"ACGT"[index] as usize] = index as u8;
        code[b"acgt"[index] as usize] = index as u8;
        index += 1;
    }
    code
};

/// Shannon entropy in bits of a k-mer with the given base counts.
fn shannon_entropy(counts: [usize; 4], k: usize) -> f64 {
    counts
        .iter()
        .filter(|count| **count > 0)
        .map(|count| {
            let share = *count as f64 / k as f64;
            -share * share.log2()
        })
        .sum()
}

/// The low `2 * k` bits, which hold a 2-bit k-mer.
fn kmer_mask(k: u32) -> u64 {
    if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

/// The reverse complement of a 2-bit k-mer of `k` bases.
fn reverse_code(forward: u64, k: usize) -> u64 {
    let mut word = !forward;
    word = ((word >> 2) & 0x3333_3333_3333_3333) | ((word & 0x3333_3333_3333_3333) << 2);
    word = ((word >> 4) & 0x0f0f_0f0f_0f0f_0f0f) | ((word & 0x0f0f_0f0f_0f0f_0f0f) << 4);
    word.swap_bytes() >> (64 - 2 * k)
}

/// One step of a scan: a k-mer, or a byte other than A, C, G or T that ends a clean run.
enum Step {
    Kmer(KmerCode),
    Break,
}

/// Calls `visit` for every k-mer of `sequence` and for every byte that breaks a run.
fn scan(sequence: &[u8], k: u8, mut visit: impl FnMut(Step)) {
    let k = u32::from(k);
    let mask = kmer_mask(k);
    let (mut forward, mut reverse, mut bases) = (0u64, 0u64, 0u32);
    for (index, &byte) in sequence.iter().enumerate() {
        let code = CODE[usize::from(byte)];
        if code == OTHER {
            bases = 0;
            visit(Step::Break);
            continue;
        }
        let code = u64::from(code);
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((3 - code) << (2 * (k - 1)));
        bases = (bases + 1).min(k);
        if bases == k {
            visit(Step::Kmer(KmerCode {
                position: index + 1 - k as usize,
                forward,
                reverse,
            }));
        }
    }
}

/// The k-mers of one clean run, in order: their hashes, and per k-mer whether the canonical
/// k-mer is the reverse complement (bit 0) and whether it reaches the entropy limit (bit 1).
#[derive(Default)]
struct RunKmers {
    hashes: Vec<u64>,
    flags: Vec<u8>,
    /// Start of the run in its contig.
    offset: usize,
    k: usize,
}

impl RunKmers {
    fn passes(&self, at: usize) -> bool {
        self.flags[at] & 2 != 0
    }

    /// The k-mer at index `at` of the run.
    fn kmer(&self, at: usize) -> Kmer {
        Kmer {
            hash: self.hashes[at],
            position: (self.offset + at) as u64,
            reverse: self.flags[at] & 1 != 0,
        }
    }
}

/// Calls `select` with the index of the smallest hash (leftmost on ties) of every window of `w`
/// hashes, once per index and in order; with fewer hashes than `w`, of all of them. The minimum
/// is searched for again only when it leaves the window.
fn minimizers(hashes: &[u64], w: usize, mut select: impl FnMut(usize)) {
    if hashes.is_empty() {
        return;
    }
    let mut best = leftmost_minimum(&hashes[..hashes.len().min(w)]);
    select(best);
    for end in w..hashes.len() {
        if best + w <= end {
            let start = end + 1 - w;
            best = start + leftmost_minimum(&hashes[start..=end]);
            select(best);
        } else if hashes[end] < hashes[best] {
            best = end;
            select(best);
        }
    }
}

/// The index of the smallest of some hashes, leftmost on ties.
fn leftmost_minimum(hashes: &[u64]) -> usize {
    let (mut best, mut lowest) = (0, u64::MAX);
    for (index, &hash) in hashes.iter().enumerate() {
        if index == 0 || hash < lowest {
            (best, lowest) = (index, hash);
        }
    }
    best
}

#[cfg(test)]
#[path = "seeds_tests.rs"]
mod tests;
