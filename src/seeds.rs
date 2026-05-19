//! Canonical k-mers, their jamhash values, and the seeds of a contig: the smallest k-mer of each
//! minimizer window, with the low-complexity filter.
//!
//! The entry side keeps seeds only; the query side looks up every query k-mer. A k-mer below the
//! entropy limit never becomes an ordinary seed. Runs of such k-mers get seeds of their own,
//! chosen the same way inside the run, with a salted hash. The query side salts the hash of every
//! query k-mer below the limit, so both sides meet.

use std::collections::VecDeque;
use thiserror::Error;

/// Errors of the seeds module.
#[derive(Debug, Error, PartialEq)]
pub enum SeedError {
    /// The seed scheme has a value outside its range.
    #[error("invalid seed scheme: {0}")]
    Scheme(&'static str),
}

/// Mixed into the hash of a seed of a low-complexity run, so it differs from the hash of the same
/// k-mer as an ordinary seed.
pub const LOW_COMPLEXITY_SALT: u64 = 0x6c6f_7720_636f_6d70;

/// How seeds are chosen: k-mer size, minimizer window and entropy limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeedScheme {
    /// Bases per k-mer, 1 to 32.
    pub k: u8,
    /// K-mers per minimizer window, at least 1.
    pub window: u16,
    /// The lowest Shannon entropy (bits, 0 to 2) of an ordinary seed.
    pub min_entropy: f64,
}

impl Default for SeedScheme {
    fn default() -> Self {
        Self {
            k: 25,
            window: 100,
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
        if self.window == 0 {
            return Err(SeedError::Scheme("the minimizer window must be at least 1"));
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

/// The hash of a canonical 2-bit k-mer.
pub fn kmer_hash(canonical: u64) -> u64 {
    jamhash::jamhash_u64(canonical)
}

const BREAK: u8 = 4;

const CODE: [u8; 256] = {
    let mut code = [BREAK; 256];
    let mut index = 0;
    while index < 4 {
        code[b"ACGT"[index] as usize] = index as u8;
        code[b"acgt"[index] as usize] = index as u8;
        index += 1;
    }
    code
};

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
    pub fn passes_entropy(&self, kmer: u64) -> bool {
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

    /// Calls `emit` for every seed of a contig. Every minimizer window of a clean run (only A, C,
    /// G and T) selects its smallest hash, the leftmost on ties; a clean run shorter than a full
    /// window selects its smallest. A selected k-mer below the entropy limit gives no seed;
    /// runs of such k-mers select seeds of their own, with salted hashes.
    pub fn seeds(&self, contig: &[u8], mut emit: impl FnMut(Kmer)) {
        let window = usize::from(self.scheme.window);
        let mut all = WindowMinimum::new(window);
        let mut low = WindowMinimum::new(window);
        let mut emit_all = |selected: Option<Candidate>| {
            if let Some(selected) = selected.filter(|selected| selected.passes) {
                emit(selected.kmer);
            }
        };
        let mut low_seeds = Vec::new();
        let mut keep_low = |selected: Option<Candidate>| {
            if let Some(selected) = selected {
                low_seeds.push(Kmer {
                    hash: selected.kmer.hash ^ LOW_COMPLEXITY_SALT,
                    ..selected.kmer
                });
            }
        };
        scan(contig, self.scheme.k, |step| match step {
            Step::Break => {
                emit_all(all.finish());
                keep_low(low.finish());
            }
            Step::Kmer {
                position,
                canonical,
                reverse,
            } => {
                let candidate = Candidate {
                    kmer: Kmer {
                        hash: kmer_hash(canonical),
                        position,
                        reverse,
                    },
                    passes: self.passes_entropy(canonical),
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
        low_seeds.into_iter().for_each(emit);
    }

    /// Calls `emit` for every k-mer of a query, in order. A k-mer below the entropy limit gets the
    /// salted hash, as the seeds of low-complexity runs do.
    pub fn query_kmers(&self, query: &[u8], mut emit: impl FnMut(Kmer)) {
        scan(query, self.scheme.k, |step| {
            if let Step::Kmer {
                position,
                canonical,
                reverse,
            } = step
            {
                let mut hash = kmer_hash(canonical);
                if !self.passes_entropy(canonical) {
                    hash ^= LOW_COMPLEXITY_SALT;
                }
                emit(Kmer {
                    hash,
                    position,
                    reverse,
                });
            }
        });
    }
}

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

fn kmer_mask(k: u32) -> u64 {
    if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    }
}

/// One step of a scan: a k-mer, or a byte other than A, C, G or T that ends a clean run.
enum Step {
    Kmer {
        position: u64,
        canonical: u64,
        reverse: bool,
    },
    Break,
}

/// Calls `visit` for every canonical k-mer of `sequence` and for every byte that breaks a run.
fn scan(sequence: &[u8], k: u8, mut visit: impl FnMut(Step)) {
    let k = u32::from(k);
    let mask = kmer_mask(k);
    let (mut forward, mut reverse, mut bases) = (0u64, 0u64, 0u32);
    for (index, &byte) in sequence.iter().enumerate() {
        let code = CODE[byte as usize];
        if code == BREAK {
            bases = 0;
            visit(Step::Break);
            continue;
        }
        let code = u64::from(code);
        forward = ((forward << 2) | code) & mask;
        reverse = (reverse >> 2) | ((3 - code) << (2 * (k - 1)));
        bases = (bases + 1).min(k);
        if bases == k {
            visit(Step::Kmer {
                position: (index + 1 - k as usize) as u64,
                canonical: forward.min(reverse),
                reverse: reverse < forward,
            });
        }
    }
}

#[derive(Clone, Copy)]
struct Candidate {
    kmer: Kmer,
    passes: bool,
}

/// The sliding minimum over the k-mers of one clean run.
struct WindowMinimum {
    window: usize,
    queue: VecDeque<Candidate>,
    kmers: usize,
    last: Option<u64>,
}

impl WindowMinimum {
    fn new(window: usize) -> Self {
        Self {
            window,
            queue: VecDeque::with_capacity(window + 1),
            kmers: 0,
            last: None,
        }
    }

    /// Adds the next k-mer of the run, and returns the minimum of a full window when it is new.
    fn push(&mut self, candidate: Candidate) -> Option<Candidate> {
        while self
            .queue
            .back()
            .is_some_and(|back| back.kmer.hash > candidate.kmer.hash)
        {
            self.queue.pop_back();
        }
        self.queue.push_back(candidate);
        self.kmers += 1;
        if self.kmers < self.window {
            return None;
        }
        let window_start = candidate.kmer.position + 1 - self.window as u64;
        while self
            .queue
            .front()
            .is_some_and(|front| front.kmer.position < window_start)
        {
            self.queue.pop_front();
        }
        let best = *self.queue.front()?;
        if self.last == Some(best.kmer.position) {
            return None;
        }
        self.last = Some(best.kmer.position);
        Some(best)
    }

    /// Ends the run: a run shorter than a full window returns its minimum.
    fn finish(&mut self) -> Option<Candidate> {
        let short = (self.kmers < self.window)
            .then(|| self.queue.front().copied())
            .flatten();
        self.queue.clear();
        self.kmers = 0;
        self.last = None;
        short
    }
}


#[cfg(test)]
mod tests;
