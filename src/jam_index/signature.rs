//! Bounded contig-stratified k21 MinHash signature selection.

use super::manifest::ScreenSelectionPolicy;
use crate::jamhash_u64_v1;
use needletail::Sequence;
use std::collections::{BTreeMap, BTreeSet};
use thiserror::Error;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContigSignature {
    pub requested_budget: u32,
    pub eligible_kmers: u64,
    pub hashes: Vec<u64>,
}

#[derive(Clone, Debug)]
pub struct MetagenomeSignatureBuilder {
    policy: ScreenSelectionPolicy,
    whole: BottomK,
    whole_sources: BTreeMap<u64, u32>,
    union: BTreeSet<u64>,
    contig_count: u64,
    total_bases: u64,
}

impl MetagenomeSignatureBuilder {
    pub fn new(policy: ScreenSelectionPolicy) -> Result<Self, SignatureSelectionError> {
        policy
            .validate()
            .map_err(|error| SignatureSelectionError::InvalidPolicy(error.to_string()))?;
        Ok(Self {
            whole: BottomK::new(policy.whole_metagenome_budget as usize)?,
            whole_sources: BTreeMap::new(),
            policy,
            union: BTreeSet::new(),
            contig_count: 0,
            total_bases: 0,
        })
    }

    pub fn add_contig(
        &mut self,
        sequence: &[u8],
    ) -> Result<ContigSignature, SignatureSelectionError> {
        let length =
            u64::try_from(sequence.len()).map_err(|_| SignatureSelectionError::Overflow)?;
        let requested_budget = self.policy.contig_budget.budget_for_bases(length);
        let contig_id =
            u32::try_from(self.contig_count).map_err(|_| SignatureSelectionError::Overflow)?;
        let mut contig = BottomK::new(requested_budget as usize)?;
        let mut eligible_kmers = 0u64;
        let normalized = sequence.normalize(false);
        for (_, kmer, _) in normalized.bit_kmers(self.policy.k, true) {
            let hash = jamhash_u64_v1(kmer.0);
            if hash == 0 {
                continue;
            }
            eligible_kmers = eligible_kmers.saturating_add(1);
            contig.insert(hash);
            if let BottomKUpdate::Retained { evicted } = self.whole.insert(hash) {
                if let Some(evicted) = evicted {
                    self.whole_sources.remove(&evicted);
                }
                self.whole_sources.entry(hash).or_insert(contig_id);
            }
        }
        let hashes = contig.into_sorted();
        self.union.extend(hashes.iter().copied());
        self.contig_count = self.contig_count.saturating_add(1);
        self.total_bases = self.total_bases.saturating_add(length);
        Ok(ContigSignature {
            requested_budget,
            eligible_kmers,
            hashes,
        })
    }

    #[must_use]
    pub fn finish(mut self) -> MetagenomeSignatures {
        let whole_metagenome_hashes = self.whole.into_sorted();
        self.union.extend(whole_metagenome_hashes.iter().copied());
        MetagenomeSignatures {
            contig_count: self.contig_count,
            total_bases: self.total_bases,
            whole_metagenome_hashes,
            whole_hash_contigs: self.whole_sources.into_iter().collect(),
            union_hashes: self.union.into_iter().collect(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeSignatures {
    pub contig_count: u64,
    pub total_bases: u64,
    pub whole_metagenome_hashes: Vec<u64>,
    pub whole_hash_contigs: Vec<(u64, u32)>,
    pub union_hashes: Vec<u64>,
}

#[derive(Clone, Debug)]
struct BottomK {
    capacity: usize,
    hashes: BTreeSet<u64>,
}

impl BottomK {
    fn new(capacity: usize) -> Result<Self, SignatureSelectionError> {
        if capacity == 0 {
            return Err(SignatureSelectionError::ZeroBudget);
        }
        Ok(Self {
            capacity,
            hashes: BTreeSet::new(),
        })
    }

    fn insert(&mut self, hash: u64) -> BottomKUpdate {
        if hash == 0 || self.hashes.contains(&hash) {
            return BottomKUpdate::Unchanged;
        }
        if self.hashes.len() < self.capacity {
            self.hashes.insert(hash);
            return BottomKUpdate::Retained { evicted: None };
        }
        let Some(&largest) = self.hashes.last() else {
            return BottomKUpdate::Unchanged;
        };
        if hash < largest {
            self.hashes.remove(&largest);
            self.hashes.insert(hash);
            return BottomKUpdate::Retained {
                evicted: Some(largest),
            };
        }
        BottomKUpdate::Unchanged
    }

    fn into_sorted(self) -> Vec<u64> {
        self.hashes.into_iter().collect()
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BottomKUpdate {
    Unchanged,
    Retained { evicted: Option<u64> },
}

#[derive(Debug, Error, Eq, PartialEq)]
pub enum SignatureSelectionError {
    #[error("invalid signature policy: {0}")]
    InvalidPolicy(String),
    #[error("signature budget must be greater than zero")]
    ZeroBudget,
    #[error("signature selection coordinate overflow")]
    Overflow,
}
