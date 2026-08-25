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
        let requested_budget = self.policy.contig_signature_budget(length);
        let contig_id =
            u32::try_from(self.contig_count).map_err(|_| SignatureSelectionError::Overflow)?;
        let mut contig = self
            .policy
            .spatial_segment_bases()
            .is_none()
            .then(|| BottomK::new(requested_budget as usize))
            .transpose()?;
        let mut spatial_minima = Vec::<(u64, u64)>::new();
        let mut active_segment = None;
        let mut active_candidates = Vec::<(u64, u64)>::new();
        let mut eligible_kmers = 0u64;
        let normalized = sequence.normalize(false);
        for (position, kmer, _) in normalized.bit_kmers(self.policy.k, true) {
            let hash = jamhash_u64_v1(kmer.0);
            if hash == 0 {
                continue;
            }
            eligible_kmers = eligible_kmers.saturating_add(1);
            if let Some(segment_bases) = self.policy.spatial_segment_bases() {
                let position =
                    u64::try_from(position).map_err(|_| SignatureSelectionError::Overflow)?;
                let segment = position / u64::from(segment_bases);
                if active_segment != Some(segment) {
                    if let Some(previous) = active_segment {
                        append_segment_signatures(
                            &mut spatial_minima,
                            previous,
                            &active_candidates,
                            self.policy.spatial_signatures_per_segment().unwrap_or(1),
                        );
                    }
                    active_segment = Some(segment);
                    active_candidates.clear();
                }
                active_candidates.push((position, hash));
            } else if let Some(contig) = contig.as_mut() {
                contig.insert(hash);
            }
            if let BottomKUpdate::Retained { evicted } = self.whole.insert(hash) {
                if let Some(evicted) = evicted {
                    self.whole_sources.remove(&evicted);
                }
                self.whole_sources.entry(hash).or_insert(contig_id);
            }
        }
        if let Some(segment) = active_segment {
            append_segment_signatures(
                &mut spatial_minima,
                segment,
                &active_candidates,
                self.policy.spatial_signatures_per_segment().unwrap_or(1),
            );
        }
        let hashes = match contig {
            Some(contig) => contig.into_sorted(),
            None => spatial_hashes(spatial_minima, self.policy.contig_budget.maximum),
        };
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

fn append_segment_signatures(
    output: &mut Vec<(u64, u64)>,
    segment: u64,
    candidates: &[(u64, u64)],
    requested: u32,
) {
    let mut by_hash = BTreeMap::<u64, Vec<u64>>::new();
    for (position, hash) in candidates {
        by_hash.entry(*hash).or_default().push(*position);
    }
    let Some((&first_hash, first_positions)) = by_hash.first_key_value() else {
        return;
    };
    output.push((segment, first_hash));
    if requested < 2 {
        return;
    }
    let separated = by_hash.iter().skip(1).find_map(|(hash, positions)| {
        positions
            .iter()
            .any(|position| {
                first_positions
                    .iter()
                    .any(|first| position.abs_diff(*first) >= 32)
            })
            .then_some(*hash)
    });
    if let Some(second_hash) = separated.or_else(|| by_hash.keys().nth(1).copied()) {
        output.push((segment, second_hash));
    }
}

fn spatial_hashes(minima: Vec<(u64, u64)>, maximum: u32) -> Vec<u64> {
    let mut seen = BTreeSet::new();
    let mut unique = minima
        .into_iter()
        .filter_map(|(segment, hash)| seen.insert(hash).then_some((segment, hash)))
        .collect::<Vec<_>>();
    let maximum = usize::try_from(maximum).unwrap_or(usize::MAX);
    if unique.len() > maximum {
        let count = unique.len();
        unique = (0..maximum)
            .map(|index| unique[index.saturating_mul(count) / maximum])
            .collect();
    }
    let mut hashes = unique.into_iter().map(|(_, hash)| hash).collect::<Vec<_>>();
    hashes.sort_unstable();
    hashes
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_signature_segment_prefers_a_separated_second_position() {
        let mut selected = Vec::new();
        append_segment_signatures(&mut selected, 0, &[(10, 1), (20, 2), (100, 3)], 2);
        assert_eq!(selected, vec![(0, 1), (0, 3)]);

        selected.clear();
        append_segment_signatures(&mut selected, 0, &[(10, 1), (20, 2)], 2);
        assert_eq!(selected, vec![(0, 1), (0, 2)]);
    }
}
