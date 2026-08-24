//! Gear fragments with bounded internal exact mini-sketches.

use super::fragments::{ExactFragment, GearConfig, GearError, fragment_bytes, fragment_sequence};
use super::tables::base_code;
use jamhash::jamhash_u64;
use std::collections::BTreeMap;

/// Parameters for an internal fragment signature.  No k smaller than 21 is
/// accepted because mini-sketch hits are later promoted to positional seeds.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MiniSketchConfig {
    pub k: u8,
    pub max_seeds: u16,
    pub scale: u64,
}

impl MiniSketchConfig {
    #[must_use]
    pub const fn k21(max_seeds: u16) -> Self {
        Self {
            k: 21,
            max_seeds,
            scale: 1,
        }
    }

    #[must_use]
    pub const fn k31(max_seeds: u16) -> Self {
        Self {
            k: 31,
            max_seeds,
            scale: 1,
        }
    }

    fn validate(self) -> Result<(), MiniSketchError> {
        if self.k < 21 || self.k > 31 {
            return Err(MiniSketchError::UnsupportedK(self.k));
        }
        if self.scale == 0 {
            return Err(MiniSketchError::ZeroScale);
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MiniSketchError {
    UnsupportedK(u8),
    ZeroScale,
    CoordinateOverflow,
    Fragment(GearError),
}

impl std::fmt::Display for MiniSketchError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::UnsupportedK(k) => {
                write!(formatter, "Gear mini-sketch k must be 21..31, got {k}")
            }
            Self::ZeroScale => formatter.write_str("Gear mini-sketch scale must be non-zero"),
            Self::CoordinateOverflow => formatter.write_str("Gear mini-sketch coordinate overflow"),
            Self::Fragment(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for MiniSketchError {}

impl From<GearError> for MiniSketchError {
    fn from(error: GearError) -> Self {
        Self::Fragment(error)
    }
}

/// One exact local mini-sketch seed.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct MiniSketchSeed {
    pub local_position: u32,
    pub k: u8,
    pub hash: u64,
    pub packed_kmer: u64,
    pub reverse: bool,
}

/// Fragment metadata plus its internal seed signature.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FragmentMiniSketch {
    pub fragment: ExactFragment,
    pub seeds: Vec<MiniSketchSeed>,
}

/// Exact local seed correspondence retained for caller-side anchor
/// verification.  The orientation flags are relative to the forward source
/// contig, not merely to each fragment's oriented byte view.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct MiniSketchLocalPair {
    pub query_local_position: u32,
    pub target_local_position: u32,
    pub query_reverse: bool,
    pub target_reverse: bool,
    pub relative_reverse: bool,
}

/// A candidate related-fragment match.  Internal exact keys have already been
/// compared, but callers must still verify the local sequence before creating
/// anchors or alignment evidence.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RelatedFragmentCandidate {
    pub fragment: ExactFragment,
    pub votes: usize,
    pub local_pairs: Vec<MiniSketchLocalPair>,
}

/// Digest and exact-key lookup index for one set of fragments.
#[derive(Clone, Debug, Default)]
pub struct FragmentMiniSketchIndex {
    fragments: Vec<FragmentMiniSketch>,
    by_digest: BTreeMap<u64, Vec<usize>>,
    by_seed_hash: BTreeMap<u64, Vec<(usize, MiniSketchSeed)>>,
}

impl FragmentMiniSketchIndex {
    pub fn insert(&mut self, sketch: FragmentMiniSketch) {
        let index = self.fragments.len();
        self.by_digest
            .entry(sketch.fragment.digest)
            .or_default()
            .push(index);
        for seed in sketch.seeds.iter().copied() {
            self.by_seed_hash
                .entry(seed.hash)
                .or_default()
                .push((index, seed));
        }
        self.fragments.push(sketch);
    }

    #[must_use]
    pub fn exact_digest_candidates(&self, digest: u64) -> Vec<&ExactFragment> {
        self.by_digest
            .get(&digest)
            .into_iter()
            .flatten()
            .filter_map(|index| self.fragments.get(*index).map(|sketch| &sketch.fragment))
            .collect()
    }

    #[must_use]
    pub fn seed_hash_candidate_count(&self, hash: u64) -> usize {
        self.by_seed_hash.get(&hash).map_or(0, Vec::len)
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.fragments.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.fragments.is_empty()
    }
}

fn canonical_kmer(sequence: &[u8]) -> Option<(u64, bool)> {
    if sequence.is_empty() || sequence.len() > 31 {
        return None;
    }
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for byte in sequence.iter().copied() {
        let code = base_code(byte)?;
        forward = (forward << 2) | u64::from(code);
        reverse = (reverse >> 2) | (u64::from(3 - code) << (2 * (sequence.len() - 1)));
    }
    if reverse < forward {
        Some((reverse, true))
    } else {
        Some((forward, false))
    }
}

/// Build a bounded exact mini-sketch for one fragment.
pub fn build_fragment_mini_sketch(
    sequence: &[u8],
    fragment: ExactFragment,
    config: MiniSketchConfig,
) -> Result<FragmentMiniSketch, MiniSketchError> {
    config.validate()?;
    let bytes = fragment_bytes(sequence, fragment)
        .ok_or(MiniSketchError::Fragment(GearError::CoordinateOverflow))?;
    let k = usize::from(config.k);
    let mut candidates = Vec::new();
    if bytes.len() >= k {
        for position in 0..=bytes.len() - k {
            let Some((packed_kmer, reverse)) = canonical_kmer(&bytes[position..position + k])
            else {
                continue;
            };
            let hash = jamhash_u64(packed_kmer);
            if hash == 0 || !hash.is_multiple_of(config.scale) {
                continue;
            }
            candidates.push(MiniSketchSeed {
                local_position: u32::try_from(position)
                    .map_err(|_| MiniSketchError::CoordinateOverflow)?,
                k: config.k,
                hash,
                packed_kmer,
                reverse,
            });
        }
    }
    candidates.sort_unstable_by_key(|seed| (seed.hash, seed.local_position, seed.packed_kmer));
    if config.max_seeds != 0 {
        candidates.truncate(usize::from(config.max_seeds));
    }
    candidates.sort_unstable_by_key(|seed| (seed.local_position, seed.hash, seed.packed_kmer));
    Ok(FragmentMiniSketch {
        fragment,
        seeds: candidates,
    })
}

/// Build sketches for all records produced by one Gear stream.
pub fn build_fragment_mini_sketches(
    sequence: &[u8],
    contig_id: u32,
    gear_config: GearConfig,
    mode: super::fragments::FragmentationMode,
    config: MiniSketchConfig,
) -> Result<Vec<FragmentMiniSketch>, MiniSketchError> {
    let fragments = fragment_sequence(sequence, contig_id, gear_config, mode)?;
    fragments
        .into_iter()
        .map(|fragment| build_fragment_mini_sketch(sequence, fragment, config))
        .collect()
}

/// Find related fragments through shared exact internal seed keys.  A digest
/// collision can only add work: the packed k-mer and scheme k are compared
/// before a vote is counted.
#[must_use]
pub fn lookup_related_fragments(
    index: &FragmentMiniSketchIndex,
    query: &FragmentMiniSketch,
    minimum_votes: u16,
) -> Vec<RelatedFragmentCandidate> {
    let mut votes: BTreeMap<usize, Vec<MiniSketchLocalPair>> = BTreeMap::new();
    for query_seed in query.seeds.iter().copied() {
        let Some(candidates) = index.by_seed_hash.get(&query_seed.hash) else {
            continue;
        };
        for (fragment_index, target_seed) in candidates {
            if target_seed.k != query_seed.k
                || target_seed.hash != query_seed.hash
                || target_seed.packed_kmer != query_seed.packed_kmer
            {
                continue;
            }
            let target_fragment = index.fragments[*fragment_index].fragment;
            let query_reverse = query_seed.reverse
                ^ matches!(
                    query.fragment.orientation,
                    super::fragments::FragmentOrientation::Reverse
                );
            let target_reverse = target_seed.reverse
                ^ matches!(
                    target_fragment.orientation,
                    super::fragments::FragmentOrientation::Reverse
                );
            votes
                .entry(*fragment_index)
                .or_default()
                .push(MiniSketchLocalPair {
                    query_local_position: query_seed.local_position,
                    target_local_position: target_seed.local_position,
                    query_reverse,
                    target_reverse,
                    relative_reverse: query_reverse ^ target_reverse,
                });
        }
    }
    votes
        .into_iter()
        .filter_map(|(fragment_index, mut local_pairs)| {
            local_pairs.sort_unstable();
            local_pairs.dedup();
            let votes = local_pairs.len();
            (votes >= usize::from(minimum_votes)).then(|| RelatedFragmentCandidate {
                fragment: index.fragments[fragment_index].fragment,
                votes,
                local_pairs,
            })
        })
        .collect()
}
