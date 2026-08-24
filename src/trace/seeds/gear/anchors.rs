//! Gear-selected exact canonical anchors.

use super::fragments::{
    FragmentOrientation, FragmentationMode, GearConfig, fragment_bytes, fragment_sequence,
};
use super::tables::base_code;
use jamhash::jamhash_u64;
use std::collections::BTreeMap;

/// Reverse-complement-symmetric 31-span/21-informative-base pattern.
///
/// The set contains the first ten bases, the centre base, and the last ten
/// bases.  Its 31-bit mask is equal to its 31-bit reversal.
pub const SPACED31_WEIGHT21_MASK: u32 = 0x7fe0_83ff;

/// Anchor scheme selected at Gear boundaries.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum GearAnchorScheme {
    CanonicalK31,
    CanonicalK21,
    Spaced31Weight21,
}

impl GearAnchorScheme {
    #[must_use]
    pub const fn scheme_id(self) -> u32 {
        match self {
            Self::CanonicalK31 => 1,
            Self::CanonicalK21 => 2,
            Self::Spaced31Weight21 => 3,
        }
    }

    #[must_use]
    pub const fn span(self) -> u8 {
        match self {
            Self::CanonicalK31 => 31,
            Self::CanonicalK21 => 21,
            Self::Spaced31Weight21 => 31,
        }
    }

    #[must_use]
    pub const fn informative_bases(self) -> u8 {
        match self {
            Self::CanonicalK31 => 31,
            Self::CanonicalK21 => 21,
            Self::Spaced31Weight21 => 21,
        }
    }

    #[must_use]
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::CanonicalK31 => "gear-anchor-k31",
            Self::CanonicalK21 => "gear-anchor-k21-rescue",
            Self::Spaced31Weight21 => "gear-anchor-spaced31-weight21",
        }
    }
}

/// Limits for one Gear-selected anchor stream.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GearAnchorConfig {
    pub scheme: GearAnchorScheme,
    /// Zero means unbounded.  A bound is applied after deterministic position
    /// ordering, so it never changes a prefix under hash-map iteration.
    pub max_anchors: u32,
}

impl GearAnchorConfig {
    #[must_use]
    pub const fn k31() -> Self {
        Self {
            scheme: GearAnchorScheme::CanonicalK31,
            max_anchors: 0,
        }
    }

    #[must_use]
    pub const fn k21_rescue() -> Self {
        Self {
            scheme: GearAnchorScheme::CanonicalK21,
            max_anchors: 0,
        }
    }

    #[must_use]
    pub const fn spaced31_weight21() -> Self {
        Self {
            scheme: GearAnchorScheme::Spaced31Weight21,
            max_anchors: 0,
        }
    }
}

/// One exact anchor selected at a Gear fragment boundary.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct GearAnchor {
    pub position: u64,
    pub span: u8,
    pub informative_bases: u8,
    pub scheme: GearAnchorScheme,
    pub hash: u64,
    /// Exact canonical packed selected bases.  Hash equality alone is never
    /// sufficient to accept this record.
    pub packed_selected: u64,
    pub reverse: bool,
    pub fragment_index: u32,
}

impl GearAnchor {
    #[must_use]
    pub fn exact_key(self) -> (u32, u8, u8, u64) {
        (
            self.scheme.scheme_id(),
            self.span,
            self.informative_bases,
            self.packed_selected,
        )
    }
}

/// Hash-only lookup index.  The returned candidates still require
/// `verify_anchor_candidate` or `verify_anchor_sequences`.
#[derive(Clone, Debug, Default)]
pub struct GearAnchorIndex {
    by_hash: BTreeMap<u64, Vec<GearAnchor>>,
}

impl GearAnchorIndex {
    pub fn insert(&mut self, anchor: GearAnchor) {
        self.by_hash.entry(anchor.hash).or_default().push(anchor);
    }

    #[must_use]
    pub fn candidates(&self, hash: u64) -> &[GearAnchor] {
        self.by_hash.get(&hash).map_or(&[], Vec::as_slice)
    }

    #[must_use]
    pub fn verified_candidates(&self, query: GearAnchor) -> Vec<GearAnchorVerification> {
        self.candidates(query.hash)
            .iter()
            .copied()
            .filter_map(|candidate| {
                verify_anchor_candidate(query, candidate)
                    .then_some(GearAnchorVerification { query, candidate })
            })
            .collect()
    }
}

/// Exact key result after a candidate hash has been checked.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct GearAnchorVerification {
    pub query: GearAnchor,
    pub candidate: GearAnchor,
}

/// Error from invalid Gear anchor configuration or sequence coordinates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GearSelectedAnchorError {
    InvalidScheme,
    CoordinateOverflow,
}

impl std::fmt::Display for GearSelectedAnchorError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidScheme => formatter.write_str("invalid Gear anchor scheme"),
            Self::CoordinateOverflow => formatter.write_str("Gear anchor coordinate overflow"),
        }
    }
}

impl std::error::Error for GearSelectedAnchorError {}

fn canonical_packed(sequence: &[u8]) -> Option<(u64, bool)> {
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

fn spaced_packed(sequence: &[u8]) -> Option<(u64, bool)> {
    if sequence.len() != 31 {
        return None;
    }
    let positions = (0usize..31).filter(|position| (SPACED31_WEIGHT21_MASK >> position) & 1 == 1);
    let mut forward = 0u64;
    for position in positions.clone() {
        forward = (forward << 2) | u64::from(base_code(sequence[position])?);
    }
    let mut reverse = 0u64;
    for position in positions {
        let source = sequence[30 - position];
        reverse = (reverse << 2) | u64::from(3 - base_code(source)?);
    }
    if reverse < forward {
        Some((reverse, true))
    } else {
        Some((forward, false))
    }
}

fn key_at(sequence: &[u8], position: usize, scheme: GearAnchorScheme) -> Option<(u64, bool)> {
    let span = usize::from(scheme.span());
    let window = sequence.get(position..position.checked_add(span)?)?;
    match scheme {
        GearAnchorScheme::CanonicalK31 | GearAnchorScheme::CanonicalK21 => canonical_packed(window),
        GearAnchorScheme::Spaced31Weight21 => spaced_packed(window),
    }
}

fn make_anchor(
    sequence: &[u8],
    position: usize,
    config: GearAnchorConfig,
    fragment_index: usize,
) -> Result<Option<GearAnchor>, GearSelectedAnchorError> {
    let Some((packed_selected, reverse)) = key_at(sequence, position, config.scheme) else {
        return Ok(None);
    };
    let hash = jamhash_u64(packed_selected);
    if hash == 0 {
        return Ok(None);
    }
    let position =
        u64::try_from(position).map_err(|_| GearSelectedAnchorError::CoordinateOverflow)?;
    let fragment_index =
        u32::try_from(fragment_index).map_err(|_| GearSelectedAnchorError::CoordinateOverflow)?;
    Ok(Some(GearAnchor {
        position,
        span: config.scheme.span(),
        informative_bases: config.scheme.informative_bases(),
        scheme: config.scheme,
        hash,
        packed_selected,
        reverse,
        fragment_index,
    }))
}

/// Build one canonical anchor at an explicit forward-sequence coordinate.
/// This helper is also used to verify that anchors selected from a reverse
/// fragment report orientation relative to the forward source.
pub fn canonical_anchor_at(
    sequence: &[u8],
    position: u64,
    scheme: GearAnchorScheme,
) -> Option<GearAnchor> {
    let position = usize::try_from(position).ok()?;
    let (packed_selected, reverse) = key_at(sequence, position, scheme)?;
    let hash = jamhash_u64(packed_selected);
    (hash != 0).then_some(GearAnchor {
        position: u64::try_from(position).ok()?,
        span: scheme.span(),
        informative_bases: scheme.informative_bases(),
        scheme,
        hash,
        packed_selected,
        reverse,
        fragment_index: 0,
    })
}

/// Select exact canonical anchors at the starts of Gear fragments.  The
/// caller can choose `BothStrands` to inspect both fragmentations; every
/// returned anchor still needs target sequence verification after lookup.
pub fn gear_selected_anchors(
    sequence: &[u8],
    config: GearConfig,
    mode: FragmentationMode,
    anchor_config: GearAnchorConfig,
) -> Result<Vec<GearAnchor>, GearSelectedAnchorError> {
    let fragments = fragment_sequence(sequence, 0, config, mode)
        .map_err(|_| GearSelectedAnchorError::CoordinateOverflow)?;
    let mut anchors = Vec::new();
    for (fragment_index, fragment) in fragments.iter().copied().enumerate() {
        let Some(bytes) = fragment_bytes(sequence, fragment) else {
            continue;
        };
        let Some(anchor) = make_anchor(&bytes, 0, anchor_config, fragment_index)? else {
            continue;
        };
        let position = match fragment.orientation {
            FragmentOrientation::Forward => fragment.start,
            FragmentOrientation::Reverse => fragment
                .start
                .checked_add(u64::from(fragment.length))
                .and_then(|end| end.checked_sub(u64::from(anchor.span)))
                .ok_or(GearSelectedAnchorError::CoordinateOverflow)?,
        };
        anchors.push(GearAnchor {
            position,
            // `anchor.reverse` is relative to the oriented fragment bytes.
            // Convert it back to the forward source coordinate orientation.
            reverse: anchor.reverse ^ (fragment.orientation == FragmentOrientation::Reverse),
            ..anchor
        });
    }
    anchors.sort_unstable_by_key(|anchor| {
        (
            anchor.position,
            anchor.scheme,
            anchor.fragment_index,
            anchor.packed_selected,
        )
    });
    if anchor_config.max_anchors != 0 {
        anchors.truncate(
            usize::try_from(anchor_config.max_anchors)
                .map_err(|_| GearSelectedAnchorError::CoordinateOverflow)?,
        );
    }
    Ok(anchors)
}

/// Verify a candidate using scheme metadata and exact packed selected bases.
/// This is the mandatory post-hash path for an indexed anchor.
#[must_use]
pub fn verify_anchor_candidate(query: GearAnchor, candidate: GearAnchor) -> bool {
    query.hash == candidate.hash && query.exact_key() == candidate.exact_key()
}

/// Verify a query anchor against target sequence bytes.  `target_position` is
/// a forward-coordinate start and `target_reverse` requests reverse-complement
/// verification of that target window.
#[must_use]
pub fn verify_anchor_sequences(
    query_sequence: &[u8],
    query: GearAnchor,
    target_sequence: &[u8],
    target_position: u64,
    target_reverse: bool,
) -> bool {
    let Ok(query_position) = usize::try_from(query.position) else {
        return false;
    };
    let Ok(target_position) = usize::try_from(target_position) else {
        return false;
    };
    let Some((query_packed, _)) = key_at(query_sequence, query_position, query.scheme) else {
        return false;
    };
    let span = usize::from(query.span);
    let Some(target_end) = target_position.checked_add(span) else {
        return false;
    };
    let Some(window) = target_sequence.get(target_position..target_end) else {
        return false;
    };
    let target_window = if target_reverse {
        super::fragments::reverse_complement(window)
    } else {
        window.to_vec()
    };
    let Some((target_packed, _)) = key_at(&target_window, 0, query.scheme) else {
        return false;
    };
    let target_hash = jamhash_u64(target_packed);
    target_hash != 0 && target_hash == query.hash && target_packed == query_packed
}
