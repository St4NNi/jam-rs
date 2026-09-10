use crate::jidx_builder::{JidxBuildConfig, JidxBuildError, select_seeds};
use needletail::Sequence;

pub const CORE_BASES: u8 = 15;
pub const HAS_CONTEXT_21: u8 = 2;
pub const HAS_CONTEXT_31: u8 = 4;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct SharedKey {
    pub core: u32,
    pub context: u32,
    pub length: u8,
}

impl SharedKey {
    pub fn core(core: u32) -> Self {
        Self {
            core,
            context: 0,
            length: CORE_BASES,
        }
    }

    pub(crate) fn context_code(self) -> Option<u64> {
        if self.core >= 1 << 30 {
            return None;
        }
        match self.length {
            15 if self.context == 0 => Some(0),
            21 if self.context < 1 << 12 => Some((1 << 62) | u64::from(self.context)),
            31 => Some((2 << 62) | u64::from(self.context)),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SharedSeed {
    pub core: u32,
    pub context: u32,
    pub flags: u8,
    pub position: u64,
}

impl SharedSeed {
    pub fn canonical_orientation(self) -> bool {
        self.flags & 1 != 0
    }

    pub fn key(self, length: u8) -> Option<SharedKey> {
        let context = match length {
            15 => 0,
            21 if self.flags & HAS_CONTEXT_21 != 0 => self.context >> 20,
            31 if self.flags & HAS_CONTEXT_31 != 0 => self.context,
            _ => return None,
        };
        Some(SharedKey {
            core: self.core,
            context,
            length,
        })
    }
}

pub fn select_shared_seeds(
    sequence: &[u8],
    window: u16,
) -> Result<Vec<SharedSeed>, JidxBuildError> {
    if window == 0 {
        return Err(JidxBuildError::Invalid("shared minimizer window"));
    }
    let normalized = sequence.normalize(false);
    select_seeds(
        &normalized,
        JidxBuildConfig {
            k: CORE_BASES,
            minimizer_window: window,
            rescue_k15: false,
        },
    )?
    .into_iter()
    .map(|seed| {
        context_seed(
            &normalized,
            seed.position as usize,
            seed.packed_key as u32,
            seed.canonical_orientation,
            false,
        )
        .ok_or(JidxBuildError::Invalid("selected shared core"))
    })
    .collect()
}

pub fn context_seed(
    sequence: &[u8],
    position: usize,
    core: u32,
    reverse: bool,
    circular: bool,
) -> Option<SharedSeed> {
    if core >= 1 << 30
        || sequence.len() < 15
        || position >= sequence.len()
        || (!circular && position.checked_add(15)? > sequence.len())
    {
        return None;
    }
    let mut seed = SharedSeed {
        core,
        context: 0,
        flags: u8::from(reverse),
        position: position as u64,
    };
    if let Some(inner) = shell(sequence, position, reverse, circular, 3, 0) {
        seed.context = inner << 20;
        seed.flags |= HAS_CONTEXT_21;
        if let Some(outer) = shell(sequence, position, reverse, circular, 8, 3) {
            seed.context |= outer;
            seed.flags |= HAS_CONTEXT_31;
        }
    }
    Some(seed)
}

fn shell(
    sequence: &[u8],
    position: usize,
    reverse: bool,
    circular: bool,
    outer: i128,
    inner: i128,
) -> Option<u32> {
    if sequence.len() < (15 + 2 * outer) as usize {
        return None;
    }
    let mut packed = 0u32;
    for offset in (-outer..-inner).chain((15 + inner)..(15 + outer)) {
        let offset = if reverse { 14 - offset } else { offset };
        let absolute = position as i128 + offset;
        let at = if circular {
            absolute.rem_euclid(sequence.len() as i128) as usize
        } else {
            usize::try_from(absolute).ok()?
        };
        let mut base = match sequence.get(at)?.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        if reverse {
            base ^= 3;
        }
        packed = (packed << 2) | base;
    }
    Some(packed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn centered_contexts_follow_the_canonical_core_on_both_strands() {
        let sequence = b"AACGTTGCAACGATCGTAGGCTAACCGTAGCTACGATTCGACCGTAGCTAACGTC";
        let forward = select_shared_seeds(sequence, 4).unwrap();
        let reverse = select_shared_seeds(&sequence.reverse_complement(), 4).unwrap();
        let mirror = forward
            .iter()
            .map(|seed| {
                (
                    seed.core,
                    seed.context,
                    seed.flags ^ 1,
                    sequence.len() as u64 - 15 - seed.position,
                )
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(
            reverse
                .iter()
                .map(|seed| (seed.core, seed.context, seed.flags, seed.position))
                .collect::<BTreeSet<_>>(),
            mirror
        );
        for seed in forward
            .iter()
            .filter(|seed| seed.flags & HAS_CONTEXT_31 != 0)
        {
            let start = seed.position as usize;
            let context = &sequence[start - 8..start + 23];
            let context = if seed.canonical_orientation() {
                context.reverse_complement()
            } else {
                context.to_vec()
            };
            let bits = |sequence: &[u8]| {
                sequence
                    .bit_kmers(sequence.len() as u8, false)
                    .next()
                    .unwrap()
                    .1
                    .0 as u32
            };
            assert_eq!(bits(&context[8..23]), seed.core);
            assert_eq!(
                (bits(&context[5..8]) << 6) | bits(&context[23..26]),
                seed.context >> 20
            );
            assert_eq!(
                (bits(&context[..5]) << 10) | bits(&context[26..]),
                seed.context & ((1 << 20) - 1)
            );
        }
    }

    #[test]
    fn missing_context_never_discards_a_valid_core() {
        let sequence = b"ACGTTGCAACGATCG";
        let seeds = select_shared_seeds(sequence, 64).unwrap();
        assert_eq!(seeds.len(), 1);
        assert!(seeds[0].key(15).is_some());
        assert!(seeds[0].key(21).is_none());
        assert!(seeds[0].key(31).is_none());
        let mut extended = b"AAANNNACGTTGCAACGATCGTTTAAA".to_vec();
        let seed = context_seed(
            &extended,
            6,
            seeds[0].core,
            seeds[0].canonical_orientation(),
            false,
        )
        .unwrap();
        assert!(seed.key(21).is_none());
        extended[..6].fill(b'A');
        let seed = context_seed(
            &extended,
            6,
            seeds[0].core,
            seeds[0].canonical_orientation(),
            false,
        )
        .unwrap();
        assert!(seed.key(21).is_some());
        assert!(seed.key(31).is_none());
        assert!(select_shared_seeds(b"ACGT", 64).unwrap().is_empty());
        assert!(select_shared_seeds(sequence, 0).is_err());
    }

    #[test]
    fn nested_keys_and_circular_contexts_are_exact() {
        let sequence = b"ACGTTGCAACGATCGTAGGCTAACCGTAGCTACGATTCGA";
        let (position, key, reverse) = sequence.bit_kmers(15, true).next().unwrap();
        let seed = context_seed(sequence, position, key.0 as u32, reverse, true).unwrap();
        assert!(seed.key(31).is_some());
        let mut linear = sequence[sequence.len() - 8..].to_vec();
        linear.extend_from_slice(sequence);
        assert_eq!(
            context_seed(&linear, 8, key.0 as u32, reverse, false)
                .unwrap()
                .context,
            seed.context
        );
        assert!(
            SharedKey {
                core: 1 << 30,
                context: 0,
                length: 15
            }
            .context_code()
            .is_none()
        );
        assert!(
            SharedKey {
                core: 3,
                context: 1 << 12,
                length: 21
            }
            .context_code()
            .is_none()
        );
        assert!(
            SharedKey {
                core: 3,
                context: 0,
                length: 17
            }
            .context_code()
            .is_none()
        );
        assert_eq!(seed.key(15).unwrap().context_code(), Some(0));
        assert_eq!(
            seed.key(21).unwrap().context_code(),
            Some((1 << 62) | u64::from(seed.context >> 20))
        );
    }
}
