//! Deterministic Gear tables and rolling state.

/// Independent table seeds used by the boundary-invariance experiments.  The
/// values are fixed so benchmark and archive builds remain reproducible while
/// still exercising uncorrelated streams.
pub const INDEPENDENT_TABLE_SEEDS: [u64; 8] = [
    0x243f_6a88_85a3_08d3,
    0x1319_8a2e_0370_7344,
    0xa409_3822_299f_31d0,
    0x082e_fa98_ec4e_6c89,
    0x4528_21e6_38d0_1377,
    0xbe54_66cf_34e9_0c6c,
    0xc0ac_29b7_c97c_50dd,
    0x3f84_d5b5_b547_0917,
];

/// The three bounded table families used by the Gear experiments.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum GearTableKind {
    /// One independently generated value for every byte value.
    SingleBase,
    /// A value derived from the previous/current byte pair.
    Dinucleotide,
    /// A value derived from the most recent four unambiguous bases packed as
    /// one byte.  An ambiguous run uses the single-byte fallback value.
    PackedFourBase,
}

impl GearTableKind {
    #[must_use]
    pub const fn identifier(self) -> &'static str {
        match self {
            Self::SingleBase => "gear-single-base-v1",
            Self::Dinucleotide => "gear-dinucleotide-v1",
            Self::PackedFourBase => "gear-packed-four-base-v1",
        }
    }
}

/// A deterministic Gear table.  The table is generated from a stored seed,
/// so archive builders and experiments do not depend on random state.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GearTable {
    kind: GearTableKind,
    seed: u64,
    values: [u64; 256],
}

impl GearTable {
    #[must_use]
    pub fn new(kind: GearTableKind, seed: u64) -> Self {
        let mut values = [0u64; 256];
        let mut index = 0usize;
        while index < values.len() {
            values[index] = splitmix64(seed.wrapping_add(index as u64));
            index += 1;
        }
        Self { kind, seed, values }
    }

    #[must_use]
    pub const fn kind(&self) -> GearTableKind {
        self.kind
    }

    #[must_use]
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    #[must_use]
    pub const fn identifier(&self) -> &'static str {
        self.kind.identifier()
    }

    #[must_use]
    pub(crate) fn value(&self, previous: Option<u8>, current: u8, packed_four: Option<u8>) -> u64 {
        match self.kind {
            GearTableKind::SingleBase => self.values[usize::from(current)],
            GearTableKind::Dinucleotide => {
                let previous = previous.unwrap_or(0);
                let pair = (u16::from(previous) << 8) | u16::from(current);
                splitmix64(self.seed ^ u64::from(pair).wrapping_mul(0x9e37_79b9_7f4a_7c15))
            }
            GearTableKind::PackedFourBase => packed_four
                .map(|packed| self.values[usize::from(packed)])
                .unwrap_or_else(|| self.values[usize::from(current)]),
        }
    }
}

/// Rolling state for one content-defined chunk.
#[derive(Clone, Debug, Default)]
pub(crate) struct GearHasher {
    hash: u64,
    previous: Option<u8>,
    packed_four: u8,
    packed_len: u8,
}

impl GearHasher {
    pub(crate) fn update(&mut self, table: &GearTable, byte: u8) -> u64 {
        let base = base_code(byte);
        let packed = match base {
            Some(code) => {
                self.packed_four = (self.packed_four << 2) | code;
                self.packed_len = self.packed_len.saturating_add(1).min(4);
                (self.packed_len == 4).then_some(self.packed_four)
            }
            None => {
                self.packed_four = 0;
                self.packed_len = 0;
                None
            }
        };
        let value = table.value(self.previous, byte, packed);
        self.hash = self.hash.rotate_left(1).wrapping_add(value);
        self.previous = Some(byte);
        self.hash
    }

    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

#[must_use]
pub(crate) fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[must_use]
pub(crate) fn base_code(byte: u8) -> Option<u8> {
    match byte.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

#[must_use]
pub(crate) fn normalized_byte(byte: u8) -> u8 {
    match byte.to_ascii_uppercase() {
        b'U' => b'T',
        value => value,
    }
}
