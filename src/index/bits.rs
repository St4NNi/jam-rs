//! Little-endian fields, bit-packed arrays and varints shared by the index files.

/// Zero bytes after every bit-packed array, so a value is always read with one 8-byte load.
pub(crate) const BIT_PADDING: usize = 8;

/// The u64 at `offset`, or `None` past the end.
pub(crate) fn u64_at(bytes: &[u8], offset: usize) -> Option<u64> {
    let field = bytes.get(offset..offset.checked_add(8)?)?;
    Some(u64::from_le_bytes(field.try_into().ok()?))
}

/// The u32 at `offset`, or `None` past the end.
pub(crate) fn u32_at(bytes: &[u8], offset: usize) -> Option<u32> {
    let field = bytes.get(offset..offset.checked_add(4)?)?;
    Some(u32::from_le_bytes(field.try_into().ok()?))
}

/// The u16 at `offset`, or `None` past the end.
pub(crate) fn u16_at(bytes: &[u8], offset: usize) -> Option<u16> {
    let field = bytes.get(offset..offset.checked_add(2)?)?;
    Some(u16::from_le_bytes(field.try_into().ok()?))
}

/// Bits needed to write `value` (0 for 0).
pub(crate) fn bits_needed(value: u64) -> u32 {
    64 - value.leading_zeros()
}

/// The value of `width` bits (at most 57) that starts at bit `bit` of a bit-packed array, the
/// first value in the low bits. Bits past the end read as zero.
pub(crate) fn read_bits(bytes: &[u8], bit: u64, width: u32) -> u64 {
    if width == 0 {
        return 0;
    }
    let start = (bit / 8) as usize;
    let word = match bytes.get(start..start.saturating_add(8)) {
        Some(word) => u64::from_le_bytes(word.try_into().unwrap_or_default()),
        None => {
            let mut word = [0u8; 8];
            let tail = bytes.get(start..).unwrap_or_default();
            word[..tail.len()].copy_from_slice(tail);
            u64::from_le_bytes(word)
        }
    };
    (word >> (bit % 8)) & ((1u64 << width) - 1)
}

/// Appends values of a fixed width to a bit-packed array.
pub(crate) struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    start: usize,
    bits: u64,
}

impl<'a> BitWriter<'a> {
    /// Starts a bit-packed array at the end of `out`.
    pub(crate) fn new(out: &'a mut Vec<u8>) -> Self {
        let start = out.len();
        Self {
            out,
            start,
            bits: 0,
        }
    }

    /// Appends the low `width` bits (at most 57) of `value`.
    pub(crate) fn push(&mut self, value: u64, width: u32) {
        if width == 0 {
            return;
        }
        let end = self.start + (self.bits + u64::from(width)).div_ceil(8) as usize;
        if self.out.len() < end {
            self.out.resize(end, 0);
        }
        let mut value = value & ((1u64 << width) - 1);
        let mut bit = self.bits;
        let mut left = width;
        while left > 0 {
            let byte = self.start + (bit / 8) as usize;
            let shift = (bit % 8) as u32;
            let take = left.min(8 - shift);
            self.out[byte] |= ((value & ((1 << take) - 1)) as u8) << shift;
            value >>= take;
            bit += u64::from(take);
            left -= take;
        }
        self.bits = bit;
    }

    /// Ends the array at its last whole byte. A file adds `BIT_PADDING` after its last array.
    pub(crate) fn finish(self) {
        let end = self.start + self.bits.div_ceil(8) as usize;
        self.out.resize(end, 0);
    }
}

/// Appends `value` as a varint: seven bits per byte, the high bit set on all but the last.
pub(crate) fn put_varint(out: &mut Vec<u8>, mut value: u64) {
    while value >= 0x80 {
        out.push(value as u8 | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

/// Reads the varint at `*at` and moves past it, or `None` when it is cut off or too long.
pub(crate) fn read_varint(bytes: &[u8], at: &mut usize) -> Option<u64> {
    let mut value = 0u64;
    for shift in (0..64).step_by(7) {
        let byte = *bytes.get(*at)?;
        *at += 1;
        value |= u64::from(byte & 0x7f) << shift;
        if byte < 0x80 {
            return Some(value);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::rng;
    use rand::RngExt;

    #[test]
    fn bit_arrays_round_trip() {
        let mut draws = rng(3);
        for width in [0u32, 1, 5, 13, 32, 41, 57] {
            let values: Vec<u64> = (0..200)
                .map(|_| draws.random::<u64>() & ((1u64 << width) - 1))
                .collect();
            let mut out = vec![0xaa];
            let mut writer = BitWriter::new(&mut out);
            for value in &values {
                writer.push(*value, width);
            }
            writer.finish();
            assert_eq!(out.len(), 1 + (200 * width as usize).div_ceil(8));
            out.resize(out.len() + BIT_PADDING, 0);
            for (index, value) in values.iter().enumerate() {
                assert_eq!(
                    read_bits(&out[1..], index as u64 * u64::from(width), width),
                    *value
                );
            }
        }
        assert_eq!(read_bits(&[0xff], 4, 8), 0x0f);
    }

    #[test]
    fn varints_round_trip() {
        let values = [0u64, 1, 127, 128, 300, u32::MAX as u64, u64::MAX];
        let mut out = Vec::new();
        for value in values {
            put_varint(&mut out, value);
        }
        let mut at = 0;
        for value in values {
            assert_eq!(read_varint(&out, &mut at), Some(value));
        }
        assert_eq!(read_varint(&out, &mut at), None);
        assert_eq!(read_varint(&[0x80; 11], &mut 0), None);
        assert_eq!(u64_at(&[0; 7], 0), None);
        assert_eq!(u32_at(&[1, 0, 0, 0], 0), Some(1));
        assert_eq!(u16_at(&[1, 2], 1), None);
        assert_eq!(bits_needed(0), 0);
        assert_eq!(bits_needed(255), 8);
    }
}
