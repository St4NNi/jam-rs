#[cfg(not(target_endian = "little"))]
compile_error!("JAM format requires a little-endian platform");

use bytemuck::{Pod, Zeroable};

pub const MAGIC: [u8; 4] = *b"JAM\0";
pub const VERSION: u32 = 3;

pub const PAGE_SIZE: usize = 4096;

#[inline]
pub const fn align_to_page(offset: usize) -> usize {
    (offset + PAGE_SIZE - 1) & !(PAGE_SIZE - 1)
}
pub const BUCKET_COUNT: usize = 256;
pub const BUCKET_BITS: u8 = 8;
pub const ENTRY_SIZE: usize = 12;
pub const HEADER_SIZE: usize = 160;
pub const BUCKET_META_SIZE: usize = 32;
pub const BUCKET_TABLE_SIZE: usize = BUCKET_COUNT * BUCKET_META_SIZE;
pub const DATA_START: usize = HEADER_SIZE + BUCKET_TABLE_SIZE;

#[inline(always)]
pub fn bucket_id(hash: u64) -> usize {
    (hash & 0xFF) as usize
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct Header {
    pub magic: [u8; 4],
    pub version: u32,
    pub flags: u64,

    pub entry_count: u64,
    pub unique_hash_count: u64,
    pub sample_count: u32,
    pub bucket_count: u16,
    pub bucket_bits: u8,
    pub entry_size: u8,

    pub hash_threshold: u64,
    pub kmer_size: u8,
    pub _param_reserved: [u8; 3],
    pub min_entropy: f32,

    pub bucket_table_offset: u64,
    pub entries_offset: u64,
    pub filters_offset: u64,
    pub bias_table_offset: u64,

    pub entries_size: u64,
    pub filters_size: u64,
    pub bias_table_size: u64,

    pub sample_names_offset: u64,
    pub sample_names_size: u64,
    pub sample_sizes_offset: u64,
    pub sample_sizes_size: u64,

    pub _padding: [u8; 16],
}

pub const FLAG_HAS_BIAS_TABLE: u64 = 1 << 0;

const _: () = assert!(std::mem::size_of::<Header>() == 160);

impl Header {
    pub fn validate(&self) -> Result<(), FormatError> {
        if self.magic != MAGIC {
            return Err(FormatError::InvalidMagic(self.magic));
        }
        if self.version != VERSION {
            return Err(FormatError::UnsupportedVersion(self.version));
        }
        if self.bucket_count != BUCKET_COUNT as u16 {
            return Err(FormatError::InvalidBucketCount(self.bucket_count));
        }
        if self.entry_size != ENTRY_SIZE as u8 {
            return Err(FormatError::InvalidEntrySize(self.entry_size));
        }
        if self.hash_threshold == 0 {
            return Err(FormatError::InvalidHashThreshold);
        }
        Ok(())
    }
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable, Default)]
pub struct BucketMeta {
    pub entry_offset: u64,
    pub entry_count: u64,
    pub filter_offset: u64,
    pub filter_size: u64,
}

const _: () = assert!(std::mem::size_of::<BucketMeta>() == 32);

#[repr(C, packed)]
#[derive(Debug, Clone, Copy, Pod, Zeroable, PartialEq, Eq, PartialOrd, Ord)]
pub struct Entry {
    pub hash: u64,
    pub sample_id: u32,
}

const _: () = assert!(std::mem::size_of::<Entry>() == 12);

impl Entry {
    #[inline]
    pub fn new(hash: u64, sample_id: u32) -> Self {
        Self { hash, sample_id }
    }

    #[inline]
    pub fn bucket_id(&self) -> usize {
        bucket_id(self.hash)
    }
}

#[derive(Debug, thiserror::Error)]
pub enum FormatError {
    #[error("Invalid magic bytes: {0:?}")]
    InvalidMagic([u8; 4]),

    #[error("Unsupported version: {0}")]
    UnsupportedVersion(u32),

    #[error("Invalid bucket count: {0}")]
    InvalidBucketCount(u16),

    #[error("Invalid entry size: {0}")]
    InvalidEntrySize(u8),

    #[error("Invalid hash threshold: must be > 0")]
    InvalidHashThreshold,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_struct_sizes() {
        assert_eq!(std::mem::size_of::<Header>(), 160);
        assert_eq!(std::mem::size_of::<BucketMeta>(), 32);
        assert_eq!(std::mem::size_of::<Entry>(), 12);
    }

    #[test]
    fn test_bucket_id() {
        assert_eq!(bucket_id(0x0000_0000_0000_0000), 0);
        assert_eq!(bucket_id(0x0000_0000_0000_00FF), 255);
        assert_eq!(bucket_id(0xFFFF_FFFF_FFFF_FF00), 0);
        assert_eq!(bucket_id(0xABCD_EF12_3456_7842), 0x42);
    }

    #[test]
    fn test_entry_ordering() {
        let e1 = Entry::new(100, 1);
        let e2 = Entry::new(100, 2);
        let e3 = Entry::new(200, 1);

        assert!(e1 < e2);
        assert!(e2 < e3);
        assert!(e1 < e3);
    }

    #[test]
    fn test_bucket_id_distribution() {
        let threshold: u64 = (u64::MAX as f64 * 0.001) as u64;
        let mut bucket_counts = [0usize; 256];

        for i in 0..100_000u64 {
            let hash = i.wrapping_mul(0x517cc1b727220a95) % threshold;
            bucket_counts[bucket_id(hash)] += 1;
        }

        let avg = 100_000 / 256;
        for (i, &count) in bucket_counts.iter().enumerate() {
            let deviation = (count as f64 - avg as f64).abs() / avg as f64;
            assert!(deviation < 0.3, "Bucket {} has skewed count: {}", i, count);
        }
    }

    #[test]
    fn test_header_validate_valid() {
        let mut header = Header::zeroed();
        header.magic = MAGIC;
        header.version = VERSION;
        header.bucket_count = BUCKET_COUNT as u16;
        header.entry_size = ENTRY_SIZE as u8;
        header.hash_threshold = u64::MAX;
        assert!(header.validate().is_ok());
    }

    #[test]
    fn test_header_validate_zero_threshold() {
        let mut header = Header::zeroed();
        header.magic = MAGIC;
        header.version = VERSION;
        header.bucket_count = BUCKET_COUNT as u16;
        header.entry_size = ENTRY_SIZE as u8;
        header.hash_threshold = 0;
        assert!(matches!(
            header.validate(),
            Err(FormatError::InvalidHashThreshold)
        ));
    }

    #[test]
    fn test_header_validate_bad_magic() {
        let mut header = Header::zeroed();
        header.magic = *b"BAD\0";
        header.version = VERSION;
        header.bucket_count = BUCKET_COUNT as u16;
        header.entry_size = ENTRY_SIZE as u8;
        header.hash_threshold = u64::MAX;
        assert!(matches!(
            header.validate(),
            Err(FormatError::InvalidMagic(_))
        ));
    }

    #[test]
    fn test_header_validate_bad_version() {
        let mut header = Header::zeroed();
        header.magic = MAGIC;
        header.version = 99;
        header.bucket_count = BUCKET_COUNT as u16;
        header.entry_size = ENTRY_SIZE as u8;
        header.hash_threshold = u64::MAX;
        assert!(matches!(
            header.validate(),
            Err(FormatError::UnsupportedVersion(99))
        ));
    }

    #[test]
    fn test_entry_bucket_id() {
        let entry = Entry::new(0xABCD_EF12_3456_7842, 5);
        assert_eq!(entry.bucket_id(), 0x42);
    }
}
