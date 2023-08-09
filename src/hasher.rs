use std::hash::Hasher;

/// Taken from finch-rs: https://github.com/onecodex/finch-rs/blob/master/lib/src/sketch_schemes/hashing.rs
/// If we're using a `HashMap` where the keys themselves are hashes, it's
/// a little silly to re-hash them. That's where the `NoHashHasher` comes in.
#[derive(Default)]
pub struct NoHashHasher(u64);

impl Hasher for NoHashHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        *self = NoHashHasher(
            (u64::from(bytes[0]) << 24)
                + (u64::from(bytes[1]) << 16)
                + (u64::from(bytes[2]) << 8)
                + u64::from(bytes[3]),
        );
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.0 = i as u64;
    }

    fn finish(&self) -> u64 {
        self.0
    }
}
