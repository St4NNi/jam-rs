//! Compact positional seed sections for JMA format 1.
//!
//! This module is deliberately independent of the legacy fixed-record seed
//! pages.  It keeps the membership filter, exact keys, posting headers, and
//! position payloads in separate ranges so a caller can stop after a header
//! lookup for repetitive seeds.

mod builder;
mod collection;
mod format;
mod reader;

pub use builder::{SeedBuildConfig, encode_seed_section};
pub use collection::{
    COLLECTION_ENTRY_SIZE, COLLECTION_HEADER_SIZE, COLLECTION_MAGIC, COLLECTION_VERSION,
    SeedCollection, SeedCollectionEntry, SeedCollectionHeader, decode_collection_directory,
    decode_collection_header, encode_seed_collection, validate_collection_entries,
};
pub use format::{
    COMPACT_HEADER_SIZE, COMPACT_SEED_MAGIC, COMPACT_SEED_VERSION, KEY_PAGE_RECORD_SIZE,
    KeyEncoding, OccurrenceEncoding, POSITION_PAGE_RECORD_SIZE, POSTING_FLAG_POSITION_BEARING,
    POSTING_HEADER_RECORD_SIZE, PREFIX_RECORD_SIZE, PostingClass, SeedBinding, SeedHeader,
    SeedKeyPage, SeedPositionPage, SeedPostingHeader, SeedPrefix,
};
pub use reader::{
    LookupOptions, SeedDirectoryPrefix, SeedIndex, SeedKeyLookupPrefix, SeedLookupMetrics,
    SeedLookupResult, SeedMatch, decode_directory_prefix, decode_key_lookup_prefix,
    decode_key_page_records, decode_position_page_record, decode_position_posting,
    decode_posting_header_record, decode_seed_header,
};

use crate::jma::{JmaError, JmaResult, SeedOccurrence, SeedQuery};

/// The immutable hash identity used by JMA seed sections.
pub const JAMHASH_ALGORITHM_ID: u32 = 1;

/// The compact section accepts the released k-mer representation but never a
/// seed shorter than k=21.
pub fn validate_query(query: SeedQuery, expected_k: u8) -> JmaResult<()> {
    if query.k != expected_k {
        return Err(JmaError::CorruptSection(format!(
            "seed query k={} does not match section k={expected_k}",
            query.k
        )));
    }
    if !(21..=32).contains(&query.k)
        || !crate::jma::format::valid_canonical_kmer(query.k, query.canonical_kmer)
    {
        return Err(JmaError::CorruptSection(
            "seed query has an invalid packed canonical k-mer".to_string(),
        ));
    }
    if query.hash == 0 || crate::jamhash_u64_v1(query.canonical_kmer) != query.hash {
        return Err(JmaError::CorruptSection(
            "seed query hash does not match its packed canonical k-mer".to_string(),
        ));
    }
    Ok(())
}

pub(crate) fn validate_occurrence(occurrence: SeedOccurrence, k: u8) -> JmaResult<()> {
    if occurrence.position > u64::MAX - u64::from(k) {
        return Err(JmaError::OffsetOverflow);
    }
    Ok(())
}

/// Check a static membership filter before requesting an exact key page.
#[must_use]
pub fn membership_contains(filter: &[u8], hash: u64) -> bool {
    let bit_count = filter.len().saturating_mul(8);
    bit_count != 0
        && builder::filter_bits(hash, bit_count)
            .into_iter()
            .all(|bit| filter[bit / 8] & (1u8 << (bit % 8)) != 0)
}
