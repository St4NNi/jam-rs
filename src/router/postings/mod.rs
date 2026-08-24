//! Immutable witness posting codecs and a small lookup adapter.
//!
//! The on-disk router reader can implement [`PostingSource`] without using
//! the in-memory implementation in this module.  Keeping the search-facing
//! contract here makes it possible to exercise routing and handoff without
//! opening a JMA object or depending on a particular page layout.

use crate::router::{RouterContractError, WitnessKey};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use thiserror::Error;

/// IDs up to this count are kept inline in a posting header.
pub const INLINE_POSTING_LIMIT: usize = 3;
/// Independent delta blocks bound the amount of work for a single lookup.
pub const DEFAULT_POSTING_BLOCK_SIZE: usize = 128;

/// Codec used for the document-ID portion of a posting.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PostingCodec {
    /// Zero to three sorted IDs are stored directly in the record.
    Inline,
    /// Each bounded block stores an absolute first ID followed by varint
    /// deltas.  Blocks can be decoded independently by a range reader.
    DeltaVarByte,
}

/// Metadata which can be inspected before decoding any posting payload.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PostingHeader {
    pub scheme_id: u32,
    pub document_frequency: u32,
    pub posting_count: u32,
    pub payload_byte_length: u32,
    pub common_or_repetitive: bool,
    pub suppressed: bool,
    pub position_bearing: bool,
    pub position_count: u32,
    pub position_payload_byte_length: u32,
    pub position_capped: bool,
    pub codec: PostingCodec,
    pub block_count: u32,
}

impl PostingHeader {
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.posting_count == 0
    }
}

/// One optional position-bearing occurrence in a collection posting.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct PositionPosting {
    pub metagenome_id: u32,
    pub contig_id: u32,
    pub position: u64,
    pub reverse: bool,
}

/// Input used by [`InMemoryPostingSource::insert`].
#[derive(Clone, Debug, Default)]
pub struct PostingInput {
    pub tier: u32,
    pub scheme_id: u32,
    pub key: Option<WitnessKey>,
    pub document_ids: Vec<u32>,
    pub positions: Vec<PositionPosting>,
    pub position_cap_per_document: Option<usize>,
    pub common_or_repetitive: bool,
    pub suppressed: bool,
}

impl PostingInput {
    #[must_use]
    pub fn new(tier: u32, key: WitnessKey, document_ids: Vec<u32>) -> Self {
        Self {
            tier,
            scheme_id: 1,
            key: Some(key),
            document_ids,
            ..Self::default()
        }
    }

    #[must_use]
    pub fn with_positions(mut self, positions: Vec<PositionPosting>) -> Self {
        self.positions = positions;
        self
    }

    #[must_use]
    pub const fn with_scheme_id(mut self, scheme_id: u32) -> Self {
        self.scheme_id = scheme_id;
        self
    }

    #[must_use]
    pub const fn with_position_cap(mut self, cap_per_document: usize) -> Self {
        self.position_cap_per_document = Some(cap_per_document);
        self
    }

    #[must_use]
    pub const fn mark_common(mut self) -> Self {
        self.common_or_repetitive = true;
        self
    }

    #[must_use]
    pub const fn mark_suppressed(mut self) -> Self {
        self.suppressed = true;
        self
    }
}

/// Search-facing access to one immutable posting source.
pub trait PostingSource {
    /// Number of documents represented by the collection catalog.
    fn collection_size(&self) -> u32;

    /// Resolve a compact document ID to its stable metagenome identifier.
    fn metagenome_id(&self, document_id: u32) -> Option<&str>;

    /// Read only the header.  A miss or a header-level suppression must not
    /// require decoding either payload.
    fn header(&self, tier: u32, key: &WitnessKey) -> Option<PostingHeader>;

    /// Decode sorted document IDs after the caller has accepted the header.
    fn decode_document_ids(&self, tier: u32, key: &WitnessKey) -> Result<Vec<u32>, PostingError>;

    /// Visit sorted IDs without requiring a caller to collect a complete
    /// posting.  A range-backed implementation should override this method
    /// and decode one bounded block at a time; the default preserves the
    /// adapter contract for simple readers.
    fn for_each_document_id(
        &self,
        tier: u32,
        key: &WitnessKey,
        visitor: &mut dyn FnMut(u32) -> Result<(), PostingError>,
    ) -> Result<(), PostingError> {
        for document_id in self.decode_document_ids(tier, key)? {
            visitor(document_id)?;
        }
        Ok(())
    }

    /// Decode position-bearing occurrences for one selected document.
    fn decode_positions(
        &self,
        tier: u32,
        key: &WitnessKey,
        document_id: u32,
        max_positions: usize,
    ) -> Result<Vec<PositionPosting>, PostingError>;
}

/// Payload decode counters useful for proving that common-key suppression is
/// header-first.  A disk reader may expose equivalent counters independently.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PostingDecodeStats {
    pub document_payload_decodes: u64,
    pub position_payload_decodes: u64,
}

#[derive(Debug, Default)]
struct DecodeCounters {
    document_payload_decodes: AtomicU64,
    position_payload_decodes: AtomicU64,
}

impl DecodeCounters {
    fn snapshot(&self) -> PostingDecodeStats {
        PostingDecodeStats {
            document_payload_decodes: self.document_payload_decodes.load(Ordering::Relaxed),
            position_payload_decodes: self.position_payload_decodes.load(Ordering::Relaxed),
        }
    }
}

#[derive(Clone, Debug)]
struct DeltaBlock {
    first_id: u32,
    count: u16,
    payload: Vec<u8>,
}

#[derive(Clone, Debug)]
enum DocumentPayload {
    Inline(Vec<u32>),
    Blocks(Vec<DeltaBlock>),
}

#[derive(Clone, Debug)]
struct EncodedPosting {
    key: WitnessKey,
    tier: u32,
    header: PostingHeader,
    documents: DocumentPayload,
    positions: Vec<PositionPosting>,
}

/// A deterministic, local posting source used by router tests and adapters.
///
/// The production reader only needs to implement [`PostingSource`].  This
/// type intentionally stores each exact `(tier, packed kmer, jamhash)` key;
/// looking up a matching digest with a different packed kmer is therefore a
/// miss and can never create candidate evidence.
#[derive(Clone, Debug)]
pub struct InMemoryPostingSource {
    document_names: Vec<String>,
    postings: BTreeMap<(u32, WitnessKey), EncodedPosting>,
    counters: Arc<DecodeCounters>,
}

impl InMemoryPostingSource {
    #[must_use]
    pub fn new(document_names: Vec<String>) -> Self {
        Self {
            document_names,
            postings: BTreeMap::new(),
            counters: Arc::new(DecodeCounters::default()),
        }
    }

    /// Insert one posting, canonicalizing IDs and sorting occurrences.
    pub fn insert(&mut self, mut input: PostingInput) -> Result<(), PostingError> {
        let key = input.key.take().ok_or(PostingError::MissingKey)?;
        let key = WitnessKey::checked(key.packed, key.jamhash)?;
        if input.tier == 0 || input.scheme_id == 0 {
            return Err(PostingError::InvalidTier);
        }

        let mut document_ids = input.document_ids;
        for occurrence in &input.positions {
            document_ids.push(occurrence.metagenome_id);
        }
        document_ids.sort_unstable();
        document_ids.dedup();
        if let Some(document_id) = document_ids
            .iter()
            .copied()
            .find(|id| (*id as usize) >= self.document_names.len())
        {
            return Err(PostingError::UnknownDocument(document_id));
        }

        let documents = encode_documents(&document_ids);
        let position_result = encode_positions(&input.positions, input.position_cap_per_document);
        let payload_byte_length = document_payload_len(&documents)?;
        let position_payload_byte_length = position_result
            .occurrences
            .len()
            .checked_mul(POSITION_RECORD_BYTES)
            .and_then(|value| u32::try_from(value).ok())
            .ok_or(PostingError::SizeOverflow)?;
        let header = PostingHeader {
            scheme_id: input.scheme_id,
            document_frequency: u32::try_from(document_ids.len())
                .map_err(|_| PostingError::SizeOverflow)?,
            posting_count: u32::try_from(document_ids.len())
                .map_err(|_| PostingError::SizeOverflow)?,
            payload_byte_length,
            common_or_repetitive: input.common_or_repetitive,
            suppressed: input.suppressed,
            position_bearing: !position_result.occurrences.is_empty(),
            position_count: u32::try_from(position_result.occurrences.len())
                .map_err(|_| PostingError::SizeOverflow)?,
            position_payload_byte_length,
            position_capped: position_result.capped,
            codec: match &documents {
                DocumentPayload::Inline(_) => PostingCodec::Inline,
                DocumentPayload::Blocks(_) => PostingCodec::DeltaVarByte,
            },
            block_count: match &documents {
                DocumentPayload::Inline(_) => 0,
                DocumentPayload::Blocks(blocks) => {
                    u32::try_from(blocks.len()).map_err(|_| PostingError::SizeOverflow)?
                }
            },
        };
        let encoded = EncodedPosting {
            key,
            tier: input.tier,
            header,
            documents,
            positions: position_result.occurrences,
        };
        if self.postings.insert((input.tier, key), encoded).is_some() {
            return Err(PostingError::DuplicatePosting {
                tier: input.tier,
                key,
            });
        }
        Ok(())
    }

    #[must_use]
    pub fn stats(&self) -> PostingDecodeStats {
        self.counters.snapshot()
    }

    pub fn reset_stats(&self) {
        self.counters
            .document_payload_decodes
            .store(0, Ordering::Relaxed);
        self.counters
            .position_payload_decodes
            .store(0, Ordering::Relaxed);
    }

    #[must_use]
    pub fn posting_count(&self) -> usize {
        self.postings.len()
    }

    #[must_use]
    pub fn document_names(&self) -> &[String] {
        &self.document_names
    }

    fn get(&self, tier: u32, key: &WitnessKey) -> Option<&EncodedPosting> {
        self.postings.get(&(tier, *key)).filter(|posting| {
            // Keep this explicit check next to the map lookup.  It protects
            // adapters which choose to index by digest only.
            posting.key == *key && posting.tier == tier
        })
    }
}

impl PostingSource for InMemoryPostingSource {
    fn collection_size(&self) -> u32 {
        u32::try_from(self.document_names.len()).unwrap_or(u32::MAX)
    }

    fn metagenome_id(&self, document_id: u32) -> Option<&str> {
        self.document_names
            .get(usize::try_from(document_id).ok()?)
            .map(String::as_str)
    }

    fn header(&self, tier: u32, key: &WitnessKey) -> Option<PostingHeader> {
        self.get(tier, key).map(|posting| posting.header)
    }

    fn decode_document_ids(&self, tier: u32, key: &WitnessKey) -> Result<Vec<u32>, PostingError> {
        let posting = self
            .get(tier, key)
            .ok_or(PostingError::MissingPosting { tier, key: *key })?;
        self.counters
            .document_payload_decodes
            .fetch_add(1, Ordering::Relaxed);
        decode_documents(&posting.documents)
    }

    fn for_each_document_id(
        &self,
        tier: u32,
        key: &WitnessKey,
        visitor: &mut dyn FnMut(u32) -> Result<(), PostingError>,
    ) -> Result<(), PostingError> {
        let posting = self
            .get(tier, key)
            .ok_or(PostingError::MissingPosting { tier, key: *key })?;
        self.counters
            .document_payload_decodes
            .fetch_add(1, Ordering::Relaxed);
        for_each_document(&posting.documents, visitor)
    }

    fn decode_positions(
        &self,
        tier: u32,
        key: &WitnessKey,
        document_id: u32,
        max_positions: usize,
    ) -> Result<Vec<PositionPosting>, PostingError> {
        let posting = self
            .get(tier, key)
            .ok_or(PostingError::MissingPosting { tier, key: *key })?;
        if !posting.header.position_bearing || max_positions == 0 {
            return Ok(Vec::new());
        }
        self.counters
            .position_payload_decodes
            .fetch_add(1, Ordering::Relaxed);
        Ok(posting
            .positions
            .iter()
            .copied()
            .filter(|occurrence| occurrence.metagenome_id == document_id)
            .take(max_positions)
            .collect())
    }
}

const POSITION_RECORD_BYTES: usize = 17;

#[derive(Clone, Debug)]
struct EncodedPositions {
    occurrences: Vec<PositionPosting>,
    capped: bool,
}

fn encode_positions(
    occurrences: &[PositionPosting],
    cap_per_document: Option<usize>,
) -> EncodedPositions {
    let mut sorted = occurrences.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    let Some(cap) = cap_per_document else {
        return EncodedPositions {
            occurrences: sorted,
            capped: false,
        };
    };

    let mut kept = Vec::with_capacity(sorted.len().min(cap));
    let mut counts = BTreeMap::<u32, usize>::new();
    let mut capped = false;
    for occurrence in sorted {
        let count = counts.entry(occurrence.metagenome_id).or_default();
        if *count < cap {
            kept.push(occurrence);
            *count += 1;
        } else {
            capped = true;
        }
    }
    EncodedPositions {
        occurrences: kept,
        capped,
    }
}

fn encode_documents(document_ids: &[u32]) -> DocumentPayload {
    if document_ids.len() <= INLINE_POSTING_LIMIT {
        return DocumentPayload::Inline(document_ids.to_vec());
    }
    let blocks = document_ids
        .chunks(DEFAULT_POSTING_BLOCK_SIZE)
        .map(|chunk| {
            let first_id = chunk[0];
            let mut payload = Vec::new();
            let mut previous = first_id;
            for &document_id in &chunk[1..] {
                encode_var_u32(document_id - previous, &mut payload);
                previous = document_id;
            }
            DeltaBlock {
                first_id,
                count: u16::try_from(chunk.len()).unwrap_or(u16::MAX),
                payload,
            }
        })
        .collect();
    DocumentPayload::Blocks(blocks)
}

fn document_payload_len(payload: &DocumentPayload) -> Result<u32, PostingError> {
    let length = match payload {
        DocumentPayload::Inline(_) => 0,
        DocumentPayload::Blocks(blocks) => blocks.iter().try_fold(0_usize, |total, block| {
            total
                .checked_add(std::mem::size_of::<u32>())
                .and_then(|sum| sum.checked_add(block.payload.len()))
                .ok_or(PostingError::SizeOverflow)
        })?,
    };
    u32::try_from(length).map_err(|_| PostingError::SizeOverflow)
}

fn decode_documents(payload: &DocumentPayload) -> Result<Vec<u32>, PostingError> {
    let mut result = Vec::new();
    for_each_document(payload, &mut |document_id| {
        result.push(document_id);
        Ok(())
    })?;
    Ok(result)
}

fn for_each_document(
    payload: &DocumentPayload,
    visitor: &mut dyn FnMut(u32) -> Result<(), PostingError>,
) -> Result<(), PostingError> {
    match payload {
        DocumentPayload::Inline(ids) => {
            for &document_id in ids {
                visitor(document_id)?;
            }
        }
        DocumentPayload::Blocks(blocks) => {
            for block in blocks {
                if block.count == 0 {
                    return Err(PostingError::CorruptPayload(
                        "zero-length delta block".to_string(),
                    ));
                }
                visitor(block.first_id)?;
                let mut previous = block.first_id;
                let mut offset = 0;
                for _ in 1..usize::from(block.count) {
                    let delta = decode_var_u32(&block.payload, &mut offset)?;
                    previous = previous
                        .checked_add(delta)
                        .ok_or(PostingError::DeltaOverflow)?;
                    visitor(previous)?;
                }
                if offset != block.payload.len() {
                    return Err(PostingError::CorruptPayload(
                        "trailing bytes in delta block".to_string(),
                    ));
                }
            }
        }
    }
    Ok(())
}

fn encode_var_u32(mut value: u32, output: &mut Vec<u8>) {
    while value >= 0x80 {
        output.push((value as u8 & 0x7f) | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
}

fn decode_var_u32(payload: &[u8], offset: &mut usize) -> Result<u32, PostingError> {
    let mut value = 0_u32;
    let mut shift = 0_u32;
    loop {
        let byte = *payload
            .get(*offset)
            .ok_or_else(|| PostingError::CorruptPayload("truncated varint".to_string()))?;
        *offset += 1;
        value |= u32::from(byte & 0x7f)
            .checked_shl(shift)
            .ok_or(PostingError::DeltaOverflow)?;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
        shift = shift.checked_add(7).ok_or(PostingError::DeltaOverflow)?;
        if shift >= 32 {
            return Err(PostingError::DeltaOverflow);
        }
    }
}

/// Errors raised while constructing or decoding a posting.
#[derive(Debug, Error, Eq, PartialEq)]
pub enum PostingError {
    #[error("posting key is required")]
    MissingKey,
    #[error("posting tier must be greater than zero")]
    InvalidTier,
    #[error("document {0} is not present in the catalog")]
    UnknownDocument(u32),
    #[error("posting already exists for tier {tier} and key {key:?}")]
    DuplicatePosting { tier: u32, key: WitnessKey },
    #[error("posting is missing for tier {tier} and key {key:?}")]
    MissingPosting { tier: u32, key: WitnessKey },
    #[error("posting payload is corrupt: {0}")]
    CorruptPayload(String),
    #[error("posting delta overflow")]
    DeltaOverflow,
    #[error("posting is too large to represent")]
    SizeOverflow,
    #[error(transparent)]
    Contract(#[from] RouterContractError),
}
