//! Search adapter from the immutable router reader to posting codecs.

use crate::router::postings::{
    PositionPosting, PostingCodec, PostingError, PostingHeader, PostingSource,
    decode_document_ids_wire, decode_positions_wire,
};
use crate::router::reader::{KeyMatch, RawPosting, RouterReader};
use crate::router::{WitnessKey, WitnessScheme};
use std::collections::BTreeMap;
use std::sync::Mutex;

pub struct RouterPostingSource {
    reader: RouterReader,
    scheme: WitnessScheme,
    key_matches: Mutex<BTreeMap<WitnessKey, KeyMatch>>,
    postings: Mutex<BTreeMap<WitnessKey, RawPosting>>,
}

impl RouterPostingSource {
    pub fn new(reader: RouterReader, scheme: WitnessScheme) -> Result<Self, PostingError> {
        scheme
            .validate()
            .map_err(|error| PostingError::Backend(error.to_string()))?;
        Ok(Self {
            reader,
            scheme,
            key_matches: Mutex::new(BTreeMap::new()),
            postings: Mutex::new(BTreeMap::new()),
        })
    }

    #[must_use]
    pub fn reader(&self) -> &RouterReader {
        &self.reader
    }

    #[must_use]
    pub fn scheme(&self) -> &WitnessScheme {
        &self.scheme
    }

    fn key_match(&self, key: &WitnessKey) -> Result<Option<KeyMatch>, PostingError> {
        if let Some(found) = self
            .key_matches
            .lock()
            .map_err(|_| PostingError::Backend("router key cache is poisoned".to_string()))?
            .get(key)
            .copied()
        {
            return Ok(Some(found));
        }
        let found = self
            .reader
            .lookup(*key)
            .map_err(|error| PostingError::Backend(error.to_string()))?;
        if let Some(found) = found {
            self.key_matches
                .lock()
                .map_err(|_| PostingError::Backend("router key cache is poisoned".to_string()))?
                .insert(*key, found);
        }
        Ok(found)
    }

    fn raw_posting(&self, key: &WitnessKey) -> Result<RawPosting, PostingError> {
        if let Some(posting) = self
            .postings
            .lock()
            .map_err(|_| PostingError::Backend("router posting cache is poisoned".to_string()))?
            .get(key)
            .cloned()
        {
            return Ok(posting);
        }
        let key_match = self
            .key_match(key)?
            .ok_or(PostingError::MissingPosting { tier: 0, key: *key })?;
        let posting = self
            .reader
            .read_posting(key_match)
            .map_err(|error| PostingError::Backend(error.to_string()))?;
        self.postings
            .lock()
            .map_err(|_| PostingError::Backend("router posting cache is poisoned".to_string()))?
            .insert(*key, posting.clone());
        Ok(posting)
    }

    fn tier_available(&self, tier: u32, key: &WitnessKey) -> Result<bool, PostingError> {
        self.scheme
            .includes_hash(key.jamhash, tier)
            .map_err(|error| PostingError::Backend(error.to_string()))
    }
}

impl PostingSource for RouterPostingSource {
    fn collection_size(&self) -> u32 {
        u32::try_from(self.reader.metagenomes().len()).unwrap_or(u32::MAX)
    }

    fn metagenome_id(&self, document_id: u32) -> Option<&str> {
        self.reader
            .metagenomes()
            .get(usize::try_from(document_id).ok()?)
            .map(|entry| entry.id.as_str())
    }

    fn header(&self, tier: u32, key: &WitnessKey) -> Option<PostingHeader> {
        self.try_header(tier, key).ok().flatten()
    }

    fn try_header(
        &self,
        tier: u32,
        key: &WitnessKey,
    ) -> Result<Option<PostingHeader>, PostingError> {
        if !self.tier_available(tier, key)? {
            return Ok(None);
        }
        let Some(found) = self.key_match(key)? else {
            return Ok(None);
        };
        let header = self
            .reader
            .posting_header(found.key_index)
            .map_err(|error| PostingError::Backend(error.to_string()))?;
        Ok(Some(PostingHeader {
            scheme_id: self.scheme.scheme_id,
            document_frequency: header.document_frequency,
            posting_count: header.posting_count,
            payload_byte_length: u32::try_from(header.posting_length)
                .map_err(|_| PostingError::SizeOverflow)?,
            common_or_repetitive: header.flags & crate::router::format::PostingHeader::FLAG_COMMON
                != 0,
            suppressed: header.flags & crate::router::format::PostingHeader::FLAG_SUPPRESSED != 0,
            position_bearing: header.has_positions(),
            position_count: u32::try_from(header.position_length / 17).unwrap_or(u32::MAX),
            position_payload_byte_length: u32::try_from(header.position_length)
                .map_err(|_| PostingError::SizeOverflow)?,
            position_capped: false,
            codec: if header.posting_count <= 3 {
                PostingCodec::Inline
            } else {
                PostingCodec::DeltaVarByte
            },
            block_count: u32::from(header.posting_count > 3),
        }))
    }

    fn decode_document_ids(&self, tier: u32, key: &WitnessKey) -> Result<Vec<u32>, PostingError> {
        if !self.tier_available(tier, key)? {
            return Err(PostingError::MissingPosting { tier, key: *key });
        }
        decode_document_ids_wire(&self.raw_posting(key)?.posting)
    }

    fn decode_positions(
        &self,
        tier: u32,
        key: &WitnessKey,
        document_id: u32,
        max_positions: usize,
    ) -> Result<Vec<PositionPosting>, PostingError> {
        if !self.tier_available(tier, key)? || max_positions == 0 {
            return Ok(Vec::new());
        }
        let Some(bytes) = self.raw_posting(key)?.positions else {
            return Ok(Vec::new());
        };
        Ok(decode_positions_wire(&bytes)?
            .into_iter()
            .filter(|position| position.metagenome_id == document_id)
            .take(max_positions)
            .collect())
    }
}
