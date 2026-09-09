use crate::jidx::{seed_length, sha256};
use crate::jidx_reader::{Contig, JidxReader, Metagenome, SeedDocument, SeedEntry, SeedOccurrence};
use crate::owner_format::{OwnerDocument, OwnerSeed};
use crate::owner_reader::OwnerReader;
use crate::trace::TraceError;

pub(crate) enum TraceIndex {
    Shard(Box<JidxReader>),
    Owner(OwnerReader),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceCacheIdentity {
    Shard([u64; 7]),
    Owner([u8; 32]),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceDocument {
    Shard(SeedDocument),
    Owner(OwnerDocument),
}

impl TraceDocument {
    pub(crate) fn metagenome_id(self) -> u32 {
        match self {
            Self::Shard(document) => document.metagenome_id,
            Self::Owner(document) => document.metagenome_id,
        }
    }

    pub(crate) fn occurrence_count(self) -> u64 {
        match self {
            Self::Shard(document) => document.occurrence_count,
            Self::Owner(document) => document.occurrence_count,
        }
    }
}

impl TraceIndex {
    pub(crate) fn shard(&self) -> Result<&JidxReader, TraceError> {
        match self {
            Self::Shard(index) => Ok(index),
            Self::Owner(_) => Err(TraceError::Invalid("expected a metagenome shard")),
        }
    }

    pub(crate) fn k(&self) -> u8 {
        match self {
            Self::Shard(index) => index.header().k,
            Self::Owner(index) => index.k(),
        }
    }

    pub(crate) fn rescue_k15(&self) -> bool {
        match self {
            Self::Shard(index) => index.header().rescue_k15,
            Self::Owner(index) => index.rescue_k15(),
        }
    }

    pub(crate) fn document_count(&self) -> u32 {
        match self {
            Self::Shard(index) => index.header().document_count,
            Self::Owner(index) => index.document_count(),
        }
    }

    pub(crate) fn header_sha256(&self) -> Result<[u8; 32], TraceError> {
        match self {
            Self::Shard(index) => Ok(sha256(&index.header().encode()?)),
            Self::Owner(index) => Ok(index.header_sha256()),
        }
    }

    pub(crate) fn manifest_sha256(&self) -> [u8; 32] {
        match self {
            Self::Shard(index) => index.header().manifest_sha256,
            Self::Owner(index) => index.manifest_sha256(),
        }
    }

    pub(crate) fn body_sha256(&self) -> [u8; 32] {
        match self {
            Self::Shard(index) => index.header().body_sha256,
            Self::Owner(index) => index.body_sha256(),
        }
    }

    pub(crate) fn cache_file_identity(&self) -> std::io::Result<Option<TraceCacheIdentity>> {
        match self {
            Self::Shard(index) => Ok(index.cache_file_identity()?.map(TraceCacheIdentity::Shard)),
            Self::Owner(index) => Ok(index.cache_identity()?.map(TraceCacheIdentity::Owner)),
        }
    }

    pub(crate) fn verify_checksum(&self) -> Result<(), TraceError> {
        match self {
            Self::Shard(index) => Ok(index.verify_checksum()?),
            Self::Owner(index) => Ok(index.verify_checksum()?),
        }
    }

    pub(crate) fn verify_query_filter_pages<'a>(
        &self,
        keys: impl IntoIterator<Item = &'a u64>,
    ) -> Result<(), TraceError> {
        match self {
            Self::Shard(index) => Ok(index.verify_query_filter_pages(keys)?),
            Self::Owner(index) => {
                for &key in keys {
                    seed_length(index.k(), index.rescue_k15(), key)?;
                }
                Ok(())
            }
        }
    }

    pub(crate) fn find_seeds_batch(
        &self,
        keys: &[u64],
    ) -> Result<Vec<Option<SeedEntry>>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.find_seeds_batch(keys)?),
            Self::Owner(index) => index
                .find_seeds_batch(keys)?
                .into_iter()
                .map(|seed| {
                    seed.map(|seed| {
                        Ok(SeedEntry {
                            packed_key: seed.packed_key,
                            document_frequency: u32::try_from(seed.document_frequency)
                                .map_err(|_| TraceError::Invalid("document frequency"))?,
                            document_offset: seed.block_ordinal,
                        })
                    })
                    .transpose()
                })
                .collect(),
        }
    }

    #[cfg(test)]
    pub(crate) fn find_seed(&self, key: u64) -> Result<Option<SeedEntry>, TraceError> {
        self.find_seeds_batch(&[key])
            .map(|mut seeds| seeds.remove(0))
    }

    pub(crate) fn seed_documents(&self, seed: SeedEntry) -> Result<Vec<TraceDocument>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index
                .seed_documents(seed)?
                .into_iter()
                .map(TraceDocument::Shard)
                .collect()),
            Self::Owner(index) => Ok(index
                .seed_documents(owner_seed(seed))?
                .into_iter()
                .map(TraceDocument::Owner)
                .collect()),
        }
    }

    pub(crate) fn seed_document_occurrences(
        &self,
        seed: SeedEntry,
        document: TraceDocument,
    ) -> Result<Vec<SeedOccurrence>, TraceError> {
        match (self, document) {
            (Self::Shard(index), TraceDocument::Shard(document)) => {
                Ok(index.seed_document_occurrences(seed, document)?)
            }
            (Self::Owner(index), TraceDocument::Owner(document)) => {
                Ok(index.seed_document_occurrences(owner_seed(seed), document)?)
            }
            _ => Err(TraceError::Invalid("seed document index")),
        }
    }

    pub(crate) fn metagenome_name(&self, id: u32) -> Result<Option<&str>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.metagenome_name(id)?),
            Self::Owner(index) => Ok(index.metagenome_name(id)?),
        }
    }

    pub(crate) fn metagenome(&self, id: u32) -> Result<Option<Metagenome<'_>>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.metagenome(id)?),
            Self::Owner(index) => Ok(index.metagenome(id)?),
        }
    }

    pub(crate) fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.contig(id)?),
            Self::Owner(index) => Ok(index.contig(id)?),
        }
    }

    pub(crate) fn advise_first_document_rows(&self, seeds: &[Option<SeedEntry>]) {
        if let Self::Shard(index) = self {
            index.advise_first_document_rows(seeds);
        }
    }

    pub(crate) fn enable_selected_front_metadata(&self) {
        if let Self::Shard(index) = self {
            index.enable_selected_front_metadata();
        }
    }
}

fn owner_seed(seed: SeedEntry) -> OwnerSeed {
    OwnerSeed {
        packed_key: seed.packed_key,
        document_frequency: u64::from(seed.document_frequency),
        block_ordinal: seed.document_offset,
    }
}
