use crate::jidx::{seed_length, sha256};
use crate::jidx_reader::{Contig, JidxReader, Metagenome, SeedDocument, SeedEntry, SeedOccurrence};
use crate::owner_format::{OwnerDocument, OwnerSeed};
use crate::owner_reader::OwnerReader;
use crate::shared_reader::{SharedGroup, SharedMember, SharedReader};
use crate::shared_seed::SharedKey;
use crate::trace::TraceError;

pub(crate) enum TraceIndex {
    Shard(Box<JidxReader>),
    Owner(OwnerReader),
    Shared(Box<SharedReader>),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceCacheIdentity {
    Shard([u64; 7]),
    Owner([u8; 32]),
    Shared([u64; 7], u64),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceSeed {
    Ordinary(SeedEntry),
    Shared(SharedGroup),
}

impl TraceSeed {
    pub(crate) fn packed_key(self) -> u64 {
        match self {
            Self::Ordinary(seed) => seed.packed_key,
            Self::Shared(group) => group.key().packed().unwrap(),
        }
    }

    pub(crate) fn document_frequency(self) -> u32 {
        match self {
            Self::Ordinary(seed) => seed.document_frequency,
            Self::Shared(group) => group.member_count(),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TraceDocument {
    Shard(SeedDocument),
    Owner(OwnerDocument),
    Shared(SharedMember),
}

impl TraceDocument {
    pub(crate) fn metagenome_id(self) -> u32 {
        match self {
            Self::Shard(document) => document.metagenome_id,
            Self::Owner(document) => document.metagenome_id,
            Self::Shared(member) => member.metagenome_id,
        }
    }

    pub(crate) fn occurrence_count(self) -> u64 {
        match self {
            Self::Shard(document) => document.occurrence_count,
            Self::Owner(document) => document.occurrence_count,
            Self::Shared(member) => member.occurrence_count(),
        }
    }
}

impl TraceIndex {
    pub(crate) fn shard(&self) -> Result<&JidxReader, TraceError> {
        match self {
            Self::Shard(index) => Ok(index),
            Self::Owner(_) | Self::Shared(_) => {
                Err(TraceError::Invalid("expected a metagenome shard"))
            }
        }
    }

    pub(crate) fn k(&self) -> u8 {
        match self {
            Self::Shard(index) => index.header().k,
            Self::Owner(index) => index.k(),
            Self::Shared(_) => 15,
        }
    }

    pub(crate) fn rescue_k15(&self) -> bool {
        match self {
            Self::Shard(index) => index.header().rescue_k15,
            Self::Owner(index) => index.rescue_k15(),
            Self::Shared(_) => false,
        }
    }

    pub(crate) fn document_count(&self) -> u32 {
        match self {
            Self::Shard(index) => index.header().document_count,
            Self::Owner(index) => index.document_count(),
            Self::Shared(index) => index.document_count(),
        }
    }

    pub(crate) fn header_sha256(&self) -> Result<[u8; 32], TraceError> {
        match self {
            Self::Shard(index) => Ok(sha256(&index.header().encode()?)),
            Self::Owner(index) => Ok(index.header_sha256()),
            Self::Shared(index) => Ok(index.header_sha256()?),
        }
    }

    pub(crate) fn manifest_sha256(&self) -> [u8; 32] {
        match self {
            Self::Shard(index) => index.header().manifest_sha256,
            Self::Owner(index) => index.manifest_sha256(),
            Self::Shared(index) => index.manifest_sha256(),
        }
    }

    pub(crate) fn body_sha256(&self) -> [u8; 32] {
        match self {
            Self::Shard(index) => index.header().body_sha256,
            Self::Owner(index) => index.body_sha256(),
            Self::Shared(index) => index.body_sha256(),
        }
    }

    pub(crate) fn cache_file_identity(&self) -> std::io::Result<Option<TraceCacheIdentity>> {
        match self {
            Self::Shard(index) => Ok(index.cache_file_identity()?.map(TraceCacheIdentity::Shard)),
            Self::Owner(index) => Ok(index.cache_identity()?.map(TraceCacheIdentity::Owner)),
            Self::Shared(index) => Ok(index
                .cache_identity()
                .map_err(std::io::Error::other)?
                .map(|identity| TraceCacheIdentity::Shared(identity, index.reader_token()))),
        }
    }

    pub(crate) fn verify_checksum(&self) -> Result<(), TraceError> {
        match self {
            Self::Shard(index) => Ok(index.verify_checksum()?),
            Self::Owner(index) => Ok(index.verify_checksum()?),
            Self::Shared(index) => Ok(index.verify_checksum()?),
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
            Self::Shared(_) => {
                for &key in keys {
                    shared_key(key)?;
                }
                Ok(())
            }
        }
    }

    pub(crate) fn find_seeds_batch(
        &self,
        keys: &[u64],
    ) -> Result<Vec<Option<TraceSeed>>, TraceError> {
        match self {
            Self::Shard(index) => {
                let seeds = index.find_seeds_batch(keys)?;
                index.advise_first_document_rows(&seeds);
                Ok(seeds
                    .into_iter()
                    .map(|seed| seed.map(TraceSeed::Ordinary))
                    .collect())
            }
            Self::Shared(index) => {
                let contexts = keys
                    .iter()
                    .map(|&key| shared_key(key))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(index
                    .find_many(&contexts)?
                    .into_iter()
                    .map(|group| group.map(TraceSeed::Shared))
                    .collect())
            }
            Self::Owner(index) => index
                .find_seeds_batch(keys)?
                .into_iter()
                .map(|seed| {
                    seed.map(|seed| {
                        Ok(TraceSeed::Ordinary(SeedEntry {
                            packed_key: seed.packed_key,
                            document_frequency: u32::try_from(seed.document_frequency)
                                .map_err(|_| TraceError::Invalid("document frequency"))?,
                            document_offset: seed.block_ordinal,
                        }))
                    })
                    .transpose()
                })
                .collect(),
        }
    }

    #[cfg(test)]
    pub(crate) fn find_seed(&self, key: u64) -> Result<Option<TraceSeed>, TraceError> {
        self.find_seeds_batch(&[key])
            .map(|mut seeds| seeds.remove(0))
    }

    pub(crate) fn find_seeds_in_cores(
        &self,
        keys: &[u64],
        cores: &crate::trace_batch::SharedCoreLookups,
    ) -> Result<Vec<Option<TraceSeed>>, TraceError> {
        let Self::Shared(reader) = self else {
            return Err(TraceError::Invalid("resolved core index"));
        };
        if self.cache_file_identity()? != Some(cores.identity) {
            return Err(TraceError::Invalid("resolved core identity"));
        }
        let mut requests = Vec::new();
        requests
            .try_reserve_exact(keys.len())
            .map_err(|_| crate::shared_format::SharedError::ResourceLimit)?;
        for (ordinal, &key) in keys.iter().enumerate() {
            requests.push((ordinal, shared_key(key)?));
        }
        let order =
            |&(ordinal, key): &(usize, SharedKey)| (key.core, key.context_code().unwrap(), ordinal);
        if !requests
            .windows(2)
            .all(|pair| order(&pair[0]) <= order(&pair[1]))
        {
            requests.sort_unstable_by_key(order);
        }
        let mut result = Vec::new();
        result
            .try_reserve_exact(keys.len())
            .map_err(|_| crate::shared_format::SharedError::ResourceLimit)?;
        result.resize(keys.len(), None);
        let mut contexts = Vec::new();
        contexts
            .try_reserve_exact(keys.len())
            .map_err(|_| crate::shared_format::SharedError::ResourceLimit)?;
        let mut groups = Vec::new();
        groups
            .try_reserve_exact(keys.len())
            .map_err(|_| crate::shared_format::SharedError::ResourceLimit)?;
        let operation = reader.posting_operation()?;
        let mut core_ordinal = 0;
        for same_core in requests.chunk_by(|left, right| left.1.core == right.1.core) {
            let core_key = same_core[0].1.core;
            while cores
                .groups
                .get(core_ordinal)
                .is_some_and(|group| group.key().core < core_key)
            {
                core_ordinal += 1;
            }
            let Some(&core) = cores
                .groups
                .get(core_ordinal)
                .filter(|group| group.key().core == core_key)
            else {
                continue;
            };
            contexts.clear();
            for &(ordinal, key) in same_core {
                if key.length == 15 {
                    result[ordinal] = Some(TraceSeed::Shared(core));
                } else {
                    contexts.push(key);
                }
            }
            groups.resize(contexts.len(), None);
            operation.find_in_core_into(core, &contexts, &mut groups)?;
            if !contexts.is_empty() {
                for ((ordinal, _), group) in same_core
                    .iter()
                    .filter(|(_, key)| key.length != 15)
                    .zip(groups.iter().copied())
                {
                    result[*ordinal] = group.map(TraceSeed::Shared);
                }
            }
        }
        operation.finish()?;
        Ok(result)
    }

    pub(crate) fn seed_documents(&self, seed: TraceSeed) -> Result<Vec<TraceDocument>, TraceError> {
        match (self, seed) {
            (Self::Shard(index), TraceSeed::Ordinary(seed)) => Ok(index
                .seed_documents(seed)?
                .into_iter()
                .map(TraceDocument::Shard)
                .collect()),
            (Self::Owner(index), TraceSeed::Ordinary(seed)) => Ok(index
                .seed_documents(owner_seed(seed))?
                .into_iter()
                .map(TraceDocument::Owner)
                .collect()),
            (Self::Shared(index), TraceSeed::Shared(group)) => Ok(index
                .members(group)?
                .into_iter()
                .map(TraceDocument::Shared)
                .collect()),
            _ => Err(TraceError::Invalid("seed index")),
        }
    }

    pub(crate) fn seed_document_occurrences(
        &self,
        seed: TraceSeed,
        document: TraceDocument,
    ) -> Result<Vec<SeedOccurrence>, TraceError> {
        match (self, seed, document) {
            (Self::Shard(index), TraceSeed::Ordinary(seed), TraceDocument::Shard(document)) => {
                Ok(index.seed_document_occurrences(seed, document)?)
            }
            (Self::Owner(index), TraceSeed::Ordinary(seed), TraceDocument::Owner(document)) => {
                Ok(index.seed_document_occurrences(owner_seed(seed), document)?)
            }
            (Self::Shared(index), TraceSeed::Shared(group), TraceDocument::Shared(member)) => {
                Ok(index.member_occurrences(group, member)?)
            }
            _ => Err(TraceError::Invalid("seed document index")),
        }
    }

    pub(crate) fn metagenome_name(&self, id: u32) -> Result<Option<&str>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.metagenome_name(id)?),
            Self::Owner(index) => Ok(index.metagenome_name(id)?),
            Self::Shared(index) => Ok(index.metagenome(id)?.map(|source| source.name)),
        }
    }

    pub(crate) fn visit_occurrences(
        &self,
        seed: TraceSeed,
        document: TraceDocument,
        mut visit: impl FnMut(&[SeedOccurrence]) -> Result<(), TraceError>,
    ) -> Result<(), TraceError> {
        if let (Self::Shared(index), TraceSeed::Shared(group), TraceDocument::Shared(member)) =
            (self, seed, document)
        {
            let occurrence_count = member.occurrence_count();
            let mut start = 0;
            while start < occurrence_count {
                let block = index.occurrence_block(group, member, start, 4096)?;
                if block.is_empty() {
                    return Err(TraceError::Invalid("shared occurrence progress"));
                }
                visit(&block)?;
                start += block.len() as u64;
            }
            Ok(())
        } else {
            visit(&self.seed_document_occurrences(seed, document)?)
        }
    }

    pub(crate) fn metagenome(&self, id: u32) -> Result<Option<Metagenome<'_>>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.metagenome(id)?),
            Self::Owner(index) => Ok(index.metagenome(id)?),
            Self::Shared(index) => Ok(index.metagenome(id)?),
        }
    }

    pub(crate) fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, TraceError> {
        match self {
            Self::Shard(index) => Ok(index.contig(id)?),
            Self::Owner(index) => Ok(index.contig(id)?),
            Self::Shared(index) => Ok(index.contig(id)?),
        }
    }

    pub(crate) fn numeric_contig(
        &self,
        id: u32,
    ) -> Result<Option<crate::shared_reader::NumericContig>, TraceError> {
        if let Self::Shared(index) = self {
            return Ok(index.numeric_contig(id)?);
        }
        Ok(self
            .contig(id)?
            .map(|contig| crate::shared_reader::NumericContig {
                id: contig.id,
                metagenome_id: contig.metagenome_id,
                length: contig.length,
            }))
    }

    pub(crate) fn enable_selected_front_metadata(&self) {
        if let Self::Shard(index) = self {
            index.enable_selected_front_metadata();
        }
    }

    pub(crate) fn is_shared(&self) -> bool {
        matches!(self, Self::Shared(_))
    }

    pub(crate) fn seed_length(&self, key: u64) -> Result<u8, TraceError> {
        if self.is_shared() {
            Ok(shared_key(key)?.length)
        } else {
            Ok(seed_length(self.k(), self.rescue_k15(), key)?)
        }
    }
}

fn shared_key(key: u64) -> Result<SharedKey, TraceError> {
    SharedKey::unpack(key).ok_or(TraceError::Invalid("shared context key"))
}

fn owner_seed(seed: SeedEntry) -> OwnerSeed {
    OwnerSeed {
        packed_key: seed.packed_key,
        document_frequency: u64::from(seed.document_frequency),
        block_ordinal: seed.document_offset,
    }
}
