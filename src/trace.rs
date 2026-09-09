use crate::alignment::{
    AlignmentConfig, AlignmentError, Interval, Strand, TraceAlignmentWorkspace,
};
use crate::bgzf::{BgzfError, BgzfReader};
use crate::jidx::{JidxError, RESCUE_K15_TAG, seed_length, sha256, sha256_reader};
use crate::jidx_reader::{
    ContigId, JidxReader, JidxReaderError, MetagenomeId, SEED_LOOKUP_BATCH_KEYS, SeedDocument,
    SeedEntry,
};
use crate::mosaic::{Fragment, Mosaic, MosaicError, build_mosaic};
use crate::query::{QueryEngine, QueryError, QuerySketch};
use crate::range_source::S3Config;
use crate::reader::ReaderError;
use needletail::Sequence;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceConfig {
    pub min_containment: f64,
    pub max_metagenomes: usize,
    pub use_sketch: bool,
    pub min_seed_hits: u32,
    pub diagonal_bin_bases: u64,
    pub flank_bases: u64,
    pub min_identity: f64,
    pub min_aligned_bases: u64,
    pub endpoint_bases: usize,
    pub circular: bool,
    pub verify_resources: bool,
    pub alignment: AlignmentConfig,
}

impl Default for TraceConfig {
    fn default() -> Self {
        Self {
            min_containment: 0.01,
            max_metagenomes: usize::MAX,
            use_sketch: true,
            min_seed_hits: 2,
            diagonal_bin_bases: 64,
            flank_bases: 256,
            min_identity: 0.8,
            min_aligned_bases: 40,
            endpoint_bases: 256,
            circular: true,
            verify_resources: false,
            alignment: AlignmentConfig {
                max_cells: 1 << 24,
                ..AlignmentConfig::default()
            },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TraceResult {
    pub query_id: String,
    pub query_length: u64,
    pub index: TraceIndexIdentity,
    pub completion: SearchCompletion,
    pub candidates_screened: u32,
    pub metagenomes: Vec<MetagenomeTrace>,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TraceIndexIdentity {
    pub manifest_sha256: String,
    pub body_sha256: String,
    pub seed_k: u8,
    pub rescue_k15: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum SearchCompletion {
    Complete,
    CandidateBudgetExceeded { candidates_omitted: u32 },
}

#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct MetagenomeTrace {
    pub metagenome_id: MetagenomeId,
    pub name: String,
    pub shared_hashes: u32,
    pub containment: f64,
    pub exact_seed_hits: u64,
    pub compressed_bytes_read: u64,
    pub range_requests: u64,
    pub bgzf_blocks_decoded: u64,
    pub contigs: Vec<TraceContig>,
    pub mosaic: Mosaic,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct TraceContig {
    pub id: ContigId,
    pub name: String,
}

pub struct TraceEngine {
    jam_path: PathBuf,
    screen: QueryEngine,
    index: JidxReader,
    sample_to_metagenome: Vec<MetagenomeId>,
    s3: Option<S3Config>,
}

pub(crate) struct PreparedQuery {
    pub(crate) query_id: String,
    pub(crate) query_length: u64,
    query: Vec<u8>,
    positions_by_key: BTreeMap<u64, Vec<QuerySeed>>,
    lookup_identity: [u8; 32],
}

pub(crate) struct TraceCensus {
    pub(crate) candidates: Vec<Candidate>,
    pub(crate) frequencies: Vec<(u64, u32)>,
    pub(crate) lookups: Option<CachedSeedLookups>,
}

pub(crate) const LOOKUP_CACHE_BYTES: usize = 256 * 1024 * 1024;
static LOOKUP_CACHE_AVAILABLE: AtomicUsize = AtomicUsize::new(LOOKUP_CACHE_BYTES);

struct CacheReservation<'a> {
    available: &'a AtomicUsize,
    bytes: usize,
}

impl<'a> CacheReservation<'a> {
    fn acquire(available: &'a AtomicUsize, bytes: usize) -> Option<Self> {
        available
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |remaining| {
                remaining.checked_sub(bytes)
            })
            .ok()
            .map(|_| Self { available, bytes })
    }
}

impl Drop for CacheReservation<'_> {
    fn drop(&mut self) {
        self.available.fetch_add(self.bytes, Ordering::Relaxed);
    }
}

pub(crate) struct CachedSeedLookups {
    header_sha256: [u8; 32],
    query_identity: [u8; 32],
    file_identity: [u64; 7],
    through: Option<u64>,
    frozen: bool,
    groups: Vec<CachedDocumentGroup>,
    documents: Vec<SeedDocument>,
    _reservation: CacheReservation<'static>,
}

#[derive(Clone, Copy)]
struct CachedDocumentGroup {
    seed: SeedEntry,
    document_start: usize,
    document_count: usize,
}

#[derive(Clone, Copy)]
struct CachedSeedLookup<'a> {
    seed: SeedEntry,
    documents: &'a [SeedDocument],
}

impl CachedSeedLookups {
    fn new(
        header_sha256: [u8; 32],
        query_identity: [u8; 32],
        file_identity: [u64; 7],
        bytes: usize,
        query_keys: usize,
        document_count: u32,
    ) -> Option<Self> {
        let fixed = std::mem::size_of::<Option<Self>>().checked_add(4096)?;
        let documents_per_group = usize::try_from(document_count).ok()?;
        let bytes_per_group = std::mem::size_of::<CachedDocumentGroup>()
            .checked_add(documents_per_group.checked_mul(std::mem::size_of::<SeedDocument>())?)?;
        let group_limit = bytes
            .checked_sub(fixed)?
            .checked_div(bytes_per_group)?
            .min(query_keys);
        if group_limit == 0 {
            return None;
        }
        let document_limit = group_limit.checked_mul(documents_per_group)?;
        let reservation = CacheReservation::acquire(&LOOKUP_CACHE_AVAILABLE, bytes)?;
        let mut groups = Vec::new();
        let mut documents = Vec::new();
        groups.try_reserve_exact(group_limit).ok()?;
        documents.try_reserve_exact(document_limit).ok()?;
        let allocated = fixed
            .checked_add(
                groups
                    .capacity()
                    .checked_mul(std::mem::size_of::<CachedDocumentGroup>())?,
            )?
            .checked_add(
                documents
                    .capacity()
                    .checked_mul(std::mem::size_of::<SeedDocument>())?,
            )?;
        if allocated > reservation.bytes {
            return None;
        }
        Some(Self {
            header_sha256,
            query_identity,
            file_identity,
            through: None,
            frozen: false,
            groups,
            documents,
            _reservation: reservation,
        })
    }

    fn get(&self, key: u64) -> Option<Option<CachedSeedLookup<'_>>> {
        self.through.is_some_and(|through| key <= through).then(|| {
            self.groups
                .binary_search_by_key(&key, |group| group.seed.packed_key)
                .ok()
                .map(|index| {
                    let group = self.groups[index];
                    CachedSeedLookup {
                        seed: group.seed,
                        documents: &self.documents
                            [group.document_start..group.document_start + group.document_count],
                    }
                })
        })
    }

    fn cache_negative(&mut self, key: u64) {
        if !self.frozen {
            self.through = Some(key);
        }
    }

    fn cache_group(&mut self, key: u64, seed: SeedEntry, documents: &[SeedDocument]) {
        if self.frozen {
            return;
        }
        let Some(document_end) = self.documents.len().checked_add(documents.len()) else {
            self.frozen = true;
            return;
        };
        if self.groups.len() == self.groups.capacity() || document_end > self.documents.capacity() {
            self.frozen = true;
            return;
        }
        let document_start = self.documents.len();
        self.documents.extend_from_slice(documents);
        self.groups.push(CachedDocumentGroup {
            seed,
            document_start,
            document_count: documents.len(),
        });
        self.through = Some(key);
    }
}

impl TraceEngine {
    pub fn open(
        jam: impl AsRef<Path>,
        jidx: impl AsRef<Path>,
        manifest: impl AsRef<Path>,
        s3: Option<S3Config>,
    ) -> Result<Self, TraceError> {
        let jam = jam.as_ref();
        let index = JidxReader::open(jidx)?;
        let manifest_sha256 = sha256_reader(BufReader::new(File::open(manifest)?))?;
        if manifest_sha256 != index.header().manifest_sha256 {
            return Err(TraceError::Invalid(
                "JIDX belongs to a different root manifest",
            ));
        }
        let screen = QueryEngine::open(jam)?;
        let mut by_name = HashMap::new();
        for id in 0..index.header().document_count {
            let name = index
                .metagenome_name(id)?
                .ok_or(TraceError::Invalid("missing JIDX metagenome"))?;
            if by_name.insert(name.to_string(), id).is_some() {
                return Err(TraceError::Invalid("duplicate JIDX metagenome"));
            }
        }
        let mut sample_to_metagenome = Vec::with_capacity(screen.reader().sample_names().len());
        for name in screen.reader().sample_names() {
            sample_to_metagenome.push(
                *by_name
                    .get(name)
                    .ok_or(TraceError::Invalid("JAM and JIDX names differ"))?,
            );
        }
        if sample_to_metagenome.len() != by_name.len() {
            return Err(TraceError::Invalid("JAM and JIDX names differ"));
        }
        Ok(Self {
            jam_path: jam.to_path_buf(),
            screen,
            index,
            sample_to_metagenome,
            s3,
        })
    }

    pub fn verify_index(&self) -> Result<(), TraceError> {
        if sha256_reader(BufReader::new(File::open(&self.jam_path)?))?
            != self.index.header().jam_sha256
        {
            return Err(TraceError::Invalid(
                "JIDX belongs to a different JAM database",
            ));
        }
        self.index.verify_checksum()?;
        Ok(())
    }

    pub(crate) fn index(&self) -> &JidxReader {
        &self.index
    }

    pub fn search(
        &self,
        query_id: impl Into<String>,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<TraceResult, TraceError> {
        let prepared = prepare_query(
            query_id,
            sequence,
            config,
            self.index.header().k,
            self.index.header().rescue_k15,
        )?;
        let cache_bytes = LOOKUP_CACHE_BYTES / rayon::current_num_threads().max(1);
        let census = self.candidate_census(&prepared, config, cache_bytes)?;
        let completion = candidate_completion(census.candidates.len(), config.max_metagenomes)?;
        let mut candidates = census.candidates;
        candidates.truncate(config.max_metagenomes);
        let mut frequencies = census.frequencies;
        frequencies.sort_unstable_by_key(|&(key, frequency)| (frequency, key));
        let key_order = frequencies
            .into_iter()
            .map(|(key, _)| key)
            .collect::<Vec<_>>();
        let candidates_screened =
            u32::try_from(candidates.len()).map_err(|_| TraceError::Invalid("candidate count"))?;
        let metagenomes = self.trace_selected(
            &prepared,
            candidates,
            &key_order,
            census.lookups.as_ref(),
            config,
        )?;
        Ok(TraceResult {
            query_id: prepared.query_id,
            query_length: prepared.query_length,
            index: TraceIndexIdentity {
                manifest_sha256: digest_hex(self.index.header().manifest_sha256),
                body_sha256: digest_hex(self.index.header().body_sha256),
                seed_k: self.index.header().k,
                rescue_k15: self.index.header().rescue_k15,
            },
            completion,
            candidates_screened,
            metagenomes,
        })
    }

    pub(crate) fn candidate_census(
        &self,
        prepared: &PreparedQuery,
        config: TraceConfig,
        cache_bytes: usize,
    ) -> Result<TraceCensus, TraceError> {
        validate_config(config)?;
        let sketch_candidates =
            self.screen_candidates(&prepared.query_id, &prepared.query, config)?;
        let mut candidates = BTreeMap::new();
        for candidate in sketch_candidates {
            candidates.insert(candidate.id, candidate);
        }
        self.index
            .verify_query_filter_pages(prepared.positions_by_key.keys())?;
        let mut frequencies = Vec::new();
        let header_sha256 = sha256(&self.index.header().encode()?);
        let mut lookups =
            self.index
                .cache_file_identity()
                .ok()
                .flatten()
                .and_then(|file_identity| {
                    CachedSeedLookups::new(
                        header_sha256,
                        prepared.lookup_identity,
                        file_identity,
                        cache_bytes,
                        prepared.positions_by_key.len(),
                        self.index.header().document_count,
                    )
                });
        let mut entries = prepared.positions_by_key.iter();
        let mut chunk = Vec::new();
        chunk
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        let mut packed_keys = Vec::new();
        packed_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        loop {
            chunk.clear();
            packed_keys.clear();
            for _ in 0..SEED_LOOKUP_BATCH_KEYS {
                let Some((&packed_key, query_seeds)) = entries.next() else {
                    break;
                };
                chunk.push((packed_key, query_seeds));
                packed_keys.push(packed_key);
            }
            if chunk.is_empty() {
                break;
            }

            let index_seeds = self.index.find_seeds_batch(&packed_keys)?;
            self.index.advise_first_document_rows(&index_seeds);
            for ((packed_key, query_seeds), index_seed) in chunk.iter().copied().zip(index_seeds) {
                let Some(index_seed) = index_seed else {
                    if let Some(lookups) = &mut lookups {
                        lookups.cache_negative(packed_key);
                    }
                    continue;
                };
                frequencies.push((packed_key, index_seed.document_frequency));
                let query_positions = u64::try_from(query_seeds.len())
                    .map_err(|_| TraceError::Invalid("query seed count"))?;
                let documents = self.index.seed_documents(index_seed)?;
                for &document in &documents {
                    let hits = document
                        .occurrence_count
                        .checked_mul(query_positions)
                        .ok_or(TraceError::Invalid("exact seed hit count"))?;
                    if let Some(candidate) = candidates.get_mut(&document.metagenome_id) {
                        candidate.exact_seed_hits = candidate
                            .exact_seed_hits
                            .checked_add(hits)
                            .ok_or(TraceError::Invalid("exact seed hit count"))?;
                        continue;
                    }
                    let name = self
                        .index
                        .metagenome_name(document.metagenome_id)?
                        .ok_or(TraceError::Invalid("JIDX metagenome ID"))?
                        .to_string();
                    candidates.insert(
                        document.metagenome_id,
                        Candidate {
                            id: document.metagenome_id,
                            name,
                            shared_hashes: 0,
                            containment: 0.0,
                            exact_seed_hits: hits,
                        },
                    );
                }
                if let Some(lookups) = &mut lookups {
                    lookups.cache_group(packed_key, index_seed, &documents);
                }
            }
        }
        let mut candidates = candidates.into_values().collect::<Vec<_>>();
        candidates.sort_by(compare_candidates);
        Ok(TraceCensus {
            candidates,
            frequencies,
            lookups,
        })
    }

    pub(crate) fn trace_selected(
        &self,
        prepared: &PreparedQuery,
        candidates: Vec<Candidate>,
        key_order: &[u64],
        lookups: Option<&CachedSeedLookups>,
        config: TraceConfig,
    ) -> Result<Vec<MetagenomeTrace>, TraceError> {
        validate_config(config)?;
        self.index.enable_selected_front_metadata();
        let lookups = if let Some(lookups) = lookups {
            if lookups.header_sha256 != sha256(&self.index.header().encode()?)
                || lookups.query_identity != prepared.lookup_identity
            {
                return Err(TraceError::Invalid("cached seed lookup identity"));
            }
            match self.index.cache_file_identity() {
                Ok(Some(identity)) if identity == lookups.file_identity => Some(lookups),
                Ok(Some(_)) => return Err(TraceError::Invalid("cached seed lookup identity")),
                _ => None,
            }
        } else {
            None
        };
        let candidate_ids = candidates
            .iter()
            .map(|candidate| candidate.id)
            .collect::<HashSet<_>>();
        let mut region_hits = BTreeMap::<RegionKey, Vec<SeedHit>>::new();
        let mut packed_keys = Vec::new();
        packed_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        let mut filter_keys = Vec::new();
        filter_keys
            .try_reserve_exact(SEED_LOOKUP_BATCH_KEYS)
            .map_err(|_| TraceError::Invalid("query seed batch"))?;
        for key_chunk in key_order.chunks(SEED_LOOKUP_BATCH_KEYS) {
            packed_keys.clear();
            for &packed_key in key_chunk {
                if prepared.positions_by_key.contains_key(&packed_key) {
                    packed_keys.push(packed_key);
                }
            }
            if packed_keys.is_empty() {
                continue;
            }
            filter_keys.clear();
            filter_keys.extend(packed_keys.iter().copied().filter(|&key| {
                lookups.is_none_or(|lookups| lookups.through.is_none_or(|through| key > through))
            }));
            filter_keys.sort_unstable();
            filter_keys.dedup();
            self.index.verify_query_filter_pages(filter_keys.iter())?;
            let uncached_seeds = self.index.find_seeds_batch(&filter_keys)?;
            self.index.advise_first_document_rows(&uncached_seeds);
            let index_seeds = packed_keys
                .iter()
                .map(|&key| match lookups.and_then(|lookups| lookups.get(key)) {
                    Some(Some(cached)) => (Some(cached.seed), Some(cached.documents)),
                    Some(None) => (None, None),
                    None => (
                        uncached_seeds
                            [filter_keys.binary_search(&key).expect("uncached query key")],
                        None,
                    ),
                })
                .collect::<Vec<_>>();

            for (&packed_key, (index_seed, cached_documents)) in packed_keys.iter().zip(index_seeds)
            {
                let query_seeds = prepared
                    .positions_by_key
                    .get(&packed_key)
                    .expect("present query key");
                let Some(index_seed) = index_seed else {
                    continue;
                };
                let seed_k = seed_length(
                    self.index.header().k,
                    self.index.header().rescue_k15,
                    packed_key,
                )?;
                let decoded_documents;
                let documents = if let Some(documents) = cached_documents {
                    documents
                } else {
                    decoded_documents = self.index.seed_documents(index_seed)?;
                    &decoded_documents
                };
                for document in documents
                    .iter()
                    .copied()
                    .filter(|document| candidate_ids.contains(&document.metagenome_id))
                {
                    let occurrences = self.index.seed_document_occurrences(index_seed, document)?;
                    for seed in query_seeds {
                        for occurrence in &occurrences {
                            let contig = self
                                .index
                                .contig(occurrence.contig_id)?
                                .ok_or(TraceError::Invalid("missing occurrence contig"))?;
                            let strand =
                                if seed.canonical_orientation == occurrence.canonical_orientation {
                                    Strand::Forward
                                } else {
                                    Strand::Reverse
                                };
                            let oriented_position = match strand {
                                Strand::Forward => occurrence.position,
                                Strand::Reverse => contig
                                    .length
                                    .checked_sub(
                                        occurrence
                                            .position
                                            .checked_add(u64::from(seed_k))
                                            .ok_or(TraceError::Invalid("occurrence position"))?,
                                    )
                                    .ok_or(TraceError::Invalid("occurrence position"))?,
                            };
                            let diagonal =
                                i128::from(oriented_position) - i128::from(seed.position);
                            region_hits
                                .entry(RegionKey {
                                    metagenome_id: contig.metagenome_id,
                                    contig_id: contig.id,
                                    strand,
                                    k: seed_k,
                                })
                                .or_default()
                                .push(SeedHit {
                                    query: seed.position,
                                    target: oriented_position,
                                    diagonal,
                                });
                        }
                    }
                }
            }
        }
        let regions = form_regions(region_hits, config.diagonal_bin_bases);

        let tasks = self.tasks(regions, prepared.query_length, config)?;
        let (loaded, reads) = self.load_ranges(&tasks, config.verify_resources)?;
        let fragments = self.align_tasks(&prepared.query, &tasks, &loaded, config)?;
        let mut by_metagenome = BTreeMap::<MetagenomeId, Vec<Fragment>>::new();
        for (metagenome_id, fragment) in fragments {
            by_metagenome
                .entry(metagenome_id)
                .or_default()
                .push(fragment);
        }

        let mut metagenomes = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            let fragments = by_metagenome.remove(&candidate.id).unwrap_or_default();
            let mosaic = build_mosaic(prepared.query_length, &fragments)?;
            let mut contig_ids = BTreeSet::new();
            for fragment in mosaic
                .primary
                .iter()
                .map(|selected| &selected.fragment)
                .chain(mosaic.alternatives.iter())
            {
                contig_ids.insert(fragment.contig_id);
            }
            let mut contigs = Vec::with_capacity(contig_ids.len());
            for id in contig_ids {
                let contig = self
                    .index
                    .contig(id)?
                    .ok_or(TraceError::Invalid("missing result contig"))?;
                contigs.push(TraceContig {
                    id,
                    name: contig.name.to_string(),
                });
            }
            let (stats, bgzf_blocks_decoded) =
                reads.get(&candidate.id).copied().unwrap_or_default();
            metagenomes.push(MetagenomeTrace {
                metagenome_id: candidate.id,
                name: candidate.name,
                shared_hashes: candidate.shared_hashes,
                containment: candidate.containment,
                exact_seed_hits: candidate.exact_seed_hits,
                compressed_bytes_read: stats.bytes_read,
                range_requests: stats.read_requests,
                bgzf_blocks_decoded,
                contigs,
                mosaic,
            });
        }
        Ok(metagenomes)
    }

    fn screen_candidates(
        &self,
        query_id: &str,
        query: &[u8],
        config: TraceConfig,
    ) -> Result<Vec<Candidate>, TraceError> {
        if !config.use_sketch {
            return Ok(Vec::new());
        }
        let sketch = QuerySketch::from_sequence(query_id, query, self.screen.reader())?;
        let result = self
            .screen
            .query_sketch(&sketch)
            .into_iter()
            .next()
            .ok_or(TraceError::Invalid("missing JAM query result"))?;
        result
            .matches
            .into_iter()
            .filter(|candidate| candidate.containment >= config.min_containment)
            .map(|candidate| {
                let id = *self
                    .sample_to_metagenome
                    .get(candidate.sample_id as usize)
                    .ok_or(TraceError::Invalid("JAM sample ID"))?;
                let name = self
                    .index
                    .metagenome_name(id)?
                    .ok_or(TraceError::Invalid("JIDX metagenome ID"))?
                    .to_string();
                Ok(Candidate {
                    id,
                    name,
                    shared_hashes: candidate.hit_count,
                    containment: candidate.containment,
                    exact_seed_hits: 0,
                })
            })
            .collect::<Result<Vec<_>, TraceError>>()
    }

    fn tasks(
        &self,
        regions: Vec<(RegionKey, RegionAccumulator)>,
        query_length: u64,
        config: TraceConfig,
    ) -> Result<Vec<AlignmentTask>, TraceError> {
        let mut tasks = Vec::new();
        for (key, region) in regions {
            if region.hits < minimum_region_hits(key.k, config.min_seed_hits) {
                continue;
            }
            let k = u64::from(key.k);
            let contig = self
                .index
                .contig(key.contig_id)?
                .ok_or(TraceError::Invalid("missing region contig"))?;
            let query_start = region.query_start.saturating_sub(config.flank_bases);
            let query_end = region
                .query_end
                .saturating_add(k)
                .saturating_add(config.flank_bases);
            let query_span = if config.circular {
                query_end.saturating_sub(query_start).min(query_length)
            } else {
                query_end.min(query_length).saturating_sub(query_start)
            };
            if query_span == 0 {
                continue;
            }
            let oriented_start = region.target_start.saturating_sub(config.flank_bases);
            let oriented_end = region
                .target_end
                .saturating_add(k)
                .saturating_add(config.flank_bases)
                .min(contig.length);
            let (target_start, target_end) = match key.strand {
                Strand::Forward => (oriented_start, oriented_end),
                Strand::Reverse => (contig.length - oriented_end, contig.length - oriented_start),
            };
            let target_relative = region.target_start - oriented_start;
            let query_relative = region.query_start - query_start;
            let diagonal_offset =
                i64::try_from(i128::from(target_relative) - i128::from(query_relative))
                    .map_err(|_| TraceError::Invalid("task diagonal"))?;
            tasks.push(AlignmentTask {
                metagenome_id: key.metagenome_id,
                contig_id: key.contig_id,
                strand: key.strand,
                query_start,
                query_span,
                target_start,
                target_end,
                diagonal_offset,
            });
        }
        tasks.sort_unstable_by_key(|task| {
            (
                task.metagenome_id,
                task.contig_id,
                task.strand,
                task.target_start,
                task.query_start,
            )
        });
        Ok(tasks)
    }

    fn load_ranges(
        &self,
        tasks: &[AlignmentTask],
        verify: bool,
    ) -> Result<LoadedRanges, TraceError> {
        let mut spans = BTreeMap::<(MetagenomeId, ContigId), Vec<(u64, u64)>>::new();
        for task in tasks {
            spans
                .entry((task.metagenome_id, task.contig_id))
                .or_default()
                .push((task.target_start, task.target_end));
        }
        for contig_spans in spans.values_mut() {
            coalesce_spans(contig_spans);
        }
        let mut loaded = BTreeMap::new();
        let mut reads = HashMap::new();
        for metagenome_id in tasks
            .iter()
            .map(|task| task.metagenome_id)
            .collect::<BTreeSet<_>>()
        {
            let source = self
                .index
                .metagenome(metagenome_id)?
                .ok_or(TraceError::Invalid("missing source metagenome"))?;
            let mut reader = BgzfReader::open(source, self.s3.as_ref(), verify)?;
            for (&(_, contig_id), contig_spans) in
                spans.range((metagenome_id, 0)..=(metagenome_id, u32::MAX))
            {
                let contig = self
                    .index
                    .contig(contig_id)?
                    .ok_or(TraceError::Invalid("missing source contig"))?;
                let loaded_ranges = loaded
                    .entry((metagenome_id, contig_id))
                    .or_insert_with(Vec::new);
                for &(start, end) in contig_spans {
                    loaded_ranges.push(LoadedRange {
                        offset: start,
                        end,
                        sequence: reader.read_contig_range(contig, start, end)?,
                    });
                }
            }
            reads.insert(
                metagenome_id,
                (reader.range_stats(), reader.blocks_decoded()),
            );
        }
        Ok((loaded, reads))
    }

    fn align_tasks(
        &self,
        query: &[u8],
        tasks: &[AlignmentTask],
        loaded: &BTreeMap<(MetagenomeId, ContigId), Vec<LoadedRange>>,
        config: TraceConfig,
    ) -> Result<Vec<(MetagenomeId, Fragment)>, TraceError> {
        let (max_query_bases, max_target_bases) =
            tasks
                .iter()
                .try_fold((0usize, 0usize), |(query_bases, target_bases), task| {
                    let query = usize::try_from(task.query_span)
                        .map_err(|_| TraceError::Invalid("query window"))?;
                    let target = task
                        .target_end
                        .checked_sub(task.target_start)
                        .and_then(|span| usize::try_from(span).ok())
                        .ok_or(TraceError::Invalid("loaded range"))?;
                    Ok::<_, TraceError>((query_bases.max(query), target_bases.max(target)))
                })?;
        let workspace = || {
            TraceAlignmentWorkspace::acquire(
                max_query_bases,
                max_target_bases,
                config.endpoint_bases,
                config.alignment,
            )
        };
        tasks
            .par_iter()
            .map_init(workspace, |workspace, task| {
                let workspace = workspace
                    .as_mut()
                    .map_err(|error| TraceError::AlignmentAdmission(error.to_string()))?
                    .workspace_mut();
                let query_window =
                    linearize_query(query, task.query_start, task.query_span, config.circular)?;
                let loaded = loaded
                    .get(&(task.metagenome_id, task.contig_id))
                    .and_then(|ranges| {
                        ranges.iter().find(|range| {
                            range.offset <= task.target_start && range.end >= task.target_end
                        })
                    })
                    .ok_or(TraceError::Invalid("missing loaded range"))?;
                let start = usize::try_from(task.target_start - loaded.offset)
                    .map_err(|_| TraceError::Invalid("loaded range"))?;
                let end = usize::try_from(task.target_end - loaded.offset)
                    .map_err(|_| TraceError::Invalid("loaded range"))?;
                let target = loaded
                    .sequence
                    .get(start..end)
                    .ok_or(TraceError::Invalid("loaded range"))?;
                let mut alignment_config = config.alignment;
                alignment_config.diagonal_offset = task.diagonal_offset;
                let core = match workspace.align_oriented(
                    &query_window,
                    target,
                    task.target_start,
                    task.strand,
                    alignment_config,
                ) {
                    Ok(alignment) => alignment,
                    Err(AlignmentError::NoAlignment) => return Ok(None),
                    Err(error) => return Err(error.into()),
                };
                let completed = workspace
                    .complete_endpoints(
                        core.clone(),
                        &query_window,
                        target,
                        task.target_start,
                        config.endpoint_bases,
                        alignment_config,
                    )?
                    .alignment;
                let alignment = if completed.identity() >= config.min_identity {
                    completed
                } else {
                    core
                };
                if alignment.identity() < config.min_identity
                    || alignment.query_interval.len() < config.min_aligned_bases
                {
                    return Ok(None);
                }
                let query_segments = query_segments(
                    task.query_start,
                    alignment.query_interval,
                    u64::try_from(query.len()).map_err(|_| TraceError::Invalid("query length"))?,
                    config.circular,
                )?;
                Ok(Some((
                    task.metagenome_id,
                    Fragment {
                        contig_id: task.contig_id,
                        query_segments,
                        alignment,
                    },
                )))
            })
            .filter_map(|result| match result {
                Ok(Some(value)) => Some(Ok(value)),
                Ok(None) => None,
                Err(error) => Some(Err(error)),
            })
            .collect()
    }
}

type LoadedRanges = (
    BTreeMap<(MetagenomeId, ContigId), Vec<LoadedRange>>,
    HashMap<MetagenomeId, (crate::range_source::RangeStats, u64)>,
);

struct LoadedRange {
    offset: u64,
    end: u64,
    sequence: Vec<u8>,
}

pub(crate) struct Candidate {
    pub(crate) id: MetagenomeId,
    pub(crate) name: String,
    pub(crate) shared_hashes: u32,
    pub(crate) containment: f64,
    pub(crate) exact_seed_hits: u64,
}

pub(crate) fn compare_candidates(left: &Candidate, right: &Candidate) -> std::cmp::Ordering {
    right
        .exact_seed_hits
        .cmp(&left.exact_seed_hits)
        .then_with(|| right.containment.total_cmp(&left.containment))
        .then_with(|| right.shared_hashes.cmp(&left.shared_hashes))
        .then_with(|| left.name.cmp(&right.name))
}

#[derive(Clone, Copy)]
struct QuerySeed {
    packed_key: u64,
    position: u64,
    canonical_orientation: bool,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct RegionKey {
    metagenome_id: MetagenomeId,
    contig_id: ContigId,
    strand: Strand,
    k: u8,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct SeedHit {
    query: u64,
    target: u64,
    diagonal: i128,
}

struct RegionAccumulator {
    query_start: u64,
    query_end: u64,
    target_start: u64,
    target_end: u64,
    diagonal_min: i128,
    diagonal_max: i128,
    hits: u32,
}

impl RegionAccumulator {
    fn new(hit: SeedHit) -> Self {
        Self {
            query_start: hit.query,
            query_end: hit.query,
            target_start: hit.target,
            target_end: hit.target,
            diagonal_min: hit.diagonal,
            diagonal_max: hit.diagonal,
            hits: 1,
        }
    }

    fn accepts(&self, hit: SeedHit, max_diagonal_drift: u64) -> bool {
        hit.query > self.query_end
            && hit.target > self.target_end
            && self.diagonal_min.min(hit.diagonal) + i128::from(max_diagonal_drift)
                >= self.diagonal_max.max(hit.diagonal)
    }

    fn add(&mut self, hit: SeedHit) {
        self.query_end = hit.query;
        self.target_end = hit.target;
        self.diagonal_min = self.diagonal_min.min(hit.diagonal);
        self.diagonal_max = self.diagonal_max.max(hit.diagonal);
        self.hits = self.hits.saturating_add(1);
    }
}

fn form_regions(
    hits_by_contig: BTreeMap<RegionKey, Vec<SeedHit>>,
    max_diagonal_drift: u64,
) -> Vec<(RegionKey, RegionAccumulator)> {
    let mut output = Vec::new();
    for (key, mut hits) in hits_by_contig {
        hits.sort_unstable_by_key(|hit| (hit.query, hit.target));
        let mut regions = Vec::<RegionAccumulator>::new();
        for hit in hits {
            if let Some(region) = regions
                .iter_mut()
                .rev()
                .find(|region| region.accepts(hit, max_diagonal_drift))
            {
                region.add(hit);
            } else {
                regions.push(RegionAccumulator::new(hit));
            }
        }
        output.extend(regions.into_iter().map(|region| (key, region)));
    }
    output
}

struct AlignmentTask {
    metagenome_id: MetagenomeId,
    contig_id: ContigId,
    strand: Strand,
    query_start: u64,
    query_span: u64,
    target_start: u64,
    target_end: u64,
    diagonal_offset: i64,
}

pub(crate) fn prepare_query(
    query_id: impl Into<String>,
    sequence: &[u8],
    config: TraceConfig,
    k: u8,
    rescue_k15: bool,
) -> Result<PreparedQuery, TraceError> {
    validate_config(config)?;
    let query_id = query_id.into();
    if query_id.is_empty()
        || query_id
            .bytes()
            .any(|byte| matches!(byte, 0 | b'\n' | b'\r'))
    {
        return Err(TraceError::Invalid("query ID"));
    }
    let query = sequence.normalize(false).into_owned();
    if query.is_empty() {
        return Err(TraceError::Invalid("empty query"));
    }
    let query_length =
        u64::try_from(query.len()).map_err(|_| TraceError::Invalid("query length"))?;
    let mut positions_by_key = BTreeMap::<u64, Vec<QuerySeed>>::new();
    for seed in query_seeds(&query, k, rescue_k15, config.circular)? {
        positions_by_key
            .entry(seed.packed_key)
            .or_default()
            .push(seed);
    }
    let mut lookup_identity = [0; 35];
    lookup_identity[..32].copy_from_slice(&sha256(&query));
    lookup_identity[32..].copy_from_slice(&[k, u8::from(rescue_k15), u8::from(config.circular)]);
    Ok(PreparedQuery {
        query_id,
        query_length,
        query,
        positions_by_key,
        lookup_identity: sha256(&lookup_identity),
    })
}

fn query_seeds(
    query: &[u8],
    k: u8,
    rescue_k15: bool,
    circular: bool,
) -> Result<Vec<QuerySeed>, TraceError> {
    let mut seeds = extract_query_seeds(query, k, circular)?;
    if rescue_k15 {
        let mut rescue = extract_query_seeds(query, 15, circular)?;
        for seed in &mut rescue {
            seed.packed_key |= RESCUE_K15_TAG;
        }
        seeds.extend(rescue);
    }
    Ok(seeds)
}

fn extract_query_seeds(query: &[u8], k: u8, circular: bool) -> Result<Vec<QuerySeed>, TraceError> {
    if query.len() < usize::from(k) {
        return Ok(Vec::new());
    }
    let mut sequence = query.to_vec();
    if circular {
        sequence.extend_from_slice(&query[..usize::from(k) - 1]);
    }
    let mut seeds = Vec::new();
    for (position, kmer, orientation) in sequence.bit_kmers(k, true) {
        if position >= query.len() {
            break;
        }
        let seed = QuerySeed {
            packed_key: kmer.0,
            position: u64::try_from(position)
                .map_err(|_| TraceError::Invalid("query seed position"))?,
            canonical_orientation: orientation,
        };
        seeds.push(seed);
    }
    Ok(seeds)
}

fn minimum_region_hits(k: u8, configured: u32) -> u32 {
    if k == 15 {
        configured.max(3)
    } else {
        configured
    }
}

pub(crate) fn digest_hex(digest: [u8; 32]) -> String {
    let mut output = String::with_capacity(64);
    for byte in digest {
        write!(&mut output, "{byte:02x}").expect("writing to a string cannot fail");
    }
    output
}

pub(crate) fn candidate_completion(
    candidate_count: usize,
    max_metagenomes: usize,
) -> Result<SearchCompletion, TraceError> {
    let omitted = candidate_count.saturating_sub(max_metagenomes);
    if omitted == 0 {
        Ok(SearchCompletion::Complete)
    } else {
        Ok(SearchCompletion::CandidateBudgetExceeded {
            candidates_omitted: u32::try_from(omitted)
                .map_err(|_| TraceError::Invalid("candidate count"))?,
        })
    }
}

fn coalesce_spans(spans: &mut Vec<(u64, u64)>) {
    spans.sort_unstable();
    let mut merged = Vec::with_capacity(spans.len());
    for &(start, end) in spans.iter() {
        if let Some((_, previous_end)) = merged.last_mut()
            && start <= *previous_end
        {
            *previous_end = (*previous_end).max(end);
        } else {
            merged.push((start, end));
        }
    }
    *spans = merged;
}

fn linearize_query(
    query: &[u8],
    start: u64,
    span: u64,
    circular: bool,
) -> Result<Vec<u8>, TraceError> {
    let start = usize::try_from(start).map_err(|_| TraceError::Invalid("query window"))?;
    let span = usize::try_from(span).map_err(|_| TraceError::Invalid("query window"))?;
    if !circular {
        return query
            .get(start..start + span)
            .map(<[u8]>::to_vec)
            .ok_or(TraceError::Invalid("query window"));
    }
    Ok((0..span)
        .map(|offset| query[(start + offset) % query.len()])
        .collect())
}

fn query_segments(
    window_start: u64,
    local: Interval,
    query_length: u64,
    circular: bool,
) -> Result<Vec<Interval>, TraceError> {
    let start = window_start
        .checked_add(local.start)
        .ok_or(TraceError::Invalid("query coordinates"))?;
    let end = window_start
        .checked_add(local.end)
        .ok_or(TraceError::Invalid("query coordinates"))?;
    if !circular || end <= query_length {
        return Ok(vec![Interval::new(start, end)?]);
    }
    if start >= query_length {
        return Ok(vec![Interval::new(
            start - query_length,
            end - query_length,
        )?]);
    }
    let wrapped = end % query_length;
    let mut segments = vec![Interval::new(start, query_length)?];
    if wrapped != 0 {
        segments.push(Interval::new(0, wrapped)?);
    }
    Ok(segments)
}

fn validate_config(config: TraceConfig) -> Result<(), TraceError> {
    if !config.min_containment.is_finite()
        || !(0.0..=1.0).contains(&config.min_containment)
        || config.max_metagenomes == 0
        || config.min_seed_hits == 0
        || config.diagonal_bin_bases == 0
        || !config.min_identity.is_finite()
        || !(0.0..=1.0).contains(&config.min_identity)
        || config.min_aligned_bases == 0
    {
        return Err(TraceError::Invalid("trace configuration"));
    }
    Ok(())
}

#[derive(Debug, Error)]
pub enum TraceError {
    #[error("trace I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Database(#[from] ReaderError),
    #[error(transparent)]
    Query(#[from] QueryError),
    #[error(transparent)]
    Jidx(#[from] JidxReaderError),
    #[error(transparent)]
    JidxFormat(#[from] JidxError),
    #[error(transparent)]
    Bgzf(#[from] BgzfError),
    #[error(transparent)]
    Alignment(#[from] AlignmentError),
    #[error("alignment workspace admission failed: {0}")]
    AlignmentAdmission(String),
    #[error(transparent)]
    Mosaic(#[from] MosaicError),
    #[error("invalid trace input: {0}")]
    Invalid(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::AlignmentWorkspace;
    use crate::cli::handlers::{TraceArgs, TraceInput, handle_trace_command};
    use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
    use crate::writer::{BuildConfig, build};
    use noodles_bgzf::{self as bgzf, gzi};
    use std::io::Write;

    #[test]
    fn lookup_cache_reservations_share_and_release_the_byte_limit() {
        let available = AtomicUsize::new(10);
        let first = CacheReservation::acquire(&available, 6).unwrap();
        assert!(CacheReservation::acquire(&available, 5).is_none());
        assert_eq!(available.load(Ordering::Relaxed), 4);
        let second = CacheReservation::acquire(&available, 4).unwrap();
        assert_eq!(available.load(Ordering::Relaxed), 0);
        drop(first);
        assert_eq!(available.load(Ordering::Relaxed), 6);
        drop(second);
        assert_eq!(available.load(Ordering::Relaxed), 10);
    }

    fn sequence() -> String {
        let mut state = 7u64;
        (0..128)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect()
    }

    #[test]
    fn trace_defaults_align_long_windows_with_explicit_caps_preserved() {
        let sequence = vec![b'A'; 64_000];
        let mut workspace = AlignmentWorkspace::default();
        let explicit = TraceConfig {
            alignment: AlignmentConfig::default(),
            ..TraceConfig::default()
        };
        assert!(matches!(
            workspace.align(&sequence, &sequence, explicit.alignment),
            Err(AlignmentError::MatrixTooLarge {
                cells,
                max_cells: 4_000_000,
            }) if cells > 4_000_000
        ));

        let config = TraceConfig::default();
        for (target, strand) in [
            (sequence.clone(), Strand::Forward),
            (vec![b'T'; sequence.len()], Strand::Reverse),
        ] {
            let alignment = workspace
                .align_oriented(&sequence, &target, 100, strand, config.alignment)
                .unwrap();
            alignment.validate_cigar().unwrap();
            assert_eq!(alignment.query_interval, Interval::new(0, 64_000).unwrap());
            assert_eq!(
                alignment.target_interval,
                Interval::new(100, 64_100).unwrap()
            );
            assert_eq!(alignment.strand, strand);
            assert_eq!(alignment.cigar, "64000=");
            assert_eq!(alignment.matches, 64_000);
            assert_eq!(alignment.identity(), 1.0);
        }
    }

    #[test]
    fn candidate_census_matches_scalar_across_batch_boundary() {
        let directory = tempfile::tempdir().unwrap();
        let mut state = 97u64;
        let sequence = (0..50_000)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect::<String>();
        let fasta = directory.path().join("sample.fa");
        std::fs::write(&fasta, format!(">sample\n{sequence}\n")).unwrap();
        let jam = directory.path().join("database.jam");
        build(
            &[fasta],
            &jam,
            &BuildConfig {
                kmer_size: 21,
                fscale: 1,
                singleton: true,
                memory: 1,
                ..BuildConfig::default()
            },
        )
        .unwrap();

        let bgzf_path = directory.path().join("sample.bgz");
        let fai_path = directory.path().join("sample.bgz.fai");
        let gzi_path = directory.path().join("sample.bgz.gzi");
        let mut raw = b">contig\n".to_vec();
        for line in sequence.as_bytes().chunks(80) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(&raw).unwrap();
        writer.finish().unwrap();
        std::fs::write(
            &fai_path,
            format!("contig\t{}\t8\t80\t81\n", sequence.len()),
        )
        .unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display()
            ),
        )
        .unwrap();
        let jidx = directory.path().join("database.jidx");
        build_local_jidx(
            &jam,
            &manifest,
            &jidx,
            JidxBuildConfig {
                k: 21,
                minimizer_window: 1,
                rescue_k15: false,
            },
        )
        .unwrap();

        let engine = TraceEngine::open(&jam, &jidx, &manifest, None).unwrap();
        let config = TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        };
        let prepared = prepare_query("batch", sequence.as_bytes(), config, 21, false).unwrap();
        assert!(prepared.positions_by_key.len() > SEED_LOOKUP_BATCH_KEYS);
        let census = engine.candidate_census(&prepared, config, 0).unwrap();

        let mut scalar_frequencies = Vec::new();
        let mut scalar_hits = BTreeMap::<MetagenomeId, u64>::new();
        for (&packed_key, query_seeds) in &prepared.positions_by_key {
            let Some(seed) = engine.index.find_seed(packed_key).unwrap() else {
                continue;
            };
            scalar_frequencies.push((packed_key, seed.document_frequency));
            let query_positions = u64::try_from(query_seeds.len()).unwrap();
            for document in engine.index.seed_documents(seed).unwrap() {
                *scalar_hits.entry(document.metagenome_id).or_default() +=
                    document.occurrence_count * query_positions;
            }
        }

        assert_eq!(census.frequencies, scalar_frequencies);
        assert_eq!(census.candidates.len(), scalar_hits.len());
        assert!(census.candidates.iter().all(|candidate| {
            scalar_hits.get(&candidate.id).copied() == Some(candidate.exact_seed_hits)
        }));

        let packed_order = census
            .frequencies
            .iter()
            .map(|&(key, _)| key)
            .collect::<Vec<_>>();
        assert!(packed_order.len() > SEED_LOOKUP_BATCH_KEYS);
        let mut rarity_order = packed_order.clone();
        rarity_order.reverse();
        assert_ne!(rarity_order, packed_order);
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &rarity_order, None, config)
                .unwrap()
                .is_empty()
        );
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &packed_order, None, config)
                .unwrap()
                .is_empty()
        );

        let missing = (0..1u64 << 42)
            .find(|key| !prepared.positions_by_key.contains_key(key))
            .unwrap();
        let mut repeated_and_missing = rarity_order.clone();
        repeated_and_missing.insert(SEED_LOOKUP_BATCH_KEYS - 1, missing);
        repeated_and_missing.insert(SEED_LOOKUP_BATCH_KEYS, rarity_order[SEED_LOOKUP_BATCH_KEYS]);
        assert!(
            engine
                .trace_selected(&prepared, Vec::new(), &repeated_and_missing, None, config)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn traces_an_arbitrary_query_deterministically() {
        let directory = tempfile::tempdir().unwrap();
        let sequence = sequence();
        let fasta = directory.path().join("sample.fa");
        std::fs::write(&fasta, format!(">sample\n{sequence}\n")).unwrap();
        let fasta2 = directory.path().join("sample2.fa");
        std::fs::write(&fasta2, format!(">sample2\n{sequence}\n")).unwrap();
        let jam = directory.path().join("database.jam");
        build(
            &[fasta, fasta2],
            &jam,
            &BuildConfig {
                kmer_size: 5,
                fscale: 1,
                singleton: true,
                memory: 1,
                ..BuildConfig::default()
            },
        )
        .unwrap();

        let bgzf_path = directory.path().join("sample.bgz");
        let fai_path = directory.path().join("sample.bgz.fai");
        let gzi_path = directory.path().join("sample.bgz.gzi");
        let mut raw = b">contig\n".to_vec();
        for line in sequence.as_bytes().chunks(16) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(&raw).unwrap();
        writer.finish().unwrap();
        std::fs::write(&fai_path, b"contig\t128\t8\t16\t17\n").unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let bgzf_path2 = directory.path().join("sample2.bgz");
        let fai_path2 = directory.path().join("sample2.bgz.fai");
        let gzi_path2 = directory.path().join("sample2.bgz.gzi");
        std::fs::copy(&bgzf_path, &bgzf_path2).unwrap();
        std::fs::copy(&fai_path, &fai_path2).unwrap();
        std::fs::copy(&gzi_path, &gzi_path2).unwrap();
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}},{{\"name\":\"sample2\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display(),
                bgzf_path2.display(),
                fai_path2.display(),
                gzi_path2.display()
            ),
        )
        .unwrap();
        let jidx = directory.path().join("database.jidx");
        build_local_jidx(
            &jam,
            &manifest,
            &jidx,
            JidxBuildConfig {
                k: 5,
                minimizer_window: 4,
                rescue_k15: false,
            },
        )
        .unwrap();

        let engine = TraceEngine::open(&jam, &jidx, &manifest, None).unwrap();
        let wrong_manifest = directory.path().join("wrong-manifest.json");
        std::fs::write(&wrong_manifest, b"{}").unwrap();
        assert!(TraceEngine::open(&jam, &jidx, wrong_manifest, None).is_err());
        engine.verify_index().unwrap();
        let config = TraceConfig {
            min_seed_hits: 2,
            diagonal_bin_bases: 32,
            flank_bases: 32,
            min_aligned_bases: 32,
            endpoint_bases: 32,
            alignment: AlignmentConfig {
                band_width: 32,
                ..AlignmentConfig::default()
            },
            ..TraceConfig::default()
        };
        let direct_config = TraceConfig {
            use_sketch: false,
            ..config
        };
        let prepared =
            prepare_query("census", sequence.as_bytes(), direct_config, 5, false).unwrap();
        let census = engine
            .candidate_census(&prepared, direct_config, 0)
            .unwrap();
        let mut explicit_hits = HashMap::<MetagenomeId, u64>::new();
        for (&packed_key, query_seeds) in &prepared.positions_by_key {
            let Some(seed) = engine.index.find_seed(packed_key).unwrap() else {
                continue;
            };
            for document in engine.index.seed_documents(seed).unwrap() {
                let occurrences = engine
                    .index
                    .seed_document_occurrences(seed, document)
                    .unwrap();
                *explicit_hits.entry(document.metagenome_id).or_default() +=
                    u64::try_from(occurrences.len() * query_seeds.len()).unwrap();
            }
        }
        assert!(census.candidates.iter().all(|candidate| {
            explicit_hits.get(&candidate.id).copied() == Some(candidate.exact_seed_hits)
        }));
        let mut rare_first = census.frequencies;
        rare_first.sort_unstable_by_key(|&(key, frequency)| (frequency, key));
        let rare_first = rare_first
            .into_iter()
            .map(|(key, _)| key)
            .collect::<Vec<_>>();
        let rare_first_result = engine
            .trace_selected(
                &prepared,
                census.candidates,
                &rare_first,
                None,
                direct_config,
            )
            .unwrap();
        let reverse_census = engine
            .candidate_census(&prepared, direct_config, 0)
            .unwrap();
        let reverse_order = rare_first.iter().rev().copied().collect::<Vec<_>>();
        let reverse_result = engine
            .trace_selected(
                &prepared,
                reverse_census.candidates,
                &reverse_order,
                None,
                direct_config,
            )
            .unwrap();
        assert_eq!(rare_first_result, reverse_result);
        let cache_fixed = std::mem::size_of::<Option<CachedSeedLookups>>() + 4096;
        let cache_per_key = std::mem::size_of::<CachedDocumentGroup>()
            + usize::try_from(engine.index.header().document_count).unwrap()
                * std::mem::size_of::<SeedDocument>();
        let one_group_cache_bytes = cache_fixed + cache_per_key;
        let full_cache_bytes = cache_fixed + cache_per_key * prepared.positions_by_key.len();
        for cache_bytes in [0, cache_fixed, one_group_cache_bytes, full_cache_bytes] {
            let mut cached = engine
                .candidate_census(&prepared, direct_config, cache_bytes)
                .unwrap();
            if let Some(lookups) = &cached.lookups {
                let allocated = cache_fixed
                    + lookups.groups.capacity() * std::mem::size_of::<CachedDocumentGroup>()
                    + lookups.documents.capacity() * std::mem::size_of::<SeedDocument>();
                assert!(allocated <= lookups._reservation.bytes);
                assert_eq!(
                    lookups.documents.len(),
                    lookups
                        .groups
                        .iter()
                        .map(|group| group.document_count)
                        .sum::<usize>()
                );
                if cache_bytes == one_group_cache_bytes {
                    assert_eq!(lookups.groups.len(), 1);
                    assert!(lookups.frozen);
                }
                if cache_bytes == full_cache_bytes {
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .all(|&key| lookups.get(key).is_some())
                    );
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .any(|&key| matches!(lookups.get(key), Some(None)))
                    );
                } else {
                    assert!(
                        prepared
                            .positions_by_key
                            .keys()
                            .any(|&key| lookups.get(key).is_none())
                    );
                }
                for &key in prepared.positions_by_key.keys() {
                    match lookups.get(key) {
                        Some(Some(cached)) => {
                            assert_eq!(Some(cached.seed), engine.index.find_seed(key).unwrap());
                            assert_eq!(
                                cached.documents,
                                engine.index.seed_documents(cached.seed).unwrap()
                            );
                        }
                        Some(None) => assert_eq!(None, engine.index.find_seed(key).unwrap()),
                        None => {}
                    }
                }
            } else {
                assert!(cache_bytes <= cache_fixed || cfg!(not(unix)));
            }
            let result = engine
                .trace_selected(
                    &prepared,
                    cached.candidates,
                    &reverse_order,
                    cached.lookups.as_ref(),
                    direct_config,
                )
                .unwrap();
            assert_eq!(rare_first_result, result);
            if let Some(lookups) = &mut cached.lookups {
                lookups.header_sha256[0] ^= 1;
                assert!(
                    engine
                        .trace_selected(
                            &prepared,
                            Vec::new(),
                            &reverse_order,
                            Some(lookups),
                            direct_config
                        )
                        .is_err()
                );
                lookups.header_sha256[0] ^= 1;
                lookups.query_identity[0] ^= 1;
                assert!(
                    engine
                        .trace_selected(
                            &prepared,
                            Vec::new(),
                            &reverse_order,
                            Some(lookups),
                            direct_config
                        )
                        .is_err()
                );
            }
        }
        let first = engine
            .search("plasmid", sequence.as_bytes(), config)
            .unwrap();
        let padded_query = [
            b"N".repeat(64),
            sequence.as_bytes().to_vec(),
            b"N".repeat(64),
        ]
        .concat();
        let padded_target = [
            b"A".repeat(64),
            sequence.as_bytes().to_vec(),
            b"T".repeat(64),
        ]
        .concat();
        let fragments = engine
            .align_tasks(
                &padded_query,
                &[AlignmentTask {
                    metagenome_id: 0,
                    contig_id: 0,
                    strand: Strand::Forward,
                    query_start: 0,
                    query_span: padded_query.len() as u64,
                    target_start: 0,
                    target_end: padded_target.len() as u64,
                    diagonal_offset: 0,
                }],
                &BTreeMap::from([(
                    (0, 0),
                    vec![LoadedRange {
                        offset: 0,
                        end: padded_target.len() as u64,
                        sequence: padded_target,
                    }],
                )]),
                TraceConfig {
                    endpoint_bases: 64,
                    circular: false,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(fragments.len(), 1);
        assert_eq!(fragments[0].1.alignment.identity(), 1.0);
        assert_eq!(
            fragments[0].1.alignment.query_interval,
            Interval::new(64, 192).unwrap()
        );
        let second = engine
            .search("plasmid", sequence.as_bytes(), config)
            .unwrap();
        assert_eq!(first, second);
        assert_eq!(first.completion, SearchCompletion::Complete);
        assert_eq!(first.index.seed_k, 5);
        assert!(!first.index.rescue_k15);
        assert_eq!(first.index.manifest_sha256.len(), 64);
        assert_eq!(first.index.body_sha256.len(), 64);
        assert_eq!(first.metagenomes.len(), 2);
        assert_eq!(first.metagenomes[0].mosaic.covered_bases, 128);
        assert_eq!(first.metagenomes[0].contigs[0].name, "contig");
        assert!(first.metagenomes[0].bgzf_blocks_decoded > 0);
        serde_json::to_vec(&first).unwrap();

        let query = directory.path().join("query.fa");
        std::fs::write(&query, format!(">sample description\n{sequence}\n")).unwrap();
        let output = directory.path().join("trace.jsonl");
        handle_trace_command(TraceArgs {
            query: query.clone(),
            input: TraceInput::Shard {
                database: jam.clone(),
                index: jidx.clone(),
                manifest: manifest.clone(),
            },
            audit_index: false,
            output: output.clone(),
            query_id: None,
            config,
            s3: None,
            force: false,
        })
        .unwrap();
        let published = std::fs::read_to_string(output).unwrap();
        assert!(published.ends_with('\n'));
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&published).unwrap()["query_id"],
            "sample"
        );

        std::fs::write(&query, format!(">second\n{sequence}\n>first\n{sequence}\n")).unwrap();
        let batch_output = directory.path().join("batch.jsonl");
        handle_trace_command(TraceArgs {
            query,
            input: TraceInput::Shard {
                database: jam,
                index: jidx.clone(),
                manifest,
            },
            audit_index: false,
            output: batch_output.clone(),
            query_id: None,
            config,
            s3: None,
            force: false,
        })
        .unwrap();
        let batch = std::fs::read_to_string(batch_output)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
            .collect::<Vec<_>>();
        let singles = ["second", "first"].map(|id| {
            serde_json::to_value(engine.search(id, sequence.as_bytes(), config).unwrap()).unwrap()
        });
        assert_eq!(batch, singles);

        let direct = engine
            .search(
                "direct",
                sequence.as_bytes(),
                TraceConfig {
                    use_sketch: false,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(direct.completion, SearchCompletion::Complete);
        assert_eq!(direct.metagenomes.len(), 2);
        assert_eq!(direct.metagenomes[0].shared_hashes, 0);
        assert_eq!(direct.metagenomes[0].mosaic.covered_bases, 128);

        let absent = (0..1 << 10)
            .map(|packed| {
                (0..5)
                    .rev()
                    .map(|shift| b"ACGT"[(packed >> (shift * 2)) & 3] as char)
                    .collect::<String>()
            })
            .find(|word| {
                let reverse_complement = word
                    .bytes()
                    .rev()
                    .map(|base| match base {
                        b'A' => 'T',
                        b'C' => 'G',
                        b'G' => 'C',
                        b'T' => 'A',
                        _ => unreachable!(),
                    })
                    .collect::<String>();
                !sequence.contains(word) && !sequence.contains(&reverse_complement)
            })
            .unwrap();
        let partial = format!("{}{}", &sequence[..64], absent.repeat(13));
        let sketch_miss = engine
            .search(
                "sketch-miss",
                partial.as_bytes(),
                TraceConfig {
                    min_containment: 1.0,
                    circular: false,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(sketch_miss.metagenomes.len(), 2);
        assert!(
            sketch_miss
                .metagenomes
                .iter()
                .all(|metagenome| metagenome.shared_hashes == 0)
        );

        std::fs::remove_file(bgzf_path2).unwrap();
        let capped = engine
            .search(
                "capped",
                sequence.as_bytes(),
                TraceConfig {
                    max_metagenomes: 1,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(capped.metagenomes.len(), 1);
        assert_eq!(capped.metagenomes[0].name, "sample");
        assert_eq!(
            capped.completion,
            SearchCompletion::CandidateBudgetExceeded {
                candidates_omitted: 1
            }
        );
        #[cfg(unix)]
        {
            let cached = engine
                .candidate_census(&prepared, direct_config, full_cache_bytes)
                .unwrap();
            let times = std::fs::FileTimes::new().set_modified(
                std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(42),
            );
            std::fs::OpenOptions::new()
                .write(true)
                .open(jidx)
                .unwrap()
                .set_times(times)
                .unwrap();
            assert!(matches!(
                engine.trace_selected(
                    &prepared,
                    cached.candidates,
                    &rare_first,
                    cached.lookups.as_ref(),
                    direct_config,
                ),
                Err(TraceError::Invalid("cached seed lookup identity"))
            ));
        }
    }

    #[test]
    fn maps_fully_wrapped_query_intervals() {
        assert_eq!(
            query_segments(120, Interval::new(10, 20).unwrap(), 128, true).unwrap(),
            vec![Interval::new(2, 12).unwrap()]
        );
    }

    #[test]
    fn retains_every_repeated_query_seed_position() {
        let seeds = query_seeds(b"AAAAAA", 3, false, false).unwrap();
        assert_eq!(
            seeds.iter().map(|seed| seed.position).collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );
        assert!(
            seeds
                .iter()
                .all(|seed| seed.packed_key == seeds[0].packed_key)
        );
    }

    #[test]
    fn groups_hits_across_diagonal_boundaries_and_small_indels() {
        let key = RegionKey {
            metagenome_id: 0,
            contig_id: 0,
            strand: Strand::Forward,
            k: 21,
        };
        let grouped = form_regions(
            BTreeMap::from([(
                key,
                vec![
                    SeedHit {
                        query: 0,
                        target: 63,
                        diagonal: 63,
                    },
                    SeedHit {
                        query: 20,
                        target: 85,
                        diagonal: 65,
                    },
                ],
            )]),
            4,
        );
        assert_eq!(grouped.len(), 1);
        assert_eq!(grouped[0].0, key);
        assert_eq!(grouped[0].1.hits, 2);
    }

    #[test]
    fn separates_hits_beyond_diagonal_drift() {
        let key = RegionKey {
            metagenome_id: 0,
            contig_id: 0,
            strand: Strand::Forward,
            k: 21,
        };
        let grouped = form_regions(
            BTreeMap::from([(
                key,
                vec![
                    SeedHit {
                        query: 0,
                        target: 10,
                        diagonal: 10,
                    },
                    SeedHit {
                        query: 20,
                        target: 40,
                        diagonal: 20,
                    },
                ],
            )]),
            4,
        );
        assert_eq!(grouped.len(), 2);
    }

    #[test]
    fn coalesces_only_overlapping_or_adjacent_target_spans() {
        let mut spans = vec![(100, 110), (15, 20), (0, 10), (8, 15), (200, 205)];
        coalesce_spans(&mut spans);
        assert_eq!(spans, vec![(0, 20), (100, 110), (200, 205)]);
    }

    #[test]
    fn extracts_every_tagged_k15_rescue_position() {
        let seeds = query_seeds(b"ACGTACGTACGTACGTACGTA", 21, true, false).unwrap();
        assert_eq!(
            seeds
                .iter()
                .filter(|seed| seed.packed_key & RESCUE_K15_TAG != 0)
                .map(|seed| seed.position)
                .collect::<Vec<_>>(),
            (0..7).collect::<Vec<_>>()
        );
        assert_eq!(
            seeds
                .iter()
                .filter(|seed| seed.packed_key & RESCUE_K15_TAG == 0)
                .count(),
            1
        );
        assert_eq!(minimum_region_hits(15, 2), 3);
        assert_eq!(minimum_region_hits(21, 2), 2);
    }
}
