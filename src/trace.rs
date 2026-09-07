use crate::alignment::{AlignmentConfig, AlignmentError, AlignmentWorkspace, Interval, Strand};
use crate::bgzf::{BgzfError, BgzfReader};
use crate::jidx::{JidxError, RESCUE_K15_TAG, seed_length, sha256_reader};
use crate::jidx_reader::{ContigId, JidxReader, JidxReaderError, MetagenomeId};
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
            alignment: AlignmentConfig::default(),
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

    pub fn search(
        &self,
        query_id: impl Into<String>,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<TraceResult, TraceError> {
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
        let sketch_candidates = self.screen_candidates(&query_id, &query, config)?;
        let query_seeds = query_seeds(
            &query,
            self.index.header().k,
            self.index.header().rescue_k15,
            config.circular,
        )?;
        let mut query_seeds_by_key = BTreeMap::<u64, Vec<QuerySeed>>::new();
        for seed in query_seeds {
            query_seeds_by_key
                .entry(seed.packed_key)
                .or_default()
                .push(seed);
        }
        let mut region_hits = BTreeMap::<RegionKey, Vec<SeedHit>>::new();
        let mut seed_hits = HashMap::<MetagenomeId, u64>::new();
        for (packed_key, query_seeds) in query_seeds_by_key {
            let Some(index_seed) = self.index.find_seed(packed_key)? else {
                continue;
            };
            let seed_k = seed_length(
                self.index.header().k,
                self.index.header().rescue_k15,
                packed_key,
            )?;
            let occurrences = self.index.seed_occurrences(index_seed)?;
            for seed in query_seeds {
                for occurrence in &occurrences {
                    let contig = self
                        .index
                        .contig(occurrence.contig_id)?
                        .ok_or(TraceError::Invalid("missing occurrence contig"))?;
                    let strand = if seed.canonical_orientation == occurrence.canonical_orientation {
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
                    let diagonal = i128::from(oriented_position) - i128::from(seed.position);
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
                    let hits = seed_hits.entry(contig.metagenome_id).or_default();
                    *hits = hits.saturating_add(1);
                }
            }
        }

        let (candidates, completion) =
            self.rank_candidates(sketch_candidates, &seed_hits, config.max_metagenomes)?;
        let candidate_ids: HashSet<_> = candidates.iter().map(|candidate| candidate.id).collect();
        region_hits.retain(|key, _| candidate_ids.contains(&key.metagenome_id));
        let regions = form_regions(region_hits, config.diagonal_bin_bases);

        let tasks = self.tasks(regions, query_length, config)?;
        let (loaded, reads) = self.load_ranges(&tasks, config.verify_resources)?;
        let fragments = self.align_tasks(&query, &tasks, &loaded, config)?;
        let mut by_metagenome = BTreeMap::<MetagenomeId, Vec<Fragment>>::new();
        for (metagenome_id, fragment) in fragments {
            by_metagenome
                .entry(metagenome_id)
                .or_default()
                .push(fragment);
        }

        let candidates_screened =
            u32::try_from(candidates.len()).map_err(|_| TraceError::Invalid("candidate count"))?;
        let mut metagenomes = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            let fragments = by_metagenome.remove(&candidate.id).unwrap_or_default();
            let mosaic = build_mosaic(query_length, &fragments)?;
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
                exact_seed_hits: seed_hits.get(&candidate.id).copied().unwrap_or(0),
                compressed_bytes_read: stats.bytes_read,
                range_requests: stats.read_requests,
                bgzf_blocks_decoded,
                contigs,
                mosaic,
            });
        }
        Ok(TraceResult {
            query_id,
            query_length,
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

    fn rank_candidates(
        &self,
        sketch_candidates: Vec<Candidate>,
        seed_hits: &HashMap<MetagenomeId, u64>,
        max_metagenomes: usize,
    ) -> Result<(Vec<Candidate>, SearchCompletion), TraceError> {
        let mut candidates = BTreeMap::new();
        for candidate in sketch_candidates {
            candidates.insert(candidate.id, candidate);
        }
        for (&id, &exact_seed_hits) in seed_hits {
            if let Some(candidate) = candidates.get_mut(&id) {
                candidate.exact_seed_hits = exact_seed_hits;
                continue;
            }
            let name = self
                .index
                .metagenome_name(id)?
                .ok_or(TraceError::Invalid("JIDX metagenome ID"))?
                .to_string();
            candidates.insert(
                id,
                Candidate {
                    id,
                    name,
                    shared_hashes: 0,
                    containment: 0.0,
                    exact_seed_hits,
                },
            );
        }
        let mut candidates = candidates.into_values().collect::<Vec<_>>();
        candidates.sort_by(|left, right| {
            right
                .exact_seed_hits
                .cmp(&left.exact_seed_hits)
                .then_with(|| right.containment.total_cmp(&left.containment))
                .then_with(|| right.shared_hashes.cmp(&left.shared_hashes))
                .then_with(|| left.name.cmp(&right.name))
        });
        let completion = candidate_completion(candidates.len(), max_metagenomes)?;
        candidates.truncate(max_metagenomes);
        Ok((candidates, completion))
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
        tasks
            .par_iter()
            .map_init(AlignmentWorkspace::default, |workspace, task| {
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

struct Candidate {
    id: MetagenomeId,
    name: String,
    shared_hashes: u32,
    containment: f64,
    exact_seed_hits: u64,
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

fn digest_hex(digest: [u8; 32]) -> String {
    let mut output = String::with_capacity(64);
    for byte in digest {
        write!(&mut output, "{byte:02x}").expect("writing to a string cannot fail");
    }
    output
}

fn candidate_completion(
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
    #[error(transparent)]
    Mosaic(#[from] MosaicError),
    #[error("invalid trace input: {0}")]
    Invalid(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::handlers::{TraceArgs, handle_trace_command};
    use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
    use crate::writer::{BuildConfig, build};
    use noodles_bgzf::{self as bgzf, gzi};
    use std::io::Write;

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
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}},{{\"name\":\"sample2\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display(),
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
            query,
            database: jam,
            index: jidx,
            manifest,
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
        assert_eq!(
            capped.completion,
            SearchCompletion::CandidateBudgetExceeded {
                candidates_omitted: 1
            }
        );

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
