use crate::jidx::sha256;
use crate::range_source::S3Config;
use crate::trace::{
    CachedSeedLookups, Candidate, LOOKUP_CACHE_BYTES, MetagenomeTrace, PreparedQuery,
    SearchCompletion, TraceConfig, TraceEngine, TraceError, candidate_completion,
    compare_candidates, digest_hex, prepare_query,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::path::{Path, PathBuf};

#[derive(Deserialize)]
struct CollectionManifest {
    version: u32,
    shards: Vec<CollectionShard>,
}

#[derive(Deserialize)]
struct CollectionShard {
    database: PathBuf,
    index: PathBuf,
    manifest: PathBuf,
    index_header_sha256: String,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct CollectionTraceResult {
    pub root_sha256: String,
    pub query_id: String,
    pub query_length: u64,
    pub completion: SearchCompletion,
    pub candidates_screened: u32,
    pub metagenomes: Vec<CollectionMetagenomeTrace>,
}

#[derive(Debug, PartialEq, Serialize)]
pub struct CollectionMetagenomeTrace {
    pub shard_ordinal: u32,
    pub trace: MetagenomeTrace,
}

pub struct CollectionTraceEngine {
    root_sha256: String,
    shards: Vec<CollectionShard>,
    k: u8,
    rescue_k15: bool,
    s3: Option<S3Config>,
}

struct PlannedQuery {
    prepared: PreparedQuery,
    completion: SearchCompletion,
    candidates_screened: u32,
    key_order: Vec<u64>,
    selected: BTreeMap<usize, Vec<Candidate>>,
    lookups: Vec<Option<CachedSeedLookups>>,
    metagenomes: Vec<CollectionMetagenomeTrace>,
}

impl CollectionTraceEngine {
    pub fn open(root: impl AsRef<Path>, s3: Option<S3Config>) -> Result<Self, TraceError> {
        let root = root.as_ref();
        let bytes = std::fs::read(root)?;
        let mut manifest: CollectionManifest = serde_json::from_slice(&bytes)
            .map_err(|_| TraceError::Invalid("collection manifest JSON"))?;
        if manifest.version != 1 || manifest.shards.is_empty() {
            return Err(TraceError::Invalid("collection version or shard count"));
        }
        u32::try_from(manifest.shards.len())
            .map_err(|_| TraceError::Invalid("collection shard count"))?;
        let parent = root.parent().unwrap_or_else(|| Path::new("."));
        let mut names = HashSet::new();
        let mut seed_contract = None;
        for shard in &mut manifest.shards {
            shard.database = parent.join(&shard.database);
            shard.index = parent.join(&shard.index);
            shard.manifest = parent.join(&shard.manifest);
            let engine = open_shard(shard, s3.clone())?;
            let header = engine.index().header();
            let contract = (
                header.k,
                header.rescue_k15,
                header.minimizer_window,
                header.seed_scheme,
            );
            if seed_contract.is_some_and(|expected| expected != contract) {
                return Err(TraceError::Invalid("collection seed contracts differ"));
            }
            seed_contract = Some(contract);
            for id in 0..header.document_count {
                let name = engine
                    .index()
                    .metagenome_name(id)?
                    .ok_or(TraceError::Invalid("missing collection metagenome"))?;
                if !names.insert(name.to_string()) {
                    return Err(TraceError::Invalid("duplicate collection metagenome"));
                }
            }
        }
        let (k, rescue_k15, _, _) = seed_contract.expect("nonempty collection");
        Ok(Self {
            root_sha256: digest_hex(sha256(&bytes)),
            shards: manifest.shards,
            k,
            rescue_k15,
            s3,
        })
    }

    pub fn verify_index(&self) -> Result<(), TraceError> {
        self.shards
            .par_iter()
            .try_for_each(|shard| open_shard(shard, self.s3.clone())?.verify_index())
    }

    pub fn search(
        &self,
        query_id: impl Into<String>,
        sequence: &[u8],
        config: TraceConfig,
    ) -> Result<CollectionTraceResult, TraceError> {
        let prepared = prepare_query(query_id, sequence, config, self.k, self.rescue_k15)?;
        self.search_prepared(vec![prepared], config)?
            .pop()
            .ok_or(TraceError::Invalid("empty collection query batch"))
    }

    pub(crate) fn search_batch(
        &self,
        queries: &[(String, Vec<u8>)],
        config: TraceConfig,
    ) -> Result<Vec<CollectionTraceResult>, TraceError> {
        if queries.is_empty() {
            return Ok(Vec::new());
        }
        let prepared = queries
            .par_iter()
            .map(|(id, sequence)| {
                prepare_query(id.as_str(), sequence, config, self.k, self.rescue_k15)
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.search_prepared(prepared, config)
    }

    fn search_prepared(
        &self,
        prepared: Vec<PreparedQuery>,
        config: TraceConfig,
    ) -> Result<Vec<CollectionTraceResult>, TraceError> {
        let worker_limit = rayon::current_num_threads().max(1);
        let live_shards = worker_limit.div_ceil(prepared.len()).max(1);
        let cache_bytes = LOOKUP_CACHE_BYTES / prepared.len() / self.shards.len();
        let mut lookups = (0..prepared.len())
            .map(|_| (0..self.shards.len()).map(|_| None).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let mut frequencies = (0..prepared.len())
            .map(|_| BTreeMap::<u64, u64>::new())
            .collect::<Vec<_>>();
        let mut candidates = (0..prepared.len()).map(|_| Vec::new()).collect::<Vec<_>>();
        for (chunk, shards) in self.shards.chunks(live_shards).enumerate() {
            let shard_censuses = shards
                .par_iter()
                .map(|shard| {
                    let engine = open_shard(shard, self.s3.clone())?;
                    prepared
                        .par_iter()
                        .map(|query| engine.candidate_census(query, config, cache_bytes))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Vec<_>, TraceError>>()?;
            for (offset, shard) in shard_censuses.into_iter().enumerate() {
                let ordinal = chunk * live_shards + offset;
                for (query, census) in shard.into_iter().enumerate() {
                    lookups[query][ordinal] = census.lookups;
                    for (key, frequency) in census.frequencies {
                        let total = frequencies[query].entry(key).or_default();
                        *total = total
                            .checked_add(u64::from(frequency))
                            .ok_or(TraceError::Invalid("collection document frequency"))?;
                    }
                    candidates[query].extend(
                        census
                            .candidates
                            .into_iter()
                            .map(|candidate| (ordinal, candidate)),
                    );
                }
            }
        }

        let mut plans = prepared
            .into_iter()
            .zip(frequencies)
            .zip(candidates)
            .zip(lookups)
            .map(|(((query, frequencies), candidates), lookups)| {
                plan_query(query, frequencies, candidates, lookups, config)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let selected_shards = plans
            .iter()
            .flat_map(|plan| plan.selected.keys().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        for ordinals in selected_shards.chunks(live_shards) {
            let work = ordinals
                .iter()
                .map(|&ordinal| {
                    let queries = plans
                        .iter_mut()
                        .enumerate()
                        .filter_map(|(query, plan)| {
                            plan.selected
                                .remove(&ordinal)
                                .map(|candidates| (query, candidates, plan.lookups[ordinal].take()))
                        })
                        .collect::<Vec<_>>();
                    (ordinal, queries)
                })
                .collect::<Vec<_>>();
            let traces = work
                .into_par_iter()
                .map(|(ordinal, queries)| {
                    let engine = open_shard(&self.shards[ordinal], self.s3.clone())?;
                    queries
                        .into_par_iter()
                        .map(|(query, candidates, lookups)| {
                            let traces = engine.trace_selected(
                                &plans[query].prepared,
                                candidates,
                                &plans[query].key_order,
                                lookups.as_ref(),
                                config,
                            )?;
                            Ok((
                                query,
                                traces
                                    .into_iter()
                                    .map(|trace| CollectionMetagenomeTrace {
                                        shard_ordinal: ordinal as u32,
                                        trace,
                                    })
                                    .collect::<Vec<_>>(),
                            ))
                        })
                        .collect::<Result<Vec<_>, TraceError>>()
                })
                .collect::<Result<Vec<_>, TraceError>>()?;
            for (query, traces) in traces.into_iter().flatten() {
                plans[query].metagenomes.extend(traces);
            }
        }

        Ok(plans
            .into_iter()
            .map(|mut plan| {
                sort_metagenomes(&mut plan.metagenomes);
                CollectionTraceResult {
                    root_sha256: self.root_sha256.clone(),
                    query_id: plan.prepared.query_id,
                    query_length: plan.prepared.query_length,
                    completion: plan.completion,
                    candidates_screened: plan.candidates_screened,
                    metagenomes: plan.metagenomes,
                }
            })
            .collect())
    }
}

fn plan_query(
    prepared: PreparedQuery,
    frequencies: BTreeMap<u64, u64>,
    mut candidates: Vec<(usize, Candidate)>,
    mut lookups: Vec<Option<CachedSeedLookups>>,
    config: TraceConfig,
) -> Result<PlannedQuery, TraceError> {
    candidates.sort_by(|(left_shard, left), (right_shard, right)| {
        compare_candidates(left, right)
            .then_with(|| left_shard.cmp(right_shard))
            .then_with(|| left.id.cmp(&right.id))
    });
    let completion = candidate_completion(candidates.len(), config.max_metagenomes)?;
    candidates.truncate(config.max_metagenomes);
    let candidates_screened = u32::try_from(candidates.len())
        .map_err(|_| TraceError::Invalid("collection candidate count"))?;
    let mut selected = BTreeMap::<usize, Vec<Candidate>>::new();
    for (ordinal, candidate) in candidates {
        selected.entry(ordinal).or_default().push(candidate);
    }
    for (ordinal, lookups) in lookups.iter_mut().enumerate() {
        if !selected.contains_key(&ordinal) {
            *lookups = None;
        }
    }
    let mut key_order = frequencies.into_iter().collect::<Vec<_>>();
    key_order.sort_unstable_by_key(|&(key, frequency)| (frequency, key));
    Ok(PlannedQuery {
        prepared,
        completion,
        candidates_screened,
        key_order: key_order.into_iter().map(|(key, _)| key).collect(),
        selected,
        lookups,
        metagenomes: Vec::new(),
    })
}

fn sort_metagenomes(metagenomes: &mut [CollectionMetagenomeTrace]) {
    metagenomes.sort_by(|left, right| {
        right
            .trace
            .exact_seed_hits
            .cmp(&left.trace.exact_seed_hits)
            .then_with(|| right.trace.containment.total_cmp(&left.trace.containment))
            .then_with(|| right.trace.shared_hashes.cmp(&left.trace.shared_hashes))
            .then_with(|| left.trace.name.cmp(&right.trace.name))
            .then_with(|| left.shard_ordinal.cmp(&right.shard_ordinal))
            .then_with(|| left.trace.metagenome_id.cmp(&right.trace.metagenome_id))
    });
}

fn open_shard(shard: &CollectionShard, s3: Option<S3Config>) -> Result<TraceEngine, TraceError> {
    let engine = TraceEngine::open(&shard.database, &shard.index, &shard.manifest, s3)?;
    // Header decoding requires canonical reserved bytes, so re-encoding preserves its exact bytes.
    if digest_hex(sha256(&engine.index().header().encode()?)) != shard.index_header_sha256 {
        return Err(TraceError::Invalid("collection JIDX header checksum"));
    }
    Ok(engine)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Cli;
    use crate::cli::handlers::{TraceArgs, TraceInput, handle_trace_command};
    use crate::jidx_builder::{JidxBuildConfig, build_local_jidx};
    use crate::jidx_reader::JidxReader;
    use crate::writer::{BuildConfig, build};
    use clap::Parser;
    use noodles_bgzf::{self as bgzf, gzi};
    use serde_json::{Value, json};
    use std::fs::File;
    use std::io::Write;

    fn shard(directory: &Path, names: &[&str], sequence: &str, window: u16) -> Value {
        let mut inputs = Vec::new();
        let mut metagenomes = Vec::new();
        for name in names {
            let fasta = directory.join(format!("{name}.fa"));
            let raw = format!(">contig\n{sequence}\n");
            std::fs::write(&fasta, format!(">{name}\n{sequence}\n")).unwrap();
            inputs.push(fasta);
            let bgzf_path = directory.join(format!("{name}.bgz"));
            let fai_path = directory.join(format!("{name}.bgz.fai"));
            let gzi_path = directory.join(format!("{name}.bgz.gzi"));
            let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
            writer.write_all(raw.as_bytes()).unwrap();
            writer.finish().unwrap();
            std::fs::write(
                &fai_path,
                format!(
                    "contig\t{}\t8\t{}\t{}\n",
                    sequence.len(),
                    sequence.len(),
                    sequence.len() + 1
                ),
            )
            .unwrap();
            gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
            metagenomes
                .push(json!({"name": name, "bgzf": bgzf_path, "fai": fai_path, "gzi": gzi_path}));
        }
        let jam = directory.join(format!("{}.jam", names[0]));
        build(
            &inputs,
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
        let manifest = directory.join(format!("{}.json", names[0]));
        std::fs::write(
            &manifest,
            serde_json::to_vec(&json!({"metagenomes": metagenomes})).unwrap(),
        )
        .unwrap();
        let index = directory.join(format!("{}.jidx", names[0]));
        build_local_jidx(
            &jam,
            &manifest,
            &index,
            JidxBuildConfig {
                k: 5,
                minimizer_window: window,
                rescue_k15: false,
            },
        )
        .unwrap();
        let header = JidxReader::open(&index).unwrap().header().encode().unwrap();
        json!({
            "database": jam.file_name().unwrap().to_str().unwrap(),
            "index": index.file_name().unwrap().to_str().unwrap(),
            "manifest": manifest.file_name().unwrap().to_str().unwrap(),
            "index_header_sha256": digest_hex(sha256(&header)),
        })
    }

    #[test]
    fn collection_merges_complete_shard_searches_and_caps_once() {
        let directory = tempfile::tempdir().unwrap();
        let mut state = 127u64;
        let sequence = (0..192)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect::<String>();
        let shards = [
            shard(directory.path(), &["zulu", "yankee"], &sequence, 4),
            shard(directory.path(), &["alpha", "bravo"], &sequence, 4),
            shard(directory.path(), &["charlie"], &sequence.repeat(2), 4),
        ];
        let root = directory.path().join("collection.json");
        let root_bytes = serde_json::to_vec(&json!({"version": 1, "shards": shards})).unwrap();
        std::fs::write(&root, &root_bytes).unwrap();
        let engine = CollectionTraceEngine::open(&root, None).unwrap();
        engine.verify_index().unwrap();
        let config = TraceConfig {
            use_sketch: false,
            circular: false,
            ..TraceConfig::default()
        };
        let first = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| engine.search("query", sequence.as_bytes(), config).unwrap());
        assert_eq!(first.root_sha256, digest_hex(sha256(&root_bytes)));
        assert_eq!(first.completion, SearchCompletion::Complete);
        assert_eq!(first.candidates_screened, 5);
        assert_eq!(
            first
                .metagenomes
                .iter()
                .map(|result| result.trace.name.as_str())
                .collect::<Vec<_>>(),
            ["charlie", "alpha", "bravo", "yankee", "zulu"]
        );
        for (ordinal, shard) in engine.shards.iter().enumerate() {
            let single = open_shard(shard, None)
                .unwrap()
                .search("query", sequence.as_bytes(), config)
                .unwrap();
            let merged = first
                .metagenomes
                .iter()
                .filter(|result| result.shard_ordinal == ordinal as u32)
                .map(|result| result.trace.clone())
                .collect::<Vec<_>>();
            assert_eq!(single.metagenomes, merged);
            assert!(
                merged
                    .iter()
                    .all(|trace| trace.mosaic.covered_bases == sequence.len() as u64)
            );
        }
        let parallel = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap()
            .install(|| engine.search("query", sequence.as_bytes(), config).unwrap());
        assert_eq!(first, parallel);
        let capped = engine
            .search(
                "query",
                sequence.as_bytes(),
                TraceConfig {
                    max_metagenomes: 1,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(
            capped.completion,
            SearchCompletion::CandidateBudgetExceeded {
                candidates_omitted: 4
            }
        );
        assert_eq!(capped.metagenomes.len(), 1);
        assert_eq!(capped.metagenomes[0].trace.name, "charlie");
        assert_eq!(capped.metagenomes[0].shard_ordinal, 2);
        let partial_query = format!("{sequence}TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        let sketch_miss = engine
            .search(
                "partial",
                partial_query.as_bytes(),
                TraceConfig {
                    use_sketch: true,
                    min_containment: 1.0,
                    ..config
                },
            )
            .unwrap();
        assert_eq!(sketch_miss.candidates_screened, 5);
        assert!(
            sketch_miss
                .metagenomes
                .iter()
                .all(|result| result.trace.shared_hashes == 0 && result.trace.exact_seed_hits > 0)
        );
        let negative = engine
            .search(
                "negative",
                b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
                config,
            )
            .unwrap();
        assert_eq!(negative.completion, SearchCompletion::Complete);
        assert_eq!(negative.candidates_screened, 0);

        let batch_queries = vec![
            ("query".to_string(), sequence.as_bytes().to_vec()),
            ("partial".to_string(), partial_query.as_bytes().to_vec()),
            (
                "negative".to_string(),
                b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_vec(),
            ),
        ];
        let singles = batch_queries
            .iter()
            .map(|(id, sequence)| engine.search(id.as_str(), sequence, config).unwrap())
            .collect::<Vec<_>>();
        for threads in [1, 2, 4] {
            let batch = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| engine.search_batch(&batch_queries, config).unwrap());
            assert_eq!(batch, singles);
        }

        let query = directory.path().join("query.fa");
        std::fs::write(
            &query,
            format!(">query\n{sequence}\n>negative\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"),
        )
        .unwrap();
        let output = directory.path().join("results.jsonl");
        handle_trace_command(TraceArgs {
            query: query.clone(),
            input: TraceInput::Collection(root.clone()),
            audit_index: false,
            output: output.clone(),
            query_id: None,
            config,
            s3: None,
            force: false,
        })
        .unwrap();
        let records = std::fs::read_to_string(&output)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str::<Value>(line).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(
            records,
            vec![
                serde_json::to_value(&first).unwrap(),
                serde_json::to_value(&negative).unwrap()
            ]
        );
        std::fs::remove_file(&engine.shards[1].index).unwrap();
        assert!(
            engine
                .search(
                    "negative",
                    b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
                    config
                )
                .is_err()
        );
        assert!(
            handle_trace_command(TraceArgs {
                query,
                input: TraceInput::Collection(root),
                audit_index: false,
                output: output.clone(),
                query_id: None,
                config,
                s3: None,
                force: true
            })
            .is_err()
        );
        assert_eq!(std::fs::read_to_string(output).unwrap().lines().count(), 2);
    }

    #[test]
    fn batch_applies_global_caps_independently_per_query() {
        let directory = tempfile::tempdir().unwrap();
        let sequence = |mut state: u64| {
            (0..192)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[(state >> 62) as usize] as char
                })
                .collect::<String>()
        };
        let first = sequence(17);
        let second = sequence(91);
        let shards = [
            shard(directory.path(), &["alpha"], &first, 4),
            shard(directory.path(), &["bravo"], &second, 4),
            shard(directory.path(), &["zulu"], &first, 4),
            shard(directory.path(), &["yankee"], &second, 4),
        ];
        let root = directory.path().join("collection.json");
        std::fs::write(
            &root,
            serde_json::to_vec(&json!({"version": 1, "shards": shards})).unwrap(),
        )
        .unwrap();
        let engine = CollectionTraceEngine::open(root, None).unwrap();
        let config = TraceConfig {
            use_sketch: false,
            circular: false,
            max_metagenomes: 1,
            ..TraceConfig::default()
        };
        let queries = vec![
            ("first".to_string(), first.into_bytes()),
            ("second".to_string(), second.into_bytes()),
        ];
        let batch = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| engine.search_batch(&queries, config).unwrap());
        let singles = queries
            .iter()
            .map(|(id, sequence)| engine.search(id.as_str(), sequence, config).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(batch, singles);
        assert_eq!(batch[0].metagenomes[0].trace.name, "alpha");
        assert_eq!(batch[1].metagenomes[0].trace.name, "bravo");
        assert!(batch.iter().all(|result| matches!(
            result.completion,
            SearchCompletion::CandidateBudgetExceeded {
                candidates_omitted: 3
            }
        )));
    }

    #[test]
    fn collection_rejects_duplicate_members_seed_mismatch_and_changed_header() {
        let directory = tempfile::tempdir().unwrap();
        let first = shard(directory.path(), &["first"], "ACGTTGCAACGT", 4);
        let second = shard(directory.path(), &["second"], "ACGTTGCAACGT", 8);
        let root = directory.path().join("collection.json");
        for shards in [
            vec![first.clone(), first.clone()],
            vec![first.clone(), second],
        ] {
            std::fs::write(
                &root,
                serde_json::to_vec(&json!({"version": 1, "shards": shards})).unwrap(),
            )
            .unwrap();
            assert!(CollectionTraceEngine::open(&root, None).is_err());
        }
        let mut changed = first;
        changed["index_header_sha256"] = json!("00".repeat(32));
        std::fs::write(
            &root,
            serde_json::to_vec(&json!({"version": 1, "shards": [changed]})).unwrap(),
        )
        .unwrap();
        assert!(matches!(
            CollectionTraceEngine::open(&root, None),
            Err(TraceError::Invalid("collection JIDX header checksum"))
        ));
    }

    #[test]
    fn trace_cli_requires_one_complete_input() {
        for arguments in [
            vec![
                "jam",
                "trace",
                "-q",
                "q.fa",
                "-o",
                "o.jsonl",
                "--collection",
                "root.json",
            ],
            vec![
                "jam", "trace", "-q", "q.fa", "-o", "o.jsonl", "-d", "d.jam", "-i", "d.jidx", "-M",
                "m.json",
            ],
        ] {
            assert!(Cli::try_parse_from(arguments).is_ok());
        }
        for arguments in [
            vec!["jam", "trace", "-q", "q.fa", "-o", "o.jsonl"],
            vec!["jam", "trace", "-q", "q.fa", "-o", "o.jsonl", "-d", "d.jam"],
            vec![
                "jam",
                "trace",
                "-q",
                "q.fa",
                "-o",
                "o.jsonl",
                "--collection",
                "root.json",
                "-i",
                "d.jidx",
            ],
        ] {
            assert!(Cli::try_parse_from(arguments).is_err());
        }
    }
}
