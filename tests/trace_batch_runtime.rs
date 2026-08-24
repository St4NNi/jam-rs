use jam_rs::archive::{
    ArchiveContig, ArchiveError, ArchiveMetadata, ArchiveMetrics, ArchiveResult, SeedKey,
    SeedLookupResult, SeedMatch, SeedOccurrence, SeedSchemeDescriptor, SeedSchemeId,
    SequenceRequest, SequenceSlice, TraceArchive,
};
use jam_rs::router::{
    CandidateWindowEvidence, HashAlgorithmId, PositionalWitnessOccurrence, RoutedCandidate,
    SharedWitness, WitnessClass, WitnessKey, WitnessScheme,
};
use jam_rs::trace::batch::runtime::{BatchArchiveFactory, JmaBatchRuntime, RuntimeTraceResult};
use jam_rs::trace::batch::{
    BatchCandidateReadPlan, BatchExecutionOptions, BatchKeyPageRequest, BatchOutputOrder,
    BatchQueryRecord, BatchResultRecord, BatchSequenceBlockRequest, execute_batch,
    prepare_batch_queries,
};
use jam_rs::trace::contracts::{BatchCandidateWork, WitnessPlanningRequest};
use jam_rs::trace::model::{QueryDescriptor, QueryKind, TopologyRequested};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

#[derive(Default)]
struct FakeState {
    seed_calls: Mutex<Vec<Vec<SeedKey>>>,
    sequence_calls: Mutex<Vec<Vec<SequenceRequest>>>,
    opens: Mutex<Vec<String>>,
}

#[derive(Clone)]
struct FakeArchive {
    state: Arc<FakeState>,
    matches: BTreeMap<SeedKey, Vec<SeedOccurrence>>,
}

impl TraceArchive for FakeArchive {
    fn metadata(&self) -> ArchiveResult<ArchiveMetadata> {
        Ok(ArchiveMetadata {
            format_identifier: "fake-jma".to_string(),
            format_version: 1,
            layout_identifier: 1,
            source_assembly_sha256: [0; 32],
            archive_sha256: None,
            builder_version: "test".to_string(),
            source_commit: "test".to_string(),
            hash_algorithm: "jamhash_u64_v1".to_string(),
            total_bases: 10_000,
            contigs: vec![ArchiveContig {
                id: 0,
                name: "contig-0".to_string(),
                length: 10_000,
            }],
        })
    }

    fn available_seed_schemes(&self) -> ArchiveResult<Vec<SeedSchemeDescriptor>> {
        Ok(vec![SeedSchemeDescriptor {
            scheme_id: 1,
            algorithm_id: 1,
            span: 21,
            informative_bases: 21,
            density_parameter: 20,
            bucket_bits: 16,
            key_encoding: 1,
            occurrence_encoding: 1,
            flags: 0,
        }])
    }

    fn lookup_seeds(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<SeedLookupResult> {
        self.lookup_seeds_bounded(scheme, keys, None)
    }

    fn lookup_seeds_bounded(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
        max_occurrences: Option<u32>,
    ) -> ArchiveResult<SeedLookupResult> {
        if scheme != SeedSchemeId(1) {
            return Err(ArchiveError::UnknownSeedScheme(scheme.0));
        }
        self.state.seed_calls.lock().unwrap().push(keys.to_vec());
        let matches = keys
            .iter()
            .enumerate()
            .filter_map(|(key_index, key)| {
                let mut occurrences = self.matches.get(key)?.clone();
                if let Some(limit) = max_occurrences {
                    occurrences.truncate(limit as usize);
                }
                Some(SeedMatch {
                    key_index,
                    occurrences,
                })
            })
            .collect();
        Ok(SeedLookupResult {
            matches,
            metrics: Default::default(),
        })
    }

    fn read_sequences(&self, requests: &[SequenceRequest]) -> ArchiveResult<Vec<SequenceSlice>> {
        self.state
            .sequence_calls
            .lock()
            .unwrap()
            .push(requests.to_vec());
        requests
            .iter()
            .copied()
            .map(|request| {
                Ok(SequenceSlice {
                    request,
                    bases: vec![b'A'; request.len() as usize],
                })
            })
            .collect()
    }

    fn metrics(&self) -> ArchiveMetrics {
        ArchiveMetrics::default()
    }
}

#[derive(Clone)]
struct FakeFactory {
    state: Arc<FakeState>,
    matches: BTreeMap<SeedKey, Vec<SeedOccurrence>>,
}

impl BatchArchiveFactory for FakeFactory {
    type Archive = FakeArchive;

    fn open(
        &mut self,
        metagenome_id: &str,
    ) -> Result<Self::Archive, jam_rs::trace::batch::runtime::RuntimeError> {
        self.state
            .opens
            .lock()
            .unwrap()
            .push(metagenome_id.to_string());
        Ok(FakeArchive {
            state: Arc::clone(&self.state),
            matches: self.matches.clone(),
        })
    }
}

fn scheme() -> WitnessScheme {
    WitnessScheme {
        scheme_id: 1,
        k: 21,
        base_scale: 20,
        available_scales: vec![20, 50, 100, 200, 500],
        hash_id: HashAlgorithmId::JamhashU64V1,
        zero_excluded: true,
    }
}

fn record(query_id: &str) -> BatchQueryRecord {
    let sequence = vec![b'G'; 256];
    BatchQueryRecord::new(
        QueryDescriptor {
            query_id: query_id.to_string(),
            query_length: sequence.len() as u64,
            query_kind: QueryKind::Plasmid,
            topology_requested: TopologyRequested::Linear,
        },
        sequence,
    )
}

fn seed_key(key: WitnessKey) -> SeedKey {
    let bytes = key.packed.to_be_bytes();
    SeedKey {
        digest: key.jamhash,
        verification: bytes[2..].to_vec(),
    }
}

fn shared(key: WitnessKey) -> SharedWitness {
    SharedWitness {
        key,
        query_position: 4,
        query_reverse: false,
        query_window_id: 0,
        document_frequency: 1,
        witness_tier: 20,
        witness_class: WitnessClass::Rare,
    }
}

fn work(query_index: u32, metagenome_id: &str, key: WitnessKey) -> BatchCandidateWork {
    let witness = shared(key);
    BatchCandidateWork {
        query_index,
        metagenome_id: metagenome_id.to_string(),
        candidate: RoutedCandidate {
            metagenome_id: metagenome_id.to_string(),
            rare_shared_witnesses: 1,
            common_shared_witnesses: 0,
            supported_query_windows: 1,
            weighted_witness_sum: 1.0,
            estimated_query_containment: 1.0,
            candidate_tier: 20,
            window_evidence: CandidateWindowEvidence {
                matched_query_windows: 1,
                total_eligible_query_windows: 1,
                longest_supported_window_run: 1,
                rare_witness_count: 1,
                common_witness_count: 0,
                total_shared_witness_count: 1,
            },
            shared_witnesses: vec![witness],
            positional_witnesses: vec![PositionalWitnessOccurrence {
                witness,
                contig_id: 0,
                position: 10,
                reverse: false,
                scheme_id: 1,
            }],
        },
    }
}

fn make_runtime(
    state: Arc<FakeState>,
    matches: BTreeMap<SeedKey, Vec<SeedOccurrence>>,
) -> JmaBatchRuntime<FakeFactory> {
    JmaBatchRuntime::new(FakeFactory { state, matches }, SeedSchemeId(1), Some(8)).unwrap()
}

fn setup() -> (
    jam_rs::trace::batch::BatchPlan,
    Vec<BatchCandidateWork>,
    WitnessKey,
    BTreeMap<SeedKey, Vec<SeedOccurrence>>,
) {
    let records = vec![record("q-00"), record("q-01"), record("q-02")];
    let plan = prepare_batch_queries(
        &records,
        &scheme(),
        jam_rs::trace::batch::BatchPrepareOptions::new(
            WitnessPlanningRequest {
                min_trace_length: 160,
                target_identity: 0.99,
                max_zero_witness_risk: 0.01,
                strict: false,
            },
            128,
        ),
    )
    .unwrap();
    let key = WitnessKey::from_packed(0x12345).unwrap();
    let occurrence = SeedOccurrence {
        contig_id: 0,
        position: 10,
        span: 21,
        reverse: false,
    };
    let mut matches = BTreeMap::new();
    matches.insert(seed_key(key), vec![occurrence]);
    let works = vec![
        work(0, "mg-a", key),
        work(1, "mg-a", key),
        work(2, "mg-b", key),
    ];
    (plan, works, key, matches)
}

#[test]
fn grouped_jma_runtime_opens_once_and_unions_public_reads() {
    let (plan, works, _key, matches) = setup();
    let state = Arc::new(FakeState::default());
    let mut runtime = make_runtime(Arc::clone(&state), matches);
    let output = execute_batch(
        &plan,
        &works,
        &mut runtime,
        BatchExecutionOptions {
            output_order: BatchOutputOrder::Deterministic,
        },
    )
    .unwrap();

    assert_eq!(output.records.len(), 3);
    assert_eq!(output.metrics.metagenome_open_count, 2);
    assert_eq!(output.metrics.exact_key_requests, 2);
    assert_eq!(output.metrics.key_page_requests, 2);
    assert_eq!(output.metrics.key_pages_decoded, 2);
    assert_eq!(output.metrics.sequence_block_requests, 2);
    assert_eq!(output.metrics.sequence_blocks_decoded, 2);
    assert_eq!(
        state.opens.lock().unwrap().as_slice(),
        vec!["mg-a".to_string(), "mg-b".to_string()].as_slice()
    );
    assert_eq!(state.seed_calls.lock().unwrap().len(), 2);
    assert!(
        state
            .seed_calls
            .lock()
            .unwrap()
            .iter()
            .all(|call| call.len() == 1)
    );
    assert_eq!(state.sequence_calls.lock().unwrap().len(), 2);
    assert!(
        state
            .sequence_calls
            .lock()
            .unwrap()
            .iter()
            .all(|call| call.len() == 1)
    );
    assert!(output.records.iter().all(|record| {
        record.result.exact_key_matches == 1
            && record.result.occurrence_count == 1
            && record.result.sequence_bases == 21
            && record.key.query_id == record.result.query_id
            && record.key.metagenome_id == record.result.metagenome_id
    }));
}

#[test]
fn grouped_runtime_results_equal_separate_query_results() {
    let (plan, works, _key, matches) = setup();
    let grouped_state = Arc::new(FakeState::default());
    let mut grouped_runtime = make_runtime(Arc::clone(&grouped_state), matches.clone());
    let grouped = execute_batch(
        &plan,
        &works,
        &mut grouped_runtime,
        BatchExecutionOptions::default(),
    )
    .unwrap()
    .records;

    let separate_state = Arc::new(FakeState::default());
    let mut separate_records = Vec::<BatchResultRecord<RuntimeTraceResult>>::new();
    for selected in &works {
        let mut runtime = make_runtime(Arc::clone(&separate_state), matches.clone());
        separate_records.extend(
            execute_batch(
                &plan,
                std::slice::from_ref(selected),
                &mut runtime,
                BatchExecutionOptions::default(),
            )
            .unwrap()
            .records,
        );
    }
    separate_records.sort_by(|left, right| left.key.cmp(&right.key));
    assert_eq!(grouped, separate_records);
}

#[test]
fn explicit_plans_union_two_exact_keys_into_one_public_lookup() {
    let records = vec![record("q-00"), record("q-01")];
    let plan = prepare_batch_queries(
        &records,
        &scheme(),
        jam_rs::trace::batch::BatchPrepareOptions::new(
            WitnessPlanningRequest {
                min_trace_length: 160,
                target_identity: 0.99,
                max_zero_witness_risk: 0.01,
                strict: false,
            },
            128,
        ),
    )
    .unwrap();
    let first = WitnessKey::from_packed(0x23456).unwrap();
    let second = WitnessKey::from_packed(0x34567).unwrap();
    let works = vec![work(0, "mg-a", first), work(1, "mg-a", second)];
    let occurrence = SeedOccurrence {
        contig_id: 0,
        position: 10,
        span: 21,
        reverse: false,
    };
    let matches = BTreeMap::from([
        (seed_key(first), vec![occurrence]),
        (seed_key(second), vec![occurrence]),
    ]);
    let state = Arc::new(FakeState::default());
    let mut runtime = make_runtime(Arc::clone(&state), matches);
    let shared_plan = BatchCandidateReadPlan {
        key_pages: vec![
            BatchKeyPageRequest::new(17, first),
            BatchKeyPageRequest::new(17, second),
        ],
        sequence_blocks: vec![BatchSequenceBlockRequest::new(0, 10, 31)],
    };
    runtime.set_candidate_plan(0, "mg-a", shared_plan.clone());
    runtime.set_candidate_plan(1, "mg-a", shared_plan);

    let output = execute_batch(
        &plan,
        &works,
        &mut runtime,
        BatchExecutionOptions::default(),
    )
    .unwrap();
    assert_eq!(output.metrics.exact_key_requests, 2);
    assert_eq!(output.metrics.key_page_requests, 1);
    assert_eq!(output.metrics.key_pages_decoded, 1);
    assert_eq!(state.seed_calls.lock().unwrap().len(), 1);
    assert_eq!(state.seed_calls.lock().unwrap()[0].len(), 2);
    assert!(
        output
            .records
            .iter()
            .all(|record| record.result.exact_key_matches == 1)
    );
}
