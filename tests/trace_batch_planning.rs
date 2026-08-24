use jam_rs::router::{
    CandidateWindowEvidence, HashAlgorithmId, RoutedCandidate, WitnessKey, WitnessScheme,
};
use jam_rs::trace::batch::fixture::{FixtureTraceResult, InMemoryBatchAdapter};
use jam_rs::trace::batch::{
    BatchCandidateReadPlan, BatchExecutionOptions, BatchKeyPageRequest, BatchOutputOrder,
    BatchPlan, BatchPrepareOptions, BatchQueryRecord, BatchResultRecord, BatchRouter,
    BatchSequenceBlockRequest, execute_batch, group_candidate_work, prepare_batch_queries,
    route_batch, validate_batch_queries,
};
use jam_rs::trace::contracts::{BatchCandidateWork, WitnessPlanningRequest};
use jam_rs::trace::model::{QueryDescriptor, QueryKind, TopologyRequested};
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{self, Display, Formatter};

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

fn sequence(length: usize, seed: u64) -> Vec<u8> {
    let alphabet = b"ACGT";
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            alphabet[(state as usize) % alphabet.len()]
        })
        .collect()
}

fn record(query_id: &str, length: usize, seed: u64) -> BatchQueryRecord {
    let sequence = sequence(length, seed);
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

fn options() -> BatchPrepareOptions {
    BatchPrepareOptions::new(
        WitnessPlanningRequest {
            min_trace_length: 160,
            target_identity: 0.99,
            max_zero_witness_risk: 0.01,
            strict: false,
        },
        128,
    )
}

fn candidate(query_index: u32, metagenome_id: &str, key: WitnessKey) -> BatchCandidateWork {
    let candidate = RoutedCandidate {
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
        shared_witnesses: Vec::new(),
        positional_witnesses: Vec::new(),
    };
    let _ = key;
    BatchCandidateWork {
        query_index,
        metagenome_id: metagenome_id.to_string(),
        candidate,
    }
}

#[test]
fn batch_input_validation_happens_before_witness_work() {
    let duplicate = vec![record("same", 512, 1), record("same", 512, 2)];
    assert!(validate_batch_queries(&duplicate).is_err());

    let mut wrong = record("wrong", 512, 3);
    wrong.descriptor.query_length += 1;
    assert!(validate_batch_queries(&[wrong]).is_err());

    let plan = prepare_batch_queries(
        &[record("q0", 512, 4), record("q1", 512, 5)],
        &scheme(),
        options(),
    )
    .unwrap();
    assert_eq!(plan.query_count(), 2);
    assert_eq!(plan.queries[0].query_index, 0);
    assert_eq!(plan.queries[1].query_index, 1);
    assert!(!plan.queries[0].witnesses.is_empty());
    // The query witness vector is prepared once at the dense base scale;
    // expansion is only a view over those exact occurrences.
    assert!(plan.queries[0].witnesses.len() < plan.queries[0].sequence.len());
}

#[derive(Debug)]
struct RecordingRouterError;

impl Display for RecordingRouterError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str("recording router error")
    }
}

impl Error for RecordingRouterError {}

#[derive(Default)]
struct RecordingRouter {
    calls: Vec<(u32, usize)>,
    candidates: BTreeMap<u32, Vec<RoutedCandidate>>,
}

impl BatchRouter for RecordingRouter {
    type Error = RecordingRouterError;

    fn route_query(
        &mut self,
        query: &jam_rs::trace::contracts::PreparedBatchQuery,
        witnesses: &[jam_rs::router::search::TieredQueryWitness],
    ) -> Result<Vec<RoutedCandidate>, Self::Error> {
        self.calls.push((query.query_index, witnesses.len()));
        Ok(self
            .candidates
            .get(&query.query_index)
            .cloned()
            .unwrap_or_default())
    }
}

#[test]
fn route_batch_prepares_witnesses_once_and_attaches_query_indices() {
    let plan = prepare_batch_queries(
        &[record("q0", 512, 6), record("q1", 512, 7)],
        &scheme(),
        options(),
    )
    .unwrap();
    let key = WitnessKey::from_packed(101).unwrap();
    let mut router = RecordingRouter::default();
    router
        .candidates
        .insert(0, vec![candidate(0, "mg-a", key).candidate]);
    router
        .candidates
        .insert(1, vec![candidate(1, "mg-a", key).candidate]);

    let works = route_batch(&plan, &mut router).unwrap();
    assert_eq!(works.len(), 2);
    assert_eq!(works[0].query_index, 0);
    assert_eq!(works[1].query_index, 1);
    assert_eq!(router.calls.len(), 2);
    assert!(router.calls.iter().all(|(_, count)| *count > 0));
}

fn read_plan(key: WitnessKey, block: BatchSequenceBlockRequest) -> BatchCandidateReadPlan {
    BatchCandidateReadPlan {
        key_pages: vec![BatchKeyPageRequest::new(7, key)],
        sequence_blocks: vec![block],
    }
}

fn fixture_plan() -> (
    BatchPlan,
    Vec<BatchCandidateWork>,
    InMemoryBatchAdapter,
    WitnessKey,
) {
    let records = vec![
        record("q-b", 512, 8),
        record("q-a", 512, 9),
        record("q-c", 512, 10),
    ];
    let plan = prepare_batch_queries(&records, &scheme(), options()).unwrap();
    let key = WitnessKey::from_packed(201).unwrap();
    let block = BatchSequenceBlockRequest::new(2, 100, 140);
    let works = vec![
        candidate(0, "mg-a", key),
        candidate(1, "mg-a", key),
        candidate(2, "mg-b", key),
    ];
    let mut adapter = InMemoryBatchAdapter::new();
    adapter.insert_key_page("mg-a", 7, [key]);
    adapter.insert_key_page("mg-b", 7, [key]);
    adapter.insert_sequence_block("mg-a", block, vec![b'A'; 40]);
    adapter.insert_sequence_block("mg-b", block, vec![b'C'; 40]);
    for work in &works {
        adapter.insert_candidate_plan(work.query_index, &work.metagenome_id, read_plan(key, block));
    }
    (plan, works, adapter, key)
}

#[test]
fn grouped_execution_reuses_jma_open_pages_and_blocks_without_changing_results() {
    let (plan, works, mut adapter, _key) = fixture_plan();
    let grouped = execute_batch(
        &plan,
        &works,
        &mut adapter,
        BatchExecutionOptions {
            output_order: BatchOutputOrder::Deterministic,
        },
    )
    .unwrap();
    assert_eq!(grouped.records.len(), 3);
    assert_eq!(grouped.metrics.metagenome_open_count, 2);
    assert_eq!(grouped.metrics.exact_key_requests, 2);
    assert_eq!(grouped.metrics.key_page_requests, 2);
    assert_eq!(grouped.metrics.key_pages_decoded, 2);
    assert_eq!(grouped.metrics.sequence_block_requests, 2);
    assert_eq!(grouped.metrics.sequence_blocks_decoded, 2);
    assert_eq!(adapter.metrics().opened_metagenomes, vec!["mg-a", "mg-b"]);
    assert_eq!(adapter.metrics().key_page_reads.len(), 2);
    assert_eq!(adapter.metrics().sequence_block_reads.len(), 2);

    let grouped_records = grouped.records.clone();
    adapter.reset_metrics();
    let mut separate_records = Vec::<BatchResultRecord<FixtureTraceResult>>::new();
    for work in &works {
        let output = execute_batch(
            &plan,
            std::slice::from_ref(work),
            &mut adapter,
            BatchExecutionOptions {
                output_order: BatchOutputOrder::Deterministic,
            },
        )
        .unwrap();
        separate_records.extend(output.records);
    }
    separate_records.sort_by(|left, right| left.key.cmp(&right.key));
    assert_eq!(grouped_records, separate_records);
}

#[test]
fn grouped_work_is_stable_and_keeps_result_keys_distinct() {
    let key = WitnessKey::from_packed(303).unwrap();
    let works = vec![
        candidate(2, "mg-b", key),
        candidate(0, "mg-a", key),
        candidate(1, "mg-a", key),
    ];
    let grouped = group_candidate_work(&works);
    assert_eq!(
        grouped.keys().cloned().collect::<Vec<_>>(),
        vec!["mg-a", "mg-b"]
    );
    assert_eq!(
        grouped["mg-a"]
            .iter()
            .map(|work| work.query_index)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
}

#[test]
fn hundred_query_fixture_has_no_lost_records_or_unbounded_page_reads() {
    let records = (0_u32..100)
        .map(|index| record(&format!("q-{index:03}"), 256, 1000 + u64::from(index)))
        .collect::<Vec<_>>();
    let plan = prepare_batch_queries(&records, &scheme(), options()).unwrap();
    let key = WitnessKey::from_packed(401).unwrap();
    let block = BatchSequenceBlockRequest::new(1, 0, 32);
    let works = (0..100)
        .map(|index| candidate(index, "mg-single", key))
        .collect::<Vec<_>>();
    let mut adapter = InMemoryBatchAdapter::new();
    adapter.insert_key_page("mg-single", 7, [key]);
    adapter.insert_sequence_block("mg-single", block, vec![b'G'; 32]);
    for work in &works {
        adapter.insert_candidate_plan(work.query_index, "mg-single", read_plan(key, block));
    }

    let output = execute_batch(
        &plan,
        &works,
        &mut adapter,
        BatchExecutionOptions::default(),
    )
    .unwrap();
    assert_eq!(output.records.len(), 100);
    assert_eq!(output.metrics.metagenome_open_count, 1);
    assert_eq!(output.metrics.key_pages_decoded, 1);
    assert_eq!(output.metrics.sequence_blocks_decoded, 1);
    assert!(
        output
            .records
            .windows(2)
            .all(|pair| pair[0].key < pair[1].key)
    );
}

#[test]
fn routing_error_retains_query_identity() {
    let plan = prepare_batch_queries(&[record("q-fail", 256, 42)], &scheme(), options()).unwrap();
    struct FailingRouter;
    impl BatchRouter for FailingRouter {
        type Error = RecordingRouterError;
        fn route_query(
            &mut self,
            _query: &jam_rs::trace::contracts::PreparedBatchQuery,
            _witnesses: &[jam_rs::router::search::TieredQueryWitness],
        ) -> Result<Vec<RoutedCandidate>, Self::Error> {
            Err(RecordingRouterError)
        }
    }
    let error = route_batch(&plan, &mut FailingRouter).unwrap_err();
    assert_eq!(error.query_id, "q-fail");
}
