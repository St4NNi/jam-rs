use jam_rs::resource::ByteRange;
use jam_rs::resource::read_plan::{
    PhysicalReadObservation, RangePlanLimits, RangePlanner, RangeRequest, ReadPlanError,
};

fn range(offset: u64, length: u64) -> ByteRange {
    ByteRange::new(offset, length).unwrap()
}

#[test]
fn nearby_ranges_coalesce_with_stable_logical_mapping_and_accounting() {
    let planner = RangePlanner::new(RangePlanLimits {
        max_gap_bytes: 5,
        max_overread_bytes: 5,
        max_physical_range_bytes: 64,
        max_concurrency: 2,
    })
    .unwrap();
    let plan = planner
        .plan([
            RangeRequest::new(30, range(50, 5)),
            RangeRequest::new(10, range(10, 10)),
            RangeRequest::new(20, range(25, 10)),
        ])
        .unwrap();

    assert_eq!(plan.physical_ranges().len(), 2);
    assert_eq!(plan.physical_ranges()[0].range, range(10, 25));
    assert_eq!(plan.physical_ranges()[0].logical_request_ids, vec![10, 20]);
    assert_eq!(plan.physical_ranges()[0].useful_logical_bytes, 20);
    assert_eq!(plan.physical_ranges()[0].overread_bytes, 5);
    assert_eq!(plan.mappings()[0].request_id, 30);
    assert_eq!(plan.mappings()[0].physical_index, Some(1));
    assert_eq!(plan.mappings()[1].request_id, 10);
    assert_eq!(plan.mappings()[1].physical_offset, Some(0));
    assert_eq!(plan.mappings()[2].request_id, 20);
    assert_eq!(plan.mappings()[2].physical_offset, Some(15));

    let counters = plan.counters();
    assert_eq!(counters.logical_bytes, 25);
    assert_eq!(counters.physical_requested_bytes, 30);
    assert_eq!(counters.physical_returned_bytes, 30);
    assert_eq!(counters.useful_decoded_bytes, 25);
    assert_eq!(counters.overread_bytes, 5);
    assert_eq!(counters.physical_requests, 2);
}

#[test]
fn cache_hits_are_mapped_without_physical_requests() {
    let planner = RangePlanner::new(RangePlanLimits {
        max_gap_bytes: 16,
        max_overread_bytes: 16,
        max_physical_range_bytes: 64,
        max_concurrency: 4,
    })
    .unwrap();
    let plan = planner
        .plan([
            RangeRequest::cached(1, range(0, 4)),
            RangeRequest::new(2, range(8, 4)),
        ])
        .unwrap();
    assert_eq!(plan.physical_ranges().len(), 1);
    assert_eq!(plan.mappings()[0].physical_index, None);
    assert!(plan.mappings()[0].cache_hit);
    assert_eq!(plan.mappings()[1].physical_index, Some(0));
    let counters = plan.counters();
    assert_eq!(counters.cache_hits, 1);
    assert_eq!(counters.cache_misses, 1);
    assert_eq!(counters.cache_bytes, 4);
    assert_eq!(counters.physical_returned_bytes, 4);
    assert_eq!(counters.useful_decoded_bytes, 8);
}

#[test]
fn overread_and_physical_caps_stop_coalescing() {
    let planner = RangePlanner::new(RangePlanLimits {
        max_gap_bytes: 10,
        max_overread_bytes: 2,
        max_physical_range_bytes: 12,
        max_concurrency: 1,
    })
    .unwrap();
    let plan = planner
        .plan([
            RangeRequest::new(1, range(0, 4)),
            RangeRequest::new(2, range(7, 4)),
            RangeRequest::new(3, range(20, 4)),
        ])
        .unwrap();
    assert_eq!(plan.physical_ranges().len(), 3);
    assert_eq!(plan.batches().len(), 3);
    assert_eq!(plan.batches()[1].physical_indices, vec![1]);
}

#[test]
fn planner_rejects_overflow_and_ranges_above_cap() {
    let planner = RangePlanner::new(RangePlanLimits {
        max_gap_bytes: 1,
        max_overread_bytes: 1,
        max_physical_range_bytes: 8,
        max_concurrency: 1,
    })
    .unwrap();
    assert_eq!(
        planner
            .plan([RangeRequest::new(
                9,
                ByteRange {
                    offset: u64::MAX,
                    length: 1,
                },
            )])
            .unwrap_err(),
        ReadPlanError::RangeOverflow { request_id: 9 }
    );
    assert_eq!(
        planner
            .plan([RangeRequest::new(8, range(0, 9))])
            .unwrap_err(),
        ReadPlanError::RangeExceedsPhysicalCap { request_id: 8 }
    );
}

#[test]
fn physical_observations_replace_estimates_without_losing_cache_accounting() {
    let planner = RangePlanner::new(RangePlanLimits {
        max_gap_bytes: 0,
        max_overread_bytes: 0,
        max_physical_range_bytes: 32,
        max_concurrency: 2,
    })
    .unwrap();
    let plan = planner
        .plan([
            RangeRequest::cached(1, range(0, 3)),
            RangeRequest::new(2, range(10, 5)),
        ])
        .unwrap();
    let counters = plan
        .accounting_with_returns([PhysicalReadObservation {
            physical_index: 0,
            returned_bytes: 7,
            useful_decoded_bytes: 5,
        }])
        .unwrap();
    assert_eq!(counters.cache_bytes, 3);
    assert_eq!(counters.physical_returned_bytes, 7);
    assert_eq!(counters.useful_decoded_bytes, 8);
    assert_eq!(counters.overread_bytes, 2);
}
