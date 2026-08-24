use jam_rs::trace::memory::{
    ScratchCapacities, ScratchPoolError, TraceWorkerScratch, TraceWorkerScratchPool,
};

#[test]
fn estimated_bytes_matches_every_scratch_buffer_element_type() {
    let capacities = ScratchCapacities {
        key_decode: 2,
        posting_decode: 3,
        anchors: 5,
        chains: 7,
        sequence: 11,
        score_rows: 13,
        traceback: 17,
        cigar: 19,
        edit_script: 23,
        output: 29,
    };
    // u8 fields: key/posting/sequence/cigar/edit/output; u64 fields:
    // anchors/chains; i32 fields: score rows/traceback.
    assert_eq!(capacities.estimated_bytes(), 303);
}

#[test]
fn reset_reuses_all_worker_buffer_allocations() {
    let capacities = ScratchCapacities {
        key_decode: 8,
        posting_decode: 16,
        anchors: 4,
        chains: 3,
        sequence: 32,
        score_rows: 8,
        traceback: 8,
        cigar: 16,
        edit_script: 16,
        output: 32,
    };
    let mut scratch = TraceWorkerScratch::with_capacities(capacities);
    scratch.key_decode().extend_from_slice(&[1, 2, 3]);
    scratch.posting_decode().extend_from_slice(&[4, 5]);
    scratch.anchors().extend_from_slice(&[6, 7]);
    scratch.chains().push(8);
    scratch.sequence().extend_from_slice(&[9, 10]);
    scratch.score_rows().push(11);
    scratch.traceback().push(12);
    scratch.cigar().push(13);
    scratch.edit_script().push(14);
    scratch.output().push(15);
    let before = scratch.capacities();
    scratch.reset();
    assert_eq!(scratch.capacities(), before);
    assert!(scratch.key_decode().is_empty());
    assert!(scratch.anchors().is_empty());
    scratch.key_decode().extend_from_slice(&[16, 17]);
    assert_eq!(scratch.key_decode().as_slice(), &[16, 17]);
}

#[test]
fn scratch_pool_has_one_reusable_buffer_per_worker_and_closes_cleanly() {
    let pool = TraceWorkerScratchPool::new(2, ScratchCapacities::default())
        .expect("positive worker count");
    let first = pool.acquire().expect("first worker scratch");
    let second = pool.acquire().expect("second worker scratch");
    assert_eq!(pool.available(), 0);
    drop(first);
    assert_eq!(pool.available(), 1);
    drop(second);
    assert_eq!(pool.available(), 2);
    pool.close();
    assert!(matches!(pool.try_acquire(), Err(ScratchPoolError::Closed)));
}

#[test]
fn scratch_pool_rejects_zero_workers() {
    assert!(matches!(
        TraceWorkerScratchPool::new(0, ScratchCapacities::default()),
        Err(ScratchPoolError::InvalidWorkerCount)
    ));
}
