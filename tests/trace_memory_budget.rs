use jam_rs::trace::memory::{HeapRank, TraceMemorySemaphore, retain_top_k};
use jam_rs::trace::{MemoryBudgetError, TraceMemoryEstimate};
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

fn estimate(bytes: u64) -> TraceMemoryEstimate {
    TraceMemoryEstimate {
        output_bytes: bytes,
        ..TraceMemoryEstimate::default()
    }
}

#[test]
fn reservation_releases_after_drop_and_rejects_one_oversized_candidate() {
    let budget = TraceMemorySemaphore::new(10);
    let reservation = budget.acquire(estimate(7)).expect("7 bytes fit");
    assert_eq!(reservation.bytes(), 7);
    assert_eq!(budget.in_use_bytes(), 7);
    assert!(matches!(
        budget.acquire(estimate(11)),
        Err(MemoryBudgetError::CandidateTooLarge {
            requested_bytes: 11,
            budget_bytes: 10,
        })
    ));
    drop(reservation);
    assert_eq!(budget.in_use_bytes(), 0);
    assert_eq!(budget.available_bytes(), 10);
}

#[test]
fn fifo_admission_prevents_small_candidates_from_bypassing_a_large_waiter() {
    let budget = Arc::new(TraceMemorySemaphore::new(3));
    let held = budget.acquire(estimate(3)).expect("initial reservation");
    let first_budget = Arc::clone(&budget);
    let (first_started_tx, first_started_rx) = mpsc::channel();
    let first = thread::spawn(move || {
        let reservation = first_budget.acquire(estimate(3)).expect("first waiter");
        first_started_tx.send(()).expect("first progress signal");
        thread::sleep(Duration::from_millis(20));
        drop(reservation);
    });
    while budget.waiting_candidates() != 1 {
        thread::yield_now();
    }

    let second_budget = Arc::clone(&budget);
    let (second_started_tx, second_started_rx) = mpsc::channel();
    let second = thread::spawn(move || {
        let reservation = second_budget.acquire(estimate(1)).expect("second waiter");
        second_started_tx.send(()).expect("second progress signal");
        drop(reservation);
    });
    while budget.waiting_candidates() != 2 {
        thread::yield_now();
    }

    drop(held);
    first_started_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("large waiter must make progress");
    assert!(second_started_rx.try_recv().is_err());
    first.join().expect("first waiter thread");
    second_started_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("small waiter must make progress after first release");
    second.join().expect("second waiter thread");
    assert_eq!(budget.in_use_bytes(), 0);
}

#[test]
fn close_releases_blocked_waiters_without_leaking_existing_reservations() {
    let budget = Arc::new(TraceMemorySemaphore::new(1));
    let held = budget.acquire(estimate(1)).expect("initial reservation");
    let waiting_budget = Arc::clone(&budget);
    let waiter = thread::spawn(move || waiting_budget.acquire(estimate(1)));
    while budget.waiting_candidates() != 1 {
        thread::yield_now();
    }
    budget.close();
    assert!(matches!(
        waiter.join().expect("waiter thread"),
        Err(MemoryBudgetError::Closed)
    ));
    assert!(budget.is_closed());
    assert_eq!(budget.in_use_bytes(), 1);
    drop(held);
    assert_eq!(budget.in_use_bytes(), 0);
    assert!(matches!(
        budget.acquire(estimate(1)),
        Err(MemoryBudgetError::Closed)
    ));
}

#[test]
fn with_candidate_releases_only_after_serialization_closure_returns() {
    let budget = TraceMemorySemaphore::new(32);
    let result = budget
        .with_candidate_result(estimate(12), |reservation| {
            assert_eq!(reservation.bytes(), 12);
            assert_eq!(budget.in_use_bytes(), 12);
            42_u32
        })
        .expect("candidate closure");
    assert_eq!(result, 42);
    assert_eq!(budget.in_use_bytes(), 0);
}

#[test]
fn top_k_streaming_helper_never_retains_more_than_the_cap() {
    let mut rank_calls = 0;
    let values = retain_top_k(0_u64..10_000, 3, |value| {
        rank_calls += 1;
        HeapRank::new(*value as i64, *value)
    });
    assert_eq!(rank_calls, 10_000);
    assert_eq!(values, vec![9_999, 9_998, 9_997]);
}
