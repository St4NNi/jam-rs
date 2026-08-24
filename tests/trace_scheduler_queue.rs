use jam_rs::trace::scheduler::{CandidateScheduler, CandidateWorkQueue, QueueError};
use std::sync::Arc;
use std::sync::mpsc;
use std::thread;
use std::time::Duration;

#[test]
fn queue_preserves_fifo_order_and_hard_capacity() {
    let queue = CandidateWorkQueue::new(2).expect("positive capacity");
    assert_eq!(queue.capacity(), 2);
    queue.push(1_u32).expect("first item");
    queue.push(2_u32).expect("second item");
    assert_eq!(queue.len(), 2);
    assert_eq!(queue.try_push(3), Err(QueueError::Full));
    assert_eq!(queue.pop().expect("first pop"), Some(1));
    assert_eq!(queue.pop().expect("second pop"), Some(2));
    assert_eq!(queue.try_pop(), None);
}

#[test]
fn blocking_producer_makes_progress_after_worker_pops_a_slot() {
    let queue = Arc::new(CandidateWorkQueue::new(1).expect("positive capacity"));
    queue.push(1_u32).expect("initial item");
    let producer_queue = Arc::clone(&queue);
    let (submitted_tx, submitted_rx) = mpsc::channel();
    let producer = thread::spawn(move || {
        producer_queue.push(2_u32).expect("blocked producer");
        submitted_tx.send(()).expect("producer progress signal");
    });
    thread::sleep(Duration::from_millis(10));
    assert!(submitted_rx.try_recv().is_err());
    assert_eq!(queue.pop().expect("free one slot"), Some(1));
    submitted_rx
        .recv_timeout(Duration::from_secs(1))
        .expect("producer must make progress");
    producer.join().expect("producer thread");
    assert_eq!(queue.pop().expect("second item"), Some(2));
}

#[test]
fn close_wakes_blocked_producers_and_drained_consumers() {
    let queue = Arc::new(CandidateWorkQueue::new(1).expect("positive capacity"));
    queue.push(1_u32).expect("initial item");
    let producer_queue = Arc::clone(&queue);
    let producer = thread::spawn(move || producer_queue.push(2_u32));
    thread::sleep(Duration::from_millis(10));
    queue.close();
    assert_eq!(
        producer.join().expect("producer thread"),
        Err(QueueError::Closed)
    );
    assert_eq!(queue.pop().expect("drain initial item"), Some(1));
    assert_eq!(queue.pop().expect("closed and empty"), None);
    assert_eq!(queue.try_push(3), Err(QueueError::Closed));
}

#[test]
fn multiple_workers_drain_each_candidate_once_with_bounded_memory() {
    let queue = Arc::new(CandidateWorkQueue::new(3).expect("positive capacity"));
    let mut workers = Vec::new();
    for _ in 0..4 {
        let worker_queue = Arc::clone(&queue);
        workers.push(thread::spawn(move || {
            let mut items = Vec::new();
            while let Some(item) = worker_queue.pop().expect("queue pop") {
                items.push(item);
            }
            items
        }));
    }
    for item in 0_u32..100 {
        queue.push(item).expect("bounded producer");
        assert!(queue.len() <= queue.capacity());
    }
    queue.close();
    let mut all = workers
        .into_iter()
        .flat_map(|worker| worker.join().expect("worker thread"))
        .collect::<Vec<_>>();
    all.sort_unstable();
    assert_eq!(all, (0_u32..100).collect::<Vec<_>>());
}

#[test]
fn scheduler_facade_closes_and_drains_its_queue() {
    let scheduler = CandidateScheduler::new(2).expect("positive capacity");
    scheduler.submit("a").expect("first candidate");
    scheduler.submit("b").expect("second candidate");
    scheduler.finish();
    assert_eq!(scheduler.next().expect("first candidate"), Some("a"));
    assert_eq!(scheduler.next().expect("second candidate"), Some("b"));
    assert_eq!(scheduler.next().expect("closed queue"), None);
}

#[test]
fn zero_capacity_is_rejected_before_work_can_be_queued() {
    assert!(matches!(
        CandidateWorkQueue::<u32>::new(0),
        Err(QueueError::InvalidCapacity)
    ));
}
