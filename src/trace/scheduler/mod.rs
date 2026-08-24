//! Bounded candidate work scheduling.

mod queue;

pub use queue::{
    BoundedCandidateQueue, CandidateScheduler, CandidateWorkQueue, QueueError, QueueState,
};
