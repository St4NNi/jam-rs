//! The shard a server answers for: opened once, then searched by many requests at once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ServiceError;
use super::protocol::{Reply, SearchRequest};
use crate::index::Shard;
use crate::loss::MissModel;
use crate::report::{Record, ReportError};
use crate::search::{EntryEvent, SearchConfig, SearchControl, Searcher};
use std::sync::atomic::Ordering;
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, Mutex, PoisonError};
use tokio::sync::mpsc;

/// Entry events held between the search and their messages.
const EVENT_QUEUE: usize = 64;

/// Opens the shard again; `None` when it cannot be opened.
pub type Reopen = Box<dyn Fn() -> Option<Shard> + Send + Sync>;

/// An open shard with its name and miss models.
pub struct ServedShard {
    name: String,
    shard: Mutex<Arc<Shard>>,
    /// The miss models for seeds hit exactly and for seeds hit with one mismatch.
    models: [MissModel; 2],
    /// How to open the shard again once an S3 object of it was replaced.
    reopen: Option<Reopen>,
}

impl ServedShard {
    /// Serves `shard` under `name`, the name its records carry. Fails when the miss models
    /// (exact, then one mismatch) do not fit the shard's k and minimizer window.
    pub fn new(name: String, shard: Shard, models: [MissModel; 2]) -> Result<Self, ServiceError> {
        Searcher::new(&shard, &models[0], &models[1])?;
        Ok(Self {
            name,
            shard: Mutex::new(Arc::new(shard)),
            models,
            reopen: None,
        })
    }

    /// The same shard, opened again with `reopen` after an S3 object of it was replaced, so a
    /// shard uploaded again is served without a restart.
    pub fn with_reopen(mut self, reopen: Reopen) -> Self {
        self.reopen = Some(reopen);
        self
    }

    /// The shard that searches use now.
    fn current(&self) -> Arc<Shard> {
        let shard = self.shard.lock().unwrap_or_else(PoisonError::into_inner);
        Arc::clone(&shard)
    }

    /// Opens the shard again if `searched` is stale and no search did it yet. Searches that
    /// fail for the replaced object end with a retryable error; the next one finds the new shard.
    fn refresh(&self, searched: &Arc<Shard>) {
        let Some(reopen) = &self.reopen else {
            return;
        };
        if !searched.is_stale() {
            return;
        }
        let mut shard = self.shard.lock().unwrap_or_else(PoisonError::into_inner);
        if Arc::ptr_eq(&shard, searched)
            && let Some(opened) = reopen()
        {
            *shard = Arc::new(opened);
        }
    }

    /// The name the records of this shard carry.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Searches the shard and sends one message per entry event, then the summary; a search
    /// that cannot run sends a failed summary. It blocks, so it runs off the async threads. When
    /// `replies` is closed (the client left), it sets the stop flag and the search ends.
    pub fn search(
        &self,
        request: &SearchRequest,
        config: &SearchConfig,
        control: SearchControl<'_>,
        replies: &mpsc::Sender<Reply>,
    ) {
        let shard = self.current();
        let summary = match self.search_events(&shard, request, config, control, replies) {
            Ok(summary) => Reply::record(self.name.clone(), summary),
            Err(error) => Reply::failed(self.name.clone(), request, &error),
        };
        self.refresh(&shard);
        if replies.blocking_send(summary).is_err() {
            control.stop.store(true, Ordering::Relaxed);
        }
    }

    /// Runs the search on a scoped thread and turns its events into messages. Leaving the loop
    /// drops the event receiver, which also stops the search.
    fn search_events(
        &self,
        shard: &Shard,
        request: &SearchRequest,
        config: &SearchConfig,
        control: SearchControl<'_>,
        replies: &mpsc::Sender<Reply>,
    ) -> Result<Record, ServiceError> {
        // A searcher borrows the shard, so each search builds its own; that is cheap.
        let searcher = Searcher::new(shard, &self.models[0], &self.models[1])?;
        let (events, received) = sync_channel(EVENT_QUEUE);
        let query = &request.query_name;
        std::thread::scope(|scope| {
            let search = scope.spawn(move || {
                searcher.search(&request.bases, request.circular, config, control, &events)
            });
            for event in received {
                let reply = Reply::record(self.name.clone(), self.record(shard, query, &event)?);
                if replies.blocking_send(reply).is_err() {
                    control.stop.store(true, Ordering::Relaxed);
                    break;
                }
            }
            let summary = search.join().map_err(|_| ServiceError::Panic)??;
            Ok(Record::summary(query, &summary))
        })
    }

    fn record(
        &self,
        shard: &Shard,
        query: &str,
        event: &EntryEvent,
    ) -> Result<Record, ServiceError> {
        let index = event.entry();
        let entry = (shard.manifest().entries.get(index as usize))
            .ok_or(ReportError::UnknownEntry(index))?;
        Ok(Record::from_event(query, &entry.name, event)?)
    }
}
