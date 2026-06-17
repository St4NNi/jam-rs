//! The shard a server answers for: opened once at start, then searched by many requests at once.

use super::ServiceError;
use super::frame::{Frame, ShardIdentity};
use super::records::{EventRecord, SummaryRecord};
use crate::index::Shard;
use crate::index::manifest::{MANIFEST_FILE, sha256_file};
use crate::loss::MissModel;
use crate::search::{SearchConfig, SearchControl, Searcher};
use std::path::PathBuf;
use std::time::{Duration, Instant};

/// Where a server finds its shard.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShardOptions {
    /// The shard directory.
    pub directory: PathBuf,
    /// The directory of the packed sequences, when not `directory/packed` (joined shards may
    /// share one).
    pub packed_directory: Option<PathBuf>,
}

/// An open shard with its miss model and identity. The sketch stays in RAM while it is open.
pub struct ServedShard {
    shard: Shard,
    model: MissModel,
    identity: ShardIdentity,
    opened_in: Duration,
}

impl ServedShard {
    /// Opens the shard of `options`. The miss model must fit the shard's k and minimizer window;
    /// it is built once by the caller, as it takes a while.
    pub fn open(options: &ShardOptions, model: MissModel) -> Result<Self, ServiceError> {
        let started = Instant::now();
        let directory = &options.directory;
        let shard = match &options.packed_directory {
            Some(packed) => Shard::open_with_packed(directory, packed)?,
            None => Shard::open(directory)?,
        };
        // Fails early when the model does not fit the shard.
        Searcher::new(&shard, &model)?;
        let name = directory
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| directory.display().to_string());
        let (_, digest) = sha256_file(&directory.join(MANIFEST_FILE))?;
        let identity = ShardIdentity::of(&shard, name, digest);
        Ok(Self {
            shard,
            model,
            identity,
            opened_in: started.elapsed(),
        })
    }

    /// What the server tells clients about this shard.
    pub fn identity(&self) -> &ShardIdentity {
        &self.identity
    }

    /// The open shard.
    pub fn shard(&self) -> &Shard {
        &self.shard
    }

    /// How long opening took, so a slow start can be told apart from a hung one.
    pub fn opened_in(&self) -> Duration {
        self.opened_in
    }

    /// Searches the shard and sends its frames, tagged with `tag`: one per entry event, then the
    /// summary, or `ShardFailed` when the search cannot start.
    pub fn search(
        &self,
        query: &[u8],
        circular: bool,
        config: &SearchConfig,
        control: SearchControl<'_>,
        tag: u32,
        send: &(dyn Fn(Frame) + Sync),
    ) {
        let on_event = |event| {
            send(Frame::Event {
                shard: tag,
                record: EventRecord::new(&self.shard, event),
            })
        };
        let outcome = Searcher::new(&self.shard, &self.model)
            .and_then(|searcher| searcher.search(query, circular, config, control, &on_event));
        send(match outcome {
            Ok(summary) => Frame::Summary {
                shard: tag,
                summary: SummaryRecord::new(&summary),
            },
            Err(error) => Frame::ShardFailed {
                shard: tag,
                message: error.to_string(),
            },
        });
    }
}
