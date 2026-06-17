//! Small generated shards and running servers for the service tests.

use super::frame::{Frame, SearchRequest};
use super::local::{ServedShard, ShardOptions};
use super::records::{EventRecord, SummaryRecord};
use super::server::{Server, ServerConfig};
use super::settings::SearchSettings;
use crate::decide::QueryCoverage;
use crate::index::{BuildOptions, EntryInput, Parameters, build_shard};
use crate::loss::MissModel;
use crate::search::{SearchConfig, SearchControl};
use std::net::TcpListener;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, OnceLock};
use std::thread::JoinHandle;

/// The query length of the tests.
pub const QUERY_LENGTH: usize = 5_000;

/// Random bases from a xorshift generator seeded with `seed`.
pub fn bases(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed.max(1);
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state % 4) as usize]
        })
        .collect()
}

/// The query of the tests.
pub fn query() -> Vec<u8> {
    bases(0x5eed, QUERY_LENGTH)
}

/// Builds shard `name` in `directory` with three entries: the query planted in random bases with
/// one substitution every 500 bases, and two unrelated entries.
pub fn build(directory: &Path, name: &str, seed: u64) -> PathBuf {
    let mut planted = bases(seed, 20_000);
    let mut copy = query();
    for base in copy.iter_mut().step_by(500) {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    planted.splice(6_000..6_000 + QUERY_LENGTH, copy);
    let sequences = [planted, bases(seed + 1, 15_000), bases(seed + 2, 12_000)];
    let inputs: Vec<EntryInput> = sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| {
            let path = directory.join(format!("{name}_entry{index}.fa"));
            let mut text = format!(">{name}_contig{index}\n").into_bytes();
            text.extend_from_slice(sequence);
            std::fs::write(&path, text).unwrap();
            EntryInput {
                name: format!("{name}_entry{index}"),
                path,
            }
        })
        .collect();
    let output = directory.join(name);
    let options = BuildOptions {
        parameters: Parameters {
            k: 21,
            window: 20,
            bin_bits: 8,
            ..Parameters::default()
        },
        bucket_bits: None,
        threads: 2,
        temp_dir: directory.to_path_buf(),
        entry_per_record: false,
    };
    build_shard(&inputs, &output, &options).unwrap();
    output
}

/// Opens the shard at `directory` with the test miss model.
pub fn open(directory: &Path) -> ServedShard {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    let model = MODEL.get_or_init(|| MissModel::with_layouts(21, 20, 200).unwrap());
    let options = ShardOptions {
        directory: directory.to_path_buf(),
        packed_directory: None,
    };
    ServedShard::open(&options, model.clone()).unwrap()
}

/// The search settings of the tests: 90% query coverage at 95% identity.
pub fn config() -> SearchConfig {
    let mut config = SearchConfig {
        threads: 2,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config.rule.min_query_coverage = QueryCoverage::Share {
        numerator: 9,
        decimal_places: 1,
    };
    config
}

/// A search for the test query.
pub fn request(deadline_ms: Option<u64>) -> SearchRequest {
    SearchRequest {
        query_name: "query".into(),
        query: String::from_utf8(query()).unwrap(),
        circular: false,
        settings: SearchSettings::of(&config()),
        deadline_ms,
    }
}

/// The frames of a search of `shard` in this process, tagged with `tag`.
pub fn local_frames(shard: &ServedShard, tag: u32) -> Vec<Frame> {
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let (sender, receiver) = sync_channel(16);
    std::thread::scope(|scope| {
        scope.spawn(move || shard.search(&query(), false, &config(), control, tag, &sender));
        receiver.into_iter().collect()
    })
}

/// The event records of `shard` as sorted JSON text, as the order between entries may change.
pub fn events_of(frames: &[Frame], shard: u32) -> Vec<String> {
    let mut events: Vec<String> = frames
        .iter()
        .filter_map(|frame| match frame {
            Frame::Event { shard: tag, record } if *tag == shard => Some(record),
            _ => None,
        })
        .map(|record: &EventRecord| serde_json::to_string(record).unwrap())
        .collect();
    events.sort();
    events
}

/// The summary of `shard`, with its times cleared.
pub fn summary_of(frames: &[Frame], shard: u32) -> SummaryRecord {
    let mut summary = frames
        .iter()
        .find_map(|frame| match frame {
            Frame::Summary {
                shard: tag,
                summary,
            } if *tag == shard => Some(summary.clone()),
            _ => None,
        })
        .unwrap();
    summary.sketch_seconds = 0.0;
    summary.locate_seconds = 0.0;
    summary.align_seconds = 0.0;
    summary
}

/// A server that serves on a localhost port until dropped.
pub struct Running {
    pub address: String,
    pub server: Arc<Server>,
    stop: Arc<AtomicBool>,
    thread: Option<JoinHandle<()>>,
}

impl Running {
    /// Serves `server` on a free localhost port.
    pub fn start(server: Server) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let address = listener.local_addr().unwrap().to_string();
        let (server, stop) = (Arc::new(server), Arc::new(AtomicBool::new(false)));
        let (serving, stopping) = (server.clone(), stop.clone());
        let thread = std::thread::spawn(move || serving.serve(listener, &stopping).unwrap());
        Self {
            address,
            server,
            stop,
            thread: Some(thread),
        }
    }

    /// Serves `shard` alone with the given configuration.
    pub fn alone(shard: ServedShard, config: ServerConfig) -> Self {
        Self::start(Server::new(shard, config))
    }
}

impl Drop for Running {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}
