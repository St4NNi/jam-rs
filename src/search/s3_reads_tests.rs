//! Searches of shards on a slow test store: the results equal a local search, cold and warm, and
//! the store counts each request by the part of the shard it reads and how many run at once.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::control::SearchControl;
use super::driver::{SearchConfig, SearchSummary, Searcher};
use super::events::EntryEvent;
use crate::index::{Location, Part, Plan, Shard, ShardOptions, build_shard, warm};
use crate::storage::{BlockCache, CacheConfig, S3Config, S3Limits};
use crate::test_support::{
    Reply, build_options, diverged, etag_of, no_stop, random_bases, rng, runtime, serve,
    write_entries,
};
use rand::RngExt;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

const SHARDS: usize = 2;
const ENTRIES: usize = 40;
const QUERY_LENGTH: usize = 1_300;
/// The wait of the store per request, so requests overlap as on S3.
const DELAY: Duration = Duration::from_millis(3);

/// The bytes and ETag of a file of the test store.
type StoredFile = Arc<(Vec<u8>, String)>;

/// How fast the test store answers: a wait per request, then the bytes at the rate of one
/// stream, and all answers share one link of `link_rate` bytes per second (0: no rate).
#[derive(Clone, Copy)]
struct Speed {
    latency: Duration,
    stream_rate: u64,
    link_rate: u64,
}

/// A store that serves the files below a directory at a `Speed`, and logs each request as
/// (key, first byte, end), with the most requests it answered at once.
struct Store {
    speed: Speed,
    log: Mutex<Vec<(String, u64, u64)>>,
    active: AtomicUsize,
    peak: AtomicUsize,
    /// When the shared link is free again.
    link_free: Mutex<Instant>,
    /// The bytes and ETag of each file read so far.
    files: Mutex<HashMap<String, StoredFile>>,
}

impl Store {
    /// Forgets the requests so far.
    fn reset(&self) {
        self.log.lock().unwrap().clear();
        self.peak.store(0, Ordering::SeqCst);
    }

    /// The bytes and ETag of the file `key` below `root`, read once.
    fn file(&self, root: &Path, key: &str) -> Option<StoredFile> {
        let mut files = self.files.lock().unwrap();
        if let Some(file) = files.get(key) {
            return Some(Arc::clone(file));
        }
        let data = std::fs::read(root.join(key)).ok()?;
        let etag = etag_of(&data);
        let file = Arc::new((data, etag));
        files.insert(key.to_string(), Arc::clone(&file));
        Some(file)
    }

    /// Waits as long as an answer of `bytes` bytes takes at the store's speed.
    fn wait(&self, bytes: u64) {
        let Speed {
            latency,
            stream_rate,
            link_rate,
        } = self.speed;
        let time = |rate: u64| match rate {
            0 => Duration::ZERO,
            rate => Duration::from_secs_f64(bytes as f64 / rate as f64),
        };
        let now = Instant::now();
        let link_done = {
            let mut free = self.link_free.lock().unwrap();
            *free = (*free).max(now) + time(link_rate);
            *free
        };
        let done = link_done.max(now + latency + time(stream_rate));
        std::thread::sleep(done.saturating_duration_since(Instant::now()));
    }
}

fn serve_slow(root: &Path) -> (crate::test_support::Server, Arc<Store>) {
    let speed = Speed {
        latency: DELAY,
        stream_rate: 0,
        link_rate: 0,
    };
    serve_at(root, speed)
}

fn serve_at(root: &Path, speed: Speed) -> (crate::test_support::Server, Arc<Store>) {
    let store = Arc::new(Store {
        speed,
        log: Mutex::default(),
        active: AtomicUsize::new(0),
        peak: AtomicUsize::new(0),
        link_free: Mutex::new(Instant::now()),
        files: Mutex::default(),
    });
    let (root, seen) = (root.to_path_buf(), Arc::clone(&store));
    let server = serve(move |_, request| {
        let now = seen.active.fetch_add(1, Ordering::SeqCst) + 1;
        seen.peak.fetch_max(now, Ordering::SeqCst);
        let key = request.path.strip_prefix("/bucket/").unwrap_or_default();
        let reply = match seen.file(&root, key) {
            Some(file) => {
                let (data, etag) = (&file.0, &file.1);
                let (start, last) = request.range().unwrap_or((0, data.len().saturating_sub(1)));
                seen.wait((last + 1 - start) as u64);
                let range = (key.to_string(), start as u64, last as u64 + 1);
                seen.log.lock().unwrap().push(range);
                Reply::object_tagged(data, request, etag.clone())
            }
            None => {
                seen.wait(0);
                Reply::status(404)
            }
        };
        seen.active.fetch_sub(1, Ordering::SeqCst);
        reply
    });
    (server, store)
}

/// `shards` shards of 40 entries of three contigs of 20 to 100 kb, for a query of `length` bases. The
/// query lies with 1% changes in entries 0 and 1 of each shard, a long one cut into a piece per
/// contig, as a plasmid in a metagenome; stretches of 60 to 400 of its bases lie in most other
/// entries, so they become candidates with few seed hits, as most candidates of a large index
/// are.
fn build_shards(directory: &Path, length: usize, shards: usize) -> Vec<u8> {
    let mut random = rng(0x53_ead);
    let query = random_bases(&mut random, length);
    for shard in 0..shards {
        let mut entries = Vec::new();
        for entry in 0..ENTRIES {
            let mut contigs: Vec<Vec<u8>> = (0..3)
                .map(|_| {
                    let length = random.random_range(20_000..100_000);
                    random_bases(&mut random, length)
                })
                .collect();
            let planted = match entry {
                0 | 1 => {
                    let copy = diverged(&mut random, &query, 0.99);
                    let pieces = if length > 10_000 { 3 } else { 1 };
                    copy.chunks(length.div_ceil(pieces))
                        .map(<[u8]>::to_vec)
                        .collect()
                }
                _ if entry % 8 == 7 => Vec::new(),
                _ => {
                    let piece = random.random_range(60..400);
                    let start = random.random_range(0..length - piece);
                    vec![query[start..start + piece].to_vec()]
                }
            };
            for (at, piece) in planted.into_iter().enumerate() {
                let contig = &mut contigs[(entry + at) % 3];
                let at = random.random_range(0..contig.len());
                contig.splice(at..at, piece);
            }
            entries.push(contigs);
        }
        let inputs = directory.join(format!("input{shard}"));
        std::fs::create_dir_all(&inputs).unwrap();
        let inputs = write_entries(&inputs, &entries);
        let output = directory.join(format!("shard{shard}.jam"));
        build_shard(&inputs, &output, &build_options(directory, 2)).unwrap();
    }
    query
}

fn config(threads: usize) -> SearchConfig {
    let mut config = SearchConfig {
        threads,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config
}

/// The summary of a search without its times, and the events of each entry in order.
fn search_shard(
    shard: &Shard,
    query: &[u8],
    threads: usize,
) -> (String, BTreeMap<u32, Vec<EntryEvent>>) {
    let searcher = Searcher::new(shard).unwrap();
    let (sender, receiver) = sync_channel(1 << 16);
    let control: SearchControl<'_> = no_stop();
    let summary = searcher.search(query, false, &config(threads), control, &sender);
    drop(sender);
    let mut events = BTreeMap::<u32, Vec<EntryEvent>>::new();
    for event in receiver {
        events.entry(event.entry()).or_default().push(event);
    }
    let mut summary: SearchSummary = summary.unwrap();
    summary.sketch.seconds = 0.0;
    summary.locate.seconds = 0.0;
    summary.align.seconds = 0.0;
    summary.times = Default::default();
    (format!("{summary:?}"), events)
}

/// The results of searching every shard in turn, as one text to compare.
fn search_all(shards: &[Shard], query: &[u8], threads: usize) -> String {
    let results: Vec<_> = shards
        .iter()
        .map(|shard| search_shard(shard, query, threads))
        .collect();
    format!("{results:?}")
}

/// The part of a shard that a request reads; `local` are the shards opened locally.
fn part_of(local: &[Shard], key: &str, start: u64, end: u64) -> &'static str {
    if end - start == 1 {
        return "check";
    }
    let number: usize = key
        .trim_start_matches("shard")
        .split('.')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    let manifest = local[number].manifest();
    if key.ends_with(".jam") {
        let tables = local[number].seed_file().contig_tables();
        let locators = manifest.locator_part.clone().unwrap_or_default();
        return match start {
            0 if end <= 256 => "seeds header",
            _ if start < manifest.sketch_part.end => "sketch",
            _ if tables.contains(&start) => "contig tables",
            _ if start < locators.end => "locators",
            _ => "manifest",
        };
    }
    let sequence = (manifest.entries.iter()).find(|entry| {
        entry
            .sequence
            .as_ref()
            .is_some_and(|s| s.range().contains(&start))
    });
    let Some(entry) = sequence else {
        return "packed gap";
    };
    let offset = entry.sequence.as_ref().unwrap().offset;
    let data = entry.bases.div_ceil(4);
    let names_at = offset + 4096 + data + 4 * data.div_ceil(4092);
    let header = start < offset + 4096;
    let bases = start < names_at && end > offset + 4096;
    match (header, bases, end > names_at) {
        (true, false, false) => "packed header",
        (false, true, false) => "packed bases",
        (false, false, true) => "packed names",
        (false, true, true) => "packed bases and names",
        _ => "packed header and more",
    }
}

/// Requests by part, and the most at once.
fn counts(local: &[Shard], store: &Store) -> (BTreeMap<&'static str, usize>, usize) {
    let mut counts = BTreeMap::new();
    for (key, start, end) in store.log.lock().unwrap().iter() {
        *counts.entry(part_of(local, key, *start, *end)).or_default() += 1;
    }
    (counts, store.peak.load(Ordering::SeqCst))
}

/// Opens the shards on the store, through a new empty block cache, with the pins in `pins`.
fn open_s3(config: &S3Config, pins: Option<&Path>, shards: usize) -> Vec<Shard> {
    let cache = CacheConfig {
        read_ahead: S3Limits::default().max_in_flight,
        ..CacheConfig::new(64 << 20, 64 << 10)
    };
    let options = ShardOptions {
        cache: Some(BlockCache::open(&cache, runtime()).unwrap()),
        pins: pins.map(Path::to_path_buf),
    };
    (0..shards)
        .map(|shard| {
            let location = Location::s3(&format!("s3://bucket/shard{shard}.jam"), config).unwrap();
            Shard::open_with(&location, &options).unwrap()
        })
        .collect()
}

/// Warms `cache` with `parts` of every shard, in RAM up to `memory` bytes, else on disk.
fn warm_all(config: &S3Config, cache: &Path, parts: &[Part], memory: u64, shards: usize) {
    for shard in 0..shards {
        let location = Location::s3(&format!("s3://bucket/shard{shard}.jam"), config).unwrap();
        let head = location.head().unwrap();
        let plan = Plan::choose(&head, memory, 1 << 30, parts).unwrap();
        warm(&location, &head, &plan, cache).unwrap();
    }
}

/// A search of two shards on the store gives the results of a local search, from an empty block
/// cache and with pinned parts. Prints the requests by part for each way.
#[test]
fn store_search_same() {
    store_search(QUERY_LENGTH);
}

/// The same for a query as long as a plasmid.
#[test]
fn store_plasmid_same() {
    store_search(30_000);
}

fn store_search(length: usize) {
    let directory = tempfile::tempdir().unwrap();
    let query = build_shards(directory.path(), length, SHARDS);
    let local: Vec<Shard> = (0..SHARDS)
        .map(|shard| Shard::open(&directory.path().join(format!("shard{shard}.jam"))).unwrap())
        .collect();
    let expected = search_all(&local, &query, 2);
    for shard in &local {
        let (summary, _) = search_shard(shard, &query, 2);
        let align = summary.split("align: ").nth(1).unwrap_or_default();
        eprintln!("{length} bases: align {}", &align[..align.len().min(400)]);
    }
    let (server, store) = serve_slow(directory.path());
    let s3 = S3Config::new(
        "local",
        Some(&server.address),
        true,
        S3Limits::default(),
        runtime(),
    );
    let s3 = s3.unwrap().unsigned();
    let cold = |pins: Option<&Path>, name: &str| {
        store.reset();
        let shards = open_s3(&s3, pins, SHARDS);
        let (opened, _) = counts(&local, &store);
        store.reset();
        assert_eq!(search_all(&shards, &query, 2), expected, "{name}");
        let (searched, peak) = counts(&local, &store);
        eprintln!("{name}: open {opened:?}; search {searched:?}, at most {peak} at once");
        // A second search of the open shards, warm, gives the same results.
        assert_eq!(search_all(&shards, &query, 2), expected, "{name}, warm");
        // With many candidates, the entries ahead keep many more requests in flight than the two
        // workers. A long query has only its two hits per shard.
        assert!(peak > 8 || length > 10_000, "{name}: {peak}");
        searched
    };
    let searched = cold(None, "cold");
    assert!(searched.keys().any(|part| part.contains("bases")));
    let cache = directory.path().join("cache");
    let in_ram = directory.path().join("cache-ram");
    warm_all(
        &s3,
        &in_ram,
        &[Part::Sketch, Part::Locators],
        1 << 30,
        SHARDS,
    );
    cold(Some(&in_ram), "sketch and locators pinned in RAM");
    warm_all(&s3, &cache, &[Part::Sketch, Part::Locators], 0, SHARDS);
    let searched = cold(Some(&cache), "sketch and locators pinned");
    // The seeds file is read from the pins alone.
    for part in [
        "seeds header",
        "sketch",
        "locators",
        "contig tables",
        "manifest",
    ] {
        assert!(!searched.contains_key(part), "{part}: {searched:?}");
    }
    let parts = [Part::Sketch, Part::Locators, Part::Headers];
    warm_all(&s3, &cache, &parts, 0, SHARDS);
    let searched = cold(Some(&cache), "sketch, locators and headers pinned");
    assert!(
        !searched.keys().any(|part| part.contains("header")),
        "{searched:?}"
    );
}

/// Ten shards on a store with the latency of S3 (50 ms per request) and a bandwidth, searched
/// with 8 threads, cold and with sketch, locators and headers pinned: the results equal a local
/// search, and the requests, bytes, the most at once and the wall time are printed. Run it with
/// `--ignored` for the numbers.
#[test]
#[ignore]
fn ten_shards_timed() {
    for length in [QUERY_LENGTH, 30_000] {
        timed_search(length, 10);
    }
}

fn timed_search(length: usize, shards: usize) {
    let directory = tempfile::tempdir().unwrap();
    let query = build_shards(directory.path(), length, shards);
    let local: Vec<Shard> = (0..shards)
        .map(|shard| Shard::open(&directory.path().join(format!("shard{shard}.jam"))).unwrap())
        .collect();
    let expected = search_all(&local, &query, 8);
    let speed = Speed {
        latency: Duration::from_millis(50),
        stream_rate: 100 << 20,
        link_rate: 1 << 30,
    };
    let (server, store) = serve_at(directory.path(), speed);
    let s3 = S3Config::new(
        "local",
        Some(&server.address),
        true,
        S3Limits::default(),
        runtime(),
    );
    let s3 = s3.unwrap().unsigned();
    let run = |pins: Option<&Path>, name: &str| {
        store.reset();
        let opened = Instant::now();
        let shards = open_s3(&s3, pins, shards);
        let (open_requests, open_seconds) = (store.log.lock().unwrap().len(), opened.elapsed());
        store.reset();
        let start = Instant::now();
        assert_eq!(search_all(&shards, &query, 8), expected, "{name}");
        let seconds = start.elapsed().as_secs_f64();
        let (parts, peak) = counts(&local, &store);
        let log = store.log.lock().unwrap();
        let bytes: u64 = log.iter().map(|(_, start, end)| end - start).sum();
        eprintln!(
            "{length} bases, {name}: open {open_requests} requests in {:.2} s; search {} \
             requests, {bytes} bytes, at most {peak} at once, {seconds:.2} s; {parts:?}",
            open_seconds.as_secs_f64(),
            log.len()
        );
    };
    run(None, "cold");
    let cache = directory.path().join("cache");
    warm_all(
        &s3,
        &cache,
        &[Part::Sketch, Part::Locators, Part::Headers],
        0,
        shards,
    );
    run(Some(&cache), "warm");
}
