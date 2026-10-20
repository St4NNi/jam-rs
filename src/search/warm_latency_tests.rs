//! Timing of warm searches of a plasmid-length query: a local shard, and a shard on a test store
//! with its sketch and locators pinned on disk or in RAM. Run by hand with `--ignored`.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::align_stage::{AlignStageConfig, PreparedQuery};
use super::control::SearchControl;
use super::driver::{SearchConfig, SearchSummary, Searcher};
use super::events::EntryEvent;
use super::query::QueryKmers;
use super::seed_mismatches::lookup_plan;
use super::sketch_stage::keep_found;
use crate::index::{
    BuildOptions, Location, Parameters, Part, Plan, Shard, ShardOptions, build_shard, warm,
};
use crate::loss::MissModel;
use crate::seeds::Seeder;
use crate::storage::{BlockCache, CacheConfig, S3Config, S3Limits};
use crate::test_support::{
    Reply, diverged, no_stop, random_bases, rng, runtime, serve, write_entries,
};
use rand::RngExt;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::mpsc::sync_channel;
use std::time::Instant;

const ENTRIES: usize = 150;
const QUERY_LENGTH: usize = 100_000;

/// One shard with the default seed scheme: 150 entries of three contigs of 20 to 100 kb. The
/// query lies with 1% changes in entries 0 and 1, cut into three pieces, and stretches of 60 to
/// 400 of its bases lie in most other entries.
fn build(directory: &Path) -> Vec<u8> {
    let mut random = rng(0x1a7e);
    let query = random_bases(&mut random, QUERY_LENGTH);
    let mut entries = Vec::new();
    for entry in 0..ENTRIES {
        let mut contigs: Vec<Vec<u8>> = (0..3)
            .map(|_| {
                let length = random.random_range(20_000..100_000);
                random_bases(&mut random, length)
            })
            .collect();
        let planted: Vec<Vec<u8>> = match entry {
            0 | 1 => (diverged(&mut random, &query, 0.99).chunks(QUERY_LENGTH.div_ceil(3)))
                .map(<[u8]>::to_vec)
                .collect(),
            _ if entry % 8 == 7 => Vec::new(),
            _ => {
                let piece = random.random_range(60..400);
                let start = random.random_range(0..QUERY_LENGTH - piece);
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
    let options = BuildOptions {
        parameters: Parameters::default(),
        bucket_bits: None,
        threads: 8,
        temporary_directory: directory.to_path_buf(),
        singleton: false,
        sketch_only: false,
        chunks: 1,
    };
    let inputs = directory.join("input");
    std::fs::create_dir_all(&inputs).unwrap();
    let inputs = write_entries(&inputs, &entries);
    build_shard(&inputs, &directory.join("shard.jam"), &options).unwrap();
    query
}

fn config() -> SearchConfig {
    let mut config = SearchConfig {
        threads: 8,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config
}

/// Searches once: the summary and the events of each entry, in order, as text.
fn search(shard: &Shard, query: &[u8]) -> (SearchSummary, String) {
    let searcher = Searcher::new(shard).unwrap();
    let (sender, receiver) = sync_channel(1 << 16);
    let control: SearchControl<'_> = no_stop();
    let summary = searcher.search(query, false, &config(), control, &sender);
    drop(sender);
    let mut events = BTreeMap::<u32, Vec<EntryEvent>>::new();
    for event in receiver {
        events.entry(event.entry()).or_default().push(event);
    }
    (summary.unwrap(), format!("{events:?}"))
}

/// Prints the parts of three warm searches after one first search, and returns the events.
fn time_searches(name: &str, shard: &Shard, query: &[u8]) -> String {
    let (_, events) = search(shard, query);
    for round in 0..3 {
        let (summary, again) = search(shard, query);
        assert_eq!(again, events, "{name}");
        let times = summary.times;
        let part = |time: &super::PartTime| {
            let reads = time.storage.pinned_reads;
            format!(
                "{:.4} s ({:.4} cpu, {reads} pinned)",
                time.wall_seconds, time.cpu_seconds
            )
        };
        eprintln!(
            "{name} {round}: open {:.4}; prepare {}; sketch {}; entries {}; {} possible",
            times.open_shard_seconds,
            part(&times.prepare_query),
            part(&times.sketch),
            part(&times.entries),
            summary.sketch.possible_entries,
        );
    }
    events
}

/// Times the pieces of the query preparation: the lookups with the sketch filter, and the
/// tables of stage 3.
fn time_prepare(name: &str, shard: &Shard, query: &[u8]) {
    let scheme = shard.manifest().parameters.scheme;
    let seeder = Seeder::new(scheme).unwrap();
    let models = [0, 1, 2].map(|level| MissModel::new(scheme, level).unwrap());
    let config = config();
    let length = query.len() as u64;
    let plan = lookup_plan(&models, &config.rule, config.max_miss, length, false);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    pool.install(|| {
        let started = Instant::now();
        let all = QueryKmers::planned(&seeder, query, false, &plan);
        let made = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let keep = |lookups: &mut Vec<(u64, u32)>| keep_found(shard.seed_file(), lookups);
        let kept = QueryKmers::filtered(&seeder, query, false, &plan, keep).unwrap();
        let filtered = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let prepared = PreparedQuery::new(query, false, &AlignStageConfig::default()).unwrap();
        let tables = started.elapsed().as_secs_f64();
        eprintln!(
            "{name}: level {}; {} lookups, {} kept; made {made:.4} s; filtered {filtered:.4} s; \
             stage 3 tables {tables:.4} s ({} bases)",
            plan.level(),
            all.lookups(),
            kept.distinct().len(),
            prepared.length(),
        );
    });
}

/// Warm searches of a 100 kb query on a local shard and on a test store with pinned parts give
/// the same events; prints the time of each part.
#[test]
#[ignore = "timing, run by hand"]
fn warm_plasmid_timing() {
    let directory = tempfile::tempdir().unwrap();
    let started = Instant::now();
    let query = build(directory.path());
    eprintln!("built in {:.1} s", started.elapsed().as_secs_f64());
    let local = Shard::open(&directory.path().join("shard.jam")).unwrap();
    let sketch = local.seed_file().sketch_part();
    eprintln!("sketch part {} bytes", sketch.end - sketch.start);
    time_prepare("local", &local, &query);
    let expected = time_searches("local", &local, &query);
    let root = directory.path().to_path_buf();
    let server = serve(move |_, request| {
        let key = request.path.strip_prefix("/bucket/").unwrap_or_default();
        match std::fs::read(root.join(key)) {
            Ok(data) => Reply::object(&data, request),
            Err(_) => Reply::status(404),
        }
    });
    let limits = S3Limits::default();
    let s3 = S3Config::new("local", Some(&server.address), true, limits, runtime());
    let s3 = s3.unwrap().unsigned();
    let location = Location::s3("s3://bucket/shard.jam", &s3).unwrap();
    for (name, memory) in [("pinned on disk", 0), ("pinned in RAM", 1u64 << 30)] {
        let cache = directory.path().join(format!("cache-{memory}"));
        let head = location.head().unwrap();
        let plan = Plan::choose(&head, memory, 1 << 30, &[Part::Sketch, Part::Locators]);
        warm(&location, &head, &plan.unwrap(), &cache).unwrap();
        let options = ShardOptions {
            cache: Some(
                BlockCache::open(&CacheConfig::new(256 << 20, 64 << 10), runtime()).unwrap(),
            ),
            pins: Some(cache),
        };
        let shard = Shard::open_with(&location, &options).unwrap();
        time_prepare(name, &shard, &query);
        assert_eq!(time_searches(name, &shard, &query), expected, "{name}");
    }
}
