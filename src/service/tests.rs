//! Tests over local iroh endpoints: a server, a coordinator with peers, and clients.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::*;
use crate::decide::QueryCoverage;
use crate::index::{EntryInput, Shard, build_shard};
use crate::report::{Dropped, Record, ShardRecord, Unfinished};
use crate::search::{DropReason, SearchConfig, SearchControl, SearchStage, SeedMismatches};
use crate::test_support::{
    build_options, random_bases, rng, with_substitutions,
};
use iroh::protocol::Router;
use iroh::{Endpoint, EndpointAddr, SecretKey};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;
use tokio::sync::mpsc;

const QUERY_LENGTH: usize = 5_000;
const WAIT: Duration = Duration::from_secs(10);

fn query() -> Vec<u8> {
    random_bases(&mut rng(0x5eed), QUERY_LENGTH)
}

/// Shard `name` with three entries: the query planted at 99% identity, and two unrelated ones.
fn build(directory: &Path, name: &str, seed: u64) -> ServedShard {
    let mut random = rng(seed);
    let mut planted = random_bases(&mut random, 20_000);
    let copy = with_substitutions(&mut random, &query(), 0.99);
    planted.splice(6_000..6_000 + QUERY_LENGTH, copy);
    let sequences = [
        planted,
        random_bases(&mut random, 15_000),
        random_bases(&mut random, 12_000),
    ];
    let inputs: Vec<EntryInput> = (sequences.iter().enumerate())
        .map(|(index, sequence)| {
            let path = directory.join(format!("{name}_entry{index}.fa"));
            let mut text = format!(">{name}_contig{index}\n").into_bytes();
            text.extend_from_slice(sequence);
            std::fs::write(&path, text).unwrap();
            let name = format!("{name}_entry{index}");
            EntryInput { name, path }
        })
        .collect();
    let options = build_options(directory, 2);
    let output = directory.join(format!("{name}.jam"));
    build_shard(&inputs, &output, &options).unwrap();
    ServedShard::new(name.into(), Shard::open(&output).unwrap()).unwrap()
}

/// 90% query coverage at 95% identity, on two threads.
fn config() -> SearchConfig {
    let mut config = SearchConfig {
        threads: 2,
        ..SearchConfig::default()
    };
    config.rule.min_identity = 0.95;
    config.rule.min_fragment_identity = 0.0;
    config.rule.min_query_coverage = QueryCoverage::Share {
        numerator: 9,
        decimal_places: 1,
    };
    config
}

fn request(deadline_ms: Option<u64>) -> SearchRequest {
    SearchRequest {
        query_name: "query".into(),
        bases: query(),
        circular: false,
        settings: SearchSettings::of(&config()),
        deadline_ms,
        local_only: false,
    }
}

/// The JSON line of a record message, with the times of a summary cleared.
fn line(reply: Reply) -> String {
    let mut record = reply.into_record().unwrap();
    if let Record::Summary(summary) = &mut record.record {
        summary.sketch.seconds = 0.0;
        summary.locate.seconds = 0.0;
        summary.align.seconds = 0.0;
        summary.times = Default::default();
    }
    serde_json::to_string(&record).unwrap()
}

/// The sorted lines of a search of `shard` in this process, as the order between entries may
/// change.
fn local_lines(shard: &ServedShard) -> Vec<String> {
    let (replies, mut received) = mpsc::channel(1_024);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    shard.search(&request(None), &config(), control, &replies);
    drop(replies);
    let mut lines = Vec::new();
    while let Some(reply) = received.blocking_recv() {
        lines.push(line(reply));
    }
    lines.sort();
    lines
}

async fn endpoint() -> Endpoint {
    let local = "127.0.0.1:0".parse().unwrap();
    bind(SecretKey::generate(), &[local]).await.unwrap()
}

/// The address to dial `endpoint` at.
fn address(endpoint: &Endpoint) -> EndpointAddr {
    let socket = endpoint.bound_sockets().into_iter().find(|s| s.is_ipv4());
    EndpointAddr::new(endpoint.id()).with_ip_addr(socket.unwrap())
}

async fn serve(shard: ServedShard, peers: Vec<EndpointAddr>, limits: ServerLimits) -> Router {
    let server = Server::new(endpoint().await, shard, peers, config(), limits);
    server.spawn()
}

async fn remote_lines(client: &Endpoint, router: &Router, request: &SearchRequest) -> Vec<String> {
    let mut answer = search(client, address(router.endpoint()), request)
        .await
        .unwrap();
    let mut lines = Vec::new();
    while let Some(reply) = answer.next().await.unwrap() {
        lines.push(line(reply));
    }
    lines.sort();
    lines
}

fn confirms(lines: &[String], shard: &str) -> bool {
    let needle = format!("\"shard\":\"{shard}\",\"kind\":\"confirmed\"");
    lines
        .iter()
        .any(|line| line.starts_with(&format!("{{{needle}")))
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn remote_matches_local() {
    let directory = tempfile::tempdir().unwrap();
    let shard = build(directory.path(), "alpha", 11);
    let expected = tokio::task::block_in_place(|| local_lines(&shard));
    assert!(confirms(&expected, "alpha"), "{expected:?}");
    let router = serve(shard, Vec::new(), ServerLimits::default()).await;
    let client = endpoint().await;
    assert_eq!(
        remote_lines(&client, &router, &request(None)).await,
        expected
    );
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn request_sets_mismatches() {
    let directory = tempfile::tempdir().unwrap();
    let shard = build(directory.path(), "alpha", 11);
    let router = serve(shard, Vec::new(), ServerLimits::default()).await;
    let client = endpoint().await;
    for (mode, mismatches) in [(SeedMismatches::Zero, 0), (SeedMismatches::One, 1)] {
        let mut request = request(None);
        request.settings.seed_mismatches = mode;
        let lines = remote_lines(&client, &router, &request).await;
        assert!(confirms(&lines, "alpha"), "{lines:?}");
        let summary = lines
            .iter()
            .find(|line| line.contains("\"kind\":\"summary\""));
        let summary: serde_json::Value = serde_json::from_str(summary.unwrap()).unwrap();
        assert_eq!(summary["seed_mismatches"], mismatches);
    }
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn coordinator_merges_shards() {
    let directory = tempfile::tempdir().unwrap();
    let names = [("alpha", 11), ("beta", 22), ("gamma", 33)];
    let shards: Vec<_> = (names.iter())
        .map(|(name, seed)| build(directory.path(), name, *seed))
        .collect();
    let mut expected: Vec<String> =
        tokio::task::block_in_place(|| shards.iter().flat_map(local_lines).collect());
    expected.sort();
    let mut shards = shards.into_iter();
    let own = shards.next().unwrap();
    let mut peers = Vec::new();
    for shard in shards {
        peers.push(serve(shard, Vec::new(), ServerLimits::default()).await);
    }
    let addresses = peers.iter().map(|peer| address(peer.endpoint())).collect();
    let coordinator = serve(own, addresses, ServerLimits::default()).await;
    let lines = remote_lines(&endpoint().await, &coordinator, &request(None)).await;
    assert_eq!(lines, expected);
    for (name, _) in names {
        assert!(confirms(&lines, name), "{name}");
    }
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn failing_peer_fails() {
    let directory = tempfile::tempdir().unwrap();
    // An endpoint that does not speak the protocol, so the handshake fails.
    let silent = endpoint().await;
    let peers = vec![address(&silent)];
    let coordinator = serve(
        build(directory.path(), "alpha", 11),
        peers,
        Default::default(),
    )
    .await;
    let lines = remote_lines(&endpoint().await, &coordinator, &request(None)).await;
    assert!(confirms(&lines, "alpha"));
    let failed = format!("{{\"shard\":\"{}\",\"kind\":\"summary\"", silent.id());
    let failed: Vec<_> = lines.iter().filter(|l| l.starts_with(&failed)).collect();
    assert_eq!(failed.len(), 1, "{lines:?}");
    assert!(failed[0].contains("\"status\":\"failed\""), "{}", failed[0]);
}

/// A server that runs no search and lets one wait, so the second search is busy.
async fn full_server(directory: &Path) -> Router {
    let limits = ServerLimits {
        searches: 0,
        queue: 1,
    };
    serve(build(directory, "alpha", 11), Vec::new(), limits).await
}

/// Whether a new search is refused as busy within a short wait.
async fn refused(client: &Endpoint, router: &Router) -> bool {
    let mut answer = (search(client, address(router.endpoint()), &request(None)).await).unwrap();
    let next = tokio::time::timeout(Duration::from_millis(500), answer.next()).await;
    matches!(next, Ok(Err(ServiceError::Refused(reason))) if reason == server::BUSY)
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn busy_server_refuses() {
    let directory = tempfile::tempdir().unwrap();
    let router = full_server(directory.path()).await;
    let client = endpoint().await;
    let server_address = address(router.endpoint());
    let mut first = search(&client, server_address.clone(), &request(None))
        .await
        .unwrap();
    let mut second = search(&client, server_address, &request(None))
        .await
        .unwrap();
    // Whichever search arrives first takes the only place and waits; the other is refused.
    let answer = tokio::select! {
        answer = first.next() => answer,
        answer = second.next() => answer,
    };
    assert!(matches!(answer, Err(ServiceError::Refused(reason)) if reason == server::BUSY));
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn disconnect_frees_place() {
    let directory = tempfile::tempdir().unwrap();
    let router = full_server(directory.path()).await;
    let client = endpoint().await;
    let waiting = search(&client, address(router.endpoint()), &request(None)).await;
    let started = tokio::time::Instant::now();
    while !refused(&client, &router).await {
        assert!(started.elapsed() < WAIT, "no search was refused");
    }
    drop(waiting);
    // The place of the dropped search frees up, so a new search waits instead of being refused.
    while refused(&client, &router).await {
        assert!(started.elapsed() < 2 * WAIT, "the place never came free");
    }
}

/// A search that waits for a place is refused at its deadline instead of waiting on.
#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn queued_search_expires() {
    let directory = tempfile::tempdir().unwrap();
    let router = full_server(directory.path()).await;
    let client = endpoint().await;
    let mut answer =
        (search(&client, address(router.endpoint()), &request(Some(200))).await).unwrap();
    let next = tokio::time::timeout(WAIT, answer.next()).await;
    assert!(
        matches!(next, Ok(Err(ServiceError::Refused(ref reason))) if reason == server::EXPIRED),
        "{next:?}"
    );
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn deadline_ends_search() {
    let directory = tempfile::tempdir().unwrap();
    let router = serve(
        build(directory.path(), "alpha", 11),
        Vec::new(),
        Default::default(),
    )
    .await;
    let lines = remote_lines(&endpoint().await, &router, &request(Some(0))).await;
    let summary = lines
        .iter()
        .find(|line| line.contains("\"kind\":\"summary\""));
    assert!(
        summary.unwrap().contains("\"reason\":\"deadline\""),
        "{lines:?}"
    );
}

#[test]
fn closed_channel_stops() {
    let directory = tempfile::tempdir().unwrap();
    let shard = build(directory.path(), "alpha", 11);
    let (replies, received) = mpsc::channel(1);
    drop(received);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    shard.search(&request(None), &config(), control, &replies);
    assert!(stop.load(Ordering::Relaxed));
}

#[test]
fn records_cross_postcard() {
    let directory = tempfile::tempdir().unwrap();
    let shard = build(directory.path(), "alpha", 11);
    let (replies, mut received) = mpsc::channel(1_024);
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    shard.search(&request(None), &config(), control, &replies);
    drop(replies);
    let mut records = Vec::new();
    while let Some(reply) = received.blocking_recv() {
        records.push(reply.into_record().unwrap());
    }
    let (query, entry) = (String::from("query"), String::from("alpha_entry1"));
    let stage = SearchStage::Align;
    let reason = DropReason::Decision;
    for record in [
        Record::Dropped(Dropped {
            query: query.clone(),
            entry: entry.clone(),
            stage,
            reason,
        }),
        Record::Unfinished(Unfinished {
            query,
            entry,
            stage,
        }),
    ] {
        let shard = "alpha".into();
        records.push(ShardRecord { shard, record });
    }
    let mut kinds = Vec::new();
    for record in records {
        let json = serde_json::to_string(&record).unwrap();
        let bytes = postcard::to_allocvec(&Reply::record(record.shard, record.record)).unwrap();
        let back: Reply = postcard::from_bytes(&bytes).unwrap();
        let back = back.into_record().unwrap();
        assert_eq!(serde_json::to_string(&back).unwrap(), json);
        let value: serde_json::Value = serde_json::from_str(&json).unwrap();
        kinds.push(value["kind"].as_str().unwrap().to_string());
    }
    kinds.sort();
    kinds.dedup();
    let all = [
        "confirmed",
        "dropped",
        "located",
        "possible",
        "summary",
        "unfinished",
    ];
    assert_eq!(kinds, all);
    for message in [
        Reply::End,
        Reply::Refused {
            reason: "busy".into(),
        },
    ] {
        let bytes = postcard::to_allocvec(&message).unwrap();
        assert_eq!(postcard::from_bytes::<Reply>(&bytes).unwrap(), message);
    }
}

/// A shard on S3 whose objects are uploaded again: the search that meets the change fails with
/// the retryable error, the server opens the shard again, and the next search gives the results
/// of the new shard.
#[test]
fn replaced_shard_reopened() {
    use crate::index::{Location, ShardOptions};
    use crate::storage::{S3Config, S3Limits};
    use crate::test_support::{runtime, serve_seen};
    let (served, new) = (tempfile::tempdir().unwrap(), tempfile::tempdir().unwrap());
    build(served.path(), "shard", 1);
    let expected = local_lines(&build(new.path(), "shard", 2));
    let (server, _) = serve_seen(served.path());
    let endpoint = Some(server.address.as_str());
    let config = S3Config::new("local", endpoint, true, S3Limits::default(), runtime());
    let location = Location::s3("s3://bucket/shard.jam", &config.unwrap().unsigned()).unwrap();
    let open = move || Shard::open_with(&location, &ShardOptions::default());
    let shard = ServedShard::new("shard".into(), open().unwrap()).unwrap();
    let shard = shard.with_reopen(Box::new(move || open().ok()));
    let before = local_lines(&shard);
    assert_ne!(before, expected);
    assert_eq!(local_lines(&shard), before);
    for name in ["shard.jam", "shard.jamseq"] {
        std::fs::copy(new.path().join(name), served.path().join(name)).unwrap();
    }
    let failed = local_lines(&shard).concat();
    assert!(failed.contains("retry"), "{failed}");
    assert_eq!(local_lines(&shard), expected);
}
