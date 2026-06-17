//! End-to-end tests over localhost sockets: a server, a client and a coordinator of peers.

use super::client::Client;
use super::coordinator::Peer;
use super::frame::{Frame, PROTOCOL_VERSION, SearchRequest};
use super::records::{CompletionRecord, EntryState};
use super::server::{Server, ServerConfig};
use super::test_shards::{
    Running, build, config, events_of, local_frames, open, request, summary_of,
};
use super::{ServiceError, connection::Connection};
use std::net::TcpListener;
use std::time::Duration;

const TIMEOUT: Duration = Duration::from_secs(10);

fn server_config() -> ServerConfig {
    ServerConfig {
        search: config(),
        ..ServerConfig::default()
    }
}

/// Every frame of one search through a new client of `address`.
fn remote_frames(address: &str, request: &SearchRequest) -> Vec<Frame> {
    let mut client = Client::connect(address, TIMEOUT).unwrap();
    let frames = client.search(request, None).unwrap();
    frames.collect::<Result<_, _>>().unwrap()
}

/// Whether shard `shard` has a confirmed event of entry 0, the planted one.
fn confirms_planted_entry(frames: &[Frame], shard: u32) -> bool {
    frames.iter().any(|frame| {
        matches!(frame, Frame::Event { shard: tag, record }
            if *tag == shard && record.entry == 0
                && matches!(record.state, EntryState::Confirmed(_)))
    })
}

#[test]
fn a_server_answers_with_the_events_of_a_local_search() {
    let directory = tempfile::tempdir().unwrap();
    let path = build(directory.path(), "alpha", 11);
    let expected = local_frames(&open(&path), 0);
    assert!(confirms_planted_entry(&expected, 0));
    let running = Running::alone(open(&path), server_config());

    let frames = remote_frames(&running.address, &request(None));
    assert_eq!(events_of(&frames, 0), events_of(&expected, 0));
    assert_eq!(summary_of(&frames, 0), summary_of(&expected, 0));
    assert_eq!(frames.last(), Some(&Frame::Complete { complete: true }));

    // A connection answers one search after another.
    let mut client = Client::connect(&running.address, TIMEOUT).unwrap();
    assert_eq!(client.shards()[0].name, "alpha");
    for _ in 0..2 {
        let frames: Vec<Frame> = client
            .search(&request(None), None)
            .unwrap()
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(events_of(&frames, 0), events_of(&expected, 0));
    }
}

#[test]
fn a_coordinator_merges_the_events_and_summaries_of_every_shard() {
    let directory = tempfile::tempdir().unwrap();
    let paths: Vec<_> = ["alpha", "beta", "gamma"]
        .iter()
        .zip([11, 22, 33])
        .map(|(name, seed)| build(directory.path(), name, seed))
        .collect();
    let mut servers: Vec<Running> = paths[1..]
        .iter()
        .map(|path| Running::alone(open(path), server_config()))
        .collect();
    let peers: Vec<Peer> = servers
        .iter()
        .map(|running| Peer::connect(&running.address, TIMEOUT).unwrap())
        .collect();
    let coordinator = Server::with_peers(open(&paths[0]), server_config(), peers).unwrap();
    let coordinator = Running::start(coordinator);

    let frames = remote_frames(&coordinator.address, &request(None));
    for (tag, path) in paths.iter().enumerate() {
        let tag = tag as u32;
        let expected = local_frames(&open(path), tag);
        assert!(confirms_planted_entry(&frames, tag));
        assert_eq!(events_of(&frames, tag), events_of(&expected, tag));
        assert_eq!(summary_of(&frames, tag), summary_of(&expected, tag));
    }
    assert_eq!(frames.last(), Some(&Frame::Complete { complete: true }));
    let client = Client::connect(&coordinator.address, TIMEOUT).unwrap();
    let names: Vec<&str> = client.shards().iter().map(|s| s.name.as_str()).collect();
    assert_eq!(names, ["alpha", "beta", "gamma"]);
    drop(client);

    // A lost peer gives `ShardFailed`; the other shards still answer.
    drop(servers.pop());
    let frames = remote_frames(&coordinator.address, &request(None));
    assert!(frames.iter().any(|frame| matches!(
        frame,
        Frame::ShardFailed { shard: 2, message } if message.contains("127.0.0.1")
    )));
    assert!(confirms_planted_entry(&frames, 0));
    assert!(confirms_planted_entry(&frames, 1));
    assert_eq!(frames.last(), Some(&Frame::Complete { complete: false }));
}

#[test]
fn a_deadline_stops_the_search_and_still_ends_every_shard() {
    let directory = tempfile::tempdir().unwrap();
    let running = Running::alone(open(&build(directory.path(), "alpha", 11)), server_config());
    let frames = remote_frames(&running.address, &request(Some(0)));
    let summary = summary_of(&frames, 0);
    assert_eq!(
        summary.completion,
        CompletionRecord::Stopped { deadline: true }
    );
    assert_eq!(frames.last(), Some(&Frame::Complete { complete: false }));
}

#[test]
fn another_protocol_version_is_refused_with_a_clear_error() {
    let directory = tempfile::tempdir().unwrap();
    let running = Running::alone(open(&build(directory.path(), "alpha", 11)), server_config());
    let mut connection = Connection::open(&running.address, TIMEOUT).unwrap();
    assert!(matches!(
        connection.read().unwrap(),
        Some(Frame::Hello { .. })
    ));
    connection
        .write(&Frame::Search {
            protocol_version: PROTOCOL_VERSION + 1,
            request: request(None),
        })
        .unwrap();
    let Some(Frame::Refused { message }) = connection.read().unwrap() else {
        panic!("no refusal");
    };
    assert!(message.contains("protocol version"), "{message}");

    // A server of another version is named as such by the client.
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap().to_string();
    std::thread::scope(|scope| {
        scope.spawn(|| {
            let mut server = Connection::new(listener.accept().unwrap().0).unwrap();
            let hello = Frame::Hello {
                protocol_version: PROTOCOL_VERSION + 1,
                shards: Vec::new(),
            };
            server.write(&hello).unwrap();
        });
        let error = Client::connect(&address, TIMEOUT).err().unwrap();
        assert!(matches!(error, ServiceError::Version { .. }), "{error}");
    });
}

#[test]
fn a_peer_must_answer_for_one_shard_with_the_same_seeds() {
    let directory = tempfile::tempdir().unwrap();
    let path = build(directory.path(), "alpha", 11);
    let running = Running::alone(open(&path), server_config());
    let peer = Peer::connect(&running.address, TIMEOUT).unwrap();
    // The same shard twice is refused.
    let error = Server::with_peers(open(&path), server_config(), vec![peer]).err();
    assert!(matches!(error, Some(ServiceError::Peer { .. })));
    // An address where nobody listens is a peer error.
    let free = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = free.local_addr().unwrap().to_string();
    drop(free);
    assert!(Peer::connect(&address, TIMEOUT).is_err());
}
