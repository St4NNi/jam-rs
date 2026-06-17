use super::*;
use crate::service::client::Client;
use crate::service::test_shards::{Running, build, config, open, request};
use std::sync::atomic::AtomicBool;

const TIMEOUT: Duration = Duration::from_secs(10);

/// Waits until `check` holds, for at most ten seconds.
fn wait_until(check: impl Fn() -> bool) {
    let started = Instant::now();
    while !check() {
        assert!(started.elapsed() < TIMEOUT, "the condition never held");
        std::thread::sleep(Duration::from_millis(10));
    }
}

fn limited(searches: usize, queue: usize, connections: usize) -> ServerConfig {
    ServerConfig {
        limits: ServerLimits {
            connections,
            searches,
            queue,
        },
        search: config(),
    }
}

#[test]
fn searches_over_the_limit_wait_or_are_refused_and_a_closed_client_leaves_the_queue() {
    let directory = tempfile::tempdir().unwrap();
    let path = build(directory.path(), "alpha", 11);
    let running = Running::alone(open(&path), limited(1, 1, 8));
    let admission = &running.server.admission;
    let stop = AtomicBool::new(false);
    let control = SearchControl {
        deadline: None,
        stop: &stop,
    };
    let place = admission.search(control).unwrap();

    let mut waiting = Client::connect(&running.address, TIMEOUT).unwrap();
    let frames = waiting.search(&request(None), None).unwrap();
    wait_until(|| admission.searches() == (1, 1));

    let mut refused = Client::connect(&running.address, TIMEOUT).unwrap();
    let error = refused.search(&request(None), None).unwrap().next();
    assert!(matches!(error, Some(Err(ServiceError::Refused(message))) if message.contains("busy")));

    drop(frames);
    drop(waiting);
    wait_until(|| admission.searches() == (1, 0));
    drop(place);
    let frames: Vec<Frame> = refused
        .search(&request(None), None)
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(frames.last(), Some(&Frame::Complete { complete: true }));
}

#[test]
fn connections_over_the_limit_are_refused() {
    let directory = tempfile::tempdir().unwrap();
    let running = Running::alone(
        open(&build(directory.path(), "alpha", 11)),
        limited(1, 1, 1),
    );
    let first = Client::connect(&running.address, TIMEOUT).unwrap();
    let error = Client::connect(&running.address, TIMEOUT).err();
    assert!(matches!(error, Some(ServiceError::Refused(_))));
    drop(first);
    wait_until(|| Client::connect(&running.address, TIMEOUT).is_ok());
}

/// A peer that says hello for `identity` and then never answers a search. It reports whether the
/// coordinator closed the stalled search connection.
fn stalled_peer(listener: TcpListener, identity: ShardIdentity) -> bool {
    let hello = Frame::Hello {
        protocol_version: PROTOCOL_VERSION,
        shards: vec![identity],
    };
    // The first connection only pins the shard.
    let mut pinning = Connection::new(listener.accept().unwrap().0).unwrap();
    pinning.write(&hello).unwrap();
    drop(pinning);
    let mut search = Connection::new(listener.accept().unwrap().0).unwrap();
    search.write(&hello).unwrap();
    assert!(matches!(search.read(), Ok(Some(Frame::Search { .. }))));
    search.set_read_deadline(Some(Instant::now() + TIMEOUT));
    matches!(search.read(), Ok(None))
}

#[test]
fn a_stalled_peer_holds_back_only_its_shard_and_a_closed_client_stops_it() {
    let directory = tempfile::tempdir().unwrap();
    let local = open(&build(directory.path(), "alpha", 11));
    let identity = open(&build(directory.path(), "beta", 22))
        .identity()
        .clone();
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let address = listener.local_addr().unwrap().to_string();
    std::thread::scope(|scope| {
        let peer = scope.spawn(|| stalled_peer(listener, identity));
        let peers = vec![Peer::connect(&address, TIMEOUT).unwrap()];
        let server = Server::with_peers(local, limited(1, 1, 8), peers).unwrap();
        let running = Running::start(server);

        let mut client = Client::connect(&running.address, TIMEOUT).unwrap();
        let mut frames = client.search(&request(None), None).unwrap();
        // The local shard ends while the peer stalls.
        let local_summary = frames.find(|frame| {
            matches!(frame, Ok(Frame::Summary { shard: 0, summary }) if summary.is_complete())
        });
        assert!(local_summary.is_some());
        drop(frames);
        drop(client);
        assert!(
            peer.join().unwrap(),
            "the stalled peer connection stayed open"
        );
        wait_until(|| running.server.admission.searches() == (0, 0));
    });
}
