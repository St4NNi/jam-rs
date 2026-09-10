//! Test data shared by all tests: seeded random bases and bytes, mutated copies, test shards, and
//! a small HTTP server on a local port that plays an S3 store.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use crate::index::{BuildOptions, EntryInput, Parameters, Shard, build_shard};
use crate::loss::MissModel;
use crate::search::SearchControl;
use crate::seeds::SeedScheme;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::net::TcpListener;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};
use tokio::runtime::{Handle, Runtime};

/// A random generator that gives the same values for the same seed on every run.
pub fn rng(seed: u64) -> StdRng {
    StdRng::seed_from_u64(seed)
}

/// `length` random bases, each of `A`, `C`, `G` and `T` with the same chance.
pub fn random_bases(rng: &mut StdRng, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| b"ACGT"[rng.random_range(0..4)])
        .collect()
}

/// `length` random bytes.
pub fn random_bytes(rng: &mut StdRng, length: usize) -> Vec<u8> {
    let mut bytes = vec![0; length];
    rng.fill(&mut bytes[..]);
    bytes
}

/// The handle of one tokio runtime shared by all tests.
pub fn runtime() -> Handle {
    static RUNTIME: OnceLock<Runtime> = OnceLock::new();
    RUNTIME
        .get_or_init(|| Runtime::new().unwrap())
        .handle()
        .clone()
}

/// A copy of `bases` in which each base becomes another base with chance `1 - identity`.
pub fn with_substitutions(rng: &mut StdRng, bases: &[u8], identity: f64) -> Vec<u8> {
    bases
        .iter()
        .map(|&base| {
            if rng.random::<f64>() < identity {
                return base;
            }
            let others: Vec<u8> = b"ACGT"
                .iter()
                .copied()
                .filter(|&other| other != base)
                .collect();
            others[rng.random_range(0..others.len())]
        })
        .collect()
}

/// The FASTA text of the contigs of entry `entry`, named `entry{entry}_contig{index}`.
pub fn fasta_text(entry: usize, contigs: &[Vec<u8>]) -> Vec<u8> {
    let mut text = Vec::new();
    for (index, bases) in contigs.iter().enumerate() {
        text.extend_from_slice(format!(">entry{entry}_contig{index} description\n").as_bytes());
        for line in bases.chunks(70) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    text
}

/// Writes each entry (a list of contigs) to `directory/entry{index}.fa` and returns the inputs,
/// named `entry{index}`.
pub fn write_entries(directory: &Path, entries: &[Vec<Vec<u8>>]) -> Vec<EntryInput> {
    let inputs = entries.iter().enumerate().map(|(index, contigs)| {
        let path = directory.join(format!("entry{index}.fa"));
        std::fs::write(&path, fasta_text(index, contigs)).unwrap();
        EntryInput::from_path(path)
    });
    inputs.collect()
}

/// The parameters of the test shards: k 21, w 20 and bins of 256 bases, so short queries have
/// enough seeds.
pub fn parameters() -> Parameters {
    Parameters {
        scheme: SeedScheme {
            k: 21,
            w: 20,
            ..SeedScheme::default()
        },
        bin_bits: 8,
        ..Parameters::default()
    }
}

/// The build options of the test shards, with `parameters`.
pub fn build_options(directory: &Path, threads: usize) -> BuildOptions {
    BuildOptions {
        parameters: parameters(),
        bucket_bits: None,
        threads,
        temporary_directory: directory.to_path_buf(),
        entry_per_record: false,
    }
}

/// Builds and opens the shard `directory/shard` of `entries`, each a list of contigs.
pub fn test_shard(directory: &Path, entries: &[Vec<Vec<u8>>]) -> Shard {
    let inputs = write_entries(directory, entries);
    let output = directory.join("shard");
    build_shard(&inputs, &output, &build_options(directory, 2)).unwrap();
    Shard::open(&output).unwrap()
}

/// A copy of `bases` with exactly `1 - identity` of its bases changed, at random places.
pub fn diverged(rng: &mut StdRng, bases: &[u8], identity: f64) -> Vec<u8> {
    let mut copy = bases.to_vec();
    let changes = ((1.0 - identity) * bases.len() as f64).round() as usize;
    let mut places: Vec<usize> = (0..bases.len()).collect();
    for index in 0..changes {
        places.swap(index, rng.random_range(index..bases.len()));
        let at = places[index];
        let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != copy[at]).collect();
        copy[at] = others[rng.random_range(0..others.len())];
    }
    copy
}

/// The miss model of the test shards, built once per test process. 50 layouts per length keep
/// the build short and are enough for tests.
pub fn miss_model() -> &'static MissModel {
    static MODEL: OnceLock<MissModel> = OnceLock::new();
    MODEL.get_or_init(|| MissModel::with_layouts(parameters().scheme, 50).unwrap())
}

/// A search control without a deadline whose stop flag is never set.
pub fn no_stop() -> SearchControl<'static> {
    static NEVER: AtomicBool = AtomicBool::new(false);
    SearchControl {
        deadline: None,
        stop: &NEVER,
    }
}

/// A request as the server saw it.
pub struct Request {
    pub method: String,
    pub path: String,
    pub headers: HashMap<String, String>,
}

impl Request {
    /// The first and last byte of the requested range.
    pub fn range(&self) -> (usize, usize) {
        let value = self.headers["range"].strip_prefix("bytes=").unwrap();
        let (start, last) = value.split_once('-').unwrap();
        (start.parse().unwrap(), last.parse().unwrap())
    }
}

/// An answer: status, headers, the length it declares and the body.
pub struct Reply {
    pub status: u16,
    pub headers: Vec<(String, String)>,
    pub declared: usize,
    pub body: Vec<u8>,
}

impl Reply {
    pub fn status(status: u16) -> Self {
        let headers = Vec::new();
        let body = Vec::new();
        Self {
            status,
            headers,
            declared: 0,
            body,
        }
    }

    /// The answer to a range request for `data` with `etag`, which honours `If-Match`.
    pub fn range_of(data: &[u8], etag: &str, request: &Request) -> Self {
        if request
            .headers
            .get("if-match")
            .is_some_and(|wanted| wanted != etag)
        {
            return Self::status(412);
        }
        let (start, last) = request.range();
        let body = data[start..=last].to_vec();
        let range = format!("bytes {start}-{last}/{}", data.len());
        Self {
            status: 206,
            headers: vec![
                ("Content-Range".into(), range),
                ("ETag".into(), etag.into()),
            ],
            declared: body.len(),
            body,
        }
    }
}

/// A server that answers each request with `handler(index, request)`, where `index` counts
/// requests from 0. It stops when dropped.
pub struct Server {
    pub address: String,
    pub requests: Arc<AtomicUsize>,
    stop: Arc<AtomicBool>,
}

impl Drop for Server {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::SeqCst);
        // One more connection wakes the accepting thread, which then sees the stop flag.
        let _ = std::net::TcpStream::connect(self.address.trim_start_matches("http://"));
    }
}

/// Starts a server on a free local port.
pub fn serve(handler: impl Fn(usize, &Request) -> Reply + Send + Sync + 'static) -> Server {
    let listener = TcpListener::bind("127.0.0.1:0").unwrap();
    let server = Server {
        address: format!("http://{}", listener.local_addr().unwrap()),
        requests: Arc::new(AtomicUsize::new(0)),
        stop: Arc::new(AtomicBool::new(false)),
    };
    let (requests, stop) = (Arc::clone(&server.requests), Arc::clone(&server.stop));
    let handler = Arc::new(handler);
    std::thread::spawn(move || {
        for stream in listener.incoming() {
            if stop.load(Ordering::SeqCst) {
                break;
            }
            let Ok(stream) = stream else {
                continue;
            };
            let index = requests.fetch_add(1, Ordering::SeqCst);
            let handler = Arc::clone(&handler);
            std::thread::spawn(move || answer(stream, index, &*handler));
        }
    });
    server
}

fn answer(
    mut stream: std::net::TcpStream,
    index: usize,
    handler: &dyn Fn(usize, &Request) -> Reply,
) {
    let mut reader = BufReader::new(stream.try_clone().unwrap());
    let mut first = String::new();
    reader.read_line(&mut first).unwrap();
    let mut parts = first.split_whitespace();
    let method = parts.next().unwrap().to_string();
    let path = parts.next().unwrap().to_string();
    let mut headers = HashMap::new();
    loop {
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        let Some((name, value)) = line.trim_end().split_once(':') else {
            break;
        };
        headers.insert(name.to_ascii_lowercase(), value.trim().to_string());
    }
    let request = Request {
        method,
        path,
        headers,
    };
    let reply = handler(index, &request);
    let mut head = format!(
        "HTTP/1.1 {} reply\r\nContent-Length: {}\r\nConnection: close\r\n",
        reply.status, reply.declared
    );
    for (name, value) in &reply.headers {
        head.push_str(&format!("{name}: {value}\r\n"));
    }
    // The client may give up and close the connection, so write errors are fine.
    let _ = stream.write_all(format!("{head}\r\n").as_bytes());
    let _ = stream.write_all(&reply.body);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeds_repeat_bases() {
        assert_eq!(random_bases(&mut rng(7), 50), random_bases(&mut rng(7), 50));
        assert_ne!(random_bases(&mut rng(7), 50), random_bases(&mut rng(8), 50));
    }

    #[test]
    fn substitutions_follow_identity() {
        let mut rng = rng(1);
        let bases = random_bases(&mut rng, 100_000);
        let copy = with_substitutions(&mut rng, &bases, 0.9);
        let same = bases.iter().zip(&copy).filter(|(a, b)| a == b).count();
        assert!((89_000..91_000).contains(&same), "{same}");
        assert_eq!(with_substitutions(&mut rng, &bases, 1.0), bases);
    }

    #[test]
    fn diverged_changes_exactly() {
        let mut rng = rng(2);
        let bases = random_bases(&mut rng, 1_000);
        let copy = diverged(&mut rng, &bases, 0.95);
        assert_eq!(bases.iter().zip(&copy).filter(|(a, b)| a != b).count(), 50);
    }
}
