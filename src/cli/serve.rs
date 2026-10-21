//! `jam serve`: serve one shard over iroh until SIGINT or SIGTERM, as the coordinator of its peers
//! when it has some.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::CliError;
use super::storage::{CacheArgs, ShardStorage, shard_name};
use crate::search::SearchConfig;
use crate::service::{
    Reopen, ServedShard, Server, ServerLimits, ServiceError, bind, parse_address,
};
use clap::Args;
use iroh::{EndpointAddr, SecretKey};
use std::fs::OpenOptions;
use std::io::{ErrorKind, Write};
use std::net::SocketAddr;
use std::os::unix::fs::OpenOptionsExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use tokio::signal::unix::{SignalKind, signal};

/// The options of `jam serve`.
#[derive(Debug, Args)]
pub struct ServeArgs {
    /// The shard: NAME.jam or s3://bucket/prefix/NAME.jam
    shard: PathBuf,
    /// The secret key file; its public key is the node id. It is created on first start
    #[arg(long)]
    key: PathBuf,
    /// A peer server as NODE_ID@ADDRESS; every search also runs there. Repeat for more peers
    #[arg(long, value_parser = parse_address)]
    peer: Vec<EndpointAddr>,
    /// A UDP address to listen on; repeat to add an IPv6 one
    #[arg(long, default_value = "0.0.0.0:4100")]
    bind: Vec<SocketAddr>,
    /// Searches that run at once
    #[arg(long, default_value_t = 4)]
    searches: usize,
    /// Searches that wait for a place; more are refused as busy
    #[arg(long, default_value_t = 32)]
    queue: usize,
    /// Threads of each search [default: all cores]
    #[arg(short, long)]
    threads: Option<usize>,
    #[command(flatten)]
    cache: CacheArgs,
}

/// Opens the shard, prints `NODE_ID@ADDRESS` for each bound socket on one line, and serves
/// until SIGINT or SIGTERM. Running searches are cut at the stop.
pub fn run(args: &ServeArgs) -> Result<(), CliError> {
    if args.searches == 0 || args.threads == Some(0) {
        return Err(CliError::Usage(
            "--searches and --threads must be at least 1".into(),
        ));
    }
    let key = secret_key(&args.key)?;
    let mut storage = ShardStorage::new()?;
    storage.set_searches(args.searches);
    let storage = Arc::new(Mutex::new(storage));
    let mut shards =
        lock(&storage).open_all(std::slice::from_ref(&args.shard), &args.cache, true)?;
    let shard = shards
        .pop()
        .ok_or_else(|| CliError::Usage("no shard".into()))?;
    let name = shard_name(&args.shard);
    let shard = ServedShard::new(name, shard)?.with_reopen(reopen(&storage, args));
    let config = SearchConfig {
        threads: args.threads.unwrap_or(0),
        ..SearchConfig::default()
    };
    let limits = ServerLimits {
        searches: args.searches,
        queue: args.queue,
    };
    let runtime = tokio::runtime::Runtime::new().map_err(ServiceError::from)?;
    runtime.block_on(async {
        let endpoint = bind(key, &args.bind).await?;
        let id = endpoint.id();
        let sockets = endpoint.bound_sockets();
        let addresses: Vec<String> = sockets.iter().map(|at| format!("{id}@{at}")).collect();
        println!("{}", addresses.join(" "));
        let router = Server::new(endpoint, shard, args.peer.clone(), config, limits).spawn();
        let stopped = stop_signal().await;
        let shutdown = router.shutdown().await;
        stopped.map_err(ServiceError::from)?;
        shutdown.map_err(|error| ServiceError::Io(std::io::Error::other(error)))?;
        Ok(())
    })
}

fn lock(storage: &Mutex<ShardStorage>) -> std::sync::MutexGuard<'_, ShardStorage> {
    storage
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// Opens the shard again, with a fresh read of its manifest, for when an object of it on S3 was
/// replaced.
fn reopen(storage: &Arc<Mutex<ShardStorage>>, args: &ServeArgs) -> Reopen {
    let (storage, shard, cache) = (
        Arc::downgrade(storage),
        args.shard.clone(),
        args.cache.clone(),
    );
    Box::new(move || {
        let storage = storage.upgrade()?;
        let opened = lock(&storage).open_all(std::slice::from_ref(&shard), &cache, true);
        opened
            .map(|mut shards| shards.pop())
            .unwrap_or_else(|error| {
                eprintln!("jam serve: the shard could not be opened again: {error}");
                None
            })
    })
}

/// The secret key in `path` (32 bytes). On first start a new key is written there, readable by
/// its owner only.
fn secret_key(path: &Path) -> Result<SecretKey, CliError> {
    let file_error = |source| CliError::File {
        path: path.to_path_buf(),
        source,
    };
    match std::fs::read(path) {
        Ok(bytes) => {
            let bytes: [u8; 32] = bytes.try_into().map_err(|_| {
                CliError::Usage(format!("{} must hold a 32-byte key", path.display()))
            })?;
            Ok(SecretKey::from_bytes(&bytes))
        }
        Err(error) if error.kind() == ErrorKind::NotFound => {
            let key = SecretKey::generate();
            let mut file = (OpenOptions::new().write(true).create_new(true).mode(0o600))
                .open(path)
                .map_err(file_error)?;
            file.write_all(&key.to_bytes()).map_err(file_error)?;
            Ok(key)
        }
        Err(error) => Err(file_error(error)),
    }
}

/// Waits for SIGINT or SIGTERM.
async fn stop_signal() -> std::io::Result<()> {
    let mut terminate = signal(SignalKind::terminate())?;
    tokio::select! {
        result = tokio::signal::ctrl_c() => result,
        _ = terminate.recv() => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn key_file_persists() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("node.key");
        let first = secret_key(&path).unwrap();
        let mode = std::fs::metadata(&path).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o600);
        assert_eq!(secret_key(&path).unwrap().public(), first.public());
        std::fs::write(&path, b"short").unwrap();
        assert!(matches!(secret_key(&path), Err(CliError::Usage(_))));
    }
}
