//! Sequence files: the packed sequences of the entries of an index one after another, each on a
//! 4 KiB boundary, split into chunks of about equal size; each entry lies in one chunk.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::manifest::FileRecord;
use super::packed::BLOCK_BYTES;
use super::{AtPath, IndexError};
use sha2::{Digest, Sha256};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

/// The extension of an index file.
pub(crate) const INDEX_EXTENSION: &str = "jam";
/// The extension of a sequence file.
pub(crate) const SEQUENCE_EXTENSION: &str = "jamseq";
/// Chunks written at once.
const WRITERS: usize = 8;

/// The names of the sequence files of the index `name`: `NAME.jamseq` for one chunk, else
/// `NAME.1.jamseq` to `NAME.N.jamseq`.
pub(crate) fn sequence_file_names(name: &str, chunks: usize) -> Vec<String> {
    match chunks {
        1 => vec![format!("{name}.{SEQUENCE_EXTENSION}")],
        _ => (1..=chunks)
            .map(|chunk| format!("{name}.{chunk}.{SEQUENCE_EXTENSION}"))
            .collect(),
    }
}

/// Where one entry lies: its chunk and its offset there.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Place {
    pub(crate) chunk: usize,
    pub(crate) offset: u64,
}

/// Places entries of these sizes into chunks, the largest first (the first on ties), each into
/// the chunk with the fewest bytes so far (the first on ties); in a chunk they keep their order,
/// each on a 4 KiB boundary. Also gives the size of each chunk.
pub(crate) fn place_entries(
    sizes: &[u64],
    chunks: usize,
) -> Result<(Vec<Place>, Vec<u64>), IndexError> {
    if chunks == 0 || chunks > sizes.len().max(1) {
        return Err(IndexError::Input(format!(
            "the chunks must be 1 to {}, the entries",
            sizes.len().max(1)
        )));
    }
    let mut order: Vec<usize> = (0..sizes.len()).collect();
    order.sort_by_key(|&entry| (Reverse(sizes[entry]), entry));
    let mut smallest: BinaryHeap<Reverse<(u64, usize)>> =
        (0..chunks).map(|chunk| Reverse((0, chunk))).collect();
    let mut chunk_of = vec![0; sizes.len()];
    for entry in order {
        let Some(Reverse((bytes, chunk))) = smallest.pop() else {
            break;
        };
        chunk_of[entry] = chunk;
        smallest.push(Reverse((bytes + aligned(sizes[entry]), chunk)));
    }
    let mut ends = vec![0u64; chunks];
    let places = (sizes.iter().zip(&chunk_of))
        .map(|(&size, &chunk)| {
            let offset = aligned(ends[chunk]);
            ends[chunk] = offset + size;
            Place { chunk, offset }
        })
        .collect();
    Ok((places, ends))
}

/// `bytes` rounded up to the next 4 KiB boundary.
fn aligned(bytes: u64) -> u64 {
    bytes.next_multiple_of(BLOCK_BYTES)
}

/// The packed sequence of one entry: `bytes` bytes of the file `path` from `offset`.
#[derive(Clone, Debug)]
pub(crate) struct Piece {
    pub(crate) path: PathBuf,
    pub(crate) offset: u64,
    pub(crate) bytes: u64,
}

/// Writes the chunk files `paths` from `pieces` placed by `places`, and returns the record of
/// each chunk under the names `names`. Each chunk is hashed while it is written. With `remove`,
/// each piece's file is removed once it is copied, so the disk holds each byte about once.
pub(crate) fn write_sequence_files(
    paths: &[PathBuf],
    names: &[String],
    pieces: &[Piece],
    places: &[Place],
    remove: bool,
) -> Result<Vec<FileRecord>, IndexError> {
    let next = AtomicUsize::new(0);
    let records = Mutex::new(vec![None; paths.len()]);
    let write = || -> Result<(), IndexError> {
        loop {
            let chunk = next.fetch_add(1, Ordering::Relaxed);
            let Some(path) = paths.get(chunk) else {
                return Ok(());
            };
            let members = (pieces.iter().zip(places)).filter(|(_, place)| place.chunk == chunk);
            let (bytes, sha256) = write_chunk(path, members, remove)?;
            let path = names[chunk].clone();
            let record = FileRecord {
                path,
                bytes,
                sha256,
            };
            crate::storage::lock(&records)[chunk] = Some(record);
        }
    };
    std::thread::scope(|scope| {
        let writers: Vec<_> = (0..WRITERS.min(paths.len()))
            .map(|_| scope.spawn(write))
            .collect();
        (writers.into_iter())
            .map(|writer| {
                writer
                    .join()
                    .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
            })
            .collect::<Result<Vec<()>, IndexError>>()
    })?;
    let records = records.into_inner().unwrap_or_else(|e| e.into_inner());
    Ok(records.into_iter().flatten().collect())
}

/// Writes the pieces of one chunk at their places into a new file at `path`; gives its size
/// and SHA-256 digest.
fn write_chunk<'a>(
    path: &Path,
    members: impl Iterator<Item = (&'a Piece, &'a Place)>,
    remove: bool,
) -> Result<(u64, String), IndexError> {
    let file = std::fs::File::create_new(path).at(path)?;
    let mut out = std::io::BufWriter::with_capacity(1 << 20, file);
    let mut digest = Sha256::new();
    let mut buffer = vec![0u8; 1 << 20];
    let mut written = 0u64;
    for (piece, place) in members {
        let padding = vec![0u8; (place.offset - written) as usize];
        digest.update(&padding);
        out.write_all(&padding).at(path)?;
        let mut source = std::fs::File::open(&piece.path).at(&piece.path)?;
        source.seek(SeekFrom::Start(piece.offset)).at(&piece.path)?;
        let mut left = piece.bytes;
        while left > 0 {
            let length = left.min(buffer.len() as u64) as usize;
            source.read_exact(&mut buffer[..length]).at(&piece.path)?;
            digest.update(&buffer[..length]);
            out.write_all(&buffer[..length]).at(path)?;
            left -= length as u64;
        }
        written = place.offset + piece.bytes;
        if remove {
            std::fs::remove_file(&piece.path).at(&piece.path)?;
        }
    }
    let file = (out.into_inner())
        .map_err(|error| error.into_error())
        .at(path)?;
    file.sync_all().at(path)?;
    Ok((written, hex::encode(digest.finalize())))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::manifest::sha256_file;

    #[test]
    fn chunk_names() {
        assert_eq!(sequence_file_names("a", 1), ["a.jamseq"]);
        assert_eq!(
            sequence_file_names("a", 3),
            ["a.1.jamseq", "a.2.jamseq", "a.3.jamseq"]
        );
    }

    /// Chunks differ by at most one entry, each entry lies in one chunk, and entries keep their
    /// order on 4 KiB boundaries without overlap, whatever the input order.
    #[test]
    fn chunks_balanced() {
        let mut rng = crate::test_support::rng(5);
        for chunks in [1, 2, 3, 7] {
            let sizes: Vec<u64> = (0..40)
                .map(|_| rand::RngExt::random_range(&mut rng, 1..2_000_000))
                .collect();
            let (places, ends) = place_entries(&sizes, chunks).unwrap();
            let largest = *sizes.iter().max().unwrap() + BLOCK_BYTES;
            let (low, high) = (ends.iter().min().unwrap(), ends.iter().max().unwrap());
            assert!(high - low <= largest, "{chunks}: {ends:?}");
            for (chunk, &chunk_end) in ends.iter().enumerate() {
                let mut end = 0;
                for (place, size) in places.iter().zip(&sizes) {
                    if place.chunk == chunk {
                        assert!(place.offset >= end && place.offset % BLOCK_BYTES == 0);
                        end = place.offset + size;
                    }
                }
                assert_eq!(end, chunk_end);
            }
        }
        let (places, _) = place_entries(&[10_000, 30_000, 20_000, 5_000], 2).unwrap();
        let chunks: Vec<usize> = places.iter().map(|place| place.chunk).collect();
        assert_eq!(chunks, [1, 0, 1, 0]);
        assert!(place_entries(&[1, 2], 3).is_err());
        assert!(place_entries(&[1, 2], 0).is_err());
    }

    #[test]
    fn chunks_written() {
        let directory = tempfile::tempdir().unwrap();
        let sizes = [5_000u64, 100, 9_000];
        let pieces: Vec<Piece> = (sizes.iter().enumerate())
            .map(|(entry, &bytes)| {
                let path = directory.path().join(format!("piece{entry}"));
                std::fs::write(&path, vec![entry as u8 + 1; bytes as usize]).unwrap();
                Piece {
                    path,
                    offset: 0,
                    bytes,
                }
            })
            .collect();
        let (places, ends) = place_entries(&sizes, 2).unwrap();
        let names = sequence_file_names("x", 2);
        let paths: Vec<PathBuf> = names.iter().map(|n| directory.path().join(n)).collect();
        let records = write_sequence_files(&paths, &names, &pieces, &places, true).unwrap();
        for (chunk, record) in records.iter().enumerate() {
            let digest = sha256_file(&paths[chunk]).unwrap();
            assert_eq!(
                (record.bytes, &record.sha256),
                (ends[chunk], &digest.sha256)
            );
            assert_eq!(record.path, names[chunk]);
        }
        let first = std::fs::read(&paths[places[1].chunk]).unwrap();
        let at = places[1].offset as usize;
        assert_eq!(first[at..at + 100], [2; 100]);
        assert!(pieces.iter().all(|piece| !piece.path.exists()));
    }
}
