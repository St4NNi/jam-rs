//! Reading an assembly and joining its contigs into batches, which are searched like queries.
// Copyright (c) 2023-2026 Sebastian Beyvers
// SPDX-License-Identifier: MIT

use super::ScreenError;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Bases between two contigs of a batch, so that no seed or alignment joins them.
pub(super) const SPACER: usize = 64;

/// Bytes of the read buffers; a FASTA parser alone reads a gzip file in tiny pieces, which is
/// very slow on a network file system.
const READ_BUFFER: usize = 1 << 20;

/// One piece of a contig inside a batch: its start in the batch, its length, the contig and the
/// position of the piece in the contig. Pieces of a long contig overlap; each owns the positions
/// from the middle of one overlap to the middle of the next, so each alignment is kept once.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Piece {
    pub start: u64,
    pub length: u64,
    pub contig: u32,
    pub offset: u64,
    /// The contig positions the piece owns.
    pub owned: (u64, u64),
    /// Where the previous piece of the contig ends in the contig, or `offset` for the first.
    pub previous_end: u64,
}

impl Piece {
    /// The position after the last base of the piece in the batch.
    pub fn end(&self) -> u64 {
        self.start + self.length
    }
}

/// The bases of a batch, with spacers of N between the pieces, and the pieces.
#[derive(Debug, Default)]
pub(super) struct Batch {
    pub bases: Vec<u8>,
    pub pieces: Vec<Piece>,
}

impl Batch {
    /// Adds `piece` with the bases `sequence`; its start becomes the end of the batch.
    fn add(&mut self, piece: Piece, sequence: &[u8]) {
        if !self.bases.is_empty() {
            self.bases.extend(std::iter::repeat_n(b'N', SPACER));
        }
        self.pieces.push(Piece {
            start: self.bases.len() as u64,
            length: sequence.len() as u64,
            ..piece
        });
        self.bases.extend_from_slice(sequence);
    }

    /// A batch of only the pieces with the given indices, ascending, joined again with spacers.
    pub fn subset(&self, pieces: &[u32]) -> Batch {
        let mut subset = Batch::default();
        for &index in pieces {
            let piece = self.pieces[index as usize];
            subset.add(
                piece,
                &self.bases[piece.start as usize..piece.end() as usize],
            );
        }
        subset
    }
}

/// The pieces of a contig of `length` bases: at most `batch_bases` long, each overlapping the
/// next by `overlap` bases, so that seeds and alignments near a cut have context on both sides.
fn cut_contig(contig: u32, length: u64, batch_bases: u64, overlap: u64) -> Vec<Piece> {
    let overlap = overlap.min(batch_bases / 2);
    let mut pieces = Vec::new();
    let (mut offset, mut previous_end) = (0, 0);
    loop {
        let end = (offset + batch_bases).min(length);
        let owned_start = if offset == 0 { 0 } else { offset + overlap / 2 };
        let owned_end = if end == length {
            length
        } else {
            end - overlap + overlap / 2
        };
        pieces.push(Piece {
            start: 0,
            length: end - offset,
            contig,
            offset,
            owned: (owned_start, owned_end),
            previous_end: if offset == 0 { 0 } else { previous_end },
        });
        if end == length {
            return pieces;
        }
        previous_end = end;
        offset = end - overlap;
    }
}

/// A contig of the assembly that was searched: its name (the first word of the header) and its
/// length.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) struct AssemblyContig {
    pub name: String,
    pub length: u64,
}

/// The contigs of an assembly that were searched, by id, and those left out.
#[derive(Debug, Default)]
pub(super) struct AssemblyContigs {
    pub list: Vec<AssemblyContig>,
    pub bases: u64,
    pub skipped: u64,
}

/// Receives every searched contig: its id, the header line of its record and its bases.
pub(super) type ContigSink<'a> = dyn FnMut(u32, &[u8], &[u8]) -> Result<(), ScreenError> + 'a;

/// Receives every batch with its index.
pub(super) type BatchVisit<'a> = dyn FnMut(usize, &Batch) -> Result<(), ScreenError> + 'a;

/// Passes the contigs of an assembly (plain or gzip FASTA) to `visit` in batches of at most
/// `batch_bases` bases, the same on every read, cutting longer ones in pieces overlapping by
/// `overlap`; skips contigs of fewer than `min_bases` known bases. `None` for an empty file.
pub(super) fn read_batches(
    path: &Path,
    batch_bases: usize,
    overlap: usize,
    min_bases: usize,
    visit: &mut BatchVisit<'_>,
    mut sink: Option<&mut ContigSink<'_>>,
) -> Result<Option<AssemblyContigs>, ScreenError> {
    let failed = |error: &dyn std::fmt::Display| input_error(path, error);
    let Some(mut reader) = open_fasta(path)? else {
        return Ok(None);
    };
    let mut contigs = AssemblyContigs::default();
    let (mut index, mut batch) = (0, Batch::default());
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| failed(&error))?;
        let sequence = record.seq();
        let known = sequence
            .iter()
            .filter(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        if sequence.len() < min_bases || known < min_bases {
            contigs.skipped += 1;
            continue;
        }
        let contig = contigs.list.len() as u32;
        if let Some(sink) = sink.as_mut() {
            sink(contig, record.id(), &sequence)?;
        }
        let name = String::from_utf8_lossy(record.id())
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .to_string();
        contigs.list.push(AssemblyContig {
            name,
            length: sequence.len() as u64,
        });
        contigs.bases += sequence.len() as u64;
        let length = sequence.len() as u64;
        for piece in cut_contig(contig, length, batch_bases as u64, overlap as u64) {
            let bases = &sequence[piece.offset as usize..(piece.offset + piece.length) as usize];
            if !batch.bases.is_empty() && batch.bases.len() + SPACER + bases.len() > batch_bases {
                visit(index, &batch)?;
                index += 1;
                batch = Batch::default();
            }
            batch.add(piece, bases);
        }
    }
    if !batch.bases.is_empty() {
        visit(index, &batch)?;
    }
    Ok(Some(contigs))
}

fn input_error(path: &Path, error: &dyn std::fmt::Display) -> ScreenError {
    ScreenError::Input(format!("{}: {error}", path.display()))
}

/// A FASTA reader of a plain or gzip file through a large buffer; `None` for an empty file.
fn open_fasta(
    path: &Path,
) -> Result<Option<Box<dyn needletail::parser::FastxReader>>, ScreenError> {
    let failed = |error: &dyn std::fmt::Display| input_error(path, error);
    let file = std::fs::File::open(path).map_err(|error| failed(&error))?;
    let mut source = BufReader::with_capacity(READ_BUFFER, file);
    let magic = source.fill_buf().map_err(|error| failed(&error))?;
    let decoded: Box<dyn Read + Send> = if magic.starts_with(&[0x1f, 0x8b]) {
        Box::new(BufReader::with_capacity(
            READ_BUFFER,
            flate2::bufread::MultiGzDecoder::new(source),
        ))
    } else {
        Box::new(source)
    };
    match needletail::parse_fastx_reader(decoded) {
        Ok(reader) => Ok(Some(reader)),
        Err(error) if error.kind == needletail::errors::ParseErrorKind::EmptyFile => Ok(None),
        Err(error) => Err(failed(&error)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Seen = Vec<(Vec<u8>, Vec<Piece>)>;

    fn batches(fasta: &[u8], batch_bases: usize) -> (Seen, AssemblyContigs) {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("assembly.fa");
        std::fs::write(&path, fasta).unwrap();
        let mut seen = Vec::new();
        let contigs = read_batches(
            &path,
            batch_bases,
            0,
            4,
            &mut |index, batch| {
                assert_eq!(index, seen.len());
                seen.push((batch.bases.clone(), batch.pieces.clone()));
                Ok(())
            },
            None,
        )
        .unwrap()
        .unwrap();
        (seen, contigs)
    }

    #[test]
    fn batches_join_contigs() {
        let mut fasta = b">a one\nACGTACGT\n>n\nNNNNNN\n>b\nGGGGCCCC\n>c\n".to_vec();
        fasta.extend_from_slice(&[b'T'; 88]);
        fasta.push(b'\n');
        let (seen, contigs) = batches(&fasta, 90);
        let names: Vec<&str> = contigs.list.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["a", "b", "c"]);
        assert_eq!((contigs.skipped, contigs.bases), (1, 104));
        assert_eq!(seen.len(), 2);
        let (bases, pieces) = &seen[0];
        assert_eq!(bases.len(), 8 + SPACER + 8);
        assert_eq!(&bases[8 + SPACER..], b"GGGGCCCC");
        assert_eq!(pieces[1].start, 8 + SPACER as u64);
        // The contig of 88 bases fits no batch with others, so it starts a new one.
        assert_eq!(
            seen[1].1,
            vec![Piece {
                start: 0,
                length: 88,
                contig: 2,
                offset: 0,
                owned: (0, 88),
                previous_end: 0,
            }]
        );
        let (seen, _) = batches(&fasta, 50);
        let last: Vec<Piece> = seen.iter().flat_map(|(_, pieces)| pieces.clone()).collect();
        assert_eq!(
            last[2..]
                .iter()
                .map(|p| (p.offset, p.length))
                .collect::<Vec<_>>(),
            [(0, 50), (50, 38)]
        );
    }

    /// Pieces overlap, and the parts they own join without a gap or a double.
    #[test]
    fn pieces_overlap() {
        let pieces = cut_contig(3, 1_000_000, 500_000, 16_384);
        let spans: Vec<(u64, u64)> = pieces.iter().map(|p| (p.offset, p.length)).collect();
        assert_eq!(spans, [(0, 500_000), (483_616, 500_000), (967_232, 32_768)]);
        let owned: Vec<(u64, u64)> = pieces.iter().map(|p| p.owned).collect();
        let middle = [491_808, 975_424];
        assert_eq!(
            owned,
            [
                (0, middle[0]),
                (middle[0], middle[1]),
                (middle[1], 1_000_000)
            ]
        );
        let ends: Vec<u64> = pieces.iter().map(|p| p.previous_end).collect();
        assert_eq!(ends, [0, 500_000, 983_616]);
        assert_eq!(cut_contig(0, 100, 500_000, 16_384)[0].owned, (0, 100));
    }

    #[test]
    fn subset_keeps_offsets() {
        let mut batch = Batch::default();
        let piece = |contig, offset| Piece {
            start: 0,
            length: 0,
            contig,
            offset,
            owned: (offset, offset + 4),
            previous_end: offset,
        };
        batch.add(piece(0, 0), b"ACGT");
        batch.add(piece(1, 0), b"ACGT");
        batch.add(piece(2, 7), b"GGCC");
        let subset = batch.subset(&[0, 2]);
        assert_eq!(subset.bases.len(), 8 + SPACER);
        assert_eq!(&subset.bases[4 + SPACER..], b"GGCC");
        assert_eq!((subset.pieces[1].contig, subset.pieces[1].offset), (2, 7));
    }
}
