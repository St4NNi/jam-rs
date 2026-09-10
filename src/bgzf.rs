use crate::bgzf_cache::{BgzfBlockCache, BgzfBlockIdentity, MAX_BGZF_BLOCK_BYTES};
use crate::jidx::sha256;
use crate::jidx_reader::{Contig, Metagenome, MetagenomeId};
use crate::range_source::{RangeSource, RangeSourceError, RangeStats, S3Config};
use noodles_bgzf::{self as bgzf, gzi};
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use thiserror::Error;

const MAX_UNCOMPRESSED_BLOCK_BYTES: usize = MAX_BGZF_BLOCK_BYTES;

pub struct BgzfReader {
    metagenome_id: MetagenomeId,
    reader: bgzf::io::Reader<RangeSource>,
    index: gzi::Index,
    cached_block: Option<CachedBlock>,
    shared_cache: Option<Arc<BgzfBlockCache>>,
    cache_source: CacheSource,
    local_path: Option<PathBuf>,
    blocks_decoded: u64,
    decode_nanoseconds: Option<u64>,
}

#[derive(Clone, Copy)]
struct CacheSource {
    source_sha256: [u8; 32],
    source_bytes: u64,
    locator_sha256: [u8; 32],
    local_file: Option<[u64; 7]>,
}

struct CachedBlock {
    compressed_offset: u64,
    uncompressed_offset: u64,
    data: Vec<u8>,
}

impl BgzfReader {
    pub fn open(
        source: Metagenome<'_>,
        s3: Option<&S3Config>,
        verify: bool,
    ) -> Result<Self, BgzfError> {
        Self::from_sources(
            source,
            RangeSource::open(source.bgzf_uri, source.bgzf_bytes, s3)?,
            verify,
            None,
        )
    }

    pub fn open_with_cache(
        source: Metagenome<'_>,
        s3: Option<&S3Config>,
        verify: bool,
        cache: Arc<BgzfBlockCache>,
    ) -> Result<Self, BgzfError> {
        Self::from_sources(
            source,
            RangeSource::open(source.bgzf_uri, source.bgzf_bytes, s3)?,
            verify,
            Some(cache),
        )
    }

    fn from_sources(
        source: Metagenome<'_>,
        mut bgzf_source: RangeSource,
        verify: bool,
        shared_cache: Option<Arc<BgzfBlockCache>>,
    ) -> Result<Self, BgzfError> {
        if bgzf_source.len() != source.bgzf_bytes {
            return Err(BgzfError::SizeMismatch);
        }
        if verify && !bgzf_source.verify_sha256(source.bgzf_sha256)? {
            return Err(BgzfError::ChecksumMismatch);
        }
        let mut index_reader = gzi::io::Reader::new(source.gzi);
        let index = index_reader.read_index()?;
        validate_gzi(&index, source.bgzf_bytes)?;
        let (local_path, local_file, locator_sha256) = if shared_cache.is_some() {
            let path = local_path(source.bgzf_uri);
            let identity = path.as_deref().map(file_identity).transpose()?;
            (path, identity, sha256(source.bgzf_uri.as_bytes()))
        } else {
            (None, None, [0; 32])
        };
        Ok(Self {
            metagenome_id: source.id,
            reader: bgzf::io::Reader::new(bgzf_source),
            index,
            cached_block: None,
            shared_cache,
            cache_source: CacheSource {
                source_sha256: source.bgzf_sha256,
                source_bytes: source.bgzf_bytes,
                locator_sha256,
                local_file,
            },
            local_path,
            blocks_decoded: 0,
            decode_nanoseconds: None,
        })
    }

    pub fn range_stats(&self) -> RangeStats {
        self.reader.get_ref().stats()
    }

    pub fn blocks_decoded(&self) -> u64 {
        self.blocks_decoded
    }

    pub(crate) fn enable_timing(&mut self) {
        self.decode_nanoseconds = Some(0);
        self.reader.get_mut().enable_timing();
    }

    pub(crate) fn decompression_nanoseconds(&self) -> Option<u64> {
        Some(
            self.decode_nanoseconds?
                .saturating_sub(self.range_stats().read_nanoseconds?),
        )
    }

    pub fn read_contig_range(
        &mut self,
        contig: Contig<'_>,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, BgzfError> {
        if contig.metagenome_id != self.metagenome_id || start > end || end > contig.length {
            return Err(BgzfError::InvalidRange);
        }
        if start == end {
            return Ok(Vec::new());
        }
        if self.shared_cache.is_some() && self.local_file_changed()? {
            return Err(BgzfError::SourceChanged);
        }
        let line_bases = u64::from(contig.line_bases);
        let line_width = u64::from(contig.line_width);
        if line_bases == 0 || line_width < line_bases || line_width > line_bases.saturating_add(2) {
            return Err(BgzfError::InvalidFasta);
        }
        let uncompressed_offset = contig
            .fasta_offset
            .checked_add(
                (start / line_bases)
                    .checked_mul(line_width)
                    .ok_or(BgzfError::InvalidRange)?,
            )
            .and_then(|offset| offset.checked_add(start % line_bases))
            .ok_or(BgzfError::InvalidRange)?;
        let length = usize::try_from(end - start).map_err(|_| BgzfError::InvalidRange)?;
        let mut sequence = Vec::with_capacity(length);
        let mut position = start;
        let mut file_offset = uncompressed_offset;
        while position < end {
            let count = (end - position).min(line_bases - position % line_bases);
            let old_len = sequence.len();
            sequence.resize(
                old_len
                    .checked_add(usize::try_from(count).map_err(|_| BgzfError::InvalidRange)?)
                    .ok_or(BgzfError::InvalidRange)?,
                0,
            );
            self.read_uncompressed_exact(file_offset, &mut sequence[old_len..])?;
            file_offset = file_offset
                .checked_add(count)
                .ok_or(BgzfError::InvalidRange)?;
            position += count;
            if position < end && position.is_multiple_of(line_bases) {
                let newline_len = usize::try_from(line_width - line_bases)
                    .map_err(|_| BgzfError::InvalidRange)?;
                let mut newline = [0; 2];
                self.read_uncompressed_exact(file_offset, &mut newline[..newline_len])?;
                file_offset = file_offset
                    .checked_add(newline_len as u64)
                    .ok_or(BgzfError::InvalidRange)?;
                if newline[..newline_len]
                    .iter()
                    .any(|byte| !matches!(byte, b'\n' | b'\r'))
                {
                    return Err(BgzfError::InvalidFasta);
                }
            }
        }
        if sequence.iter().any(|base| !is_iupac(*base)) {
            return Err(BgzfError::InvalidFasta);
        }
        sequence.make_ascii_uppercase();
        Ok(sequence)
    }

    fn read_uncompressed_exact(
        &mut self,
        mut uncompressed_offset: u64,
        mut output: &mut [u8],
    ) -> Result<(), BgzfError> {
        while !output.is_empty() {
            let virtual_position = self.index.query(uncompressed_offset)?;
            let compressed_offset = virtual_position.compressed();
            let block_offset = usize::from(virtual_position.uncompressed());
            let block_start = uncompressed_offset
                .checked_sub(block_offset as u64)
                .ok_or(BgzfError::InvalidGzi)?;
            if self.shared_cache.is_some() {
                let count = self.copy_from_shared_block(
                    compressed_offset,
                    block_start,
                    block_offset,
                    output,
                )?;
                uncompressed_offset = uncompressed_offset
                    .checked_add(count as u64)
                    .ok_or(BgzfError::InvalidRange)?;
                output = &mut output[count..];
                continue;
            }
            self.load_block(compressed_offset, block_start)?;
            let block = self.cached_block.as_ref().ok_or(BgzfError::InvalidGzi)?;
            let available = block
                .data
                .get(block_offset..)
                .ok_or(BgzfError::InvalidGzi)?;
            if available.is_empty() {
                return Err(BgzfError::InvalidGzi);
            }
            let count = output.len().min(available.len());
            output[..count].copy_from_slice(&available[..count]);
            uncompressed_offset = uncompressed_offset
                .checked_add(count as u64)
                .ok_or(BgzfError::InvalidRange)?;
            output = &mut output[count..];
        }
        Ok(())
    }

    fn load_block(
        &mut self,
        compressed_offset: u64,
        uncompressed_offset: u64,
    ) -> Result<(), BgzfError> {
        if self.cached_block.as_ref().is_some_and(|block| {
            block.compressed_offset == compressed_offset
                && block.uncompressed_offset == uncompressed_offset
        }) {
            return Ok(());
        }
        let data = self.decode_block(compressed_offset)?;
        self.cached_block = Some(CachedBlock {
            compressed_offset,
            uncompressed_offset,
            data,
        });
        self.blocks_decoded = self.blocks_decoded.saturating_add(1);
        Ok(())
    }

    fn copy_from_shared_block(
        &mut self,
        compressed_offset: u64,
        uncompressed_offset: u64,
        block_offset: usize,
        output: &mut [u8],
    ) -> Result<usize, BgzfError> {
        let cache = Arc::clone(self.shared_cache.as_ref().expect("shared BGZF cache"));
        let source = self.cache_source;
        let identity = BgzfBlockIdentity {
            source_sha256: source.source_sha256,
            source_bytes: source.source_bytes,
            locator_sha256: source.locator_sha256,
            local_file: source.local_file,
            compressed_offset,
            uncompressed_offset,
        };
        cache.with_block(
            identity,
            || {
                let data = self.decode_block(compressed_offset)?;
                self.blocks_decoded = self.blocks_decoded.saturating_add(1);
                Ok(data)
            },
            |data| {
                let available = data.get(block_offset..).ok_or(BgzfError::InvalidGzi)?;
                if available.is_empty() {
                    return Err(BgzfError::InvalidGzi);
                }
                let count = output.len().min(available.len());
                output[..count].copy_from_slice(&available[..count]);
                Ok(count)
            },
        )
    }

    fn decode_block(&mut self, compressed_offset: u64) -> Result<Vec<u8>, BgzfError> {
        let started = self.decode_nanoseconds.map(|_| std::time::Instant::now());
        let virtual_position = bgzf::VirtualPosition::try_from((compressed_offset, 0))
            .map_err(|_| BgzfError::InvalidGzi)?;
        self.reader.seek(virtual_position)?;
        let data = self.reader.fill_buf()?;
        if data.is_empty() || data.len() > MAX_UNCOMPRESSED_BLOCK_BYTES {
            return Err(BgzfError::InvalidGzi);
        }
        let result = data.to_vec();
        if let (Some(started), Some(elapsed)) = (started, &mut self.decode_nanoseconds) {
            *elapsed = elapsed.saturating_add(started.elapsed().as_nanos() as u64);
        }
        Ok(result)
    }

    fn local_file_changed(&self) -> Result<bool, BgzfError> {
        match (&self.local_path, self.cache_source.local_file) {
            (Some(path), Some(expected)) => Ok(file_identity(path)? != expected),
            (Some(_), None) => Ok(true),
            (None, _) => Ok(false),
        }
    }
}

fn local_path(uri: &str) -> Option<PathBuf> {
    if let Some(path) = uri.strip_prefix("file://") {
        return path.starts_with('/').then(|| PathBuf::from(path));
    }
    (!uri.contains("://")).then(|| PathBuf::from(uri))
}

#[cfg(unix)]
fn file_identity(path: &Path) -> io::Result<[u64; 7]> {
    use std::os::unix::fs::MetadataExt;
    let metadata = path.metadata()?;
    Ok([
        metadata.dev(),
        metadata.ino(),
        metadata.len(),
        metadata.mtime() as u64,
        metadata.mtime_nsec() as u64,
        metadata.ctime() as u64,
        metadata.ctime_nsec() as u64,
    ])
}

#[cfg(not(unix))]
fn file_identity(_path: &Path) -> io::Result<[u64; 7]> {
    Err(io::Error::new(
        io::ErrorKind::Unsupported,
        "shared BGZF cache requires local file identity",
    ))
}

fn validate_gzi(index: &gzi::Index, bgzf_bytes: u64) -> Result<(), BgzfError> {
    let mut previous = (0, 0);
    for &(compressed_offset, uncompressed_offset) in index.as_ref() {
        if compressed_offset <= previous.0
            || compressed_offset >= bgzf_bytes
            || bgzf::VirtualPosition::new(compressed_offset, 0).is_none()
            || uncompressed_offset <= previous.1
            || uncompressed_offset - previous.1 > MAX_UNCOMPRESSED_BLOCK_BYTES as u64
        {
            return Err(BgzfError::InvalidGzi);
        }
        previous = (compressed_offset, uncompressed_offset);
    }
    Ok(())
}

fn is_iupac(base: u8) -> bool {
    matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'S'
            | b'W'
            | b'K'
            | b'M'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
    )
}

#[derive(Debug, Error)]
pub enum BgzfError {
    #[error("BGZF I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Range(#[from] RangeSourceError),
    #[error("BGZF resource size differs from JIDX metadata")]
    SizeMismatch,
    #[error("BGZF resource checksum differs from JIDX metadata")]
    ChecksumMismatch,
    #[error("BGZF resource identity changed during cached access")]
    SourceChanged,
    #[error("embedded BGZF index is invalid")]
    InvalidGzi,
    #[error("BGZF contig range is invalid")]
    InvalidRange,
    #[error("BGZF FASTA layout differs from JIDX metadata")]
    InvalidFasta,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx::sha256_reader;
    use crate::jidx_reader::{ContigId, MetagenomeId};
    use std::fs::File;
    use std::io::{BufReader, Write};
    use std::path::Path;

    fn digest(path: &Path) -> [u8; 32] {
        sha256_reader(BufReader::new(File::open(path).unwrap())).unwrap()
    }

    #[test]
    fn reads_checked_ranges_across_fasta_lines() {
        let directory = tempfile::tempdir().unwrap();
        let bgzf_path = directory.path().join("sequence.bgz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(b">ctg\nACGT\nTGCA\nAAAA\n").unwrap();
        writer.finish().unwrap();
        let mut gzi_writer = gzi::io::Writer::new(Vec::new());
        gzi_writer.write_index(&gzi::Index::default()).unwrap();
        let gzi = gzi_writer.into_inner();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let source = Metagenome {
            id: MetagenomeId::default(),
            name: "sample",
            bgzf_uri: &bgzf_uri,
            bgzf_bytes: bgzf_path.metadata().unwrap().len(),
            bgzf_sha256: digest(&bgzf_path),
            contig_start: ContigId::default(),
            contig_count: 1,
            gzi: &gzi,
        };
        let mut reader = BgzfReader::open(source, None, true).unwrap();
        let mut mismatched = source;
        mismatched.bgzf_sha256 = [9; 32];
        assert!(BgzfReader::open(mismatched, None, true).is_err());
        let contig = Contig {
            id: 0,
            metagenome_id: 0,
            name: "ctg",
            length: 12,
            fasta_offset: 5,
            line_bases: 4,
            line_width: 5,
        };
        assert_eq!(
            reader.read_contig_range(contig, 2, 10).unwrap(),
            b"GTTGCAAA"
        );
        assert_eq!(reader.blocks_decoded(), 1);
        assert!(reader.range_stats().bytes_read > 0);
        assert!(reader.read_contig_range(contig, 10, 13).is_err());
    }

    #[test]
    fn reuses_one_block_across_contigs() {
        let directory = tempfile::tempdir().unwrap();
        let bgzf_path = directory.path().join("sequence.bgz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(b">a\nACGT\n>b\nrysw\n").unwrap();
        writer.finish().unwrap();
        let mut gzi_writer = gzi::io::Writer::new(Vec::new());
        gzi_writer.write_index(&gzi::Index::default()).unwrap();
        let gzi = gzi_writer.into_inner();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let source = Metagenome {
            id: 0,
            name: "sample",
            bgzf_uri: &bgzf_uri,
            bgzf_bytes: bgzf_path.metadata().unwrap().len(),
            bgzf_sha256: digest(&bgzf_path),
            contig_start: 0,
            contig_count: 2,
            gzi: &gzi,
        };
        let mut reader = BgzfReader::open(source, None, false).unwrap();
        let first = Contig {
            id: 0,
            metagenome_id: 0,
            name: "a",
            length: 4,
            fasta_offset: 3,
            line_bases: 4,
            line_width: 5,
        };
        let second = Contig {
            id: 1,
            metagenome_id: 0,
            name: "b",
            length: 4,
            fasta_offset: 11,
            line_bases: 4,
            line_width: 5,
        };

        assert_eq!(reader.read_contig_range(first, 0, 4).unwrap(), b"ACGT");
        assert_eq!(reader.read_contig_range(second, 0, 4).unwrap(), b"RYSW");
        assert_eq!(reader.blocks_decoded(), 1);
    }

    #[test]
    fn batch_cache_shares_one_decode_and_preserves_iupac() {
        let directory = tempfile::tempdir().unwrap();
        let bgzf_path = directory.path().join("sequence.bgz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(b">ctg\nacgtryswkmbdhvn\n").unwrap();
        writer.finish().unwrap();
        let mut gzi_writer = gzi::io::Writer::new(Vec::new());
        gzi_writer.write_index(&gzi::Index::default()).unwrap();
        let gzi = gzi_writer.into_inner();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let source = Metagenome {
            id: 0,
            name: "sample",
            bgzf_uri: &bgzf_uri,
            bgzf_bytes: bgzf_path.metadata().unwrap().len(),
            bgzf_sha256: digest(&bgzf_path),
            contig_start: 0,
            contig_count: 1,
            gzi: &gzi,
        };
        let contig = Contig {
            id: 0,
            metagenome_id: 0,
            name: "ctg",
            length: 15,
            fasta_offset: 5,
            line_bases: 15,
            line_width: 16,
        };
        let cache = Arc::new(BgzfBlockCache::new(32 * 1024 * 1024).unwrap());
        let mut first =
            BgzfReader::open_with_cache(source, None, false, Arc::clone(&cache)).unwrap();
        let mut second =
            BgzfReader::open_with_cache(source, None, false, Arc::clone(&cache)).unwrap();

        assert_eq!(
            first.read_contig_range(contig, 0, 15).unwrap(),
            b"ACGTRYSWKMBDHVN"
        );
        assert_eq!(
            second.read_contig_range(contig, 3, 12).unwrap(),
            b"TRYSWKMBD"
        );
        assert_eq!(first.blocks_decoded(), 1);
        assert_eq!(second.blocks_decoded(), 0);
        let stats = cache.stats();
        assert_eq!(stats.blocks_decoded, 1);
        assert_eq!(stats.resident_blocks, 1);
        assert!(stats.accounted_bytes <= 32 * 1024 * 1024);
    }

    #[test]
    fn changed_local_source_cannot_reuse_declared_identity() {
        let directory = tempfile::tempdir().unwrap();
        let bgzf_path = directory.path().join("sequence.bgz");
        let write_sequence = |sequence: &[u8]| {
            let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
            writer.write_all(b">ctg\n").unwrap();
            writer.write_all(sequence).unwrap();
            writer.write_all(b"\n").unwrap();
            writer.finish().unwrap();
        };
        write_sequence(b"AAAA");
        let initial_bytes = bgzf_path.metadata().unwrap().len();
        let mut gzi_writer = gzi::io::Writer::new(Vec::new());
        gzi_writer.write_index(&gzi::Index::default()).unwrap();
        let gzi = gzi_writer.into_inner();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let source = Metagenome {
            id: 0,
            name: "sample",
            bgzf_uri: &bgzf_uri,
            bgzf_bytes: initial_bytes,
            bgzf_sha256: digest(&bgzf_path),
            contig_start: 0,
            contig_count: 1,
            gzi: &gzi,
        };
        let contig = Contig {
            id: 0,
            metagenome_id: 0,
            name: "ctg",
            length: 4,
            fasta_offset: 5,
            line_bases: 4,
            line_width: 5,
        };
        let cache = Arc::new(BgzfBlockCache::new(32 * 1024 * 1024).unwrap());
        let mut stale =
            BgzfReader::open_with_cache(source, None, false, Arc::clone(&cache)).unwrap();
        assert_eq!(stale.read_contig_range(contig, 0, 4).unwrap(), b"AAAA");

        write_sequence(b"CCCC");
        assert_eq!(bgzf_path.metadata().unwrap().len(), initial_bytes);
        assert!(matches!(
            stale.read_contig_range(contig, 0, 4),
            Err(BgzfError::SourceChanged)
        ));
        let mut current =
            BgzfReader::open_with_cache(source, None, false, Arc::clone(&cache)).unwrap();
        assert_eq!(current.read_contig_range(contig, 0, 4).unwrap(), b"CCCC");
        assert_eq!(cache.stats().blocks_decoded, 2);
    }

    #[test]
    fn reads_across_embedded_gzi_block_boundary() {
        let directory = tempfile::tempdir().unwrap();
        let bgzf_path = directory.path().join("sequence.bgz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(b">ctg\nACGT").unwrap();
        writer.flush().unwrap();
        let second_compressed_offset = writer.position();
        writer.write_all(b"TGCA\n").unwrap();
        writer.finish().unwrap();
        let index = gzi::Index::from(vec![(second_compressed_offset, 9)]);
        let mut gzi_writer = gzi::io::Writer::new(Vec::new());
        gzi_writer.write_index(&index).unwrap();
        let gzi = gzi_writer.into_inner();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let source = Metagenome {
            id: 0,
            name: "sample",
            bgzf_uri: &bgzf_uri,
            bgzf_bytes: bgzf_path.metadata().unwrap().len(),
            bgzf_sha256: digest(&bgzf_path),
            contig_start: 0,
            contig_count: 1,
            gzi: &gzi,
        };
        let mut reader = BgzfReader::open(source, None, false).unwrap();
        let contig = Contig {
            id: 0,
            metagenome_id: 0,
            name: "ctg",
            length: 8,
            fasta_offset: 5,
            line_bases: 8,
            line_width: 9,
        };

        assert_eq!(reader.read_contig_range(contig, 2, 8).unwrap(), b"GTTGCA");
        assert_eq!(reader.blocks_decoded(), 2);
        assert_eq!(reader.read_contig_range(contig, 6, 8).unwrap(), b"CA");
        assert_eq!(reader.blocks_decoded(), 2);
    }

    #[test]
    fn rejects_invalid_gzi_offsets() {
        let index = gzi::Index::from(vec![(32, 64), (16, 128)]);
        assert!(matches!(
            validate_gzi(&index, 256),
            Err(BgzfError::InvalidGzi)
        ));
    }
}
