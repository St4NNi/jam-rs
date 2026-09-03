use crate::jidx_reader::{Contig, Metagenome, MetagenomeId};
use crate::range_source::{RangeSource, RangeSourceError, RangeStats, S3Config};
use noodles_bgzf::{self as bgzf, gzi};
use std::io::{self, Read};
use thiserror::Error;

pub struct BgzfReader {
    metagenome_id: MetagenomeId,
    reader: bgzf::io::Reader<RangeSource>,
    index: gzi::Index,
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
            RangeSource::open(source.fai_uri, source.fai_bytes, s3)?,
            RangeSource::open(source.gzi_uri, source.gzi_bytes, s3)?,
            verify,
        )
    }

    fn from_sources(
        source: Metagenome<'_>,
        mut bgzf_source: RangeSource,
        mut fai_source: RangeSource,
        mut gzi_source: RangeSource,
        verify: bool,
    ) -> Result<Self, BgzfError> {
        if bgzf_source.len() != source.bgzf_bytes
            || fai_source.len() != source.fai_bytes
            || gzi_source.len() != source.gzi_bytes
        {
            return Err(BgzfError::SizeMismatch);
        }
        if verify
            && (!bgzf_source.verify_sha256(source.bgzf_sha256)?
                || !fai_source.verify_sha256(source.fai_sha256)?
                || !gzi_source.verify_sha256(source.gzi_sha256)?)
        {
            return Err(BgzfError::ChecksumMismatch);
        }
        let mut index_reader = gzi::io::Reader::new(gzi_source);
        let index = index_reader.read_index()?;
        Ok(Self {
            metagenome_id: source.id,
            reader: bgzf::io::Reader::new(bgzf_source),
            index,
        })
    }

    pub fn range_stats(&self) -> RangeStats {
        self.reader.get_ref().stats()
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
        let line_bases = u64::from(contig.line_bases);
        let line_width = u64::from(contig.line_width);
        let uncompressed_offset = contig
            .fasta_offset
            .checked_add(
                (start / line_bases)
                    .checked_mul(line_width)
                    .ok_or(BgzfError::InvalidRange)?,
            )
            .and_then(|offset| offset.checked_add(start % line_bases))
            .ok_or(BgzfError::InvalidRange)?;
        self.reader
            .seek_by_uncompressed_position(&self.index, uncompressed_offset)?;

        let length = usize::try_from(end - start).map_err(|_| BgzfError::InvalidRange)?;
        let mut sequence = Vec::with_capacity(length);
        let mut position = start;
        while position < end {
            let count = (end - position).min(line_bases - position % line_bases);
            let old_len = sequence.len();
            sequence.resize(
                old_len
                    .checked_add(usize::try_from(count).map_err(|_| BgzfError::InvalidRange)?)
                    .ok_or(BgzfError::InvalidRange)?,
                0,
            );
            self.reader.read_exact(&mut sequence[old_len..])?;
            position += count;
            if position < end && position.is_multiple_of(line_bases) {
                let newline_len = usize::try_from(line_width - line_bases)
                    .map_err(|_| BgzfError::InvalidRange)?;
                let mut newline = [0; 2];
                self.reader.read_exact(&mut newline[..newline_len])?;
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
        let fai_path = directory.path().join("sequence.bgz.fai");
        let gzi_path = directory.path().join("sequence.bgz.gzi");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(b">ctg\nACGT\nTGCA\nAAAA\n").unwrap();
        writer.finish().unwrap();
        std::fs::write(&fai_path, b"ctg\t12\t5\t4\t5\n").unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let bgzf_uri = bgzf_path.to_string_lossy();
        let fai_uri = fai_path.to_string_lossy();
        let gzi_uri = gzi_path.to_string_lossy();
        let source = Metagenome {
            id: MetagenomeId::default(),
            name: "sample",
            bgzf_uri: &bgzf_uri,
            fai_uri: &fai_uri,
            gzi_uri: &gzi_uri,
            bgzf_bytes: bgzf_path.metadata().unwrap().len(),
            fai_bytes: fai_path.metadata().unwrap().len(),
            gzi_bytes: gzi_path.metadata().unwrap().len(),
            bgzf_sha256: digest(&bgzf_path),
            fai_sha256: digest(&fai_path),
            gzi_sha256: digest(&gzi_path),
            contig_start: ContigId::default(),
            contig_count: 1,
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
        assert!(reader.range_stats().bytes_read > 0);
        assert!(reader.read_contig_range(contig, 10, 13).is_err());
    }
}
