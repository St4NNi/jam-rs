use crate::bgzf::BgzfReader;
use crate::jidx::{sha256, sha256_reader};
use crate::jidx_reader::{Contig, Metagenome};
use crate::jidx_writer::{
    ContigInput, JidxInput, JidxWriteError, JidxWriteStats, MetagenomeInput, SelectedSeed,
    write_jidx,
};
use crate::reader::{JamReader, ReaderError};
use jamhash::jamhash_u64;
use needletail::Sequence;
use serde::Deserialize;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JidxBuildConfig {
    pub k: u8,
    pub segment_bases: u32,
    pub seeds_per_segment: u16,
}

impl Default for JidxBuildConfig {
    fn default() -> Self {
        Self {
            k: 21,
            segment_bases: 256,
            seeds_per_segment: 2,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JidxBuildStats {
    pub source_bases: u64,
    pub written: JidxWriteStats,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    metagenomes: Vec<ManifestMetagenome>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct ManifestMetagenome {
    name: String,
    bgzf: String,
    fai: String,
    gzi: String,
}

pub fn build_local_jidx(
    jam_path: impl AsRef<Path>,
    manifest_path: impl AsRef<Path>,
    output: impl AsRef<Path>,
    config: JidxBuildConfig,
) -> Result<JidxBuildStats, JidxBuildError> {
    validate_config(config)?;
    let jam_path = jam_path.as_ref();
    let manifest_path = manifest_path.as_ref();
    let manifest_bytes = std::fs::read(manifest_path)?;
    let manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
    let database = JamReader::open(jam_path)?;
    validate_names(&database, &manifest)?;
    let jam_sha256 = file_digest(jam_path)?.1;
    let manifest_sha256 = sha256(&manifest_bytes);
    let base = manifest_path.parent().unwrap_or_else(|| Path::new("."));

    let mut metagenomes = Vec::with_capacity(manifest.metagenomes.len());
    let mut source_bases = 0u64;
    for entry in manifest.metagenomes {
        let bgzf_path = resolve_local(base, &entry.bgzf)?;
        let fai_path = resolve_local(base, &entry.fai)?;
        let gzi_path = resolve_local(base, &entry.gzi)?;
        let (bgzf_bytes, bgzf_sha256) = file_digest(&bgzf_path)?;
        let (fai_bytes, fai_sha256) = file_digest(&fai_path)?;
        let (gzi_bytes, gzi_sha256) = file_digest(&gzi_path)?;
        let fai = parse_fai(&std::fs::read(&fai_path)?)?;
        let bgzf_uri = path_text(&bgzf_path)?;
        let fai_uri = path_text(&fai_path)?;
        let gzi_uri = path_text(&gzi_path)?;
        let source = Metagenome {
            id: 0,
            name: &entry.name,
            bgzf_uri: &bgzf_uri,
            fai_uri: &fai_uri,
            gzi_uri: &gzi_uri,
            bgzf_bytes,
            fai_bytes,
            gzi_bytes,
            bgzf_sha256,
            fai_sha256,
            gzi_sha256,
            contig_start: 0,
            contig_count: u32::try_from(fai.len())
                .map_err(|_| JidxBuildError::Invalid("contig count"))?,
        };
        let mut reader = BgzfReader::open(source, None, false)?;
        let mut contigs = Vec::with_capacity(fai.len());
        for (id, record) in fai.into_iter().enumerate() {
            let contig = Contig {
                id: u32::try_from(id).map_err(|_| JidxBuildError::Invalid("contig count"))?,
                metagenome_id: 0,
                name: &record.name,
                length: record.length,
                fasta_offset: record.offset,
                line_bases: record.line_bases,
                line_width: record.line_width,
            };
            let sequence = reader.read_contig_range(contig, 0, record.length)?;
            let seeds = select_seeds(&sequence, config)?;
            source_bases = source_bases
                .checked_add(record.length)
                .ok_or(JidxBuildError::Invalid("source bases"))?;
            contigs.push(ContigInput {
                name: record.name,
                length: record.length,
                fasta_offset: record.offset,
                line_bases: record.line_bases,
                line_width: record.line_width,
                seeds,
            });
        }
        metagenomes.push(MetagenomeInput {
            name: entry.name,
            bgzf_uri,
            fai_uri,
            gzi_uri,
            bgzf_bytes,
            fai_bytes,
            gzi_bytes,
            bgzf_sha256,
            fai_sha256,
            gzi_sha256,
            contigs,
        });
    }
    let written = write_jidx(
        output,
        &JidxInput {
            k: config.k,
            segment_bases: config.segment_bases,
            seeds_per_segment: config.seeds_per_segment,
            jam_sha256,
            manifest_sha256,
            metagenomes,
        },
    )?;
    Ok(JidxBuildStats {
        source_bases,
        written,
    })
}

fn validate_config(config: JidxBuildConfig) -> Result<(), JidxBuildError> {
    if !(1..=32).contains(&config.k)
        || config.segment_bases < u32::from(config.k)
        || config.seeds_per_segment == 0
        || u32::from(config.seeds_per_segment) > config.segment_bases - u32::from(config.k) + 1
    {
        return Err(JidxBuildError::Invalid("seed selection"));
    }
    Ok(())
}

fn validate_names(database: &JamReader, manifest: &Manifest) -> Result<(), JidxBuildError> {
    let database_names: BTreeSet<_> = database.sample_names().iter().map(String::as_str).collect();
    let manifest_names: BTreeSet<_> = manifest
        .metagenomes
        .iter()
        .map(|entry| entry.name.as_str())
        .collect();
    if database_names.len() != database.sample_names().len()
        || manifest_names.len() != manifest.metagenomes.len()
        || database_names != manifest_names
    {
        return Err(JidxBuildError::Invalid("manifest metagenome names"));
    }
    Ok(())
}

fn select_seeds(
    sequence: &[u8],
    config: JidxBuildConfig,
) -> Result<Vec<SelectedSeed>, JidxBuildError> {
    let normalized = sequence.normalize(false);
    let mut output = Vec::new();
    let mut candidates = Vec::new();
    let mut current_segment = None;
    for (position, kmer, orientation) in normalized.bit_kmers(config.k, true) {
        let position =
            u64::try_from(position).map_err(|_| JidxBuildError::Invalid("seed position"))?;
        let segment = position / u64::from(config.segment_bases);
        if current_segment.is_some_and(|current| current != segment) {
            append_minima(&mut output, &mut candidates, config.seeds_per_segment);
        }
        current_segment = Some(segment);
        candidates.push((
            jamhash_u64(kmer.0),
            SelectedSeed {
                packed_key: kmer.0,
                position,
                canonical_orientation: orientation,
            },
        ));
    }
    append_minima(&mut output, &mut candidates, config.seeds_per_segment);
    output.sort_unstable_by_key(|seed| seed.position);
    Ok(output)
}

fn append_minima(
    output: &mut Vec<SelectedSeed>,
    candidates: &mut Vec<(u64, SelectedSeed)>,
    count: u16,
) {
    candidates.sort_unstable_by(|left, right| {
        (
            left.0,
            left.1.packed_key,
            left.1.position,
            left.1.canonical_orientation,
        )
            .cmp(&(
                right.0,
                right.1.packed_key,
                right.1.position,
                right.1.canonical_orientation,
            ))
    });
    output.extend(
        candidates
            .iter()
            .take(usize::from(count))
            .map(|(_, seed)| *seed),
    );
    candidates.clear();
}

struct FaiRecord {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u32,
    line_width: u32,
}

fn parse_fai(bytes: &[u8]) -> Result<Vec<FaiRecord>, JidxBuildError> {
    let text = std::str::from_utf8(bytes).map_err(|_| JidxBuildError::Invalid("FAI encoding"))?;
    let mut names = BTreeSet::new();
    let mut records = Vec::new();
    for line in text.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5 || fields[0].is_empty() || !names.insert(fields[0]) {
            return Err(JidxBuildError::Invalid("FAI record"));
        }
        let record = FaiRecord {
            name: fields[0].to_string(),
            length: fields[1]
                .parse()
                .map_err(|_| JidxBuildError::Invalid("FAI length"))?,
            offset: fields[2]
                .parse()
                .map_err(|_| JidxBuildError::Invalid("FAI offset"))?,
            line_bases: fields[3]
                .parse()
                .map_err(|_| JidxBuildError::Invalid("FAI line bases"))?,
            line_width: fields[4]
                .parse()
                .map_err(|_| JidxBuildError::Invalid("FAI line width"))?,
        };
        if record.length == 0
            || record.line_bases == 0
            || record.line_width < record.line_bases
            || record.line_width > record.line_bases.saturating_add(2)
        {
            return Err(JidxBuildError::Invalid("FAI record"));
        }
        records.push(record);
    }
    if records.is_empty() {
        return Err(JidxBuildError::Invalid("FAI records"));
    }
    Ok(records)
}

fn resolve_local(base: &Path, value: &str) -> Result<PathBuf, JidxBuildError> {
    let path = if let Some(path) = value.strip_prefix("file://") {
        if !path.starts_with('/') {
            return Err(JidxBuildError::Invalid("local URI"));
        }
        PathBuf::from(path)
    } else {
        if value.contains("://") {
            return Err(JidxBuildError::Invalid("local URI"));
        }
        let path = PathBuf::from(value);
        if path.is_absolute() {
            path
        } else {
            base.join(path)
        }
    };
    Ok(path.canonicalize()?)
}

fn path_text(path: &Path) -> Result<String, JidxBuildError> {
    path.to_str()
        .map(str::to_owned)
        .ok_or(JidxBuildError::Invalid("non-UTF-8 resource path"))
}

fn file_digest(path: &Path) -> Result<(u64, [u8; 32]), JidxBuildError> {
    let file = File::open(path)?;
    let bytes = file.metadata()?.len();
    Ok((bytes, sha256_reader(BufReader::new(file))?))
}

#[derive(Debug, Error)]
pub enum JidxBuildError {
    #[error("JIDX build I/O failed: {0}")]
    Io(#[from] io::Error),
    #[error("JIDX manifest is invalid: {0}")]
    Json(#[from] serde_json::Error),
    #[error(transparent)]
    Jam(#[from] ReaderError),
    #[error(transparent)]
    Bgzf(#[from] crate::bgzf::BgzfError),
    #[error(transparent)]
    Write(#[from] JidxWriteError),
    #[error("invalid JIDX build input: {0}")]
    Invalid(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx_reader::JidxReader;
    use crate::writer::{BuildConfig, build};
    use noodles_bgzf::{self as bgzf, gzi};
    use std::io::Write;

    #[test]
    fn builds_query_independent_index_from_local_bgzf() {
        let directory = tempfile::tempdir().unwrap();
        let source = directory.path().join("sample.fa");
        let sequence = "ACGTTGCAACGATCGTACGTTGCAACGATCGTACGTTGCAACGATCGTACGTTGCAACGATCGT";
        std::fs::write(&source, format!(">sample\n{sequence}\n")).unwrap();
        let jam = directory.path().join("database.jam");
        build(
            &[source],
            &jam,
            &BuildConfig {
                kmer_size: 5,
                fscale: 1,
                singleton: true,
                memory: 1,
                ..BuildConfig::default()
            },
        )
        .unwrap();

        let bgzf_path = directory.path().join("sample.bgz");
        let fai_path = directory.path().join("sample.bgz.fai");
        let gzi_path = directory.path().join("sample.bgz.gzi");
        let mut raw = b">contig\n".to_vec();
        for line in sequence.as_bytes().chunks(16) {
            raw.extend_from_slice(line);
            raw.push(b'\n');
        }
        let mut writer = bgzf::io::Writer::new(File::create(&bgzf_path).unwrap());
        writer.write_all(&raw).unwrap();
        writer.finish().unwrap();
        std::fs::write(
            &fai_path,
            format!("contig\t{}\t8\t16\t17\n", sequence.len()),
        )
        .unwrap();
        gzi::fs::write(&gzi_path, &gzi::Index::default()).unwrap();
        let manifest = directory.path().join("manifest.json");
        std::fs::write(
            &manifest,
            format!(
                "{{\"metagenomes\":[{{\"name\":\"sample\",\"bgzf\":\"{}\",\"fai\":\"{}\",\"gzi\":\"{}\"}}]}}",
                bgzf_path.display(),
                fai_path.display(),
                gzi_path.display()
            ),
        )
        .unwrap();
        let output = directory.path().join("index.jidx");
        let stats = build_local_jidx(
            &jam,
            &manifest,
            &output,
            JidxBuildConfig {
                k: 5,
                segment_bases: 16,
                seeds_per_segment: 2,
            },
        )
        .unwrap();
        assert_eq!(stats.source_bases, sequence.len() as u64);
        assert!(stats.written.seeds > 0);
        assert!(stats.written.occurrences <= 10);
        JidxReader::open(output).unwrap().verify_checksum().unwrap();
    }

    #[test]
    fn window_selection_is_bounded_and_deterministic() {
        let sequence = b"ACGTTGCAACGATCGTACGTTGCAACGATCGT";
        let config = JidxBuildConfig {
            k: 5,
            segment_bases: 16,
            seeds_per_segment: 2,
        };
        let first = select_seeds(sequence, config).unwrap();
        assert_eq!(select_seeds(sequence, config).unwrap(), first);
        assert!(first.len() <= 4);
        assert!(first.iter().all(|seed| seed.packed_key < 1 << 10));
    }
}
