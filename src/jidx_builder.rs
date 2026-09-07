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
    pub minimizer_window: u16,
}

impl Default for JidxBuildConfig {
    fn default() -> Self {
        Self {
            k: 21,
            minimizer_window: 16,
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
        let fai = parse_fai(&std::fs::read(&fai_path)?)?;
        let gzi = std::fs::read(&gzi_path)?;
        let bgzf_uri = path_text(&bgzf_path)?;
        let source = Metagenome {
            id: 0,
            name: &entry.name,
            bgzf_uri: &bgzf_uri,
            bgzf_bytes,
            bgzf_sha256,
            contig_start: 0,
            contig_count: u32::try_from(fai.len())
                .map_err(|_| JidxBuildError::Invalid("contig count"))?,
            gzi: &gzi,
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
            bgzf_bytes,
            bgzf_sha256,
            gzi,
            contigs,
        });
    }
    let written = write_jidx(
        output,
        &JidxInput {
            k: config.k,
            minimizer_window: config.minimizer_window,
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
    if !(1..=32).contains(&config.k) || config.minimizer_window == 0 {
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
    let mut previous_position = None;
    for (position, kmer, orientation) in normalized.bit_kmers(config.k, true) {
        let position =
            u64::try_from(position).map_err(|_| JidxBuildError::Invalid("seed position"))?;
        if previous_position.is_some_and(|previous| position != previous + 1) {
            append_minimizers(&mut output, &mut candidates, config.minimizer_window);
        }
        previous_position = Some(position);
        candidates.push((
            jamhash_u64(kmer.0),
            SelectedSeed {
                packed_key: kmer.0,
                position,
                canonical_orientation: orientation,
            },
        ));
    }
    append_minimizers(&mut output, &mut candidates, config.minimizer_window);
    Ok(output)
}

fn append_minimizers(
    output: &mut Vec<SelectedSeed>,
    candidates: &mut Vec<(u64, SelectedSeed)>,
    window: u16,
) {
    if candidates.is_empty() {
        return;
    }
    let window = usize::from(window).min(candidates.len());
    let mut selected = vec![false; candidates.len()];
    for (offset, slice) in candidates.windows(window).enumerate() {
        let minimum = slice
            .iter()
            .map(|(hash, seed)| (*hash, seed.packed_key))
            .min()
            .expect("a minimizer window is nonempty");
        for (index, (hash, seed)) in slice.iter().enumerate() {
            if (*hash, seed.packed_key) == minimum {
                selected[offset + index] = true;
            }
        }
    }
    output.extend(
        candidates
            .iter()
            .zip(selected)
            .filter_map(|((_, seed), selected)| selected.then_some(*seed)),
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
                minimizer_window: 16,
            },
        )
        .unwrap();
        assert_eq!(stats.source_bases, sequence.len() as u64);
        assert!(stats.written.seeds > 0);
        let index = JidxReader::open(output).unwrap();
        index.verify_checksum().unwrap();
        std::fs::remove_file(fai_path).unwrap();
        std::fs::remove_file(gzi_path).unwrap();
        let metagenome = index.metagenome(0).unwrap().unwrap();
        let contig = index.contig(0).unwrap().unwrap();
        let mut reader = BgzfReader::open(metagenome, None, false).unwrap();
        assert_eq!(
            reader.read_contig_range(contig, 3, 17).unwrap(),
            sequence.as_bytes()[3..17]
        );
    }

    #[test]
    fn sliding_windows_are_covered_deterministically() {
        let sequence = b"ACGTTGCAACGATCGTACGTTGCAACGATCGT";
        let config = JidxBuildConfig {
            k: 5,
            minimizer_window: 4,
        };
        let first = select_seeds(sequence, config).unwrap();
        assert_eq!(select_seeds(sequence, config).unwrap(), first);
        assert!(first.iter().all(|seed| seed.packed_key < 1 << 10));

        let candidates: Vec<_> = sequence
            .bit_kmers(config.k, true)
            .map(|(position, kmer, _)| (position, jamhash_u64(kmer.0), kmer.0))
            .collect();
        for window in candidates.windows(usize::from(config.minimizer_window)) {
            let minimum = window
                .iter()
                .map(|(_, hash, key)| (*hash, *key))
                .min()
                .unwrap();
            assert!(window.iter().any(|(position, hash, key)| {
                (*hash, *key) == minimum
                    && first
                        .iter()
                        .any(|selected| selected.position == *position as u64)
            }));
        }
    }

    #[test]
    fn tied_minima_and_short_runs_are_retained() {
        let config = JidxBuildConfig {
            k: 3,
            minimizer_window: 4,
        };
        let tied = select_seeds(b"AAAAAA", config).unwrap();
        assert_eq!(
            tied.iter().map(|seed| seed.position).collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );

        let split = select_seeds(b"AAAAANCCCCC", config).unwrap();
        assert_eq!(
            split.iter().map(|seed| seed.position).collect::<Vec<_>>(),
            vec![0, 1, 2, 6, 7, 8]
        );
        assert_eq!(select_seeds(b"ACG", config).unwrap()[0].position, 0);
        assert!(select_seeds(b"AC", config).unwrap().is_empty());
    }

    #[test]
    fn reverse_complement_selection_is_equivalent() {
        let sequence = b"ACGTTGCAACGATCGTAGGCTAACCGT";
        let config = JidxBuildConfig {
            k: 5,
            minimizer_window: 4,
        };
        let forward = select_seeds(sequence, config).unwrap();
        let reverse = select_seeds(&sequence.reverse_complement(), config).unwrap();
        let last_start = sequence.len() as u64 - u64::from(config.k);
        let mirrored: BTreeSet<_> = forward
            .iter()
            .map(|seed| {
                (
                    seed.packed_key,
                    last_start - seed.position,
                    !seed.canonical_orientation,
                )
            })
            .collect();
        assert_eq!(
            reverse
                .iter()
                .map(|seed| (seed.packed_key, seed.position, seed.canonical_orientation))
                .collect::<BTreeSet<_>>(),
            mirrored
        );
    }

    #[test]
    fn zero_minimizer_window_is_rejected() {
        assert!(
            validate_config(JidxBuildConfig {
                k: 21,
                minimizer_window: 0,
            })
            .is_err()
        );
    }
}
