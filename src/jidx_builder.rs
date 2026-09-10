use crate::bgzf::BgzfReader;
use crate::jidx::{RESCUE_K15_TAG, sha256, sha256_reader};
use crate::jidx_reader::{Contig, Metagenome};
use crate::jidx_writer::{
    ContigInput, JidxInput, JidxWriteError, JidxWriteStats, JidxWriter, MetagenomeInput,
    SelectedSeed,
};
use crate::reader::{JamReader, ReaderError};
use jamhash::jamhash_u64;
use needletail::Sequence;
use serde::Deserialize;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use thiserror::Error;

const SEED_CORE_BYTES: u64 = 1024 * 1024;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JidxBuildConfig {
    pub k: u8,
    pub minimizer_window: u16,
    pub rescue_k15: bool,
}

impl Default for JidxBuildConfig {
    fn default() -> Self {
        Self {
            k: 21,
            minimizer_window: 16,
            rescue_k15: false,
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
    let mut manifest: Manifest = serde_json::from_slice(&manifest_bytes)?;
    let database = JamReader::open(jam_path)?;
    validate_names(&database, &manifest)?;
    let jam_sha256 = file_digest(jam_path)?.1;
    let manifest_sha256 = sha256(&manifest_bytes);
    let base = manifest_path.parent().unwrap_or_else(|| Path::new("."));
    manifest
        .metagenomes
        .sort_unstable_by(|left, right| left.name.cmp(&right.name));

    let mut writer = JidxWriter::new(
        output,
        &JidxInput {
            k: config.k,
            minimizer_window: config.minimizer_window,
            rescue_k15: config.rescue_k15,
            jam_sha256,
            manifest_sha256,
        },
    )?;
    let mut source_bases = 0u64;
    for entry in manifest.metagenomes {
        let bgzf_path = resolve_local(base, &entry.bgzf)?;
        let fai_path = resolve_local(base, &entry.fai)?;
        let gzi_path = resolve_local(base, &entry.gzi)?;
        let (bgzf_bytes, bgzf_sha256) = file_digest(&bgzf_path)?;
        let contig_count = visit_fai(&fai_path, true, |_| Ok(()))?;
        let gzi = std::fs::read(&gzi_path)?;
        let bgzf_uri = path_text(&bgzf_path)?;
        let source = Metagenome {
            id: 0,
            name: &entry.name,
            bgzf_uri: &bgzf_uri,
            bgzf_bytes,
            bgzf_sha256,
            contig_start: 0,
            contig_count,
            gzi: &gzi,
        };
        let mut reader = BgzfReader::open(source, None, false)?;
        writer.begin_metagenome(MetagenomeInput {
            name: entry.name,
            bgzf_uri,
            bgzf_bytes,
            bgzf_sha256,
            gzi,
        })?;
        let streamed_count = visit_fai(&fai_path, false, |record| {
            let contig_id = writer.begin_contig(ContigInput {
                name: record.name.clone(),
                length: record.length,
                fasta_offset: record.offset,
                line_bases: record.line_bases,
                line_width: record.line_width,
            })?;
            write_contig_seeds(&mut writer, &mut reader, contig_id, &record, config)?;
            source_bases = source_bases
                .checked_add(record.length)
                .ok_or(JidxBuildError::Invalid("source bases"))?;
            Ok(())
        })?;
        if streamed_count != contig_count {
            return Err(JidxBuildError::Invalid("contig count"));
        }
    }
    let written = writer.finish()?;
    Ok(JidxBuildStats {
        source_bases,
        written,
    })
}

fn validate_config(config: JidxBuildConfig) -> Result<(), JidxBuildError> {
    if !(1..=32).contains(&config.k)
        || config.minimizer_window == 0
        || (config.rescue_k15 && config.k != 21)
    {
        return Err(JidxBuildError::Invalid("seed selection"));
    }
    Ok(())
}

fn write_contig_seeds(
    writer: &mut JidxWriter,
    reader: &mut BgzfReader,
    contig_id: u32,
    record: &FaiRecord,
    config: JidxBuildConfig,
) -> Result<(), JidxBuildError> {
    let overlap = seed_overlap(config)?;
    let mut core_start = 0u64;
    while core_start < record.length {
        let core_end = core_start
            .saturating_add(SEED_CORE_BYTES)
            .min(record.length);
        let read_start = core_start.saturating_sub(overlap);
        let read_end = core_end
            .checked_add(overlap)
            .ok_or(JidxBuildError::Invalid("seed chunk range"))?
            .min(record.length);
        let sequence = reader.read_contig_range(
            Contig {
                id: contig_id,
                metagenome_id: 0,
                name: &record.name,
                length: record.length,
                fasta_offset: record.offset,
                line_bases: record.line_bases,
                line_width: record.line_width,
            },
            read_start,
            read_end,
        )?;
        let seeds = select_core_seeds(&sequence, read_start, core_start, core_end, config)?;
        if !seeds.is_empty() {
            writer.add_seeds(contig_id, &seeds)?;
        }
        core_start = core_end;
    }
    Ok(())
}

fn seed_overlap(config: JidxBuildConfig) -> Result<u64, JidxBuildError> {
    u64::from(config.k)
        .checked_add(u64::from(config.minimizer_window))
        .and_then(|value| value.checked_sub(2))
        .ok_or(JidxBuildError::Invalid("seed chunk overlap"))
}

fn select_core_seeds(
    sequence: &[u8],
    read_start: u64,
    core_start: u64,
    core_end: u64,
    config: JidxBuildConfig,
) -> Result<Vec<SelectedSeed>, JidxBuildError> {
    let mut seeds = select_index_seeds(sequence, config)?;
    for seed in &mut seeds {
        seed.position = seed
            .position
            .checked_add(read_start)
            .ok_or(JidxBuildError::Invalid("seed position"))?;
    }
    seeds.retain(|seed| (core_start..core_end).contains(&seed.position));
    Ok(seeds)
}

fn select_index_seeds(
    sequence: &[u8],
    config: JidxBuildConfig,
) -> Result<Vec<SelectedSeed>, JidxBuildError> {
    let mut seeds = select_seeds(sequence, config)?;
    if config.rescue_k15 {
        let mut rescue = select_seeds(
            sequence,
            JidxBuildConfig {
                k: 15,
                minimizer_window: config.minimizer_window,
                rescue_k15: false,
            },
        )?;
        for seed in &mut rescue {
            seed.packed_key |= RESCUE_K15_TAG;
        }
        seeds.extend(rescue);
        seeds.sort_unstable_by_key(|seed| (seed.position, seed.packed_key));
    }
    Ok(seeds)
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

pub(crate) fn select_seeds(
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

fn visit_fai(
    path: &Path,
    validate_unique_names: bool,
    mut visit: impl FnMut(FaiRecord) -> Result<(), JidxBuildError>,
) -> Result<u32, JidxBuildError> {
    let mut names = validate_unique_names.then(BTreeSet::new);
    let mut count = 0u32;
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 5
            || fields[0].is_empty()
            || names
                .as_mut()
                .is_some_and(|names| !names.insert(fields[0].to_owned()))
        {
            return Err(JidxBuildError::Invalid("FAI record"));
        }
        let record = FaiRecord {
            name: fields[0].to_owned(),
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
        visit(record)?;
        count = count
            .checked_add(1)
            .ok_or(JidxBuildError::Invalid("contig count"))?;
    }
    if count == 0 {
        return Err(JidxBuildError::Invalid("FAI records"));
    }
    Ok(count)
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

    fn collect_chunked_seeds(
        sequence: &[u8],
        config: JidxBuildConfig,
        core_bytes: u64,
    ) -> Vec<SelectedSeed> {
        let length = sequence.len() as u64;
        let overlap = seed_overlap(config).unwrap();
        let mut output = Vec::new();
        let mut core_start = 0u64;
        while core_start < length {
            let core_end = core_start.saturating_add(core_bytes).min(length);
            let read_start = core_start.saturating_sub(overlap);
            let read_end = core_end.saturating_add(overlap).min(length);
            output.extend(
                select_core_seeds(
                    &sequence[read_start as usize..read_end as usize],
                    read_start,
                    core_start,
                    core_end,
                    config,
                )
                .unwrap(),
            );
            core_start = core_end;
        }
        output
    }

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
                rescue_k15: false,
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
            rescue_k15: false,
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
    fn chunked_selection_matches_whole_contigs() {
        let mut sequence = b"ACGT".repeat(64);
        sequence[45] = b'N';
        sequence[70] = b'N';
        sequence[91] = b'N';
        sequence[115..170].fill(b'A');
        for minimizer_window in [16, 32, 64] {
            let config = JidxBuildConfig {
                k: 21,
                minimizer_window,
                rescue_k15: true,
            };
            for sequence in [sequence.clone(), sequence.reverse_complement()] {
                assert_eq!(
                    collect_chunked_seeds(&sequence, config, 23),
                    select_index_seeds(&sequence, config).unwrap()
                );
            }
        }
    }

    #[test]
    fn tied_minima_and_short_runs_are_retained() {
        let config = JidxBuildConfig {
            k: 3,
            minimizer_window: 4,
            rescue_k15: false,
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
    fn rescue_keys_preserve_reverse_complement_positions() {
        let sequence = b"ACGTTGCAACGATCGTAGGCTAACCGTAGCTACGATTCGA";
        let config = JidxBuildConfig {
            k: 21,
            minimizer_window: 4,
            rescue_k15: true,
        };
        let forward = select_index_seeds(sequence, config).unwrap();
        let reverse = select_index_seeds(&sequence.reverse_complement(), config).unwrap();
        assert!(
            forward
                .iter()
                .any(|seed| seed.packed_key & RESCUE_K15_TAG == 0)
        );
        assert!(forward.iter().any(|seed| {
            seed.packed_key & RESCUE_K15_TAG != 0 && (seed.packed_key & !RESCUE_K15_TAG) < 1 << 30
        }));
        let mirrored: BTreeSet<_> = forward
            .iter()
            .map(|seed| {
                let k =
                    crate::jidx::seed_length(config.k, config.rescue_k15, seed.packed_key).unwrap();
                (
                    seed.packed_key,
                    sequence.len() as u64 - u64::from(k) - seed.position,
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
    fn invalid_selection_configs_are_rejected() {
        assert!(
            validate_config(JidxBuildConfig {
                k: 21,
                minimizer_window: 0,
                rescue_k15: false,
            })
            .is_err()
        );
        assert!(
            validate_config(JidxBuildConfig {
                k: 20,
                minimizer_window: 16,
                rescue_k15: true,
            })
            .is_err()
        );
        assert!(
            validate_config(JidxBuildConfig {
                k: 21,
                minimizer_window: 16,
                rescue_k15: true,
            })
            .is_ok()
        );
    }
}
