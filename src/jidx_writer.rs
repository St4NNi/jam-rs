use crate::jidx::{
    CONTIG_POSTING_SIZE, CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_POSTING_SIZE,
    DOCUMENT_RECORD_SIZE, DocumentRecord, FilterKind, HEADER_SIZE, Header, JidxError, PostingCodec,
    SECTION_COUNT, SEED_RECORD_SIZE, SectionDescriptor, SectionKind, SeedScheme, StringRef,
    put_u32, sha256_reader,
};
use crate::jidx_postings::{
    SeedEntry, SeedOccurrence, encode_occurrence, encode_seed, validate_packed_key,
};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Seek, SeekFrom, Write};
use std::path::Path;
use thiserror::Error;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SelectedSeed {
    pub packed_key: u64,
    pub position: u64,
    pub canonical_orientation: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ContigInput {
    pub name: String,
    pub length: u64,
    pub fasta_offset: u64,
    pub line_bases: u32,
    pub line_width: u32,
    pub seeds: Vec<SelectedSeed>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeInput {
    pub name: String,
    pub bgzf_uri: String,
    pub fai_uri: String,
    pub gzi_uri: String,
    pub bgzf_bytes: u64,
    pub fai_bytes: u64,
    pub gzi_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub fai_sha256: [u8; 32],
    pub gzi_sha256: [u8; 32],
    pub contigs: Vec<ContigInput>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct JidxInput {
    pub k: u8,
    pub segment_bases: u32,
    pub seeds_per_segment: u16,
    pub jam_sha256: [u8; 32],
    pub manifest_sha256: [u8; 32],
    pub metagenomes: Vec<MetagenomeInput>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct JidxWriteStats {
    pub metagenomes: u32,
    pub contigs: u32,
    pub seeds: u64,
    pub occurrences: u64,
    pub file_bytes: u64,
    pub file_sha256: [u8; 32],
}

pub fn write_jidx(
    path: impl AsRef<Path>,
    input: &JidxInput,
) -> Result<JidxWriteStats, JidxWriteError> {
    let path = path.as_ref();
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if path.file_name().is_none() {
        return Err(JidxWriteError::Invalid("output path"));
    }
    if path.try_exists()? {
        return Err(JidxWriteError::Invalid("output exists"));
    }
    let prepared = prepare(input)?;
    let mut temporary = tempfile::Builder::new()
        .prefix(".jidx-")
        .tempfile_in(parent)?;
    let file = temporary.as_file_mut();
    file.write_all(&[0; HEADER_SIZE])?;

    write_padding(file, prepared.sections[0].offset)?;
    file.write_all(&prepared.strings)?;
    write_padding(file, prepared.sections[1].offset)?;
    for record in prepared.documents {
        file.write_all(&record.encode())?;
    }
    write_padding(file, prepared.sections[2].offset)?;
    for record in prepared.contigs {
        file.write_all(&record.encode())?;
    }
    write_padding(file, prepared.sections[3].offset)?;
    let mut document_offset = 0;
    let mut occurrence_offset = 0;
    for (key, occurrences) in &prepared.seeds {
        let documents = documents_for(occurrences, &prepared.contig_documents);
        file.write_all(&encode_seed(SeedEntry {
            packed_key: *key,
            document_frequency: u32::try_from(documents.len())
                .map_err(|_| JidxWriteError::Invalid("document frequency"))?,
            occurrence_count: u64::try_from(occurrences.len())
                .map_err(|_| JidxWriteError::Invalid("occurrence count"))?,
            document_offset,
            occurrence_offset,
        }))?;
        document_offset = document_offset
            .checked_add(byte_len(documents.len(), DOCUMENT_POSTING_SIZE)?)
            .ok_or(JidxWriteError::Invalid("document postings"))?;
        occurrence_offset = occurrence_offset
            .checked_add(byte_len(occurrences.len(), CONTIG_POSTING_SIZE)?)
            .ok_or(JidxWriteError::Invalid("contig postings"))?;
    }
    write_padding(file, prepared.sections[4].offset)?;
    for occurrences in prepared.seeds.values() {
        for document in documents_for(occurrences, &prepared.contig_documents) {
            let mut bytes = [0; DOCUMENT_POSTING_SIZE as usize];
            put_u32(&mut bytes, 0, document);
            file.write_all(&bytes)?;
        }
    }
    write_padding(file, prepared.sections[5].offset)?;
    for occurrences in prepared.seeds.values() {
        for occurrence in occurrences {
            file.write_all(&encode_occurrence(*occurrence))?;
        }
    }
    let file_bytes = file.stream_position()?;
    if file_bytes != prepared.file_bytes {
        return Err(JidxWriteError::Invalid("written length"));
    }
    file.flush()?;
    file.seek(SeekFrom::Start(HEADER_SIZE as u64))?;
    let body_sha256 = sha256_reader(&mut *file)?;
    let header = Header {
        k: input.k,
        seed_scheme: SeedScheme::WindowMinHash,
        posting_codec: PostingCodec::Raw,
        filter: FilterKind::None,
        document_count: u32::try_from(input.metagenomes.len())
            .map_err(|_| JidxWriteError::Invalid("metagenome count"))?,
        contig_count: u32::try_from(prepared.contig_documents.len())
            .map_err(|_| JidxWriteError::Invalid("contig count"))?,
        seed_count: u64::try_from(prepared.seeds.len())
            .map_err(|_| JidxWriteError::Invalid("seed count"))?,
        occurrence_count: prepared.occurrences,
        segment_bases: input.segment_bases,
        seeds_per_segment: input.seeds_per_segment,
        jam_sha256: input.jam_sha256,
        manifest_sha256: input.manifest_sha256,
        body_sha256,
        sections: prepared.sections,
    };
    file.seek(SeekFrom::Start(0))?;
    file.write_all(&header.encode()?)?;
    file.flush()?;
    file.sync_all()?;
    file.seek(SeekFrom::Start(0))?;
    let file_sha256 = sha256_reader(&mut *file)?;
    temporary
        .persist_noclobber(path)
        .map_err(|error| error.error)?;
    sync_directory(parent)?;
    Ok(JidxWriteStats {
        metagenomes: header.document_count,
        contigs: header.contig_count,
        seeds: header.seed_count,
        occurrences: header.occurrence_count,
        file_bytes,
        file_sha256,
    })
}

struct Prepared {
    strings: Vec<u8>,
    documents: Vec<DocumentRecord>,
    contigs: Vec<ContigRecord>,
    contig_documents: Vec<u32>,
    seeds: BTreeMap<u64, Vec<SeedOccurrence>>,
    occurrences: u64,
    sections: [SectionDescriptor; SECTION_COUNT],
    file_bytes: u64,
}

fn prepare(input: &JidxInput) -> Result<Prepared, JidxWriteError> {
    if input.metagenomes.is_empty()
        || !(1..=32).contains(&input.k)
        || input.segment_bases < u32::from(input.k)
        || input.seeds_per_segment == 0
        || u32::from(input.seeds_per_segment) > input.segment_bases - u32::from(input.k) + 1
        || input.jam_sha256 == [0; 32]
        || input.manifest_sha256 == [0; 32]
    {
        return Err(JidxWriteError::Invalid("index metadata"));
    }
    let mut metagenomes: Vec<_> = input.metagenomes.iter().collect();
    metagenomes.sort_unstable_by(|left, right| left.name.cmp(&right.name));
    if metagenomes
        .windows(2)
        .any(|pair| pair[0].name == pair[1].name)
    {
        return Err(JidxWriteError::Invalid("duplicate metagenome"));
    }

    let mut strings = Vec::new();
    let mut documents = Vec::with_capacity(metagenomes.len());
    let mut contigs = Vec::new();
    let mut contig_documents = Vec::new();
    let mut seeds = BTreeMap::<u64, Vec<SeedOccurrence>>::new();
    for (document_id, metagenome) in metagenomes.into_iter().enumerate() {
        validate_metagenome(metagenome)?;
        let document_id =
            u32::try_from(document_id).map_err(|_| JidxWriteError::Invalid("metagenome count"))?;
        let mut source_contigs: Vec<_> = metagenome.contigs.iter().collect();
        source_contigs.sort_unstable_by(|left, right| left.name.cmp(&right.name));
        if source_contigs.is_empty()
            || source_contigs
                .windows(2)
                .any(|pair| pair[0].name == pair[1].name)
        {
            return Err(JidxWriteError::Invalid("metagenome contigs"));
        }
        let contig_start =
            u32::try_from(contigs.len()).map_err(|_| JidxWriteError::Invalid("contig count"))?;
        for contig in source_contigs {
            validate_contig(contig, input)?;
            let contig_id = u32::try_from(contigs.len())
                .map_err(|_| JidxWriteError::Invalid("contig count"))?;
            contigs.push(ContigRecord {
                document_id,
                name: push_string(&mut strings, &contig.name)?,
                length: contig.length,
                fasta_offset: contig.fasta_offset,
                line_bases: contig.line_bases,
                line_width: contig.line_width,
            });
            contig_documents.push(document_id);
            for seed in &contig.seeds {
                seeds
                    .entry(seed.packed_key)
                    .or_default()
                    .push(SeedOccurrence {
                        contig_id,
                        position: seed.position,
                        canonical_orientation: seed.canonical_orientation,
                    });
            }
        }
        documents.push(DocumentRecord {
            name: push_string(&mut strings, &metagenome.name)?,
            bgzf_uri: push_string(&mut strings, &metagenome.bgzf_uri)?,
            fai_uri: push_string(&mut strings, &metagenome.fai_uri)?,
            gzi_uri: push_string(&mut strings, &metagenome.gzi_uri)?,
            bgzf_bytes: metagenome.bgzf_bytes,
            fai_bytes: metagenome.fai_bytes,
            gzi_bytes: metagenome.gzi_bytes,
            contig_start,
            contig_count: u32::try_from(contigs.len())
                .map_err(|_| JidxWriteError::Invalid("contig count"))?
                - contig_start,
            bgzf_sha256: metagenome.bgzf_sha256,
            fai_sha256: metagenome.fai_sha256,
            gzi_sha256: metagenome.gzi_sha256,
        });
    }
    for occurrences in seeds.values_mut() {
        occurrences.sort_unstable_by_key(|occurrence| (occurrence.contig_id, occurrence.position));
        if occurrences.windows(2).any(|pair| {
            (pair[0].contig_id, pair[0].position) == (pair[1].contig_id, pair[1].position)
        }) {
            return Err(JidxWriteError::Invalid("duplicate seed occurrence"));
        }
    }
    let occurrences = seeds.values().try_fold(0u64, |total, values| {
        total
            .checked_add(
                u64::try_from(values.len())
                    .map_err(|_| JidxWriteError::Invalid("occurrence count"))?,
            )
            .ok_or(JidxWriteError::Invalid("occurrence count"))
    })?;
    let document_postings = seeds.values().try_fold(0u64, |total, values| {
        total
            .checked_add(
                u64::try_from(documents_for(values, &contig_documents).len())
                    .map_err(|_| JidxWriteError::Invalid("document postings"))?,
            )
            .ok_or(JidxWriteError::Invalid("document postings"))
    })?;
    let lengths = [
        u64::try_from(strings.len()).map_err(|_| JidxWriteError::Invalid("string table"))?,
        byte_len(documents.len(), DOCUMENT_RECORD_SIZE)?,
        byte_len(contigs.len(), CONTIG_RECORD_SIZE)?,
        byte_len(seeds.len(), SEED_RECORD_SIZE)?,
        document_postings
            .checked_mul(u64::from(DOCUMENT_POSTING_SIZE))
            .ok_or(JidxWriteError::Invalid("document postings"))?,
        occurrences
            .checked_mul(u64::from(CONTIG_POSTING_SIZE))
            .ok_or(JidxWriteError::Invalid("contig postings"))?,
    ];
    let (sections, file_bytes) = section_layout(lengths)?;
    Ok(Prepared {
        strings,
        documents,
        contigs,
        contig_documents,
        seeds,
        occurrences,
        sections,
        file_bytes,
    })
}

fn validate_metagenome(input: &MetagenomeInput) -> Result<(), JidxWriteError> {
    if input.bgzf_bytes == 0
        || input.fai_bytes == 0
        || input.gzi_bytes == 0
        || input.bgzf_sha256 == [0; 32]
        || input.fai_sha256 == [0; 32]
        || input.gzi_sha256 == [0; 32]
    {
        return Err(JidxWriteError::Invalid("metagenome metadata"));
    }
    Ok(())
}

fn validate_contig(input: &ContigInput, index: &JidxInput) -> Result<(), JidxWriteError> {
    if input.length == 0 || input.line_bases == 0 || input.line_width < input.line_bases {
        return Err(JidxWriteError::Invalid("contig metadata"));
    }
    let mut seeds = input.seeds.iter().collect::<Vec<_>>();
    seeds.sort_unstable_by_key(|seed| (seed.position, seed.packed_key));
    if seeds
        .windows(2)
        .any(|pair| pair[0].position == pair[1].position)
    {
        return Err(JidxWriteError::Invalid("duplicate seed position"));
    }
    let mut segment = None;
    let mut count = 0u16;
    for seed in seeds {
        validate_packed_key(seed.packed_key, index.k)?;
        if seed
            .position
            .checked_add(u64::from(index.k))
            .is_none_or(|end| end > input.length)
        {
            return Err(JidxWriteError::Invalid("seed position"));
        }
        let current = seed.position / u64::from(index.segment_bases);
        if segment != Some(current) {
            segment = Some(current);
            count = 0;
        }
        count = count
            .checked_add(1)
            .ok_or(JidxWriteError::Invalid("segment seed count"))?;
        if count > index.seeds_per_segment {
            return Err(JidxWriteError::Invalid("segment seed count"));
        }
    }
    Ok(())
}

fn push_string(strings: &mut Vec<u8>, value: &str) -> Result<StringRef, JidxWriteError> {
    if value.is_empty() || value.bytes().any(|byte| matches!(byte, 0 | b'\n' | b'\r')) {
        return Err(JidxWriteError::Invalid("string value"));
    }
    let reference = StringRef {
        offset: u32::try_from(strings.len())
            .map_err(|_| JidxWriteError::Invalid("string table"))?,
        length: u32::try_from(value.len()).map_err(|_| JidxWriteError::Invalid("string value"))?,
    };
    strings.extend_from_slice(value.as_bytes());
    Ok(reference)
}

fn documents_for(occurrences: &[SeedOccurrence], contig_documents: &[u32]) -> Vec<u32> {
    let mut documents = Vec::new();
    for occurrence in occurrences {
        let document = contig_documents[occurrence.contig_id as usize];
        if documents.last().copied() != Some(document) {
            documents.push(document);
        }
    }
    documents
}

fn byte_len(count: usize, record_size: u32) -> Result<u64, JidxWriteError> {
    u64::try_from(count)
        .ok()
        .and_then(|count| count.checked_mul(u64::from(record_size)))
        .ok_or(JidxWriteError::Invalid("section length"))
}

fn section_layout(
    lengths: [u64; SECTION_COUNT],
) -> Result<([SectionDescriptor; SECTION_COUNT], u64), JidxWriteError> {
    let mut offset = HEADER_SIZE as u64;
    let mut sections = Vec::with_capacity(SECTION_COUNT);
    for (kind, length) in SectionKind::ALL.into_iter().zip(lengths) {
        offset = offset
            .checked_add(7)
            .ok_or(JidxWriteError::Invalid("file length"))?
            & !7;
        sections.push(SectionDescriptor {
            kind,
            record_size: kind.record_size(),
            offset,
            length,
        });
        offset = offset
            .checked_add(length)
            .ok_or(JidxWriteError::Invalid("file length"))?;
    }
    Ok((
        sections.try_into().expect("fixed JIDX section count"),
        offset,
    ))
}

fn write_padding(file: &mut File, target: u64) -> io::Result<()> {
    let position = file.stream_position()?;
    let padding = target
        .checked_sub(position)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "JIDX section overlap"))?;
    let padding = usize::try_from(padding)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "JIDX padding"))?;
    file.write_all(&vec![0; padding])
}

#[cfg(unix)]
fn sync_directory(path: &Path) -> io::Result<()> {
    File::open(path)?.sync_all()
}

#[cfg(not(unix))]
fn sync_directory(_path: &Path) -> io::Result<()> {
    Ok(())
}

#[derive(Debug, Error)]
pub enum JidxWriteError {
    #[error("JIDX write failed: {0}")]
    Io(#[from] io::Error),
    #[error(transparent)]
    Format(#[from] JidxError),
    #[error("invalid JIDX input: {0}")]
    Invalid(&'static str),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jidx_reader::JidxReader;

    fn metagenome(name: &str, contig: &str, seeds: &[(u64, u64)]) -> MetagenomeInput {
        let marker = name.as_bytes()[0];
        MetagenomeInput {
            name: name.into(),
            bgzf_uri: format!("{name}.bgz"),
            fai_uri: format!("{name}.bgz.fai"),
            gzi_uri: format!("{name}.bgz.gzi"),
            bgzf_bytes: 100,
            fai_bytes: 20,
            gzi_bytes: 16,
            bgzf_sha256: [marker; 32],
            fai_sha256: [marker.wrapping_add(1); 32],
            gzi_sha256: [marker.wrapping_add(2); 32],
            contigs: vec![ContigInput {
                name: contig.into(),
                length: 100,
                fasta_offset: 4,
                line_bases: 100,
                line_width: 101,
                seeds: seeds
                    .iter()
                    .map(|(packed_key, position)| SelectedSeed {
                        packed_key: *packed_key,
                        position: *position,
                        canonical_orientation: position % 2 == 0,
                    })
                    .collect(),
            }],
        }
    }

    fn input() -> JidxInput {
        JidxInput {
            k: 5,
            segment_bases: 32,
            seeds_per_segment: 2,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
            metagenomes: vec![
                metagenome("z", "z-contig", &[(7, 3), (9, 40)]),
                metagenome("a", "a-contig", &[(7, 2)]),
            ],
        }
    }

    #[test]
    fn writes_deterministic_complete_index() {
        let directory = tempfile::tempdir().unwrap();
        let first = directory.path().join("first.jidx");
        let second = directory.path().join("second.jidx");
        let input = input();
        let stats = write_jidx(&first, &input).unwrap();
        let mut reversed = input.clone();
        reversed.metagenomes.reverse();
        write_jidx(&second, &reversed).unwrap();
        assert_eq!(
            std::fs::read(&first).unwrap(),
            std::fs::read(&second).unwrap()
        );
        assert_eq!((stats.metagenomes, stats.contigs, stats.seeds), (2, 2, 2));

        let reader = JidxReader::open(&first).unwrap();
        reader.verify_checksum().unwrap();
        assert_eq!(reader.metagenome(0).unwrap().unwrap().name, "a");
        let seed = reader.find_seed(7).unwrap().unwrap();
        assert_eq!(reader.seed_metagenomes(seed).unwrap(), [0, 1]);
        assert_eq!(reader.seed_occurrences(seed).unwrap().len(), 2);
    }

    #[test]
    fn publication_does_not_replace_existing_output() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("index.jidx");
        write_jidx(&path, &input()).unwrap();
        let before = std::fs::read(&path).unwrap();
        assert!(write_jidx(&path, &input()).is_err());
        assert_eq!(std::fs::read(path).unwrap(), before);
    }

    #[test]
    fn rejects_more_seeds_than_the_recorded_window_policy() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("index.jidx");
        let mut input = input();
        input.metagenomes[0].contigs[0].seeds = vec![
            SelectedSeed {
                packed_key: 1,
                position: 0,
                canonical_orientation: false,
            },
            SelectedSeed {
                packed_key: 2,
                position: 1,
                canonical_orientation: false,
            },
            SelectedSeed {
                packed_key: 3,
                position: 2,
                canonical_orientation: false,
            },
        ];
        assert!(write_jidx(path, &input).is_err());
    }
}
