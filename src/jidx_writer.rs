use crate::jidx::{
    CONTIG_RECORD_SIZE, ContigRecord, DOCUMENT_POSTING_SIZE, DOCUMENT_RECORD_SIZE, DocumentRecord,
    FilterKind, HEADER_SIZE, Header, JidxError, PAGE_SIZE, PostingCodec, SECTION_COUNT,
    SectionDescriptor, SectionKind, SeedScheme, StringRef, seed_length, sha256, sha256_reader,
};
use crate::jidx_postings::{SeedEntry, SeedOccurrence, encode_seed};
use crate::jidx_runs::{RunRecord, Runs};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
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
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MetagenomeInput {
    pub name: String,
    pub bgzf_uri: String,
    pub bgzf_bytes: u64,
    pub bgzf_sha256: [u8; 32],
    pub gzi: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct JidxInput {
    pub k: u8,
    pub rescue_k15: bool,
    pub minimizer_window: u16,
    pub jam_sha256: [u8; 32],
    pub manifest_sha256: [u8; 32],
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

pub struct JidxWriter {
    path: PathBuf,
    input: JidxInput,
    scratch: tempfile::TempDir,
    runs: Runs,
    strings: Vec<u8>,
    gzi: Vec<u8>,
    documents: Vec<DocumentRecord>,
    contigs: Vec<ContigRecord>,
    current_metagenome: Option<MetagenomeInput>,
    contig_start: u32,
    contig_names: HashSet<String>,
    previous_seed: Option<(u64, u64, u8)>,
}

impl JidxWriter {
    pub fn new(path: impl AsRef<Path>, input: &JidxInput) -> Result<Self, JidxWriteError> {
        Self::with_run_records(
            path.as_ref(),
            input,
            64 * 1024 * 1024 / std::mem::size_of::<RunRecord>(),
        )
    }

    fn with_run_records(
        path: &Path,
        input: &JidxInput,
        max_records: usize,
    ) -> Result<Self, JidxWriteError> {
        if path.file_name().is_none() || path.try_exists()? {
            return Err(JidxWriteError::Invalid("output path or existing output"));
        }
        if !(1..=32).contains(&input.k)
            || (input.rescue_k15 && input.k != 21)
            || input.minimizer_window == 0
            || input.jam_sha256 == [0; 32]
            || input.manifest_sha256 == [0; 32]
        {
            return Err(JidxWriteError::Invalid("index metadata"));
        }
        let scratch = tempfile::Builder::new()
            .prefix(".jidx-runs-")
            .tempdir_in(output_parent(path))?;
        let runs = Runs::new(scratch.path(), max_records)?;
        Ok(Self {
            path: path.to_owned(),
            input: input.clone(),
            scratch,
            runs,
            strings: Vec::new(),
            gzi: Vec::new(),
            documents: Vec::new(),
            contigs: Vec::new(),
            current_metagenome: None,
            contig_start: 0,
            contig_names: HashSet::new(),
            previous_seed: None,
        })
    }

    pub fn begin_metagenome(&mut self, input: MetagenomeInput) -> Result<(), JidxWriteError> {
        validate_metagenome(&input)?;
        if self
            .current_metagenome
            .as_ref()
            .is_some_and(|previous| previous.name >= input.name)
        {
            return Err(JidxWriteError::Invalid("metagenome order"));
        }
        self.finish_metagenome()?;
        self.contig_start = u32::try_from(self.contigs.len())
            .map_err(|_| JidxWriteError::Invalid("contig count"))?;
        self.contig_names.clear();
        self.current_metagenome = Some(input);
        Ok(())
    }

    fn finish_metagenome(&mut self) -> Result<(), JidxWriteError> {
        let Some(input) = self.current_metagenome.take() else {
            return Ok(());
        };
        let contig_count = u32::try_from(self.contigs.len())
            .map_err(|_| JidxWriteError::Invalid("contig count"))?
            - self.contig_start;
        if contig_count == 0 {
            return Err(JidxWriteError::Invalid("metagenome contigs"));
        }
        self.documents.push(DocumentRecord {
            name: push_string(&mut self.strings, &input.name)?,
            bgzf_uri: push_string(&mut self.strings, &input.bgzf_uri)?,
            bgzf_bytes: input.bgzf_bytes,
            contig_start: self.contig_start,
            contig_count,
            bgzf_sha256: input.bgzf_sha256,
            gzi_offset: self.gzi.len() as u64,
            gzi_length: input.gzi.len() as u64,
        });
        self.gzi.extend_from_slice(&input.gzi);
        Ok(())
    }

    pub fn begin_contig(&mut self, input: ContigInput) -> Result<u32, JidxWriteError> {
        if self.current_metagenome.is_none() || !self.contig_names.insert(input.name.clone()) {
            return Err(JidxWriteError::Invalid("metagenome contigs"));
        }
        validate_contig(&input)?;
        let contig_id = u32::try_from(self.contigs.len())
            .map_err(|_| JidxWriteError::Invalid("contig count"))?;
        self.contigs.push(ContigRecord {
            document_id: u32::try_from(self.documents.len())
                .map_err(|_| JidxWriteError::Invalid("metagenome count"))?,
            name: push_string(&mut self.strings, &input.name)?,
            length: input.length,
            fasta_offset: input.fasta_offset,
            line_bases: input.line_bases,
            line_width: input.line_width,
        });
        self.previous_seed = None;
        Ok(contig_id)
    }

    pub fn add_seeds(
        &mut self,
        contig_id: u32,
        seeds: &[SelectedSeed],
    ) -> Result<(), JidxWriteError> {
        if self.current_metagenome.is_none()
            || contig_id < self.contig_start
            || contig_id as usize + 1 != self.contigs.len()
        {
            return Err(JidxWriteError::Invalid("active contig"));
        }
        let length = self.contigs[contig_id as usize].length;
        for seed in seeds {
            let k = seed_length(self.input.k, self.input.rescue_k15, seed.packed_key)?;
            if seed
                .position
                .checked_add(u64::from(k))
                .is_none_or(|end| end > length)
            {
                return Err(JidxWriteError::Invalid("seed position"));
            }
            if self
                .previous_seed
                .is_some_and(|(position, key, previous_k)| {
                    (position, key) >= (seed.position, seed.packed_key)
                        || (position == seed.position && previous_k == k)
                })
            {
                return Err(JidxWriteError::Invalid(
                    "duplicate or unordered seed position",
                ));
            }
            self.runs.push(RunRecord {
                packed_key: seed.packed_key,
                contig_id,
                position: seed.position,
                canonical_orientation: seed.canonical_orientation,
            })?;
            self.previous_seed = Some((seed.position, seed.packed_key, k));
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<JidxWriteStats, JidxWriteError> {
        self.finish_metagenome()?;
        if self.documents.is_empty() {
            return Err(JidxWriteError::Invalid("index metadata"));
        }
        self.contig_names.clear();
        let mut merged = self.runs.finish()?;
        let mut seed_file = BufWriter::with_capacity(
            1024 * 1024,
            File::create(self.scratch.path().join("seeds"))?,
        );
        let mut document_file = BufWriter::with_capacity(
            1024 * 1024,
            File::create(self.scratch.path().join("documents"))?,
        );
        let mut occurrence_file = BufWriter::with_capacity(
            1024 * 1024,
            File::create(self.scratch.path().join("occurrences"))?,
        );
        let mut entry: Option<SeedEntry> = None;
        let mut group: Option<DocumentGroup> = None;
        let mut previous = None;
        let mut seed_count = 0u64;
        let mut occurrence_count = 0u64;
        let mut document_bytes = 0u64;
        let mut occurrence_bytes = 0u64;
        while let Some(row) = merged.next_record()? {
            let document = self
                .contigs
                .get(row.contig_id as usize)
                .ok_or(JidxWriteError::Invalid("contig ordinal"))?
                .document_id;
            let new_key = entry.is_some_and(|entry| entry.packed_key != row.packed_key);
            if new_key
                || group
                    .as_ref()
                    .is_some_and(|group| group.document_id != document)
            {
                document_bytes = document_bytes
                    .checked_add(group.take().expect("active document").finish(
                        &mut document_file,
                        &mut occurrence_file,
                        &mut occurrence_bytes,
                    )?)
                    .ok_or(JidxWriteError::Invalid("document postings"))?;
            }
            if new_key {
                seed_file.write_all(&encode_seed(entry.take().expect("active seed")))?;
            }
            if entry.is_none() {
                seed_count = seed_count
                    .checked_add(1)
                    .ok_or(JidxWriteError::Invalid("seed count"))?;
                entry = Some(SeedEntry {
                    packed_key: row.packed_key,
                    document_frequency: 0,
                    document_offset: document_bytes,
                });
            }
            if previous == Some((row.packed_key, row.contig_id, row.position)) {
                return Err(JidxWriteError::Invalid("duplicate seed occurrence"));
            }
            previous = Some((row.packed_key, row.contig_id, row.position));
            let occurrence = SeedOccurrence {
                contig_id: row.contig_id,
                position: row.position,
                canonical_orientation: row.canonical_orientation,
            };
            if let Some(group) = &mut group {
                group.push(occurrence, &mut occurrence_file, &mut occurrence_bytes)?;
            } else {
                let entry = entry.as_mut().expect("active seed");
                entry.document_frequency = entry
                    .document_frequency
                    .checked_add(1)
                    .ok_or(JidxWriteError::Invalid("document frequency"))?;
                group = Some(DocumentGroup {
                    document_id: document,
                    contig_start: self.documents[document as usize].contig_start,
                    first: occurrence,
                    previous: occurrence,
                    count: 1,
                    offset: occurrence_bytes,
                });
            }
            occurrence_count = occurrence_count
                .checked_add(1)
                .ok_or(JidxWriteError::Invalid("occurrence count"))?;
        }
        if let Some(group) = group {
            group.finish(
                &mut document_file,
                &mut occurrence_file,
                &mut occurrence_bytes,
            )?;
        }
        if let Some(entry) = entry {
            seed_file.write_all(&encode_seed(entry))?;
        }
        seed_file.flush()?;
        document_file.flush()?;
        occurrence_file.flush()?;
        let lengths = [
            self.strings.len() as u64,
            byte_len(self.documents.len(), DOCUMENT_RECORD_SIZE)?,
            byte_len(self.contigs.len(), CONTIG_RECORD_SIZE)?,
            seed_file.get_ref().metadata()?.len(),
            document_file.get_ref().metadata()?.len(),
            occurrence_file.get_ref().metadata()?.len(),
            self.gzi.len() as u64,
            0,
        ];
        drop((merged, seed_file, document_file, occurrence_file));
        let (sections, expected_file_bytes) = section_layout(lengths)?;
        let parent = output_parent(&self.path);
        let mut temporary = tempfile::Builder::new()
            .prefix(".jidx-")
            .tempfile_in(parent)?;
        let checksum_input = temporary.reopen()?;
        let mut file = BufWriter::with_capacity(1024 * 1024, temporary.as_file_mut());
        file.write_all(&[0; HEADER_SIZE])?;
        write_padding(&mut file, sections[0].offset)?;
        file.write_all(&self.strings)?;
        write_padding(&mut file, sections[1].offset)?;
        for record in &self.documents {
            file.write_all(&record.encode())?;
        }
        write_padding(&mut file, sections[2].offset)?;
        for record in &self.contigs {
            file.write_all(&record.encode())?;
        }
        for (section, name) in sections[3..6]
            .iter()
            .zip(["seeds", "documents", "occurrences"])
        {
            write_padding(&mut file, section.offset)?;
            io::copy(&mut File::open(self.scratch.path().join(name))?, &mut file)?;
        }
        write_padding(&mut file, sections[6].offset)?;
        file.write_all(&self.gzi)?;
        write_padding(&mut file, sections[7].offset)?;
        file.flush()?;
        let mut checksum_input = BufReader::new(checksum_input);
        checksum_input.seek(SeekFrom::Start(PAGE_SIZE))?;
        let mut page = [0; PAGE_SIZE as usize];
        for _ in 0..sections[7].length / 32 {
            checksum_input.read_exact(&mut page)?;
            file.write_all(&sha256(&page))?;
        }
        let file_bytes = file.stream_position()?;
        if file_bytes != expected_file_bytes {
            return Err(JidxWriteError::Invalid("written length"));
        }
        file.flush()?;
        file.seek(SeekFrom::Start(HEADER_SIZE as u64))?;
        let body_sha256 = sha256_reader(&mut **file.get_mut())?;
        let header = Header {
            k: self.input.k,
            rescue_k15: self.input.rescue_k15,
            seed_scheme: SeedScheme::SlidingMinimizer,
            posting_codec: PostingCodec::DeltaVarint,
            filter: FilterKind::None,
            document_count: u32::try_from(self.documents.len())
                .map_err(|_| JidxWriteError::Invalid("metagenome count"))?,
            contig_count: u32::try_from(self.contigs.len())
                .map_err(|_| JidxWriteError::Invalid("contig count"))?,
            seed_count,
            occurrence_count,
            minimizer_window: self.input.minimizer_window,
            jam_sha256: self.input.jam_sha256,
            manifest_sha256: self.input.manifest_sha256,
            body_sha256,
            sections,
        };
        file.seek(SeekFrom::Start(0))?;
        file.write_all(&header.encode()?)?;
        file.flush()?;
        file.get_ref().sync_all()?;
        file.seek(SeekFrom::Start(0))?;
        let file_sha256 = sha256_reader(&mut **file.get_mut())?;
        drop(file);
        temporary
            .persist_noclobber(&self.path)
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
}

struct DocumentGroup {
    document_id: u32,
    contig_start: u32,
    first: SeedOccurrence,
    previous: SeedOccurrence,
    count: u64,
    offset: u64,
}

impl DocumentGroup {
    fn push(
        &mut self,
        occurrence: SeedOccurrence,
        file: &mut BufWriter<File>,
        bytes: &mut u64,
    ) -> Result<(), JidxWriteError> {
        if self.count == 1 {
            write_delta(file, bytes, self.contig_start, None, self.first)?;
        }
        write_delta(
            file,
            bytes,
            self.contig_start,
            Some(self.previous),
            occurrence,
        )?;
        self.previous = occurrence;
        self.count = self
            .count
            .checked_add(1)
            .ok_or(JidxWriteError::Invalid("occurrence count"))?;
        Ok(())
    }

    fn finish(
        self,
        documents: &mut BufWriter<File>,
        occurrences: &mut BufWriter<File>,
        bytes: &mut u64,
    ) -> Result<u64, JidxWriteError> {
        let local_contig = self
            .first
            .contig_id
            .checked_sub(self.contig_start)
            .ok_or(JidxWriteError::Invalid("contig ordinal"))?;
        let mut row = [0; 32];
        crate::jidx::put_u32(&mut row, 0, self.document_id);
        if self.count == 1 && local_contig != u32::MAX && self.first.position <= u64::MAX >> 1 {
            crate::jidx::put_u32(&mut row, 4, local_contig);
            crate::jidx::put_u64(
                &mut row,
                8,
                (self.first.position << 1) | u64::from(self.first.canonical_orientation),
            );
            documents.write_all(&row[..DOCUMENT_POSTING_SIZE as usize])?;
            Ok(u64::from(DOCUMENT_POSTING_SIZE))
        } else {
            if self.count == 1 {
                write_delta(occurrences, bytes, self.contig_start, None, self.first)?;
            }
            crate::jidx::put_u32(&mut row, 4, u32::MAX);
            crate::jidx::put_u64(&mut row, 8, self.offset);
            crate::jidx::put_u64(&mut row, 16, self.count);
            crate::jidx::put_u64(
                &mut row,
                24,
                bytes
                    .checked_sub(self.offset)
                    .ok_or(JidxWriteError::Invalid("occurrence length"))?,
            );
            documents.write_all(&row)?;
            Ok(row.len() as u64)
        }
    }
}

fn write_delta(
    file: &mut BufWriter<File>,
    bytes: &mut u64,
    contig_start: u32,
    previous: Option<SeedOccurrence>,
    occurrence: SeedOccurrence,
) -> Result<(), JidxWriteError> {
    let delta_contig = occurrence
        .contig_id
        .checked_sub(previous.map_or(contig_start, |value| value.contig_id))
        .ok_or(JidxWriteError::Invalid("contig order"))?;
    let position = match previous {
        Some(previous) if previous.contig_id == occurrence.contig_id => occurrence
            .position
            .checked_sub(previous.position)
            .filter(|delta| *delta != 0)
            .ok_or(JidxWriteError::Invalid("position order"))?,
        _ => occurrence.position,
    };
    let mut buffer = unsigned_varint::encode::u64_buffer();
    for value in [
        (u64::from(delta_contig) << 1) | u64::from(occurrence.canonical_orientation),
        position,
    ] {
        let encoded = unsigned_varint::encode::u64(value, &mut buffer);
        file.write_all(encoded)?;
        *bytes = bytes
            .checked_add(encoded.len() as u64)
            .ok_or(JidxWriteError::Invalid("occurrence bytes"))?;
    }
    Ok(())
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn validate_metagenome(input: &MetagenomeInput) -> Result<(), JidxWriteError> {
    if input.bgzf_bytes == 0 || input.bgzf_sha256 == [0; 32] || input.gzi.len() < 8 {
        return Err(JidxWriteError::Invalid("metagenome metadata"));
    }
    let mut gzi = noodles_bgzf::gzi::io::Reader::new(input.gzi.as_slice());
    gzi.read_index()?;
    Ok(())
}

fn validate_contig(input: &ContigInput) -> Result<(), JidxWriteError> {
    if input.length == 0
        || input.line_bases == 0
        || input.line_width < input.line_bases
        || input.line_width > input.line_bases.saturating_add(2)
    {
        return Err(JidxWriteError::Invalid("contig metadata"));
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
    for (kind, mut length) in SectionKind::ALL.into_iter().zip(lengths) {
        offset = offset
            .checked_add(PAGE_SIZE - 1)
            .ok_or(JidxWriteError::Invalid("file length"))?
            & !(PAGE_SIZE - 1);
        if kind == SectionKind::BlockChecksums {
            length = (offset / PAGE_SIZE - 1)
                .checked_mul(32)
                .ok_or(JidxWriteError::Invalid("checksum table length"))?;
        }
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

fn write_padding(file: &mut BufWriter<&mut File>, target: u64) -> io::Result<()> {
    let position = file.stream_position()?;
    let padding = target
        .checked_sub(position)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "JIDX section overlap"))?;
    let padding = usize::try_from(padding)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "JIDX padding"))?;
    file.write_all(&vec![0; padding])
}

#[cfg(unix)]
pub(crate) fn sync_directory(path: &Path) -> io::Result<()> {
    File::open(path)?.sync_all()
}

#[cfg(not(unix))]
pub(crate) fn sync_directory(_path: &Path) -> io::Result<()> {
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

    fn input() -> JidxInput {
        JidxInput {
            k: 5,
            rescue_k15: false,
            minimizer_window: 16,
            jam_sha256: [1; 32],
            manifest_sha256: [2; 32],
        }
    }

    fn metagenome(name: &str) -> MetagenomeInput {
        MetagenomeInput {
            name: name.into(),
            bgzf_uri: format!("{name}.bgz"),
            bgzf_bytes: 100,
            bgzf_sha256: [name.as_bytes()[0]; 32],
            gzi: vec![0; 8],
        }
    }

    fn contig(name: &str) -> ContigInput {
        ContigInput {
            name: name.into(),
            length: 100,
            fasta_offset: 4,
            line_bases: 100,
            line_width: 101,
        }
    }

    fn fixture(path: &Path, max_records: usize) -> Result<JidxWriteStats, JidxWriteError> {
        let mut writer = JidxWriter::with_run_records(path, &input(), max_records)?;
        for (name, positions) in [("a", &[(7, 2)][..]), ("z", &[(7, 3), (9, 40)][..])] {
            writer.begin_metagenome(metagenome(name))?;
            let id = writer.begin_contig(contig(&format!("{name}-contig")))?;
            for (packed_key, position) in positions {
                writer.add_seeds(
                    id,
                    &[SelectedSeed {
                        packed_key: *packed_key,
                        position: *position,
                        canonical_orientation: position % 2 == 0,
                    }],
                )?;
            }
        }
        writer.finish()
    }

    #[test]
    fn writes_deterministic_complete_index_across_spills() {
        let directory = tempfile::tempdir().unwrap();
        let first = directory.path().join("first.jidx");
        let second = directory.path().join("second.jidx");
        let stats = fixture(&first, 1).unwrap();
        fixture(&second, 100).unwrap();
        assert_eq!(
            std::fs::read(&first).unwrap(),
            std::fs::read(&second).unwrap()
        );
        assert_eq!(
            (
                stats.metagenomes,
                stats.contigs,
                stats.seeds,
                stats.occurrences
            ),
            (2, 2, 2, 3)
        );
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
        fixture(&path, 2).unwrap();
        let before = std::fs::read(&path).unwrap();
        assert!(fixture(&path, 2).is_err());
        assert_eq!(std::fs::read(path).unwrap(), before);
    }

    #[test]
    fn packs_document_groups_and_large_positions() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("packed.jidx");
        let mut writer = JidxWriter::with_run_records(&path, &input(), 1).unwrap();
        writer.begin_metagenome(metagenome("a")).unwrap();
        for (name, length, seeds) in [
            ("z", 100, vec![(7, 2, false), (7, 10, true)]),
            ("a", 100, vec![(7, 1, false), (9, 3, true)]),
            ("huge", u64::MAX, vec![(11, 1u64 << 63, true)]),
        ] {
            let id = writer
                .begin_contig(ContigInput {
                    length,
                    ..contig(name)
                })
                .unwrap();
            for (packed_key, position, canonical_orientation) in seeds {
                writer
                    .add_seeds(
                        id,
                        &[SelectedSeed {
                            packed_key,
                            position,
                            canonical_orientation,
                        }],
                    )
                    .unwrap();
            }
        }
        let stats = writer.finish().unwrap();
        assert_eq!((stats.seeds, stats.occurrences), (3, 5));
        let reader = JidxReader::open(&path).unwrap();
        reader.verify_checksum().unwrap();
        assert_eq!(reader.contig(0).unwrap().unwrap().name, "z");
        let seven = reader.find_seed(7).unwrap().unwrap();
        assert_eq!(
            reader.seed_occurrences(seven).unwrap(),
            [
                SeedOccurrence {
                    contig_id: 0,
                    position: 2,
                    canonical_orientation: false
                },
                SeedOccurrence {
                    contig_id: 0,
                    position: 10,
                    canonical_orientation: true
                },
                SeedOccurrence {
                    contig_id: 1,
                    position: 1,
                    canonical_orientation: false
                },
            ]
        );
        let eleven = reader.find_seed(11).unwrap().unwrap();
        assert_eq!(
            reader.seed_occurrences(eleven).unwrap()[0].position,
            1u64 << 63
        );
        let bytes = std::fs::read(&path).unwrap();
        let positions = reader.header().section(SectionKind::ContigPostings);
        let start = positions.offset as usize;
        assert_eq!(
            &bytes[start..start + positions.length as usize],
            &[
                0, 2, 1, 8, 2, 1, 5, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 1,
            ]
        );
        let documents = reader.header().section(SectionKind::DocumentPostings);
        assert_eq!(documents.length, 80);
        let start = documents.offset as usize;
        assert_eq!(&bytes[start + 4..start + 8], &u32::MAX.to_le_bytes());
        assert_eq!(&bytes[start + 16..start + 24], &3u64.to_le_bytes());
        assert_eq!(&bytes[start + 24..start + 32], &6u64.to_le_bytes());
        assert_eq!(&bytes[start + 36..start + 40], &1u32.to_le_bytes());
        assert_eq!(&bytes[start + 40..start + 48], &7u64.to_le_bytes());
        assert_eq!(&bytes[start + 56..start + 64], &6u64.to_le_bytes());
    }

    #[test]
    fn accepts_ties_and_rejects_duplicate_positions_across_batches() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("index.jidx");
        let mut writer = JidxWriter::with_run_records(&path, &input(), 1).unwrap();
        writer.begin_metagenome(metagenome("a")).unwrap();
        let id = writer.begin_contig(contig("contig")).unwrap();
        for position in 0..3 {
            writer
                .add_seeds(
                    id,
                    &[SelectedSeed {
                        packed_key: 1,
                        position,
                        canonical_orientation: false,
                    }],
                )
                .unwrap();
        }
        assert!(
            writer
                .add_seeds(
                    id,
                    &[SelectedSeed {
                        packed_key: 2,
                        position: 2,
                        canonical_orientation: false,
                    }]
                )
                .is_err()
        );
        assert!(!path.exists());
    }

    #[test]
    fn rejects_empty_or_repeated_metagenomes_and_contigs() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("index.jidx");
        let mut writer = JidxWriter::with_run_records(&path, &input(), 1).unwrap();
        writer.begin_metagenome(metagenome("a")).unwrap();
        assert!(writer.begin_metagenome(metagenome("z")).is_err());
        let mut writer = JidxWriter::with_run_records(&path, &input(), 1).unwrap();
        writer.begin_metagenome(metagenome("a")).unwrap();
        writer.begin_contig(contig("contig")).unwrap();
        assert!(writer.begin_contig(contig("contig")).is_err());
        assert!(writer.begin_metagenome(metagenome("a")).is_err());
        assert!(!path.exists());
    }
}
