use crate::jidx::{seed_length, sha256};
use crate::jidx_reader::{Contig, Metagenome};
use crate::owner_file::OwnerFile;
use crate::owner_format::{
    OWNER_HEADER_SIZE, OwnerDocument, OwnerHeader, OwnerReaderError, OwnerSeed,
};
use crate::owner_observer::{OwnerReadObserver, OwnerReadSnapshot};
use serde::Deserialize;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct OwnerManifest {
    version: u32,
    complete: bool,
    metadata_owner: u32,
    owners: Vec<OwnerManifestFile>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct OwnerManifestFile {
    path: PathBuf,
    header_sha256: String,
}

pub(crate) struct OwnerReader {
    root_sha256: [u8; 32],
    files: Vec<OwnerFile>,
    ranges: Vec<(u64, u64, usize)>,
    metadata_owner: usize,
    complete: bool,
    document_widths: Vec<u8>,
    observer: Arc<OwnerReadObserver>,
    report_on_drop: bool,
}

impl OwnerReader {
    pub(crate) fn open(path: impl AsRef<Path>) -> Result<Self, OwnerReaderError> {
        let report = std::env::var_os("JAM_OWNER_OBSERVE").is_some_and(|value| value == "1");
        let observer = if report {
            OwnerReadObserver::enabled(4096).map_err(OwnerReaderError::Invalid)?
        } else {
            OwnerReadObserver::disabled()
        };
        let mut reader = Self::open_with_observer(path, Arc::new(observer))?;
        reader.report_on_drop = report;
        Ok(reader)
    }

    pub(crate) fn open_with_observer(
        path: impl AsRef<Path>,
        observer: Arc<OwnerReadObserver>,
    ) -> Result<Self, OwnerReaderError> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)?;
        observer.record_open(bytes.len() as u64, 0, 0, 0, 0);
        let manifest: OwnerManifest = serde_json::from_slice(&bytes)
            .map_err(|_| OwnerReaderError::Invalid("owner manifest JSON"))?;
        if manifest.version != 1 || manifest.owners.is_empty() {
            return Err(OwnerReaderError::Invalid("owner manifest"));
        }
        let expected_count = u32::try_from(manifest.owners.len())
            .map_err(|_| OwnerReaderError::Invalid("owner count"))?;
        if manifest.metadata_owner >= expected_count {
            return Err(OwnerReaderError::Invalid("metadata owner"));
        }
        let parent = path.parent().unwrap_or_else(|| Path::new("."));
        let mut files = Vec::with_capacity(manifest.owners.len());
        let hot_cache_bytes = Arc::new(AtomicUsize::new(0));
        for (ordinal, entry) in manifest.owners.into_iter().enumerate() {
            let file = OwnerFile::open(
                parent.join(entry.path),
                Arc::clone(&hot_cache_bytes),
                Arc::clone(&observer),
            )?;
            if file.header.owner_ordinal != ordinal as u32
                || file.header.owner_count != expected_count
                || hex(&sha256(&file.mmap[..OWNER_HEADER_SIZE])) != entry.header_sha256
            {
                return Err(OwnerReaderError::Invalid("owner identity"));
            }
            files.push(file);
        }
        let contract = &files[0].header;
        if files.iter().any(|file| {
            let header = &file.header;
            header.generation_id != contract.generation_id
                || (header.k, header.rescue_k15, header.minimizer_window)
                    != (contract.k, contract.rescue_k15, contract.minimizer_window)
                || (header.document_count, header.contig_count)
                    != (contract.document_count, contract.contig_count)
                || header.metadata_sha256 != contract.metadata_sha256
        }) {
            return Err(OwnerReaderError::Invalid("owner generation contract"));
        }
        if files.iter().enumerate().any(|(index, file)| {
            file.header.has_metadata() != (index == manifest.metadata_owner as usize)
        }) {
            return Err(OwnerReaderError::Invalid("metadata owner"));
        }
        let mut ranges = files
            .iter()
            .enumerate()
            .map(|(index, file)| (file.header.first_key, file.header.last_key, index))
            .collect::<Vec<_>>();
        ranges.sort_unstable();
        for pair in ranges.windows(2) {
            if pair[0].1 >= pair[1].0 {
                return Err(OwnerReaderError::Invalid("overlapping owner ranges"));
            }
        }
        if manifest.complete
            && (ranges[0].0 != 0
                || ranges.last().expect("nonempty ranges").1 != u64::MAX
                || ranges
                    .iter()
                    .any(|&(_, _, index)| !files[index].header.complete_range())
                || ranges
                    .windows(2)
                    .any(|pair| pair[0].1.checked_add(1) != Some(pair[1].0)))
        {
            return Err(OwnerReaderError::Invalid("incomplete owner generation"));
        }
        let document_widths = (0..contract.document_count)
            .map(|id| {
                files[manifest.metadata_owner as usize]
                    .document_record(id)
                    .map(|record| record.locus_bits)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            root_sha256: sha256(&bytes),
            files,
            ranges,
            metadata_owner: manifest.metadata_owner as usize,
            complete: manifest.complete,
            document_widths,
            observer,
            report_on_drop: false,
        })
    }

    pub(crate) fn read_snapshot(&self) -> OwnerReadSnapshot {
        self.observer.snapshot()
    }

    pub(crate) fn header(&self) -> &OwnerHeader {
        &self.files[self.metadata_owner].header
    }

    pub(crate) fn cache_identity(&self) -> io::Result<Option<[u8; 32]>> {
        if self.files.iter().any(|file| file.file_identity.is_none()) {
            return Ok(None);
        }
        let mut bytes = Vec::with_capacity(32 + self.files.len() * 7 * 8);
        bytes.extend_from_slice(&self.root_sha256);
        for file in &self.files {
            file.verify_unchanged()?;
            for value in file.file_identity.expect("checked owner file identity") {
                bytes.extend_from_slice(&value.to_le_bytes());
            }
        }
        Ok(Some(sha256(&bytes)))
    }

    pub(crate) fn k(&self) -> u8 {
        self.header().k
    }

    pub(crate) fn rescue_k15(&self) -> bool {
        self.header().rescue_k15
    }

    pub(crate) fn document_count(&self) -> u32 {
        self.header().document_count
    }

    pub(crate) fn header_sha256(&self) -> [u8; 32] {
        let mut bytes = Vec::with_capacity(self.files.len() * 32);
        for file in &self.files {
            bytes.extend_from_slice(&sha256(&file.mmap[..OWNER_HEADER_SIZE]));
        }
        sha256(&bytes)
    }

    pub(crate) fn manifest_sha256(&self) -> [u8; 32] {
        self.root_sha256
    }

    pub(crate) fn body_sha256(&self) -> [u8; 32] {
        let mut bytes = Vec::with_capacity(self.files.len() * 32);
        for file in &self.files {
            bytes.extend_from_slice(&file.header.body_sha256);
        }
        sha256(&bytes)
    }

    pub(crate) fn is_complete(&self) -> bool {
        self.complete
    }

    pub(crate) fn find_seeds_batch(
        &self,
        keys: &[u64],
    ) -> Result<Vec<Option<OwnerSeed>>, OwnerReaderError> {
        let mut out = vec![None; keys.len()];
        let mut grouped = std::collections::BTreeMap::<usize, Vec<(usize, u64)>>::new();
        for (position, &key) in keys.iter().enumerate() {
            let owner = self.owner_for(key)?;
            grouped.entry(owner).or_default().push((position, key));
        }
        let touched = grouped.len() as u64;
        for (owner, requests) in grouped {
            let found = self.files[owner]
                .find_seeds_batch(&requests.iter().map(|&(_, key)| key).collect::<Vec<_>>())?;
            for ((position, _), seed) in requests.into_iter().zip(found) {
                out[position] = seed;
            }
        }
        let present = out.iter().filter(|seed| seed.is_some()).count() as u64;
        self.observer.record_route(
            keys.len() as u64,
            present,
            keys.len() as u64 - present,
            touched,
        );
        Ok(out)
    }

    pub(crate) fn seed_documents(
        &self,
        seed: OwnerSeed,
    ) -> Result<Vec<OwnerDocument>, OwnerReaderError> {
        let owner = self.owner_for(seed.packed_key)?;
        self.files[owner].seed_documents(seed)
    }

    pub(crate) fn seed_document_occurrences(
        &self,
        seed: OwnerSeed,
        document: OwnerDocument,
    ) -> Result<Vec<crate::jidx_reader::SeedOccurrence>, OwnerReaderError> {
        if document.seed_key != seed.packed_key {
            return Err(OwnerReaderError::Invalid("seed document"));
        }
        let owner = usize::try_from(document.owner_ordinal)
            .map_err(|_| OwnerReaderError::Invalid("owner ordinal"))?;
        let local = self
            .files
            .get(owner)
            .ok_or(OwnerReaderError::Invalid("owner ordinal"))?
            .document_occurrences(document, &self.document_widths)?;
        let k = seed_length(self.k(), self.rescue_k15(), seed.packed_key)
            .map_err(|_| OwnerReaderError::Invalid("seed key"))?;
        self.files[self.metadata_owner].locus_occurrences(document.metagenome_id, &local, k)
    }

    pub(crate) fn metagenome_name(&self, id: u32) -> Result<Option<&str>, OwnerReaderError> {
        Ok(self.metagenome(id)?.map(|record| record.name))
    }

    pub(crate) fn metagenome(&self, id: u32) -> Result<Option<Metagenome<'_>>, OwnerReaderError> {
        self.files[self.metadata_owner].metagenome(id)
    }

    pub(crate) fn contig(&self, id: u32) -> Result<Option<Contig<'_>>, OwnerReaderError> {
        self.files[self.metadata_owner].contig(id)
    }

    pub(crate) fn verify_checksum(&self) -> Result<(), OwnerReaderError> {
        for file in &self.files {
            file.verify_checksum()?;
        }
        self.files[self.metadata_owner].audit_metadata()?;
        for (owner, file) in self.files.iter().enumerate() {
            file.validate_directory()?;
            let mut keys = 0u64;
            let mut occurrences = 0u64;
            for ordinal in 0..file.block_count()? {
                let record = file.block_record(ordinal)?;
                let decoded_hot = file.hot_block(ordinal, record)?;
                keys = keys
                    .checked_add(decoded_hot.keys.len() as u64)
                    .ok_or(OwnerReaderError::Invalid("key count"))?;
                for entry in &decoded_hot.keys {
                    let seed = OwnerSeed {
                        packed_key: entry.key,
                        document_frequency: entry.document_frequency,
                        block_ordinal: ordinal,
                    };
                    for member in crate::owner_postings::key_members(&decoded_hot, entry)? {
                        let decoded = self.seed_document_occurrences(
                            seed,
                            OwnerDocument {
                                metagenome_id: member.document_id,
                                occurrence_count: member.occurrence_count,
                                seed_key: entry.key,
                                owner_ordinal: owner as u32,
                                block_ordinal: ordinal,
                                member_ordinal: member.member_ordinal,
                            },
                        )?;
                        if decoded.len() as u64 != member.occurrence_count {
                            return Err(OwnerReaderError::Invalid("occurrence count"));
                        }
                        occurrences = occurrences
                            .checked_add(member.occurrence_count)
                            .ok_or(OwnerReaderError::Invalid("occurrence count"))?;
                    }
                }
            }
            if keys != file.header.key_count || occurrences != file.header.occurrence_count {
                return Err(OwnerReaderError::Invalid("owner totals"));
            }
        }
        Ok(())
    }

    fn owner_for(&self, key: u64) -> Result<usize, OwnerReaderError> {
        let position = self.ranges.partition_point(|range| range.0 <= key);
        if position == 0 || key > self.ranges[position - 1].1 {
            return Err(OwnerReaderError::KeyNotCovered(key));
        }
        Ok(self.ranges[position - 1].2)
    }
}

impl Drop for OwnerReader {
    fn drop(&mut self) {
        if self.report_on_drop {
            if let Ok(json) = serde_json::to_string(&self.read_snapshot()) {
                eprintln!("owner_read_stats {json}");
            }
        }
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}
