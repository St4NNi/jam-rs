//! Position-free archive adapter that generates dense seeds from selected contigs.

use super::part::JamIndexPartReader;
use crate::archive::{
    ArchiveContig, ArchiveError, ArchiveMetadata, ArchiveMetrics, ArchiveResult, SeedKey,
    SeedLookupMetrics, SeedLookupResult, SeedMatch, SeedOccurrence, SeedSchemeDescriptor,
    SeedSchemeId, SequenceRequest, SequenceSlice, TraceArchive,
};
use needletail::Sequence;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, Mutex};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct ExactContigMatch {
    pub contig_id: String,
    pub strand: crate::trace::model::Strand,
    pub query_start: u64,
    pub target_start: u64,
    pub length: u64,
}

pub const K21_SCHEME: u32 = 0x4a49_0015;
pub const K31_SCHEME: u32 = 0x4a49_001f;

#[derive(Clone, Debug)]
struct LoadedContig {
    id: u32,
    name: String,
    sequence: Arc<[u8]>,
}

pub struct JamIndexArchive {
    metadata: ArchiveMetadata,
    contigs: BTreeMap<u32, LoadedContig>,
    metrics: Mutex<ArchiveMetrics>,
}

impl JamIndexArchive {
    pub fn load(
        reader: &JamIndexPartReader,
        metagenome_id: u32,
        contig_ids: impl IntoIterator<Item = u32>,
    ) -> ArchiveResult<Self> {
        let allowed = reader
            .metagenome_contigs(metagenome_id)
            .map_err(part_error)?;
        let mut selected = BTreeSet::new();
        for contig_id in contig_ids {
            if contig_id < allowed.start || contig_id >= allowed.end {
                return Err(ArchiveError::UnknownContig(contig_id));
            }
            selected.insert(contig_id);
        }
        if selected.is_empty() {
            return Err(ArchiveError::CorruptMetadata(
                "Jam Index candidate has no selected contigs".to_string(),
            ));
        }
        let selected_ids = selected.iter().copied().collect::<Vec<_>>();
        let loaded = reader
            .read_contigs(metagenome_id, &selected_ids)
            .map_err(part_error)?;
        let loaded_sequences = loaded
            .contigs
            .into_iter()
            .map(|(contig_id, sequence)| (contig_id, Arc::<[u8]>::from(sequence)))
            .collect::<BTreeMap<_, _>>();
        Self::from_loaded(
            reader,
            metagenome_id,
            &loaded_sequences,
            selected,
            loaded.source_bytes,
        )
    }

    pub(crate) fn from_loaded(
        reader: &JamIndexPartReader,
        metagenome_id: u32,
        loaded: &BTreeMap<u32, Arc<[u8]>>,
        contig_ids: impl IntoIterator<Item = u32>,
        source_bytes: u64,
    ) -> ArchiveResult<Self> {
        let metagenome =
            reader
                .metagenomes()
                .get(usize::try_from(metagenome_id).map_err(|_| {
                    ArchiveError::CorruptMetadata("metagenome ID overflow".to_string())
                })?)
                .ok_or_else(|| {
                    ArchiveError::CorruptMetadata(format!(
                        "unknown Jam Index metagenome {metagenome_id}"
                    ))
                })?;
        let allowed = reader
            .metagenome_contigs(metagenome_id)
            .map_err(part_error)?;
        let mut selected = BTreeSet::new();
        for contig_id in contig_ids {
            if !allowed.contains(&contig_id) {
                return Err(ArchiveError::UnknownContig(contig_id));
            }
            selected.insert(contig_id);
        }
        if selected.is_empty() {
            return Err(ArchiveError::CorruptMetadata(
                "Jam Index candidate has no selected contigs".to_string(),
            ));
        }
        let mut contigs = BTreeMap::new();
        for contig_id in selected {
            let descriptor =
                reader
                    .contigs()
                    .get(usize::try_from(contig_id).map_err(|_| {
                        ArchiveError::CorruptMetadata("contig ID overflow".to_string())
                    })?)
                    .ok_or(ArchiveError::UnknownContig(contig_id))?;
            let sequence = loaded
                .get(&contig_id)
                .cloned()
                .ok_or(ArchiveError::UnknownContig(contig_id))?;
            contigs.insert(
                contig_id,
                LoadedContig {
                    id: contig_id,
                    name: descriptor.name.clone(),
                    sequence,
                },
            );
        }
        let archive_contigs = contigs
            .values()
            .map(|contig| ArchiveContig {
                id: contig.id,
                name: contig.name.clone(),
                length: u64::try_from(contig.sequence.len()).unwrap_or(u64::MAX),
            })
            .collect::<Vec<_>>();
        let total_bases = archive_contigs
            .iter()
            .map(|contig| contig.length)
            .fold(0u64, u64::saturating_add);
        Ok(Self {
            metadata: ArchiveMetadata {
                format_identifier: "jam-index-part".to_string(),
                format_version: 1,
                layout_identifier: 1,
                source_assembly_sha256: metagenome.source_sha256,
                archive_sha256: None,
                builder_version: env!("CARGO_PKG_VERSION").to_string(),
                source_commit: option_env!("JAM_RS_SOURCE_COMMIT")
                    .unwrap_or("unknown")
                    .to_string(),
                hash_algorithm: crate::archive::JAMHASH_ALGORITHM_ID.to_string(),
                total_bases,
                contigs: archive_contigs,
            },
            contigs,
            metrics: Mutex::new(ArchiveMetrics {
                sequence_bytes_read: source_bytes,
                decoded_sequence_bases: total_bases,
                ..ArchiveMetrics::default()
            }),
        })
    }

    pub(crate) fn exact_matches(&self, query: &[u8], allow_wrap: bool) -> Vec<ExactContigMatch> {
        let _profile = crate::profiling::scope("exact_sequence_scan");
        if query.is_empty() {
            return Vec::new();
        }
        let reverse_query = query.reverse_complement();
        let mut doubled_query = Vec::new();
        if allow_wrap {
            doubled_query.reserve(query.len().saturating_mul(2).saturating_sub(1));
            doubled_query.extend_from_slice(query);
            doubled_query.extend_from_slice(&query[..query.len().saturating_sub(1)]);
        }
        let mut matches = Vec::new();
        for contig in self.contigs.values() {
            let found = find_exact(contig.sequence.as_ref(), query)
                .map(|target_start| (crate::trace::model::Strand::Forward, 0, target_start));
            let found = found.or_else(|| {
                find_exact(contig.sequence.as_ref(), &reverse_query)
                    .map(|target_start| (crate::trace::model::Strand::Reverse, 0, target_start))
            });
            let found = found.or_else(|| {
                (allow_wrap && contig.sequence.len() == query.len())
                    .then(|| {
                        find_exact(&doubled_query, contig.sequence.as_ref()).map(|query_start| {
                            (crate::trace::model::Strand::Forward, query_start, 0)
                        })
                    })
                    .flatten()
            });
            let found = found.or_else(|| {
                (allow_wrap && contig.sequence.len() == query.len())
                    .then(|| {
                        let reverse_contig = contig.sequence.reverse_complement();
                        find_exact(&doubled_query, &reverse_contig).map(|query_start| {
                            (crate::trace::model::Strand::Reverse, query_start, 0)
                        })
                    })
                    .flatten()
            });
            if let Some((strand, query_start, target_start)) = found {
                matches.push(ExactContigMatch {
                    contig_id: contig.name.clone(),
                    strand,
                    query_start: u64::try_from(query_start).unwrap_or(u64::MAX),
                    target_start: u64::try_from(target_start).unwrap_or(u64::MAX),
                    length: u64::try_from(query.len()).unwrap_or(u64::MAX),
                });
            }
        }
        matches
    }

    fn scheme(&self, scheme: SeedSchemeId) -> ArchiveResult<u8> {
        match scheme.0 {
            K21_SCHEME => Ok(21),
            K31_SCHEME => Ok(31),
            other => Err(ArchiveError::UnknownSeedScheme(other)),
        }
    }

    fn lookup(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
        max_occurrences: Option<u32>,
    ) -> ArchiveResult<SeedLookupResult> {
        let k = self.scheme(scheme)?;
        let _profile = crate::profiling::scope(if k == 21 {
            "dense_k21_generation"
        } else {
            "dense_k31_generation"
        });
        let width = usize::from(k).div_ceil(4);
        let mut requested = BTreeMap::<u64, Vec<(usize, u64)>>::new();
        for (index, key) in keys.iter().enumerate() {
            if key.verification.len() != width {
                return Err(ArchiveError::CorruptMetadata(format!(
                    "seed verification width does not match k={k}"
                )));
            }
            let packed = key.verification.iter().fold(0u64, |value, byte| {
                value.checked_shl(8).unwrap_or(0) | u64::from(*byte)
            });
            requested
                .entry(key.digest)
                .or_default()
                .push((index, packed));
        }
        let mut occurrences = vec![Vec::<SeedOccurrence>::new(); keys.len()];
        let mut counts = vec![0u64; keys.len()];
        let mut tested = 0u64;
        for contig in self.contigs.values() {
            let normalized = contig.sequence.as_ref().normalize(false);
            for (position, kmer, reverse) in normalized.bit_kmers(k, true) {
                tested = tested.saturating_add(1);
                let hash = crate::jamhash_u64_v1(kmer.0);
                let Some(candidates) = requested.get(&hash) else {
                    continue;
                };
                for &(key_index, packed) in candidates {
                    if packed != kmer.0 {
                        continue;
                    }
                    counts[key_index] = counts[key_index].saturating_add(1);
                    if max_occurrences.is_none_or(|limit| counts[key_index] <= u64::from(limit)) {
                        occurrences[key_index].push(SeedOccurrence {
                            contig_id: contig.id,
                            position: u64::try_from(position).map_err(|_| {
                                ArchiveError::CorruptMetadata("seed position overflow".to_string())
                            })?,
                            span: u16::from(k),
                            reverse,
                        });
                    } else {
                        occurrences[key_index].clear();
                    }
                }
            }
        }
        crate::profiling::add_counter("dense_signature_comparisons", tested);
        let before = counts.iter().copied().sum();
        let after = occurrences
            .iter()
            .map(|items| u64::try_from(items.len()).unwrap_or(u64::MAX))
            .sum();
        let matches = occurrences
            .into_iter()
            .enumerate()
            .filter_map(|(key_index, occurrences)| {
                (!occurrences.is_empty()).then_some(SeedMatch {
                    key_index,
                    occurrences,
                })
            })
            .collect();
        Ok(SeedLookupResult {
            matches,
            metrics: SeedLookupMetrics {
                pages_read: 0,
                keys_tested: tested,
                occurrences_before_limits: before,
                occurrences_after_limits: after,
                bytes_read: 0,
            },
        })
    }
}

impl TraceArchive for JamIndexArchive {
    fn metadata(&self) -> ArchiveResult<ArchiveMetadata> {
        Ok(self.metadata.clone())
    }

    fn available_seed_schemes(&self) -> ArchiveResult<Vec<SeedSchemeDescriptor>> {
        Ok([K21_SCHEME, K31_SCHEME]
            .into_iter()
            .map(|scheme_id| {
                let k = if scheme_id == K21_SCHEME { 21 } else { 31 };
                SeedSchemeDescriptor {
                    scheme_id,
                    algorithm_id: 1,
                    span: k,
                    informative_bases: k,
                    density_parameter: 1,
                    bucket_bits: 0,
                    key_encoding: 1,
                    occurrence_encoding: 1,
                    flags: 0,
                }
            })
            .collect())
    }

    fn lookup_seeds(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<SeedLookupResult> {
        self.lookup(scheme, keys, None)
    }

    fn lookup_seeds_bounded(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
        max_occurrences: Option<u32>,
    ) -> ArchiveResult<SeedLookupResult> {
        self.lookup(scheme, keys, max_occurrences)
    }

    fn read_sequences(&self, requests: &[SequenceRequest]) -> ArchiveResult<Vec<SequenceSlice>> {
        let mut slices = Vec::with_capacity(requests.len());
        let mut decoded = 0u64;
        for request in requests {
            let contig = self
                .contigs
                .get(&request.contig_id)
                .ok_or(ArchiveError::UnknownContig(request.contig_id))?;
            let start = usize::try_from(request.start).map_err(|_| {
                ArchiveError::CorruptMetadata("sequence start overflow".to_string())
            })?;
            let end = usize::try_from(request.end)
                .map_err(|_| ArchiveError::CorruptMetadata("sequence end overflow".to_string()))?;
            let mut bases = contig
                .sequence
                .get(start..end)
                .ok_or(ArchiveError::InvalidSequenceRange {
                    start: request.start,
                    end: request.end,
                })?
                .to_vec();
            if request.reverse_complement {
                bases.reverse();
                for base in &mut bases {
                    *base = crate::sequence::complement_base(*base)
                        .map_err(|error| ArchiveError::CorruptMetadata(error.to_string()))?;
                }
            }
            decoded = decoded.saturating_add(u64::try_from(bases.len()).unwrap_or(u64::MAX));
            slices.push(SequenceSlice {
                request: *request,
                bases,
            });
        }
        let mut metrics = self
            .metrics
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        metrics.decoded_sequence_bases = metrics.decoded_sequence_bases.saturating_add(decoded);
        Ok(slices)
    }

    fn metrics(&self) -> ArchiveMetrics {
        *self
            .metrics
            .lock()
            .unwrap_or_else(|error| error.into_inner())
    }
}

fn find_exact(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() {
        return Some(0);
    }
    let mut prefix = vec![0usize; needle.len()];
    let mut matched = 0usize;
    for index in 1..needle.len() {
        while matched > 0 && !base_equal(needle[index], needle[matched]) {
            matched = prefix[matched - 1];
        }
        if base_equal(needle[index], needle[matched]) {
            matched += 1;
        }
        prefix[index] = matched;
    }
    matched = 0;
    for (index, base) in haystack.iter().copied().enumerate() {
        while matched > 0 && !base_equal(base, needle[matched]) {
            matched = prefix[matched - 1];
        }
        if base_equal(base, needle[matched]) {
            matched += 1;
        }
        if matched == needle.len() {
            return Some(index + 1 - needle.len());
        }
    }
    None
}

fn base_equal(left: u8, right: u8) -> bool {
    let normalize = |base: u8| match base.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    };
    normalize(left) == normalize(right)
}

fn part_error(error: super::part::JamIndexPartError) -> ArchiveError {
    ArchiveError::Backend(error.to_string())
}
