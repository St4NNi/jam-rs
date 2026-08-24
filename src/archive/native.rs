//! Native JMA implementation of the backend-neutral trace archive contract.

use super::{
    ArchiveContig, ArchiveError, ArchiveMetadata, ArchiveMetrics, ArchiveResult, SeedKey,
    SeedLookupResult, SeedMatch, SeedOccurrence, SeedSchemeDescriptor, SeedSchemeId,
    SequenceRequest, SequenceSlice, TraceArchive,
};
use crate::jma::reader::JmaReader;
use crate::jma::{ArchiveReader, SeedQuery, SequenceRange};
use crate::resource::RangeReader;

pub struct NativeJmaArchive<R> {
    reader: JmaReader<R>,
}

impl<R: RangeReader> NativeJmaArchive<R> {
    pub fn from_resource(resource: R) -> ArchiveResult<Self> {
        Ok(Self {
            reader: JmaReader::from_resource(resource).map_err(backend_error)?,
        })
    }

    #[must_use]
    pub fn reader(&self) -> &JmaReader<R> {
        &self.reader
    }

    pub fn into_reader(self) -> JmaReader<R> {
        self.reader
    }
}

impl<R: RangeReader> TraceArchive for NativeJmaArchive<R> {
    fn metadata(&self) -> ArchiveResult<ArchiveMetadata> {
        let metadata = self.reader.metadata_fields();
        Ok(ArchiveMetadata {
            format_identifier: metadata.format_identifier.clone(),
            format_version: metadata.format_version,
            layout_identifier: metadata.layout_identifier,
            source_assembly_sha256: metadata.source_assembly_sha256,
            archive_sha256: metadata.archive_sha256,
            builder_version: metadata.builder_version.clone(),
            source_commit: metadata.source_commit.clone(),
            hash_algorithm: metadata.hash_algorithm.clone(),
            total_bases: self.reader.header().total_bases,
            contigs: self
                .reader
                .contigs()
                .iter()
                .map(|contig| ArchiveContig {
                    id: contig.id,
                    name: contig.name.clone(),
                    length: contig.length,
                })
                .collect(),
        })
    }

    fn available_seed_schemes(&self) -> ArchiveResult<Vec<SeedSchemeDescriptor>> {
        Ok(self.reader.seed_schemes().copied().collect())
    }

    fn lookup_seeds(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<SeedLookupResult> {
        let descriptor = self
            .reader
            .seed_schemes()
            .find(|descriptor| descriptor.scheme_id == scheme.0)
            .copied()
            .ok_or(ArchiveError::UnknownSeedScheme(scheme.0))?;
        let k = u8::try_from(descriptor.span).map_err(|_| {
            ArchiveError::CorruptMetadata("native seed span exceeds 8 bits".to_string())
        })?;
        let expected_width = usize::from(k).div_ceil(4);
        if k == 0 || k > 32 || descriptor.key_encoding != 1 {
            return Err(ArchiveError::CorruptMetadata(format!(
                "unsupported native seed descriptor {}",
                descriptor.scheme_id
            )));
        }
        let before = self.reader.archive_metrics();
        let mut matches = Vec::new();
        for (key_index, key) in keys.iter().enumerate() {
            if key.verification.len() != expected_width {
                return Err(ArchiveError::CorruptMetadata(format!(
                    "seed verification width {} does not match k={k}",
                    key.verification.len()
                )));
            }
            let canonical_kmer = key
                .verification
                .iter()
                .try_fold(0u64, |value, byte| {
                    value.checked_shl(8).map(|value| value | u64::from(*byte))
                })
                .ok_or_else(|| {
                    ArchiveError::CorruptMetadata("seed verification overflows u64".to_string())
                })?;
            if !crate::jma::format::valid_canonical_kmer(k, canonical_kmer) {
                return Err(ArchiveError::CorruptMetadata(
                    "seed verification contains non-zero padding bits".to_string(),
                ));
            }
            let occurrences = self
                .reader
                .seed_occurrences_for_scheme(
                    scheme,
                    SeedQuery {
                        k,
                        hash: key.digest,
                        canonical_kmer,
                    },
                )
                .map_err(backend_error)?
                .into_iter()
                .map(|occurrence| SeedOccurrence {
                    contig_id: occurrence.contig_id,
                    position: occurrence.position,
                    span: descriptor.span,
                    reverse: occurrence.reverse,
                })
                .collect::<Vec<_>>();
            if !occurrences.is_empty() {
                matches.push(SeedMatch {
                    key_index,
                    occurrences,
                });
            }
        }
        let after = self.reader.archive_metrics();
        let occurrence_count = matches.iter().try_fold(0u64, |count, seed_match| {
            count.checked_add(u64::try_from(seed_match.occurrences.len()).ok()?)
        });
        let occurrence_count = occurrence_count
            .ok_or_else(|| ArchiveError::Backend("seed occurrence count overflow".to_string()))?;
        Ok(SeedLookupResult {
            matches,
            metrics: super::SeedLookupMetrics {
                pages_read: after
                    .resource
                    .seed_buckets_read
                    .saturating_sub(before.resource.seed_buckets_read),
                keys_tested: u64::try_from(keys.len()).unwrap_or(u64::MAX),
                occurrences_before_limits: occurrence_count,
                occurrences_after_limits: occurrence_count,
                bytes_read: after.seed_bytes_read.saturating_sub(before.seed_bytes_read),
            },
        })
    }

    fn read_sequences(&self, requests: &[SequenceRequest]) -> ArchiveResult<Vec<SequenceSlice>> {
        if let Some(request) = requests.iter().find(|request| request.start > request.end) {
            return Err(ArchiveError::InvalidSequenceRange {
                start: request.start,
                end: request.end,
            });
        }
        let mut order = (0..requests.len()).collect::<Vec<_>>();
        order.sort_by_key(|index| {
            let request = requests[*index];
            (request.contig_id, request.start, request.end, *index)
        });
        let mut output = vec![None; requests.len()];
        let mut cursor = 0usize;
        while cursor < order.len() {
            let first = requests[order[cursor]];
            let mut end = first.end;
            let mut group_end = cursor + 1;
            while group_end < order.len() {
                let next = requests[order[group_end]];
                if next.contig_id != first.contig_id || next.start > end {
                    break;
                }
                end = end.max(next.end);
                group_end += 1;
            }
            let bases = self
                .reader
                .read_sequence(
                    first.contig_id,
                    SequenceRange::new(first.start, end).map_err(backend_error)?,
                )
                .map_err(backend_error)?;
            for &index in &order[cursor..group_end] {
                let request = requests[index];
                let start = usize::try_from(request.start - first.start).map_err(|_| {
                    ArchiveError::InvalidSequenceRange {
                        start: request.start,
                        end: request.end,
                    }
                })?;
                let end = usize::try_from(request.end - first.start).map_err(|_| {
                    ArchiveError::InvalidSequenceRange {
                        start: request.start,
                        end: request.end,
                    }
                })?;
                let mut selected = bases
                    .get(start..end)
                    .ok_or(ArchiveError::InvalidSequenceRange {
                        start: request.start,
                        end: request.end,
                    })?
                    .to_vec();
                if request.reverse_complement {
                    selected.reverse();
                    for base in &mut selected {
                        *base = crate::sequence::complement_base(*base)
                            .map_err(|error| ArchiveError::Backend(error.to_string()))?;
                    }
                }
                output[index] = Some(SequenceSlice {
                    request,
                    bases: selected,
                });
            }
            cursor = group_end;
        }
        output
            .into_iter()
            .map(|slice| {
                slice.ok_or_else(|| {
                    ArchiveError::Backend("sequence request was not processed".to_string())
                })
            })
            .collect()
    }

    fn metrics(&self) -> ArchiveMetrics {
        self.reader.archive_metrics()
    }
}

fn backend_error(error: crate::jma::JmaError) -> ArchiveError {
    match error {
        crate::jma::JmaError::InvalidSequenceRange { start, end } => {
            ArchiveError::InvalidSequenceRange { start, end }
        }
        crate::jma::JmaError::UnknownContig(contig) => ArchiveError::UnknownContig(contig),
        crate::jma::JmaError::ChecksumMismatch(section) => ArchiveError::ChecksumMismatch(section),
        crate::jma::JmaError::Resource(error) => ArchiveError::Resource(error),
        error => ArchiveError::Backend(error.to_string()),
    }
}
