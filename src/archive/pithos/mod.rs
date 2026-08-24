//! Backend-neutral adapter for the Pithos biosequence profile.
//!
//! Pithos stores the logical trace members (metadata, contig/seed indexes, and
//! sequence slabs) separately.  This module deliberately does not know about
//! the Pithos container bytes: the source implementation supplies profile
//! records and performs member/range reads.  That keeps the trace path free of
//! a JMA decoder and lets the same adapter be used by mmap and HTTP range
//! sources.

use crate::archive::{
    ArchiveContig, ArchiveError, ArchiveMetadata, ArchiveMetrics, ArchiveResult, SeedKey,
    SeedLookupMetrics, SeedMatch, SeedOccurrence, SeedSchemeDescriptor, SeedSchemeId,
    SequenceRequest, SequenceSlice, TraceArchive,
};
use crate::resource::ResourceMetrics;
use crate::sequence::{complement_base, encode_base, normalize_ambiguity};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};

/// Physical sequence organization used by a Pithos biosequence profile.
///
/// Both organizations expose the same logical contig-range API.  A source is
/// responsible for selecting only the member or packed slab ranges needed by
/// a request.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PithosSequenceOrganization {
    OneMemberPerContig,
    PackedSlabs,
}

/// Metadata and index descriptors loaded from the Pithos profile directory.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PithosProfile {
    pub metadata: ArchiveMetadata,
    pub seed_schemes: Vec<SeedSchemeDescriptor>,
    pub sequence_organization: PithosSequenceOrganization,
}

impl PithosProfile {
    fn validate(&self) -> ArchiveResult<()> {
        let mut contig_ids = self
            .metadata
            .contigs
            .iter()
            .map(|contig| contig.id)
            .collect::<Vec<_>>();
        contig_ids.sort_unstable();
        if contig_ids.windows(2).any(|window| window[0] == window[1]) {
            return Err(ArchiveError::CorruptMetadata(
                "Pithos profile contains duplicate contig IDs".to_string(),
            ));
        }
        let total_bases = self
            .metadata
            .contigs
            .iter()
            .try_fold(0u64, |total, contig| total.checked_add(contig.length))
            .ok_or_else(|| {
                ArchiveError::CorruptMetadata(
                    "Pithos contig lengths overflow total base count".to_string(),
                )
            })?;
        if total_bases != self.metadata.total_bases {
            return Err(ArchiveError::CorruptMetadata(format!(
                "Pithos total base count {} does not match contigs {}",
                self.metadata.total_bases, total_bases
            )));
        }
        let mut scheme_ids = self
            .seed_schemes
            .iter()
            .map(|scheme| scheme.scheme_id)
            .collect::<Vec<_>>();
        scheme_ids.sort_unstable();
        if scheme_ids.windows(2).any(|window| window[0] == window[1]) {
            return Err(ArchiveError::CorruptMetadata(
                "Pithos profile contains duplicate seed scheme IDs".to_string(),
            ));
        }
        Ok(())
    }

    fn contig(&self, id: u32) -> Option<&ArchiveContig> {
        self.metadata.contigs.iter().find(|contig| contig.id == id)
    }

    fn scheme(&self, id: SeedSchemeId) -> Option<&SeedSchemeDescriptor> {
        self.seed_schemes
            .iter()
            .find(|scheme| scheme.scheme_id == id.0)
    }
}

/// A source-side request.  Coordinates are zero-based half-open contig
/// coordinates; reverse-complement requests retain the forward coordinate
/// range and reverse only the returned bases.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct PithosSequenceRequest {
    pub contig_id: u32,
    pub start: u64,
    pub end: u64,
    pub reverse_complement: bool,
}

impl PithosSequenceRequest {
    #[must_use]
    pub const fn len(self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start == self.end
    }
}

impl From<SequenceRequest> for PithosSequenceRequest {
    fn from(request: SequenceRequest) -> Self {
        Self {
            contig_id: request.contig_id,
            start: request.start,
            end: request.end,
            reverse_complement: request.reverse_complement,
        }
    }
}

/// A sequence range returned by a Pithos source.
///
/// `checksum` is SHA-256 over the decoded bases in the returned orientation.
/// It is intentionally part of the source contract so mmap and remote
/// implementations cannot silently skip block integrity checks.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PithosSequenceSlice {
    pub request: PithosSequenceRequest,
    pub bases: Vec<u8>,
    pub checksum: [u8; 32],
}

impl PithosSequenceSlice {
    #[must_use]
    pub fn from_bases(request: PithosSequenceRequest, bases: Vec<u8>) -> Self {
        let checksum = sha256(&bases);
        Self {
            request,
            bases,
            checksum,
        }
    }
}

/// A seed occurrence returned by the Pithos index after a digest candidate.
/// The source must return the exact packed verification material with each
/// hit.  The adapter checks it again because a digest collision is allowed to
/// increase lookup work but must never become an anchor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PithosSeedMatch {
    pub key_index: usize,
    pub digest: u64,
    pub verification: Vec<u8>,
    pub occurrences: Vec<SeedOccurrence>,
}

/// Result of a source-side seed lookup.  The byte/page counters are optional
/// source measurements; the adapter also accepts sources that expose the same
/// values through [`PithosSourceMetrics`].
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct PithosSeedLookup {
    pub matches: Vec<PithosSeedMatch>,
    pub pages_read: u64,
    pub bytes_read: u64,
}

/// Cumulative counters supplied by the Pithos mmap or remote-range source.
/// `resource` is copied into the backend-neutral metrics without changing its
/// identity or checksum fields.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PithosSourceMetrics {
    pub resource: ResourceMetrics,
    pub mapped_bytes: u64,
    pub resident_bytes: u64,
    pub metadata_bytes_read: u64,
    pub seed_bytes_read: u64,
    pub sequence_bytes_read: u64,
    pub decoded_sequence_bases: u64,
    pub coalesced_range_requests: u64,
}

/// Source contract implemented by the Pithos profile reader.
///
/// The implementation in the Pithos checkout is expected to bind these
/// operations to `/metadata/jam.cbor`, `/index/*`, and either per-contig or
/// packed sequence members.  Keeping this trait in jam-rs avoids making the
/// backend-neutral trace contract depend on Pithos' container model.
pub trait PithosBiosequenceSource: Send + Sync {
    fn profile(&self) -> ArchiveResult<PithosProfile>;

    fn lookup_seed_keys(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<PithosSeedLookup>;

    fn read_sequence_ranges(
        &self,
        requests: &[PithosSequenceRequest],
    ) -> ArchiveResult<Vec<PithosSequenceSlice>>;

    fn metrics(&self) -> PithosSourceMetrics;
}

/// Controls integrity checks at the Jam/Pithos boundary.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PithosArchiveOptions {
    /// Require every returned range to carry and pass a decoded-byte checksum.
    pub require_checksums: bool,
}

impl Default for PithosArchiveOptions {
    fn default() -> Self {
        Self {
            require_checksums: true,
        }
    }
}

#[derive(Debug, Default)]
struct AdapterCounters {
    fallback_seed_bytes: AtomicU64,
    coalesced_ranges: AtomicU64,
}

/// `TraceArchive` implementation backed by a Pithos biosequence profile.
pub struct PithosBiosequenceArchive<S> {
    source: S,
    profile: PithosProfile,
    options: PithosArchiveOptions,
    counters: AdapterCounters,
}

impl<S: PithosBiosequenceSource> PithosBiosequenceArchive<S> {
    pub fn new(source: S) -> ArchiveResult<Self> {
        Self::with_options(source, PithosArchiveOptions::default())
    }

    pub fn with_options(source: S, options: PithosArchiveOptions) -> ArchiveResult<Self> {
        let profile = source.profile()?;
        profile.validate()?;
        Ok(Self {
            source,
            profile,
            options,
            counters: AdapterCounters::default(),
        })
    }

    #[must_use]
    pub fn profile(&self) -> &PithosProfile {
        &self.profile
    }

    #[must_use]
    pub fn sequence_organization(&self) -> PithosSequenceOrganization {
        self.profile.sequence_organization
    }

    #[must_use]
    pub fn source(&self) -> &S {
        &self.source
    }

    fn validate_request(&self, request: SequenceRequest) -> ArchiveResult<()> {
        if request.start > request.end {
            return Err(ArchiveError::InvalidSequenceRange {
                start: request.start,
                end: request.end,
            });
        }
        let contig = self
            .profile
            .contig(request.contig_id)
            .ok_or(ArchiveError::UnknownContig(request.contig_id))?;
        if request.end > contig.length {
            return Err(ArchiveError::InvalidSequenceRange {
                start: request.start,
                end: request.end,
            });
        }
        Ok(())
    }

    fn validate_occurrence(&self, occurrence: SeedOccurrence) -> ArchiveResult<()> {
        let contig = self
            .profile
            .contig(occurrence.contig_id)
            .ok_or(ArchiveError::UnknownContig(occurrence.contig_id))?;
        if occurrence.span == 0
            || occurrence
                .position
                .checked_add(u64::from(occurrence.span))
                .is_none_or(|end| end > contig.length)
        {
            return Err(ArchiveError::CorruptMetadata(format!(
                "Pithos occurrence exceeds contig {}",
                occurrence.contig_id
            )));
        }
        Ok(())
    }

    fn archive_metrics(&self) -> ArchiveMetrics {
        let source = self.source.metrics();
        let fallback_seed_bytes = self.counters.fallback_seed_bytes.load(Ordering::Relaxed);
        let coalesced_ranges = self.counters.coalesced_ranges.load(Ordering::Relaxed);
        ArchiveMetrics {
            resource: source.resource,
            mapped_bytes: source.mapped_bytes,
            resident_bytes: source.resident_bytes,
            metadata_bytes_read: source.metadata_bytes_read,
            seed_bytes_read: source.seed_bytes_read.saturating_add(fallback_seed_bytes),
            sequence_bytes_read: source.sequence_bytes_read,
            decoded_sequence_bases: source.decoded_sequence_bases,
            coalesced_range_requests: source
                .coalesced_range_requests
                .saturating_add(coalesced_ranges),
        }
    }
}

impl<S: PithosBiosequenceSource> TraceArchive for PithosBiosequenceArchive<S> {
    fn metadata(&self) -> ArchiveResult<ArchiveMetadata> {
        Ok(self.profile.metadata.clone())
    }

    fn available_seed_schemes(&self) -> ArchiveResult<Vec<SeedSchemeDescriptor>> {
        Ok(self.profile.seed_schemes.clone())
    }

    fn lookup_seeds(
        &self,
        scheme: SeedSchemeId,
        keys: &[SeedKey],
    ) -> ArchiveResult<crate::archive::SeedLookupResult> {
        if self.profile.scheme(scheme).is_none() {
            return Err(ArchiveError::UnknownSeedScheme(scheme.0));
        }
        let before = self.source.metrics();
        let lookup = self.source.lookup_seed_keys(scheme, keys)?;
        let after = self.source.metrics();
        let source_seed_delta = after.seed_bytes_read.saturating_sub(before.seed_bytes_read);
        if source_seed_delta == 0 && lookup.bytes_read != 0 {
            self.counters
                .fallback_seed_bytes
                .fetch_add(lookup.bytes_read, Ordering::Relaxed);
        }
        let mut occurrences_by_key = HashMap::<usize, Vec<SeedOccurrence>>::new();
        let mut occurrences_before = 0u64;
        let mut occurrences_after = 0u64;
        for hit in lookup.matches {
            if hit.key_index >= keys.len() {
                return Err(ArchiveError::CorruptMetadata(format!(
                    "Pithos seed hit key index {} exceeds {} query keys",
                    hit.key_index,
                    keys.len()
                )));
            }
            let key = &keys[hit.key_index];
            if hit.digest != key.digest || hit.verification != key.verification {
                // A digest collision or a malformed source record is never an
                // exact key match.  It is intentionally ignored rather than
                // promoted to an anchor.
                continue;
            }
            occurrences_before = occurrences_before
                .checked_add(u64::try_from(hit.occurrences.len()).map_err(|_| {
                    ArchiveError::Backend("Pithos occurrence count overflows u64".to_string())
                })?)
                .ok_or_else(|| {
                    ArchiveError::Backend("Pithos occurrence count overflows u64".to_string())
                })?;
            let occurrences = occurrences_by_key.entry(hit.key_index).or_default();
            for occurrence in hit.occurrences {
                self.validate_occurrence(occurrence)?;
                occurrences.push(occurrence);
            }
        }
        let mut matches = Vec::with_capacity(occurrences_by_key.len());
        for (key_index, mut occurrences) in occurrences_by_key {
            occurrences.sort_unstable_by_key(|occurrence| {
                (
                    occurrence.contig_id,
                    occurrence.position,
                    occurrence.span,
                    occurrence.reverse,
                )
            });
            occurrences.dedup();
            occurrences_after = occurrences_after
                .checked_add(u64::try_from(occurrences.len()).map_err(|_| {
                    ArchiveError::Backend("Pithos occurrence count overflows u64".to_string())
                })?)
                .ok_or_else(|| {
                    ArchiveError::Backend("Pithos occurrence count overflows u64".to_string())
                })?;
            if !occurrences.is_empty() {
                matches.push(SeedMatch {
                    key_index,
                    occurrences,
                });
            }
        }
        matches.sort_unstable_by_key(|seed_match| seed_match.key_index);
        Ok(crate::archive::SeedLookupResult {
            matches,
            metrics: SeedLookupMetrics {
                pages_read: lookup.pages_read,
                keys_tested: u64::try_from(keys.len()).unwrap_or(u64::MAX),
                occurrences_before_limits: occurrences_before,
                occurrences_after_limits: occurrences_after,
                bytes_read: lookup
                    .bytes_read
                    .max(after.seed_bytes_read.saturating_sub(before.seed_bytes_read)),
            },
        })
    }

    fn read_sequences(&self, requests: &[SequenceRequest]) -> ArchiveResult<Vec<SequenceSlice>> {
        for request in requests {
            self.validate_request(*request)?;
        }
        let mut output = requests
            .iter()
            .copied()
            .map(|request| {
                request.is_empty().then_some(SequenceSlice {
                    request,
                    bases: Vec::new(),
                })
            })
            .collect::<Vec<_>>();
        let mut order = requests
            .iter()
            .enumerate()
            .filter(|(_, request)| !request.is_empty())
            .map(|(index, request)| (index, *request))
            .collect::<Vec<_>>();
        order.sort_unstable_by_key(|(index, request)| {
            (
                request.contig_id,
                request.reverse_complement,
                request.start,
                request.end,
                *index,
            )
        });

        let mut merged = Vec::new();
        let mut members = Vec::<Vec<usize>>::new();
        for (index, request) in order {
            let can_extend = merged
                .last()
                .is_some_and(|merged_request: &PithosSequenceRequest| {
                    merged_request.contig_id == request.contig_id
                        && merged_request.reverse_complement == request.reverse_complement
                        && request.start <= merged_request.end
                });
            if can_extend {
                let current = merged
                    .last_mut()
                    .expect("can_extend implies a merged request exists");
                current.end = current.end.max(request.end);
                members
                    .last_mut()
                    .expect("can_extend implies a member group exists")
                    .push(index);
            } else {
                merged.push(PithosSequenceRequest::from(request));
                members.push(vec![index]);
            }
        }
        if merged.is_empty() {
            return output
                .into_iter()
                .map(|slice| {
                    slice.ok_or_else(|| {
                        ArchiveError::Backend(
                            "Pithos empty sequence request was not processed".to_string(),
                        )
                    })
                })
                .collect();
        }

        let before = self.source.metrics();
        let slices = self.source.read_sequence_ranges(&merged)?;
        let after = self.source.metrics();
        if after.coalesced_range_requests == before.coalesced_range_requests {
            self.counters.coalesced_ranges.fetch_add(
                u64::try_from(merged.len()).unwrap_or(u64::MAX),
                Ordering::Relaxed,
            );
        }
        let mut by_request = HashMap::with_capacity(slices.len());
        for slice in slices {
            if by_request.insert(slice.request, slice).is_some() {
                return Err(ArchiveError::CorruptMetadata(
                    "Pithos source returned duplicate sequence ranges".to_string(),
                ));
            }
        }
        for (merged_request, group) in merged.iter().zip(members.iter()) {
            let slice = by_request.remove(merged_request).ok_or_else(|| {
                ArchiveError::Backend(format!(
                    "Pithos source omitted sequence range contig={} start={} end={}",
                    merged_request.contig_id, merged_request.start, merged_request.end
                ))
            })?;
            validate_sequence_slice(&self.profile, &slice, self.options.require_checksums)?;
            for &index in group {
                let request = requests[index];
                let start = if request.reverse_complement {
                    merged_request.end.checked_sub(request.end).ok_or({
                        ArchiveError::InvalidSequenceRange {
                            start: request.start,
                            end: request.end,
                        }
                    })?
                } else {
                    request.start.checked_sub(merged_request.start).ok_or(
                        ArchiveError::InvalidSequenceRange {
                            start: request.start,
                            end: request.end,
                        },
                    )?
                };
                let end = if request.reverse_complement {
                    merged_request.end.checked_sub(request.start).ok_or(
                        ArchiveError::InvalidSequenceRange {
                            start: request.start,
                            end: request.end,
                        },
                    )?
                } else {
                    request.end.checked_sub(merged_request.start).ok_or(
                        ArchiveError::InvalidSequenceRange {
                            start: request.start,
                            end: request.end,
                        },
                    )?
                };
                let start = usize::try_from(start).map_err(|_| {
                    ArchiveError::Backend(
                        "Pithos sequence offset exceeds platform usize".to_string(),
                    )
                })?;
                let end = usize::try_from(end).map_err(|_| {
                    ArchiveError::Backend(
                        "Pithos sequence offset exceeds platform usize".to_string(),
                    )
                })?;
                let bases = slice.bases.get(start..end).ok_or_else(|| {
                    ArchiveError::CorruptMetadata(
                        "Pithos sequence source returned an invalid range".to_string(),
                    )
                })?;
                output[index] = Some(SequenceSlice {
                    request,
                    bases: bases.to_vec(),
                });
            }
        }
        if let Some(extra) = by_request.keys().next().copied() {
            return Err(ArchiveError::CorruptMetadata(format!(
                "Pithos source returned unexpected sequence range contig={} start={} end={}",
                extra.contig_id, extra.start, extra.end
            )));
        }
        output
            .into_iter()
            .map(|slice| {
                slice.ok_or_else(|| {
                    ArchiveError::Backend("Pithos sequence request was not processed".to_string())
                })
            })
            .collect()
    }

    fn metrics(&self) -> ArchiveMetrics {
        self.archive_metrics()
    }
}

fn validate_sequence_slice(
    profile: &PithosProfile,
    slice: &PithosSequenceSlice,
    require_checksum: bool,
) -> ArchiveResult<()> {
    let contig = profile
        .contig(slice.request.contig_id)
        .ok_or(ArchiveError::UnknownContig(slice.request.contig_id))?;
    if slice.request.start > slice.request.end || slice.request.end > contig.length {
        return Err(ArchiveError::InvalidSequenceRange {
            start: slice.request.start,
            end: slice.request.end,
        });
    }
    let expected_len = usize::try_from(slice.request.len()).map_err(|_| {
        ArchiveError::Backend("Pithos sequence length exceeds platform usize".to_string())
    })?;
    if slice.bases.len() != expected_len {
        return Err(ArchiveError::CorruptMetadata(format!(
            "Pithos sequence range has {} decoded bases, expected {}",
            slice.bases.len(),
            expected_len
        )));
    }
    for &base in &slice.bases {
        if encode_base(base).is_none() && normalize_ambiguity(base).is_err() {
            return Err(ArchiveError::CorruptMetadata(format!(
                "Pithos sequence contains unsupported base byte {}",
                base
            )));
        }
    }
    if require_checksum {
        let actual = sha256(&slice.bases);
        if actual != slice.checksum {
            return Err(ArchiveError::ChecksumMismatch(format!(
                "Pithos sequence contig {}:{}..{}",
                slice.request.contig_id, slice.request.start, slice.request.end
            )));
        }
    }
    Ok(())
}

fn sha256(bytes: &[u8]) -> [u8; 32] {
    let digest = Sha256::digest(bytes);
    digest.into()
}

/// Reverse-complement helper for source implementations that expose forward
/// slab bytes but not a direct reverse-complement read.
pub fn reverse_complement_in_place(bases: &mut [u8]) -> ArchiveResult<()> {
    bases.reverse();
    for base in bases {
        *base = complement_base(*base).map_err(|error| ArchiveError::Backend(error.to_string()))?;
    }
    Ok(())
}
