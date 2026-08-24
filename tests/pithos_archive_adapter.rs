use jam_rs::archive::pithos::{
    PithosBiosequenceArchive, PithosBiosequenceSource, PithosProfile, PithosSeedLookup,
    PithosSeedMatch, PithosSequenceOrganization, PithosSequenceRequest, PithosSequenceSlice,
    PithosSourceMetrics,
};
use jam_rs::archive::{
    ArchiveContig, ArchiveMetadata, SeedKey, SeedOccurrence, SeedSchemeDescriptor, SeedSchemeId,
    SequenceRequest, TraceArchive,
};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

const SCHEME: SeedSchemeDescriptor = SeedSchemeDescriptor {
    scheme_id: 7,
    algorithm_id: 1,
    span: 31,
    informative_bases: 31,
    density_parameter: 1,
    bucket_bits: 8,
    key_encoding: 1,
    occurrence_encoding: 1,
    flags: 0,
};

#[derive(Clone)]
struct MockSource {
    profile: PithosProfile,
    contigs: HashMap<u32, Vec<u8>>,
    hits: Vec<PithosSeedMatch>,
    calls: Arc<Mutex<Vec<Vec<PithosSequenceRequest>>>>,
    metrics: Arc<Mutex<PithosSourceMetrics>>,
    corrupt_checksum: bool,
    extra_ranges: Vec<PithosSequenceRequest>,
}

impl MockSource {
    fn new(organization: PithosSequenceOrganization) -> Self {
        let first = b"ACGTNACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let second = b"TTGCAACCGGTTNNAACCGGTTCCGGAA".to_vec();
        let contigs = HashMap::from([(0, first), (1, second)]);
        let total_bases = contigs.values().map(|sequence| sequence.len() as u64).sum();
        let metadata = ArchiveMetadata {
            format_identifier: "pithos-biosequence".to_string(),
            format_version: 1,
            layout_identifier: 0x5042_494f,
            source_assembly_sha256: [3; 32],
            archive_sha256: Some([4; 32]),
            builder_version: "test".to_string(),
            source_commit: "test".to_string(),
            hash_algorithm: "jamhash_u64_v1".to_string(),
            total_bases,
            contigs: vec![
                ArchiveContig {
                    id: 0,
                    name: "one".to_string(),
                    length: contigs[&0].len() as u64,
                },
                ArchiveContig {
                    id: 1,
                    name: "two".to_string(),
                    length: contigs[&1].len() as u64,
                },
            ],
        };
        Self {
            profile: PithosProfile {
                metadata,
                seed_schemes: vec![SCHEME],
                sequence_organization: organization,
            },
            contigs,
            hits: Vec::new(),
            calls: Arc::new(Mutex::new(Vec::new())),
            metrics: Arc::new(Mutex::new(PithosSourceMetrics::default())),
            corrupt_checksum: false,
            extra_ranges: Vec::new(),
        }
    }

    fn with_hit(mut self, hit: PithosSeedMatch) -> Self {
        self.hits.push(hit);
        self
    }

    fn with_corrupt_checksum(mut self) -> Self {
        self.corrupt_checksum = true;
        self
    }

    fn with_extra_range(mut self, request: PithosSequenceRequest) -> Self {
        self.extra_ranges.push(request);
        self
    }
}

impl PithosBiosequenceSource for MockSource {
    fn profile(&self) -> jam_rs::archive::ArchiveResult<PithosProfile> {
        Ok(self.profile.clone())
    }

    fn lookup_seed_keys(
        &self,
        _scheme: SeedSchemeId,
        _keys: &[SeedKey],
    ) -> jam_rs::archive::ArchiveResult<PithosSeedLookup> {
        Ok(PithosSeedLookup {
            matches: self.hits.clone(),
            pages_read: self.hits.len() as u64,
            bytes_read: (self.hits.len() * 64) as u64,
        })
    }

    fn read_sequence_ranges(
        &self,
        requests: &[PithosSequenceRequest],
    ) -> jam_rs::archive::ArchiveResult<Vec<PithosSequenceSlice>> {
        self.calls.lock().unwrap().push(requests.to_vec());
        let mut metrics = self.metrics.lock().unwrap();
        metrics.resource.range_requests += (requests.len() + self.extra_ranges.len()) as u64;
        let mut output = Vec::with_capacity(requests.len() + self.extra_ranges.len());
        for request in requests.iter().chain(self.extra_ranges.iter()) {
            let sequence = &self.contigs[&request.contig_id];
            let mut bases = sequence[request.start as usize..request.end as usize].to_vec();
            if request.reverse_complement {
                jam_rs::archive::pithos::reverse_complement_in_place(&mut bases)?;
            }
            metrics.sequence_bytes_read += bases.len() as u64;
            metrics.decoded_sequence_bases += bases.len() as u64;
            metrics.resource.requested_bytes += bases.len() as u64;
            metrics.resource.returned_bytes += bases.len() as u64;
            let mut slice = PithosSequenceSlice::from_bases(*request, bases);
            if self.corrupt_checksum {
                slice.checksum[0] ^= 1;
            }
            output.push(slice);
        }
        Ok(output)
    }

    fn metrics(&self) -> PithosSourceMetrics {
        *self.metrics.lock().unwrap()
    }
}

fn key() -> SeedKey {
    SeedKey {
        digest: 0x1234,
        verification: vec![1, 2, 3, 4],
    }
}

fn matching_hit(occurrence: SeedOccurrence) -> PithosSeedMatch {
    PithosSeedMatch {
        key_index: 0,
        digest: key().digest,
        verification: key().verification,
        occurrences: vec![occurrence],
    }
}

#[test]
fn per_contig_and_packed_slab_sources_return_the_same_logical_results() {
    let occurrence = SeedOccurrence {
        contig_id: 0,
        position: 2,
        span: 4,
        reverse: false,
    };
    let per_contig = MockSource::new(PithosSequenceOrganization::OneMemberPerContig)
        .with_hit(matching_hit(occurrence));
    let packed =
        MockSource::new(PithosSequenceOrganization::PackedSlabs).with_hit(matching_hit(occurrence));
    let per_contig_archive = PithosBiosequenceArchive::new(per_contig.clone()).unwrap();
    let packed_archive = PithosBiosequenceArchive::new(packed.clone()).unwrap();

    let per_lookup = per_contig_archive
        .lookup_seeds(SeedSchemeId(SCHEME.scheme_id), &[key()])
        .unwrap();
    let packed_lookup = packed_archive
        .lookup_seeds(SeedSchemeId(SCHEME.scheme_id), &[key()])
        .unwrap();
    assert_eq!(per_lookup.matches, packed_lookup.matches);

    let requests = [
        SequenceRequest::new(0, 4, 15, false).unwrap(),
        SequenceRequest::new(0, 6, 13, true).unwrap(),
    ];
    let per_slices = per_contig_archive.read_sequences(&requests).unwrap();
    let packed_slices = packed_archive.read_sequences(&requests).unwrap();
    assert_eq!(per_slices, packed_slices);
    assert_eq!(
        per_contig_archive.profile().sequence_organization,
        PithosSequenceOrganization::OneMemberPerContig
    );
    assert_eq!(
        packed_archive.profile().sequence_organization,
        PithosSequenceOrganization::PackedSlabs
    );
}

#[test]
fn digest_collisions_never_create_seed_matches() {
    let wrong = PithosSeedMatch {
        key_index: 0,
        digest: key().digest,
        verification: vec![9, 9, 9, 9],
        occurrences: vec![SeedOccurrence {
            contig_id: 0,
            position: 0,
            span: 4,
            reverse: false,
        }],
    };
    let source = MockSource::new(PithosSequenceOrganization::PackedSlabs)
        .with_hit(wrong)
        .with_hit(matching_hit(SeedOccurrence {
            contig_id: 0,
            position: 4,
            span: 4,
            reverse: false,
        }));
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let result = archive
        .lookup_seeds(SeedSchemeId(SCHEME.scheme_id), &[key()])
        .unwrap();
    assert_eq!(result.matches.len(), 1);
    assert_eq!(result.matches[0].occurrences[0].position, 4);
    assert_eq!(result.metrics.occurrences_before_limits, 1);
}

#[test]
fn verified_hits_for_one_key_are_aggregated_and_deduplicated() {
    let source = MockSource::new(PithosSequenceOrganization::PackedSlabs)
        .with_hit(matching_hit(SeedOccurrence {
            contig_id: 0,
            position: 4,
            span: 4,
            reverse: false,
        }))
        .with_hit(PithosSeedMatch {
            key_index: 0,
            digest: key().digest,
            verification: key().verification,
            occurrences: vec![
                SeedOccurrence {
                    contig_id: 0,
                    position: 4,
                    span: 4,
                    reverse: false,
                },
                SeedOccurrence {
                    contig_id: 0,
                    position: 8,
                    span: 4,
                    reverse: false,
                },
            ],
        });
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let result = archive
        .lookup_seeds(SeedSchemeId(SCHEME.scheme_id), &[key()])
        .unwrap();
    assert_eq!(result.matches.len(), 1);
    assert_eq!(
        result.matches[0]
            .occurrences
            .iter()
            .map(|occurrence| occurrence.position)
            .collect::<Vec<_>>(),
        vec![4, 8]
    );
    assert_eq!(result.metrics.occurrences_before_limits, 3);
    assert_eq!(result.metrics.occurrences_after_limits, 2);
}

#[test]
fn sequence_ranges_are_selected_coalesced_and_never_full_object_reads() {
    let source = MockSource::new(PithosSequenceOrganization::PackedSlabs);
    let calls = source.calls.clone();
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let requests = [
        SequenceRequest::new(0, 5, 11, false).unwrap(),
        SequenceRequest::new(0, 10, 17, false).unwrap(),
    ];
    let slices = archive.read_sequences(&requests).unwrap();
    assert_eq!(slices[0].bases.len(), 6);
    assert_eq!(slices[1].bases.len(), 7);
    let calls = calls.lock().unwrap();
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].len(), 1);
    assert_eq!(calls[0][0].start, 5);
    assert_eq!(calls[0][0].end, 17);
    assert!(calls[0][0].end < 38);
    assert_eq!(archive.metrics().resource.range_requests, 1);
}

#[test]
fn mixed_empty_and_non_empty_sequence_requests_preserve_empty_outputs() {
    let source = MockSource::new(PithosSequenceOrganization::PackedSlabs);
    let calls = source.calls.clone();
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let requests = [
        SequenceRequest::new(0, 7, 7, false).unwrap(),
        SequenceRequest::new(0, 5, 11, false).unwrap(),
    ];
    let slices = archive.read_sequences(&requests).unwrap();
    assert_eq!(slices.len(), requests.len());
    assert_eq!(slices[0].request, requests[0]);
    assert!(slices[0].bases.is_empty());
    assert_eq!(slices[1].bases, b"ACGTAC");
    let calls = calls.lock().unwrap();
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0][0].start, 5);
    assert_eq!(calls[0][0].end, 11);
}

#[test]
fn unexpected_extra_sequence_ranges_are_rejected() {
    let source = MockSource::new(PithosSequenceOrganization::PackedSlabs).with_extra_range(
        PithosSequenceRequest {
            contig_id: 1,
            start: 0,
            end: 1,
            reverse_complement: false,
        },
    );
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let request = SequenceRequest::new(0, 5, 11, false).unwrap();
    let error = archive.read_sequences(&[request]).unwrap_err();
    assert!(matches!(
        error,
        jam_rs::archive::ArchiveError::CorruptMetadata(message)
            if message.contains("unexpected sequence range")
    ));
}

#[test]
fn checksum_failures_are_rejected_before_sequence_projection() {
    let source =
        MockSource::new(PithosSequenceOrganization::OneMemberPerContig).with_corrupt_checksum();
    let archive = PithosBiosequenceArchive::new(source).unwrap();
    let request = SequenceRequest::new(0, 0, 10, false).unwrap();
    assert!(matches!(
        archive.read_sequences(&[request]),
        Err(jam_rs::archive::ArchiveError::ChecksumMismatch(_))
    ));
}
