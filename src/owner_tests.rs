use crate::owner_format::{
    HAS_METADATA, OWNER_CONTIG_SIZE, OWNER_DOCUMENT_SIZE, OWNER_HEADER_SIZE, OWNER_PAGE_SIZE,
    OWNER_SECTION_COUNT, OwnerHeader, OwnerSection, OwnerSectionDescriptor,
};

#[test]
fn owner_header_supports_more_than_two_billion_contigs_without_allocation() {
    let contig_count = 3_000_000_000u32;
    let lengths = [
        0,
        u64::from(OWNER_DOCUMENT_SIZE),
        u64::from(contig_count) * u64::from(OWNER_CONTIG_SIZE),
        0,
        0,
        0,
        0,
    ];
    let mut offset = OWNER_HEADER_SIZE as u64;
    let mut sections = Vec::with_capacity(OWNER_SECTION_COUNT);
    for (kind, length) in OwnerSection::ALL.into_iter().take(7).zip(lengths) {
        offset = offset.next_multiple_of(OWNER_PAGE_SIZE);
        sections.push(OwnerSectionDescriptor {
            kind,
            record_size: kind.record_size(),
            offset,
            length,
        });
        offset = offset.checked_add(length).unwrap();
    }
    offset = offset.next_multiple_of(OWNER_PAGE_SIZE);
    let checksum_length = crate::owner_format::checksum_layout(offset / OWNER_PAGE_SIZE - 1)
        .unwrap()
        .iter()
        .map(|level| level.page_count * OWNER_PAGE_SIZE)
        .sum();
    sections.push(OwnerSectionDescriptor {
        kind: OwnerSection::PageChecksums,
        record_size: OwnerSection::PageChecksums.record_size(),
        offset,
        length: checksum_length,
    });
    let header = OwnerHeader {
        flags: HAS_METADATA,
        k: 15,
        rescue_k15: false,
        minimizer_window: 4,
        owner_ordinal: 0,
        owner_count: 1,
        first_key: 0,
        last_key: u64::MAX,
        key_count: 0,
        occurrence_count: 0,
        document_count: 1,
        contig_count,
        generation_id: [1; 32],
        body_sha256: [2; 32],
        checksum_root_sha256: [3; 32],
        sections: sections.try_into().unwrap(),
    };
    let bytes = header.encode().unwrap();
    let decoded = OwnerHeader::decode(&bytes, offset + checksum_length).unwrap();
    assert_eq!(decoded.contig_count, contig_count);
}
pub(crate) fn random_sequence(mut state: u64, length: usize) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

