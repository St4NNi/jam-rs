use crate::jidx::{
    CONTIG_POSTING_SIZE, DOCUMENT_POSTING_SIZE, Header, JidxError, SEED_RECORD_SIZE, SectionKind,
    read_u32, read_u64,
};
use crate::jidx_reader::JidxReader;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedEntry {
    pub packed_key: u64,
    pub document_frequency: u32,
    pub occurrence_count: u64,
    pub(crate) document_offset: u64,
    pub(crate) occurrence_offset: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedOccurrence {
    pub contig_id: u32,
    pub position: u64,
    pub canonical_orientation: bool,
}

pub(crate) fn lookup(reader: &JidxReader, packed_key: u64) -> Result<Option<SeedEntry>, JidxError> {
    let header = reader.header();
    validate_packed_key(packed_key, header.k)?;
    let mut low = 0;
    let mut high = header.seed_count;
    while low < high {
        let middle = low + (high - low) / 2;
        let record = seed_record(reader, middle)?;
        match record.packed_key.cmp(&packed_key) {
            std::cmp::Ordering::Less => low = middle + 1,
            std::cmp::Ordering::Greater => high = middle,
            std::cmp::Ordering::Equal => {
                validate_record(header, record)?;
                return Ok(Some(record.into()));
            }
        }
    }
    Ok(None)
}

pub(crate) fn validate_table(reader: &JidxReader) -> Result<(), JidxError> {
    let header = reader.header();
    let mut previous = None;
    let mut document_offset = 0;
    let mut occurrence_offset = 0;
    for index in 0..header.seed_count {
        let record = seed_record(reader, index)?;
        validate_record(header, record)?;
        if record.document_offset != document_offset
            || record.occurrence_offset != occurrence_offset
        {
            return Err(JidxError::Invalid("posting order"));
        }
        validate_document_bytes(document_bytes(reader, record.into())?, header)?;
        validate_occurrence_bytes(occurrence_bytes(reader, record.into())?, header)?;
        if previous.is_some_and(|key| key >= record.packed_key) {
            return Err(JidxError::Invalid("seed order"));
        }
        document_offset = document_offset
            .checked_add(u64::from(record.document_count) * u64::from(DOCUMENT_POSTING_SIZE))
            .ok_or(JidxError::Invalid("document posting length"))?;
        occurrence_offset = occurrence_offset
            .checked_add(record.occurrence_count * u64::from(CONTIG_POSTING_SIZE))
            .ok_or(JidxError::Invalid("contig posting length"))?;
        previous = Some(record.packed_key);
    }
    if document_offset != header.section(SectionKind::DocumentPostings).length
        || occurrence_offset != header.section(SectionKind::ContigPostings).length
    {
        return Err(JidxError::Invalid("posting coverage"));
    }
    Ok(())
}

pub(crate) fn entry(reader: &JidxReader, index: u64) -> Result<SeedEntry, JidxError> {
    let header = reader.header();
    let record = seed_record(reader, index)?;
    validate_record(header, record)?;
    Ok(record.into())
}

pub(crate) fn documents(reader: &JidxReader, seed: SeedEntry) -> Result<Vec<u32>, JidxError> {
    let header = reader.header();
    let bytes = document_bytes(reader, seed)?;
    validate_document_bytes(bytes, header)?;
    Ok(bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|bytes| u32::from_le_bytes(*bytes))
        .collect())
}

pub(crate) fn occurrences(
    reader: &JidxReader,
    seed: SeedEntry,
) -> Result<Vec<SeedOccurrence>, JidxError> {
    let header = reader.header();
    let bytes = occurrence_bytes(reader, seed)?;
    validate_occurrence_bytes(bytes, header)?;
    Ok(bytes
        .as_chunks::<16>()
        .0
        .iter()
        .map(|bytes| SeedOccurrence {
            contig_id: read_u32(bytes, 0),
            canonical_orientation: bytes[4] == 1,
            position: read_u64(bytes, 8),
        })
        .collect())
}

#[derive(Clone, Copy)]
struct SeedRecord {
    packed_key: u64,
    document_offset: u64,
    document_count: u32,
    occurrence_offset: u64,
    occurrence_count: u64,
}

impl From<SeedRecord> for SeedEntry {
    fn from(record: SeedRecord) -> Self {
        Self {
            packed_key: record.packed_key,
            document_frequency: record.document_count,
            occurrence_count: record.occurrence_count,
            document_offset: record.document_offset,
            occurrence_offset: record.occurrence_offset,
        }
    }
}

fn seed_record(reader: &JidxReader, index: u64) -> Result<SeedRecord, JidxError> {
    let header = reader.header();
    let section = header.section(SectionKind::Seeds);
    let start = section
        .offset
        .checked_add(
            index
                .checked_mul(u64::from(SEED_RECORD_SIZE))
                .ok_or(JidxError::Invalid("seed offset"))?,
        )
        .ok_or(JidxError::Invalid("seed offset"))?;
    let end = start
        .checked_add(u64::from(SEED_RECORD_SIZE))
        .ok_or(JidxError::Invalid("seed range"))?;
    let bytes = reader.checked_bytes(start, end)?;
    if read_u32(bytes, 20) != 0 {
        return Err(JidxError::Invalid("seed reservation"));
    }
    Ok(SeedRecord {
        packed_key: read_u64(bytes, 0),
        document_offset: read_u64(bytes, 8),
        document_count: read_u32(bytes, 16),
        occurrence_offset: read_u64(bytes, 24),
        occurrence_count: read_u64(bytes, 32),
    })
}

fn validate_record(header: &Header, record: SeedRecord) -> Result<(), JidxError> {
    validate_packed_key(record.packed_key, header.k)?;
    if record.document_count == 0
        || record.document_count > header.document_count
        || record.occurrence_count == 0
        || record.occurrence_count > header.occurrence_count
        || !record
            .document_offset
            .is_multiple_of(u64::from(DOCUMENT_POSTING_SIZE))
        || !record
            .occurrence_offset
            .is_multiple_of(u64::from(CONTIG_POSTING_SIZE))
    {
        return Err(JidxError::Invalid("seed record"));
    }
    let document_end = record
        .document_offset
        .checked_add(
            u64::from(record.document_count)
                .checked_mul(u64::from(DOCUMENT_POSTING_SIZE))
                .ok_or(JidxError::Invalid("document posting range"))?,
        )
        .ok_or(JidxError::Invalid("document posting range"))?;
    let occurrence_end = record
        .occurrence_offset
        .checked_add(
            record
                .occurrence_count
                .checked_mul(u64::from(CONTIG_POSTING_SIZE))
                .ok_or(JidxError::Invalid("contig posting range"))?,
        )
        .ok_or(JidxError::Invalid("contig posting range"))?;
    if document_end > header.section(SectionKind::DocumentPostings).length
        || occurrence_end > header.section(SectionKind::ContigPostings).length
    {
        return Err(JidxError::Invalid("seed posting range"));
    }
    Ok(())
}

fn document_bytes(reader: &JidxReader, seed: SeedEntry) -> Result<&[u8], JidxError> {
    let header = reader.header();
    let section = header.section(SectionKind::DocumentPostings);
    let start = section
        .offset
        .checked_add(seed.document_offset)
        .ok_or(JidxError::Invalid("document posting offset"))?;
    let length = u64::from(seed.document_frequency)
        .checked_mul(u64::from(DOCUMENT_POSTING_SIZE))
        .ok_or(JidxError::Invalid("document posting length"))?;
    let end = start
        .checked_add(length)
        .ok_or(JidxError::Invalid("document posting range"))?;
    reader.checked_bytes(start, end)
}

fn occurrence_bytes(reader: &JidxReader, seed: SeedEntry) -> Result<&[u8], JidxError> {
    let header = reader.header();
    let section = header.section(SectionKind::ContigPostings);
    let start = section
        .offset
        .checked_add(seed.occurrence_offset)
        .ok_or(JidxError::Invalid("contig posting offset"))?;
    let length = seed
        .occurrence_count
        .checked_mul(u64::from(CONTIG_POSTING_SIZE))
        .ok_or(JidxError::Invalid("contig posting length"))?;
    let end = start
        .checked_add(length)
        .ok_or(JidxError::Invalid("contig posting range"))?;
    reader.checked_bytes(start, end)
}

fn validate_document_bytes(bytes: &[u8], header: &Header) -> Result<(), JidxError> {
    let mut previous = None;
    for bytes in bytes.as_chunks::<4>().0 {
        let id = u32::from_le_bytes(*bytes);
        if id >= header.document_count || previous.is_some_and(|previous| previous >= id) {
            return Err(JidxError::Invalid("document postings"));
        }
        previous = Some(id);
    }
    Ok(())
}

fn validate_occurrence_bytes(bytes: &[u8], header: &Header) -> Result<(), JidxError> {
    let mut previous = None;
    for bytes in bytes.as_chunks::<16>().0 {
        let contig_id = read_u32(bytes, 0);
        let position = read_u64(bytes, 8);
        if contig_id >= header.contig_count
            || bytes[4] > 1
            || bytes[5..8].iter().any(|byte| *byte != 0)
            || previous.is_some_and(|previous| previous >= (contig_id, position))
        {
            return Err(JidxError::Invalid("contig postings"));
        }
        previous = Some((contig_id, position));
    }
    Ok(())
}

pub(crate) fn validate_packed_key(key: u64, k: u8) -> Result<(), JidxError> {
    if k < 32 && key >= 1u64 << (2 * k) {
        return Err(JidxError::Invalid("packed seed"));
    }
    Ok(())
}

pub(crate) fn encode_seed(seed: SeedEntry) -> [u8; SEED_RECORD_SIZE as usize] {
    let mut bytes = [0; SEED_RECORD_SIZE as usize];
    crate::jidx::put_u64(&mut bytes, 0, seed.packed_key);
    crate::jidx::put_u64(&mut bytes, 8, seed.document_offset);
    crate::jidx::put_u32(&mut bytes, 16, seed.document_frequency);
    crate::jidx::put_u64(&mut bytes, 24, seed.occurrence_offset);
    crate::jidx::put_u64(&mut bytes, 32, seed.occurrence_count);
    bytes
}

pub(crate) fn encode_occurrence(occurrence: SeedOccurrence) -> [u8; CONTIG_POSTING_SIZE as usize] {
    let mut bytes = [0; CONTIG_POSTING_SIZE as usize];
    crate::jidx::put_u32(&mut bytes, 0, occurrence.contig_id);
    bytes[4] = u8::from(occurrence.canonical_orientation);
    crate::jidx::put_u64(&mut bytes, 8, occurrence.position);
    bytes
}
