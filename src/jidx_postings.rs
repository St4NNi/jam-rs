use crate::jidx::{
    CONTIG_POSTING_SIZE, DOCUMENT_POSTING_SIZE, Header, JidxError, SEED_RECORD_SIZE, SectionKind,
    read_u32, read_u64,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedEntry {
    pub packed_key: u64,
    pub document_frequency: u32,
    pub occurrence_count: u64,
    pub(crate) document_offset: u64,
    pub(crate) occurrence_offset: u64,
}

pub(crate) fn lookup(
    file: &[u8],
    header: &Header,
    packed_key: u64,
) -> Result<Option<SeedEntry>, JidxError> {
    validate_key(packed_key, header.k)?;
    let mut low = 0;
    let mut high = header.seed_count;
    while low < high {
        let middle = low + (high - low) / 2;
        let record = seed_record(file, header, middle)?;
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

pub(crate) fn validate_table(file: &[u8], header: &Header) -> Result<(), JidxError> {
    let mut previous = None;
    for index in 0..header.seed_count {
        let record = seed_record(file, header, index)?;
        validate_record(header, record)?;
        validate_document_bytes(document_bytes(file, header, record.into())?, header)?;
        if previous.is_some_and(|key| key >= record.packed_key) {
            return Err(JidxError::Invalid("seed order"));
        }
        previous = Some(record.packed_key);
    }
    Ok(())
}

pub(crate) fn documents(
    file: &[u8],
    header: &Header,
    seed: SeedEntry,
) -> Result<Vec<u32>, JidxError> {
    let bytes = document_bytes(file, header, seed)?;
    validate_document_bytes(bytes, header)?;
    Ok(bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|bytes| u32::from_le_bytes(*bytes))
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

fn seed_record(file: &[u8], header: &Header, index: u64) -> Result<SeedRecord, JidxError> {
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
    let start = usize::try_from(start).map_err(|_| JidxError::Invalid("seed offset"))?;
    let end = usize::try_from(end).map_err(|_| JidxError::Invalid("seed range"))?;
    let bytes = file
        .get(start..end)
        .ok_or(JidxError::Invalid("seed range"))?;
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
    validate_key(record.packed_key, header.k)?;
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

fn document_bytes<'a>(
    file: &'a [u8],
    header: &Header,
    seed: SeedEntry,
) -> Result<&'a [u8], JidxError> {
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
    let start =
        usize::try_from(start).map_err(|_| JidxError::Invalid("document posting offset"))?;
    let end = usize::try_from(end).map_err(|_| JidxError::Invalid("document posting range"))?;
    file.get(start..end)
        .ok_or(JidxError::Invalid("document posting range"))
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

fn validate_key(key: u64, k: u8) -> Result<(), JidxError> {
    if k < 32 && key >= 1u64 << (2 * k) {
        return Err(JidxError::Invalid("packed seed"));
    }
    Ok(())
}
