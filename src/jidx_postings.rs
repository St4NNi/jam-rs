use crate::jidx::{
    DOCUMENT_POSTING_SIZE, Header, JidxError, SEED_RECORD_SIZE, SectionKind, read_u32, read_u64,
    seed_length,
};
use crate::jidx_reader::JidxReader;

const EXTERNAL_SENTINEL: u32 = u32::MAX;
const EXTERNAL_DOCUMENT_POSTING_SIZE: u64 = 32;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedEntry {
    pub packed_key: u64,
    pub document_frequency: u32,
    pub(crate) document_offset: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedOccurrence {
    pub contig_id: u32,
    pub position: u64,
    pub canonical_orientation: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SeedDocument {
    pub metagenome_id: u32,
    pub occurrence_count: u64,
    seed_key: u64,
    document_id: u32,
    location: DocumentLocation,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum DocumentLocation {
    Inline {
        local_contig: u32,
        payload: u64,
    },
    Delta {
        data_start: u64,
        data_end: u64,
        occurrence_count: u64,
    },
}

struct BatchLookupState {
    low: u64,
    high: u64,
    result: Option<SeedEntry>,
}

pub(crate) fn lookup(reader: &JidxReader, packed_key: u64) -> Result<Option<SeedEntry>, JidxError> {
    let header = reader.header();
    seed_length(header.k, header.rescue_k15, packed_key)?;
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

pub(crate) fn lookup_batch(
    reader: &JidxReader,
    packed_keys: &[u64],
) -> Result<Vec<Option<SeedEntry>>, JidxError> {
    let header = reader.header();
    let mut states = Vec::new();
    states
        .try_reserve_exact(packed_keys.len())
        .map_err(|_| JidxError::Invalid("seed lookup batch"))?;
    for &packed_key in packed_keys {
        seed_length(header.k, header.rescue_k15, packed_key)?;
        states.push(BatchLookupState {
            low: 0,
            high: header.seed_count,
            result: None,
        });
    }

    let mut midpoints = Vec::<(u64, usize)>::new();
    midpoints
        .try_reserve_exact(packed_keys.len())
        .map_err(|_| JidxError::Invalid("seed lookup batch"))?;
    loop {
        midpoints.clear();
        for (index, state) in states.iter().enumerate() {
            if state.low < state.high {
                midpoints.push((state.low + (state.high - state.low) / 2, index));
            }
        }
        if midpoints.is_empty() {
            break;
        }
        midpoints.sort_unstable();

        let mut start = 0usize;
        while start < midpoints.len() {
            let midpoint = midpoints[start].0;
            let record = seed_record(reader, midpoint)?;
            let mut end = start + 1;
            while end < midpoints.len() && midpoints[end].0 == midpoint {
                end += 1;
            }
            for &(_, index) in &midpoints[start..end] {
                let state = &mut states[index];
                match record.packed_key.cmp(&packed_keys[index]) {
                    std::cmp::Ordering::Less => state.low = midpoint + 1,
                    std::cmp::Ordering::Greater => state.high = midpoint,
                    std::cmp::Ordering::Equal => {
                        validate_record(header, record)?;
                        state.result = Some(record.into());
                        state.low = state.high;
                    }
                }
            }
            start = end;
        }
    }

    let mut output = Vec::new();
    output
        .try_reserve_exact(states.len())
        .map_err(|_| JidxError::Invalid("seed lookup batch"))?;
    output.extend(states.into_iter().map(|state| state.result));
    Ok(output)
}

pub(crate) fn validate_table(reader: &JidxReader) -> Result<(), JidxError> {
    let header = reader.header();
    let cold = header.section(SectionKind::ContigPostings);
    let mut previous_key = None;
    let mut document_offset = 0u64;
    let mut cold_offset = 0u64;
    let mut occurrence_count = 0u64;
    for index in 0..header.seed_count {
        let seed = entry(reader, index)?;
        if previous_key.is_some_and(|key| key >= seed.packed_key)
            || seed.document_offset != document_offset
        {
            return Err(JidxError::Invalid("posting order"));
        }
        let (documents, next_document_offset) = documents_with_end(reader, seed)?;
        document_offset = next_document_offset;
        for document in documents {
            occurrence_count = occurrence_count
                .checked_add(document.occurrence_count)
                .ok_or(JidxError::Invalid("occurrence count"))?;
            if let DocumentLocation::Delta {
                data_start,
                data_end,
                ..
            } = document.location
            {
                let expected = cold
                    .offset
                    .checked_add(cold_offset)
                    .ok_or(JidxError::Invalid("contig posting range"))?;
                if data_start != expected {
                    return Err(JidxError::Invalid("contig posting order"));
                }
                cold_offset = data_end
                    .checked_sub(cold.offset)
                    .ok_or(JidxError::Invalid("contig posting range"))?;
            }
        }
        previous_key = Some(seed.packed_key);
    }
    if document_offset != header.section(SectionKind::DocumentPostings).length
        || cold_offset != cold.length
        || occurrence_count != header.occurrence_count
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

pub(crate) fn documents(
    reader: &JidxReader,
    seed: SeedEntry,
) -> Result<Vec<SeedDocument>, JidxError> {
    Ok(documents_with_end(reader, seed)?.0)
}

fn documents_with_end(
    reader: &JidxReader,
    seed: SeedEntry,
) -> Result<(Vec<SeedDocument>, u64), JidxError> {
    let section = reader.header().section(SectionKind::DocumentPostings);
    let mut cursor = section
        .offset
        .checked_add(seed.document_offset)
        .ok_or(JidxError::Invalid("document posting offset"))?;
    let section_end = section
        .offset
        .checked_add(section.length)
        .ok_or(JidxError::Invalid("document posting range"))?;
    let capacity = usize::try_from(seed.document_frequency)
        .map_err(|_| JidxError::Invalid("document frequency"))?;
    let mut documents = Vec::new();
    documents
        .try_reserve_exact(capacity)
        .map_err(|_| JidxError::Invalid("document frequency"))?;
    let mut previous = None;
    let mut seed_occurrence_count = 0u64;
    for _ in 0..seed.document_frequency {
        let inline_end = cursor
            .checked_add(u64::from(DOCUMENT_POSTING_SIZE))
            .ok_or(JidxError::Invalid("document posting range"))?;
        if inline_end > section_end {
            return Err(JidxError::Invalid("document posting range"));
        }
        let bytes = reader.checked_bytes(cursor, inline_end)?;
        let document_id = read_u32(bytes, 0);
        if document_id >= reader.header().document_count
            || previous.is_some_and(|previous| previous >= document_id)
        {
            return Err(JidxError::Invalid("document postings"));
        }
        let document_record = reader.document_record(document_id)?;
        let local_contig = read_u32(bytes, 4);
        let payload = read_u64(bytes, 8);
        let (occurrence_count, location) = if local_contig == EXTERNAL_SENTINEL {
            let row_end = cursor
                .checked_add(EXTERNAL_DOCUMENT_POSTING_SIZE)
                .ok_or(JidxError::Invalid("document posting range"))?;
            if row_end > section_end {
                return Err(JidxError::Invalid("document posting range"));
            }
            let bytes = reader.checked_bytes(cursor, row_end)?;
            let cold = reader.header().section(SectionKind::ContigPostings);
            let data_start = cold
                .offset
                .checked_add(payload)
                .ok_or(JidxError::Invalid("contig posting offset"))?;
            let occurrence_count = read_u64(bytes, 16);
            let encoded_length = read_u64(bytes, 24);
            let data_end = data_start
                .checked_add(encoded_length)
                .ok_or(JidxError::Invalid("contig posting range"))?;
            let cold_end = cold
                .offset
                .checked_add(cold.length)
                .ok_or(JidxError::Invalid("contig posting range"))?;
            let minimum_length = occurrence_count
                .checked_mul(2)
                .ok_or(JidxError::Invalid("contig posting block"))?;
            if occurrence_count == 0
                || occurrence_count > reader.header().occurrence_count
                || encoded_length < minimum_length
                || data_end > cold_end
            {
                return Err(JidxError::Invalid("contig posting block"));
            }
            (
                occurrence_count,
                DocumentLocation::Delta {
                    data_start,
                    data_end,
                    occurrence_count,
                },
            )
        } else {
            if local_contig >= document_record.contig_count {
                return Err(JidxError::Invalid("contig posting ordinal"));
            }
            (
                1,
                DocumentLocation::Inline {
                    local_contig,
                    payload,
                },
            )
        };
        if occurrence_count > reader.header().occurrence_count {
            return Err(JidxError::Invalid("occurrence count"));
        }
        seed_occurrence_count = seed_occurrence_count
            .checked_add(occurrence_count)
            .ok_or(JidxError::Invalid("occurrence count"))?;
        if seed_occurrence_count > reader.header().occurrence_count {
            return Err(JidxError::Invalid("occurrence count"));
        }
        documents.push(SeedDocument {
            metagenome_id: document_id,
            occurrence_count,
            seed_key: seed.packed_key,
            document_id,
            location,
        });
        previous = Some(document_id);
        cursor = if local_contig == EXTERNAL_SENTINEL {
            cursor
                .checked_add(EXTERNAL_DOCUMENT_POSTING_SIZE)
                .ok_or(JidxError::Invalid("document posting range"))?
        } else {
            inline_end
        };
    }
    let relative_end = cursor
        .checked_sub(section.offset)
        .ok_or(JidxError::Invalid("document posting range"))?;
    Ok((documents, relative_end))
}

pub(crate) fn document_occurrences(
    reader: &JidxReader,
    seed: SeedEntry,
    document: SeedDocument,
) -> Result<Vec<SeedOccurrence>, JidxError> {
    if document.seed_key != seed.packed_key
        || document.metagenome_id != document.document_id
        || document.occurrence_count != location_count(document.location)
    {
        return Err(JidxError::Invalid("seed document"));
    }
    let record = reader.document_record(document.document_id)?;
    match document.location {
        DocumentLocation::Inline {
            local_contig,
            payload,
        } => Ok(vec![occurrence(
            record.contig_start,
            record.contig_count,
            local_contig,
            payload >> 1,
            payload & 1 == 1,
        )?]),
        DocumentLocation::Delta {
            data_start,
            data_end,
            occurrence_count,
            ..
        } => decode_delta_occurrences(
            reader.checked_bytes(data_start, data_end)?,
            record.contig_start,
            record.contig_count,
            occurrence_count,
        ),
    }
}

fn location_count(location: DocumentLocation) -> u64 {
    match location {
        DocumentLocation::Inline { .. } => 1,
        DocumentLocation::Delta {
            occurrence_count, ..
        } => occurrence_count,
    }
}

fn decode_delta_occurrences(
    mut bytes: &[u8],
    contig_start: u32,
    contig_count: u32,
    count: u64,
) -> Result<Vec<SeedOccurrence>, JidxError> {
    let capacity = usize::try_from(count).map_err(|_| JidxError::Invalid("occurrence count"))?;
    let mut occurrences = Vec::new();
    occurrences
        .try_reserve_exact(capacity)
        .map_err(|_| JidxError::Invalid("occurrence count"))?;
    let mut previous_local = 0u32;
    let mut previous_position = 0u64;
    for index in 0..count {
        let packed_contig = take_varint(&mut bytes)?;
        let delta_contig = packed_contig >> 1;
        let orientation = packed_contig & 1 == 1;
        let encoded_position = take_varint(&mut bytes)?;
        let local_contig = if index == 0 {
            u32::try_from(delta_contig).map_err(|_| JidxError::Invalid("contig posting ordinal"))?
        } else {
            let delta_contig = u32::try_from(delta_contig)
                .map_err(|_| JidxError::Invalid("contig posting ordinal"))?;
            previous_local
                .checked_add(delta_contig)
                .ok_or(JidxError::Invalid("contig posting ordinal"))?
        };
        let position = if index != 0 && local_contig == previous_local {
            if encoded_position == 0 {
                return Err(JidxError::Invalid("contig posting position"));
            }
            previous_position
                .checked_add(encoded_position)
                .ok_or(JidxError::Invalid("contig posting position"))?
        } else {
            encoded_position
        };
        occurrences.push(occurrence(
            contig_start,
            contig_count,
            local_contig,
            position,
            orientation,
        )?);
        previous_local = local_contig;
        previous_position = position;
    }
    if !bytes.is_empty() {
        return Err(JidxError::Invalid("contig posting length"));
    }
    Ok(occurrences)
}

fn occurrence(
    contig_start: u32,
    contig_count: u32,
    local_contig: u32,
    position: u64,
    canonical_orientation: bool,
) -> Result<SeedOccurrence, JidxError> {
    if local_contig >= contig_count {
        return Err(JidxError::Invalid("contig posting ordinal"));
    }
    Ok(SeedOccurrence {
        contig_id: contig_start
            .checked_add(local_contig)
            .ok_or(JidxError::Invalid("contig posting ordinal"))?,
        position,
        canonical_orientation,
    })
}

fn take_varint(bytes: &mut &[u8]) -> Result<u64, JidxError> {
    let input = *bytes;
    let (value, remaining) = unsigned_varint::decode::u64(input)
        .map_err(|_| JidxError::Invalid("contig posting varint"))?;
    let consumed = input.len() - remaining.len();
    let mut buffer = unsigned_varint::encode::u64_buffer();
    if unsigned_varint::encode::u64(value, &mut buffer) != &input[..consumed] {
        return Err(JidxError::Invalid("contig posting varint"));
    }
    *bytes = remaining;
    Ok(value)
}

#[derive(Clone, Copy)]
struct SeedRecord {
    packed_key: u64,
    document_offset: u64,
    document_count: u32,
}

impl From<SeedRecord> for SeedEntry {
    fn from(record: SeedRecord) -> Self {
        Self {
            packed_key: record.packed_key,
            document_frequency: record.document_count,
            document_offset: record.document_offset,
        }
    }
}

fn seed_record(reader: &JidxReader, index: u64) -> Result<SeedRecord, JidxError> {
    let section = reader.header().section(SectionKind::Seeds);
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
    })
}

fn validate_record(header: &Header, record: SeedRecord) -> Result<(), JidxError> {
    seed_length(header.k, header.rescue_k15, record.packed_key)?;
    if record.document_count == 0
        || record.document_count > header.document_count
        || !record
            .document_offset
            .is_multiple_of(u64::from(DOCUMENT_POSTING_SIZE))
    {
        return Err(JidxError::Invalid("seed record"));
    }
    if record.document_offset > header.section(SectionKind::DocumentPostings).length {
        return Err(JidxError::Invalid("seed posting range"));
    }
    Ok(())
}

pub(crate) fn encode_seed(seed: SeedEntry) -> [u8; SEED_RECORD_SIZE as usize] {
    let mut bytes = [0; SEED_RECORD_SIZE as usize];
    crate::jidx::put_u64(&mut bytes, 0, seed.packed_key);
    crate::jidx::put_u64(&mut bytes, 8, seed.document_offset);
    crate::jidx::put_u32(&mut bytes, 16, seed.document_frequency);
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn varints(values: &[u64]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for value in values {
            let mut buffer = unsigned_varint::encode::u64_buffer();
            bytes.extend_from_slice(unsigned_varint::encode::u64(*value, &mut buffer));
        }
        bytes
    }

    #[test]
    fn delta_decoder_rejects_truncation_duplicate_and_trailing_bytes() {
        assert!(decode_delta_occurrences(&[0], 0, 1, 1).is_err());
        assert!(decode_delta_occurrences(&[0, 5, 0, 0], 0, 1, 2).is_err());
        assert!(decode_delta_occurrences(&[0, 5, 0], 0, 1, 1).is_err());
    }

    #[test]
    fn delta_decoder_checks_contig_and_position_overflow() {
        let contig_overflow = varints(&[((u64::from(u32::MAX) + 1) << 1), 0]);
        assert!(decode_delta_occurrences(&contig_overflow, 0, u32::MAX, 1).is_err());

        let position_overflow = varints(&[0, u64::MAX, 0, 1]);
        assert!(decode_delta_occurrences(&position_overflow, 0, 1, 2).is_err());
    }
}
