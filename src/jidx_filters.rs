use crate::jidx::{JidxError, PAGE_SIZE, SEED_RECORD_SIZE, SectionKind};
use crate::jidx_postings;
use crate::jidx_reader::JidxReader;
use std::fs::OpenOptions;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use xorf::{BinaryFuse8, BinaryFuse8Ref, DmaSerializable, Filter, FilterRef};

const GROUP_KEYS: u64 = 1_000_000;
const DIRECTORY_RECORD_SIZE: u64 = 40;
const DESCRIPTOR_SIZE: usize = 20;

pub(crate) struct FilterDirectory {
    records: Box<[FilterRecord]>,
    seed_count: u64,
    section_offset: u64,
    section_length: u64,
}

#[derive(Clone, Copy)]
struct FilterRecord {
    first_key: u64,
    last_key: u64,
    filter_start: u64,
    filter_end: u64,
    key_count: u32,
}

pub(crate) fn load(reader: &JidxReader) -> Result<FilterDirectory, JidxError> {
    let header = reader.header();
    let section = header.section(SectionKind::SeedFilters);
    let group_count = group_count(header.seed_count)?;
    if group_count == 0 {
        if section.length != 0 {
            return Err(JidxError::Invalid("seed filter section"));
        }
        return Ok(FilterDirectory {
            records: Box::new([]),
            seed_count: 0,
            section_offset: section.offset,
            section_length: 0,
        });
    }

    let directory_length = group_count
        .checked_mul(DIRECTORY_RECORD_SIZE)
        .ok_or(JidxError::Invalid("seed filter directory"))?;
    let directory_end = section
        .offset
        .checked_add(directory_length)
        .ok_or(JidxError::Invalid("seed filter directory"))?;
    let section_end = section
        .offset
        .checked_add(section.length)
        .ok_or(JidxError::Invalid("seed filter section"))?;
    if directory_end > section_end {
        return Err(JidxError::Invalid("seed filter directory"));
    }
    let capacity =
        usize::try_from(group_count).map_err(|_| JidxError::Invalid("seed filter directory"))?;
    let mut records = Vec::new();
    records
        .try_reserve_exact(capacity)
        .map_err(|_| JidxError::Invalid("seed filter directory"))?;
    let mut expected_filter_offset = align_page(directory_length)?;
    let mut previous_last = None;
    for group in 0..group_count {
        let relative = group
            .checked_mul(DIRECTORY_RECORD_SIZE)
            .ok_or(JidxError::Invalid("seed filter directory"))?;
        let start = section
            .offset
            .checked_add(relative)
            .ok_or(JidxError::Invalid("seed filter directory"))?;
        let end = start
            .checked_add(DIRECTORY_RECORD_SIZE)
            .ok_or(JidxError::Invalid("seed filter directory"))?;
        let bytes = reader.checked_bytes(start, end)?;
        let first_key = read_u64(bytes, 0);
        let last_key = read_u64(bytes, 8);
        let filter_offset = read_u64(bytes, 16);
        let filter_length = read_u64(bytes, 24);
        let key_count = read_u32(bytes, 32);
        let descriptor_length = read_u16(bytes, 36);
        if read_u16(bytes, 38) != 0
            || descriptor_length != DESCRIPTOR_SIZE as u16
            || first_key > last_key
            || previous_last.is_some_and(|previous| previous >= first_key)
            || filter_offset != expected_filter_offset
            || !filter_offset.is_multiple_of(PAGE_SIZE)
            || filter_length < DESCRIPTOR_SIZE as u64 + 1
        {
            return Err(JidxError::Invalid("seed filter record"));
        }
        let expected_count = (header.seed_count - group * GROUP_KEYS).min(GROUP_KEYS);
        if u64::from(key_count) != expected_count {
            return Err(JidxError::Invalid("seed filter key count"));
        }
        let first_ordinal = group
            .checked_mul(GROUP_KEYS)
            .ok_or(JidxError::Invalid("seed filter ordinal"))?;
        let last_ordinal = first_ordinal
            .checked_add(expected_count - 1)
            .ok_or(JidxError::Invalid("seed filter ordinal"))?;
        if jidx_postings::entry(reader, first_ordinal)?.packed_key != first_key
            || jidx_postings::entry(reader, last_ordinal)?.packed_key != last_key
        {
            return Err(JidxError::Invalid("seed filter key range"));
        }
        let filter_start = section
            .offset
            .checked_add(filter_offset)
            .ok_or(JidxError::Invalid("seed filter range"))?;
        let filter_end = filter_start
            .checked_add(filter_length)
            .ok_or(JidxError::Invalid("seed filter range"))?;
        if filter_end > section_end {
            return Err(JidxError::Invalid("seed filter range"));
        }
        records.push(FilterRecord {
            first_key,
            last_key,
            filter_start,
            filter_end,
            key_count,
        });
        expected_filter_offset = align_page(
            filter_offset
                .checked_add(filter_length)
                .ok_or(JidxError::Invalid("seed filter range"))?,
        )?;
        previous_last = Some(last_key);
    }
    if records
        .last()
        .is_none_or(|record| record.filter_end != section_end)
    {
        return Err(JidxError::Invalid("seed filter coverage"));
    }
    Ok(FilterDirectory {
        records: records.into_boxed_slice(),
        seed_count: header.seed_count,
        section_offset: section.offset,
        section_length: section.length,
    })
}

pub(crate) fn contains(
    reader: &JidxReader,
    directory: &FilterDirectory,
    key: u64,
) -> Result<bool, JidxError> {
    validate_identity(reader, directory)?;
    let index = directory
        .records
        .partition_point(|record| record.last_key < key);
    let Some(record) = directory.records.get(index) else {
        return Ok(false);
    };
    if key < record.first_key {
        return Ok(false);
    }
    page_local_contains(reader, *record, key)
}

#[cfg(test)]
pub(crate) fn reference_contains(
    reader: &JidxReader,
    directory: &FilterDirectory,
    key: u64,
) -> Result<bool, JidxError> {
    validate_identity(reader, directory)?;
    let index = directory
        .records
        .partition_point(|record| record.last_key < key);
    let Some(record) = directory.records.get(index) else {
        return Ok(false);
    };
    if key < record.first_key {
        return Ok(false);
    }
    Ok(filter(reader, *record)?.contains(&key))
}

#[cfg(test)]
pub(crate) fn fingerprint_offsets(
    reader: &JidxReader,
    directory: &FilterDirectory,
    key: u64,
) -> Result<[u64; 3], JidxError> {
    validate_identity(reader, directory)?;
    let index = directory
        .records
        .partition_point(|record| record.last_key < key);
    let record = *directory
        .records
        .get(index)
        .filter(|record| key >= record.first_key)
        .ok_or(JidxError::Invalid("seed filter test key"))?;
    let descriptor_end = record
        .filter_start
        .checked_add(DESCRIPTOR_SIZE as u64)
        .ok_or(JidxError::Invalid("seed filter descriptor"))?;
    let descriptor_bytes = reader.checked_bytes(record.filter_start, descriptor_end)?;
    let fingerprint_length = usize::try_from(record.filter_end - descriptor_end)
        .map_err(|_| JidxError::Invalid("seed filter fingerprints"))?;
    validate_descriptor(descriptor_bytes, fingerprint_length)?;
    let (_, indexes) = fingerprint_indices(descriptor_bytes, key);
    Ok(indexes.map(|index| descriptor_end + u64::from(index)))
}

fn page_local_contains(
    reader: &JidxReader,
    record: FilterRecord,
    key: u64,
) -> Result<bool, JidxError> {
    let descriptor_end = record
        .filter_start
        .checked_add(DESCRIPTOR_SIZE as u64)
        .ok_or(JidxError::Invalid("seed filter descriptor"))?;
    let descriptor_bytes = reader.checked_bytes(record.filter_start, descriptor_end)?;
    let fingerprint_length = record
        .filter_end
        .checked_sub(descriptor_end)
        .and_then(|length| usize::try_from(length).ok())
        .ok_or(JidxError::Invalid("seed filter fingerprints"))?;
    validate_descriptor(descriptor_bytes, fingerprint_length)?;
    let (fingerprint, indexes) = fingerprint_indices(descriptor_bytes, key);
    let mut value = fingerprint;
    for index in indexes {
        let start = descriptor_end
            .checked_add(u64::from(index))
            .ok_or(JidxError::Invalid("seed filter fingerprint"))?;
        let end = start
            .checked_add(1)
            .ok_or(JidxError::Invalid("seed filter fingerprint"))?;
        let byte = reader
            .checked_bytes(start, end)?
            .first()
            .copied()
            .ok_or(JidxError::Invalid("seed filter fingerprint"))?;
        value ^= byte;
    }
    Ok(value == 0)
}

// Pinned xorf 0.12 BinaryFuse8 DMA addressing; descriptors are validated first.
// The avalanche is Austin Appleby's public-domain MurmurHash3 fmix64.
fn fingerprint_indices(descriptor: &[u8], key: u64) -> (u8, [u32; 3]) {
    let mut hash = key.wrapping_add(read_u64(descriptor, 0));
    hash ^= hash >> 33;
    hash = hash.wrapping_mul(0xff51_afd7_ed55_8ccd);
    hash ^= hash >> 33;
    hash = hash.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    hash ^= hash >> 33;
    let segment_length = read_u32(descriptor, 8);
    let segment_mask = read_u32(descriptor, 12);
    let segment_count_length = read_u32(descriptor, 16);
    let first = ((u128::from(hash) * u128::from(segment_count_length)) >> 64) as u32;
    let second = (first + segment_length) ^ ((hash >> 18) as u32 & segment_mask);
    let third = (first + 2 * segment_length) ^ (hash as u32 & segment_mask);
    (xorf::fingerprint!(hash) as u8, [first, second, third])
}

pub(crate) fn audit(reader: &JidxReader, directory: &FilterDirectory) -> Result<(), JidxError> {
    validate_identity(reader, directory)?;
    let mut previous_end = directory
        .records
        .len()
        .checked_mul(DIRECTORY_RECORD_SIZE as usize)
        .and_then(|length| u64::try_from(length).ok())
        .ok_or(JidxError::Invalid("seed filter directory"))?;
    let mut ordinal = 0u64;
    for record in &directory.records {
        let filter_offset = record
            .filter_start
            .checked_sub(directory.section_offset)
            .ok_or(JidxError::Invalid("seed filter range"))?;
        validate_zeroes(
            reader,
            directory.section_offset,
            previous_end,
            filter_offset,
        )?;
        let filter = filter(reader, *record)?;
        for _ in 0..record.key_count {
            let key = jidx_postings::entry(reader, ordinal)?.packed_key;
            if !filter.contains(&key) {
                return Err(JidxError::Invalid("seed filter false negative"));
            }
            ordinal = ordinal
                .checked_add(1)
                .ok_or(JidxError::Invalid("seed filter ordinal"))?;
        }
        previous_end = record
            .filter_end
            .checked_sub(directory.section_offset)
            .ok_or(JidxError::Invalid("seed filter range"))?;
    }
    if ordinal != directory.seed_count || previous_end != directory.section_length {
        return Err(JidxError::Invalid("seed filter coverage"));
    }
    Ok(())
}

pub(crate) fn build(seeds: impl Read, output: &Path, seed_count: u64) -> io::Result<()> {
    if BinaryFuse8::DESCRIPTOR_LEN != DESCRIPTOR_SIZE {
        return Err(invalid_data("unexpected BinaryFuse8 descriptor size"));
    }
    let group_count = group_count(seed_count).map_err(format_error)?;
    let directory_length = group_count
        .checked_mul(DIRECTORY_RECORD_SIZE)
        .ok_or_else(|| invalid_data("seed filter directory overflow"))?;
    let mut output = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)?,
    );
    write_zeroes(&mut output, directory_length)?;
    write_zeroes(
        &mut output,
        align_page(directory_length).map_err(format_error)? - directory_length,
    )?;

    let capacity =
        usize::try_from(GROUP_KEYS).map_err(|_| invalid_data("seed filter group size"))?;
    let mut keys = Vec::new();
    keys.try_reserve_exact(capacity)
        .map_err(|_| invalid_data("seed filter key allocation"))?;
    let mut seeds = BufReader::new(seeds);
    let mut records = Vec::new();
    records
        .try_reserve_exact(
            usize::try_from(group_count)
                .map_err(|_| invalid_data("seed filter directory allocation"))?,
        )
        .map_err(|_| invalid_data("seed filter directory allocation"))?;
    let mut previous_key = None;
    let mut remaining = seed_count;
    while remaining != 0 {
        keys.clear();
        let count = remaining.min(GROUP_KEYS);
        for _ in 0..count {
            let mut row = [0; SEED_RECORD_SIZE as usize];
            seeds.read_exact(&mut row)?;
            let key = read_u64(&row, 0);
            if row[20..].iter().any(|byte| *byte != 0)
                || previous_key.is_some_and(|previous| previous >= key)
            {
                return Err(invalid_data("invalid seed filter input"));
            }
            keys.push(key);
            previous_key = Some(key);
        }
        let filter = BinaryFuse8::try_from(keys.as_slice())
            .map_err(|_| invalid_data("BinaryFuse8 construction failed"))?;
        if keys.iter().any(|key| !filter.contains(key)) {
            return Err(invalid_data("BinaryFuse8 construction false negative"));
        }
        let position = output.stream_position()?;
        let filter_offset = align_page(position).map_err(format_error)?;
        write_zeroes(&mut output, filter_offset - position)?;
        let mut descriptor = [0; DESCRIPTOR_SIZE];
        filter.dma_copy_descriptor_to(&mut descriptor);
        output.write_all(&descriptor)?;
        output.write_all(filter.dma_fingerprints())?;
        let filter_length = u64::try_from(DESCRIPTOR_SIZE)
            .ok()
            .and_then(|length| {
                u64::try_from(filter.dma_fingerprints().len())
                    .ok()
                    .and_then(|fingerprints| length.checked_add(fingerprints))
            })
            .ok_or_else(|| invalid_data("seed filter length overflow"))?;
        records.push(EncodedRecord {
            first_key: keys[0],
            last_key: keys[keys.len() - 1],
            filter_offset,
            filter_length,
            key_count: u32::try_from(keys.len())
                .map_err(|_| invalid_data("seed filter key count"))?,
        });
        remaining -= count;
    }
    let mut trailing = [0; 1];
    if seeds.read(&mut trailing)? != 0 {
        return Err(invalid_data("trailing seed filter input"));
    }
    output.flush()?;
    output.seek(SeekFrom::Start(0))?;
    for record in records {
        output.write_all(&record.encode())?;
    }
    output.flush()
}

#[derive(Clone, Copy)]
struct EncodedRecord {
    first_key: u64,
    last_key: u64,
    filter_offset: u64,
    filter_length: u64,
    key_count: u32,
}

impl EncodedRecord {
    fn encode(self) -> [u8; DIRECTORY_RECORD_SIZE as usize] {
        let mut bytes = [0; DIRECTORY_RECORD_SIZE as usize];
        bytes[0..8].copy_from_slice(&self.first_key.to_le_bytes());
        bytes[8..16].copy_from_slice(&self.last_key.to_le_bytes());
        bytes[16..24].copy_from_slice(&self.filter_offset.to_le_bytes());
        bytes[24..32].copy_from_slice(&self.filter_length.to_le_bytes());
        bytes[32..36].copy_from_slice(&self.key_count.to_le_bytes());
        bytes[36..38].copy_from_slice(&(DESCRIPTOR_SIZE as u16).to_le_bytes());
        bytes
    }
}

fn filter<'a>(
    reader: &'a JidxReader,
    record: FilterRecord,
) -> Result<BinaryFuse8Ref<'a>, JidxError> {
    let bytes = reader.checked_bytes(record.filter_start, record.filter_end)?;
    let (descriptor, fingerprints) = bytes.split_at(DESCRIPTOR_SIZE);
    validate_descriptor(descriptor, fingerprints.len())?;
    Ok(BinaryFuse8Ref::from_dma(descriptor, fingerprints))
}

fn validate_descriptor(descriptor: &[u8], fingerprint_length: usize) -> Result<(), JidxError> {
    if descriptor.len() != DESCRIPTOR_SIZE {
        return Err(JidxError::Invalid("seed filter descriptor"));
    }
    let segment_length = read_u32(descriptor, 8);
    let segment_mask = read_u32(descriptor, 12);
    let segment_count_length = read_u32(descriptor, 16);
    if !(4..=262_144).contains(&segment_length)
        || !segment_length.is_power_of_two()
        || segment_mask != segment_length - 1
        || segment_count_length == 0
        || !segment_count_length.is_multiple_of(segment_length)
    {
        return Err(JidxError::Invalid("seed filter descriptor"));
    }
    let expected = segment_length
        .checked_mul(2)
        .and_then(|tail| segment_count_length.checked_add(tail))
        .and_then(|length| usize::try_from(length).ok())
        .ok_or(JidxError::Invalid("seed filter descriptor"))?;
    if fingerprint_length != expected {
        return Err(JidxError::Invalid("seed filter fingerprints"));
    }
    Ok(())
}

fn validate_identity(reader: &JidxReader, directory: &FilterDirectory) -> Result<(), JidxError> {
    let section = reader.header().section(SectionKind::SeedFilters);
    if directory.seed_count != reader.header().seed_count
        || directory.section_offset != section.offset
        || directory.section_length != section.length
    {
        return Err(JidxError::Invalid("seed filter identity"));
    }
    Ok(())
}

fn validate_zeroes(
    reader: &JidxReader,
    section_offset: u64,
    start: u64,
    end: u64,
) -> Result<(), JidxError> {
    if start > end {
        return Err(JidxError::Invalid("seed filter padding"));
    }
    let start = section_offset
        .checked_add(start)
        .ok_or(JidxError::Invalid("seed filter padding"))?;
    let end = section_offset
        .checked_add(end)
        .ok_or(JidxError::Invalid("seed filter padding"))?;
    if reader
        .checked_bytes(start, end)?
        .iter()
        .any(|byte| *byte != 0)
    {
        return Err(JidxError::Invalid("seed filter padding"));
    }
    Ok(())
}

fn group_count(seed_count: u64) -> Result<u64, JidxError> {
    seed_count
        .checked_add(GROUP_KEYS - 1)
        .map(|count| count / GROUP_KEYS)
        .ok_or(JidxError::Invalid("seed filter group count"))
}

fn align_page(value: u64) -> Result<u64, JidxError> {
    value
        .checked_add(PAGE_SIZE - 1)
        .map(|value| value & !(PAGE_SIZE - 1))
        .ok_or(JidxError::Invalid("seed filter alignment"))
}

fn read_u16(bytes: &[u8], offset: usize) -> u16 {
    u16::from_le_bytes(bytes[offset..offset + 2].try_into().expect("filter u16"))
}

fn read_u32(bytes: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(bytes[offset..offset + 4].try_into().expect("filter u32"))
}

fn read_u64(bytes: &[u8], offset: usize) -> u64 {
    u64::from_le_bytes(bytes[offset..offset + 8].try_into().expect("filter u64"))
}

fn write_zeroes(writer: &mut impl Write, mut length: u64) -> io::Result<()> {
    let zeroes = [0; 4096];
    while length != 0 {
        let count = usize::try_from(length.min(zeroes.len() as u64))
            .map_err(|_| invalid_data("seed filter padding"))?;
        writer.write_all(&zeroes[..count])?;
        length -= count as u64;
    }
    Ok(())
}

fn format_error(error: JidxError) -> io::Error {
    invalid_data(error.to_string())
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    struct SeedRows {
        next: u64,
        end: u64,
        row: [u8; SEED_RECORD_SIZE as usize],
        offset: usize,
    }

    impl SeedRows {
        fn new(end: u64) -> Self {
            Self {
                next: 0,
                end,
                row: [0; SEED_RECORD_SIZE as usize],
                offset: SEED_RECORD_SIZE as usize,
            }
        }
    }

    impl Read for SeedRows {
        fn read(&mut self, mut output: &mut [u8]) -> io::Result<usize> {
            let requested = output.len();
            while !output.is_empty() {
                if self.offset == self.row.len() {
                    if self.next == self.end {
                        break;
                    }
                    self.row = [0; SEED_RECORD_SIZE as usize];
                    self.row[..8].copy_from_slice(&self.next.to_le_bytes());
                    self.next += 1;
                    self.offset = 0;
                }
                let count = output.len().min(self.row.len() - self.offset);
                output[..count].copy_from_slice(&self.row[self.offset..self.offset + count]);
                self.offset += count;
                output = &mut output[count..];
            }
            Ok(requested - output.len())
        }
    }

    #[test]
    fn directory_record_has_explicit_little_endian_bytes() {
        let bytes = EncodedRecord {
            first_key: 1,
            last_key: 2,
            filter_offset: 4096,
            filter_length: 32,
            key_count: 3,
        }
        .encode();
        assert_eq!(&bytes[0..8], &1u64.to_le_bytes());
        assert_eq!(&bytes[8..16], &2u64.to_le_bytes());
        assert_eq!(&bytes[16..24], &4096u64.to_le_bytes());
        assert_eq!(&bytes[24..32], &32u64.to_le_bytes());
        assert_eq!(&bytes[32..36], &3u32.to_le_bytes());
        assert_eq!(&bytes[36..38], &20u16.to_le_bytes());
        assert_eq!(&bytes[38..40], &[0, 0]);
    }

    #[test]
    fn empty_seed_stream_builds_empty_section() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("filters");
        build(Cursor::new([]), &path, 0).unwrap();
        assert_eq!(std::fs::metadata(path).unwrap().len(), 0);
    }

    #[test]
    fn small_filter_build_is_deterministic_and_self_contained() {
        let directory = tempfile::tempdir().unwrap();
        let mut seeds = Vec::new();
        for key in [3u64, 7, 11] {
            let mut row = [0; SEED_RECORD_SIZE as usize];
            row[..8].copy_from_slice(&key.to_le_bytes());
            seeds.extend_from_slice(&row);
        }
        let first = directory.path().join("first");
        let second = directory.path().join("second");
        build(Cursor::new(&seeds), &first, 3).unwrap();
        build(Cursor::new(&seeds), &second, 3).unwrap();
        let bytes = std::fs::read(first).unwrap();
        assert_eq!(bytes, std::fs::read(second).unwrap());
        assert_eq!(&bytes[0..8], &3u64.to_le_bytes());
        assert_eq!(&bytes[8..16], &11u64.to_le_bytes());
        assert_eq!(&bytes[16..24], &PAGE_SIZE.to_le_bytes());
        assert_eq!(&bytes[32..36], &3u32.to_le_bytes());
        assert_eq!(&bytes[36..38], &20u16.to_le_bytes());
        let length = u64::from_le_bytes(bytes[24..32].try_into().unwrap()) as usize;
        let block = &bytes[PAGE_SIZE as usize..PAGE_SIZE as usize + length];
        let (descriptor, fingerprints) = block.split_at(DESCRIPTOR_SIZE);
        validate_descriptor(descriptor, fingerprints.len()).unwrap();
        let filter = BinaryFuse8Ref::from_dma(descriptor, fingerprints);
        assert!([3, 7, 11].iter().all(|key| filter.contains(key)));
    }

    #[test]
    fn million_key_boundary_builds_full_and_singleton_groups() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("filters");
        build(SeedRows::new(GROUP_KEYS + 1), &path, GROUP_KEYS + 1).unwrap();
        let bytes = std::fs::read(path).unwrap();

        let first_offset = read_u64(&bytes, 16);
        let first_length = read_u64(&bytes, 24);
        let second_offset = read_u64(&bytes, DIRECTORY_RECORD_SIZE as usize + 16);
        let second_length = read_u64(&bytes, DIRECTORY_RECORD_SIZE as usize + 24);
        assert_eq!(read_u64(&bytes, 0), 0);
        assert_eq!(read_u64(&bytes, 8), GROUP_KEYS - 1);
        assert_eq!(read_u64(&bytes, DIRECTORY_RECORD_SIZE as usize), GROUP_KEYS);
        assert_eq!(
            read_u64(&bytes, DIRECTORY_RECORD_SIZE as usize + 8),
            GROUP_KEYS
        );
        assert_eq!(read_u32(&bytes, 32), GROUP_KEYS as u32);
        assert_eq!(read_u32(&bytes, DIRECTORY_RECORD_SIZE as usize + 32), 1);
        assert_eq!(first_offset, PAGE_SIZE);
        assert_eq!(
            second_offset,
            align_page(first_offset + first_length).unwrap()
        );
        assert_eq!(bytes.len() as u64, second_offset + second_length);

        for (start, count, offset, length) in [
            (0, GROUP_KEYS, first_offset, first_length),
            (GROUP_KEYS, 1, second_offset, second_length),
        ] {
            let block = &bytes[offset as usize..(offset + length) as usize];
            let (descriptor, fingerprints) = block.split_at(DESCRIPTOR_SIZE);
            validate_descriptor(descriptor, fingerprints.len()).unwrap();
            let filter = BinaryFuse8Ref::from_dma(descriptor, fingerprints);
            assert!((start..start + count).all(|key| filter.contains(&key)));
        }
    }
}
