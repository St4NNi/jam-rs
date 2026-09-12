use crate::jidx::{put_u32, put_u64};
use crate::jidx_filters::{build_binary_fuse, validate_descriptor};
use crate::shared_file::SharedFile;
use crate::shared_format::{Section, SharedError, read_u32, read_u64};
use crate::shared_reader::CoreRow;
use xorf::{BinaryFuse8Ref, DmaSerializable, Filter, FilterRef};

const METADATA_BYTES: usize = 128;
const MAX_BYTES: usize = 8 * 1024 * 1024;
const PREFIX_END: u32 = 65_536;

#[cfg(test)]
thread_local! {
    pub(crate) static FILTER_TEST_LIMIT: std::cell::Cell<usize> = const { std::cell::Cell::new(MAX_BYTES) };
}

pub(crate) fn core_input(core: u32) -> u64 {
    u64::from(core)
}

pub(crate) fn build(source: &SharedFile, maximum_keys: usize) -> Result<Vec<u8>, SharedError> {
    source.verify_unchanged()?;
    if source.header.version != 3 {
        return Err(SharedError::Invalid("filter requires split cores"));
    }
    let count =
        usize::try_from(source.header.core_count).map_err(|_| SharedError::ResourceLimit)?;
    let mut keys = Vec::new();
    keys.try_reserve_exact(count.min(maximum_keys))
        .map_err(|_| SharedError::ResourceLimit)?;
    let mut end_prefix = if maximum_keys == 0 { 0 } else { PREFIX_END };
    let mut previous = None;
    if maximum_keys != 0 {
        'pages: for first in (0..count).step_by(1024) {
            let page = source.section(
                Section::Cores,
                first as u64 * 4,
                (count - first).min(1024) as u64 * 4,
            )?;
            for bytes in page.as_chunks::<4>().0 {
                let core = CoreRow::decode_key(bytes)?;
                if previous.is_some_and(|before| before >= core) {
                    return Err(SharedError::Invalid("filter core order"));
                }
                previous = Some(core);
                if keys.len() == maximum_keys {
                    end_prefix = core >> 14;
                    keys.truncate(keys.partition_point(|key| *key < u64::from(end_prefix) << 14));
                    break 'pages;
                }
                keys.push(core_input(core));
            }
        }
    }
    source.verify_unchanged()?;
    encode(
        &keys,
        end_prefix,
        source.header.core_count,
        &source.header.manifest_sha256,
        &source.header.body_sha256,
    )
}

fn encode(
    keys: &[u64],
    end_prefix: u32,
    dictionary_count: u64,
    manifest: &[u8; 32],
    source_body: &[u8; 32],
) -> Result<Vec<u8>, SharedError> {
    let filter = if keys.is_empty() {
        None
    } else {
        Some(build_binary_fuse(keys)?)
    };
    let fingerprint_len = filter.as_ref().map_or(0, |f| f.dma_fingerprints().len());
    let descriptor_len = if filter.is_some() { 20 } else { 0 };
    let length = METADATA_BYTES + descriptor_len + fingerprint_len;
    if length > MAX_BYTES {
        return Err(SharedError::ResourceLimit);
    }
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(length)
        .map_err(|_| SharedError::ResourceLimit)?;
    bytes.resize(METADATA_BYTES + descriptor_len, 0);
    bytes[..8].copy_from_slice(b"JCFUSE1\0");
    put_u32(&mut bytes, 8, 1);
    put_u32(&mut bytes, 12, 1);
    put_u32(&mut bytes, 16, end_prefix);
    put_u32(&mut bytes, 20, descriptor_len as u32);
    put_u64(&mut bytes, 24, keys.len() as u64);
    put_u64(&mut bytes, 32, dictionary_count);
    put_u64(&mut bytes, 40, fingerprint_len as u64);
    bytes[48..80].copy_from_slice(manifest);
    bytes[80..112].copy_from_slice(source_body);
    if let Some(filter) = filter {
        filter.dma_copy_descriptor_to(&mut bytes[METADATA_BYTES..]);
        bytes.extend_from_slice(filter.dma_fingerprints());
    }
    validate(&bytes, dictionary_count, manifest, source_body)?;
    if !keys.is_empty() {
        let borrowed = BinaryFuse8Ref::from_dma(&bytes[128..148], &bytes[148..]);
        if keys.iter().any(|key| !borrowed.contains(key)) {
            return Err(SharedError::Invalid("filter borrowed false negative"));
        }
    }
    Ok(bytes)
}

fn validate(
    bytes: &[u8],
    dictionary_count: u64,
    manifest: &[u8; 32],
    source_body: &[u8; 32],
) -> Result<(u32, u64), SharedError> {
    if bytes.len() < METADATA_BYTES
        || &bytes[..8] != b"JCFUSE1\0"
        || read_u32(bytes, 8) != 1
        || read_u32(bytes, 12) != 1
        || read_u64(bytes, 32) != dictionary_count
        || &bytes[48..80] != manifest
        || &bytes[80..112] != source_body
        || bytes[112..128].iter().any(|&byte| byte != 0)
    {
        return Err(SharedError::Invalid("core filter metadata"));
    }
    let end_prefix = read_u32(bytes, 16);
    let covered_count = read_u64(bytes, 24);
    let descriptor_len = read_u32(bytes, 20);
    let fingerprint_len = read_u64(bytes, 40);
    if end_prefix > PREFIX_END
        || covered_count > dictionary_count
        || covered_count > u64::from(end_prefix) << 14
        || dictionary_count > 1 << 30
        || (end_prefix == PREFIX_END && covered_count != dictionary_count)
        || u64::from(descriptor_len).checked_add(fingerprint_len)
            != Some((bytes.len() - METADATA_BYTES) as u64)
    {
        return Err(SharedError::Invalid("core filter extent"));
    }
    if covered_count == 0 {
        if descriptor_len != 0 || fingerprint_len != 0 {
            return Err(SharedError::Invalid("empty core filter"));
        }
    } else {
        if descriptor_len != 20 {
            return Err(SharedError::Invalid("core filter descriptor"));
        }
        validate_descriptor(&bytes[128..148], bytes.len() - 148)
            .map_err(|_| SharedError::Invalid("core filter descriptor"))?;
    }
    Ok((end_prefix, covered_count))
}

pub(crate) struct CoreFilter {
    bytes: Vec<u8>,
    end_prefix: u32,
    covered_count: u64,
    setup_ns: [u64; 3],
}

pub(crate) fn load(file: &SharedFile) -> Result<Option<CoreFilter>, SharedError> {
    if file.header.version != 4 {
        return Ok(None);
    }
    file.verify_unchanged()?;
    let length = file.header.section(Section::CoreFilter).length;
    if length > MAX_BYTES as u64 {
        return Err(SharedError::Invalid("core filter size"));
    }
    let observed = file.stats().observed;
    let started = observed.then(std::time::Instant::now);
    let bytes = file.section(
        Section::CoreFilter,
        0,
        file.header.section(Section::CoreFilter).length,
    )?;
    let authenticated = observed.then(std::time::Instant::now);
    let (end_prefix, covered_count) = validate(
        bytes,
        file.header.core_count,
        &file.header.manifest_sha256,
        &file.header.filter_source_sha256,
    )?;
    let boundary = file.section(Section::CorePrefixes, u64::from(end_prefix) * 4, 4)?;
    if u64::from(read_u32(boundary, 0)) != covered_count {
        return Err(SharedError::Invalid("core filter coverage"));
    }
    file.verify_unchanged()?;
    #[cfg(test)]
    let limit = FILTER_TEST_LIMIT.with(|limit| limit.get().min(MAX_BYTES));
    #[cfg(not(test))]
    let limit = MAX_BYTES;
    if bytes.len() > limit {
        return Err(SharedError::ResourceLimit);
    }
    let validated = observed.then(std::time::Instant::now);
    let mut owned = Vec::new();
    owned
        .try_reserve_exact(bytes.len())
        .map_err(|_| SharedError::ResourceLimit)?;
    owned.extend_from_slice(bytes);
    Ok(Some(CoreFilter {
        bytes: owned,
        end_prefix,
        covered_count,
        setup_ns: match (started, authenticated, validated) {
            (Some(start), Some(auth), Some(valid)) => [
                (auth - start).as_nanos() as u64,
                (valid - auth).as_nanos() as u64,
                valid.elapsed().as_nanos() as u64,
            ],
            _ => [0; 3],
        },
    }))
}

impl CoreFilter {
    pub(crate) fn setup_ns(&self) -> [u64; 3] {
        self.setup_ns
    }

    #[cfg(test)]
    pub(crate) fn allocation(&self) -> usize {
        self.bytes.as_ptr() as usize
    }

    pub(crate) fn view(&self) -> Option<BinaryFuse8Ref<'_>> {
        (self.covered_count != 0)
            .then(|| BinaryFuse8Ref::from_dma(&self.bytes[128..148], &self.bytes[148..]))
    }

    pub(crate) fn end_prefix(&self) -> u32 {
        self.end_prefix
    }

    pub(crate) fn bytes(&self) -> usize {
        self.bytes.capacity()
    }

    pub(crate) fn covered_count(&self) -> u64 {
        self.covered_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(bytes: &[u8], count: u64) -> Result<(u32, u64), SharedError> {
        validate(bytes, count, &[1; 32], &[2; 32])
    }

    #[test]
    fn empty_tiny_zero_and_borrowed_equality() {
        for keys in [vec![], vec![0], vec![0, 1, (1 << 30) - 1]] {
            let bytes = encode(&keys, PREFIX_END, keys.len() as u64, &[1; 32], &[2; 32]).unwrap();
            let (end_prefix, covered_count) = check(&bytes, keys.len() as u64).unwrap();
            let filter = CoreFilter {
                bytes,
                end_prefix,
                covered_count,
                setup_ns: [0; 3],
            };
            assert_eq!(filter.end_prefix(), PREFIX_END);
            assert_eq!(filter.covered_count(), keys.len() as u64);
            assert!(filter.bytes() <= MAX_BYTES);
            if keys.is_empty() {
                assert!(filter.view().is_none());
            } else {
                let owned = build_binary_fuse(&keys).unwrap();
                let borrowed = filter.view().unwrap();
                for key in keys.iter().copied().chain(0..10_000) {
                    assert_eq!(owned.contains(&key), borrowed.contains(&key));
                }
                for &key in &keys {
                    assert!(borrowed.contains(&key));
                }
            }
        }
        assert_eq!(core_input(0), 0);
        assert_eq!(core_input((1 << 30) - 1), (1 << 30) - 1);
    }

    #[test]
    fn metadata_generation_descriptor_and_truncation_fail_closed() {
        let valid = encode(&[0], PREFIX_END, 1, &[1; 32], &[2; 32]).unwrap();
        for offset in [0, 8, 12, 20, 24, 32, 40, 48, 80, 112, 136, 140, 144] {
            let mut invalid = valid.clone();
            invalid[offset] ^= 1;
            assert!(check(&invalid, 1).is_err(), "offset {offset}");
        }
        for length in 0..valid.len() {
            assert!(check(&valid[..length], 1).is_err(), "length {length}");
        }
        let mut extra = valid.clone();
        extra.push(0);
        assert!(check(&extra, 1).is_err());
        assert!(check(&valid, 2).is_err());
    }

    #[test]
    fn coverage_range_edges_are_checked() {
        for end in [0, 1, PREFIX_END - 1, PREFIX_END] {
            let bytes = encode(&[], end, 0, &[1; 32], &[2; 32]).unwrap();
            assert_eq!(check(&bytes, 0).unwrap(), (end, 0));
        }
        let bytes = encode(&[0, (1 << 14) - 1], 1, 3, &[1; 32], &[2; 32]).unwrap();
        assert_eq!(check(&bytes, 3).unwrap(), (1, 2));
        for end in [0, PREFIX_END, PREFIX_END + 1] {
            let mut invalid = bytes.clone();
            put_u32(&mut invalid, 16, end);
            assert!(check(&invalid, 3).is_err());
        }
    }
}
