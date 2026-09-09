use std::fmt;

pub const MAX_KEYS_PER_BLOCK: usize = 256;
const MAGIC: [u8; 8] = *b"JOWNBLK\0";
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 24;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EncodedOwnerBlock {
    pub hot: Vec<u8>,
    pub cold: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerKey {
    pub key: u64,
    pub members: Vec<OwnerMember>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerMember {
    pub document_id: u32,
    pub occurrences: Vec<OwnerOccurrence>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerOccurrence {
    pub local_contig: u32,
    pub position: u64,
    pub canonical_orientation: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerHotKey {
    pub key: u64,
    pub document_frequency: u64,
    pub members: Vec<OwnerHotMember>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerHotMember {
    pub document_id: u32,
    pub occurrence_count: u64,
    pub cold_offset: u64,
    pub cold_length: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerPostingsError(&'static str);

impl fmt::Display for OwnerPostingsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "invalid owner posting block: {}", self.0)
    }
}

impl std::error::Error for OwnerPostingsError {}

pub fn encode_block(keys: &[OwnerKey]) -> Result<EncodedOwnerBlock, OwnerPostingsError> {
    if keys.is_empty() || keys.len() > MAX_KEYS_PER_BLOCK {
        return Err(invalid("key count"));
    }
    let mut directory = Vec::new();
    let mut cold = Vec::new();
    let mut previous_key = None;
    for entry in keys {
        let key_code = delta_u64(previous_key, entry.key, "key order")?;
        put_varint(&mut directory, key_code);
        let document_frequency =
            u64::try_from(entry.members.len()).map_err(|_| invalid("document frequency"))?;
        if document_frequency == 0 {
            return Err(invalid("document frequency"));
        }
        put_varint(&mut directory, document_frequency);
        let mut previous_document = None;
        for member in &entry.members {
            let document_code = delta_u32(previous_document, member.document_id, "document order")?;
            put_varint(&mut directory, u64::from(document_code));
            let count =
                u64::try_from(member.occurrences.len()).map_err(|_| invalid("occurrence count"))?;
            if count == 0 {
                return Err(invalid("occurrence count"));
            }
            let start = cold.len();
            encode_occurrences(&member.occurrences, &mut cold)?;
            let length = u64::try_from(cold.len() - start).map_err(|_| invalid("cold length"))?;
            put_varint(&mut directory, count);
            put_varint(&mut directory, length);
            previous_document = Some(member.document_id);
        }
        previous_key = Some(entry.key);
    }

    let directory_length = u64::try_from(directory.len()).map_err(|_| invalid("directory"))?;
    let cold_length = u64::try_from(cold.len()).map_err(|_| invalid("cold data"))?;
    let capacity = HEADER_SIZE
        .checked_add(directory.len())
        .ok_or_else(|| invalid("block length"))?;
    let mut hot = Vec::with_capacity(capacity);
    hot.extend_from_slice(&MAGIC);
    hot.extend_from_slice(&VERSION.to_le_bytes());
    hot.extend_from_slice(&(keys.len() as u16).to_le_bytes());
    hot.extend_from_slice(&0u32.to_le_bytes());
    hot.extend_from_slice(&cold_length.to_le_bytes());
    hot.extend_from_slice(&directory);
    debug_assert_eq!(directory_length as usize, hot.len() - HEADER_SIZE);
    Ok(EncodedOwnerBlock { hot, cold })
}

pub fn parse_hot(hot: &[u8]) -> Result<Vec<OwnerHotKey>, OwnerPostingsError> {
    let mut directory = directory(hot)?;
    let key_count = usize::from(u16::from_le_bytes(
        hot[10..12].try_into().expect("owner key count"),
    ));
    let cold_length = declared_cold_length(hot)?;
    let mut cold_offset = 0u64;
    let mut keys = Vec::with_capacity(key_count);
    let mut previous_key = None;
    for _ in 0..key_count {
        let key = apply_delta_u64(previous_key, take_varint(&mut directory)?, "key order")?;
        let document_frequency = take_varint(&mut directory)?;
        let member_count =
            usize::try_from(document_frequency).map_err(|_| invalid("document frequency"))?;
        if member_count == 0 || member_count > directory.len() / 3 {
            return Err(invalid("document frequency"));
        }
        let mut members = Vec::with_capacity(member_count);
        let mut previous_document = None;
        for _ in 0..member_count {
            let document_code = take_varint(&mut directory)?;
            let document_id = apply_delta_u32(previous_document, document_code, "document order")?;
            let occurrence_count = take_varint(&mut directory)?;
            let cold_length_for_member = take_varint(&mut directory)?;
            if occurrence_count == 0 || cold_length_for_member < occurrence_count.saturating_mul(2)
            {
                return Err(invalid("member length"));
            }
            let next_cold = cold_offset
                .checked_add(cold_length_for_member)
                .ok_or_else(|| invalid("cold range"))?;
            if next_cold > cold_length {
                return Err(invalid("cold range"));
            }
            members.push(OwnerHotMember {
                document_id,
                occurrence_count,
                cold_offset,
                cold_length: cold_length_for_member,
            });
            cold_offset = next_cold;
            previous_document = Some(document_id);
        }
        keys.push(OwnerHotKey {
            key,
            document_frequency,
            members,
        });
        previous_key = Some(key);
    }
    if !directory.is_empty() || cold_offset != cold_length {
        return Err(invalid("block coverage"));
    }
    Ok(keys)
}

pub fn decode_block(hot: &[u8], cold: &[u8]) -> Result<Vec<OwnerKey>, OwnerPostingsError> {
    let directory = parse_hot(hot)?;
    let declared_cold = declared_cold_length(hot)?;
    if u64::try_from(cold.len()).ok() != Some(declared_cold) {
        return Err(invalid("cold data"));
    }
    let mut keys = Vec::with_capacity(directory.len());
    for entry in directory {
        let mut members = Vec::with_capacity(entry.members.len());
        for member in entry.members {
            members.push(OwnerMember {
                document_id: member.document_id,
                occurrences: decode_member(cold, member)?,
            });
        }
        keys.push(OwnerKey {
            key: entry.key,
            members,
        });
    }
    Ok(keys)
}

fn delta_u64(
    previous: Option<u64>,
    value: u64,
    field: &'static str,
) -> Result<u64, OwnerPostingsError> {
    match previous {
        None => Ok(value),
        Some(previous) => value
            .checked_sub(previous)
            .filter(|delta| *delta != 0)
            .ok_or_else(|| invalid(field)),
    }
}

fn delta_u32(
    previous: Option<u32>,
    value: u32,
    field: &'static str,
) -> Result<u32, OwnerPostingsError> {
    match previous {
        None => Ok(value),
        Some(previous) => value
            .checked_sub(previous)
            .filter(|delta| *delta != 0)
            .ok_or_else(|| invalid(field)),
    }
}

fn apply_delta_u64(
    previous: Option<u64>,
    delta: u64,
    field: &'static str,
) -> Result<u64, OwnerPostingsError> {
    match previous {
        None => Ok(delta),
        Some(previous) if delta != 0 => previous.checked_add(delta).ok_or_else(|| invalid(field)),
        Some(_) => Err(invalid(field)),
    }
}

fn apply_delta_u32(
    previous: Option<u32>,
    delta: u64,
    field: &'static str,
) -> Result<u32, OwnerPostingsError> {
    let delta = u32::try_from(delta).map_err(|_| invalid(field))?;
    match previous {
        None => Ok(delta),
        Some(previous) if delta != 0 => previous.checked_add(delta).ok_or_else(|| invalid(field)),
        Some(_) => Err(invalid(field)),
    }
}

fn put_varint(output: &mut Vec<u8>, value: u64) {
    let mut buffer = unsigned_varint::encode::u64_buffer();
    output.extend_from_slice(unsigned_varint::encode::u64(value, &mut buffer));
}

fn take_varint(input: &mut &[u8]) -> Result<u64, OwnerPostingsError> {
    let before = *input;
    let (value, remaining) = unsigned_varint::decode::u64(before).map_err(|_| invalid("varint"))?;
    let consumed = before.len() - remaining.len();
    let mut buffer = unsigned_varint::encode::u64_buffer();
    if unsigned_varint::encode::u64(value, &mut buffer) != &before[..consumed] {
        return Err(invalid("varint"));
    }
    *input = remaining;
    Ok(value)
}

const fn invalid(field: &'static str) -> OwnerPostingsError {
    OwnerPostingsError(field)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> Vec<OwnerKey> {
        vec![
            OwnerKey {
                key: 7,
                members: vec![OwnerMember {
                    document_id: 3,
                    occurrences: vec![OwnerOccurrence {
                        local_contig: 4,
                        position: u64::from(u32::MAX) + 9,
                        canonical_orientation: true,
                    }],
                }],
            },
            OwnerKey {
                key: u64::MAX,
                members: vec![
                    OwnerMember {
                        document_id: 0,
                        occurrences: vec![
                            OwnerOccurrence {
                                local_contig: 1,
                                position: 5,
                                canonical_orientation: false,
                            },
                            OwnerOccurrence {
                                local_contig: 1,
                                position: 5,
                                canonical_orientation: true,
                            },
                            OwnerOccurrence {
                                local_contig: u32::MAX,
                                position: u64::MAX,
                                canonical_orientation: false,
                            },
                        ],
                    },
                    OwnerMember {
                        document_id: u32::MAX,
                        occurrences: vec![OwnerOccurrence {
                            local_contig: 0,
                            position: 0,
                            canonical_orientation: false,
                        }],
                    },
                ],
            },
        ]
    }

    #[test]
    fn round_trip_preserves_wide_values_and_duplicate_positions() {
        let expected = fixture();
        let encoded = encode_block(&expected).unwrap();
        assert_eq!(decode_block(&encoded.hot, &encoded.cold).unwrap(), expected);
    }

    #[test]
    fn hot_lookup_does_not_require_valid_cold_varints() {
        let encoded = encode_block(&fixture()).unwrap();
        let hot = lookup_hot(&encoded.hot, u64::MAX).unwrap().unwrap();
        assert_eq!(hot.document_frequency, 2);
        assert_eq!(hot.members[0].occurrence_count, 3);
        let mut broken = encoded.cold;
        *broken.last_mut().unwrap() = 0x80;
        assert!(lookup_hot(&encoded.hot, u64::MAX).unwrap().is_some());
        assert!(decode_block(&encoded.hot, &broken).is_err());
    }

    #[test]
    fn rejects_bad_order_noncanonical_varints_and_ranges() {
        let mut keys = fixture();
        keys.swap(0, 1);
        assert!(encode_block(&keys).is_err());
        let mut duplicate_document = fixture();
        duplicate_document[1].members[1].document_id = 0;
        assert!(encode_block(&duplicate_document).is_err());

        let mut encoded = encode_block(&fixture()).unwrap();
        encoded.hot[16..24].copy_from_slice(&u64::MAX.to_le_bytes());
        assert!(decode_block(&encoded.hot, &encoded.cold).is_err());

        let mut encoded = encode_block(&fixture()).unwrap();
        let directory_start = HEADER_SIZE;
        encoded.hot[directory_start] = 0x87;
        encoded.hot.insert(directory_start + 1, 0);
        assert!(decode_block(&encoded.hot, &encoded.cold).is_err());
    }

    #[test]
    fn rejects_too_many_keys_and_descending_positions() {
        let keys = (0..=MAX_KEYS_PER_BLOCK)
            .map(|key| OwnerKey {
                key: key as u64,
                members: vec![OwnerMember {
                    document_id: 0,
                    occurrences: vec![OwnerOccurrence {
                        local_contig: 0,
                        position: 0,
                        canonical_orientation: false,
                    }],
                }],
            })
            .collect::<Vec<_>>();
        assert!(encode_block(&keys).is_err());
        let mut descending = fixture();
        descending[1].members[0].occurrences[1].position = 4;
        assert!(encode_block(&descending).is_err());
    }
}
