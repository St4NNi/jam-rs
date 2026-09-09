use std::fmt;

pub const MAX_KEYS_PER_BLOCK: usize = 256;
pub const MEMBER_ANCHOR_STRIDE: u64 = 16;
pub const LONG_MEMBER_OCCURRENCES: u64 = 256;
const MAGIC: [u8; 8] = *b"JOWNBLK\0";
const VERSION: u16 = 2;
const HEADER_SIZE: usize = 40;
const MAX_HOT_DECODE_BYTES: usize = 64 * 1024 * 1024;

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
pub struct OwnerHotBlock {
    pub keys: Vec<OwnerHotKey>,
    pub members: Vec<OwnerHotMember>,
    pub anchors: Vec<OwnerAnchor>,
    pub cold_bits: u64,
    pub document_count: u32,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OwnerHotKey {
    pub key: u64,
    pub document_frequency: u64,
    pub member_start: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerHotMember {
    pub document_id: u32,
    pub occurrence_count: u64,
    pub member_ordinal: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerAnchor {
    pub member_ordinal: u64,
    pub cold_bit_offset: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerMemberWindow {
    pub first_member: u64,
    pub target_member: u64,
    pub start_bit: u64,
    pub end_bit: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DecodedOwnerMember {
    pub loci: Vec<(u64, bool)>,
    pub skipped_members: u64,
    pub skipped_occurrences: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct OwnerPostingsError(pub &'static str);

impl fmt::Display for OwnerPostingsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "invalid owner posting block: {}", self.0)
    }
}

impl std::error::Error for OwnerPostingsError {}

pub fn encode_block<F>(
    keys: &[OwnerKey],
    document_widths: &[u8],
    mut to_locus: F,
) -> Result<EncodedOwnerBlock, OwnerPostingsError>
where
    F: FnMut(u32, OwnerOccurrence) -> Result<u64, OwnerPostingsError>,
{
    if keys.is_empty()
        || keys.len() > MAX_KEYS_PER_BLOCK
        || document_widths.is_empty()
        || document_widths.len() > u32::MAX as usize
        || document_widths.iter().any(|width| *width > 64)
    {
        return Err(invalid("block input"));
    }
    let mut body = Vec::new();
    let mut cold = BitWriter::default();
    let mut anchors = Vec::new();
    let mut member_ordinal = 0u64;
    let mut previous_key = None;
    for entry in keys {
        put_varint(&mut body, delta_u64(previous_key, entry.key, "key order")?);
        let df = u64::try_from(entry.members.len()).map_err(|_| invalid("document frequency"))?;
        if df == 0 || df > document_widths.len() as u64 {
            return Err(invalid("document frequency"));
        }
        put_varint(&mut body, df);
        let (document_ids, dense) = encode_document_ids(&entry.members, document_widths.len())?;
        let all_singleton = entry
            .members
            .iter()
            .all(|member| member.occurrences.len() == 1);
        body.push(u8::from(dense) | (u8::from(all_singleton) << 1));
        body.extend_from_slice(&document_ids);
        if !all_singleton {
            let mut non_singletons = vec![0; entry.members.len().div_ceil(8)];
            for (index, member) in entry.members.iter().enumerate() {
                if member.occurrences.len() != 1 {
                    set_bit(&mut non_singletons, index);
                }
            }
            body.extend_from_slice(&non_singletons);
            for member in &entry.members {
                if member.occurrences.len() != 1 {
                    put_varint(
                        &mut body,
                        u64::try_from(member.occurrences.len())
                            .map_err(|_| invalid("occurrence count"))?,
                    );
                }
            }
        }
        for member in &entry.members {
            let count =
                u64::try_from(member.occurrences.len()).map_err(|_| invalid("occurrence count"))?;
            if count == 0 {
                return Err(invalid("occurrence count"));
            }
            if member_ordinal.is_multiple_of(MEMBER_ANCHOR_STRIDE)
                || count > LONG_MEMBER_OCCURRENCES
            {
                push_anchor(&mut anchors, member_ordinal, cold.len())?;
            }
            let width = document_width(document_widths, member.document_id)?;
            encode_loci(member, width, &mut to_locus, &mut cold)?;
            member_ordinal = member_ordinal
                .checked_add(1)
                .ok_or_else(|| invalid("member count"))?;
            if count > LONG_MEMBER_OCCURRENCES {
                push_anchor(&mut anchors, member_ordinal, cold.len())?;
            }
        }
        previous_key = Some(entry.key);
    }
    push_anchor(&mut anchors, member_ordinal, cold.len())?;
    let anchor_count = u32::try_from(anchors.len()).map_err(|_| invalid("anchor count"))?;
    let mut previous_member = None;
    let mut previous_bit = None;
    for anchor in &anchors {
        put_varint(
            &mut body,
            delta_u64(previous_member, anchor.member_ordinal, "anchor order")?,
        );
        put_varint(
            &mut body,
            delta_u64(previous_bit, anchor.cold_bit_offset, "anchor order")?,
        );
        previous_member = Some(anchor.member_ordinal);
        previous_bit = Some(anchor.cold_bit_offset);
    }
    let mut hot = vec![0; HEADER_SIZE];
    hot[..8].copy_from_slice(&MAGIC);
    put_u16(&mut hot, 8, VERSION);
    put_u16(&mut hot, 10, keys.len() as u16);
    put_u32(&mut hot, 12, document_widths.len() as u32);
    put_u64(&mut hot, 16, member_ordinal);
    put_u64(&mut hot, 24, cold.len());
    put_u32(&mut hot, 32, anchor_count);
    hot.extend_from_slice(&body);
    Ok(EncodedOwnerBlock {
        hot,
        cold: cold.finish(),
    })
}

pub fn parse_hot(hot: &[u8]) -> Result<OwnerHotBlock, OwnerPostingsError> {
    let (key_count, document_count, member_count, cold_bits, anchor_count, mut input) =
        decode_header(hot)?;
    decoded_hot_bound(key_count, member_count, anchor_count)?;
    let mut keys = Vec::new();
    keys.try_reserve_exact(key_count)
        .map_err(|_| invalid("decoded hot size"))?;
    let member_capacity = usize::try_from(member_count).map_err(|_| invalid("member count"))?;
    let mut members = Vec::new();
    members
        .try_reserve_exact(member_capacity)
        .map_err(|_| invalid("decoded hot size"))?;
    let mut previous_key = None;
    for _ in 0..key_count {
        let key = apply_delta_u64(previous_key, take_varint(&mut input)?, "key order")?;
        let df = take_varint(&mut input)?;
        let remaining = member_count
            .checked_sub(members.len() as u64)
            .ok_or_else(|| invalid("member count"))?;
        if df == 0 || df > u64::from(document_count) || df > remaining {
            return Err(invalid("document frequency"));
        }
        decoded_hot_transient_bound(key_count, member_count, anchor_count, df)?;
        let flags = take_byte(&mut input)?;
        if flags & !3 != 0 {
            return Err(invalid("member flags"));
        }
        let ids = decode_document_ids(&mut input, df, document_count, flags & 1 != 0)?;
        let counts = decode_counts(&mut input, df, flags & 2 != 0)?;
        let member_start = members.len() as u64;
        for (document_id, occurrence_count) in ids.into_iter().zip(counts) {
            members.push(OwnerHotMember {
                document_id,
                occurrence_count,
                member_ordinal: members.len() as u64,
            });
        }
        keys.push(OwnerHotKey {
            key,
            document_frequency: df,
            member_start,
        });
        previous_key = Some(key);
    }
    if members.len() as u64 != member_count {
        return Err(invalid("member count"));
    }
    let mut anchors = Vec::new();
    anchors
        .try_reserve_exact(anchor_count)
        .map_err(|_| invalid("decoded hot size"))?;
    let mut previous_member = None;
    let mut previous_bit = None;
    for _ in 0..anchor_count {
        let member_ordinal =
            apply_delta_u64(previous_member, take_varint(&mut input)?, "anchor order")?;
        let cold_bit_offset =
            apply_delta_u64(previous_bit, take_varint(&mut input)?, "anchor order")?;
        anchors.push(OwnerAnchor {
            member_ordinal,
            cold_bit_offset,
        });
        previous_member = Some(member_ordinal);
        previous_bit = Some(cold_bit_offset);
    }
    if !input.is_empty()
        || anchors.first().copied()
            != Some(OwnerAnchor {
                member_ordinal: 0,
                cold_bit_offset: 0,
            })
        || anchors.last().copied()
            != Some(OwnerAnchor {
                member_ordinal: member_count,
                cold_bit_offset: cold_bits,
            })
        || anchors.iter().any(|anchor| {
            anchor.member_ordinal > member_count || anchor.cold_bit_offset > cold_bits
        })
    {
        return Err(invalid("anchor coverage"));
    }
    if !validate_anchor_ordinals(&anchors, &members)? {
        return Err(invalid("anchor set"));
    }
    Ok(OwnerHotBlock {
        keys,
        members,
        anchors,
        cold_bits,
        document_count,
    })
}

pub fn find_key(hot: &OwnerHotBlock, key: u64) -> Option<&OwnerHotKey> {
    hot.keys
        .binary_search_by_key(&key, |entry| entry.key)
        .ok()
        .map(|index| &hot.keys[index])
}

pub fn key_members<'a>(
    hot: &'a OwnerHotBlock,
    key: &OwnerHotKey,
) -> Result<&'a [OwnerHotMember], OwnerPostingsError> {
    let start = usize::try_from(key.member_start).map_err(|_| invalid("member range"))?;
    let length = usize::try_from(key.document_frequency).map_err(|_| invalid("member range"))?;
    hot.members
        .get(
            start
                ..start
                    .checked_add(length)
                    .ok_or_else(|| invalid("member range"))?,
        )
        .ok_or_else(|| invalid("member range"))
}

pub fn locate_member(
    hot: &OwnerHotBlock,
    member_ordinal: u64,
) -> Result<OwnerMemberWindow, OwnerPostingsError> {
    if member_ordinal >= hot.members.len() as u64 {
        return Err(invalid("member ordinal"));
    }
    let after = hot
        .anchors
        .partition_point(|anchor| anchor.member_ordinal <= member_ordinal);
    let before = after
        .checked_sub(1)
        .ok_or_else(|| invalid("anchor range"))?;
    let first = hot
        .anchors
        .get(before)
        .ok_or_else(|| invalid("anchor range"))?;
    let end = hot
        .anchors
        .get(after)
        .ok_or_else(|| invalid("anchor range"))?;
    Ok(OwnerMemberWindow {
        first_member: first.member_ordinal,
        target_member: member_ordinal,
        start_bit: first.cold_bit_offset,
        end_bit: end.cold_bit_offset,
    })
}

pub fn decode_member_window(
    cold_window: &[u8],
    window_byte_start: u64,
    window: OwnerMemberWindow,
    hot: &OwnerHotBlock,
    document_widths: &[u8],
    max_decoded_bytes: usize,
) -> Result<DecodedOwnerMember, OwnerPostingsError> {
    if locate_member(hot, window.target_member)? != window
        || document_widths.len() != hot.document_count as usize
    {
        return Err(invalid("member window"));
    }
    let expected_start = window.start_bit / 8;
    let expected_end = window.end_bit.div_ceil(8);
    if window_byte_start != expected_start
        || u64::try_from(cold_window.len()).ok() != Some(expected_end - expected_start)
    {
        return Err(invalid("cold window"));
    }
    let mut reader = BitReader::window(
        cold_window,
        window.start_bit - expected_start * 8,
        window.end_bit - expected_start * 8,
    )?;
    let mut skipped_occurrences = 0u64;
    let mut loci = Vec::new();
    for ordinal in window.first_member..=window.target_member {
        let member = *hot
            .members
            .get(ordinal as usize)
            .ok_or_else(|| invalid("member ordinal"))?;
        let width = document_width(document_widths, member.document_id)?;
        if ordinal == window.target_member {
            loci = decode_loci(
                &mut reader,
                member.occurrence_count,
                width,
                max_decoded_bytes,
            )?;
        } else {
            if member.occurrence_count > LONG_MEMBER_OCCURRENCES {
                return Err(invalid("unanchored long member"));
            }
            skip_loci(
                &mut reader,
                member.occurrence_count,
                width,
                max_decoded_bytes,
            )?;
            skipped_occurrences = skipped_occurrences
                .checked_add(member.occurrence_count)
                .ok_or_else(|| invalid("skipped occurrences"))?;
        }
    }
    let after = hot
        .anchors
        .partition_point(|anchor| anchor.member_ordinal <= window.target_member);
    let target_end = window
        .target_member
        .checked_add(1)
        .ok_or_else(|| invalid("member ordinal"))?;
    if hot.anchors[after].member_ordinal == target_end
        && reader.position() != window.end_bit - expected_start * 8
    {
        return Err(invalid("member boundary"));
    }
    if window.end_bit == hot.cold_bits {
        validate_tail_padding(cold_window, window.end_bit - expected_start * 8)?;
    }
    Ok(DecodedOwnerMember {
        loci,
        skipped_members: window.target_member - window.first_member,
        skipped_occurrences,
    })
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
