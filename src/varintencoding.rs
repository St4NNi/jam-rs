use std::borrow::Cow;

use heed::BoxedError;
use integer_encoding::{VarInt, VarIntReader};

pub struct VarIntEncoder;

impl heed::BytesEncode<'_> for VarIntEncoder {
    type EItem = Vec<u32>;

    fn bytes_encode(item: &Self::EItem) -> Result<Cow<'_, [u8]>, BoxedError> {
        let mut vec = Vec::new();
        for integer in item {
            vec.extend_from_slice(&integer.encode_var_vec());
        }
        Ok(Cow::Owned(vec))
    }
}

impl heed::BytesDecode<'_> for VarIntEncoder {
    type DItem = Vec<u32>;
    fn bytes_decode(bytes: &[u8]) -> Result<Self::DItem, BoxedError> {
        let mut vec = Vec::new();
        let mut bytes = bytes;
        while !bytes.is_empty() {
            vec.push(VarIntReader::read_varint(&mut bytes).map_err(|e| e.to_string())?);
        }
        Ok(vec)
    }
}
