use super::*;
use crate::index::MappedFile;
use std::collections::BTreeMap;

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

const KEY_BITS: u32 = 30;

/// Random seed locations per entry: key to sorted seed location values. Entry 2 is empty.
fn random_entries() -> Vec<BTreeMap<u64, Vec<u32>>> {
    let mut state = 0xfeed_u64;
    [3000usize, 1, 0, 40_000, 20]
        .iter()
        .map(|&keys| {
            let mut entry: BTreeMap<u64, Vec<u32>> = BTreeMap::new();
            for _ in 0..keys {
                let key = random(&mut state) >> (64 - KEY_BITS);
                let locations = entry.entry(key).or_default();
                for _ in 0..1 + random(&mut state).is_multiple_of(8) as usize * 3 {
                    locations.push((random(&mut state) % 500_000) as u32);
                }
                locations.sort_unstable();
                locations.dedup();
            }
            entry
        })
        .collect()
}

fn records(entry: &BTreeMap<u64, Vec<u32>>) -> Vec<(u64, u32)> {
    entry
        .iter()
        .flat_map(|(key, locations)| locations.iter().map(|location| (*key, *location)))
        .collect()
}

fn write(path: &Path, entries: &[BTreeMap<u64, Vec<u32>>]) {
    let mut writer = PositionsWriter::create(path, KEY_BITS, 10, entries.len()).unwrap();
    for entry in entries {
        let (section, bytes) = encode_section(&records(entry), KEY_BITS).unwrap();
        writer.push_section(&section, &bytes).unwrap();
    }
    writer.finish().unwrap();
}

fn open(path: &Path) -> Positions {
    Positions::open(Arc::new(MappedFile::open(path, false).unwrap())).unwrap()
}

#[test]
fn positions_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join(POSITIONS_FILE);
    let entries = random_entries();
    write(&path, &entries);
    let positions = open(&path);
    positions.verify_checksums().unwrap();
    assert_eq!(positions.entries(), entries.len() as u32);
    assert_eq!((positions.key_bits(), positions.bin_bits()), (KEY_BITS, 10));
    let mut state = 7u64;
    let mut out = Vec::new();
    for (index, entry) in entries.iter().enumerate() {
        let index = index as u32;
        let stored = records(entry);
        assert_eq!(
            positions.entry_locations(index).unwrap(),
            stored.len() as u64
        );
        for (key, locations) in entry {
            let hash = (key << (64 - KEY_BITS)) | (random(&mut state) >> KEY_BITS);
            out.clear();
            positions.locations(index, hash, &mut out).unwrap();
            let values: Vec<u32> = out.iter().map(|location| location.value()).collect();
            assert_eq!(&values, locations);
        }
        for _ in 0..2000 {
            let hash = random(&mut state);
            out.clear();
            positions.locations(index, hash, &mut out).unwrap();
            assert_eq!(
                out.is_empty(),
                !entry.contains_key(&(hash >> (64 - KEY_BITS)))
            );
        }
        let section = *positions.section(index).unwrap();
        let bytes = positions.section_bytes(index).unwrap();
        let mut decoded = Vec::new();
        decode_section(&section, &bytes, |key, location| {
            decoded.push((key, location))
        })
        .unwrap();
        assert_eq!(decoded, stored);
    }
    assert!(positions.locations(9, 0, &mut out).is_err());
    let location = SeedLocation {
        bin: 12345,
        reverse: true,
    };
    assert_eq!(SeedLocation::from_value(location.value()), location);
}

#[test]
fn sections_copied_into_another_file_read_the_same() {
    let directory = tempfile::tempdir().unwrap();
    let entries = random_entries();
    let first = directory.path().join("first");
    write(&first, &entries);
    let source = open(&first);
    let copy = directory.path().join("copy");
    let mut writer = PositionsWriter::create(&copy, KEY_BITS, 10, 2).unwrap();
    for entry in [3, 0] {
        let bytes = source.section_bytes(entry).unwrap();
        writer
            .push_section(source.section(entry).unwrap(), &bytes)
            .unwrap();
    }
    writer.finish().unwrap();
    let copied = open(&copy);
    let mut out = Vec::new();
    for (key, locations) in &entries[0] {
        out.clear();
        copied
            .locations(1, key << (64 - KEY_BITS), &mut out)
            .unwrap();
        assert_eq!(
            out.iter().map(|p| p.value()).collect::<Vec<_>>(),
            *locations
        );
    }
}

#[test]
fn damage_is_detected() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join(POSITIONS_FILE);
    write(&path, &random_entries());
    let original = std::fs::read(&path).unwrap();
    let damaged = directory.path().join("damaged");
    let table_end = HEADER_BYTES as usize + 5 * ROW_BYTES;
    for (at, opens) in [
        (2, false),
        (25, false),
        (140, false),
        (table_end + 100, true),
    ] {
        let mut bytes = original.clone();
        bytes[at] ^= 0x04;
        std::fs::remove_file(&damaged).ok();
        std::fs::write(&damaged, &bytes).unwrap();
        let source = Arc::new(MappedFile::open(&damaged, false).unwrap());
        match Positions::open(source) {
            Ok(positions) => {
                assert!(opens, "byte {at}");
                assert!(positions.verify_checksums().is_err());
            }
            Err(_) => assert!(!opens, "byte {at}"),
        }
    }
    std::fs::write(&damaged, &original[..original.len() - 8]).unwrap();
    let source = Arc::new(MappedFile::open(&damaged, false).unwrap());
    assert!(Positions::open(source).is_err());
}
