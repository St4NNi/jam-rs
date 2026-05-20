use super::*;
use std::collections::BTreeMap;

fn random(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// Writes a sketch of random seeds (key to color set id) and random color sets.
fn write_random(path: &Path, bucket_bits: u32) -> (BTreeMap<u64, u32>, Vec<Vec<u32>>) {
    let mut state = 0x5eed_u64 + u64::from(bucket_bits);
    let sets: Vec<Vec<u32>> = (0..300)
        .map(|_| {
            let mut set: Vec<u32> = (0..1 + random(&mut state) % 20)
                .map(|_| (random(&mut state) % 5000) as u32)
                .collect();
            set.sort_unstable();
            set.dedup();
            set
        })
        .collect();
    let key_bits = bucket_bits + 16;
    let mut seeds = BTreeMap::new();
    for _ in 0..20_000 {
        let key = random(&mut state) >> (64 - key_bits);
        // Mostly small ids, a few rare ones, as frequency ranking gives.
        let id = match random(&mut state) % 10 {
            0 => random(&mut state) % 300,
            _ => random(&mut state) % 4,
        };
        seeds.insert(key, id as u32);
    }
    let mut writer = SketchWriter::create(path, bucket_bits, 5000).unwrap();
    let mut block = Vec::new();
    let mut bucket_seeds = Vec::new();
    for bucket in 0..1usize << bucket_bits {
        bucket_seeds.clear();
        let first = (bucket as u64) << 16;
        for (key, id) in seeds.range(first..first + (1 << 16)) {
            bucket_seeds.push((*key as u16, *id));
        }
        block.clear();
        encode_bucket(&bucket_seeds, &mut block);
        if !block.is_empty() {
            writer.push_bucket(bucket, &block).unwrap();
        }
    }
    let encoded: Vec<Vec<u8>> = sets
        .iter()
        .map(|set| {
            let mut out = Vec::new();
            encode_color_set(set, &mut out);
            out
        })
        .collect();
    writer.finish(&encoded).unwrap();
    (seeds, sets)
}

#[test]
fn sketch_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    for bucket_bits in [1, 6, 12] {
        let path = directory.path().join(format!("sketch-{bucket_bits}"));
        let (seeds, sets) = write_random(&path, bucket_bits);
        let sketch = Sketch::open(&path).unwrap();
        sketch.verify_checksums().unwrap();
        assert_eq!(sketch.seeds(), seeds.len() as u64);
        assert_eq!(sketch.entries(), 5000);
        assert_eq!(sketch.color_set_count(), sets.len() as u64);
        let key_bits = bucket_bits + 16;
        let mut state = 99u64;
        for (key, id) in &seeds {
            // Any hash with the stored key finds it, whatever its lower bits.
            let hash = (key << (64 - key_bits)) | (random(&mut state) >> key_bits);
            assert_eq!(sketch.color_set(hash), Some(ColorSetId(*id)));
        }
        for _ in 0..20_000 {
            let hash = random(&mut state);
            let stored = seeds.get(&(hash >> (64 - key_bits)));
            assert_eq!(sketch.color_set(hash), stored.map(|id| ColorSetId(*id)));
        }
        let color_sets = sketch.color_sets();
        for (id, set) in sets.iter().enumerate() {
            let entries = color_sets.entries(ColorSetId(id as u32)).unwrap();
            assert_eq!(entries.len(), set.len() as u64);
            assert_eq!(entries.collect::<Vec<_>>(), *set);
        }
        assert!(color_sets.entries(ColorSetId(sets.len() as u32)).is_none());
    }
}

#[test]
fn rare_ids_are_escaped_instead_of_widening_every_value() {
    let mut ids = vec![1u32; 100];
    ids[17] = 70_000_000;
    ids[60] = 3;
    let (width, escapes) = choose_width(&ids);
    assert_eq!((width, escapes), (2, 2));
    assert_eq!(choose_width(&[0, 0, 0]), (0, 0));
    assert_eq!(choose_width(&[u32::MAX]), (32, 0));
    let seeds: Vec<(u16, u32)> = ids
        .iter()
        .enumerate()
        .map(|(i, id)| (i as u16 * 3, *id))
        .collect();
    let mut block = Vec::new();
    encode_bucket(&seeds, &mut block);
    block.extend_from_slice(&[0; BIT_PADDING]);
    for (fingerprint, id) in &seeds {
        assert_eq!(find_in_block(&block, *fingerprint), Some(ColorSetId(*id)));
        assert_eq!(find_in_block(&block, fingerprint + 1), None);
    }
}

#[test]
fn damage_is_detected() {
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join("sketch");
    write_random(&path, 6);
    let original = std::fs::read(&path).unwrap();
    let damaged = directory.path().join("damaged");
    for (at, expect_open) in [
        (3, false),
        (30, false),
        (200, true),
        (original.len() - 2, true),
    ] {
        let mut bytes = original.clone();
        bytes[at] ^= 0x10;
        std::fs::remove_file(&damaged).ok();
        std::fs::write(&damaged, &bytes).unwrap();
        match Sketch::open(&damaged) {
            Ok(sketch) => {
                assert!(expect_open, "byte {at}");
                assert!(sketch.verify_checksums().is_err(), "byte {at}");
            }
            Err(_) => assert!(!expect_open, "byte {at}"),
        }
    }
    std::fs::write(&damaged, &original[..original.len() - 1]).unwrap();
    assert!(Sketch::open(&damaged).is_err());
}
