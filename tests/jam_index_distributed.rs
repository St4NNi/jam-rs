use jam_rs::jam_index::{
    IndexPlanSource, MetagenomeSource, ScreenSelectionPolicy, build_fragment, finalize_index,
    load_manifest, merge_part, plan_index, write_part,
};
use jam_rs::provenance::sha256_file;
use jam_rs::reader::JamReader;
use jam_rs::writer::{HashSampleInput, build_hash_samples};
use std::collections::BTreeMap;

fn sequence(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 7;
            state ^= state >> 9;
            state ^= state << 8;
            b"ACGT"[(state as usize) & 3]
        })
        .collect()
}

#[test]
fn fragment_merge_matches_reference_part_bytes_and_finalizes_once() {
    let source_root = tempfile::Builder::new()
        .prefix("jam-index-distributed-sources-")
        .tempdir_in("target")
        .unwrap();
    let output_root = tempfile::Builder::new()
        .prefix("jam-index-distributed-output-")
        .tempdir_in("target")
        .unwrap();
    let staged_root = tempfile::Builder::new()
        .prefix("jam-index-distributed-staged-")
        .tempdir_in("target")
        .unwrap();
    let mut sources = Vec::new();
    let mut planned = Vec::new();
    let mut staged = BTreeMap::new();
    let mut published = BTreeMap::new();
    for index in 0..6usize {
        let id = format!("mg-{index:02}");
        let path = source_root.path().join(format!("{id}.fasta"));
        let first = sequence(400 + index * 31, 0x1234_0000 + index as u64);
        let second = sequence(240 + index * 17, 0xabcd_0000 + index as u64);
        std::fs::write(
            &path,
            format!(
                ">first\n{}\n>second\n{}\n",
                String::from_utf8_lossy(&first),
                String::from_utf8_lossy(&second)
            ),
        )
        .unwrap();
        let source = MetagenomeSource {
            metagenome_id: id.clone(),
            sequence_path: path.clone(),
        };
        let staged_path = staged_root.path().join(format!("{id}.fasta"));
        std::fs::copy(&path, &staged_path).unwrap();
        staged.insert(
            id.clone(),
            MetagenomeSource {
                metagenome_id: id.clone(),
                sequence_path: staged_path,
            },
        );
        published.insert(id.clone(), source.clone());
        sources.push(source);
        planned.push(IndexPlanSource {
            metagenome_id: id,
            source_size: std::fs::metadata(&path).unwrap().len(),
            source_sha256: sha256_file(&path).unwrap(),
            source_path: path,
            estimated_bases: 0,
            estimated_signatures: 0,
        });
    }
    let policy = ScreenSelectionPolicy::default_signatures();
    let plan = plan_index(planned, "b".repeat(64), policy.clone(), 2, 2, 4).unwrap();
    let fragments_root = output_root.path().join("fragments");
    std::fs::create_dir(&fragments_root).unwrap();
    for fragment in plan.parts.iter().flat_map(|part| &part.fragments) {
        build_fragment(
            &plan,
            fragment.fragment_id,
            &staged,
            fragments_root.join(format!("fragment-{:06}", fragment.fragment_id)),
        )
        .unwrap();
    }
    let index = output_root.path().join("index");
    std::fs::create_dir(&index).unwrap();
    std::fs::create_dir(index.join("parts")).unwrap();
    for part in &plan.parts {
        let reference = output_root
            .path()
            .join(format!("reference-{:06}", part.part_id));
        std::fs::create_dir(&reference).unwrap();
        let ordered = part
            .fragments
            .iter()
            .flat_map(|fragment| &fragment.sources)
            .map(|source| published.get(&source.metagenome_id).unwrap().clone())
            .collect::<Vec<_>>();
        let reference_data = reference.join("part.bin");
        let reference_screen = reference.join("screen.jam");
        let repeated_screen = reference.join("screen-repeat.jam");
        let result = write_part(&reference_data, &ordered, &policy).unwrap();
        let samples = result
            .screen_samples
            .iter()
            .map(|sample| HashSampleInput {
                sample_name: sample.metagenome_id.clone(),
                hashes: sample.hashes.clone(),
            })
            .collect::<Vec<_>>();
        build_hash_samples(&reference_screen, &samples, 21, 1).unwrap();
        build_hash_samples(&repeated_screen, &samples, 21, 1).unwrap();
        assert_eq!(
            sha256_file(&reference_screen).unwrap(),
            sha256_file(&repeated_screen).unwrap(),
            "reference screen construction must be byte deterministic"
        );
        let output = index.join(format!("parts/part-{:06}", part.part_id));
        let merged = merge_part(&plan, part.part_id, &fragments_root, &output).unwrap();
        assert_eq!(
            sha256_file(&reference_data).unwrap(),
            sha256_file(&output.join("part.bin")).unwrap()
        );
        let reference_reader = JamReader::open(&reference_screen).unwrap();
        let merged_reader = JamReader::open(output.join("screen.jam")).unwrap();
        assert_eq!(
            reference_reader.sample_names(),
            merged_reader.sample_names()
        );
        assert_eq!(
            reference_reader.stats().entry_count,
            merged_reader.stats().entry_count
        );
        for bucket in 0..256 {
            assert_eq!(
                reference_reader.bucket_entries(bucket),
                merged_reader.bucket_entries(bucket)
            );
        }
        let reference_bytes = std::fs::read(&reference_screen).unwrap();
        let merged_bytes = std::fs::read(output.join("screen.jam")).unwrap();
        if reference_bytes != merged_bytes {
            let first = reference_bytes
                .iter()
                .zip(&merged_bytes)
                .position(|(left, right)| left != right)
                .unwrap_or(reference_bytes.len().min(merged_bytes.len()));
            panic!(
                "screen bytes differ at {first}; reference={} merged={}; sample_sizes={:?}/{:?}",
                reference_bytes.len(),
                merged_bytes.len(),
                reference_reader.sample_sizes(),
                merged_reader.sample_sizes(),
            );
        }
        assert_eq!(
            merged.part.estimated_signature_count,
            result.estimated_signature_count
        );
        assert_eq!(
            JamReader::open(output.join("screen.jam"))
                .unwrap()
                .sample_names()
                .len(),
            ordered.len()
        );
    }
    let stats = finalize_index(&plan, &index).unwrap();
    assert_eq!(stats.total_metagenomes, sources.len() as u64);
    assert_eq!(load_manifest(&index).unwrap().parts.len(), 2);
    assert!(finalize_index(&plan, &index).is_err());
    let first_fragment = &plan.parts[0].fragments[0];
    assert!(
        build_fragment(
            &plan,
            first_fragment.fragment_id,
            &staged,
            fragments_root.join(format!("fragment-{:06}", first_fragment.fragment_id)),
        )
        .is_err()
    );
}
