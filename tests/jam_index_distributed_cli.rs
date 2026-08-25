use jam_rs::jam_index::{IndexBuildPlan, load_manifest};
use jam_rs::provenance::sha256_file;
use std::process::Command;

#[test]
fn distributed_index_commands_publish_only_after_all_parts_exist() {
    let root = tempfile::Builder::new()
        .prefix("jam-index-distributed-cli-")
        .tempdir_in("target")
        .unwrap();
    let mut catalog = String::from("metagenome_id\tresource_uri\tsha256\n");
    for index in 0..4usize {
        let id = format!("mg-{index}");
        let path = root.path().join(format!("{id}.fasta"));
        let sequence = (0..600 + index * 31)
            .map(|position| b"ACGT"[(position * 17 + position / 11 + index) & 3] as char)
            .collect::<String>();
        std::fs::write(&path, format!(">contig-{index}\n{sequence}\n")).unwrap();
        catalog.push_str(&format!(
            "{id}\t{}\t{}\n",
            path.display(),
            sha256_file(&path).unwrap()
        ));
    }
    let catalog_path = root.path().join("catalog.tsv");
    std::fs::write(&catalog_path, catalog).unwrap();
    let index = root.path().join("index");
    let fragments = index.join("fragments");
    let parts = index.join("parts");
    std::fs::create_dir_all(&fragments).unwrap();
    std::fs::create_dir_all(&parts).unwrap();
    let plan_path = index.join("plan.json");
    let plan = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "index", "plan", "--metagenomes"])
        .arg(&catalog_path)
        .args(["--output"])
        .arg(&plan_path)
        .args([
            "--parts",
            "2",
            "--fragments-per-part",
            "2",
            "--screen-policy",
            "spatial256-one",
            "--adaptive-second-minimum-threshold",
            "768",
            "--whole-metagenome-hashes",
            "512",
        ])
        .output()
        .unwrap();
    assert!(
        plan.status.success(),
        "{}",
        String::from_utf8_lossy(&plan.stderr)
    );
    let plan: IndexBuildPlan = serde_json::from_reader(std::io::BufReader::new(
        std::fs::File::open(&plan_path).unwrap(),
    ))
    .unwrap();
    assert_eq!(
        plan.selection_policy.adaptive_second_minimum_bases,
        Some(768)
    );
    for fragment in plan.parts.iter().flat_map(|part| &part.fragments) {
        let run = Command::new(env!("CARGO_BIN_EXE_jam"))
            .args(["--silent", "index", "build-fragment", "--plan"])
            .arg(&plan_path)
            .args(["--fragment-id", &fragment.fragment_id.to_string()])
            .args(["--staged-metagenomes"])
            .arg(&catalog_path)
            .args(["--output"])
            .arg(fragments.join(format!("fragment-{:06}", fragment.fragment_id)))
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    for part in &plan.parts {
        let run = Command::new(env!("CARGO_BIN_EXE_jam"))
            .args(["--silent", "index", "merge-part", "--plan"])
            .arg(&plan_path)
            .args(["--part-id", &part.part_id.to_string()])
            .args(["--fragments-root"])
            .arg(&fragments)
            .args(["--output"])
            .arg(parts.join(format!("part-{:06}", part.part_id)))
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    assert!(!index.join("manifest.json").exists());
    let finalize = Command::new(env!("CARGO_BIN_EXE_jam"))
        .args(["--silent", "index", "finalize", "--plan"])
        .arg(&plan_path)
        .args(["--output"])
        .arg(&index)
        .output()
        .unwrap();
    assert!(
        finalize.status.success(),
        "{}",
        String::from_utf8_lossy(&finalize.stderr)
    );
    let manifest = load_manifest(&index).unwrap();
    assert_eq!(manifest.parts.len(), 2);
    assert_eq!(manifest.total_metagenomes, 4);
}
