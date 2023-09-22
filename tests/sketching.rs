use std::path::{self, PathBuf};

use jam_rs::{compare::MultiComp, file_io::FileHandler};
use sourmash::sketch::Sketch;

fn get_hashes_sketch(sketch: &Sketch) -> Vec<u64> {
    if let Sketch::MinHash(minhash) = sketch {
        minhash.mins()
    } else {
        panic!("Sketch is not a MinHash sketch");
    }
}

#[test]
fn test_file_sketching() {
    let input_file = "tests/testfiles/test.small.fa";
    FileHandler::sketch_files(
        jam_rs::cli::Commands::Sketch {
            input: vec![PathBuf::from(input_file)],
            output: Some(PathBuf::from("out.test")),
            kmer_size: 33,
            fscale: None,
            kscale: None,
            nmin: None,
            nmax: None,
            format: jam_rs::cli::OutputFormats::Sourmash,
            algorithm: jam_rs::cli::HashAlgorithms::Murmur3,
            singleton: false,
        },
        None,
    )
    .unwrap();

    let created_sketch = sourmash::signature::Signature::from_path(path::Path::new("out.test"))
        .unwrap()
        .pop()
        .unwrap()
        .sketches()
        .pop()
        .unwrap();

    let expected_sketch = sourmash::signature::Signature::from_path(path::Path::new(
        "tests/testfiles/test.small.fasta.sourmash_k33.sig",
    ))
    .unwrap()
    .pop()
    .unwrap()
    .sketches()
    .pop()
    .unwrap();

    assert_eq!(
        get_hashes_sketch(&created_sketch),
        get_hashes_sketch(&expected_sketch)
    );
}

// #[test]
// fn test_multi_comp() {

//     let comp = MultiComp::new
// }
