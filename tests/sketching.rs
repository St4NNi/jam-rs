use jam_rs::{file_io::FileHandler, signature::Signature};
use sourmash::sketch::Sketch;
use std::path::{self, PathBuf};

fn get_hashes_sketch(sketch: &Sketch) -> Vec<u64> {
    if let Sketch::MinHash(minhash) = sketch {
        minhash.mins()
    } else {
        panic!("Sketch is not a MinHash sketch");
    }
}

#[test]
fn test_file_sketching_basic() {
    let input_file = "tests/testfiles/test.small.fa";
    FileHandler::sketch_files(
        jam_rs::cli::Commands::Sketch {
            input: vec![PathBuf::from(input_file)],
            output: Some(PathBuf::from("test.small.fa.test")),
            kmer_size: 33,
            fscale: None,
            kscale: None,
            nmin: None,
            nmax: None,
            stats: false,
            format: jam_rs::cli::OutputFormats::Sourmash,
            algorithm: jam_rs::cli::HashAlgorithms::Murmur3,
            singleton: false,
        },
        None,
    )
    .unwrap();

    let created_sketch =
        sourmash::signature::Signature::from_path(path::Path::new("test.small.fa.test"))
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

#[test]
fn test_file_sketching_comp() {
    let input_file = "tests/testfiles/test.small.fa";
    FileHandler::sketch_files(
        jam_rs::cli::Commands::Sketch {
            input: vec![PathBuf::from(input_file)],
            output: Some(PathBuf::from("test.small.fa.test.bin")),
            kmer_size: 33,
            fscale: None,
            kscale: None,
            nmin: None,
            nmax: None,
            stats: false,
            format: jam_rs::cli::OutputFormats::Bin,
            algorithm: jam_rs::cli::HashAlgorithms::Murmur3,
            singleton: false,
        },
        None,
    )
    .unwrap();

    let read_to_bytes = std::fs::read("test.small.fa.test.bin").unwrap();
    let mut signature: Vec<Signature> =
        bincode::deserialize_from(read_to_bytes.as_slice()).unwrap();
    let signature = signature.pop().unwrap();

    let expected_signature = sourmash::signature::Signature::from_path(path::Path::new(
        "tests/testfiles/test.small.fasta.sourmash_k33.sig",
    ))
    .unwrap()
    .pop()
    .unwrap();

    let expected_signature = jam_rs::signature::Signature::from(expected_signature);

    assert_eq!(signature.max_hash, expected_signature.max_hash);
    assert_eq!(signature.kmer_size, expected_signature.kmer_size);
    assert_eq!(signature.sketches.len(), expected_signature.sketches.len());
    assert_eq!(
        signature.sketches[0].hashes,
        expected_signature.sketches[0].hashes
    );
}
