#[path = "support/remote/mod.rs"]
mod remote;

use jam_rs::jma::ArchiveReader;
use jam_rs::jma::reader::JmaReader;
use jam_rs::jma::writer::{ArchiveParts, encode_archive};
use jam_rs::jma::{ContigMetadata, SequenceRange};
use jam_rs::resource::ResourceOpenOptions;
use jam_rs::resource::object::ObjectResource;
use jam_rs::sequence::{BlockCodec, encode_sequence_block};
use remote::HttpFixture;

#[test]
fn remote_reader_uses_ranges_for_metadata_and_selected_sequence() {
    let sequence = b"ACGTACGTACGTACGT";
    let sequence_block = encode_sequence_block(0, 0, sequence, BlockCodec::Raw2Bit).unwrap();
    let bytes = encode_archive(&ArchiveParts {
        flags: 0,
        source_sha256: [5; 32],
        contigs: vec![ContigMetadata {
            id: 0,
            name: "remote".into(),
            length: sequence.len() as u64,
        }],
        sequence_blocks: vec![sequence_block],
        seed_sections: Vec::new(),
    })
    .unwrap();
    let fixture = HttpFixture::new(bytes);
    let resource = ObjectResource::open(
        fixture.url("/remote.jma"),
        ResourceOpenOptions {
            cache_block_bytes: 2,
            max_retries: 0,
            ..ResourceOpenOptions::default()
        },
    )
    .unwrap();
    let reader = JmaReader::open(resource).unwrap();
    assert_eq!(
        reader
            .read_sequence(0, SequenceRange::new(2, 8).unwrap())
            .unwrap(),
        b"GTACGT"
    );
    let requests = fixture.requests();
    let mut gets = requests.iter().filter(|request| request.method == "GET");
    assert!(gets.clone().all(|request| request.range.is_some()));
    assert!(gets.any(|request| request.requested_bytes < 256));
}
