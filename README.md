[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)

# jam-rs

`jam-rs` is a high-speed candidate finder for known and near-known plasmid traces in metagenomic assemblies. It reports contig-level and assembly-level sketch evidence against a fixed plasmid catalog. Optional target-aware sampling may improve retrieval of represented targets, but final presence calls require independent confirmation.

The intended workflow is:

```text
metagenomic assembly FASTA
    -> jam-rs screening against a versioned plasmid reference catalog
    -> candidate plasmid references and supporting contigs
    -> read mapping, marker checks, assembly-graph checks, or other confirmation
    -> final surveillance report
```

A sketch hit is an `alert`, not proof of plasmid presence. Shared mobile elements, repeats, contamination, conserved sequence, or chromosomal integration can all create genuine sketch overlap without an intact plasmid being present.

## Install

Install the current crates.io release:

```bash
cargo install jam-rs --locked
```

Or build this source snapshot:

```bash
cargo build --release --locked
install -m 0755 target/release/jam "$HOME/.local/bin/jam"
```

The container build is also tested:

```bash
docker build -t jam-rs:0.10.0 .
docker run --rm jam-rs:0.10.0 --version
```

## Quick start

Build one reference sketch per catalog FASTA record. Use stable, whitespace-free FASTA identifiers so candidate IDs can be passed to downstream tools.

```bash
jam sketch plasmid_catalog.fasta \
  --output plasmids.jam \
  --singleton \
  --kmer-size 21 \
  --fscale 100 \
  --threads 8 \
  --memory 4
```

This writes `plasmids.jam` and the provenance sidecar `plasmids.jam.json`.

Screen an assembly, treating each contig as one query:

```bash
jam screen \
  --input assembly.fasta.gz \
  --database plasmids.jam \
  --output contig_hits.tsv \
  --summary assembly_hits.tsv \
  --metadata screen_run.json \
  --min-shared 3 \
  --min-query-containment 0.0 \
  --min-reference-containment 0.0 \
  --top-per-contig 10 \
  --top-references 100 \
  --threads 8 \
  --memory 4
```

The filters are independent. Query containment is `shared_hashes / query_hashes`; reference containment is `shared_hashes / reference_hashes`. The assembly summary computes the union of shared hashes for each reference, so overlapping contigs cannot count the same reference hash more than once.

Contig output begins with:

```text
schema_version  assembly  query_contig  reference  shared_hashes  query_hashes  reference_hashes  query_containment  reference_containment  ...
```

Assembly output begins with:

```text
schema_version  assembly  reference  supporting_contigs  shared_hashes_union  reference_hashes  aggregate_reference_containment  ...
```

Both files always contain their stable header, including for an empty assembly
or a no-hit result. Screening remains a candidate stage; downstream trace,
mapping, or contextual evidence determines whether a candidate is supported.

Inspect a database without parsing human-readable text:

```bash
jam stats --input plasmids.jam --json
```

## Positional follow-up with `jam trace`

`jam trace` is an optional follow-up for one plasmid query across many candidate
metagenomes. It uses an existing `.jam` metagenome index, then retrieves
positional evidence from JMA v1 archives or raw assembly resources:

```bash
jam archive \
  --input assembly.fasta \
  --output assembly.jma \
  --primary-scale 200 \
  --rescue-scale 500

jam trace \
  --plasmid plasmid.fasta \
  --database metagenomes.jam \
  --catalog metagenomes.tsv \
  --output plasmid.trace.jsonl.zst \
  --sensitivity balanced \
  --min-shared 3 \
  --threads 8 \
  --memory 4
```

The primary positional seed is k=31 for longer exact anchors in a repetitive
assembly; k=21 is an optional rescue for shorter or more divergent represented
sequence. This is a fixed engineering contract, not evidence that k=31 is a
universally optimal k-mer size. Candidate `.jam` sketch parameters remain
independent and configurable.

Catalog JMA/raw resources may be local, `file://`, HTTP(S), or S3-compatible
object locators. A remote candidate `.jam` is downloaded into an explicit
identity-checked cache before mmap:

```bash
jam trace \
  --plasmid plasmid.fasta \
  --database https://objects.example.org/metagenomes.jam \
  --cache-dir cache/jam \
  --catalog metagenomes.tsv \
  --output plasmid.trace.jsonl
```

Alignment and nonredundant plasmid-coordinate coverage are `supported`
evidence, not confirmation of plasmid presence or fragment linkage. See
[docs/TRACE.md](docs/TRACE.md), [docs/REMOTE_RESOURCES.md](docs/REMOTE_RESOURCES.md),
and [examples/trace/run_trace.sh](examples/trace/run_trace.sh).

## Bias-assisted candidate retrieval

Uniform FracMinHash is the baseline. A bias table is an optional, catalog-specific sampler tied to its positive set, chromosome background, k-mer size, and sampling configuration. It is not a classifier or a probability that a k-mer came from a plasmid.

```bash
jam bias create \
  --positive catalog_training.fasta \
  --negative chromosome_background.fasta \
  --output catalog.bias \
  --kmer-size 21 \
  --fscale 100 \
  --min-positive-retention 0.25

jam sketch plasmid_catalog.fasta \
  --output plasmids.bias.jam \
  --singleton \
  --kmer-size 21 \
  --bias-table catalog.bias
```

Bias output uses separate `retained_*` and `bias_weighted_*` fields, includes
`score_mode=bias` and the bias-table checksum, and reports
`uniform_hash_e_value=NA`. It never filters on an E-value approximation. A
recommended use is bias-assisted candidate retrieval followed by uniform
sketch evidence and independent mapping or contextual confirmation.

## Evaluation

The deterministic harness under [evaluation/trace/](evaluation/trace/) is a
small reproducibility and smoke-test workload. It does not establish
production-scale accuracy, remote-storage economics, or a general speed
advantage. Release measurements must record complete inputs, truth, commands,
versions, hardware and storage, raw measurements, checksums, and summary
generation.

## Compatibility and provenance

- `.jam` binary format version remains 3. Existing version-3 files remain readable; the repository tests a fixture created by jam-rs 0.9.11.
- New builds add `<database>.jam.json`; this sidecar does not alter the binary database.
- The released hash identity is `jamhash_u64_v1`, provided by the exactly pinned crates.io dependency `jamhash = 0.1.2`.
- Hash zero is excluded consistently during catalog construction and query construction.
- Screening TSV and metadata schemas are version `1.0.0`.

See [docs/TRACE.md](docs/TRACE.md) before interpreting trace results. The
pipeline example is [examples/trace/run_trace.sh](examples/trace/run_trace.sh).

## Citation

Citation metadata are in [CITATION.cff](CITATION.cff). Until a version-specific archive DOI is available, cite the jam-rs release/tag and commit used, plus jamhash and any comparator tools used in the analysis.

## License

jam-rs is available under the [MIT license](LICENSE).
