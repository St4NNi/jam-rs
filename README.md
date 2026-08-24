[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)

# jam-rs

`jam-rs` provides reference-guided sketch screening and query-centered positional follow-up for metagenomic assemblies. `jam screen` reports contig-level and assembly-level evidence against a fixed plasmid catalog. A local Jam Index searches one plasmid, phage, other mobile element, synthetic construct, or unknown nucleotide reference across a collection without a positional seed archive. Optional target-aware sampling may improve retrieval of represented catalog targets, but final presence calls require independent confirmation.

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
  --memory-target 4
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
  --memory-target 4
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

## Local collection search with Jam Index

A Jam Index is one logical local dataset with a manifest and independently
searchable, append-only parts. Each part contains a pure version-3 `.jam`
screening shard, position-free contig signatures, a compact contig directory,
and the complete two-bit contig sequences.

The source catalog binds each metagenome ID to one local FASTA/FASTQ and its
checksum:

```text
metagenome_id	resource_uri	sha256
sample_001	/data/assemblies/sample_001.fasta	<sha256>
sample_002	/data/assemblies/sample_002.fasta	<sha256>
```

Build the standard index and split parts by both bases and estimated signature
count:

```bash
jam index build \
  --metagenomes sources.tsv \
  --output metagenomes-index \
  --max-part-bases 100000000 \
  --max-part-signatures 250000 \
  --parallel-parts 4
```

New metagenomes are added as new independent parts; existing parts are not
rewritten:

```bash
jam index append \
  --metagenomes new-sources.tsv \
  --output metagenomes-index \
  --parallel-parts 4
```

Search every part in parallel and align only globally selected metagenomes and
contigs:

```bash
jam trace \
  --query element.fasta \
  --query-kind plasmid \
  --topology auto \
  --index metagenomes-index \
  --output element.trace.jsonl.zst \
  --sensitivity sensitive \
  --min-shared 1 \
  --threads 4 \
  --memory-target 2
```

The standard policy uses canonical k=21 `jamhash_u64_v1` signatures with a
length-dependent budget of 16 to 256 hashes per contig, one hash per 1,024
bases within those bounds, a 512-hash whole-metagenome sketch, and 256-base
query windows. Exact shared hashes route selected candidates to bounded contig
sets without storing nucleotide positions. Selected complete contigs generate
dense k=31 and k=21 positions in memory, verify packed k-mers, and enter the
existing mixed chaining, corridor alignment, and fragment-mosaic path. An
exact complete-contig match emits a direct `=` alignment before seed chaining.

Jam Index operation is local-only. It does not open JMA objects, store
collection-wide positional postings, or change ordinary `.jam` format or
commands.

## Existing `.jam` and JMA trace path

`jam trace` is an optional follow-up for one query sequence across many candidate
metagenomes. It uses an existing `.jam` metagenome index, then retrieves
positional evidence from JMA v1 archives or raw assembly resources:

```bash
jam archive \
  --input assembly.fasta \
  --output assembly.jma \
  --primary-scale 200 \
  --rescue-scale 500

jam trace \
  --query element.fasta \
  --query-kind plasmid \
  --topology auto \
  --database metagenomes.jam \
  --metagenomes metagenomes.tsv \
  --output element.trace.jsonl.zst \
  --sensitivity balanced \
  --min-shared 3 \
  --threads 8 \
  --io-concurrency 8 \
  --memory-target 4
```

`--query` accepts exactly one FASTA record. `--plasmid` remains a deprecated
one-release alias that implies `--query-kind plasmid`. Linear, circular, auto,
and unknown topology settings control coordinate handling. A wrapped coordinate
model does not prove that the biological molecule is circular.

The archive command writes `assembly.jma` and its checksum-bound range index
`assembly.jma.idx.json`. Production catalog rows identify both objects; strict
trace runs reject a JMA row without its range index.

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
  --query phage.fasta \
  --query-kind phage \
  --topology linear \
  --database https://objects.example.org/metagenomes.jam \
  --cache-dir cache/jam \
  --metagenomes metagenomes.tsv \
  --output phage.trace.jsonl
```

Alignment and nonredundant query-coordinate fragment mosaics are `supported`
evidence, not confirmation of autonomous element presence or physical linkage
between separate contigs. See
[docs/TRACE.md](docs/TRACE.md), [docs/ALGORITHM.md](docs/ALGORITHM.md),
[docs/REMOTE_RESOURCES.md](docs/REMOTE_RESOURCES.md), and
[examples/trace/run_trace.sh](examples/trace/run_trace.sh).

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

The bounded release evidence and exact limitations are under
[evaluation/trace-production/](evaluation/trace-production/). On the selected
anonymous 40-case matrix, the standard Jam Index recovered 40/40 candidates
with base precision 1.0, base recall 0.981834375, and interval precision and
recall 1.0. Three paired process-cold four-thread runs had a 4.364 s median and
313,300 KiB median peak RSS. The 605,224,880-base collection occupied
183,144,518 bytes, or 0.302606 index bytes per source base. The exact standalone
160-base trace, reverse-complement, origin-crossing, separate-contig, overlap,
repeat-shared, integrated, and unrelated controls were rerun.

These are measurements of one local task, not a general speed or accuracy
claim. The anonymous edited-fragment matrix missed both 80-base/90%-identity
orientations because neither retained an exact k=21 witness. Long related-query
cases remain CPU-bound, and read-derived assemblies, independently supported
natural positives, and the 1,000-assembly/100-query release scale remain
unmeasured.

## Compatibility and provenance

- `.jam` binary format version remains 3. Existing version-3 files remain readable; the repository tests a fixture created by jam-rs 0.9.11.
- Jam Index format 1 is a local manifest-and-parts dataset. Its screening shards are ordinary version-3 `.jam` files; adding a part does not rewrite earlier parts.
- New builds add `<database>.jam.json`; this sidecar does not alter the binary database.
- The released hash identity is `jamhash_u64_v1`, provided by the exactly pinned crates.io dependency `jamhash = 0.1.2`.
- Hash zero is excluded consistently during catalog construction and query construction.
- Screening TSV and metadata schemas are version `1.0.0`.
- Trace JSONL uses schema `2.0.0`; schema `1.x` remains a documented older
  output contract rather than being silently reinterpreted.

See [docs/TRACE.md](docs/TRACE.md) before interpreting trace results. The
pipeline example is [examples/trace/run_trace.sh](examples/trace/run_trace.sh).

## Citation

Citation metadata are in [CITATION.cff](CITATION.cff). Until a version-specific archive DOI is available, cite the jam-rs release/tag and commit used, plus jamhash and any comparator tools used in the analysis.

## License

jam-rs is available under the [MIT license](LICENSE).
