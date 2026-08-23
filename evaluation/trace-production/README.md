# `jam trace` production evaluation

This directory contains the reproducible evaluation harness for the
one-plasmid-to-many-metagenomes `jam trace` workflow. The benchmark measures
candidate search, sequence evidence, coverage, resource I/O, and memory. A
sketch candidate or an alignment is not, by itself, evidence that an autonomous
plasmid is present; confirmation requires independent mapping, assembly-graph,
long-read, isolate, or contextual evidence.

The checked-in manifests and scripts are intentionally separate from the
input sequences. Large source FASTA files and generated assemblies are
referenced by checksum in the raw records and are not repository data.

## Environment and measured scope

The recorded host was an AMD Ryzen 9 7945HX with 16 physical cores and 32
logical CPUs, Linux 7.1.8-arch1-3, glibc 2.44, `x86_64-unknown-linux-gnu`,
Rust/Cargo 1.98.0, and local Btrfs on SSD-backed NVMe storage. Release runs
used `cargo build --release --locked` with the repository's release profile
(LTO, one code-generation unit, `panic=abort`). See
`results/raw/environment.json` and each raw measurement for the exact source
identity and storage state.

Track A used 148 controlled derivatives containing 158,749 contigs and
2,239,344,218 bases. These derivatives reuse one checksum-verified,
accession-backed chromosome background; they are not 148 independent
metagenomes. The matrix requests query lengths from 2 kb through 200 kb,
identities from 100% through 85%, 5% through 100% query-coordinate coverage,
one, two, five, and twenty fragments, both orientations, and the declared
short/long insertion and deletion cases. Source-length combinations that
cannot be materialized remain recorded as unsupported rather than being
silently treated as negatives.

The accession-backed source set contains the catalog queries J01749.1 (pBR322),
AP000342.1 (R100), and NC_008272.1 (pKJK5), near-known records M77789.2 and
KY749247.1, catalog distractors AP001918.1 and NC_001735.1, and chromosome
sources NC_000913.3, NC_000964.3, and NC_002516.2. Every source record used in
the generated dataset is verified against
`manifests/accession-sources-v1.tsv`.

## Recreate the accession-backed source dataset

The source manifest is the authority for accessions, lengths, roles, and
checksums. The following fetches the accession-versioned records from NCBI
nuccore and then lets `accession_spike.py` verify them before generating the
four controlled assemblies:

```bash
set -euo pipefail
SOURCE_ROOT=tools/out/trace-accession-production-v1
mkdir -p "$SOURCE_ROOT"
awk -F '\t' 'NR > 1 {print $1}' \
  evaluation/trace-production/manifests/accession-sources-v1.tsv \
  | while IFS= read -r accession; do
      curl --fail --location --silent --show-error \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text"
    done > "$SOURCE_ROOT/sources.fasta"
python3 evaluation/trace-production/scripts/accession_spike.py \
  --sources-fasta "$SOURCE_ROOT/sources.fasta" \
  --source-manifest evaluation/trace-production/manifests/accession-sources-v1.tsv \
  --output "$SOURCE_ROOT/dataset"
```

The generator is deterministic for a given verified source FASTA. It fragments
the three chromosome records using the declared accession-derived seed and
adds exact, reverse-complement, partial, split-overlap, and near-known spike
records. It does not simulate reads or run an assembler. The generated
`dataset.json`, `source_verification.tsv`, and `truth.tsv` are part of the
provenance chain.

## Track A controlled truth

The controlled generator consumes a checksum-verified source root and the
matrix in `manifests/controlled-matrix-v1.json`. A full matrix run is:

```bash
SOURCE_ROOT=tools/out/spike-in/accession-v1-release-d4e0df2-20260823/dataset
SOURCE_MANIFEST=tools/out/trace-production-controlled-source/source.json
python3 evaluation/trace-production/scripts/controlled.py \
  --source-root "$SOURCE_ROOT" \
  --source-manifest "$PWD/$SOURCE_MANIFEST" \
  --matrix evaluation/trace-production/manifests/controlled-matrix-v1.json \
  --output tools/out/trace-production-controlled-v1
python3 evaluation/trace-production/scripts/controlled.py \
  --validate-output tools/out/trace-production-controlled-v1
```

For a bounded run, use explicit axes such as:

```bash
python3 evaluation/trace-production/scripts/controlled.py \
  --source-root tools/out/spike-in/accession-v1-release-d4e0df2-20260823/dataset \
  --source-manifest "$PWD/tools/out/trace-production-controlled-source/source.json" \
  --matrix evaluation/trace-production/manifests/controlled-matrix-v1.json \
  --output tools/out/trace-production-controlled-smoke \
  --query-length 2000 --identity 100 --coverage 100 --fragments 1 \
  --orientation forward --indel none --label standard
```

The generator writes `manifest.json`, `status.tsv`, `truth.tsv`,
`truth.jsonl`, and `checksums.tsv`. Truth coordinates are zero-based,
half-open query coordinates; coverage is a union, so overlapping fragments
cannot increase the true covered length twice. The materialized release
measurement set consists of the 108 identity/coverage cases, 40 fragment and
indel cases, and the separately recorded origin-crossing case. The 148-case
scale above is the count used for the production-evaluation input summary.

## Track B and Track C prerequisites

Track B (read-derived assemblies) requires a licensed, checksum-pinned read
manifest, abundance design, pinned read mixer, pinned assembler, and an
extractor that records the coordinates that actually assembled. It must score
assembled sequence recovery, not sequence that was present only in the input
reads. Run the fail-closed validator with an explicit manifest:

```bash
python3 evaluation/trace-production/scripts/read_derived.py \
  --manifest evaluation/trace-production/manifests/read_derived.example.json \
  --read-manifest evaluation/trace-production/manifests/read_derived.reads.example.json \
  --output-dir tools/out/trace-production-read-derived-v1
```

Track C requires licensed natural assemblies, query sequences, and independent
support records with checksums. It must remain separate from controlled truth:

```bash
python3 evaluation/trace-production/scripts/natural.py \
  --manifest evaluation/trace-production/manifests/natural.example.json \
  --support-manifest evaluation/trace-production/manifests/natural.support.example.json \
  --output-dir tools/out/trace-production-natural-v1
```

Neither track was available for the recorded run. There was no actual S3
endpoint or credentialed object store, no checksummed read set and pinned
assembler, and no independently supported natural-positive manifest. The
requested scale of 1,000 independent real assemblies and 100 query plasmids
was therefore not run. These are unavailable measurements, not negative
results; see `results/raw/unavailable.json`.

## Pinned comparison environment

The comparator image is defined by `environment/Dockerfile` and
`environment/versions.env`. It pins:

| Tool | Version / release | Source archive SHA-256 | Additional checksum |
| --- | --- | --- | --- |
| BLASTn | NCBI BLAST+ 2.17.0 | `3888112d8207831aa47371d93583c601f058f88b5db22dc782438b039a3a411b` | official MD5 `bdec166721de3b55f90a3badc83538e8` |
| minimap2 | 2.31, tag `v2.31`, commit prefix `3c28777` | `300bc287f05eb890c6211fa7db043ce98320a401621fadd1cfdbeabd1a6e4ab5` | — |

Build it outside the timed search commands:

```bash
evaluation/trace-production/environment/build.sh
```

The Dockerfile verifies both archive checksums, installs the comparator
binaries in a minimal image, and runs the final image as an unprivileged
`comparator` user. Installation, image creation, and comparator setup are not
included in search timings.

## Execution semantics

The Jam controlled comparator uses the release binary and one `jam trace`
invocation over the selected assembly directory/JAM index. For the fragment
comparison its exact command parameters were:

```text
jam --threads 8 --memory-target 4 --silent trace
  --plasmid q10000.fasta --database fragments.jam
  --catalog catalog.fragments.tsv --sensitivity {fast|balanced|sensitive}
  --min-shared 3 --min-plasmid-containment 0
  --min-metagenome-containment 0 --top-candidates 250
  --max-alignments 256 --cache-block-bytes 1048576
```

The three profiles are separate measured runs. `--memory-target` is a
concurrency/cache target, not a hard process-RSS ceiling.

BLASTn and minimap2 use the same query, the same 40 assembly files, the same
truth coordinates, and the same interval-normalization code. BLASTn runs
`-task megablast -evalue 1e-10 -max_hsps 10 -num_threads 8` and emits tabular
query/subject coordinates. minimap2 runs `-x asm5 --secondary=no -t 8 -c
--eqx` and emits PAF. Both comparator modes start one tool process per
assembly, so their timings include the driver and process startup for 40
assemblies. They exclude installation and setup. The normalized outputs do
not require identical rows between tools; they compare the resulting query
coordinate union against the same truth manifest.

## Controlled fragment results

The following table is one measured 10 kb, 97% identity, 80% coverage,
one-to-twenty-fragment set over 40 assemblies. `base recall` is the weighted
truth-base intersection divided by 320,000 truth bases. `interval recall` is
the mean per-assembly truth-interval recall. Precision was 1.0 for every row.
Wall/CPU are seconds and RSS is the measured peak child RSS in KiB. Values are
from the raw measurement and normalized files; rerun the summary generator to
recalculate derived values.

| Method | Assemblies with records | Base recall | Interval recall | Wall s | CPU s | Peak RSS KiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Jam fast | 20/40 | 0.191121875 | 0.400 | 0.176063 | 0.510490 | 173160 |
| Jam balanced | 38/40 | 0.592781250 | 0.8375 | 0.533666 | 1.851799 | 274356 |
| Jam sensitive | 40/40 | 0.791131250 | 0.9500 | 2.058892 | 3.984025 | 243172 |
| BLASTn 2.17.0 | 40/40 | 0.982890625 | 1.0000 | 14.844675 | 14.714881 | 74848 |
| minimap2 2.31 | 40/40 | 0.869306250 | 0.9250 | 19.264993 | 27.389398 | 90344 |

These direct comparator runs show why no matched-accuracy speed claim is made:
the Jam sensitive run does not reach the BLASTn or minimap2 recall on this
controlled set. The table is a bounded task result, not a general statement
that Jam is faster or more accurate than either comparator.
The recorded `min-shared=1` sweep produced the same candidate and coverage
counts on this track; both settings remain in the raw results.

The accession-backed exact/partial/split outcomes are also recorded directly
in the per-query `summary.tsv` files:

| Query | Query length | Balanced t8 wall s | Mean coverage | Exact candidate/alignment | Partial candidate/alignment | Split-overlap candidate/alignment | Short-fragment candidate/alignment |
| --- | ---: | ---: | ---: | --- | --- | --- | --- |
| J01749.1 (pBR322) | 4,361 | 0.073643 | 0.820454 | 1/1 | 0/0 | 0/0 | 1/1 |
| AP000342.1 (R100) | 94,281 | 1.803527 | 0.865484 | 1/1 | 0/0 | 1/1 | 1/1 |
| NC_008272.1 (pKJK5) | 54,383 | 0.775608 | 0.440052 | 1/1 | 1/1 | 0/0 | 1/1 |

The exact and short-fragment rows are candidate/alignment recovery counts,
not final biological calls. Near-known source truth is retained for later
retrieval evaluation, but it is not converted into an independent presence
claim.

## Range, block, and memory measurements

The block study has 48 complete cases: four archive block sizes (16 KiB,
64 KiB, 256 KiB, 1 MiB) crossed with four cache block sizes and three query
shapes (small, large, repeat-rich). Raw rows are in
`results/raw/resource-scaling-v1/block/summary.tsv` and include archive build,
wall/CPU/RSS, request counts, requested/returned/decoded bytes, cache events,
seed buckets, and sequence blocks.

The memory matrix has 288 complete cases: three query shapes, three
sensitivity profiles, four thread counts (1, 4, 8, 16), two candidate limits
(1, 100), and four cache block sizes. RSS is reported as an observation, not as
a hard cap. The 100-candidate large-catalog boundary produced 83 selected
candidate metagenomes, 85 trace records, and 226 alignments. At eight threads
it took 14.792359 seconds, used 741,435,014 requested bytes and peaked at
802,240 KiB RSS. The 1/4/8/16-thread raw files show why concurrency should be
chosen from a production workload rather than inferred from the tiny fixture.

The local HTTP range fixture used the large 50 kb query with sensitivity
`sensitive`, eight threads, a 64 KiB cache block, three top candidates, and
`--no-full-download-fallback`. Both cold and warm runs served 6,210,430 bytes
with 95 GET/range responses carrying HTTP 206 and two HEAD requests. Cold wall
time was 3.339435 s and warm wall time was 3.313176 s; no full-object fallback
occurred. This is an HTTP 206 test, not an actual S3 cold/warm result.

## Raw artifacts and regeneration

Raw records are grouped by evidence source:

```text
results/raw/environment.json
results/raw/unavailable.json
results/raw/accession-spike-v1/{source,<accession-version>}/
results/raw/controlled-v1/{input,measurements,normalized,truth}/
results/raw/resource-scaling-v1/{block,memory,http,large-catalog}/
```

Each measured case keeps its command, software identity, input/truth
references, timings, RSS, candidate/alignment counts, and resource metrics.
`controlled.py` writes source and derivative checksums; the run harness writes
raw-case checksums. Do not edit a raw record by hand. To regenerate summary
TSV/JSON and the output checksum manifest from the raw directory, use the
repository's summary tool:

```bash
python3 evaluation/trace-production/scripts/summarize_evidence.py \
  --raw-root evaluation/trace-production/results/raw \
  --output-dir tools/out/trace-production-summary-check
```

The command refuses an existing output directory and emits `summary.json`,
`summary.tsv`, and `raw-checksums.json`. The raw inputs, truth manifests,
commands, and source manifests remain the evidence of record; summaries are
derived views.

## Cache policy and limitations

`process-cold` starts a fresh process and case-local cache but does not evict
the Linux kernel page cache. `warm` performs one unmeasured warm-up and then
times the recorded run; Jam block caches remain process-local. These labels do
not claim cold physical storage or cold object storage.

The current evidence is sufficient to exercise deterministic truth generation,
exact/reverse/partial/split retrieval, interval normalization, comparator
semantics, range accounting, and bounded-memory behavior on the recorded
inputs. It is not sufficient to establish biological sensitivity, specificity,
absence, natural-metagenome performance, S3 behavior, or a production-scale
speed advantage. In particular, the controlled assemblies reuse one
chromosome background, the read-derived and natural-positive tracks are
unavailable, and the 1,000-assembly/100-query release-scale run was not
performed.
