# `jam trace` MGE evaluation

This directory records a reproducible, reference-guided evaluation of the
one-query-to-many-assembly `jam trace` workflow. The measured inputs contain
declared plasmid and phage sequence constructions on a reused chromosome
background. They test candidate retrieval, positional alignment, query
coverage, range I/O, and memory accounting. A sketch candidate or alignment is
sequence evidence only; it is not a biological confirmation of an autonomous
plasmid or phage.

## Scope and source manifest

The current MGE evaluation is Track A only. It uses the exact, versioned
manifest [mge-accession-sources-v1.tsv](manifests/mge-accession-sources-v1.tsv), whose SHA-256 is:

```text
46f48bd06c3e7b91785d9b30277366c94d1c2970efbc5feffa80aa82d2af3ba8
```

The verified source FASTA used for the recorded dataset has SHA-256
`7ee9b69d4dee59b600321e885799a0fd22ab8760f42a08ac06df858b1f44f42c`.
The manifest contains seven plasmid queries and three phage queries, plus the
three declared chromosome backgrounds:

| Stratum | Queries | Declared topology | Materialized derivatives | Unsupported planned cases |
| --- | ---: | --- | ---: | ---: |
| plasmid | 7 | circular | 21 | 14 |
| phage | 3 | linear | 7 | 8 |

The generated dataset therefore contains 10 query FASTA records and 28
controlled assembly derivatives. The assembly files contain 112 contigs and
425,138,705 bases; the query files contain 10 contigs and 663,117 bases.
Every derivative reuses the same declared `chromosome-v1` background group:
NC_000913.3, NC_000964.3, and NC_002516.2. They are not independent natural
metagenomes. These counts and source checksums are recorded in
`results/raw/mge-accession-v2/input/dataset.json`, `stats.json`,
`source_verification.tsv`, `status.tsv`, and `truth.tsv`.

The constructions are explicit sequence transformations:

- circular plasmids have forward, reverse-complement, and origin-crossing
  cases;
- linear lambda and T4 phages have forward and reverse cases;
- T7 additionally has a declared 160 bp terminal-repeat case at query
  coordinates `[0,160)`;
- cases requiring undeclared terminal-repeat or shared-region coordinates are
  recorded as `unsupported`, never as negative evidence.

The generator does not simulate reads, run an assembler, infer topology, or
claim biological origin. The `query_kind` and `topology` values are copied from
the manifest.

## Recreate the accession-backed dataset

The following command fetches exactly the accession-versioned records named by
the committed manifest and generates the verified subset used by the raw
results. The generator checks every declared length and sequence SHA-256 before
writing output.

```bash
set -euo pipefail
SOURCE_ROOT=tools/out/mge-accession-spike-20260823
mkdir -p "$SOURCE_ROOT"
awk -F '\t' 'NR > 1 {print $1}' \
  evaluation/trace-production/manifests/mge-accession-sources-v1.tsv \
  | while IFS= read -r accession; do
      curl --fail --location --silent --show-error \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=text"
    done > "$SOURCE_ROOT/sources.fasta"
python3 evaluation/trace-production/scripts/mge_accession_spike.py \
  --sources-fasta "$SOURCE_ROOT/sources.fasta" \
  --source-manifest evaluation/trace-production/manifests/mge-accession-sources-v1.tsv \
  --output "$SOURCE_ROOT/dataset-v2" \
  --verified-only
```

The output `status.tsv` is part of the truth contract. Omitted or unsupported
rows are not negatives. The generated truth uses zero-based, half-open query
coordinates and stores the source accession checksum, background group,
assembly checksum, orientation, origin-crossing flag, terminal-repeat flag,
and the declared query-coordinate union.

## Recorded environment

The MGE measurements were recorded on an AMD Ryzen 9 7945HX with 16 physical
cores and 32 logical CPUs, Linux 7.1.8-arch1-3, glibc 2.44,
`x86_64-unknown-linux-gnu`, Rust/Cargo 1.98.0, and local SSD-backed storage.
Release measurements used the repository release profile with optimization,
LTO, one code-generation unit, and `panic=abort`. The MGE environment record
is `results/raw/mge-accession-v2/environment.json`; each measured command also
records its source commit. Compilation, installation, environment creation,
and comparator setup are excluded from search timings.

## Jam trace command and candidate semantics

The recorded Jam command used one query, a local `.jam` database, a metagenome
catalog, balanced sensitivity, one CPU thread, one I/O task, no full-download
fallback, and explicit zero containment cutoffs:

```text
jam trace \
  --query query.fasta --query-kind {plasmid|phage} \
  --topology {circular|linear} \
  --database metagenomes.jam --metagenomes metagenomes.tsv \
  --sensitivity balanced --threads 1 --io-concurrency 1 \
  --memory-target 4 --min-shared 3 \
  --min-query-containment 0 --min-metagenome-containment 0 \
  --top-candidates 28 --max-alignments 256 \
  --no-full-download-fallback --silent
```

The T7 scale sweep varied the Jam index scale while retaining the same query
and evidence filters:

```text
jam sketch assemblies/ --kmer-size 31 --fscale 20  --output metagenomes-fscale20.jam
jam sketch assemblies/ --kmer-size 31 --fscale 200 --output metagenomes-fscale200.jam
```

For candidate comparator runs, BLASTn, minimap2, and LexicMap received exactly
the assembly IDs selected by the corresponding Jam JSONL candidate manifest.
The pBR322 candidate set contains 9 assemblies; T7 scale 20 contains 3 and
T7 scale 200 contains 2. Direct comparator runs use all 28 assembly files and
are reported separately. Candidate runs do not silently substitute a
different ranking, database, or assembly set.

## Pinned comparator versions

The exact comparator metadata and archive checksums are in
`environment/manifest.json`:

| Tool | Recorded version | Release identity | Archive SHA-256 |
| --- | --- | --- | --- |
| BLASTn | NCBI BLAST+ 2.17.0+ | 2.17.0 | `3888112d8207831aa47371d93583c601f058f88b5db22dc782438b039a3a411b` |
| minimap2 | 2.31-r1302 | tag `v2.31`, commit prefix `3c28777` | `300bc287f05eb890c6211fa7db043ce98320a401621fadd1cfdbeabd1a6e4ab5` |
| LexicMap | v0.9.0 | tag `v0.9.0`, commit prefix `e8dec8d` | `54dea6a35e0c1a25025ad649da47d7cc7f6190ae69e44624db9685bf3c0dd662` |

BLASTn was invoked per selected assembly with `-task megablast`, tabular
coordinates, `-evalue 1e-10`, `-max_hsps 10`, and one tool thread. Minimap2
was invoked per selected assembly with `-x asm5 --secondary=no -c --eqx` and
one tool thread. LexicMap used a separately measured index build followed by
search with the declared identity, query-coverage, and minimum-alignment
length thresholds. LexicMap index construction is never included in its
search time.

The comparator driver starts one BLASTn or minimap2 process per selected
assembly. All outputs are normalized through the shared query-coordinate
calculator before truth comparison. LexicMap's search output is normalized by
the same interval semantics. Tool rows need not be identical; they are
compared against the same truth intervals.

## Comparator results

The following values come from
`results/summary/mge-accession-v2-comparators/summary.tsv`. Wall and CPU are
seconds; RSS is peak KiB. Precision and recall are construction-truth interval
metrics, not biological specificity or sensitivity.

Jam wall time is the complete `jam trace` command, including sketch candidate
search and indexed alignment. BLASTn, minimap2, and LexicMap candidate wall
times start from the already selected assembly set. Their candidate rows are
therefore alignment-path comparisons, not end-to-end replacements for Jam's
candidate stage. LexicMap index construction is reported separately.

### pBR322 candidate set

All methods used the same 9 Jam-selected assemblies. Index timing is shown
only for LexicMap.

| Method | Base precision | Base recall | Interval precision | Interval recall | Search wall s | Search CPU s | Search RSS KiB | Index wall s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Jam 0.10.0 | 0.589803 | 1.000000 | 0.333333 | 1.000000 | 0.328650 | 0.325185 | 69,996 | — |
| BLASTn 2.17.0+ | 0.557530 | 1.000000 | 0.333333 | 1.000000 | 3.212296 | 3.014993 | 80,968 | — |
| minimap2 2.31 | 0.575203 | 1.000000 | 0.333333 | 1.000000 | 5.340259 | 5.295830 | 109,676 | — |
| LexicMap v0.9.0 | 0.557530 | 1.000000 | 0.333333 | 1.000000 | 0.849389 | 0.811883 | 819,144 | 68.623549 |

### T7 scale-20 and scale-200 candidate boundaries

These are candidate-scoped runs. Scale 20 selected 3 assemblies; scale 200
selected 2. The two scales are different Jam databases, not post-hoc timing
labels. At scale 200 the candidate stage missed the isolated 160 bp terminal-
repeat construction, so its T7 candidate and final recall were 2/3. The table
below evaluates only the assemblies actually selected at each scale.

| Scale | Method | Selected | Base precision | Base recall | Interval precision | Interval recall | Search wall s | Search RSS KiB | LexicMap index wall s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | Jam | 3 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.976459 | 335,372 | — |
| 20 | BLASTn | 3 | 0.998005 | 1.000000 | 0.833333 | 1.000000 | 1.310717 | 80,036 | — |
| 20 | minimap2 | 3 | 1.000000 | 0.998001 | 1.000000 | 0.666667 | 2.132898 | 102,772 | — |
| 20 | LexicMap | 3 | 0.998005 | 1.000000 | 0.833333 | 1.000000 | 0.661587 | 1,296,616 | 23.742383 |
| 200 | Jam | 2 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.715025 | 199,716 | — |
| 200 | BLASTn | 2 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.814500 | 80,192 | — |
| 200 | minimap2 | 2 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.307141 | 102,680 | — |
| 200 | LexicMap | 2 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 0.988508 | 808,432 | 15.523464 |

The direct T7 run compared all 28 assemblies. BLASTn recovered 3 assemblies
in 9.791788 s; minimap2 recovered 2 in 16.640964 s; LexicMap recovered 3
with interval precision 0.833333. LexicMap direct search was measured with
minimum alignment lengths 17, 100, and 160 bp; the 17 bp setting returned
records for all 28 assemblies but interval precision fell to 0.078571. The
100 and 160 bp settings returned the same 3 truth-supported assemblies. Direct
Jam timing is not claimed here because Jam is the candidate selector for this
comparison.

These small controlled results do not establish a general speed or accuracy
ordering. They are not a matched-accuracy speed claim, and they do not justify
the statement “jam-rs is faster than” any comparator outside these measured
commands, versions, hardware, and assembly sets.

## Jam strata, threads, and RSS

The MGE summary is in
`results/summary/mge-accession-v2/summary.tsv`. The baseline and thread rows
below are complete raw measurements. A row's p90 can be much larger than its
median because the assembly set includes long R100/T4-derived files.

Three changed query orders were executed. Two were strictly sequential and
form the timing summary; the third ran two query processes at a time and is
retained under `jam/t1/paired` but excluded from the sequential timing rows.

| Stratum | Sweep | Raw cases | Median wall s | Wall p90 s | Median CPU s | Median RSS KiB | Candidate recall | Final recall |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| phage | baseline, 1 thread | 6 | 2.649003 | 19.063631 | 2.617466 | 452,088 | 1.000000 | 1.000000 |
| phage | thread sweep, 8 threads | 3 | 1.079668 | 4.263013 | 2.080454 | 1,128,260 | 1.000000 | 1.000000 |
| plasmid | baseline, 1 thread | 14 | 8.872488 | 174.419355 | 8.779787 | 510,648 | 1.000000 | 1.000000 |
| plasmid | thread sweep, 4 threads | 1 | 59.186756 | 59.186756 | 220.728188 | 2,192,908 | 1.000000 | 1.000000 |
| plasmid | thread sweep, 8 threads | 7 | 2.344089 | 33.793347 | 7.959011 | 1,241,460 | 1.000000 | 1.000000 |
| plasmid | thread sweep, 16 threads | 1 | 59.183681 | 59.183681 | 218.989200 | 2,255,904 | 1.000000 | 1.000000 |

RSS is an observation, not a hard process maximum. The thread rows are
separate measured cases, not a smooth scaling curve; the wide plasmid p90 and
the one-case 4/16-thread rows are reasons to profile the production workload
before selecting concurrency.

The measured index artifacts were:

| Artifact | Size bytes | Build wall s | Build CPU s | Build RSS KiB |
| --- | ---: | ---: | ---: | ---: |
| Jam scale 20 | 253,338,645 | 4.962489 | 13.624409 | 961,868 |
| Jam scale 200 | 25,723,925 | 2.428672 | 7.760746 | 964,416 |
| JMA archive | 363,345,582 | 18.025189 | 17.882331 | 124,120 |
| JMA sidecar | 16,256,395 | — | — | — |

## Additional local resource measurements

The separate resource-scaling fixture remains useful for I/O and memory
regression checks, but it is not part of the 28-derivative MGE scale:

- 48 block-size cases completed across 16 KiB, 64 KiB, 256 KiB, and 1 MiB
  archive/cache blocks and three query shapes.
- 288 memory cases completed across 1/4/8/16 threads, fast/balanced/sensitive
  profiles, candidate limits 1/100, and four cache block sizes.
- The 100-candidate local boundary processed 83 candidates, 85 trace records,
  and 226 alignments. At 8 threads it took 14.792359 s, requested
  741,435,014 bytes, and peaked at 802,240 KiB RSS.
- The local HTTP range fixture served 6,210,430 bytes in both cold and warm
  runs, with 95 GET/range responses carrying HTTP 206 and two HEAD requests.
  Cold wall time was 3.339435 s; warm wall time was 3.313176 s. No
  full-object fallback occurred. This is an HTTP range test, not an actual S3
  cold/warm result.

## Track B and Track C gates

Track B, read-derived assemblies, requires a checksummed read manifest, an
abundance design, pinned mixing and assembly tools, and an extractor that
records the coordinates that actually assembled. The raw tree has no such
inputs. Track C, independently supported natural positives, requires licensed
assemblies and independent long-read, isolate-match, assembly-graph, or other
support manifests. The raw tree has none.

The following required evidence is explicitly unavailable, not negative:

- no actual S3 endpoint or AWS/S3 credential environment;
- no read-derived assembly track;
- no natural or independently supported positive track;
- no 1,000 independent real assembly backgrounds;
- no 100-query release-scale collection.

## Raw artifacts, summaries, and checksums

The MGE evidence tree is:

```text
results/raw/mge-accession-v2/input/
results/raw/mge-accession-v2/index/
results/raw/mge-accession-v2/jam/
results/raw/mge-accession-v2/comparators/
results/raw/mge-accession-v2/raw-manifest.json
results/summary/mge-accession-v2/
results/summary/mge-accession-v2-comparators/
```

`raw-manifest.json` records SHA-256 values for the raw MGE measurements and
outputs. The comparator summary records its own raw-file hashes.

Regenerate the MGE timing/accuracy summary and its raw manifest:

```bash
python3 evaluation/trace-production/scripts/summarize.py \
  --raw-dir evaluation/trace-production/results/raw/mge-accession-v2/jam \
  --output-dir tools/out/mge-accession-v2-summary \
  --matched-recall 1.0
```

Regenerate the comparator summary, keeping LexicMap index construction and
search timing separate:

```bash
python3 evaluation/trace-production/scripts/summarize_comparators.py \
  --raw-dir evaluation/trace-production/results/raw/mge-accession-v2/comparators \
  --output-dir tools/out/mge-accession-v2-comparator-summary
```

Each summary command writes to a new output directory. The raw files, source
manifest, truth manifest, commands, software versions, hardware record, and
checksums remain the evidence of record; summaries are derived views.

## Cache policy and limits

`process-cold` starts a fresh process and case-local cache but does not evict
the Linux kernel page cache. `warm` performs one unmeasured warm-up and then
times the recorded invocation. Jam block caches remain process-local. These
labels do not claim cold physical storage or cold object storage.

The measured MGE data demonstrate deterministic accession verification,
forward/reverse/origin coordinate handling, declared phage terminal-repeat
construction, candidate selection, comparator normalization, and bounded
resource accounting on this small reused-background set. They do not establish
biological sensitivity, specificity, absence, natural-metagenome behavior,
S3 behavior, or a general speed advantage. `jam trace` remains a candidate
finder and sequence-evidence stage; final presence calls require independent
confirmation.
