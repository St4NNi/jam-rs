<!-- What jam does and how to use its commands. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam

jam (just another minimizer) searches for query sequences (plasmids, genes, contigs; 1 kb to a
few Mb) in very large collections of assemblies. It computes full base-level alignments with
CIGAR strings and decides for each assembly whether it contains the query.

## Install

Pre-built binaries will be linked here once there is a release.

Until then, build jam from source. jam is written in Rust (edition 2024) and requires a recent
stable toolchain. Install Rust with [rustup](https://rustup.rs), then run:

    cargo build --release --locked

This puts the binary in `target/release/jam`. Run `jam --help` or `jam COMMAND --help` to see
all options.

The block cache uses a patched copy of the `foyer` crate. `Cargo.toml` fetches this crate from
a git repository (`[patch.crates-io]`), so the first build requires network access to GitHub.
To build offline, run `cargo vendor` on a machine with network access, then build with
`--offline`.

The `Dockerfile` builds the same binary into a small container image.

## Quick start

Build an index from your assembly files (plain or compressed FASTA). Each file becomes one
entry in the index:

    jam index build -o assemblies -t 8 assemblies/*.fa.gz

This writes the index file `assemblies.jam` and the sequence file `assemblies.jamseq`.

Search for the sequences of a query FASTA file in the index:

    jam search assemblies.jam plasmid.fa > hits.jsonl

Each query record is searched separately. For each query, jam writes one JSON line per assembly
that contains the query, followed by one summary line.

## How it works

Each assembly in the index is stored as an *entry*. With `--singleton`, each FASTA record
becomes a separate entry instead. A search runs in three stages and answers one question per
stage:

1. **Sketch**: which entries may contain the query? jam looks up the query's seeds (minimizers)
   in the index. Each entry gets a containment estimate: the share of the query that the entry
   appears to contain, estimated from shared seeds alone.
2. **Locate**: where do the shared seeds lie within those entries? This stage finds candidate
   alignment windows on the contigs of each entry.
3. **Align**: how does the query align within those windows? jam computes base-level alignments
   in each window and reports the entry if the alignments pass the decision rule.

Each stage reads more data than the stage before it, and processes only the entries kept by the
previous stage. You can stop a search after any stage with `--level`:

- `--level sketch` runs a containment search, similar to `sourmash search --containment`. It is
  fast and reads only the index file, but mutations lower the containment estimate and no
  alignments are computed.
- `--level locate` reports candidate regions (the windows on each contig where the query may
  align) without computing base-level alignments.
- `--level align` (the default) computes base-level alignments and applies the decision rule.

Between stages, jam drops entries that can no longer reach the required query coverage. This
pruning has a small, modeled probability of missing a true match; see
[The miss chance](#the-miss-chance).

## Commands

| Command | Purpose |
| --- | --- |
| `jam index build` | Build an index file and its sequence file from FASTA files. |
| `jam index reshard` | Merge and split index files into a directory of new index files. |
| `jam index warm` | Pin parts of an index in a local cache directory. |
| `jam index copy` | Copy an index between local directories and S3. |
| `jam stats` | Show the parameters and contents of an index file. |
| `jam search` | Search for query sequences in an index and align them. |
| `jam serve` | Serve an index over the network for `jam search --server`. |
| `jam screen` | Stream large assemblies against an index of small sequences. |

### jam index build

    jam index build -o NAME [OPTIONS] FASTA...

Builds the index file `NAME.jam` and the sequence file `NAME.jamseq`. Neither file may already
exist.

Normal options:

- `-o, --output NAME`: base name of the index to create.
- `-l, --list FILE`: text file with one FASTA path per line, used in addition to any paths on
  the command line.
- `--singleton`: treat each FASTA record as a separate entry, instead of treating each file as
  one entry.
- `-t, --threads N` (default: all cores): number of threads to use.

Advanced options:

- `-k, --k K` (default 25): k-mer length in bases (1 to 32).
- `-w, --window W` (default 100): number of k-mers per minimizer window.
- `--sketch-only`: write only the index file containing the sketch data. This index is much
  smaller, but it supports only `--level sketch` searches.
- `--chunks N` (default 1): split the sequence data into `NAME.1.jamseq` through
  `NAME.N.jamseq` of roughly equal size (for example, to stay within object size limits or to
  parallelize file transfers). Each entry is stored completely within one chunk.
- `--temp-dir DIR` (default: the output directory): directory for temporary files.

### jam index reshard

    jam index reshard -o DIR (--shards N | --max-seed-locations N) [OPTIONS] INDEX...

Merges and splits existing indexes into a new directory `DIR` containing shard files named
`shard-0000.jam` and so on. Output shards are balanced so that each holds a similar number of
seed locations. Inputs can be individual index files or directories produced by an earlier
reshard. Entries stay in the same order as the inputs.

- `-o, --output DIR`: the output directory to create. This directory must not already exist.
- `--shards N`: number of output shards to create.
- `--max-seed-locations N`: maximum number of seed locations per output shard, which limits the
  size of each index file. Specify either this option or `--shards`.
- `--copy-sequences`: write a new sequence file for each output shard. Without this option, the
  new index files point to the original sequence files without copying them.
- `-t, --threads N` (default: all cores): number of threads to use.
- `--temp-dir DIR` (default: the output directory): directory for temporary files.

Example:

    jam index reshard -o all --shards 4 part1.jam part2.jam

To search across multiple shards, pass all of their index files to `jam search` before the
query file.

### jam index warm

    jam index warm INDEX --cache DIR --memory SIZE --disk SIZE [--pin PARTS]

Copies byte ranges of an S3 index into a local cache directory so that later searches read them
from RAM or local disk. See [S3 and local caches](#s3-and-local-caches).

- `--cache DIR`: local cache directory to populate. `jam search` and `jam serve` read pinned
  data from this directory.
- `--memory SIZE`: RAM budget for pinned data and the block cache (for example, `16GiB`).
- `--disk SIZE`: disk budget for pinned data and the block cache (for example, `200GiB`).
- `--pin PARTS`: index sections to pin instead of the default plan. Give a comma-separated list
  of `sketch`, `locators`, `packed`, or `packed=NAME,NAME` (to pin sequences for specific
  entries). Because every item after `packed=` is read as an entry name, that item must come
  last in the list. You can pass this flag multiple times.

### jam index copy

    jam index copy SOURCE TARGET [--cache DIR]

Copies an index file and all of its sequence files. Both source and target may be a local path
or `s3://bucket/prefix`. `SOURCE` must be the index file `NAME.jam`. `TARGET` may be a
directory, `s3://bucket/prefix`, or the full path to the new index file `NAME.jam`.

- `--cache DIR`: local cache directory in which to record the S3 ETag of each uploaded file.
  Searches that use this cache directory will then verify from their first read that the files
  on S3 match the uploaded version.

### jam stats

    jam stats INDEX [--json]

Shows the parameters and summary counts of an index file: entries, contigs, bases, seeds, seed
locations, and the size in bytes of each index section.

- `--json`: print JSON instead of a summary table.

### jam search

    jam search [OPTIONS] INDEX... QUERY

Searches for each sequence of the FASTA file `QUERY` across all given index files. The query
name is the first word of its FASTA header. When using `--server`, pass only the query file.

Normal options:

- `-o, --output FILE`: write results to this file instead of standard output.
- `--format json|tsv|paf` (default `json`): output format; see [Output](#output). The `tsv` and
  `paf` formats require `--level align`.
- `-t, --threads N` (default: all cores): number of threads to use.
- `-c, --coverage C` (default 0.5): minimum share of the query that the alignments must cover.
  Give a fraction from 0 to 1 (such as `0.8`) or a base count (such as `800bp`).
- `-i, --identity I` (default 0.95): minimum sequence identity across the covered query bases.
- `--circular`: treat query sequences as circular. An alignment fragment that wraps across the
  end of the query is reported with an end coordinate beyond the query length. By default,
  queries are treated as linear.
- `-v, --verbose`: write one progress line per query to standard error. Pass twice (`-vv`) to
  also write per-assembly debug lines.

Advanced options:

- `--stretch-identity I` (default: the value of `--identity`): minimum identity required for a
  local stretch of an alignment to count toward query coverage.
- `--stretch-length N` (default 40): minimum length (in alignment columns) for a local stretch
  of an alignment to count toward query coverage.
- `--contig-coverage C` (default 0, off): count a contig only if the alignments cover at least
  this share of the contig itself.
- `--contig-identity I` (default 0, off): count a contig only if the alignments on that contig
  reach at least this identity.
- `--level sketch|locate|align` (default `align`): stop the search after this stage; see
  [How it works](#how-it-works).
- `--max-miss P` (default 0.0001): maximum allowed probability of missing a qualifying
  assembly; see [The miss chance](#the-miss-chance).
- `--seed-mismatches auto|0|1|2` (default `auto`): number of substitutions allowed when
  matching a query k-mer to an indexed seed. Setting `1` also looks up every k-mer one
  substitution away, and `2` looks up every k-mer two substitutions away. Allowing mismatches
  finds assemblies at lower identity, at the cost of extra CPU time. `auto` chooses a setting
  for each query; see [The miss chance](#the-miss-chance).
- `--deadline SECONDS`: stop searching a query after this many seconds. The summary line then
  reports the status as `stopped`.
- `--server NODE_ID@ADDRESS`: run the search on a `jam serve` server instead of local index
  files.
- `--cache DIR`, `--memory SIZE`, `--disk SIZE`: local cache directory and memory or disk
  budgets for index files on S3; see [S3 and local caches](#s3-and-local-caches). `--memory`
  defaults to the budget recorded in the cache directory, or 1GiB if there is none. `--disk`
  defaults to the recorded budget, or 0 if there is none.

Example:

    jam search assemblies.jam plasmid.fa -c 0.8 -i 0.9 --format paf -o hits.paf

### jam serve

    jam serve INDEX --key FILE [OPTIONS]

Keeps an index open and answers search requests from other machines; see [Server](#server).

- `--key FILE`: path to the secret key file. Its public key is the server's node ID. If the
  file does not exist, it is created on first start.
- `--bind ADDRESS` (default `0.0.0.0:4100`): UDP address to listen on. Pass this flag again to
  add an IPv6 address.
- `--peer NODE_ID@ADDRESS`: address of a peer server. Every search also runs on each peer. Pass
  this flag multiple times to add more peers.
- `--searches N` (default 4): maximum number of searches to run at the same time.
- `--queue N` (default 32): maximum number of searches waiting for a slot. Additional requests
  are rejected as busy.
- `-t, --threads N` (default: all cores): number of threads used per search.
- `--cache DIR`, `--memory SIZE`, `--disk SIZE`: same as for `jam search`.

### jam screen

    jam screen INDEX --assemblies LIST [OPTIONS]

Works in the opposite direction: the index holds many small sequences (such as plasmids or
genes), and large assemblies stream past it one by one. See
[Screening assemblies](#screening-assemblies).

Normal options:

- `--assemblies LIST`: text file with one assembly path per line (plain or gzipped FASTA). Each
  assembly is named after its file name.
- `-c, --coverage C` and `-i, --identity I`: minimum coverage and identity required for each
  indexed small sequence, as in `jam search`.
- `-t, --threads N` (default: all cores): number of worker threads. Each thread screens one
  assembly at a time.
- `-o, --output FILE`: write results to this file instead of standard output.

Advanced options:

- `--stretch-identity`, `--stretch-length`, `--contig-coverage`, `--contig-identity`: same as
  for `jam search`.
- `--closest-margin M` (0 to 1): for each contig, score every matching entry as the aligned
  share of the contig multiplied by the alignment identity. A contig then counts only toward
  entries whose score is within margin M of the highest-scoring entry on that contig. Without
  this option, a contig counts toward every entry it qualifies for.
- `--batch-bases N` (default 500000): number of contig bases searched together in one batch.
- `--contig-dir DIR`: write contigs that share enough with an entry into this directory, as one
  gzipped FASTA file per assembly.
- `--survey-only`: stop after the survey stage. No alignments are computed and no entries are
  reported; only the files in `--contig-dir` are written.

## Output

`jam search` writes only search results to standard output (or to `-o FILE`). Progress and
debug lines go to standard error.

### JSON lines (the default)

jam writes one JSON line per reported assembly, according to the search level. Every result
line has these keys:

- `q`: the query name.
- `a`: the assembly (entry) name.
- `s`: the stage that the assembly reached (1 at `--level sketch`, 2 at `--level locate`, or 3
  when aligned and confirmed).

At `--level sketch`, jam writes one line per candidate assembly:

    {"q":"plasmid1","a":"SAMN0001","s":1,"c":0.83}

- `c`: the containment estimate (0 to 1), which is the share of the query that the assembly
  appears to contain based on shared seeds alone. Mutations lower this estimate.

At `--level locate`, jam writes one line per located candidate:

    {"q":"plasmid1","a":"SAMN0001","s":2,"c":0.81,"w":[["contig_7",120000,148000]]}

- `c`: the containment estimate, recomputed from the seeds whose positions were located.
- `w`: the candidate alignment windows, each given as `[contig, start, end]` on the contig.

At `--level align`, jam writes one line per confirmed assembly:

    {"q":"plasmid1","a":"SAMN0001","s":3,"cov":0.982,"id":0.9953,"bases":8243,
     "f":[["contig_7",0,8243,123456,131710,"+","4000=1X4242="]]}

- `cov`: query coverage, the share of query bases covered by the alignments.
- `id`: sequence identity across the covered query bases.
- `bases`: number of covered query bases.
- `f`: aligned fragments, ordered from best to worst. Each fragment is `[contig, query start, query end, contig start,
  contig end, strand, CIGAR]`. Coordinates are 0-based and half-open (the end coordinate is one
  past the last aligned base). On a circular query, a fragment that wraps across the end of the
  query has an end coordinate past the query length. The strand is `+` or `-`. The CIGAR string
  uses `=`, `X`, `I`, and `D` as in SAM, with the contig as the reference and operations
  ordered along the contig.
- `alt`: alternative fragments in the same format (included only when present). These fragments
  cover no new query bases (for example, additional copies of a repeated region).

Share and identity values are rounded to 4 decimal places.

After all result lines for a query, jam writes one summary line:

    {"q":"plasmid1","status":"complete","n":[120,14,9,3],"miss":0.0001,"t":[0.21,0.35,1.8]}

- `q`: the query name.
- `status`: how the search ended; see the table below.
- `reason`: why the search did not finish normally (included only when the status is not
  `complete`).
- `error`: the error message (included only when the status is `failed`).
- `n`: assembly counts at each step, as `[seeded, possible, located, confirmed]`. Seeded
  assemblies share at least one seed with the query. Possible assemblies passed stage 1,
  located assemblies passed stage 2, and confirmed assemblies passed stage 3 and were reported.
- `miss`: modeled probability that the search missed a qualifying assembly; see
  [The miss chance](#the-miss-chance).
- `prune`: included only when `false`. Indicates that pruning was turned off for this query
  because seed matching alone could not stay within `--max-miss`.
- `t`: time in seconds spent in `[sketch, locate, align]`. The sketch time is wall-clock time,
  including query preparation. The locate and align times are summed across threads.

When searching multiple index files, the summary combines all of them: the counts are summed,
`miss` is the largest value across indexes, and the status is the worst outcome.

The status values:

| `status` | `reason` | Meaning |
| --- | --- | --- |
| `complete` | | Every candidate assembly was evaluated. |
| `incomplete` | `resource_limit` | An assembly required more memory or work than a limit allows (for example, when a periodic query matches many periodic regions) and could not be evaluated. |
| `stopped` | `deadline` | The search exceeded `--deadline`. |
| `stopped` | `asked` | The search was asked to stop (for example, when the program reading standard output closed the pipe). |
| `failed` | `error` | A search stage or the server failed. `error` gives the reason. |

Whenever a query does not end with `complete`, jam writes all remaining output and then exits
with an error. This also applies to `--format tsv` and `paf`, even though those formats have no
summary line.

### TSV

`--format tsv` writes a header line followed by one tab-separated line per confirmed assembly:

| Column | Meaning |
| --- | --- |
| `query` | The query name. |
| `assembly` | The assembly (entry) name. |
| `coverage` | The query coverage, as in `cov` above. |
| `identity` | The sequence identity, as in `id` above. |
| `bases` | The number of covered query bases. |
| `fragments` | The number of aligned fragments. |

### PAF

`--format paf` writes one line for each fragment of a confirmed assembly. The query is the PAF
query, and the contig is the PAF target. Columns 1 to 12 follow the standard PAF format: query
name, query length, query start, query end, strand, contig name, contig length, contig start,
contig end, matching bases, alignment columns, and mapping quality (always 255). These columns
are followed by `tp:A:P` (`tp:A:S` for an alternative fragment), `NM:i:` with the edit count,
and `cg:Z:` with the CIGAR string described above.

### Progress and debug lines

`-v` writes one progress line per query to standard error when the query's search finishes, for
example:

    plasmid1: complete; 120 seeded, 14 possible, 9 located, 3 confirmed; sketch 0.210 s, locate 0.350 s, align 1.800 s

If pruning was turned off for the query, the line ends with `; pruning off`.

`-vv` also writes JSON debug lines to standard error. Each line names its index file in `shard`
and its event type in `kind`:

- `possible`: stage 1 kept the assembly; includes its coverage upper bound and containment
  estimate.
- `located`: stage 2 found alignment windows; includes their count and total bases.
- `dropped`: the assembly was dropped from the search; includes the `stage` and the `reason`:
  `coverage_bound` (it cannot reach the required coverage), `no_seed_hits`, `decision` (its
  alignments failed the decision rule), or `resource_limit`.
- `unfinished`: the search stopped before the assembly could be evaluated.
- `summary`: full summary for the index file, including the work done in each stage and the
  stage times.

Debug lines are useful for finding out why a specific assembly was or was not reported. Their
JSON keys may change between versions.

## The decision rule

An assembly is reported when its alignments cover enough of the query at a high enough
identity (`-c` and `-i`).

Only the stretches of an alignment that reach `--stretch-identity` and are at least
`--stretch-length` columns long count toward coverage. A poor stretch inside a fragment adds
nothing. Fragments are chosen so that no query base is counted twice.

The contig rule is off by default. It keeps out contigs that only share a gene or a repeat with
the query: with `--contig-coverage`, a contig counts only if the alignments cover at least that
share of the contig itself, and with `--contig-identity`, only if they reach that identity.

Low-complexity stretches of the query (repeats of a few bases, as found by DUST masking) count
toward no coverage. A query that is mostly low-complexity cannot reach the required coverage
with its other bases. It gets no answer, rather than false hits.

## The miss chance

jam never aligns every assembly. It may miss an assembly that qualifies in two ways:

1. **Seeds.** The copy in the assembly may have no seed in common with the query, because
   every seed in it carries a substitution. Then stage 1 never sees it.
2. **Pruning.** Between the stages, jam drops assemblies whose seed hits leave too much of the
   query uncovered. A long stretch without seed hits may still hold the query, with bad luck.

`--max-miss` (default 0.0001) is the largest accepted total chance of both. jam computes a
proven upper bound U for the first part: the chance that a copy of the required length shows
no seed hit. The rest, `--max-miss` minus U, goes to the pruning. jam picks the shortest
allowance for seedless stretches that keeps the pruning within that rest. When U alone reaches
`--max-miss` (a short query, or a low `--identity`), jam turns the pruning off for that query
and reports `"prune":false`. The summary key `miss` holds the total for the query.

The bound holds under these assumptions:

- the qualifying copy is one continuous piece of at least the required length, without
  insertions or deletions;
- its substitutions fall at random places, independently of each other.

Repeats, indels or clustered changes can miss more often. Under these assumptions, jam checks
blocks of w + k - 1 bases, each of which holds a whole minimizer window and so at least one
seed. It evaluates the bound at an identity slightly below `--identity`, so that it covers any
copy that qualifies.

`--seed-mismatches auto` uses the same bound. It takes 1 for a query when U with exact seeds
is above `--max-miss`, one mismatch lowers it, and the extra work stays within about 10 CPU
seconds (estimated as 1 ms per query k-mer). It takes 2 instead when also the identity is below
0.9, the query has at most 1,800 k-mers, U with one mismatch is still above `--max-miss`, and
two mismatches lower it.

## Index files

An index is one or more shards. A shard is two kinds of file:

- `NAME.jam`, the index file: the seeds of all entries, then the manifest. The sketch part of the
  seeds answers stage 1, the locator part stage 2. The manifest holds the parameters, the
  entries, and the size and SHA-256 digest of each sequence file.
- `NAME.jamseq`, the sequence file: the sequences of all entries, 2 bits per base, one entry
  after another. Stage 3 reads them. The index file finds it next to itself. With `--chunks N`
  there are several, `NAME.1.jamseq` to `NAME.N.jamseq`.

## S3 and local caches

An index can lie on S3. Copy it there with `jam index copy`:

    jam index copy assemblies.jam s3://bucket/indexes --cache /scratch/cache

This uploads the sequence files first and the index file last, each in parts that go up at
once, so a reader never sees an index without its sequence files. `jam index copy` works the
other way too, and between two places of the same kind; a download is checked against the sizes
and digests in the index file. Any other S3 tool works as well, as long as the files keep their
names.

Then give `s3://bucket/prefix/NAME.jam` instead of a local index file to `jam search`, `jam
screen`, `jam serve`, `jam stats` and `jam index warm`. Credentials come from the usual `AWS_*`
environment variables, and the endpoint from `AWS_ENDPOINT_URL`. jam never sends HEAD requests:
it reads the size of the index file from its header, and the sizes of the sequence files from
its manifest.

jam learns the ETag of each file from its first read and sends it as `If-Match` with every later
read. Cached blocks and pins belong to one ETag. With `--cache`, the ETags are kept in the cache
directory, and the first use of each file in a run checks it with one byte. `jam index copy
--cache` keeps the ETags of the files it uploads there too. If a file was uploaded again, the
search ends with `the index changed in S3 while it was read; retry`. Run it again, and it reads
the new files. `jam serve` opens the index again by itself. A pin made for the old file is
removed; run `jam index warm` again to pin the new file.

Reading from S3 is slow for the first requests. Two caches help:

- **Pins.** `jam index warm` copies byte ranges of the index file and the sequence files into a
  local cache directory. A search or a server reads them from there, from RAM or from local
  disk. Pins are never evicted.
- **Block cache.** Everything not pinned is read in blocks. The blocks stay in RAM, and on local
  disk when you give a disk budget. They are kept across runs, so a second search is faster.

Example:

    jam index warm s3://bucket/indexes/assemblies.jam --cache /scratch/cache --memory 16GiB \
        --disk 200GiB
    jam search s3://bucket/indexes/assemblies.jam plasmid.fa --cache /scratch/cache \
        --memory 16GiB --disk 200GiB

By default `jam index warm` pins the sketch part, then the locator part, while the disk budget
lasts. Choose parts yourself with `--pin`, for example `--pin sketch,locators,packed` or
`--pin sketch,packed=NAME1,NAME2`.

## Screening assemblies

`jam screen` looks for many small sequences in large assemblies. Build the index with one entry
per record, then stream the assemblies past it:

    jam index build --singleton -o plasmids plasmids.fa
    jam screen plasmids.jam --assemblies assemblies.txt -c 0.8 -o screen.jsonl

It writes JSON lines. The key `type` names the line:

- `entry`: a small sequence (`entry`) found in an assembly (`assembly`). The line holds its
  length (`entry_length`), the covered bases and their share (`covered_bases`,
  `query_coverage`), the identity and the base counts, the assembly contigs it lies on
  (`contigs`, each with its name, length, aligned bases, identity, `score`, `best_score` and the
  number of entries it counts for), the fragments (`primary`, `alternatives`) and the uncovered
  stretches (`gaps`).
- `assembly`: the summary of one assembly: its contigs and bases, the work counts and the
  seconds of the survey and the alignments.
- `failed`: an assembly that could not be read, with the `error`. The other assemblies go on,
  and jam exits with an error at the end.

## Server

`jam serve` keeps an index open and answers searches from other machines. The machines connect
over [iroh](https://www.iroh.computer), with a secret key file for each server:

    jam serve assemblies.jam --key server.key

The key file is created on first start. Its public key is the node id. `jam serve` prints the
address as `NODE_ID@ADDRESS`. Search on the server with:

    jam search --server NODE_ID@ADDRESS plasmid.fa

The decision options work as for a local search, and the output is the same.

With `--peer`, the server is a coordinator. Every search also runs on each peer, and the
coordinator merges the results. This spreads a large index over several machines, one shard
each.

## Speed and size

These numbers are from one test on 2026-10-22, one run per case, so treat them as rough. The
collection had 1,000 SPIRE assemblies (241.8 Gbp). Searches ran with 8 threads and a cold file
cache.

- Building the index took about 10 minutes with 32 threads. The index was 108 GB (68 GB packed
  sequences, 28 GB locators, 12 GB sketch), measured with an earlier layout of the same data.
- On local NVMe, queries of 1 kb to 10 kb finished in about 1 to 10 seconds. A query of 2 Mb
  with no hit took about 40 seconds.
- On network storage (Ceph), the same searches took about 30 seconds to 3 minutes.

## License

jam is released under the MIT license. See `LICENSE`.
