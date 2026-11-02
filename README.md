<!-- What jam does and how to use its commands. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam

jam (just another minimizer) searches for query sequences (plasmids, genes, contigs; 1 kb to a
few Mb) in very large collections of assemblies. It computes full base-level alignments with
CIGAR strings. jam decides each candidate assembly and reports those that pass.

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
that contains the query. Standard output holds only these matches.

## How it works

Each assembly in the index is stored as an *entry*. With `--singleton`, each FASTA record
becomes a separate entry instead.

The two commands name things differently. `jam search` calls every indexed entry an assembly.
`jam screen` calls an indexed sequence an entry, and a streamed file an assembly.

A search runs in three stages and answers one question per
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
  one entry. Without `-k` and `-w`, such an index uses k 21 and w 10, the dense seeds that
  `jam screen` needs.
- `--filter` (needs `--singleton`): also write the filter file `NAME.jamfilter` for
  `jam screen`. It is meant for indexes of small sequences and warns above 8 GiB.
- `-t, --threads N` (default: all cores): number of threads to use.

Advanced options:

- `-k, --k K` (default 25, or 21 with `--singleton`): k-mer length in bases (1 to 32).
- `-w, --window W` (default 100, or 10 with `--singleton`): number of k-mers per minimizer
  window.
- `--sketch-only`: write only the index file containing the sketch data. This index is much
  smaller, but it supports only `--level sketch` searches. It writes no `.jamseq` file, so
  `--chunks` has no effect with it.
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
query file, for example `all/shard-*.jam`. `jam search` needs the `.jam` files. It does not
accept the output directory itself. Only `jam index reshard` reads the `index.json` in it.

### jam index warm

    jam index warm INDEX --cache DIR --disk SIZE [--memory SIZE] [--pin PARTS]

Copies byte ranges of an S3 index into a local cache directory so that later searches read them
from RAM or local disk. See [S3 and local caches](#s3-and-local-caches).

- `--cache DIR`: local cache directory to populate. `jam search` and `jam serve` read pinned
  data from this directory.
- `--memory SIZE`: RAM budget for pinned data and the block cache (for example, `16GiB`). By
  default, half of the RAM available to jam, less 8 GiB for the process and the open index
  files and 1.5 GiB for the search, and at least 1 GiB. The RAM available is the memory limit of jam's cgroup (as under Slurm or in a container), or
  the available RAM of the machine if that is less.
- `--disk SIZE`: disk budget for pinned data and the block cache (for example, `200GiB`).

A size is a number with an optional unit, such as `500`, `1.5GiB` or `200 GB`. A number
without a unit means bytes. Units are not case sensitive. `KB`, `MB`, `GB`, `TB`, `PB` and
`EB` (also `K`, `M`, `G`, `T`, `P` and `E`) are powers of 1000. `KiB`, `MiB`, `GiB`, `TiB`,
`PiB` and `EiB` (also `Ki`, `Mi`, `Gi`, `Ti`, `Pi` and `Ei`) are powers of 1024. The same
format applies to `--memory` and `--disk` in `jam search` and `jam serve` too.
- `--pin PARTS`: index sections to pin instead of the default plan. Give a comma-separated list
  of `sketch`, `locators`, `headers`, `packed`, or `packed=NAME,NAME` (to pin sequences for
  specific entries). `headers` pins only the first 4 KiB page of each packed sequence, not its
  bases. Because every item after `packed=` is read as an entry name, that item must come last
  in the list. You can pass this flag multiple times.

### jam index copy

    jam index copy SOURCE TARGET [--cache DIR]

Copies an index file and all of its sequence files. Both source and target may be a local path
or `s3://bucket/prefix`. `SOURCE` must be the index file `NAME.jam`. `TARGET` may be a
directory, `s3://bucket/prefix`, or the full path to the new index file `NAME.jam`.

- `--cache DIR`: local cache directory in which to record the S3 ETag of each uploaded file.
  Searches that use this cache directory will then verify from their first read that the files
  on S3 match the uploaded version.

### jam index filter

    jam index filter NAME.jam [-t N]

Writes the filter file `NAME.jamfilter` of a local index file (see [Index files](#index-files)),
as `jam index build --filter` does, for an index built without it.

- `-t, --threads N` (default: all cores): number of threads to use.

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
- `-v, --verbose`: write one progress line per query to standard error, with its status, counts
  and miss probability. Pass twice (`-vv`) to also write per-assembly debug lines.

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
- `--max-miss P` (default 0.001): the miss target, the total chance, per query and matching
  copy, that jam may lose a copy; see [The miss chance](#the-miss-chance).
- `--deadline SECONDS`: stop searching a query after this many seconds. jam then reports the
  query as `stopped` on standard error.
- `--server NODE_ID@ADDRESS`: run the search on a `jam serve` server instead of local index
  files.
- `--cache DIR`, `--memory SIZE`, `--disk SIZE`: local cache directory and memory or disk
  budgets for index files on S3; see [S3 and local caches](#s3-and-local-caches). `--memory`
  defaults to the budget recorded in the cache directory. Without one, it is half of the RAM
  available to jam, less 8 GiB for the process and the open index files and 1.5 GiB for each
  index file searched at the same time, and at least 1 GiB. The pins, the block cache and the buffers of S3 reads in flight share this budget;
  `-v` prints how. `--disk` defaults to the recorded budget, or 0 if there is none.
- `--s3-requests N` (default 512, at most 4096): number of S3 requests in flight at the same
  time. S3 answers many parallel requests best; lower it for a store that fails under many.
  Each request needs its own connection, so jam lowers it if the open files limit is too low.
- `--shards-at-once N` (default 4 when an index file is on S3, else 1): number of index files
  searched at the same time, so the store gets more requests at once. Each search has its own
  `--threads`; on S3 they mostly wait for reads. Results are still written in the order of the
  index files.

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
- `--cache DIR`, `--memory SIZE`, `--disk SIZE`, `--s3-requests N`: same as for `jam search`.

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
- `-t, --threads N` (default: all cores): number of threads. With 16 or fewer, all threads form
  one pool, and up to N assemblies are open at once (one per thread, at most the number of
  assemblies), so `-t 2` screens two assemblies at once. With more than 16, the threads are split
  into ceil(N/16) pools of 16 (the last one smaller), and each pool screens one assembly at a
  time, so `-t 32` screens two assemblies at once and `-t 64` four. In both cases a pool works on
  the batches of its assemblies in parallel.
- `-o, --output FILE`: write results to this file instead of standard output.

Advanced options:

- `--stretch-identity`, `--stretch-length`, `--contig-coverage`, `--contig-identity`: same as
  for `jam search`.
- `--closest-margin M` (0 to 1): for each contig, score every matching entry as the aligned
  share of the contig multiplied by the alignment identity. A contig then counts only toward
  entries whose score is at least the highest score on that contig minus M. M is a plain
  difference, not a ratio: with M 0.05 and a best score of 0.90, entries that score 0.85 or
  more count. Without this option, a contig counts toward every entry it qualifies for.
- `--batch-bases N` (default 500000): number of contig bases searched together in one batch.
- `--contig-dir DIR`: write contigs that share enough with an entry into this directory, as one
  gzipped FASTA file per assembly.
- `--survey-only`: stop after the survey stage. No alignments are computed and no entries are
  reported; only the files in `--contig-dir` are written.
- `--low-complexity`: also build walks from low-complexity seeds (such as repeats of one or two
  bases), as `jam search` always does. Off by default, because such matches mean little and can
  make one assembly take days. Without it, results can differ from earlier builds: alignments and
  their coverage, the per-contig `entries` and `best_score`, and which entries qualify (records can
  disappear or appear). With it, the screen gives exactly the output of earlier builds.
- `--compact`: write one `pair` line per small sequence and assembly instead of the `entry`
  line, without CIGAR strings and gaps, for very large result sets. See [Screening assemblies](#screening-assemblies).

## Output

`jam search` writes only matches (query x assembly) to standard output (or to `-o FILE`). A query
with no match writes nothing. Progress, problem and debug lines go to standard error.

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

- `c`: the containment estimate, recomputed from the seeds whose positions were located. It
  can be lower than at stage 1. Stage 1 can count a seed hit that the entry does not have,
  because the index keeps only 44 bits of each seed hash, so two k-mers can share a key.
  Stage 2 counts a seed only if the index lists a location for it in the entry.
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

### Problems

jam always reports a query that does not end `complete` on standard error, even without `-v`. It
writes one line with the query name, the status and the reason:

    plasmid1: stopped (deadline)

For `failed`, the line ends with `: ` and the error message. With several index files, the worst
status wins. With `-v`, the progress line below already holds this information, so there is no
extra line.

The status values:

| `status` | `reason` | Meaning |
| --- | --- | --- |
| `complete` | | Every candidate assembly was evaluated. |
| `incomplete` | `resource_limit` | An assembly required more memory or work than a limit allows (for example, when a periodic query matches many periodic regions) and could not be evaluated. |
| `stopped` | `deadline` | The search exceeded `--deadline`. |
| `stopped` | `asked` | The search was asked to stop (for example, when the program reading standard output closed the pipe). |
| `failed` | `error` | A search stage or the server failed. The message gives the reason. |

Whenever a query does not end with `complete`, jam writes all remaining output and then exits
with an error. This applies to every format.

### TSV

`--format tsv` writes a header line followed by one tab-separated line per confirmed assembly:

| Column | Meaning |
| --- | --- |
| `query` | The query name. |
| `assembly` | The assembly (entry) name. |
| `coverage` | The query coverage, as in `cov` above. |
| `identity` | The sequence identity, as in `id` above. |
| `bases` | The number of covered query bases. |
| `fragments` | The number of aligned fragments in `f`. Alternative fragments (`alt`) are not counted. |

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

    plasmid1: complete; 120 seeded, 14 possible, 9 located, 3 confirmed; miss 1e-4; sketch 0.210 s, locate 0.350 s, align 1.800 s

The counts are the assemblies at each step. Seeded assemblies share at least one seed with the
query. Possible assemblies passed stage 1, located assemblies passed stage 2, and confirmed
assemblies passed stage 3 and were reported. `miss` is the modeled probability that the search
missed a qualifying assembly; see [The miss chance](#the-miss-chance). The times are seconds. The
sketch time is wall-clock time, including query preparation. The locate and align times are
summed across threads. With several index files, the counts are summed and `miss` is the largest.

If pruning was turned off for the query, the line ends with `; pruning off`.

When `--identity` or `--stretch-identity` is below 0.8, `jam search` and `jam screen` write
one note to standard error at the start, and `jam serve` logs it once: such settings are
allowed, but slower and not validated, as the benchmarks and the tuning cover 0.8 and above.

After the last query, jam writes one hint line to standard error, even without `-v`, if some
queries could not reach `--max-miss`. It never goes to standard output. For example:

    jam: 3 of 10 queries are above the miss target --max-miss 0.001 at -i 0.8 (worst-case bound up to 0.39); at this identity the bound reaches it from about 83 kb. -v shows each query.

The length is the smallest one from which every longer query reaches `--max-miss` with the
lookups jam chooses. A search with `--server` writes the line without the length.

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
sequence identity (`-c` and `-i`).

Within an alignment, only local stretches that reach `--stretch-identity` and span at least
`--stretch-length` alignment columns count toward query coverage. A low-identity stretch inside
a fragment adds nothing. When multiple fragments align, jam chooses fragments so that no query
base is counted twice.

The contig filters are off by default. They exclude contigs that share only a short gene or a
repeat with the query: with `--contig-coverage`, a contig counts only if alignments cover at
least that share of the contig itself, and with `--contig-identity`, only if the alignments on
that contig reach at least that identity.

Low-complexity stretches of the query (short repeats identified by DUST masking) count toward
coverage like other bases. But at least half of the required coverage must come from normal bases.
A query with fewer normal bases than that gets no answer. In that case, jam reports no hits
rather than returning false positives.

## The miss chance

jam does not align the query against every assembly. It can miss a qualifying assembly in two
ways:

1. **Seeds.** The matching copy in the assembly may share no seeds with the query if every seed
   in the copy contains a substitution. When this happens, stage 1 never sees the assembly.
2. **Pruning.** Between stages, jam drops assemblies whose seed hits leave too much of the
   query uncovered. By chance, substitutions may eliminate all seed hits across a long region
   that still contains a valid match.

`--max-miss` (default 0.001) sets the miss target, the combined probability of both events. jam
first computes a proven upper bound U for the first event: the probability that a qualifying
copy of the required length produces no seed hits at all. The remaining budget, `--max-miss`
minus U, goes to pruning. jam picks the tightest allowance for unseeded stretches that keeps
the pruning miss probability within that remaining budget. When U alone reaches `--max-miss`
(for example, with a short query or a low `--identity`), jam turns pruning off for that query
and marks the `-v` line with `pruning off`. The `miss` value in that line is the total miss
probability for the query.

This upper bound holds under two assumptions:

- the qualifying copy is one continuous segment of at least the required length, with no
  insertions or deletions;
- substitutions occur independently at random positions.

Repeats, indels, or clustered mutations can cause higher miss rates. Under these assumptions,
the bound uses non-overlapping blocks of w + k - 1 bases. Each block spans a complete minimizer
window and therefore contains at least one seed. jam evaluates the bound at a sequence identity
slightly below `--identity` so that the bound covers any copy that qualifies.

jam also chooses how many mismatches the seed lookups allow, for each query. There are four
levels, from the cheapest to the most thorough:

- level 0: exact seeds only;
- level 1: also seeds with one mismatch;
- level 2: also seeds with two mismatches;
- level 3: level 2, plus seeds one insertion or deletion away.

jam first tries to prove the miss target. It takes the cheapest level whose proven bound meets
`--max-miss`, with that level on every query k-mer. Only cheap levels count for this: levels 0
and 1 at any length, and level 2 for queries of at most 8,226 k-mers. That is the number of
k-mers at which the plan below needs level 2 on every k-mer, at 80% identity and coverage 0.5.
With one hit start (see below), the bound and the pruning must meet a tenth of `--max-miss`
(0.0001 at the default, the target of the release), so these plans keep the levels of the
release. With two hit starts they may use all of `--max-miss`.

When the level is proven, jam may also ask each candidate assembly for two seed-hit starts
instead of one. A hit start is a query position with a seed hit in the assembly, from its own
seeds or from frequent seeds whose assemblies include it. Lookups with mismatches at one query
position count once. jam asks for two when the proven chance that a qualifying copy shows fewer
than two still meets `--max-miss` with the same pruning allowance. It only tries two at a level
whose bound with one hit start meets `--max-miss`. Then U above is this chance. It is computed
exactly: for each copy length L from the required length up to the query length,
but below twice the required length (a longer copy holds a qualifying piece of such a length),
the copy holds floor(L / (w + k - 1)) blocks with one seed
each, and its floor((1 - identity) L) substitutions fall on random positions. Two hit starts
only drop candidates, and a lower level only finds fewer hits, so no other proven plan keeps
fewer candidates. Candidates with fewer hit starts are dropped in stage 1, before any of their
locations or sequences are read, but only when a query has more than 16 of them in a shard.
With 16 or fewer, jam keeps them all, as they cost little to check. For a gene of 1,299 bases
at 95% identity and coverage 0.5, level 2 with two hit starts has a bound of 0.00045, or
0.00087 with the pruning.

Otherwise jam plans the lookups by seed hits. A copy that just qualifies covers the required
part of the query (coverage times length) at the identity limit. jam takes the cheapest plan at
which such a copy is expected to show 4 seed hits on the required part from 90% identity up, and
8 below, where divergent partial copies need more. All query k-mers get
exact lookups. When that is not enough, jam adds level 1 on every query position (it costs only
76 lookups per position), then level 2 for just enough of the query positions. These positions
are spread over the query in strata of 32 positions, picked by
a hash of the position, so the same query always gets the same lookups and no period of the
query lines up with them. If even level 2 on every k-mer is not enough (a short query, or a
small or fixed `-c`), jam uses it on every k-mer. For a query of at most 1,100 k-mers, the
positions with level 2 also get the lookups one insertion or deletion away. The targets of 4
and 8 hits were chosen on a separate set of genes, plasmids and reads; they are not part of the
proof.

jam makes the lookups in batches and keeps only those whose seed keys the index holds, so memory
does not grow with the number of lookups. Before it aligns a candidate assembly with few seed
hits (at most 8 query positions), jam checks them on the assembly sequence: the bin of each seed
hit, with k bases on each side, must hold the looked-up k-mer on either strand. A candidate
whose seed hits all fail this check only matched seed keys by chance, and is dropped with the
reason `no_seed_hits`.

`--max-miss` is a target, not a promise. jam can prove it only for some combinations of query
length and identity. Short queries at low identity cannot meet it as a guarantee. With the
default 0.001 and coverage 0.5, jam proves it from these query lengths:

| Identity | Met from |
| -------- | -------- |
| 95% | 1.24 kb |
| 90% | 3.22 kb |
| 85% | 23.6 kb, and also 8.18 to 8.25 kb |
| 80% | 82.6 kb |

Below these lengths jam still searches normally and finds most copies in practice. The proven
worst case is only higher there. `-v` reports the value of each query and marks it when it is
above `--max-miss`. jam also writes one hint line to standard error after the run (see
[Progress and debug lines](#progress-and-debug-lines)). The `seed_mismatches` field of the `-vv`
summary holds the highest level used, and `sketch.sampled_kmers` the query positions expected to
get it when only some do. `sketch.min_hit_starts` holds the hit starts each candidate needed (1
or 2), `sketch.few_hit_entries` the assemblies dropped for fewer, and `sketch.miss_upper_bound`
and `sketch.total_miss` the seed bound and the total with the pruning.

## Index files

An index consists of one or more shards. Each shard consists of two kinds of files:

- `NAME.jam`, the index file: stores the seeds of all entries, followed by the manifest. The
  sketch section of the seeds answers stage 1, and the locator section answers stage 2. The
  manifest stores the index parameters, the entry list, and the byte size and SHA-256 digest of
  each sequence file.
- `NAME.jamseq`, the sequence file: stores the sequences of all entries, packed at 2 bits per
  base, one entry after another. Stage 3 reads this file from the same directory as the index
  file. With `--chunks N`, the sequences are split into several files named `NAME.1.jamseq`
  through `NAME.N.jamseq`.
- `NAME.jamfilter`, the filter file, only for `jam screen` and written only on request
  (`jam index build --filter` or `jam index filter`): filters of the stored seed keys, about
  1.5 bytes per seed, one for each range of sketch buckets. `jam screen` maps it, so processes
  on one machine share it, and builds the filters in memory when the file is missing or does not
  match the index. It is meant for indexes of small sequences; above 8 GiB both commands warn on
  standard error and go on. `jam search` never reads it.

## S3 and local caches

An index can be stored on S3. Upload it with `jam index copy`:

    jam index copy assemblies.jam s3://bucket/indexes --cache /scratch/cache

This uploads the sequence files first and the index file last, transferring each file in
parallel parts so that a reader never sees an index file before its sequence files are in
place. `jam index copy` also downloads from S3 to a local path, or copies between two locations
of the same kind. Downloads are verified against the file sizes and SHA-256 digests stored in
the index manifest. You can also copy the files with any other S3 tool, as long as the file
names stay the same.

Once uploaded, pass `s3://bucket/prefix/NAME.jam` instead of a local index file to `jam search`, `jam
screen`, `jam serve`, `jam stats`, and `jam index warm`. jam reads credentials from the
standard `AWS_*` environment variables and reads a custom endpoint from `AWS_ENDPOINT_URL`. It
never sends S3 HEAD requests: it reads the size of the index file from the file header and the
sizes of the sequence files from the manifest.

On its first read of each S3 file, jam records the file's ETag and sends that ETag as
`If-Match` on every subsequent read. Cached blocks and pinned data are tied to a specific ETag.
When `--cache` is used, jam saves the ETags in the cache directory, and the first read of each
file in a run verifies the ETag with a 1-byte request. `jam index copy
--cache` also saves the ETags of the files it uploads in the cache directory. If a file on S3
is replaced while it is being used, the search stops with the error
`the index changed in S3 while it was read; retry`. Run the command again to read the new
files. `jam serve` reopens the updated index automatically. Pinned data from the old version is
removed; run `jam index warm` again to pin the updated file.

Reading directly from S3 is slow on the first requests. Two local caches speed this up:

- **Pins.** `jam index warm` copies selected byte ranges of the index file and sequence files
  into a local cache directory. Searches and servers read those ranges from RAM or local disk.
  Pinned data is never evicted.
- **Block cache.** Any data that is not pinned is read in blocks. Those blocks stay in RAM, and
  also on local disk when a disk budget is given. Disk blocks are kept across runs, making
  subsequent searches faster.

Example:

    jam index warm s3://bucket/indexes/assemblies.jam --cache /scratch/cache --memory 16GiB \
        --disk 200GiB
    jam search s3://bucket/indexes/assemblies.jam plasmid.fa --cache /scratch/cache \
        --memory 16GiB --disk 200GiB

By default, `jam index warm` pins the parts that save the most S3 requests per byte first: the
headers, then the sketch section, then the locator section. Each part goes in RAM if it fits
the memory budget with 96 MiB left for the block cache, else on disk. A part that does not fit
the disk budget is left out. A part in RAM also needs its room on disk, as it is loaded from
the cache directory at the start. The block cache gets the RAM and disk that the pins leave.
Use `--pin` to choose which sections to pin, for example `--pin sketch,locators,packed` or
`--pin sketch,packed=NAME1,NAME2`.

The headers are the 4 KiB header page of each packed sequence. Before jam reads bases of an
entry, it reads and checks that page. With `headers` pinned, the page comes from RAM or local
disk, so a candidate entry needs one S3 request fewer. The headers take 4 KiB per entry: about
1.1 GB for 270,000 entries. Packed sequences always go on disk.

Without `--cache`, a search of index files on S3 still keeps the blocks it reads in RAM, up to
the memory budget, so later queries and later index files reuse them.

## Screening assemblies

`jam screen` searches for many small sequences across large assemblies. Build an index with one
entry per FASTA record, then stream the assemblies against it:

    jam index build --singleton -o plasmids plasmids.fa
    jam screen plasmids.jam --assemblies assemblies.txt -c 0.8 -o screen.jsonl

The default singleton index is the dense one (k 21, w 10). It holds about 7.7 times more seed
data than k 25 and w 100. An older sparse index still works, as the scheme is stored in the
manifest. The filter file `NAME.jamfilter` (`jam index build --singleton --filter`) saves each
screen the time to build its seed key filter; `jam index filter` writes it for an existing index.

A screen on large assemblies runs best with large freed buffers kept in the usual memory arenas.
Set `_RJEM_MALLOC_CONF=oversize_threshold:0` in the environment of `jam screen`. The binary does
not set it, as it held 6 GiB more at the start of a search on a large index.

A screen has two stages. The first is the *survey*. jam looks up the seeds of the contigs of an
assembly, a batch at a time, in the index. It keeps the entries that share enough seeds with a
contig and aligns only those. `--survey-only` stops after the survey.

The command writes JSON lines, where the `type` key identifies the line kind:

- `entry`: reports a small sequence (`entry`) found in an assembly (`assembly`). The line
  gives the sequence length (`entry_length`), the covered bases and their share
  (`covered_bases`, `query_coverage`), the sequence identity, and the alignment base counts. It
  also lists the assembly contigs that the sequence aligns to (`contigs`), where each contig
  item includes its name, length, aligned bases, identity, `score`, `best_score`, and the
  number of entries it counts toward. Finally, it lists the aligned fragments (`primary`,
  `alternatives`) and the uncovered query intervals (`gaps`).
- `pair` (with `--compact`, instead of `entry`): the same counts as `entry`, with each contig
  as an array (id, name, length, aligned bases, identity, `score`, `best_score`, entries) and
  each fragment as an array (contig id, contig start and end, strand `+` or `-`, start and end
  on the small sequence, identity, matches, and 1 for a `primary` or 0 for an alternative
  fragment). There are no CIGAR strings and no `gaps`.
- `assembly`: summary for one screened assembly, including its contig and base counts, work
  counts, and the time in seconds spent in the survey and alignment stages.
- `failed`: reports an assembly that could not be read, along with the `error` message.
  Screening continues with the remaining assemblies, and jam exits with an error at the end.

Example lines, one of each type (the long lines are wrapped here):

    {"type":"entry","entry":"pl1","entry_length":1200,"assembly":"asm1.fa","covered_bases":1200,
     "query_coverage":1.0,"identity":0.9975,"matches":1197,"substitutions":3,"insertions":0,
     "deletions":0,"contigs":[{"id":0,"name":"ctg_1","length":3000,"aligned":1200,
     "identity":0.9975,"score":0.399,"best_score":0.399,"entries":1}],
     "primary":[{"new_bases":1200,"fragment":{"contig_index":0,
     "query_segments":[{"start":0,"end":1200}],"alignment":{"score":2385,"strand":"forward",
     "query_interval":{"start":0,"end":1200},"contig_interval":{"start":1000,"end":2200},
     "matches":1197,"substitutions":3,"insertions":0,"deletions":0,"identity":0.9975,
     "cigar":"300=1X299=1X299=1X299="}}}],"alternatives":[],"gaps":[]}
    {"type":"assembly","assembly":"asm1.fa","contigs":1,"bases":3000,"skipped_contigs":0,
     "batches":1,"surveyed_entries":1,"aligned_entries":1,"fragments":1,"alignments":1,
     "reported_entries":1,"matched_contigs":1,"matched_bases":3000,"reported_contigs":1,
     "reported_bases":3000,"survey_seconds":0.000429435,"align_seconds":0.002707405}
    {"type":"failed","assembly":"missing.fa","error":"bad assembly: /data/missing.fa: No such file or directory (os error 2)"}

## Server

`jam serve` keeps an index open and answers searches from other machines. Clients connect over
[iroh](https://www.iroh.computer), using a secret key file for each server:

    jam serve assemblies.jam --key server.key

If the key file does not exist, jam creates it on first start. Its public key is the server's
node ID. `jam serve` prints its address as `NODE_ID@ADDRESS`. Search against the server with:

    jam search --server NODE_ID@ADDRESS plasmid.fa

All decision options work the same way as in a local search, and the output is identical.

With `--peer`, the server also acts as a coordinator. Every search runs on the coordinator and
on each peer, and the coordinator merges the results. This lets you spread a large index across
several machines, with one shard per machine.

## Speed and size

These numbers come from a single test on 2026-10-22, with one run per case, so treat them as
rough estimates. The collection held 1,000 SPIRE metagenome assemblies (241.8 Gbp). Searches ran with 8
threads and a cold file cache.

- Building the index took about 10 minutes with 32 threads. The index totaled 108 GB (68 GB
  packed sequences, 28 GB locators, and 12 GB sketch), measured on an earlier layout of the
  same data.
- On local NVMe storage, queries of 1 kb to 10 kb finished in about 1 to 10 seconds. A 2 Mb
  query with no hits took about 40 seconds.
- On network storage (Ceph), the same searches took about 30 seconds to 3 minutes.

## License

jam is released under the MIT license. See `LICENSE`.
