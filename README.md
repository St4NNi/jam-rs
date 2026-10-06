<!-- What jam does and how to use its commands. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam

jam (just another minimizer) searches and aligns query sequences (plasmids, genes, contigs; 1 kb
to a few Mb) in very large collections of assemblies. It reports full base-level alignments with
CIGAR strings and a decision for every assembly.

## Install

Pre-built binaries will be linked here once there is a release.

Until then, build jam from source. It is a Rust program (edition 2024, so a recent stable
toolchain is needed). Install Rust with [rustup](https://rustup.rs), then:

    cargo build --release --locked

The binary is `target/release/jam`. Run `jam --help` or `jam COMMAND --help` for all options.

The block cache uses a patched copy of the `foyer` crate. `Cargo.toml` takes it from a git
repository (`[patch.crates-io]`), so the first build needs network access to GitHub. For an
offline build, run `cargo vendor` on a machine with network access and build with
`--offline`.

The `Dockerfile` builds the same binary into a small container image.

## Quick start

Build an index of your assemblies (plain or compressed FASTA). Each file becomes one entry:

    jam index build -o assemblies -t 8 assemblies/*.fa.gz

This writes the index file `assemblies.jam` and the sequence file `assemblies.jamseq`.

Search the records of a query file in it:

    jam search assemblies.jam plasmid.fa > hits.jsonl

Each query record is searched on its own. jam writes one JSON line per assembly that holds the
query, then one summary line per query.

## How it works

Each assembly in the index is an *entry*. With `--singleton`, each FASTA record is an entry
instead. A search answers three questions, one per stage:

1. **Sketch**: which entries may contain the query? jam looks up the query's seeds (minimizers)
   in the index. Each entry gets a containment estimate: the share of the query it seems to
   hold, from shared seeds alone.
2. **Locate**: where in those entries do the shared seeds lie? The result is a list of
   alignment windows on the contigs of each entry.
3. **Align**: how exactly does the query align there? jam aligns the query in the windows. The
   entry is reported when the alignments pass the decision rule.

Each stage reads more data than the one before and only works on what the previous stage kept.
A search can stop after any stage with `--level`:

- `--level sketch` is a containment search, like `sourmash search --containment`. It is fast
  and only needs the index file, but the estimate falls with mutations and nothing is aligned.
- `--level locate` reports candidate regions: the windows on each contig where the query may
  lie, without aligning anything.
- `--level align` (the default) aligns and decides.

Between the stages, jam drops entries that cannot reach the required query coverage any more.
This pruning is safe up to a small modeled chance, see [The miss chance](#the-miss-chance).

## Commands

| Command | Purpose |
| --- | --- |
| `jam index build` | Build an index file and its sequence file from FASTA files. |
| `jam index reshard` | Merge and split index files into a directory of new index files. |
| `jam index warm` | Pin parts of an index in a local cache directory. |
| `jam index copy` | Copy an index between local directories and S3. |
| `jam stats` | Show what an index file holds. |
| `jam search` | Search query sequences and align them. |
| `jam serve` | Serve an index over the network for `jam search --server`. |
| `jam screen` | Stream large assemblies past an index of small sequences. |

### jam search

    jam search [OPTIONS] INDEX... QUERY

Searches each record of the FASTA file `QUERY` in every index file. The name of a query is the
first word of its FASTA header. With `--server`, give only the query file.

Normal options:

- `-o, --output FILE`: write to this file instead of the standard output.
- `--format json|tsv|paf` (default `json`): the output format, see [Output](#output). `tsv` and
  `paf` need `--level align`.
- `-t, --threads N` (default: all cores): threads to use.
- `-c, --coverage C` (default 0.5): the least share of the query that the alignments must
  cover. Give a share such as `0.8`, or a base count such as `800bp`.
- `-i, --identity I` (default 0.95): the least identity over the covered query bases.
- `--circular`: treat the queries as circular. A fragment over the query end then ends past the
  query length. The default is linear.
- `-v, --verbose`: write a progress line per query to the standard error. Give it twice (`-vv`)
  to also write the debug lines.

Advanced options:

- `--stretch-identity I` (default: the value of `--identity`): the least identity of a stretch
  of an alignment to count toward coverage.
- `--stretch-length N` (default 40): the least number of alignment columns of a stretch to
  count toward coverage.
- `--contig-coverage C` (default 0, off): a contig counts only if the alignments cover at least
  this share of the contig itself.
- `--contig-identity I` (default 0, off): a contig counts only if its alignments reach at least
  this identity.
- `--level sketch|locate|align` (default `align`): the last stage that runs, see
  [How it works](#how-it-works).
- `--max-miss P` (default 0.0001): the largest accepted chance of missing an assembly that
  qualifies, see [The miss chance](#the-miss-chance).
- `--seed-mismatches auto|0|1|2` (default `auto`): how many mismatches a query k-mer may have
  against a seed. `1` also looks up every k-mer one substitution away, and `2` every k-mer two
  substitutions away. This finds assemblies at lower identity and costs more CPU. `auto` chooses
  per query, see [The miss chance](#the-miss-chance).
- `--deadline SECONDS`: stop the search of each query after this time. The summary then says
  `stopped`.
- `--server NODE_ID@ADDRESS`: search on a `jam serve` server instead of local index files.
- `--cache DIR`, `--memory SIZE`, `--disk SIZE`: the cache directory and the budgets for index
  files on S3, see [S3 and local caches](#s3-and-local-caches). `--memory` defaults to the
  budget of the plan in the cache directory, else 1GiB; `--disk` to the budget of the plan,
  else 0.

Example:

    jam search assemblies.jam plasmid.fa -c 0.8 -i 0.9 --format paf -o hits.paf

## Output

`jam search` writes only results to the standard output (or to `-o FILE`). Progress and debug
lines go to the standard error.

### JSON lines (the default)

One line per reported assembly, at the level of the search. Every line has these keys:

- `q`: the query name.
- `a`: the assembly (entry) name.
- `s`: the stage the assembly reached: 1 at `--level sketch`, 2 at `--level locate`, 3 when
  aligned and confirmed.

At `--level sketch`, one line per candidate:

    {"q":"plasmid1","a":"SAMN0001","s":1,"c":0.83}

- `c`: the containment estimate, the share of the query that the assembly seems to hold, 0 to 1,
  from shared seeds alone. Mutations lower it.

At `--level locate`, one line per located candidate:

    {"q":"plasmid1","a":"SAMN0001","s":2,"c":0.81,"w":[["contig_7",120000,148000]]}

- `c`: the containment estimate again, from the seeds whose locations were found.
- `w`: the alignment windows, each as `[contig, start, end]` on the contig.

At `--level align`, one line per confirmed assembly:

    {"q":"plasmid1","a":"SAMN0001","s":3,"cov":0.982,"id":0.9953,"bases":8243,
     "f":[["contig_7",0,8243,123456,131710,"+","4000=1X4242="]]}

- `cov`: the query coverage, the share of query bases that the alignments cover.
- `id`: the identity over the covered query bases.
- `bases`: the covered query bases.
- `f`: the fragments, best first. Each is `[contig, query start, query end, contig start,
  contig end, strand, CIGAR]`. Positions start at 0, and the end is the position after the
  last. On a circular query, a fragment over the end ends past the query length. The strand is
  `+` or `-`. The CIGAR uses `=`, `X`, `I` and `D` as in SAM, with the contig as the reference,
  in the order of the contig.
- `alt`: only when there are any: alternative fragments in the same form. They add no new query
  base, for example other copies of the same region.

Shares are rounded to 4 decimals.

After the result lines of a query comes its summary line:

    {"q":"plasmid1","status":"complete","n":[120,14,9,3],"miss":0.0001,"t":[0.21,0.35,1.8]}

- `q`: the query name.
- `status`: how the search ended, see below.
- `reason`: only when the status is not `complete`: why.
- `error`: only when the status is `failed`: the error message.
- `n`: the assemblies at each step: `[seeded, possible, located, confirmed]`. Seeded assemblies
  share at least one seed with the query. Possible ones passed stage 1, located ones stage 2,
  and confirmed ones were reported.
- `miss`: the modeled chance that the search missed an assembly that qualifies, see
  [The miss chance](#the-miss-chance).
- `prune`: only when it is `false`: the pruning was off for this query, because the seeds
  alone could not stay within `--max-miss`.
- `t`: seconds of `[sketch, locate, align]`. Sketch is wall time, including preparing the
  query. Locate and align are summed over threads.

With several index files, the summary covers all of them: the counts are summed, `miss` is the
largest, and the status is the worst.

The status values:

| `status` | `reason` | Meaning |
| --- | --- | --- |
| `complete` | | Every candidate was decided. |
| `incomplete` | `resource_limit` | An assembly needed more memory or work than a limit allows, as when a periodic query meets many periodic regions. It was not decided. |
| `stopped` | `deadline` | The search passed `--deadline`. |
| `stopped` | `asked` | The search was asked to stop, for example when the reader of its output went away. |
| `failed` | `error` | A stage or the server failed. `error` says why. |

Whenever a query does not end `complete`, jam exits with an error after writing all lines. This
also holds for `--format tsv` and `paf`, which have no summary line.

### TSV

`--format tsv` writes a header line, then one line per confirmed assembly, with tabs between
the columns:

| Column | Meaning |
| --- | --- |
| `query` | The query name. |
| `assembly` | The assembly (entry) name. |
| `coverage` | The query coverage, as `cov` above. |
| `identity` | The identity, as `id` above. |
| `bases` | The covered query bases. |
| `fragments` | The number of fragments. |

### PAF

`--format paf` writes one line per fragment of a confirmed assembly. The query is the PAF query
and the contig is the PAF target. Columns 1 to 12 are standard: query name, query length, query
start, query end, strand, contig name, contig length, contig start, contig end, matches,
alignment columns, and mapping quality (always 255). Then come `tp:A:P` (`tp:A:S` for an
alternative fragment), `NM:i:` with the edits, and `cg:Z:` with the CIGAR as above.

### Progress and debug lines

`-v` writes one line per query to the standard error when its search ends, for example:

    plasmid1: complete; 120 seeded, 14 possible, 9 located, 3 confirmed; sketch 0.210 s, locate 0.350 s, align 1.800 s

A query without pruning ends with `; pruning off`.

`-vv` also writes debug lines to the standard error, as JSON. Each names its index file in
`shard` and its kind in `kind`:

- `possible`: stage 1 kept the assembly, with its coverage bound and containment estimate.
- `located`: stage 2 found alignment windows, with their number and bases.
- `dropped`: the assembly left the search, with the `stage` and the `reason`:
  `coverage_bound` (it cannot reach the coverage), `no_seed_hits`, `decision` (the alignments
  fail the decision rule) or `resource_limit`.
- `unfinished`: the search stopped before the assembly was decided.
- `summary`: the full summary of the index file, with the work of each stage and the times.

The debug lines are for finding out why an assembly is or is not reported. Their keys may
change between versions.

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
