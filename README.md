<!-- What jam does and how to use its commands. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam

jam (just another minimizer) searches and aligns query sequences (plasmids, genes, contigs; 1 kb
to a few Mb) in very large collections of assemblies. It reports full base-level alignments with
CIGAR strings and a decision for every entry.

## How it works

You build an index of your assemblies once. Each assembly (or each FASTA record) is an *entry*.
A search answers three questions, one per stage:

1. **Sketch**: which entries may contain the query? Each one gets a containment estimate, the
   share of the query it seems to hold, from shared seeds alone.
2. **Locate**: where in those entries do the shared seeds lie?
3. **Align**: how exactly does the query align there? The entry is reported when it passes the
   decision rule.

Each stage reads more data than the one before and only works on what the previous stage kept.
A search can stop after any stage (`--level sketch|locate|align`). The loss budget
(`--loss-budget`, default 0.0001) is the largest accepted chance of missing an entry that
qualifies.

## Install

jam is a Rust program (edition 2024, so a recent stable toolchain is needed). Install Rust with
[rustup](https://rustup.rs), then:

    cargo build --release --locked

The binary is `target/release/jam`. Run `jam --help` or `jam COMMAND --help` for all options.

The block cache uses a patched copy of the `foyer` crate. `Cargo.toml` takes it from a git
repository (`[patch.crates-io]`), so the first build needs network access to GitHub. For an
offline build, run `cargo vendor` on a machine with network access and build with
`--offline`.

The `Dockerfile` builds the same binary into a small container image.

## Quick start

Build an index from FASTA files (plain or compressed). Each file becomes one entry:

    jam index build -o plasmids.shard -t 8 plasmids/*.fa

The output directory must not exist yet. Use `-l LIST` for a file with one FASTA path per line,
or `--entry-per-record` to make each FASTA record its own entry.

Search it with the records of a query file:

    jam search plasmids.shard query.fa > hits.jsonl

Each query record is searched on its own, in every shard you name. Write to a file with
`-o FILE`.

## Output

`jam search` writes JSON lines, one per entry event, then one `summary` line per query. The field
`kind` names the line: `possible`, `located`, `confirmed`, `dropped`, `unfinished` or `summary`.

- `possible`: stage 1 kept the entry (with the containment estimate).
- `located`: stage 2 found alignment windows in the entry.
- `confirmed`: the entry passed the decision. The line holds the query coverage, the identity,
  the base counts and the fragments. Each fragment has the contig name, the query segments, the
  contig interval, the strand, the score, the identity and the CIGAR.
- `dropped`: the entry failed, with the stage and the reason.
- `unfinished`: the search stopped first (for example at the deadline).

Entries without any seed hit are counted in the summary, not listed.

The `ending` of the summary is `complete` only when every possible entry was decided. It is
`resource_incomplete` when an entry was dropped for a memory or work limit (reason
`resource_limit`), `deadline` or `asked` after a stop, and `failed` after an error. Whenever a
search does not end `complete`, jam exits with an error after writing all lines, also with
`--format paf`, which has no summary line.

`--format paf` writes one PAF line per fragment of a confirmed entry instead. The query is the
PAF query and the contig is the PAF target. Columns 1 to 12 are standard, followed by `tp:A:P`
(`S` for an alternative fragment), `NM:i:` and `cg:Z:`. The CIGAR uses `=`, `X`, `I` and `D`
as in SAM, with the contig as the reference.

## The decision rule

An entry is reported when its alignments cover enough of the query at a high enough identity.

- `-c, --coverage` (default 0.5): the least query coverage. Give a share such as `0.8` or a base
  count such as `800bp`.
- `-i, --identity` (default 0.95): the target identity over the covered query bases.
- `--fragment-identity` (default: the value of `--identity`): the least identity of a stretch
  of a fragment to count.
- `--fragment-columns` (default 40): the least number of alignment columns for a fragment to
  count.

The segment rule in plain words: only the stretches of an alignment that reach the identity
limit count toward coverage. A poor stretch inside a fragment adds nothing. Fragments are
chosen so that no query base is counted twice.

The contig rule is off by default. It keeps out contigs that only share a gene or a repeat with
the query:

- `--contig-coverage`: the least share of a contig's bases that must align (0 turns the rule off).
- `--contig-identity`: the least identity of the fragments of that contig.

Other options of `jam search`:

- `--circular`: treat the queries as circular. A fragment over the query end then has two query
  segments. The default is linear.
- `--seed-mismatches auto|0|1`: with `1`, a query k-mer may also match a seed with one
  substitution. This finds entries at lower identity and costs more CPU. `auto` (the default)
  chooses per query from the identity, the loss budget and the query length.
- `--level sketch|locate|align` (default `align`): the last stage that runs.
- `--deadline SECONDS`: stop the search of each query after this time. Entries not finished are
  reported as `unfinished`.
- `-t, --threads`: threads to use (default: all cores).

Example:

    jam search plasmids.shard query.fa --coverage 0.8 --identity 0.9 --format paf -o hits.paf

### Low-complexity queries

Low-complexity stretches (repeats of a few bases, as found by DUST masking) count toward no
coverage. A query that is mostly low-complexity cannot reach the required coverage with its
other bases. It gets no answer, rather than false hits.

## Indexes

A shard is one directory with these files:

- `manifest.json`: the parameters, the entries, and the size and SHA-256 digest of each file.
- `seeds`: the seeds of all entries. Its sketch part answers stage 1. Its locator part answers
  stage 2.
- `packed/`: the sequences of each entry, 2 bits per base. Stage 3 reads them.

`jam index build` takes `-k` (bases per k-mer, default 25) and `-w` (k-mers per minimizer
window, default 100). With `--sketch-only` it writes the sketch part alone and no packed
sequences. Such a shard answers stage 1 only (`--level sketch`) but is much smaller.

An index can have several shards. `jam index reshard` merges and splits shards into a directory
of new shards, balanced by seed locations. Use `--shards N` for a number of output shards, or
`--max-seed-locations N` to bound the size of each seeds file:

    jam index reshard -o all.index --shards 4 part1.shard part2.shard

The packed sequences are linked, not copied, unless you give `--copy-packed`. The input can also
be index directories of earlier reshards. `jam search` accepts several shards (or index
directories) before the query file.

`jam stats SHARD` shows what a shard holds (`--json` for JSON).

## S3 and local caches

A shard can lie on S3. Give `s3://bucket/prefix` instead of a directory to `jam search`, `jam
screen`, `jam serve`, `jam stats` and `jam index warm`. Credentials come from the usual `AWS_*`
environment variables, and the endpoint from `AWS_ENDPOINT_URL`. jam never sends HEAD requests:
it takes object sizes from the manifest.

Reading from S3 is slow for the first requests. Two caches help:

- **Pins.** `jam index warm` copies parts of a shard into a local cache directory. A search or
  a server reads them from there, from RAM or from local disk. Pins are never evicted.
- **Block cache.** Everything not pinned is read in blocks. The blocks stay in RAM, and on local
  disk when you give a disk budget. They are kept across runs, so a second search is faster.

Example:

    jam index warm s3://bucket/plasmids --cache /scratch/cache --memory 16GiB --disk 200GiB
    jam search s3://bucket/plasmids query.fa --cache /scratch/cache --memory 16GiB --disk 200GiB

`jam index warm` needs `--cache`, `--memory` and `--disk`. By default it pins the sketch part,
then the locator part, while the disk budget lasts. Choose parts yourself with `--pin`:

    --pin sketch,locators,packed
    --pin sketch,packed=NAME1,NAME2

Every item after `packed=` is an entry name, so it comes last. The flag may repeat.

`jam search` and `jam serve` take `--cache`, `--memory` and `--disk` too. `--memory` defaults
to the budget of the plan, else 1GiB. `--disk` defaults to the budget of the plan, else 0.

## Screening assemblies

`jam screen` works the other way round. The index holds many small sequences (plasmids, genes).
Large assemblies stream past it once:

    jam screen plasmids.shard --assemblies assemblies.txt --coverage 0.8 -o screen.jsonl

`assemblies.txt` has one path per line (plain or gzip FASTA). The file name names the assembly.
Every shard entry that an assembly holds is reported as a JSON line, followed by one `assembly`
summary line per assembly (`failed` for one that could not be read).

`jam screen` takes the same decision options as `jam search` (`--coverage`, `--identity`,
`--fragment-identity`, `--fragment-columns`, `--contig-coverage`, `--contig-identity`). More
options:

- `--closest-margin M` (0 to 1): a contig counts only for the entries within this margin of its
  best one. Without it, a contig counts for every entry it qualifies for.
- `--batch-bases N` (default 500000): bases of contigs searched together.
- `--contig-dir DIR`: write the contigs that share enough with an entry to this folder, one gzip
  FASTA per assembly.
- `--survey-only`: stop after the survey. Nothing is aligned and no entries are reported. Only
  `--contig-dir` is written.
- `-t, --threads`: each thread screens one assembly at a time.

## Server

`jam serve` keeps a shard open and answers searches from other machines. The machines connect
over [iroh](https://www.iroh.computer), with a secret key file for each server:

    jam serve plasmids.shard --key server.key

The key file is created on first start. Its public key is the node id. `jam serve` prints the
address as `NODE_ID@ADDRESS`. Search on the server with:

    jam search --server NODE_ID@ADDRESS query.fa

With `--server`, give only the query file. The decision options work as for a local search.

Options of `jam serve`:

- `--bind ADDRESS` (default `0.0.0.0:4100`): a UDP address to listen on. Repeat it to add an
  IPv6 address.
- `--searches N` (default 4): searches that run at once.
- `--queue N` (default 32): searches that wait for a place. More are refused as busy.
- `-t, --threads`: threads of each search.
- `--cache`, `--memory`, `--disk`: as above, for shards on S3.
- `--peer NODE_ID@ADDRESS`: a peer server. Repeat it for more peers.

With peers, the server is a coordinator. Every search also runs on each peer, and the coordinator
merges the results. This lets you spread a large index over several machines, one shard each.

## Speed and size

These numbers are from one test on 2026-10-22, one run per case, so treat them as rough. The
collection had 1,000 SPIRE assemblies (241.8 Gbp). Searches ran with 8 threads and a cold file
cache.

- Building the index took about 10 minutes with 32 threads. The index was 108 GB (68 GB packed
  sequences, 28 GB locators, 12 GB sketch).
- On local NVMe, queries of 1 kb to 10 kb finished in about 1 to 10 seconds. A query of 2 Mb
  with no hit took about 40 seconds.
- On network storage (Ceph), the same searches took about 30 seconds to 3 minutes.

## License

jam is released under the MIT license. See `LICENSE`.
