<!-- What jam does and how to use its commands. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam

jam (just another minimizer) finds small sequences (genes, plasmids, contigs; 1 kb to a few Mb)
in very large sets of assemblies and aligns them, base by base, with CIGAR strings.

## How it works

You build an index of your assemblies once. Each assembly (or each FASTA record) is an *entry*.
A search then answers three questions, one per stage:

1. **Sketch**: which entries share enough seeds with the query to hold it?
2. **Locate**: where in those entries do the shared seeds lie?
3. **Align**: how exactly does the query align there, and does the entry qualify?

Each stage reads more data than the one before and only works on what the previous stage kept.
A search can stop after any stage (`--level sketch|locate|align`). An entry is reported when its
alignments reach the requested query coverage and identity. The loss budget (`--loss-budget`) is
the largest accepted chance of missing a qualifying entry.

## Commands

Build a shard from FASTA files (one entry per file; `--entry-per-record` makes one per record):

    jam index build -o plasmids.shard -t 8 plasmids/*.fa

Merge or split shards into an index of new shards, balanced by seed locations:

    jam index reshard -o all.index --shards 4 part1.shard part2.shard

Search every record of a query file in one or more shards (linear by default, `--circular` for
circular queries):

    jam search plasmids.shard query.fa --coverage 0.8 --identity 0.95

Show what a shard holds (`--json` for JSON):

    jam stats plasmids.shard

`jam screen` streams large assemblies past an index of small sequences, and `jam serve` serves a
shard to other machines; both are being finished.

## Output

`jam search` writes JSON lines, one per entry event, then one `summary` line per query. The field
`kind` names the record: `possible`, `located`, `confirmed`, `dropped`, `unfinished` or `summary`.
A `confirmed` line holds the entry's query coverage, identity, base counts and its fragments,
each with contig name, query segments, contig interval, strand, score, identity and CIGAR.

`--format paf` writes one PAF line per fragment of a confirmed entry instead. The query is the
PAF query and the contig is the PAF target. Columns 1 to 12 are standard, followed by `tp:A:P`
(`S` for an alternative fragment), `NM:i:` and `cg:Z:`. The CIGAR uses `=`, `X`, `I` and `D`
as in SAM, with the contig as the reference. Use `-o FILE` to write to a file.

## Build

    cargo build --release

The binary is `target/release/jam`. Run `jam --help` for all options.
