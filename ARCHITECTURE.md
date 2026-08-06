<!-- The design of jam: terms, stages, formats, modules and code rules. -->
<!-- Copyright (c) 2023-2026 Sebastian Beyvers -->
<!-- SPDX-License-Identifier: MIT -->

# jam architecture

jam (just another minimizer) finds small sequences (genes, plasmids, contigs; 1 kb to a few Mb)
in very large sets of assemblies and aligns them, base by base, with CIGAR strings. One index
answers three questions, one per stage:

1. **Sketch**: which entries share enough seeds with the query to hold it, and roughly how much
   of it?
2. **Locate**: where in those entries do the shared seeds lie?
3. **Align**: how exactly does the query align there, and does the entry qualify?

A search can stop after any stage. Each stage reads more data than the one before, so each one
only works on what the previous stage kept.

## 1. Terms

Every module, type, option, comment and document uses these terms with exactly this meaning.

| Term | Meaning |
| --- | --- |
| entry | One unit of the index: one sequence (a plasmid) or one set of sequences (an assembly). Results are reported per entry. |
| contig | One sequence (one FASTA record) of an entry. |
| query | The sequence searched for. It is linear or circular. |
| k-mer | A substring of length k; jam always uses the canonical form (the smaller of the k-mer and its reverse complement). |
| seed | A k-mer of an entry that the index stores: the smallest k-mer (by hash) of each minimizer window, and the smallest of every clean run too short for a full window. |
| minimizer window | w consecutive k-mers of a contig; each holds at least one seed. |
| query k-mer | Any k-mer of the query. The query side uses all its k-mers; only the entry side selects seeds. |
| seed hit | A query k-mer that equals a seed of an entry. |
| sketch | The seeds of an entry. Stage 1 compares the query k-mers with the sketches of all entries. |
| color set | The set of entries that hold a seed. |
| containment estimate | The share of the query that an entry is estimated to hold, from its seed hits alone. |
| coverage bound | The highest query coverage an entry can still reach, given where its seed hits lie on the query. It is wrong at most with the chance set by the loss budget. |
| loss budget | The largest accepted chance of missing a qualifying entry (for example 0.01%). Every shortcut of a search stays within it. |
| candidate entry | An entry that stage 1 keeps. |
| bin | A block of 2^b bases (1,024 by default) of an entry's sequence. Positions in the positional index are bins. |
| positional index | For each seed, its seed locations in each entry. |
| seed location | The bin and strand where a seed occurs in an entry; one is stored per occurrence. |
| fingerprint | The few bits of a seed's hash that the index stores to tell seeds apart; two seeds can share one (a false seed hit). |
| candidate bin | A bin with seed hits of a candidate entry. |
| alignment window | The stretch of a contig that stage 3 reads around a candidate bin. |
| packed sequence | The contigs of an entry, 2 bits per base, with a table of contig names, lengths and runs of N. |
| anchor | An exact match of length a (15 by default) between the query and an alignment window. |
| chain | Anchors in order along nearly the same diagonal. |
| walk | An ungapped copy of the query in an alignment window, started where a k-mer too frequent or too simple for anchors matches, and extended along its diagonal. |
| join | An alignment that links short walks in query and contig order across the pieces between them, kept when together they reach the fragment limits. |
| piece | The query and contig bases between two linked walks of a join, aligned globally. |
| alignment | A gapped, base-level alignment of a query interval to a contig interval, with score, CIGAR and counts. CIGAR follows SAM with the contig as the reference: `I` marks bases only in the query, `D` bases only in the contig. |
| fragment | An alignment placed in query and contig coordinates; a fragment over the end of a circular query has two query segments. |
| mosaic | The fragments of one entry chosen to cover the query without counting a query base twice. |
| query coverage | The share of query bases covered by an entry's mosaic. |
| identity | Matches divided by alignment columns, over the fragments of the mosaic together. |
| decision | The rule that reports an entry: query coverage, identity, and the optional contig rule (a contig counts only when a set share of its bases aligns to the query at a set identity, so a contig that only shares a gene or a repeat is left out). |
| index | All shards of one collection: a directory of shards with `index.json`, or a single shard. |
| shard | A part of an index (a set of entries) that one process opens. An index has one or more shards. |
| frequent seed | A seed whose color set holds more entries than a limit; stage 1 does not expand it per entry. |
| coordinator | The server that sends a query to every shard and merges their results. |
| screen | The reverse direction: the small sequences form the index and large assemblies stream past it once. Results use the same decision. |
| manifest | The JSON file of a shard: parameters, entries, files and their SHA-256 digests. |
| block cache | Fixed-size blocks of packed sequences kept in RAM (and optionally on local disk) between reads. |

Words not used: sample (say entry), catalog (say index), locator (say positional index), trace,
a bare "window" (say minimizer window or alignment window), a bare "hit" (say seed hit), and
"sketch" in the MinHash sense.

All hashing uses jamhash (`jamhash_u64` of the 2-bit canonical k-mer, k at most 32).

## 2. One search, step by step

```
query
  -> query k-mers, hashed
  -> stage 1, sketch  (RAM)         seed hits per entry -> containment estimate, coverage bound
                                    -> possible entries
  -> stage 2, locate  (local disk   seed locations of the possible entries: false seed hits
                       or S3)       drop out, the coverage bound is computed again
                                    -> alignment windows
  -> stage 3, align   (S3 or local  packed sequences of the windows -> anchors -> chains
                       disk)        -> alignments (CIGAR) -> fragments -> mosaic -> decision
  -> entry events, sent through one channel as they happen
```

Every entry that stage 1 keeps gets its events in order: possible (with the containment estimate
and the coverage bound), located (with its alignment windows), and at the end exactly one of
confirmed (with the decision and its alignments), dropped (with the stage and the reason) or
unfinished (the search stopped first). Entries without any seed hit are counted, not reported.
`--level sketch|locate|align` ends a search after that stage.

The coverage bound: a stretch of the query without seed hits may still hide one present piece of
up to `a` bases, where `a` is the shortest length whose miss chance, summed over all such
stretches, fits half of the loss budget. Stage 2 computes the bound again from the confirmed seed
locations with the other half. An entry whose bound is below the requested query coverage is
dropped before any sequence is read.

Frequent seeds: a seed whose color set holds more entries than a limit is not expanded per entry;
its query positions count as possibly hit for the entries of its color set. This can only raise
their bounds, so no entry is lost by it.

The containment estimate is one-way: the share of the query that the entry holds, estimated from
the seed hits alone. It counts seeds that survived unchanged, so mutations lower it.

## 3. Index on disk

A shard is one directory. All files are laid out for memory mapping: fixed-width little-endian
fields, aligned arrays, offsets instead of pointers, a version number and a checksum per section.

| File | Stage | Holds | Read |
| --- | --- | --- | --- |
| `manifest.json` | all | parameters (k, w, bin size, bucket bits, fingerprint bits), entries (name, contigs, bases), file sizes and SHA-256 digests | at open |
| `sketch` | 1 | a directory of hash buckets; per bucket the sorted fingerprints of its seeds and the color set id of each; then every distinct color set once, entries as sorted ids | memory-mapped and kept in RAM |
| `positions` | 2 | per entry a section with its contig starts and its seed locations, grouped by seed hash, so that one candidate entry needs few reads | memory-mapped from local disk, or byte ranges from S3 |
| `packed/` | 3 | one packed sequence file per entry (or the same files in another directory or under an S3 prefix) | byte ranges, through the block cache |

`FORMATS.md` gives the byte layout of every file.

An index is a directory of shards with `index.json`, which lists the shards in order with their
entry ranges. Resharding (`jam index reshard`) turns any shards into any number of new shards,
balanced by seed locations: every entry keeps its positions section and packed sequence, and each
new shard rebuilds only its sketch. Merging and splitting are the two directions of it.

## 4. Source layout

A module uses only modules listed above it.

| Module | Job |
| --- | --- |
| `seeds` | canonical k-mers, jamhash, minimizer windows, the low-complexity filter |
| `loss` | the chance that a present query stretch shows no seed hit; the loss budget |
| `storage` | reading bytes by offset: local files, S3 byte ranges, the block cache |
| `index` | the on-disk index: manifest, packed sequences, positional index, building and joining shards |
| `align` | alignment types, CIGAR, the banded alignment kernel, anchors and chains |
| `decide` | fragments, mosaic, query coverage, identity and the decision |
| `search` | the three stages for one query against one shard, with memory limits and streaming |
| `report` | the output records of a search: JSON lines for every entry event and the summary, and PAF lines with the CIGAR in the `cg:Z` tag; the command line and the network use the same records |
| `screen` | the reverse direction over streamed assemblies |
| `service` | `serve` and remote search over TCP: protocol, server, client, coordinator |
| `cli` | the commands: `index`, `search`, `screen`, `serve`, `stats` |

## 5. Code rules

- Use the terms of section 1 exactly. No abbreviations beyond k, w, CIGAR, S3, TCP and RAM.
- One job per module. Files stay under about 800 lines and functions under about 80; split by
  step, not by helper type.
- Every public item has a doc comment of one or two plain sentences.
- One error type per module; bad input never panics; `unwrap` only in tests. No anyhow. Errors
  travel with `?` and are never stored as values; work on several threads sends its results
  through one channel, and the receiving thread handles errors.
- No wrappers that only rename or forward; no builders, traits or generics with one use.
- No dead code, no experiment flags, no counters that nothing reads.
- Tests sit next to the code they test and share one test-support module; random test data
  comes from the `rand` crate with fixed seeds. Every file format has a round-trip test; the
  search has end-to-end tests on small generated entries.
- Comments describe the code as it is, never its history.
- Moves and renames are separate commits from changes of behaviour or format.

## 6. Plan

1. **Cleanup**: the Aruna style policy, the review findings, and established libraries instead of
   hand-written infrastructure (OpenDAL for S3, foyer for the block cache, iroh for the service).
2. **Across machines**: the coordinator spreads stage 3 work over idle shards (every shard can
   read every packed sequence from S3).
3. **Validation**: accuracy and speed against jam 0.10, LexicMap and BLAST on the same test sets.
4. **Speed**: faster than LexicMap on every query length, on local disk and on network storage.
5. **PlasmidTracker** moves to jam 0.11.
