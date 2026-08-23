# `jam trace` algorithm

## Scope and identifiers

`jam trace` is a query-centered fragment search workflow. It accepts one
nucleotide query sequence and searches many metagenome resources for
positional sequence evidence. The query can be a plasmid, phage, synthetic
construct, other mobile genetic element, or an unknown sequence class.

The workflow reports sequence evidence. It is not a plasmid or phage
classifier, assembler, binning method, read mapper, or proof that a query
element exists as an autonomous biological molecule.

The source exposes these layer identifiers:

```text
screen_algorithm          = jam-fracminhash-screen-v1
local_alignment_algorithm = jam-exact-seed-chain-banded-v1
mosaic_algorithm          = jam-fragment-mosaic-v1
trace_workflow            = jam-trace-v1
```

The trace metadata and JMA compatibility check retain the version-1 workflow
identifier `jam-seed-chain-align-v1`. This is the archive/run compatibility
identifier; the four layer identifiers describe the individual stages. The
workflow version is `1`.

The implementation is distinct from BLAST and minimap2. That distinction is
descriptive and does not imply a novelty, accuracy, or general performance
claim.

## Coordinate conventions

All intervals are zero-based and half-open: `[start, end)`. Query coordinates
refer to the stored query sequence. Target coordinates refer to the forward
stored contig, including when the accepted alignment is on the reverse
strand.

The result separates coordinate handling from biological topology:

```text
query_kind          = plasmid | phage | other | unknown
topology_requested  = linear | circular | auto | unknown
coordinate_model    = linear | wrap | undetermined
topology_evidence   = linear_supported | wrap_supported |
                      both_compatible | insufficient | undetermined
```

`coordinate_model=wrap` means that query coordinates were allowed to cross the
stored origin. It does not mean that the molecule is biologically circular.
Separate contigs remain separate observations; their order in a mosaic does
not prove physical linkage.

## Pipeline

The implementation follows this bounded path:

```text
read exactly one query FASTA record
build a candidate-query sketch
search the metagenome .jam candidate index
apply independent candidate containment and shared-hash limits
rank candidates deterministically and apply the candidate limit

for each selected metagenome:
    open its indexed JMA resource, or use raw FASTA as an explicit fallback
    extract the requested positional query seeds
    request exact seed occurrences and verify packed canonical k-mers
    suppress over-repetitive occurrence groups
    build bounded chains per contig, strand, and k-mer level
    retrieve only sequence windows implicated by accepted chains
    run fixed-band local affine-gap alignment
    optionally run gap-directed dense k=31 and conditional k=21 rescue
    retain accepted alignments and their evidence
    build linear and/or wrap-aware fragment mosaics
    report coverage, gaps, alternatives, warnings, and resource counters
```

The default output contains one run header, one result per selected
metagenome, and one run footer in JSONL or JSONL.zst form. All accepted local
alignments remain in the corresponding result record, while the mosaic refers
to them by deterministic alignment IDs.

## Candidate sketch stage

Candidate retrieval uses canonical k-mer sketches and the existing `.jam`
index. The query sketch examines every canonical k-mer accepted by the FASTA
parser; it does not use minimizers. Hash identity is `jamhash_u64_v1`, and
hash zero is excluded consistently from query and index sketches.

The two containment directions have separate denominators:

```text
query_containment      = shared_hashes / query_hashes
metagenome_containment = shared_hashes / metagenome_hashes
```

Candidate filtering is sketch evidence only. A candidate is not a final
presence call, and a sketch hit is not an alignment or a topology result.

## Positional seed selection

The positional stage also examines every canonical k-mer window. For k-mer
size `k` and FracMinHash scale `s`, let:

```text
h = jamhash_u64_v1(canonical_packed_kmer)
```

The retained-seed rule is:

```text
h != 0 && h < floor(u64::MAX / s)
```

Hash zero is counted as a diagnostic skip but never emitted. Smaller scale
values retain denser nested sets. Sequence ranges are normalized for case and
U bases, and windows containing non-ACGT symbols are skipped by the k-mer
iterator.

Only two seed sizes are accepted:

```text
k=31  primary layer
k=21  rescue layer
```

No seed size below 21 is used. `max_occurrences` is not an extraction density
parameter; it is applied when exact occurrences are converted into anchors.

Gap rescue supplies unresolved half-open query intervals. Each interval is
expanded by the configured flank, clipped to the query, and merged with
touching or overlapping ranges before hashing. A complete k-mer window must
lie inside the merged range. Seed positions remain global query coordinates.

The JMA lookup key is the complete tuple:

```text
(k, hash, canonical_packed_kmer)
```

A hash match without packed canonical-k-mer equality is only a lookup
candidate. It cannot become an anchor or occurrence evidence. This explicit
verification prevents a `jamhash_u64_v1` collision between unequal k-mers from
creating alignment evidence.

Seed extraction is linear in the number of bases in the requested ranges,
apart from sorting and merging the requested intervals. Retained seed memory
is proportional to the number of emitted seeds.

## Occurrence evidence and anchors

An anchor joins a query seed position to one exact JMA occurrence. Relative
orientation is the exclusive-or of the query and target canonicalization
orientation flags. Occurrence evidence records query occurrence count,
candidate occurrence count/group count, collection document frequency when
available, repetitive state, and common-sequence state.

Collection document frequency is currently unavailable in the occurrence
lookup path and is reported as unknown. A common-sequence flag is based on the
configured candidate occurrence threshold; it is evidence annotation, not a
classifier.

If a seed occurrence group has more than the profile's
`max_occurrences_per_seed`, the complete group is marked repetitive and
contributes no anchors. Hash-zero groups are discarded defensively. Remaining
anchors are sorted deterministically, duplicate anchors are removed, and the
global anchor limit is applied.

Repetitive or common sequence is not silently treated as absent. It remains
available for reporting, but support confined to such sequence is labeled and
requires independent interpretation.

## Bounded anchor chaining

Anchors are grouped by contig, strand, and k-mer size. k=31 and k=21 anchors
are never mixed in one chain.

For a linear coordinate model, each source anchor appears once. For a wrap
coordinate model, each unused source anchor is also represented once at
`query_position + query_length`. A source-index guard prevents the two copies
from being used together in one path.

For each expanded anchor `i`, the chain dynamic program is:

```text
score[i] = 100
count[i] = 1

for at most P preceding anchors j:
    reject j unless:
        query position increases
        j is from a different source anchor
        query gap <= max_query_gap
        target order is monotone for the strand
        target gap is positive and <= max_target_gap
        score[j] is finite

    gap = abs(query_gap - target_gap)
    extension = 100 - gap_penalty * gap
    candidate = score[j] + extension
    retain the best candidate for i
```

All shipped profiles use `gap_penalty=1`; chain gaps are not affine. The
predecessor scan counts rejected candidates toward `P`, which bounds
repeat-heavy groups. Ties prefer higher score, then more anchors, then the
earlier query-coordinate predecessor. Endpoints require the configured
minimum anchor count. Endpoint ties use score, count, and deterministic query,
target, and hash keys.

The path is traced back, and a path spanning more than one query length is
rejected. A wrap path is normalized into one or two ordinary query segments;
`origin_crossing` is true only when two segments result. Linear and
undetermined chains do not duplicate query anchors and cannot emit a crossing
path.

For `m` anchors in one group and predecessor cap `P`, one bounded chain pass is
`O(mP)` after sorting. Extracting up to `C` chains gives
`O(m log m + C m P)` per group. Wrap expansion changes the constant by at most
about two. Anchor, chain, and global output limits bound retained state and can
affect sensitivity for repetitive or highly fragmented candidates. Seed-
bucket decoding and temporary anchor construction can peak before those
retention limits are applied.

## Fixed-band affine local alignment

The local alignment kernel is a three-state affine-gap dynamic program. It is
not global, semiglobal, overlap, BLAST, or minimap2 alignment. Non-positive
scores reset to zero; traceback starts at the highest positive cell and stops
at a restart cell.

For query base `q` and target base `t`, with
`goe = gap_open + gap_extend`:

```text
M(i,j) = max(0, best(M,I,D at i-1,j-1) + substitution(q,t))

I(i,j) = max(0,
             I(i,j-1) + gap_extend,
             M(i,j-1) + goe,
             D(i,j-1) + goe)

D(i,j) = max(0,
             D(i-1,j) + gap_extend,
             M(i-1,j) + goe,
             I(i-1,j) + goe)
```

The first base of a gap therefore costs `gap_open + gap_extend`. Base
comparison is ASCII-case-insensitive; ambiguity symbols are not wildcards.
The shipped scoring constants are:

```text
match        +2
mismatch     -3
gap open     -5
gap extend   -1
```

The matrix is restricted to a fixed half-band around the supplied diagonal:

```text
target columns for query row i:
[max(0, i + diagonal_offset - band_width),
 min(target_length, i + diagonal_offset + band_width)]
```

The band is clipped at target bounds and never widens. Version 1 has no
X-drop. A path outside the band returns a structured traceback/band error; it
is not retried with a wider band automatically.

The best-cell tie break is deterministic: smaller query index, then target
index, then state. Traceback emits merged edit runs and a canonical CIGAR with
only `=`, `X`, `I`, and `D`. Unaligned local prefixes and suffixes are reported
by query and target intervals rather than soft clips. The edit script and CIGAR
are sufficient to reconstruct the aligned strings.

Reverse-strand alignment reverse-complements the target window for dynamic
programming and maps its selected target interval back to forward contig
coordinates. Circular query alignment materializes at most one query-length
linearized span from an explicit `query_start`, then maps the selected span
back to one or two normalized query segments.

If `Q` is the query-window length and `B` is the half-band width, the dynamic
program uses `O((Q+1)(2B+1))` cells and the same asymptotic time. The workspace
also stores row metadata and traceback operations. `DEFAULT_MAX_CELLS` is
4,000,000 cells; runner windows set a profile-derived bound. Matrix limits,
empty inputs, a band that excludes the input, no positive alignment, and
traceback inconsistencies are explicit errors.

## Fragment mosaic and coverage

The unit of final evidence is one metagenome. Each accepted local alignment
retains its contig, query segments, target interval, strand, score, identity,
CIGAR/edit script, anchor evidence, and mosaic role.

The mosaic splits query coordinates at every accepted alignment boundary. For
each atomic interval, all covering alignments are compared deterministically:

```text
alignment score, descending
identity, descending
score per query base, descending
nonrepetitive anchor count, descending
supported interval length, descending
contig identifier, ascending
alignment identifier/canonical key, ascending
```

The first alignment becomes primary support for that interval; the remaining
alignments remain alternatives. The mosaic reports selection components
rather than hiding them in one presence score.

Coverage projection has separate meanings:

```text
base coverage       query bases paired through CIGAR = or X
aligned span        query bases consumed by =, X, or D
unsupported gap     query interval with no accepted base support
alignment deletion  query interval consumed by D inside an alignment
```

`I` consumes target only. Therefore unsupported gaps and alignment deletions
are reported separately. Supported intervals are unioned, so overlapping or
duplicate contigs cannot inflate coverage beyond the query length. Supporting
contigs are listed but are not physically joined.

The result also reports common-sequence intervals, repeat-only intervals,
nonrepetitive-supported bases, alternative alignment counts, atomic interval
assignments, and all accepted alignment records. An alignment supported only
by repetitive sequence is labeled `repeat_only`; common mobile sequence is
labeled separately when its evidence meets the common-sequence rule.

## Coordinate models and topology requests

### Linear

Linear chaining uses one copy of each query anchor. Linear mosaics use ordinary
interval unions and complements. Origin-crossing alignments are excluded from
the linear model, so terminal gaps remain separate. The runner uses a full
query span starting at zero for linear window alignment; the resulting absolute
query interval cannot cross the origin.

### Circular

Circular requests select the wrap coordinate model. Chaining permits one
query-length duplicate, and alignment can use a linearized span of at most one
query length. Final intervals are normalized into `[0, query_length)` and a
crossing alignment is represented by two ordinary segments. Circular gap
complements preserve the two terminal pieces and the result records a
wrapping gap when appropriate.

This is coordinate handling only. It does not establish biological circularity.

### Auto

Auto builds the main accepted local alignment set once, then summarizes both a
linear and a wrap model from that same set. It selects wrap only when wrap
support exceeds linear support by strictly more than
`auto_topology_margin_bases`, selects linear under the converse, reports
`both_compatible` when positive support is equal, and otherwise reports
`undetermined`.

The current selection comparison uses newly supported query bases only. The
model records also contain alignment quality, nonrepetitive anchor support,
origin-crossing count, terminal gaps, contradictions, and fragment count, but
those components do not currently affect the Auto decision. The
implementation does not perform a full second candidate search and does not
claim to infer biological topology.

The runner uses wrap-aware chaining for Auto, so origin-crossing candidates can
be present in the shared accepted alignment set before the two summaries are
constructed. This is an optimization for coordinate comparison, not a
topology call.

### Unknown

Unknown retains a linear display model as the primary coordinate result and
reports `coordinate_model=undetermined` and
`topology_evidence=undetermined`. A wrap model may be retained as an
alternative summary, but terminal gaps are not merged into the primary
undetermined display.

## Gap-directed rescue

Indexed candidates can use bounded rescue rounds. The initial round hashes the
complete query with primary k=31 seeds. Later rounds use only unresolved query
gaps, expanded by profile flanks:

```text
round 1: primary k=31
round 2: denser k=31 inside eligible gaps, when configured
round 3: conditional k=21 rescue inside still-eligible gaps, when configured
```

Each rescue round requests only implicated seed buckets, builds new chains,
fetches only newly implicated sequence blocks, and aligns only new windows.
Rounds stop when no eligible gap remains, the round limit is reached, or the
per-round bucket/block/window limits are reached. Previously requested seed
keys and sequence windows are not repeated.

For Auto, terminal gaps are eligible for later rounds only after nonrepetitive
seed evidence is found near both query ends on the same contig and strand.
This is a retrieval guard, not evidence of a circular molecule. Raw FASTA
fallback does not provide indexed gap-directed rescue and reports that limit.

Rescue metrics include target gaps, seed buckets, seed keys tested, anchors,
chains, sequence blocks, alignment windows, newly supported query bases, and
elapsed time per round.

## Sensitivity profiles

The profile names are execution settings, not calibrated probabilities:

| parameter | fast | balanced | sensitive |
| --- | ---: | ---: | ---: |
| primary k / scale | 31 / 500 | 31 / 200 | 31 / 100 |
| rescue k / scale | none | 21 / 500 | 21 / 200 |
| primary max occurrences | 64 | 128 | 256 |
| max candidates | 25 | 100 | 250 |
| max anchors/candidate | 25,000 | 100,000 | 250,000 |
| max chains/candidate | 4 | 8 | 16 |
| minimum chain anchors | 3 | 3 | 2 |
| maximum alignment window | 250,000 | 1,000,000 | 4,000,000 |
| alignment half-band | 64 | 128 | 256 |
| max candidate concurrency | 4 | 4 | 2 |
| gap-rescue rounds | 1 | 3 | 3 |
| dense k=31 scale | none | 100 | 100 |
| minimum rescue gap | 250 | 200 | 100 |
| rescue flank | 64 | 96 | 128 |
| common occurrence threshold | 8 | 8 | 16 |
| Auto comparison margin | 250 | 200 | 100 |

Balanced and sensitive profiles are the settings that can use conditional
k=21 rescue. Fast does not run rescue rounds.

## Known limits and failure boundaries

- A genuine trace can be missed when no retained seed is available, a seed is
  suppressed as repetitive, the minimum anchor count is not reached, or an
  anchor/chain/window limit is reached.
- A query or target path outside the fixed alignment band can fail. There is
  no automatic band widening or X-drop recovery in this version.
- Local alignment reports the best positive region in a retrieved window. It
  does not explain every separated homologous region without separate chains
  and windows.
- The bounded output path retains at most `max_alignments_per_candidate`
  alignments. Rescue rounds keep already retained alignments and do not remove
  earlier accepted evidence; once the cap is full, additional alignments are
  omitted rather than expanding the result. Evidence beyond that cap can
  therefore be missed, especially for repeat-rich candidates.
- Common mobile elements, resistance cassettes, repeats, chromosomal
  integrations, and contamination can produce genuine sequence evidence for
  multiple query elements. Common or repeat-only support is not proof of
  element identity or autonomy.
- Separate supporting contigs are not assumed to be connected, and coordinate
  order does not establish a molecule or assembly graph path.
- Wrapped coordinates and `topology_evidence` describe compatible coordinate
  representations, not biological topology. Independent mapping, graph,
  long-read, or contextual evidence is required for a biological conclusion.
- A JMA archive whose available seed level is sparser than the requested
  profile returns a structured `seed_level_mismatch` failure. It does not fall
  back to unanchored alignment, because that would change the evidence model.
- Raw FASTA fallback bounds each alignment-window size but not the number of
  parsed contigs/windows, retains the parsed candidate assembly in memory, and
  does not provide indexed gap-directed rescue. Remote resources can
  additionally fail through range,
  cache, checksum, timeout, or authorization errors; these are reported as
  structured resource outcomes.
- Candidate screening currently collects matching sketch rows before applying
  the deterministic top-candidate limit, and anchor generation sorts a
  temporary collection before truncation. `--memory-target` is therefore a
  cache target, not a hard RSS ceiling.

`jam trace` therefore reports candidate and fragment evidence against a query
sequence. Final interpretation remains an independent biological analysis.
