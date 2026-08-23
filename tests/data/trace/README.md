# Trace workflow truth fixtures

This directory contains deterministic, artificial inputs for the trace workflow.
The sequences are deliberately generated from fixed integer seeds; they are not
biological reference sequences and must not be used to estimate sensitivity on
real plasmids. `truth.json` is the source of expected relationships. Tests must
not infer expected outcomes from a Jam result.

The catalogue is in `catalog.fa`. `assembly_trace.fa` contains positive,
mutated, split, overlapping, ambiguous, repeat-only, and unrelated contigs.
`assembly_empty.fa`, `assembly_negative.fa`, and `assembly_malformed.fa` cover
empty, no-hit, and parser-error paths respectively.

Regenerate the files into a directory with:

```text
rustc tools/generate_trace_fixtures.rs -o /tmp/generate-trace-fixtures
/tmp/generate-trace-fixtures tests/data/trace
```

The generator uses only the Rust standard library. Its output is checked in so
integration tests can use the fixtures without running a generator first.
