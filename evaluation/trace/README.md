# Trace evaluation harness

This directory contains the small deterministic harness used to exercise the
one-plasmid-to-many-metagenome trace workflow. It measures catalog and archive
construction separately from trace search and can run local resources, a
range-capable HTTP fixture, and optional locally installed BLASTn or minimap2
candidate-alignment stages.

Build the release binary before measuring:

```bash
cargo build --release --locked
evaluation/trace/run.sh --profile smoke --repeats 3 --threads 8
```

The `smoke`, `standard`, and `large` profiles are deterministic synthetic
inputs. An existing materialized dataset can be supplied without downloading
or modifying it:

```bash
evaluation/trace/run.sh \
  --dataset-root /path/to/dataset \
  --query-id plasmid_accession \
  --work-dir tools/out/trace-benchmark/example \
  --repeats 3 \
  --threads 8
```

The dataset root must contain `catalog.fasta` and `assemblies/`. Optional
`truth.tsv` and `dataset.json` files are copied into the work directory. Every
run writes its resolved commands, tool versions, host description, dataset
metadata, raw measurements, resource counters, outputs, and `summary.tsv` to a
new work directory. `summarize.py` regenerates the TSV from the measurement
JSON files.

`process-cold` means the first measured process had no harness warm-up. It does
not mean that the operating-system page cache was evicted. Compilation,
installation, environment creation, and input generation are outside search
timings.

These fixtures are correctness and measurement smoke tests. They do not
establish production accuracy, biological sensitivity or specificity, remote
object-store behavior, or a general performance advantage. Production
measurements require versioned real backgrounds, coordinate-level truth,
complete input checksums, storage details, pinned comparator versions, and
matched accuracy thresholds.
