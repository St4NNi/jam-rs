[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)

# jam-rs

Just another minhash (jam). A high-performance FracMinHash implementation for genomic sequence similarity analysis, optimized for searching plasmids, phages, and other small genomic elements in large datasets.

jam uses a custom hash function ([jamhash](https://github.com/St4NNi/jamhash)) that provides lower collision rates, 2-10x higher speed and better uniformity than murmur3. It also includes a compact memory-mapped database format (`.jam`) for fast random access, and a bias filtering system based on Count-Min Sketches to selectively increase sensitivity for target sequences.

### Installation

From [crates.io](https://crates.io/crates/jam-rs):

```bash
cargo install jam-rs
```

From source:

```bash
cargo install --git https://github.com/St4NNi/jam-rs
```

### Key Features

- **Custom hash function**: [jamhash](https://github.com/St4NNi/jamhash) provides lower collisions, better uniformity and is faster compared to murmur3
- **Bias-aware sketching**: Count-Min Sketch based compositional filtering with automatic background extraction
- **Complexity filtering**: Shannon entropy threshold to exclude low-complexity k-mers
- **Memory-efficient**: External sorting for processing datasets larger than available RAM
- **Compact storage**: 256-bucket memory-mapped `.jam` format with binary fuse filters for fast random access
- **Parallel execution**: File-level parallelization via rayon with configurable thread count
- **Tuned for speed**: jemalloc allocator, LTO, single codegen unit, `opt-level = 3`

### Usage

```console
$ jam --help
Just another (genomic) minhasher (jam), obviously blazingly fast

Usage: jam [OPTIONS] <COMMAND>

Commands:
  sketch  Sketch one or more files and write the result to an output file
  dist    Estimate containment of a query sequence against a sketch database
  bias    Build and analyze hash bias tables for filtering
  stats   Display statistics about a JAM database
  help    Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -s, --silent             Silent mode, no (additional) output to stdout
  -m, --memory <MEMORY>    Maximum memory usage in GB [default: 2]
  -h, --help               Print help
  -V, --version            Print version
```

#### Sketching

Create `.jam` databases from FASTA/FASTQ files (plain or gzip/bzip2/xz/zstd compressed). Supports single files, multiple files, or directories.

```console
$ jam sketch --help
Sketch one or more files and write the result to an output file

Usage: jam sketch [OPTIONS] --output <OUTPUT> [INPUT]...

Arguments:
  [INPUT]...  Input file(s), directories, or file with list of files to be hashed

Options:
  -o, --output <OUTPUT>          Output file (.jam format)
  -k, --kmer-size <KMER_SIZE>    K-mer size, all sketches must have the same size to be compared and below 32 [default: 21]
      --fscale <FSCALE>          Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash)
      --complexity <COMPLEXITY>   Complexity cut-off, only hash sequences with complexity above this value [default: 0.0]
      --singleton                Create a separate sketch for each sequence record
      --temp-dir <TEMP_DIR>      Custom temporary directory for intermediate files during sorting
      --bias-table <BIAS_TABLE>  Path to a bias table file (.bias) for compositional filtering
  -h, --help                     Print help
```

Examples:
```bash
# Sketch a single file
jam sketch input.fasta -o sketch.jam

# Sketch a directory with 8 threads and FracMinHash scaling
jam sketch genomes/ -o db.jam --fscale 1000 -t 8

# Filter low-complexity k-mers by Shannon entropy
jam sketch genomes/ -o db.jam --fscale 1000 --complexity 1.5

# One sketch per sequence record
jam sketch multi.fasta -o db.jam --singleton

# Apply bias filtering during sketching
jam sketch plasmids/ -o filtered.jam --bias-table host_filter.bias
```

#### Querying

Estimate containment of query sequences against a sketch database.

```console
$ jam dist --help
Estimate containment of a query sequence against a sketch database

Usage: jam dist [OPTIONS] --input <INPUT> --database <DATABASE>

Options:
  -i, --input <INPUT>        Input FASTA/FASTQ file to query
  -d, --database <DATABASE>  Database sketch (.jam file)
  -o, --output <OUTPUT>      Output to file instead of stdout
  -c, --cutoff <CUTOFF>      Cut-off value for similarity/containment [default: 0.0]
      --singleton             Singleton mode, process each query sequence separately
  -h, --help                 Print help
```

Examples:
```bash
# Query against a database with a containment cutoff
jam dist -i query.fasta -d db.jam -c 0.1 -o results.tsv

# Per-sequence queries
jam dist -i multi_query.fasta -d db.jam --singleton -c 0.1
```

Output is tab-separated: `query`, `sample_id`, `hit_count`, `containment`.

#### Bias Table Construction

Bias tables allow compositional filtering to increase sensitivity for target sequences while suppressing background noise. They work by scoring k-mers based on their enrichment in a positive (target) set relative to a negative (background) set.

The underlying data structure is a **Count-Min Sketch (CMS)**, a probabilistic structure that approximates k-mer frequencies using multiple independent hash functions mapped to a fixed-width table. This keeps memory usage constant regardless of the number of distinct k-mers. By default, the CMS uses 1,048,576 columns and 5 hash functions (~5 MB).

**How it works:**

1. K-mer frequencies from both the positive and negative input sets are counted into separate CMS tables.
2. **Background extraction**: The positive counts are subtracted from the negative counts (floored at zero). This prevents k-mers naturally shared between target and background from being penalized.
3. A log-ratio weight is computed per CMS cell: `log((pos + alpha) / (adjusted_neg + alpha))`, where `alpha` is a smoothing parameter.
4. Weights are quantized to `i8` (-127 to +127) for compact storage.
5. **Threshold calibration**: All 255 possible thresholds are evaluated. The threshold that maximizes fold enrichment (positive retention / negative retention) is selected. If a target fold enrichment is specified, the closest achievable threshold is used instead.

```console
$ jam bias create --help
Create a bias table from positive (target) and negative (background) FASTA files.
Target signal is always subtracted from background before computing bias weights.

Usage: jam bias create [OPTIONS] --positive <POSITIVE> --negative <NEGATIVE> --output <OUTPUT>

Options:
      --positive <POSITIVE>              Positive (target) FASTA file(s)
      --negative <NEGATIVE>              Negative (background) FASTA file(s)
  -o, --output <OUTPUT>                  Output bias table file (.bias)
  -k, --kmer-size <KMER_SIZE>            K-mer size (must match sketch) [default: 21]
      --fscale <FSCALE>                  FracMinHash scale (must match sketch) [default: 1000]
      --cms-width <CMS_WIDTH>            CMS columns, power of 2 recommended [default: 1048576]
      --cms-depth <CMS_DEPTH>            CMS hash functions [default: 5]
      --alpha <ALPHA>                    Smoothing parameter for log-ratio [default: 1.0]
      --fold-enrichment <FOLD_ENRICHMENT>  Target fold enrichment (auto-maximized if not set)
      --threads <THREADS>                Number of threads
  -h, --help                             Print help
```

Examples:
```bash
# Build a bias table to filter out host sequences
jam bias create --positive plasmids.fasta --negative host_genome.fasta -o host_filter.bias

# With custom fold enrichment target
jam bias create --positive targets.fasta --negative background.fasta -o filter.bias --fold-enrichment 10.0

# Inspect a bias table
jam bias stats filter.bias
jam bias stats filter.bias -o report.json
```

#### Statistics

Display database statistics including hash counts and distribution analysis.

```console
$ jam stats --help
Display statistics about a JAM database

Usage: jam stats [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>  Input JAM database (.jam file)
      --short          Short summary only
      --full           Include the full entry statistics
  -h, --help           Print help
```

Examples:
```bash
jam stats -i db.jam --short
jam stats -i db.jam --full
```

### License

This project is licensed under the MIT license. See the [LICENSE](LICENSE) file for more info.

### Feedback & Contributions

If you have any ideas, suggestions, or issues, please don't hesitate to open an issue and/or PR. Contributions to this project are always welcome! We appreciate your help in making this project better.

### Credits

This tool is inspired by [finch-rs](https://github.com/onecodex/finch-rs) and [sourmash](https://github.com/sourmash-bio/sourmash). Check them out if you need a more mature ecosystem.