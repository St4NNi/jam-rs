[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)

# jam-rs

Just another minhash (jam) implementation. A high-performance minhashing tool for genomic sequence similarity analysis, specifically optimized for plasmids and other small genomic elements.

Implements the FracMinHash algorithm for rapid similarity comparisons with enhanced metadata tracking including GC content and sequence length categories tuned for typical plasmid ranges.

### Installation

Install the latest release from [crates.io](https://crates.io/):

```bash
cargo install jam-rs
```

Or install the development version from git:

```bash
cargo install --git https://github.com/St4NNi/jam-rs
```

### Key Features

- **Plasmid-optimized**: GC content and length categories specifically tuned for plasmid analysis (30-70% GC, 1kB-500kB lengths)
- **Fast sketching**: Entropy-filtered k-mers with optimized hash functions to exclude low-complexity regions
- **Rich metadata**: Enhanced metadata tracking with file index, GC category, and length category for each hash
- **Memory-efficient**: External sorting for processing datasets larger than available RAM
- **Compact storage**: Fast random access with bucket-based `.jam` file format
- **Parallel execution**: File-level parallelization with configurable thread count

### Scaling Methods

Multiple scaling methods for different use cases:
- **FracMinHash** (`--fscale`): Restricts hash-space to a fraction of `u64::MAX` / `fscale`
- **Complexity filtering** (`--complexity`): Only hash sequences with Shannon entropy above threshold (default: 0.0)
- **Singleton mode** (`--singleton`): Creates separate sketch per sequence record
- **Bias filtering** (`--bias-table`): Apply compositional bias filtering using a pre-built bias table

### Usage

```console
$ jam --help
Just another (genomic) minhasher (jam), obviously blazingly fast

Usage: jam [OPTIONS] <COMMAND>

Commands:
  sketch  Sketch one or more files and write the result to an output file
  dist    Estimate containment of a query sequence against a sketch database
  bias    Build a bias table from positive/negative reference sequences
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

Create sketches from FASTA/FASTQ files. Supports single files, multiple files, or directories.

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
      --complexity <COMPLEXITY>  Complexity cut-off, only hash sequences with complexity above this value [default: 0.0]
      --singleton                Create a separate sketch for each sequence record
      --temp-dir <TEMP_DIR>      Custom temporary directory for intermediate files during sorting
      --bias-table <BIAS_TABLE>  Path to a bias table file (.bias) for compositional filtering
  -h, --help                     Print help
```

Examples:
```bash
# Basic plasmid sketching
jam sketch plasmid.fasta -o sketch.jam

# Multiple plasmid files with custom k-mer size
jam sketch plasmids/ -o plasmid_db.jam -k 21 -t 8

# Large collections with FracMinHash scaling and complexity filtering
jam sketch large_collection/ -o database.jam --fscale 1000000 --complexity 1.5

# Separate sketch per plasmid sequence
jam sketch multi_plasmids.fasta -o sketches.jam --singleton

# With bias filtering
jam sketch plasmids/ -o filtered.jam --bias-table host_filter.bias
```

#### Distance Calculation

Compare query sequences against a sketch database.

```console
$ jam dist --help
Estimate containment of a query sequence against a sketch database

Usage: jam dist [OPTIONS] --input <INPUT> --database <DATABASE>

Options:
  -i, --input <INPUT>        Input FASTA/FASTQ file to query
  -d, --database <DATABASE>  Database sketch (.jam file)
  -o, --output <OUTPUT>      Output to file instead of stdout
  -c, --cutoff <CUTOFF>      Cut-off value for similarity/containment [default: 0.0]
      --singleton            Singleton mode, process each query sequence separately
      --bias-table <PATH>    Path to a bias table file (.bias) for query filtering
  -h, --help                 Print help
```

Examples:
```bash
# Query plasmid against database
jam dist -i query_plasmid.fasta -d plasmid_db.jam -c 0.1 -o results.tsv

# Process each sequence separately with singleton mode
jam dist -i multi_query.fasta -d plasmid_db.jam --singleton -c 0.1
```

#### Bias Table Construction

Build a bias table from positive (target) and negative (background) reference sequences for compositional filtering.

```console
$ jam bias --help
Build a bias table from positive/negative reference sequences

Usage: jam bias [OPTIONS] --positive <POSITIVE> --negative <NEGATIVE> --output <OUTPUT>

Options:
  -p, --positive <POSITIVE>  Positive reference sequences (FASTA) - k-mers enriched here will be preferred
  -n, --negative <NEGATIVE>  Negative reference sequences (FASTA) - background/unwanted sequences
  -o, --output <OUTPUT>      Output bias table file (.bias)
      --threshold <VALUE>    Score threshold (0.0-1.0). Higher = stricter filtering [default: 0.5]
  -h, --help                 Print help
```

Examples:
```bash
# Build a bias table to filter out host sequences
jam bias -p plasmids.fasta -n host_genome.fasta -o host_filter.bias --threshold 0.7
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
# Summary statistics
jam stats -i plasmid_db.jam --short

# Detailed distributions
jam stats -i plasmid_db.jam --full
```

### Output Format

Distance results are tab-separated with columns:
```
query	sample_id	hit_count	containment
```

- `query`: Query sequence name (or "combined" in non-singleton mode)
- `sample_id`: Matched sample identifier from the database
- `hit_count`: Number of shared hashes between query and database sample
- `containment`: Containment estimate (fraction of query hashes found in sample)

### Algorithm

JAM uses entropy-filtered k-mers to exclude low-complexity regions, stores rich metadata (file index, GC category, length category) with each hash, and employs external sorting for memory-efficient processing of large datasets. The categorization system is specifically tuned for plasmid analysis, with fine-grained bins in typical plasmid size and GC content ranges.

### License

This project is licensed under the MIT license. See the [LICENSE](LICENSE) file for more info.

### Feedback & Contributions

If you have any ideas, suggestions, or issues, please don't hesitate to open an issue and/or PR. Contributions to this project are always welcome! We appreciate your help in making this project better.

### Credits

This tool is inspired by [finch-rs](https://github.com/onecodex/finch-rs) and [sourmash](https://github.com/sourmash-bio/sourmash). Check them out if you need a more mature ecosystem with well tested hash functions and more features.