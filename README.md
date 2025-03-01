[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
[![Crates.io](https://img.shields.io/crates/v/jam-rs.svg)](https://crates.io/crates/jam-rs)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)
# jam-rs

Just another minhash (jam) implementation. A high performance minhash variant to screen extremely large (metagenomic) datasets in a very short timeframe.
Implements parts of the ScaledMinHash / FracMinHash algorithm described in [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027).

Unlike traditional implementations like [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027) or [mash](https://doi.org/10.1186/s13059-016-0997-x) this version tries to focus on estimating the containment of small sequences in large sets by (optionally) introducing an intentional bias towards smaller sequences and enforcing a certain frequence. This is intended to be used to screen terabytes of data in just a few seconds / minutes.

### Installation

A pre-release is published via [crates.io](https://crates.io/) to install it use (you need to have `cargo` and the `rust-toolchain` installed, the easiest way is via [rustup.rs](https://rustup.rs/)):

```bash
cargo install jam-rs@0.1.0-beta.2
```

If you want the bleeding edge development release you can install it via git:

```bash
cargo install --git https://github.com/St4NNi/jam-rs
```

### Comparison

- Multiple algorithms: [xxhash3](https://github.com/DoumanAsh/xxhash-rust), [ahash-fallback](https://github.com/tkaitchuck/aHash/wiki/AHash-fallback-algorithm) (for kmer < 32) and legacy [murmurhash3](https://github.com/mhallin/murmurhash3-rs)
- Additional filter and sketching options to increase for specificity and sensitivity for small sequences in collections of large assembled metagenomes
- Sketch to a memory mapped database including additional metadata

### Scaling methods

Multiple different scaling methods:
  - FracMinHash (`fscale`): Restricts the hash-space to a (lower) maximum fraction of `u64::MAX` / `fscale`
  - Bias scaling, introduce a bias towards smaller sequences, this increases the Database size but ensures that smaller sequences are better covered

If `KmerCountScaling` and `MinMaxAbsoluteScaling` are used together the minimum number of hashes (per sequence record) will be guaranteed. `FracMinHash` and `KmerCountScaling` produce similar results, the first is mainly provided for sourmash compatibility.

### Usage

```console
$ jam
Just another (genomic) minhasher (jam), obviously blazingly fast

Usage: jam [OPTIONS] <COMMAND>

Commands:
  sketch  Sketch one or more files and write result to output file (or stdout)
  merge   Merge multiple input sketches into a single sketch
  dist    Estimate distance of a (small) sketch against a subset of one or more sketches as database. Requires all sketches to have the same kmer size
  help    Print this message or the help of the given subcommand(s)

Options:
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -h, --help               Print help (see more with '--help')
  -V, --version            Print version
```

#### Sketching

The easiest way to sketch files is to use the `jam sketch` command. This accepts one or more input files (fastx / fastx.gz) or a `.list` file with a full list of input files. And sketches all inputs to a specific outpuf sketch file.

```console
$ jam sketch
Sketch one or more files and write the result to an output file (or stdout)

Usage: jam sketch [OPTIONS] [INPUT]...

Arguments:
  [INPUT]...  Input file(s), one directory or one file with list of files to be hashed

Options:
  -o, --output <OUTPUT>        Output file
  -k, --kmer-size <KMER_SIZE>  kmer size, all sketches must have the same size to be compared [default: 21]
      --fscale <FSCALE>        Scale the hash space to a minimum fraction of the maximum hash value (FracMinHash)
  -t, --threads <THREADS>      Number of threads to use [default: 1]
  -f, --force                  Overwrite output files
      --nmin <NMIN>            Minimum number of k-mers (per record) to be hashed, bottom cut-off
      --nmax <NMAX>            Maximum number of k-mers (per record) to be hashed, top cut-off
      --format <FORMAT>        Change to other output formats [default: bin] [possible values: bin, sourmash]
      --algorithm <ALGORITHM>  Change the hashing algorithm [default: default] [possible values: default, ahash, xxhash, murmur3]
      --singleton              Create a separate sketch for each sequence record
  -s, --stats                  Keep extra stats for each sequence record
  -h, --help                   Print help
```

#### Dist

Calculate the distance for one or more inputs vs. a large set of database sketches. Optionally specify a minimum cutoff in percent of matching kmers. Output is optional if not specified the result will be printed to stdout.

```console
$ jam dist
Estimate containment of a (small) sketch against a subset of one or more sketches as database. Requires all sketches to have the same kmer size

Usage: jam dist [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>        Input sketch or raw file
  -d, --database <DATABASE>  Database sketch(es)
  -o, --output <OUTPUT>      Output to file instead of stdout
  -c, --cutoff <CUTOFF>      Cut-off value for similarity [default: 0.0]
  -t, --threads <THREADS>    Number of threads to use [default: 1]
  -f, --force                Overwrite output files
      --stats                Use the Stats params for restricting results
      --gc-lower <GC_LOWER>  Use GC stats with an upper bound of x% (gc_lower and gc_upper must be set)
      --gc-upper <GC_UPPER>  Use GC stats with an lower bound of y% (gc_lower and gc_upper must be set)
  -h, --help                 Print help
```


#### Merge

Merge multiple sketches into one large one.

```console
$ jam merge
Merge multiple input sketches into a single sketch

Usage: jam merge [OPTIONS] --output <OUTPUT> [INPUTS]...

Arguments:
  [INPUTS]...  One or more input sketches

Options:
  -o, --output <OUTPUT>    Output file
  -t, --threads <THREADS>  Number of threads to use [default: 1]
  -f, --force              Overwrite output files
  -h, --help               Print help
```

### License

This project is licensed under the MIT license. See the [LICENSE](LICENSE) file for more info.

### Disclaimer

jam-rs is still in active development and not ready for production use. Use at your own risk.

### Feedback & Contributions

If you have any ideas, suggestions, or issues, please don't hesitate to open an issue and/or PR. Contributions to this project are always welcome ! We appreciate your help in making this project better 

### Credits

This tool is heavily inspired by [finch-rs](https://github.com/onecodex/finch-rs)/[License](https://github.com/onecodex/finch-rs/blob/master/LICENSE.txt) and [sourmash](https://github.com/sourmash-bio/sourmash)/[License](https://github.com/sourmash-bio/sourmash/blob/latest/LICENSE). Check them out if you need a more mature ecosystem with well tested hash functions and more features.
