[![Rust](https://img.shields.io/badge/built_with-Rust-dca282.svg)](https://www.rust-lang.org/)
[![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)](https://github.com/St4NNi/jam-rs/blob/main/LICENSE)
![CI](https://github.com/St4NNi/jam-rs/actions/workflows/push.yaml/badge.svg)
[![Codecov](https://codecov.io/github/St4NNi/jam-rs/coverage.svg?branch=main)](https://codecov.io/gh/St4NNi/jam-rs)
[![Dependency status](https://deps.rs/repo/github/St4NNi/jam-rs/status.svg)](https://deps.rs/repo/github/St4NNi/jam-rs)
___
# jam-rs

Just another minhash (jam) implementation. A high performance and lossy minhash variant that can be used to screen extremely large datasets in a very short timeframe.

Implements parts of the ScaledMinHash / FracMinHash algorithm described in [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027).

Unlike traditional implementations like [sourmash](https://joss.theoj.org/papers/10.21105/joss.00027) or [mash](https://doi.org/10.1186/s13059-016-0997-x) that focus
on accurately predicting Jaccard similarities, this implementation focuses on raw speed to get a quick (yet lossy) overview. This can be used to screen terabytes of data in just a few seconds / minutes.

### Comparison

- xxhash3 instead of murmurhash
- No jaccard similarity
- Quick estimation of hash fractions from file-size


