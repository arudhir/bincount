# bqcount

A parallel k-mer counter for [BINSEQ](https://github.com/noamteyssier/binseq) files, written in Rust.
**8–12× faster than [Jellyfish](https://github.com/gmarcais/Jellyfish)** across bacterial, fungal, and metagenomic datasets at 8 threads, with 100% accuracy on every benchmark.

## Performance

Benchmarked against Jellyfish 2 on a MacBook Pro (Apple M-series) at k=21, 8 threads, full-coverage (`c100`) datasets. Wall-clock time in seconds; lower is better.

### Primary benchmark (8 threads, k=21, coverage=100×)

| Organism | Error rate | Unique k-mers | Jellyfish | bqcount | Speedup |
|---|---|---|---|---|---|
| *M. florum* (0.79 Mb) | 0.001 | 3.5M | 4.01s | 0.51s | **7.9×** |
| *M. florum* | 0.005 | 12.2M | 5.17s | 0.55s | **9.4×** |
| *M. florum* | 0.020 | 32.4M | 7.00s | 0.66s | **10.6×** |
| *E. coli* (4.6 Mb) | 0.001 | 20.6M | 20.16s | 2.21s | **9.1×** |
| *E. coli* | 0.005 | 71.1M | 24.37s | 2.62s | **9.3×** |
| *E. coli* | 0.020 | 188M | 34.71s | 3.62s | **9.6×** |
| *S. cerevisiae* (12.1 Mb) | 0.001 | 52.9M | 64.97s | 5.58s | **11.6×** |
| *S. cerevisiae* | 0.005 | 181M | 72.77s | 6.90s | **10.5×** |
| *S. cerevisiae* | 0.020 | 481M | 112.36s | 9.60s | **11.7×** |
| Metagenome 90/10 | 0.020 | 462M | 123.94s | 12.51s | **9.9×** |
| Metagenome 50/50 | 0.020 | 641M | 135.81s | 16.05s | **8.5×** |

Accuracy: **100% on all 38 configurations** (unique k-mer counts match exactly).

### Thread scaling (*E. coli*, e=0.005, k=21, c=100×)

| Threads | Jellyfish | bqcount | Speedup vs Jellyfish |
|---|---|---|---|
| 1 | 131.50s | 12.19s | 10.8× |
| 2 | 77.51s | 6.54s | 11.9× |
| 4 | 41.53s | 3.96s | 10.5× |
| 8 | 24.37s | 2.62s | 9.3× |

bqcount scales at 4.6× wall-time reduction from 1→8 threads (58% efficiency), bottlenecked by memory bandwidth rather than lock contention.

### K-mer size (*E. coli*, e=0.005, 8 threads, c=100×)

| k | Jellyfish | bqcount | Speedup |
|---|---|---|---|
| 15 | 19.54s | 2.39s | 8.2× |
| 21 | 24.37s | 2.62s | 9.3× |
| 31 | 25.58s | 2.56s | 10.0× |

Runtime is nearly flat across k because ntHash is O(1) per k-mer regardless of k.

---

## How it works

### Lock-free CAS hash table

Each slot in the table is a single `AtomicU64` packed as:

```
[ occupied:1 | remainder:R | count:C ]
```

where `R = 64 - table_bits` and `C = table_bits - 1`. The slot index encodes
the low bits of the hash; the remainder field stores the upper bits. Together
they capture the full 64-bit hash, eliminating false collisions.

- **8 bytes per slot** — vs 12 in a two-array (key + count) layout
- **No locks, ever** — every insert is a compare-exchange loop; threads never block each other
- **Graceful saturation** — if the reprobe limit is hit, a flag is set and a clean error is returned instead of panicking

### ntHash rolling hash

K-mer hashing uses [ntHash](https://github.com/bcgsc/ntHash), a rolling hash
designed for DNA. Once the first k-mer's hash is computed, each subsequent
k-mer updates it in O(1) by XORing out the departing base and XORing in the
arriving base — no matter how large k is. This is in contrast to general-purpose
hash functions that hash all k bases per k-mer.

In `--dump` mode (which needs to recover the k-mer sequence from its hash),
bqcount falls back to a packed 2-bit encoding with splitmix64 finalization,
which is invertible.

### HyperLogLog auto-sizing

Rather than requiring the user to specify a table size or over-provisioning by
default, bqcount does a fast single-pass pre-scan with a
[HyperLogLog](https://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf) sketch
(16 KB state, ~0.8% error) to estimate the number of unique k-mers, then sizes
the table to 1.5× that estimate. This keeps the load factor around 66% while
avoiding the 2× worst-case over-allocation of using the total-k-mer upper bound.

HLL uses 14-bit registers (16,384 buckets) with 5-bit max-leading-zeros values,
plus linear-counting correction for sparse cardinalities.

### BINSEQ format advantage

Input files are in the `.bq` (BINSEQ) format: sequences stored as 2-bit packed
binary with no per-record delimiters or ASCII overhead. This eliminates FASTQ
parsing entirely — records are decoded directly into sliding-window buffers for
ntHash.

---

## Usage

```
bqcount [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Path to .bq or .cbq file

Options:
  -k <K>         K-mer size, max 31 [default: 31]
  -t <T>         Threads (0 = auto-detect) [default: 0]
  -s <S>         Hash table slots (K/M/G suffixes ok). Omit to auto-size via HLL.
  -o <O>         Output file [default: stdout]
  --no-canonical  Count forward and reverse complement separately
  --dump          Output full k-mer table instead of histogram
```

**Histogram mode** (default) — output a `frequency\tcount` TSV:

```bash
bqcount reads.bq -k 21 -t 8
```

**Dump mode** — output every k-mer and its count:

```bash
bqcount reads.bq -k 21 --dump -o kmers.tsv
```

Auto-sizing prints a diagnostic to stderr:

```
HLL pre-pass: estimated ~71075895 unique 21-mers -> auto-sizing to >= 106613842 slots
Hash table: 134217728 slots (1.1 GB, 8.0 bytes/slot, max count 67108863)
```

---

## Build

```bash
cargo build --release
./target/release/bqcount --help
```

Requires Rust 1.75+. Tested on macOS (Apple Silicon) and Linux (x86-64).

---

## Dependencies

| Crate | Role |
|---|---|
| [`binseq`](https://github.com/noamteyssier/binseq) | BINSEQ file format I/O and parallel record dispatch |
| [`nthash`](https://crates.io/crates/nthash) | O(1) rolling hash for DNA k-mers |
| [`bitnuc`](https://crates.io/crates/bitnuc) | 2-bit nucleotide encoding/decoding |
| [`clap`](https://crates.io/crates/clap) | Argument parsing |
| [`hashbrown`](https://crates.io/crates/hashbrown) | HashMap for histogram aggregation |
| [`anyhow`](https://crates.io/crates/anyhow) | Error handling |
