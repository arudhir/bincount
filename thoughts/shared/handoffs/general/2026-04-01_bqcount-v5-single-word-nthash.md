---
date: 2026-04-01
session_name: general
researcher: claude
git_commit: b496e1a
branch: master
repository: bincount
topic: "bqcount v5: single-word CAS cells + ntHash rolling hash"
tags: [optimization, nthash, rolling-hash, cas, hash-table, bioinformatics, kmer-counting]
status: complete
last_updated: 2026-04-01
last_updated_by: claude
type: implementation_strategy
---

# Handoff: bqcount v5 — single-word cells + ntHash

## Task(s)

1. **Single-word CAS cells (v4)** — COMPLETED
   - Packed each slot into one `AtomicU64`: `[occupied:1][remainder:R][count:C]` where R = 64 - table_bits, C = table_bits - 1.
   - Position encodes low `table_bits` of hash, remainder encodes high bits — full 64-bit hash stored, zero false collisions.
   - 8 bytes/slot (down from 12). Default 128M table: 1.07 GB (was 1.61 GB).
   - Parallel `dump_keys` array allocated only with `--dump` for k-mer recovery.
   - `unhash` (inverse splitmix64) included but unused — `dump_keys` approach chosen for ntHash compatibility.

2. **ntHash rolling hash (v5)** — COMPLETED
   - Added `nthash = "0.5"` dependency. Uses `NtHashIterator` (canonical) and `NtHashForwardIterator` (forward-only).
   - O(1) per k-mer instead of O(k) for encode+hash. k-independent performance confirmed: k=21 and k=31 within 3%.
   - Two code paths: `count_kmers_rolling` (ntHash, histogram mode) and `count_kmers_packed` (as_2bit + splitmix64, --dump mode).
   - Refactored insert into `probe` (shared core) + `insert` (splitmix64) + `insert_hash` (pre-computed).
   - Default k changed from 21 to 31.

3. **Correctness validation** — COMPLETED
   - Histogram output identical between ntHash and packed paths (verified on mflorum e0.001 c10 — zero diff).
   - Unique k-mer counts match v3 across all tested datasets.

## Critical References

- `src/main.rs` — entire bqcount implementation (~570 lines, single file)
- Previous handoff: `thoughts/shared/handoffs/general/2026-03-30_19-49-44_bqcount-v3-cas-table-implementation.md`
- nthash crate docs: `docs.rs/nthash/0.5.1`

## Recent changes

- `src/main.rs:10` — Added `nthash` imports
- `src/main.rs:78-87` — `CasKmerTable` struct: single `Vec<AtomicU64>` cells + optional dump_keys + bit layout fields
- `src/main.rs:147-206` — `probe`/`insert`/`insert_hash` — refactored insert with shared core
- `src/main.rs:291-320` — `count_kmers_rolling` using ntHash iterators
- `src/main.rs:322-344` — `count_kmers_packed` (old path, for --dump)
- `src/main.rs:346-367` — `process_record` dispatches based on `dump` flag
- `Cargo.toml:12` — Added `nthash = "0.5"` dependency

## Learnings

1. **ntHash eliminates the O(k) bottleneck**: With splitmix64, the per-kmer cost was `as_2bit(O(k)) + canonical(O(k)) + splitmix64(O(1)) + CAS(O(1) + memory latency)`. ntHash replaces the first three with O(1) rolling XOR+rotate. On ecoli c100, this yielded 3.6x speedup (9.5s → 2.6s at 8 threads).

2. **Thread scaling dropped from 87% to 66%**: v3 single-thread was 62s → 8T was 9s (6.95x). v5 single-thread is 13.7s → 8T is 2.6s (5.3x). This makes sense: ntHash made CPU work so cheap that memory bandwidth saturates earlier. The remaining bottleneck is memory latency (~100ns per random cache miss).

3. **ntHash canonical = min(fwd_hash, rc_hash)**, which is semantically different from our old `min(packed, revcomp_packed)` then hash. Both correctly group reverse complement pairs — just using different representatives. Histogram output is identical because frequency distributions don't depend on which representative is chosen.

4. **nthash crate (v0.5) handles N bases**: Maps N→0 hash seed, panics on other non-ACGTN bases. Since binseq decodes to pure ACGT (2-bit encoding), no N handling needed.

5. **Separating hash computation from table insertion was key**: The `probe(hash, kmer_for_dump)` refactoring cleanly supports both splitmix64 and ntHash paths without code duplication. The `kmer_for_dump` parameter is zero in the ntHash path (dump_keys not allocated), avoiding any overhead.

## Performance Summary

| Dataset | Jellyfish | v3 (splitmix) | v5 (ntHash) | vs JF |
|---------|-----------|---------------|-------------|-------|
| ecoli e0.005 c100 k21 t8 | 23.57s | 9.53s | **2.62s** | 9.0x |
| yeast e0.02 c100 k21 t8 | 104.27s | 31.50s | **9.84s** | 10.6x |
| meta_50_50 e0.02 c100 k21 t8 | 130.61s | 65.01s | **15.64s** | 8.4x |

## Post-Mortem

### What Worked
- **Using the nthash crate** rather than implementing from scratch: 3 lines of code for the rolling iterator, well-tested, handles canonical hashing correctly.
- **Two separate code paths** (rolling vs packed): Clean separation. No performance compromise in either mode.
- **probe/insert/insert_hash refactoring**: Shared CAS logic, no duplication.
- **Correctness validation via histogram diff**: Built histogram from --dump output and diffed against ntHash histogram — zero differences.

### What Failed
- Nothing significant. The implementation was straightforward because the v4 single-word cell design already anticipated ntHash integration (dump_keys array, separable hash function).

### Key Decisions
- **Decision**: Use nthash crate (v0.5) instead of implementing ntHash manually
  - Alternatives: Hand-roll ntHash (~40 lines) for full control over N handling
  - Reason: Crate is well-tested, handles canonical correctly, and binseq has no N bases. No benefit to custom implementation.
- **Decision**: Keep packed path for --dump, ntHash for histogram
  - Alternatives: Use ntHash always and maintain rolling 2-bit encoder in parallel for dump
  - Reason: Simpler. --dump is rare and not performance-critical.
- **Decision**: No splitmix64 finalizer on ntHash output
  - Alternatives: Apply splitmix64 to ntHash for extra avalanche
  - Reason: ntHash has good distribution. For 640M unique k-mers, expected hash collisions ≈ 0.01. Adding splitmix64 can't improve on a bijective transform of already-distributed hashes.

## Action Items & Next Steps

1. **Software prefetching** — MEDIUM PRIORITY
   - Hash the *next* k-mer while waiting for the *current* CAS to resolve
   - Hides ~100ns memory latency per access (now the dominant cost)
   - Could yield significant gains given that CPU work is now minimal

2. **Rerun full benchmarks with v5** — should update results.tsv
   - Current results.tsv has v3 numbers only

3. **Push to GitHub** — still pending

4. **Count overflow guard** — LOW PRIORITY
   - `fetch_add(1)` can carry into remainder bits if count exceeds `count_mask`
   - For 128M table, max count is 67M (safe for practical use)
   - Small tables (e.g. 64K) could overflow with high-coverage data
   - Consider documenting or adding a minimum table size check

## Artifacts

- `src/main.rs` — bqcount v5 (~570 lines)
- `Cargo.toml` — added nthash dependency
- `bench/results/results.tsv` — v3 benchmark data (v5 not yet recorded)
- Previous handoffs in `thoughts/shared/handoffs/general/`
