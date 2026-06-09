use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::Result;
use binseq::prelude::*;
use clap::Parser;
use hashbrown::HashMap;
use nthash::{NtHashForwardIterator, NtHashIterator};

/// A binseq-native parallel k-mer counter.
///
/// Counts k-mers from BINSEQ files (.bq, .cbq) and outputs a frequency histogram.
/// Uses a lock-free CAS hash table for concurrent k-mer counting.
#[derive(Parser)]
#[command(name = "bqcount", version)]
struct Cli {
    /// Path to .bq or .cbq file
    input: String,

    /// K-mer size (max: 31)
    #[arg(short, default_value_t = 31)]
    k: usize,

    /// Number of threads (0 = auto-detect)
    #[arg(short, default_value_t = 0)]
    t: usize,

    /// Hash table slots. Must exceed expected unique k-mers. Suffixes: K, M, G.
    /// If omitted, the table is sized automatically from the input file.
    #[arg(short)]
    s: Option<String>,

    /// Output file (default: stdout)
    #[arg(short)]
    o: Option<PathBuf>,

    /// Disable canonical k-mers (count forward and reverse complement separately)
    #[arg(long)]
    no_canonical: bool,

    /// Dump full k-mer count table instead of histogram
    #[arg(long)]
    dump: bool,
}

fn parse_size(s: &str) -> Result<usize> {
    let s = s.trim();
    let last = s.as_bytes().last().copied().unwrap_or(b'0');
    let (num_str, mult) = match last {
        b'G' | b'g' => (&s[..s.len() - 1], 1_000_000_000usize),
        b'M' | b'm' => (&s[..s.len() - 1], 1_000_000),
        b'K' | b'k' => (&s[..s.len() - 1], 1_000),
        _ => (s, 1),
    };
    let n: usize = num_str
        .parse()
        .map_err(|e| anyhow::anyhow!("Invalid size '{}': {}", s, e))?;
    Ok(n * mult)
}

/// Default table size used when the input provides no usable size hint
/// (e.g. variable-length files with slen == 0).
const DEFAULT_SLOTS: usize = 128_000_000;

/// Hard cap on auto-sized tables: 8G slots = 64 GB at 8 bytes/slot.
const MAX_AUTO_SLOTS: usize = 8usize << 30; // 8_589_934_592

/// Estimate a sensible number of hash-table slots from input metadata.
///
/// `num_records` is the number of records in the file, `slen`/`xlen` are the
/// per-record primary/extended sequence lengths from the file header (xlen == 0
/// when unpaired). `k` is the k-mer size.
///
/// An upper bound on total k-mers is `num_records * ((slen - k + 1) + (xlen - k + 1))`.
/// Unique k-mers can never exceed this. We target `2 * upper_bound` slots
/// (then rounded up to a power of two by `CasKmerTable::new`) so the table stays
/// well below the ~76% load factor that triggers reprobe saturation, with a floor
/// of `DEFAULT_SLOTS` and a cap of `MAX_AUTO_SLOTS`.
///
/// Returns `None` when no useful estimate can be made (slen too short / zero),
/// signalling the caller to fall back to the default size.
fn estimate_table_size(num_records: usize, slen: usize, xlen: usize, k: usize) -> Option<usize> {
    let per_seq = |len: usize| -> usize {
        if len >= k {
            len - k + 1
        } else {
            0
        }
    };
    let kmers_per_record = per_seq(slen).saturating_add(per_seq(xlen));
    if kmers_per_record == 0 || num_records == 0 {
        return None;
    }
    let upper_bound = (num_records as u128) * (kmers_per_record as u128);
    // Target 2x headroom over the upper bound on unique k-mers.
    let target = upper_bound.saturating_mul(2);
    let target = target.min(MAX_AUTO_SLOTS as u128) as usize;
    Some(target.max(DEFAULT_SLOTS))
}

// ---- Lock-free CAS k-mer table (single-word packed cells) ----

const OCCUPIED_BIT: u64 = 1 << 63;
const MAX_REPROBE: usize = 256;

/// Lock-free concurrent hash table for k-mer counting.
///
/// Each slot is a single `AtomicU64` packed as:
///   `[occupied:1][remainder:R][count:C]`
/// where `R = 64 - table_bits` and `C = table_bits - 1`.
///
/// The slot index encodes the low `table_bits` of the hash; the remainder
/// field stores the upper bits. Together they capture the full 64-bit hash,
/// eliminating false collisions (splitmix64 is bijective).
///
/// 8 bytes per slot (vs 12 in the two-array layout).
struct CasKmerTable {
    cells: Vec<AtomicU64>,
    /// Parallel array storing original k-mers — only allocated with --dump.
    dump_keys: Option<Vec<AtomicU64>>,
    mask: usize,
    slots: usize,
    table_bits: u32,
    count_bits: u32,
    count_mask: u64,
    /// Set by worker threads when the table saturates (reprobe limit hit)
    /// instead of panicking. Checked after parallel processing completes.
    saturated: AtomicBool,
}

impl CasKmerTable {
    fn new(min_capacity: usize, dump: bool) -> Result<Self> {
        let slots = min_capacity.next_power_of_two();
        let table_bits = slots.trailing_zeros();

        if table_bits < 5 {
            anyhow::bail!(
                "Table too small: {} slots (2^{}). Need at least 32 slots.",
                slots, table_bits
            );
        }

        let count_bits = table_bits - 1;
        let count_mask = (1u64 << count_bits) - 1;
        let cells = (0..slots).map(|_| AtomicU64::new(0)).collect();
        let dump_keys = if dump {
            Some((0..slots).map(|_| AtomicU64::new(0)).collect())
        } else {
            None
        };

        Ok(Self {
            cells,
            dump_keys,
            mask: slots - 1,
            slots,
            table_bits,
            count_bits,
            count_mask,
            saturated: AtomicBool::new(false),
        })
    }

    /// Record that the table saturated (reprobe limit hit) in a worker thread.
    #[inline]
    fn mark_saturated(&self) {
        self.saturated.store(true, Ordering::Relaxed);
    }

    /// Whether any worker thread hit the reprobe limit during insertion.
    fn is_saturated(&self) -> bool {
        self.saturated.load(Ordering::Relaxed)
    }

    /// splitmix64 finalizer — fast, good avalanche, bijective on u64.
    #[inline(always)]
    fn hash(kmer: u64) -> u64 {
        let mut h = kmer;
        h ^= h >> 30;
        h = h.wrapping_mul(0xbf58476d1ce4e5b9);
        h ^= h >> 27;
        h = h.wrapping_mul(0x94d049bb133111eb);
        h ^= h >> 31;
        h
    }

    /// Inverse splitmix64 — recovers k-mer from hash.
    #[allow(dead_code)]
    fn unhash(mut h: u64) -> u64 {
        h ^= h >> 31;
        h ^= h >> 62;
        h = h.wrapping_mul(0x319642b2d24d8ec3);
        h ^= h >> 27;
        h ^= h >> 54;
        h = h.wrapping_mul(0x96de1b173f119089);
        h ^= h >> 30;
        h ^= h >> 60;
        h
    }

    /// Probe-and-insert core. `kmer_for_dump` is stored in dump_keys when
    /// claiming a new slot (only needed in --dump mode).
    #[inline(always)]
    fn probe(&self, hash: u64, kmer_for_dump: u64) -> bool {
        let pos = hash as usize & self.mask;
        let remainder = hash >> self.table_bits;
        let new_cell = OCCUPIED_BIT | (remainder << self.count_bits) | 1;

        for reprobe in 0..MAX_REPROBE {
            let idx = (pos + reprobe) & self.mask;
            let cell = &self.cells[idx];
            let mut current = cell.load(Ordering::Relaxed);

            loop {
                if current == 0 {
                    // Empty slot — try to claim it
                    match cell.compare_exchange_weak(
                        0,
                        new_cell,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => {
                            if let Some(ref keys) = self.dump_keys {
                                keys[idx].store(kmer_for_dump, Ordering::Relaxed);
                            }
                            return true;
                        }
                        Err(actual) => {
                            current = actual;
                            continue; // Re-check: might be our key or different
                        }
                    }
                }

                // Compare remainder (strip occupied bit, shift out count)
                let cell_remainder = (current & !OCCUPIED_BIT) >> self.count_bits;
                if cell_remainder == remainder {
                    // Same hash — increment count
                    cell.fetch_add(1, Ordering::Relaxed);
                    return true;
                }

                break; // Different key — move to next probe position
            }
        }
        false
    }

    /// Insert a packed k-mer (splitmix64 hash, stores k-mer for --dump).
    #[inline]
    fn insert(&self, kmer: u64) -> bool {
        self.probe(Self::hash(kmer), kmer)
    }

    /// Insert a pre-computed hash (e.g. from ntHash). No k-mer recovery.
    #[inline]
    fn insert_hash(&self, hash: u64) -> bool {
        self.probe(hash, 0)
    }

    fn occupied(&self) -> usize {
        self.cells
            .iter()
            .filter(|c| c.load(Ordering::Relaxed) != 0)
            .count()
    }

    /// Iterate over (k-mer, count) pairs. Requires dump mode for k-mer recovery.
    fn iter_kmers(&self) -> impl Iterator<Item = (u64, u32)> + '_ {
        let keys = self
            .dump_keys
            .as_ref()
            .expect("iter_kmers requires --dump mode");
        self.cells
            .iter()
            .enumerate()
            .filter_map(move |(idx, cell)| {
                let val = cell.load(Ordering::Relaxed);
                if val == 0 {
                    return None;
                }
                let count = (val & self.count_mask) as u32;
                let kmer = keys[idx].load(Ordering::Relaxed);
                Some((kmer, count))
            })
    }

    /// Iterate over counts only (for histogram mode — no k-mer recovery needed).
    fn iter_counts(&self) -> impl Iterator<Item = u32> + '_ {
        self.cells.iter().filter_map(move |cell| {
            let val = cell.load(Ordering::Relaxed);
            if val == 0 {
                return None;
            }
            Some((val & self.count_mask) as u32)
        })
    }

    fn memory_bytes(&self) -> usize {
        let base = self.slots * 8;
        if self.dump_keys.is_some() {
            base + self.slots * 8
        } else {
            base
        }
    }

    fn max_count(&self) -> u64 {
        self.count_mask
    }
}

// ---- K-mer functions ----

/// Reverse complement of a 2-bit packed k-mer.
#[inline]
fn revcomp_2bit(packed: u64, k: usize) -> u64 {
    let comp = !packed;
    let mut rev = comp;
    rev = ((rev & 0x3333_3333_3333_3333) << 2) | ((rev >> 2) & 0x3333_3333_3333_3333);
    rev = ((rev & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((rev >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
    rev = rev.swap_bytes();
    rev >> (64 - 2 * k)
}

/// Canonical k-mer: min(forward, reverse complement)
#[inline]
fn canonical_kmer(packed: u64, k: usize) -> u64 {
    let rc = revcomp_2bit(packed, k);
    packed.min(rc)
}

// ---- Parallel processor ----

#[derive(Clone)]
struct KmerCounter {
    k: usize,
    canonical: bool,
    dump: bool,
    table: Arc<CasKmerTable>,
    seq_buf: Vec<u8>,
}

/// Count k-mers using ntHash rolling hash (O(1) per k-mer). Histogram mode only.
///
/// On reprobe-limit saturation, sets the table's `saturated` flag and returns
/// early instead of panicking; the condition is surfaced cleanly by `main`.
#[inline]
fn count_kmers_rolling(seq: &[u8], k: usize, canonical: bool, table: &CasKmerTable) {
    if seq.len() < k {
        return;
    }
    if canonical {
        if let Ok(iter) = NtHashIterator::new(seq, k) {
            for hash in iter {
                if !table.insert_hash(hash) {
                    table.mark_saturated();
                    return;
                }
            }
        }
    } else if let Ok(iter) = NtHashForwardIterator::new(seq, k) {
        for hash in iter {
            if !table.insert_hash(hash) {
                table.mark_saturated();
                return;
            }
        }
    }
}

/// Count k-mers with packed 2-bit encoding (O(k) per k-mer). Required for --dump mode.
///
/// On reprobe-limit saturation, sets the table's `saturated` flag and returns
/// early instead of panicking; the condition is surfaced cleanly by `main`.
#[inline]
fn count_kmers_packed(seq: &[u8], k: usize, canonical: bool, table: &CasKmerTable) {
    if seq.len() < k {
        return;
    }
    for window in seq.windows(k) {
        if let Ok(packed) = bitnuc::as_2bit(window) {
            let key = if canonical {
                canonical_kmer(packed, k)
            } else {
                packed
            };
            if !table.insert(key) {
                table.mark_saturated();
                return;
            }
        }
    }
}

impl ParallelProcessor for KmerCounter {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        self.seq_buf.clear();
        record.decode_s(&mut self.seq_buf)?;
        if self.dump {
            count_kmers_packed(&self.seq_buf, self.k, self.canonical, &self.table);
        } else {
            count_kmers_rolling(&self.seq_buf, self.k, self.canonical, &self.table);
        }

        if record.is_paired() {
            self.seq_buf.clear();
            record.decode_x(&mut self.seq_buf)?;
            if self.dump {
                count_kmers_packed(&self.seq_buf, self.k, self.canonical, &self.table);
            } else {
                count_kmers_rolling(&self.seq_buf, self.k, self.canonical, &self.table);
            }
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        Ok(()) // No-op: k-mers inserted directly via CAS
    }
}

/// Read the per-record primary (slen) and extended (xlen) sequence lengths
/// from a reader's file header.
///
/// Only the fixed-length BQ format stores per-record sequence lengths in its
/// header (`slen`/`xlen`). The CBQ and VBQ formats are block-based /
/// variable-length and expose no file-level fixed sequence length, so we
/// return `(0, 0)` for them — `estimate_table_size` then yields `None` and the
/// caller falls back to the default table size.
fn reader_seq_lens(reader: &BinseqReader) -> (usize, usize) {
    match reader {
        BinseqReader::Bq(r) => {
            let h = r.header();
            (h.slen as usize, h.xlen as usize)
        }
        BinseqReader::Cbq(_) | BinseqReader::Vbq(_) => (0, 0),
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.k == 0 || cli.k > 31 {
        anyhow::bail!("k must be between 1 and 31, got {}", cli.k);
    }

    let reader = BinseqReader::new(&cli.input)?;
    let canonical = !cli.no_canonical;

    // Determine table size: honor an explicit -s exactly; otherwise auto-size
    // from the input file metadata (record count + per-record sequence length).
    let table_size = match &cli.s {
        Some(s) => parse_size(s)?,
        None => {
            let num_records = reader.num_records()?;
            let (slen, xlen) = reader_seq_lens(&reader);
            // xlen only contributes when the file is actually paired.
            let xlen = if reader.is_paired() { xlen } else { 0 };
            match estimate_table_size(num_records, slen, xlen, cli.k) {
                Some(slots) => {
                    if slots >= MAX_AUTO_SLOTS {
                        eprintln!(
                            "warning: auto-sized table capped at {} slots (~{:.0} GB). \
                             If the genome has more unique {}-mers than this, the table \
                             may saturate; pass an explicit -s to override.",
                            MAX_AUTO_SLOTS,
                            (MAX_AUTO_SLOTS * 8) as f64 / 1e9,
                            cli.k,
                        );
                    }
                    eprintln!(
                        "Auto-sizing hash table from input: {} records, slen={}, xlen={} \
                         -> >= {} slots requested",
                        num_records, slen, xlen, slots,
                    );
                    slots
                }
                None => {
                    eprintln!(
                        "warning: could not estimate table size from input \
                         (slen={}, records={}); using default {} slots. \
                         Pass -s to override.",
                        slen, num_records, DEFAULT_SLOTS,
                    );
                    DEFAULT_SLOTS
                }
            }
        }
    };

    let table = Arc::new(CasKmerTable::new(table_size, cli.dump)?);

    eprintln!(
        "Hash table: {} slots ({:.1} GB, {:.0} bytes/slot, max count {})",
        table.slots,
        table.memory_bytes() as f64 / 1e9,
        table.memory_bytes() as f64 / table.slots as f64,
        table.max_count(),
    );

    let counter = KmerCounter {
        k: cli.k,
        canonical,
        dump: cli.dump,
        table: table.clone(),
        seq_buf: Vec::new(),
    };

    eprintln!(
        "Counting {}-mers{} from {} using {} thread(s)...",
        cli.k,
        if canonical { " (canonical)" } else { "" },
        cli.input,
        if cli.t == 0 {
            "auto".to_string()
        } else {
            cli.t.to_string()
        },
    );

    reader.process_parallel(counter, cli.t)?;

    // Graceful failure: if any worker thread hit the reprobe limit, the table
    // saturated. Report cleanly and exit nonzero instead of panicking.
    if table.is_saturated() {
        anyhow::bail!(
            "hash table saturated at {} slots; rerun with a larger -s value",
            table.slots
        );
    }

    let unique = table.occupied();
    eprintln!("Found {} unique k-mers", unique);

    // Output
    let out: Box<dyn Write> = match &cli.o {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    let mut out = out;

    if cli.dump {
        writeln!(out, "kmer\tcount")?;
        let mut buf = Vec::new();
        for (kmer, count) in table.iter_kmers() {
            buf.clear();
            bitnuc::from_2bit(kmer, cli.k, &mut buf)?;
            let seq = std::str::from_utf8(&buf)?;
            writeln!(out, "{}\t{}", seq, count)?;
        }
    } else {
        // Build histogram: frequency -> number of k-mers with that frequency
        let mut histogram: HashMap<u32, u64> = HashMap::new();
        for count in table.iter_counts() {
            *histogram.entry(count).or_insert(0) += 1;
        }

        let mut freqs: Vec<_> = histogram.into_iter().collect();
        freqs.sort_unstable_by_key(|&(freq, _)| freq);

        writeln!(out, "frequency\tcount")?;
        for (freq, num_kmers) in freqs {
            writeln!(out, "{}\t{}", freq, num_kmers)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp_2bit() {
        let acg = bitnuc::as_2bit(b"ACG").unwrap();
        let rc = revcomp_2bit(acg, 3);
        let cgt = bitnuc::as_2bit(b"CGT").unwrap();
        assert_eq!(rc, cgt, "revcomp(ACG) should be CGT");
    }

    #[test]
    fn test_revcomp_palindrome() {
        let acgt = bitnuc::as_2bit(b"ACGT").unwrap();
        let rc = revcomp_2bit(acgt, 4);
        assert_eq!(rc, acgt, "ACGT should be its own reverse complement");
    }

    #[test]
    fn test_canonical_kmer() {
        let acg = bitnuc::as_2bit(b"ACG").unwrap();
        let cgt = bitnuc::as_2bit(b"CGT").unwrap();
        assert_eq!(
            canonical_kmer(acg, 3),
            canonical_kmer(cgt, 3),
            "ACG and CGT (revcomp pair) should have same canonical form"
        );
    }

    #[test]
    fn test_revcomp_single_base() {
        let a = bitnuc::as_2bit(b"A").unwrap();
        let t = bitnuc::as_2bit(b"T").unwrap();
        assert_eq!(revcomp_2bit(a, 1), t, "revcomp(A) should be T");

        let c = bitnuc::as_2bit(b"C").unwrap();
        let g = bitnuc::as_2bit(b"G").unwrap();
        assert_eq!(revcomp_2bit(c, 1), g, "revcomp(C) should be G");
    }

    #[test]
    fn test_hash_unhash_roundtrip() {
        for &kmer in &[0u64, 1, 42, 99, 0xDEADBEEF, u64::MAX] {
            assert_eq!(
                CasKmerTable::unhash(CasKmerTable::hash(kmer)),
                kmer,
                "unhash(hash({})) should equal {}",
                kmer,
                kmer
            );
        }
    }

    #[test]
    fn test_cas_table_basic() {
        let table = CasKmerTable::new(1024, true).unwrap();
        assert!(table.insert(42));
        assert!(table.insert(42));
        assert!(table.insert(99));

        let mut found = HashMap::new();
        for (kmer, count) in table.iter_kmers() {
            found.insert(kmer, count);
        }
        assert_eq!(found[&42], 2);
        assert_eq!(found[&99], 1);
        assert_eq!(table.occupied(), 2);
    }

    #[test]
    fn test_cas_table_zero_kmer() {
        // k-mer 0 (all A's) must not collide with empty sentinel
        let table = CasKmerTable::new(1024, true).unwrap();
        assert!(table.insert(0));
        assert!(table.insert(0));

        let mut found = HashMap::new();
        for (kmer, count) in table.iter_kmers() {
            found.insert(kmer, count);
        }
        assert_eq!(found[&0], 2);
        assert_eq!(table.occupied(), 1);
    }

    #[test]
    fn test_cas_table_many_kmers() {
        let table = CasKmerTable::new(4096, false).unwrap();
        for i in 0..1000u64 {
            assert!(table.insert(i));
        }
        assert_eq!(table.occupied(), 1000);

        let total: u32 = table.iter_counts().sum();
        assert_eq!(total, 1000);
    }

    #[test]
    fn test_single_word_memory() {
        let table = CasKmerTable::new(1024, false).unwrap();
        assert_eq!(table.memory_bytes(), 1024 * 8);

        let table_dump = CasKmerTable::new(1024, true).unwrap();
        assert_eq!(table_dump.memory_bytes(), 1024 * 16);
    }

    #[test]
    fn test_count_bits_layout() {
        // 128M table (2^27): 26 count bits, max count 67M
        let table = CasKmerTable::new(128_000_000, false).unwrap();
        assert_eq!(table.slots, 134_217_728); // 2^27
        assert_eq!(table.table_bits, 27);
        assert_eq!(table.count_bits, 26);
        assert_eq!(table.max_count(), (1 << 26) - 1);
    }

    #[test]
    fn test_estimate_table_size_small_floors_to_default() {
        // A tiny input can't justify a table smaller than the default floor.
        let est = estimate_table_size(10, 100, 0, 31).unwrap();
        assert_eq!(est, DEFAULT_SLOTS);
    }

    #[test]
    fn test_estimate_table_size_large_genome() {
        // ~6.7 Mb genome shredded into reads that overflowed the 128M default.
        // upper bound = num_records * (slen - k + 1); target = 2 * upper_bound.
        let num_records = 1_000_000usize;
        let slen = 150usize;
        let k = 31usize;
        let kmers_per_record = slen - k + 1; // 120
        let upper = num_records * kmers_per_record; // 120,000,000
        let est = estimate_table_size(num_records, slen, 0, k).unwrap();
        assert_eq!(est, (upper * 2).max(DEFAULT_SLOTS)); // 240,000,000
        // And it must comfortably exceed the unique-kmer upper bound.
        assert!(est >= upper);
    }

    #[test]
    fn test_estimate_table_size_paired_counts_both_ends() {
        let est_single = estimate_table_size(1_000_000, 150, 0, 31).unwrap();
        let est_paired = estimate_table_size(1_000_000, 150, 150, 31).unwrap();
        assert!(est_paired > est_single, "paired input should request more slots");
    }

    #[test]
    fn test_estimate_table_size_caps() {
        // Absurdly large input must be capped at MAX_AUTO_SLOTS, not overflow.
        let est = estimate_table_size(usize::MAX, 1000, 1000, 31).unwrap();
        assert_eq!(est, MAX_AUTO_SLOTS);
    }

    #[test]
    fn test_estimate_table_size_none_when_unusable() {
        // slen shorter than k, or zero records => no estimate possible.
        assert!(estimate_table_size(100, 10, 0, 31).is_none());
        assert!(estimate_table_size(0, 150, 0, 31).is_none());
        assert!(estimate_table_size(100, 0, 0, 31).is_none());
    }
}
