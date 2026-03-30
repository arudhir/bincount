use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::sync::Arc;

use anyhow::Result;
use binseq::prelude::*;
use clap::Parser;
use hashbrown::HashMap;

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
    #[arg(short, default_value_t = 21)]
    k: usize,

    /// Number of threads (0 = auto-detect)
    #[arg(short, default_value_t = 0)]
    t: usize,

    /// Hash table slots. Must exceed expected unique k-mers. Suffixes: K, M, G.
    #[arg(short, default_value = "128M")]
    s: String,

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

// ---- Lock-free CAS k-mer table ----

const OCCUPIED_BIT: u64 = 1 << 63;
const MAX_REPROBE: usize = 256;

/// Lock-free concurrent hash table for k-mer counting.
///
/// Each slot has a key (`AtomicU64`, with bit 63 as occupied flag) and a
/// count (`AtomicU32`). Insertions use CAS on the key and `fetch_add` on
/// the count — no locks, no merges, no per-thread maps.
struct CasKmerTable {
    keys: Vec<AtomicU64>,
    counts: Vec<AtomicU32>,
    mask: usize,
    slots: usize,
}

impl CasKmerTable {
    fn new(min_capacity: usize) -> Self {
        let slots = min_capacity.next_power_of_two();
        let keys = (0..slots).map(|_| AtomicU64::new(0)).collect();
        let counts = (0..slots).map(|_| AtomicU32::new(0)).collect();
        Self {
            keys,
            counts,
            mask: slots - 1,
            slots,
        }
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

    /// Insert a k-mer (or increment its count). Returns false only if
    /// the table is full (reprobe limit exceeded).
    #[inline]
    fn insert(&self, kmer: u64) -> bool {
        let pos = Self::hash(kmer) as usize & self.mask;
        let stored_key = kmer | OCCUPIED_BIT;

        for reprobe in 0..MAX_REPROBE {
            let idx = (pos + reprobe) & self.mask;
            let cell = &self.keys[idx];
            let mut current = cell.load(Ordering::Relaxed);

            loop {
                if current == 0 {
                    // Empty slot — try to claim it
                    match cell.compare_exchange_weak(
                        0,
                        stored_key,
                        Ordering::Relaxed,
                        Ordering::Relaxed,
                    ) {
                        Ok(_) => {
                            self.counts[idx].fetch_add(1, Ordering::Relaxed);
                            return true;
                        }
                        Err(actual) => {
                            current = actual;
                            continue; // Re-check: might be our key or different
                        }
                    }
                } else if current == stored_key {
                    // Same key — increment count
                    self.counts[idx].fetch_add(1, Ordering::Relaxed);
                    return true;
                } else {
                    break; // Different key — move to next probe position
                }
            }
        }
        false
    }

    fn occupied(&self) -> usize {
        self.keys
            .iter()
            .filter(|k| k.load(Ordering::Relaxed) != 0)
            .count()
    }

    fn iter(&self) -> impl Iterator<Item = (u64, u32)> + '_ {
        self.keys
            .iter()
            .zip(self.counts.iter())
            .filter_map(|(k, c)| {
                let stored = k.load(Ordering::Relaxed);
                if stored == 0 {
                    return None;
                }
                let kmer = stored & !OCCUPIED_BIT;
                let count = c.load(Ordering::Relaxed);
                Some((kmer, count))
            })
    }

    fn memory_bytes(&self) -> usize {
        self.slots * 12 // 8 bytes key + 4 bytes count
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
    table: Arc<CasKmerTable>,
    seq_buf: Vec<u8>,
}

#[inline]
fn count_kmers_in_seq(seq: &[u8], k: usize, canonical: bool, table: &CasKmerTable) {
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
                panic!(
                    "Hash table full (reprobe limit reached at {} slots). \
                     Re-run with a larger -s value.",
                    table.slots
                );
            }
        }
    }
}

impl ParallelProcessor for KmerCounter {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        self.seq_buf.clear();
        record.decode_s(&mut self.seq_buf)?;
        count_kmers_in_seq(&self.seq_buf, self.k, self.canonical, &self.table);

        if record.is_paired() {
            self.seq_buf.clear();
            record.decode_x(&mut self.seq_buf)?;
            count_kmers_in_seq(&self.seq_buf, self.k, self.canonical, &self.table);
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        Ok(()) // No-op: k-mers inserted directly via CAS
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.k == 0 || cli.k > 31 {
        anyhow::bail!("k must be between 1 and 31, got {}", cli.k);
    }

    let table_size = parse_size(&cli.s)?;
    let table = Arc::new(CasKmerTable::new(table_size));

    eprintln!(
        "Hash table: {} slots ({:.1} GB)",
        table.slots,
        table.memory_bytes() as f64 / 1e9
    );

    let reader = BinseqReader::new(&cli.input)?;
    let canonical = !cli.no_canonical;

    let counter = KmerCounter {
        k: cli.k,
        canonical,
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
        for (kmer, count) in table.iter() {
            buf.clear();
            bitnuc::from_2bit(kmer, cli.k, &mut buf)?;
            let seq = std::str::from_utf8(&buf)?;
            writeln!(out, "{}\t{}", seq, count)?;
        }
    } else {
        // Build histogram: frequency -> number of k-mers with that frequency
        let mut histogram: HashMap<u32, u64> = HashMap::new();
        for (_, count) in table.iter() {
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
    fn test_cas_table_basic() {
        let table = CasKmerTable::new(1024);
        assert!(table.insert(42));
        assert!(table.insert(42));
        assert!(table.insert(99));

        let mut found = HashMap::new();
        for (kmer, count) in table.iter() {
            found.insert(kmer, count);
        }
        assert_eq!(found[&42], 2);
        assert_eq!(found[&99], 1);
        assert_eq!(table.occupied(), 2);
    }

    #[test]
    fn test_cas_table_zero_kmer() {
        // k-mer 0 (all A's) must not collide with empty sentinel
        let table = CasKmerTable::new(1024);
        assert!(table.insert(0));
        assert!(table.insert(0));

        let mut found = HashMap::new();
        for (kmer, count) in table.iter() {
            found.insert(kmer, count);
        }
        assert_eq!(found[&0], 2);
        assert_eq!(table.occupied(), 1);
    }

    #[test]
    fn test_cas_table_many_kmers() {
        let table = CasKmerTable::new(4096);
        for i in 0..1000u64 {
            assert!(table.insert(i));
        }
        assert_eq!(table.occupied(), 1000);
    }
}
