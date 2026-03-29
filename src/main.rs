use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::Result;
use binseq::prelude::*;
use clap::Parser;
use hashbrown::HashMap;
use parking_lot::Mutex;

/// A binseq-native parallel k-mer counter.
///
/// Counts k-mers from BINSEQ files (.bq, .cbq) and outputs a frequency histogram.
#[derive(Parser)]
#[command(name = "bqcount", version)]
struct Cli {
    /// Path to .bq or .cbq file
    input: String,

    /// K-mer size (max: 32)
    #[arg(short, default_value_t = 21)]
    k: usize,

    /// Number of threads (0 = auto-detect)
    #[arg(short, default_value_t = 0)]
    t: usize,

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

/// Reverse complement of a 2-bit packed k-mer.
///
/// Encoding: A=00, C=01, G=10, T=11 (LSB-first).
/// Complement: NOT each 2-bit pair (A<->T, C<->G).
/// Reverse: swap order of 2-bit pairs.
#[inline]
fn revcomp_2bit(packed: u64, k: usize) -> u64 {
    // Complement all bits
    let comp = !packed;
    // Reverse 2-bit pairs within the u64
    let mut rev = comp;
    // Swap adjacent 2-bit groups
    rev = ((rev & 0x3333_3333_3333_3333) << 2) | ((rev >> 2) & 0x3333_3333_3333_3333);
    // Swap adjacent 4-bit groups
    rev = ((rev & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((rev >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
    // Reverse bytes
    rev = rev.swap_bytes();
    // Shift right to align k bases at LSB
    rev >> (64 - 2 * k)
}

/// Canonical k-mer: min(forward, reverse complement)
#[inline]
fn canonical_kmer(packed: u64, k: usize) -> u64 {
    let rc = revcomp_2bit(packed, k);
    packed.min(rc)
}

#[derive(Clone)]
struct KmerCounter {
    k: usize,
    canonical: bool,
    local_counts: HashMap<u64, u64>,
    global_counts: Arc<Mutex<HashMap<u64, u64>>>,
    seq_buf: Vec<u8>,
}

impl KmerCounter {
    fn new(k: usize, canonical: bool) -> Self {
        Self {
            k,
            canonical,
            local_counts: HashMap::new(),
            global_counts: Arc::new(Mutex::new(HashMap::new())),
            seq_buf: Vec::new(),
        }
    }

}

#[inline]
fn count_kmers_in_seq(seq: &[u8], k: usize, canonical: bool, counts: &mut HashMap<u64, u64>) {
    if seq.len() < k {
        return;
    }
    for window in seq.windows(k) {
        // Skip windows containing N or other invalid bases
        if let Ok(packed) = bitnuc::as_2bit(window) {
            let key = if canonical {
                canonical_kmer(packed, k)
            } else {
                packed
            };
            *counts.entry(key).or_insert(0) += 1;
        }
    }
}

impl ParallelProcessor for KmerCounter {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        self.seq_buf.clear();
        record.decode_s(&mut self.seq_buf)?;
        count_kmers_in_seq(&self.seq_buf, self.k, self.canonical, &mut self.local_counts);

        if record.is_paired() {
            self.seq_buf.clear();
            record.decode_x(&mut self.seq_buf)?;
            count_kmers_in_seq(&self.seq_buf, self.k, self.canonical, &mut self.local_counts);
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        if self.local_counts.is_empty() {
            return Ok(());
        }
        let mut global = self.global_counts.lock();
        for (kmer, count) in self.local_counts.drain() {
            *global.entry(kmer).or_insert(0) += count;
        }
        Ok(())
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.k == 0 || cli.k > 32 {
        anyhow::bail!("k must be between 1 and 32, got {}", cli.k);
    }

    let reader = BinseqReader::new(&cli.input)?;
    let canonical = !cli.no_canonical;
    let counter = KmerCounter::new(cli.k, canonical);

    eprintln!(
        "Counting {}-mers{} from {} using {} thread(s)...",
        cli.k,
        if canonical { " (canonical)" } else { "" },
        cli.input,
        if cli.t == 0 { "auto".to_string() } else { cli.t.to_string() },
    );

    reader.process_parallel(counter.clone(), cli.t)?;

    let counts = counter.global_counts.lock();
    eprintln!("Found {} unique k-mers", counts.len());

    // Output
    let out: Box<dyn Write> = match &cli.o {
        Some(path) => Box::new(BufWriter::new(std::fs::File::create(path)?)),
        None => Box::new(BufWriter::new(io::stdout().lock())),
    };
    let mut out = out;

    if cli.dump {
        // Dump full k-mer count table
        writeln!(out, "kmer\tcount")?;
        let mut buf = Vec::new();
        for (&kmer, &count) in counts.iter() {
            buf.clear();
            bitnuc::from_2bit(kmer, cli.k, &mut buf)?;
            let seq = std::str::from_utf8(&buf)?;
            writeln!(out, "{}\t{}", seq, count)?;
        }
    } else {
        // Build histogram: frequency -> number of k-mers with that frequency
        let mut histogram: HashMap<u64, u64> = HashMap::new();
        for &count in counts.values() {
            *histogram.entry(count).or_insert(0) += 1;
        }

        // Sort by frequency and output
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
        // A=00, C=01, G=10, T=11
        // "ACG" packed: G(10) C(01) A(00) = 0b10_01_00 = 0x24
        let acg = bitnuc::as_2bit(b"ACG").unwrap();
        let rc = revcomp_2bit(acg, 3);
        // revcomp("ACG") = "CGT"
        let cgt = bitnuc::as_2bit(b"CGT").unwrap();
        assert_eq!(rc, cgt, "revcomp(ACG) should be CGT");
    }

    #[test]
    fn test_revcomp_palindrome() {
        // "ACGT" is a palindrome: revcomp("ACGT") = "ACGT"
        let acgt = bitnuc::as_2bit(b"ACGT").unwrap();
        let rc = revcomp_2bit(acgt, 4);
        assert_eq!(rc, acgt, "ACGT should be its own reverse complement");
    }

    #[test]
    fn test_canonical_kmer() {
        let acg = bitnuc::as_2bit(b"ACG").unwrap();
        let cgt = bitnuc::as_2bit(b"CGT").unwrap();
        // canonical should give the same result for a k-mer and its revcomp
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
}
