#!/bin/zsh
set -euo pipefail

# ============================================================================
# Rerun bqcount benchmarks only, keeping existing jellyfish results.
# Use after freeing memory to get clean bqcount timings.
# ============================================================================

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$BENCH_DIR")"
source "$HOME/.cargo/env"

BQCOUNT="$PROJECT_DIR/target/release/bqcount"
TIME="/usr/bin/time"

DATA_DIR="$BENCH_DIR/data"
RESULTS_DIR="$BENCH_DIR/results"
RESULTS_TSV="$RESULTS_DIR/results.tsv"
ACCURACY_TSV="$RESULTS_DIR/accuracy.tsv"

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

parse_timing() {
    local file=$1
    local timeline
    timeline=$(grep 'real.*user.*sys' "$file" || true)
    local real=$(echo "$timeline" | awk '{print $1}')
    local user=$(echo "$timeline" | awk '{print $3}')
    local sys=$(echo "$timeline" | awk '{print $5}')
    local rss=$(grep 'maximum resident set size' "$file" | awk '{print $1}' || true)
    local cpu_sec=$(python3 -c "print(round(${user:-0} + ${sys:-0}, 2))")
    local rss_mb=$(python3 -c "print(round(${rss:-0} / 1048576, 1))")
    echo "${real:-0} ${cpu_sec} ${rss_mb}"
}

extract_unique_kmers_bq() {
    grep 'unique k-mers' "$1" | grep -oE '[0-9]+' | head -1
}

extract_unique_kmers_jf() {
    awk '{s+=$2} END {print s}' "$1"
}

# Compute table size needed: 2x the unique k-mers from jellyfish, rounded up to next power of 2.
# Falls back to 256M if jellyfish data is missing.
compute_table_size() {
    local org=$1 error=$2 coverage=$3 k=$4 threads=$5
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local jf_hist="$RESULTS_DIR/${tag}_jf_hist.tsv"

    if [[ ! -f "$jf_hist" ]] || [[ ! -s "$jf_hist" ]]; then
        echo "256M"
        return
    fi

    local unique_kmers
    unique_kmers=$(awk '{s+=$2} END {print s}' "$jf_hist")
    # 2x headroom, rounded up to standard sizes
    local needed=$(python3 -c "
u = int(${unique_kmers:-0})
needed = max(u * 2, 128_000_000)
# Standard sizes in ascending order
sizes = [(128,'128M'),(256,'256M'),(512,'512M'),(1024,'1G'),(2048,'2G'),(4096,'4G')]
for m, label in sizes:
    if needed <= m * 1_000_000:
        print(label)
        break
else:
    print('4G')
")
    echo "$needed"
}

# Re-run bqcount for a config, replacing old histogram and timing
run_bqcount() {
    local sweep=$1 org=$2 error=$3 coverage=$4 k=$5 threads=$6
    local datadir="$DATA_DIR/${org}/e${error}/c${coverage}"
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local bq_hist="$RESULTS_DIR/${tag}_bq_hist.tsv"
    local timing_file="$RESULTS_DIR/${tag}_bq_timing.txt"

    rm -f "$bq_hist" "$timing_file"

    local table_size
    table_size=$(compute_table_size "$org" "$error" "$coverage" "$k" "$threads")

    log "  bqcount ${tag} (table: ${table_size})..."
    { $TIME -l $BQCOUNT "$datadir/reads.bq" -k "$k" -t "$threads" -s "$table_size" -o "$bq_hist" ; } 2> "$timing_file"

    local unique_kmers=$(extract_unique_kmers_bq "$timing_file")
    local timing=$(parse_timing "$timing_file")
    local wall_sec cpu_sec peak_rss_mb
    read -r wall_sec cpu_sec peak_rss_mb <<< "$timing"

    echo -e "${sweep}\t${org}\t${error}\t${coverage}\t${k}\t${threads}\tbqcount\t${wall_sec}\t${cpu_sec}\t${peak_rss_mb}\t${unique_kmers}" >> "$RESULTS_TSV"
}

# Emit existing jellyfish result from saved timing file
emit_jellyfish() {
    local sweep=$1 org=$2 error=$3 coverage=$4 k=$5 threads=$6
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local timing_file="$RESULTS_DIR/${tag}_jf_timing.txt"
    local jf_hist="$RESULTS_DIR/${tag}_jf_hist.tsv"

    if [[ ! -f "$timing_file" ]] || [[ ! -f "$jf_hist" ]]; then
        log "  WARNING: missing jellyfish data for ${tag}, skipping"
        return
    fi

    local unique_kmers=$(extract_unique_kmers_jf "$jf_hist")
    local timing=$(parse_timing "$timing_file")
    local wall_sec cpu_sec peak_rss_mb
    read -r wall_sec cpu_sec peak_rss_mb <<< "$timing"

    echo -e "${sweep}\t${org}\t${error}\t${coverage}\t${k}\t${threads}\tjellyfish\t${wall_sec}\t${cpu_sec}\t${peak_rss_mb}\t${unique_kmers}" >> "$RESULTS_TSV"
}

compare_accuracy() {
    local org=$1 error=$2 coverage=$3 k=$4 threads=$5
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local bq_hist="$RESULTS_DIR/${tag}_bq_hist.tsv"
    local jf_hist="$RESULTS_DIR/${tag}_jf_hist.tsv"

    if [[ ! -f "$bq_hist" ]] || [[ ! -f "$jf_hist" ]]; then
        return
    fi

    local acc
    acc=$(python3 "$BENCH_DIR/compare_histograms.py" "$bq_hist" "$jf_hist")
    echo -e "${org}\t${error}\t${coverage}\t${k}\t${threads}\t${acc}" >> "$ACCURACY_TSV"
}

run_config() {
    local sweep=$1 org=$2 error=$3 coverage=$4 k=$5 threads=$6
    emit_jellyfish "$@"
    run_bqcount "$@"
    compare_accuracy "$org" "$error" "$coverage" "$k" "$threads"
}

main() {
    log "=== bqcount Rerun (keeping jellyfish results) ==="

    # Rebuild bqcount
    log "=== Building bqcount (release) ==="
    (cd "$PROJECT_DIR" && RUSTFLAGS="-C target-cpu=native" cargo build --release 2>&1 | tail -1)

    # Reset results files
    echo -e "sweep\torganism\terror_rate\tcoverage\tk\tthreads\ttool\twall_sec\tcpu_sec\tpeak_rss_mb\tunique_kmers" > "$RESULTS_TSV"
    echo -e "organism\terror_rate\tcoverage\tk\tthreads\tmatching_bins\ttotal_bins\taccuracy_pct\tbq_unique\tjf_unique" > "$ACCURACY_TSV"

    local ORGANISMS=(mflorum ecoli yeast)
    local ERRORS=(0.001 0.005 0.02)
    local COVERAGES=(10 50 100)

    # ------------------------------------------------------------------
    # Sweep 1: Primary (organism x error x coverage)
    # ------------------------------------------------------------------
    log "=== Sweep 1: Primary ==="
    for org in "${ORGANISMS[@]}"; do
        for err in "${ERRORS[@]}"; do
            for cov in "${COVERAGES[@]}"; do
                run_config "primary" "$org" "$err" "$cov" 21 8
            done
        done
    done

    for ratio in 90_10 50_50; do
        for err in "${ERRORS[@]}"; do
            run_config "primary" "meta_${ratio}" "$err" 100 21 8
        done
    done

    # ------------------------------------------------------------------
    # Sweep 2: Thread scaling (ecoli / e0.005 / c100 / k21)
    # ------------------------------------------------------------------
    log "=== Sweep 2: Thread scaling ==="
    for threads in 1 2 4; do
        run_config "threads" "ecoli" "0.005" 100 21 "$threads"
    done

    # ------------------------------------------------------------------
    # Sweep 3: K-mer size (ecoli / e0.005 / c100 / t8)
    # ------------------------------------------------------------------
    log "=== Sweep 3: K-mer size ==="
    for k in 15 31; do
        run_config "kmer" "ecoli" "0.005" 100 "$k" 8
    done

    local n_results=$(( $(wc -l < "$RESULTS_TSV") - 1 ))
    local n_accuracy=$(( $(wc -l < "$ACCURACY_TSV") - 1 ))
    log "=== Rerun complete: $n_results results, $n_accuracy accuracy comparisons ==="
}

main "$@"
