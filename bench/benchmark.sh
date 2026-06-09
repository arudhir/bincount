#!/bin/zsh
set -euo pipefail

# ============================================================================
# bqcount vs Jellyfish Benchmark
# ============================================================================

BENCH_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$BENCH_DIR")"
source "$HOME/.cargo/env"

BQCOUNT="$PROJECT_DIR/target/release/bqcount"
BQTOOLS="$(which bqtools)"
JELLYFISH="$(which jellyfish)"
WGSIM="$(which wgsim)"
TIME="/usr/bin/time"

GENOMES_DIR="$BENCH_DIR/genomes"
DATA_DIR="$BENCH_DIR/data"
RESULTS_DIR="$BENCH_DIR/results"
RESULTS_TSV="$RESULTS_DIR/results.tsv"
ACCURACY_TSV="$RESULTS_DIR/accuracy.tsv"

READ_LEN=150
WGSIM_SEED=42
WGSIM_MUT_RATE=0.001

# Genome URLs (NCBI FTP)
typeset -A GENOME_URLS
GENOME_URLS=(
    mflorum "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/008/305/GCA_000008305.1_ASM830v1/GCA_000008305.1_ASM830v1_genomic.fna.gz"
    ecoli   "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    yeast   "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
)

# ============================================================================
# Helper Functions
# ============================================================================

log() { echo "[$(date +%H:%M:%S)] $*" >&2; }

get_genome_size() {
    grep -v "^>" "$1" | tr -d '\n' | wc -c | tr -d ' '
}

compute_read_pairs() {
    local gsize=$1 coverage=$2
    echo $(( gsize * coverage / READ_LEN ))
}

parse_timing() {
    # Parse macOS /usr/bin/time -l output
    # Returns: wall_sec cpu_sec peak_rss_mb
    local file=$1
    local real user sys rss

    # macOS time format: "  3.39 real  26.50 user  0.13 sys" (all on one line)
    local timeline
    timeline=$(grep 'real.*user.*sys' "$file" || true)
    real=$(echo "$timeline" | awk '{print $1}')
    user=$(echo "$timeline" | awk '{print $3}')
    sys=$(echo "$timeline" | awk '{print $5}')
    rss=$(grep 'maximum resident set size' "$file" | awk '{print $1}' || true)

    # cpu_sec = user + sys
    local cpu_sec
    cpu_sec=$(python3 -c "print(round(${user:-0} + ${sys:-0}, 2))")
    # rss bytes -> MB
    local rss_mb
    rss_mb=$(python3 -c "print(round(${rss:-0} / 1048576, 1))")

    echo "${real:-0} ${cpu_sec} ${rss_mb}"
}

extract_unique_kmers_bq() {
    # Extract unique k-mer count from bqcount stderr captured in timing file
    grep 'unique k-mers' "$1" | grep -oE '[0-9]+' | head -1
}

extract_unique_kmers_jf() {
    # Sum all counts in jellyfish histogram
    awk '{s+=$2} END {print s}' "$1"
}

# ============================================================================
# Download Genomes
# ============================================================================

download_genomes() {
    log "=== Downloading genomes ==="
    mkdir -p "$GENOMES_DIR"

    for org in mflorum ecoli yeast; do
        local fa="$GENOMES_DIR/${org}.fa"
        if [[ -f "$fa" ]]; then
            log "  $org: already exists ($(get_genome_size "$fa") bp)"
            continue
        fi
        log "  $org: downloading..."
        curl -sL "${GENOME_URLS[$org]}" | gunzip > "$fa"
        log "  $org: $(get_genome_size "$fa") bp"
    done
}

# ============================================================================
# Generate Dataset
# ============================================================================

generate_dataset() {
    local org=$1 error=$2 coverage=$3
    local genome="$GENOMES_DIR/${org}.fa"
    local datadir="$DATA_DIR/${org}/e${error}/c${coverage}"
    local bqfile="$datadir/reads.bq"

    if [[ -f "$bqfile" ]]; then
        log "  Dataset ${org}/e${error}/c${coverage}: .bq exists, skipping"
        return
    fi

    mkdir -p "$datadir"
    local gsize
    gsize=$(get_genome_size "$genome")
    local npairs
    npairs=$(compute_read_pairs "$gsize" "$coverage")

    log "  Generating ${org}/e${error}/c${coverage}: ${npairs} pairs..."
    $WGSIM -e "$error" -r "$WGSIM_MUT_RATE" -N "$npairs" \
        -1 "$READ_LEN" -2 "$READ_LEN" -S "$WGSIM_SEED" \
        "$genome" "$datadir/R1.fq" "$datadir/R2.fq" > /dev/null 2>&1

    log "  Encoding to .bq..."
    $BQTOOLS encode "$datadir/R1.fq" "$datadir/R2.fq" -o "$bqfile" 2>/dev/null

    log "  Compressing FASTQ..."
    gzip "$datadir/R1.fq"
    gzip "$datadir/R2.fq"
}

generate_meta_mix() {
    local ratio=$1 error=$2  # ratio like "90_10"
    local ec_frac=${ratio%%_*}
    local yt_frac=${ratio##*_}
    local datadir="$DATA_DIR/meta_${ratio}/e${error}/c100"
    local bqfile="$datadir/reads.bq"

    if [[ -f "$bqfile" ]]; then
        log "  Meta mix ${ratio}/e${error}: .bq exists, skipping"
        return
    fi

    mkdir -p "$datadir"

    local ec_gsize yt_gsize combined_size total_pairs ec_pairs yt_pairs
    ec_gsize=$(get_genome_size "$GENOMES_DIR/ecoli.fa")
    yt_gsize=$(get_genome_size "$GENOMES_DIR/yeast.fa")
    combined_size=$(( ec_gsize + yt_gsize ))
    total_pairs=$(( combined_size * 100 / READ_LEN ))
    ec_pairs=$(( total_pairs * ec_frac / (ec_frac + yt_frac) ))
    yt_pairs=$(( total_pairs * yt_frac / (ec_frac + yt_frac) ))

    log "  Meta ${ratio}/e${error}: E.coli ${ec_pairs} + yeast ${yt_pairs} pairs..."

    $WGSIM -e "$error" -r "$WGSIM_MUT_RATE" -N "$ec_pairs" \
        -1 "$READ_LEN" -2 "$READ_LEN" -S "$WGSIM_SEED" \
        "$GENOMES_DIR/ecoli.fa" "$datadir/ec_R1.fq" "$datadir/ec_R2.fq" > /dev/null 2>&1

    $WGSIM -e "$error" -r "$WGSIM_MUT_RATE" -N "$yt_pairs" \
        -1 "$READ_LEN" -2 "$READ_LEN" -S "$WGSIM_SEED" \
        "$GENOMES_DIR/yeast.fa" "$datadir/yt_R1.fq" "$datadir/yt_R2.fq" > /dev/null 2>&1

    cat "$datadir/ec_R1.fq" "$datadir/yt_R1.fq" > "$datadir/R1.fq"
    cat "$datadir/ec_R2.fq" "$datadir/yt_R2.fq" > "$datadir/R2.fq"
    rm "$datadir/ec_R1.fq" "$datadir/yt_R1.fq" "$datadir/ec_R2.fq" "$datadir/yt_R2.fq"

    $BQTOOLS encode "$datadir/R1.fq" "$datadir/R2.fq" -o "$bqfile" 2>/dev/null

    gzip "$datadir/R1.fq"
    gzip "$datadir/R2.fq"
}

# ============================================================================
# Run Tools
# ============================================================================

run_jellyfish() {
    local sweep=$1 org=$2 error=$3 coverage=$4 k=$5 threads=$6
    local datadir="$DATA_DIR/${org}/e${error}/c${coverage}"
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local jf_out="$datadir/counts_${tag}.jf"
    local jf_hist="$RESULTS_DIR/${tag}_jf_hist.tsv"
    local timing_file="$RESULTS_DIR/${tag}_jf_timing.txt"

    if [[ -f "$jf_hist" ]]; then
        log "  Jellyfish ${tag}: histogram exists, skipping"
        return
    fi

    local gsize hash_size
    if [[ "$org" == meta_* ]]; then
        local ec_size yt_size
        ec_size=$(get_genome_size "$GENOMES_DIR/ecoli.fa")
        yt_size=$(get_genome_size "$GENOMES_DIR/yeast.fa")
        gsize=$(( ec_size + yt_size ))
    else
        gsize=$(get_genome_size "$GENOMES_DIR/${org}.fa")
    fi
    hash_size=$(( gsize * 4 ))

    log "  Jellyfish count ${tag}..."
    { $TIME -l $JELLYFISH count -m "$k" -s "$hash_size" -t "$threads" -C \
        -o "$jf_out" \
        <(gunzip -c "$datadir/R1.fq.gz") <(gunzip -c "$datadir/R2.fq.gz") ; } 2> "$timing_file"

    $JELLYFISH histo -h 100000 "$jf_out" > "$jf_hist"

    local unique_kmers
    unique_kmers=$(extract_unique_kmers_jf "$jf_hist")

    local timing
    timing=$(parse_timing "$timing_file")
    local wall_sec cpu_sec peak_rss_mb
    read -r wall_sec cpu_sec peak_rss_mb <<< "$timing"

    echo -e "${sweep}\t${org}\t${error}\t${coverage}\t${k}\t${threads}\tjellyfish\t${wall_sec}\t${cpu_sec}\t${peak_rss_mb}\t${unique_kmers}" >> "$RESULTS_TSV"

    # Clean up .jf file
    rm -f "$jf_out"
}

run_bqcount() {
    local sweep=$1 org=$2 error=$3 coverage=$4 k=$5 threads=$6
    local datadir="$DATA_DIR/${org}/e${error}/c${coverage}"
    local tag="${org}_e${error}_c${coverage}_k${k}_t${threads}"
    local bq_hist="$RESULTS_DIR/${tag}_bq_hist.tsv"
    local timing_file="$RESULTS_DIR/${tag}_bq_timing.txt"

    if [[ -f "$bq_hist" ]]; then
        log "  bqcount ${tag}: histogram exists, skipping"
        return
    fi

    # Auto-size the hash table: bqcount's open-addressing table panics
    # ("Hash table full") near ~76% load on the 128M default. Size to 2x the
    # unique k-mers from the jellyfish histogram (already produced for this
    # config), rounded up to a standard size. Falls back to 256M.
    local jf_hist="$RESULTS_DIR/${tag}_jf_hist.tsv"
    local table_size="256M"
    if [[ -f "$jf_hist" && -s "$jf_hist" ]]; then
        local uniq
        uniq=$(awk '{s+=$2} END {print s}' "$jf_hist")
        table_size=$(python3 -c "
u = int(${uniq:-0})
needed = max(u * 2, 128_000_000)
for m, label in [(128,'128M'),(256,'256M'),(512,'512M'),(1024,'1G'),(2048,'2G'),(4096,'4G'),(8192,'8G')]:
    if needed <= m * 1_000_000:
        print(label); break
else:
    print('8G')
")
    fi

    log "  bqcount ${tag} (table: ${table_size})..."
    { $TIME -l $BQCOUNT "$datadir/reads.bq" -k "$k" -t "$threads" -s "$table_size" -o "$bq_hist" ; } 2> "$timing_file"

    local unique_kmers
    unique_kmers=$(extract_unique_kmers_bq "$timing_file")

    local timing
    timing=$(parse_timing "$timing_file")
    local wall_sec cpu_sec peak_rss_mb
    read -r wall_sec cpu_sec peak_rss_mb <<< "$timing"

    echo -e "${sweep}\t${org}\t${error}\t${coverage}\t${k}\t${threads}\tbqcount\t${wall_sec}\t${cpu_sec}\t${peak_rss_mb}\t${unique_kmers}" >> "$RESULTS_TSV"
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
    run_jellyfish "$sweep" "$org" "$error" "$coverage" "$k" "$threads"
    run_bqcount "$sweep" "$org" "$error" "$coverage" "$k" "$threads"
    compare_accuracy "$org" "$error" "$coverage" "$k" "$threads"
}

cleanup_fastq() {
    # Delete FASTQ.gz for a dataset (called after jellyfish is done with it)
    local datadir=$1
    rm -f "$datadir/R1.fq.gz" "$datadir/R2.fq.gz"
}

# ============================================================================
# Main
# ============================================================================

main() {
    log "=== bqcount vs Jellyfish Benchmark ==="
    log "Project: $PROJECT_DIR"
    log "Bench:   $BENCH_DIR"

    # Rebuild bqcount
    log "=== Building bqcount (release) ==="
    (cd "$PROJECT_DIR" && RUSTFLAGS="-C target-cpu=native" cargo build --release 2>&1 | tail -1)

    # Download genomes
    download_genomes

    # Init results files
    mkdir -p "$RESULTS_DIR"
    echo -e "sweep\torganism\terror_rate\tcoverage\tk\tthreads\ttool\twall_sec\tcpu_sec\tpeak_rss_mb\tunique_kmers" > "$RESULTS_TSV"
    echo -e "organism\terror_rate\tcoverage\tk\tthreads\tmatching_bins\ttotal_bins\taccuracy_pct\tbq_unique\tjf_unique" > "$ACCURACY_TSV"

    # ------------------------------------------------------------------
    # Sweep 1: Primary (organism × error × coverage)
    # ------------------------------------------------------------------
    log "=== Sweep 1: Primary ==="

    local ORGANISMS=(mflorum ecoli yeast)
    local ERRORS=(0.001 0.005 0.02)
    local COVERAGES=(10 50 100)

    # Generate all single-organism datasets
    for org in "${ORGANISMS[@]}"; do
        for err in "${ERRORS[@]}"; do
            for cov in "${COVERAGES[@]}"; do
                generate_dataset "$org" "$err" "$cov"
            done
        done
    done

    # Generate metagenomic mixes (100x only)
    for ratio in 90_10 50_50; do
        for err in "${ERRORS[@]}"; do
            generate_meta_mix "$ratio" "$err"
        done
    done

    log "=== Running primary benchmarks ==="

    # Run benchmarks for single organisms
    for org in "${ORGANISMS[@]}"; do
        for err in "${ERRORS[@]}"; do
            for cov in "${COVERAGES[@]}"; do
                local datadir="$DATA_DIR/${org}/e${err}/c${cov}"
                run_config "primary" "$org" "$err" "$cov" 21 8

                # Clean FASTQ.gz unless this is the sweep 2/3 dataset
                if [[ "$org" != "ecoli" ]] || [[ "$err" != "0.005" ]] || [[ "$cov" != "100" ]]; then
                    cleanup_fastq "$datadir"
                fi
            done
        done
    done

    # Run benchmarks for meta mixes
    for ratio in 90_10 50_50; do
        for err in "${ERRORS[@]}"; do
            local datadir="$DATA_DIR/meta_${ratio}/e${err}/c100"
            run_config "primary" "meta_${ratio}" "$err" 100 21 8
            cleanup_fastq "$datadir"
        done
    done

    # ------------------------------------------------------------------
    # Sweep 2: Thread scaling (ecoli / e0.005 / c100)
    # ------------------------------------------------------------------
    log "=== Sweep 2: Thread scaling ==="
    for threads in 1 2 4 8; do
        run_config "threads" "ecoli" "0.005" 100 21 "$threads"
    done

    # ------------------------------------------------------------------
    # Sweep 3: K-mer size (ecoli / e0.005 / c100)
    # ------------------------------------------------------------------
    log "=== Sweep 3: K-mer size ==="
    for k in 15 21 31; do
        run_config "kmer" "ecoli" "0.005" 100 "$k" 8
    done

    # Clean the sweep dataset FASTQ.gz now
    cleanup_fastq "$DATA_DIR/ecoli/e0.005/c100"

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    log "=== Benchmark complete ==="
    log "Results: $RESULTS_TSV"
    log "Accuracy: $ACCURACY_TSV"

    local n_results n_accuracy
    n_results=$(( $(wc -l < "$RESULTS_TSV") - 1 ))
    n_accuracy=$(( $(wc -l < "$ACCURACY_TSV") - 1 ))
    log "Total benchmark runs: $n_results"
    log "Total accuracy comparisons: $n_accuracy"

    # Print disk usage
    log "Disk usage:"
    du -sh "$DATA_DIR" "$RESULTS_DIR" 2>/dev/null | while read -r line; do log "  $line"; done
}

main "$@"
