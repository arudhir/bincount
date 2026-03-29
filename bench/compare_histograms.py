#!/usr/bin/env python3
"""Compare bqcount and jellyfish histogram files for accuracy."""

import sys


def load_bqcount_hist(path):
    """Load bqcount histogram (TSV with header: frequency\\tcount)."""
    hist = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            freq, count = line.strip().split("\t")
            hist[int(freq)] = int(count)
    return hist


def load_jellyfish_hist(path):
    """Load jellyfish histogram (space-delimited, no header)."""
    hist = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                hist[int(parts[0])] = int(parts[1])
    return hist


def compare(bq_path, jf_path):
    bq = load_bqcount_hist(bq_path)
    jf = load_jellyfish_hist(jf_path)

    all_freqs = set(bq.keys()) | set(jf.keys())
    matching = sum(1 for f in all_freqs if bq.get(f, 0) == jf.get(f, 0))
    total = len(all_freqs)

    bq_unique = sum(bq.values())
    jf_unique = sum(jf.values())

    accuracy = matching / total * 100 if total > 0 else 100.0

    return {
        "matching_bins": matching,
        "total_bins": total,
        "accuracy_pct": round(accuracy, 2),
        "bq_unique": bq_unique,
        "jf_unique": jf_unique,
    }


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <bqcount_hist.tsv> <jellyfish_hist.tsv>", file=sys.stderr)
        sys.exit(1)

    result = compare(sys.argv[1], sys.argv[2])
    # Print tab-separated: matching_bins total_bins accuracy_pct bq_unique jf_unique
    print(
        f"{result['matching_bins']}\t{result['total_bins']}\t{result['accuracy_pct']}\t"
        f"{result['bq_unique']}\t{result['jf_unique']}"
    )
