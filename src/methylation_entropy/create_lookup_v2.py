"""Version 2: generate minimum methylation entropy lookup tables."""

import argparse
import json
from pathlib import Path

import pandas as pd

from .entropy_simulation_v2 import calculate_methylation_entropy


def find_minimum_entropy_by_count(total_meth_cpgs, n_reads=16, n_cpgs=4):
    """Find the lowest entropy allocation for a methylated CpG count."""
    if not 0 <= total_meth_cpgs <= n_reads * n_cpgs:
        raise ValueError("methylated CpG count is outside the possible range")
    minimum, best = float("inf"), None

    def search(position, remaining_cpgs, remaining_reads, counts):
        nonlocal minimum, best
        if position == 0:
            if remaining_cpgs == 0:
                candidate = counts + [remaining_reads]
                entropy = calculate_methylation_entropy(dict(enumerate(candidate)), n_cpgs)
                if entropy < minimum:
                    minimum, best = entropy, candidate
            return
        if remaining_cpgs < 0 or remaining_cpgs > remaining_reads * position:
            return
        limit = min(remaining_reads, remaining_cpgs // position)
        for count in range(limit + 1):
            search(position - 1, remaining_cpgs - count * position,
                   remaining_reads - count, counts + [count])

    search(n_cpgs, total_meth_cpgs, n_reads, [])
    return minimum, best


def create_lookup_table(output_dir, n_reads=16, n_cpgs=4):
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    total, lookup, combinations = n_reads * n_cpgs, {}, {}
    for methylated in range(total + 1):
        key = f"{methylated / total:.3f}"
        entropy, counts = find_minimum_entropy_by_count(methylated, n_reads, n_cpgs)
        lookup[key] = entropy
        combinations[key] = {"meth_cpgs": methylated, **{
            f"n{n_cpgs - index}_reads": count for index, count in enumerate(counts)}}
    (destination / "minimum_entropy_lookup.json").write_text(json.dumps(lookup), encoding="utf-8")
    (destination / "minimum_entropy_combinations.json").write_text(
        json.dumps(combinations), encoding="utf-8")
    pd.DataFrame(lookup.items(), columns=["methylation", "min_entropy"]).to_csv(
        destination / "minimum_entropy_lookup.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description="Create a minimum entropy lookup table")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--nreads", type=int, default=16)
    parser.add_argument("--ncpgs", type=int, default=4)
    args = parser.parse_args()
    create_lookup_table(args.output_dir, args.nreads, args.ncpgs)


if __name__ == "__main__":
    main()
