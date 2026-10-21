"""Time cached and varied-depth RLM workloads using synthetic pattern counts."""

import argparse
import tempfile
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd

from methylation_entropy.rlm_entropy_v3 import (
    PATTERNS, PATTERN_WEIGHTS, REQUIRED_COLUMNS, entropy_from_counts,
    process_rlm_file,
)


def unique_depth_count_pairs(data):
    counts = data.loc[:, PATTERNS].to_numpy(dtype=np.int64)
    return len(set(zip(counts.sum(axis=1), counts @ PATTERN_WEIGHTS)))


def varied_rows(number, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for index in range(number):
        depth = int(rng.integers(10, 81))
        counts = rng.multinomial(depth, np.full(16, 1 / 16))
        frequencies = counts / depth
        row = {"#chr": "chr1", "start": 3 * index, "end": 3 * index + 2,
               "entropy": entropy_from_counts(counts),
               "epipolymorphism": 1 - float(np.sum(frequencies ** 2)),
               "mean_methylation": int(counts @ PATTERN_WEIGHTS) / (4 * depth),
               "coverage": depth}
        row.update(zip(PATTERNS, counts))
        rows.append(row)
    return pd.DataFrame(rows, columns=REQUIRED_COLUMNS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeated-rows", type=int, default=50000)
    parser.add_argument("--varied-rows", type=int, default=500)
    parser.add_argument("--nsims", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=71)
    args = parser.parse_args()
    if min(args.repeated_rows, args.varied_rows, args.nsims) < 1:
        parser.error("row counts and nsims must be positive")

    example = pd.read_csv(Path(__file__).with_name("rlm_entropy_demo.bed"), sep="\t")
    repeated = example.iloc[np.arange(args.repeated_rows) % len(example)].copy()
    repeated["start"] = 3 * np.arange(len(repeated))
    repeated["end"] = repeated["start"] + 2
    varied = varied_rows(args.varied_rows, args.seed)
    with tempfile.TemporaryDirectory() as directory:
        for name, data in (("repeated", repeated), ("varied", varied)):
            source = Path(directory) / (name + ".bed")
            output = Path(directory) / (name + ".tsv")
            data.to_csv(source, sep="\t", index=False)
            start = perf_counter()
            process_rlm_file(source, output, n_sims=args.nsims, seed=args.seed)
            elapsed = perf_counter() - start
            print(f"{name}: {len(data)} rows, {unique_depth_count_pairs(data)} "
                  f"distinct depth/count pairs, {elapsed:.3f} seconds")


if __name__ == "__main__":
    main()
