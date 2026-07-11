"""Version 2: portable and reproducible methylation entropy simulation."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from collections import Counter
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

COLUMNS = ["chrom", "start", "end", "entropy", "methylation"]
OUTPUT_COLUMNS = ["chrom", "start", "end", "observed_entropy", "methylation",
                  "allelic_entropy", "stochastic_mean", "stochastic_min", "stochastic_max"]


def calculate_methylation_entropy(pattern_counts: Mapping[object, int], n_cpgs: int) -> float:
    """Calculate Shannon entropy per CpG."""
    if n_cpgs <= 0:
        raise ValueError("n_cpgs must be positive")
    total = sum(pattern_counts.values())
    if total == 0:
        return 0.0
    probabilities = np.array([count / total for count in pattern_counts.values() if count > 0])
    return float(-(probabilities * np.log2(probabilities)).sum() / n_cpgs)


def get_minimum_entropy(meth_level: float, lookup: Mapping[str, float]) -> float:
    """Return the entropy at the nearest available methylation level."""
    if not lookup:
        raise ValueError("lookup table is empty")
    key = min(lookup, key=lambda item: abs(float(item) - meth_level))
    return float(lookup[key])


def simulate_entropy_ranges(meth_level: float, n_reads: int = 128, n_cpgs: int = 4,
                            n_sims: int = 1000, rng=None) -> dict[str, float]:
    """Simulate entropy while preserving the rounded methylated-site count."""
    if not 0 <= meth_level <= 1:
        raise ValueError("meth_level must be between 0 and 1")
    if min(n_reads, n_cpgs, n_sims) <= 0:
        raise ValueError("n_reads, n_cpgs, and n_sims must be positive")
    rng = rng or np.random.default_rng()
    total = n_reads * n_cpgs
    methylated = round(total * meth_level)
    sites = np.r_[np.ones(methylated, dtype=np.uint8),
                  np.zeros(total - methylated, dtype=np.uint8)]
    entropies = np.empty(n_sims)
    for index in range(n_sims):
        patterns = rng.permutation(sites).reshape(n_reads, n_cpgs)
        entropies[index] = calculate_methylation_entropy(Counter(map(tuple, patterns)), n_cpgs)
    return {"stochastic_mean": float(np.mean(entropies)),
            "stochastic_min": float(np.percentile(entropies, 5)),
            "stochastic_max": float(np.percentile(entropies, 95))}


def _process_chunk(chunk, lookup, n_sims, n_reads, n_cpgs, seed):
    rng, rows = np.random.default_rng(seed), []
    for row in chunk.itertuples(index=False):
        if pd.notna(row.methylation) and 0 <= row.methylation <= 1:
            result = {"chrom": row.chrom, "start": row.start, "end": row.end,
                      "observed_entropy": row.entropy, "methylation": row.methylation,
                      "allelic_entropy": get_minimum_entropy(row.methylation, lookup)}
            result.update(simulate_entropy_ranges(row.methylation, n_reads, n_cpgs, n_sims, rng))
            rows.append(result)
    return rows


def process_file(input_file, output_file, lookup_file, n_sims=1000, n_reads=128,
                 n_cpgs=4, workers=None, seed=None):
    """Process a BED-like input file and write a tab-separated result."""
    with Path(lookup_file).open(encoding="utf-8") as handle:
        lookup = json.load(handle)
    data = pd.read_csv(input_file, sep="\t", header=None, names=COLUMNS)
    if data.empty:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(output_file, sep="\t", index=False)
        return
    worker_count = max(1, min(workers or max(1, mp.cpu_count() - 1), len(data)))
    chunk_count = min(len(data), worker_count * 4)
    chunks = [data.iloc[indexes] for indexes in np.array_split(np.arange(len(data)), chunk_count)]
    seeds = [int(x.generate_state(1)[0]) for x in np.random.SeedSequence(seed).spawn(len(chunks))]
    if worker_count == 1:
        groups = [_process_chunk(chunk, lookup, n_sims, n_reads, n_cpgs, value)
                  for chunk, value in zip(chunks, seeds)]
    else:
        with mp.Pool(worker_count) as pool:
            arguments = [(chunk, lookup, n_sims, n_reads, n_cpgs, value)
                         for chunk, value in zip(chunks, seeds)]
            groups = pool.starmap(_process_chunk, arguments)
    pd.DataFrame([row for group in groups for row in group], columns=OUTPUT_COLUMNS).to_csv(
        output_file, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Calculate methylation entropy reference ranges")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--lookup", required=True)
    parser.add_argument("--nsims", type=int, default=1000)
    parser.add_argument("--nreads", type=int, default=128)
    parser.add_argument("--ncpgs", type=int, default=4)
    parser.add_argument("--workers", type=int)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    process_file(args.input, args.output, args.lookup, args.nsims, args.nreads,
                 args.ncpgs, args.workers, args.seed)


if __name__ == "__main__":
    main()
