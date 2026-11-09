"""Coverage- and CpG-marginal-matched entropy analysis for RLM four-CpG windows."""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

from .rlm_entropy_v3 import PATTERNS, REQUIRED_COLUMNS, _validate_chunk


PATTERN_BITS = np.array([[int(base == "G") for base in pattern]
                         for pattern in PATTERNS], dtype=np.int64)
OUTPUT_COLUMNS = (
    "chrom", "start", "end", "coverage", "methylation", "observed_entropy",
    "informative_cpg_count",
    *(f"cpg{index}_methylated_count" for index in range(1, 5)),
    *(f"cpg{index}_methylation" for index in range(1, 5)),
    "independence_max_entropy", "observed_total_correlation",
    "independent_null_mean", "independent_null_sd", "independent_null_p05",
    "independent_null_p95",
    "entropy_minus_independent_null", "empirical_p_coordination",
    "empirical_p_dispersion", "n_sims",
)


def binary_entropy(probabilities: np.ndarray) -> np.ndarray:
    """Binary Shannon entropy, with zero at probabilities zero and one."""
    probabilities = np.asarray(probabilities, dtype=float)
    result = np.zeros_like(probabilities)
    middle = (probabilities > 0) & (probabilities < 1)
    p = probabilities[middle]
    result[middle] = -p * np.log2(p) - (1 - p) * np.log2(1 - p)
    return result


def simulate_independent_null(
    n_reads: int, methylated_by_cpg: tuple[int, int, int, int],
    n_sims: int, rng: np.random.Generator,
) -> np.ndarray:
    """Permute each CpG across molecules while fixing its methylated-call count.

    Every draw preserves coverage and all four observed CpG marginals. Batching
    limits temporary memory for high coverage or large simulation counts.
    """
    if (n_reads < 1 or n_sims < 1 or len(methylated_by_cpg) != 4
            or any(not 0 <= count <= n_reads for count in methylated_by_cpg)):
        raise ValueError("invalid coverage, CpG counts, or simulation count")
    entropies = np.empty(n_sims, dtype=float)
    batch_limit = max(1, min(1000, 250000 // n_reads))
    for start in range(0, n_sims, batch_limit):
        batch = min(batch_limit, n_sims - start)
        patterns = _draw_independent_patterns(n_reads, methylated_by_cpg,
                                              batch, rng)
        indexed = patterns.astype(np.int64) + 16 * np.arange(batch)[:, None]
        counts = np.bincount(indexed.ravel(), minlength=16 * batch).reshape(batch, 16)
        frequencies = counts / n_reads
        logarithms = np.zeros_like(frequencies)
        np.log2(frequencies, out=logarithms, where=counts > 0)
        entropies[start:start + batch] = -(frequencies * logarithms).sum(axis=1) / 4
    return np.sort(entropies)


def _draw_independent_patterns(
    n_reads: int, methylated_by_cpg: tuple[int, int, int, int],
    batch: int, rng: np.random.Generator,
) -> np.ndarray:
    """Return `batch` simulated pattern matrices with exact CpG marginals."""
    patterns = np.zeros((batch, n_reads), dtype=np.uint8)
    for index, count in enumerate(methylated_by_cpg):
        bit = 1 << (3 - index)
        if count == n_reads:
            patterns |= bit
        elif count:
            column = np.r_[np.ones(count, dtype=np.uint8),
                           np.zeros(n_reads - count, dtype=np.uint8)]
            shuffled = rng.permuted(np.tile(column, (batch, 1)), axis=1)
            patterns |= shuffled * bit
    return patterns


def process_rlm_file(
    input_file: Union[str, Path], output_file: Union[str, Path],
    n_sims: int = 1000, seed: int = 42, chunk_size: int = 10000,
) -> None:
    """Compare RLM entropy with independent CpGs at observed coverage/marginals."""
    if n_sims < 1 or chunk_size < 1 or seed < 0:
        raise ValueError("n_sims and chunk_size must be positive; seed must be nonnegative")
    source, destination = Path(input_file), Path(output_file)
    if source.resolve() == destination.resolve():
        raise ValueError("input and output must be different files")
    with source.open(encoding="utf-8") as handle:
        header = tuple(handle.readline().rstrip("\r\n").split("\t"))
    if header != REQUIRED_COLUMNS:
        raise ValueError("input must be the unmodified 23-column RLM entropy output")

    @lru_cache(maxsize=256)
    def distribution(
        signature: tuple[int, int, int, int, int]
    ) -> tuple[np.ndarray, float, float, float, float]:
        n_reads, *marginals = signature
        rng = np.random.default_rng(np.random.SeedSequence([seed, *signature]))
        draws = simulate_independent_null(n_reads, tuple(marginals), n_sims, rng)
        return (draws, float(draws.mean()), float(draws.std()),
                float(np.quantile(draws, 0.05)), float(np.quantile(draws, 0.95)))

    with pd.read_csv(source, sep="\t", chunksize=chunk_size) as reader, \
            destination.open("w", encoding="utf-8") as handle:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(handle, sep="\t", index=False)
        for chunk in reader:
            depths, methylated, observed = _validate_chunk(chunk)
            pattern_counts = chunk.loc[:, PATTERNS].to_numpy(dtype=np.int64)
            methylated_by_cpg = pattern_counts @ PATTERN_BITS
            results = []
            for index, (n_reads, totals, entropy) in enumerate(
                    zip(depths, methylated_by_cpg, observed)):
                n_reads = int(n_reads)
                marginal_counts = tuple(int(count) for count in totals)
                null, null_mean, null_sd, null_p05, null_p95 = distribution(
                    (n_reads, *marginal_counts))
                marginals = np.asarray(marginal_counts) / n_reads
                maximum_entropy = float(binary_entropy(marginals).sum() / 4)
                total_correlation = maximum_entropy - float(entropy)
                if total_correlation < -1e-10:
                    raise ValueError("observed entropy exceeds the CpG-marginal bound")
                row = {
                    "chrom": chunk["#chr"].iat[index],
                    "start": int(chunk["start"].iat[index]),
                    "end": int(chunk["end"].iat[index]),
                    "coverage": n_reads,
                    "methylation": int(methylated[index]) / (4 * n_reads),
                    "observed_entropy": float(entropy),
                    "informative_cpg_count": int(((marginals > 0) &
                                                   (marginals < 1)).sum()),
                    "independence_max_entropy": maximum_entropy,
                    "observed_total_correlation": max(0.0, total_correlation),
                    "independent_null_mean": null_mean,
                    "independent_null_sd": null_sd,
                    "independent_null_p05": null_p05,
                    "independent_null_p95": null_p95,
                    "entropy_minus_independent_null": float(entropy - null_mean),
                    "empirical_p_coordination": float(
                        (1 + np.searchsorted(null, entropy + 1e-12, side="right"))
                        / (n_sims + 1)),
                    "empirical_p_dispersion": float(
                        (1 + n_sims - np.searchsorted(null, entropy - 1e-12,
                                                       side="left")) / (n_sims + 1)),
                    "n_sims": n_sims,
                }
                row.update({f"cpg{j}_methylated_count": count
                            for j, count in enumerate(marginal_counts, 1)})
                row.update({f"cpg{j}_methylation": float(rate)
                            for j, rate in enumerate(marginals, 1)})
                results.append(row)
            pd.DataFrame(results, columns=OUTPUT_COLUMNS).to_csv(
                handle, sep="\t", index=False, header=False)
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare RLM entropy with a coverage- and CpG-marginal-matched null")
    parser.add_argument("--input", required=True, help="Unmodified RLM output_entropy.bed")
    parser.add_argument("--output", required=True, help="Tab-separated results")
    parser.add_argument("--nsims", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--chunk-size", type=int, default=10000)
    args = parser.parse_args()
    process_rlm_file(args.input, args.output, args.nsims, args.seed,
                     args.chunk_size)


if __name__ == "__main__":
    main()
