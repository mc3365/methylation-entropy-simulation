"""Depth-matched entropy references for the full RLM entropy output."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd


PATTERNS = tuple(format(index, "04b").replace("0", "g").replace("1", "G")
                 for index in range(16))
PATTERN_WEIGHTS = np.array([pattern.count("G") for pattern in PATTERNS])
REQUIRED_COLUMNS = ("#chr", "start", "end", "entropy", "epipolymorphism",
                    *PATTERNS, "mean_methylation", "coverage")
OUTPUT_COLUMNS = ("chrom", "start", "end", "observed_entropy", "methylation",
                  "coverage", "minimum_possible_entropy", "matched_null_mean",
                  "matched_null_p05", "matched_null_p95", "entropy_minus_matched_null",
                  "empirical_p_lower", "empirical_p_upper", "reference_nreads",
                  "reference_mean", "reference_p05", "reference_p95")


def minimum_possible_entropy(methylation: float) -> float:
    """Population lower bound for four-CpG entropy at a pooled mean.

    Entropy is concave in pattern frequencies. Its minimum subject to a fixed
    pooled mean therefore occurs at a vertex with at most two patterns. This
    bound is not always attainable at a given finite read depth.
    """
    if not math.isfinite(methylation) or not 0 <= methylation <= 1:
        raise ValueError("methylation must be finite and between 0 and 1")
    target = 4 * methylation
    if any(math.isclose(target, weight, abs_tol=1e-12) for weight in range(5)):
        return 0.0
    candidates = []
    for low in range(5):
        for high in range(low + 1, 5):
            if low <= target <= high:
                fraction = (target - low) / (high - low)
                entropy = -(fraction * math.log2(fraction)
                            + (1 - fraction) * math.log2(1 - fraction)) / 4
                candidates.append(entropy)
    return min(candidates)


def entropy_from_counts(counts: np.ndarray) -> float:
    """Use the same per-CpG Shannon convention as RLM."""
    n_reads = int(counts.sum())
    if n_reads == 0:
        raise ValueError("pattern counts must include at least one read")
    positive = counts[counts > 0] / n_reads
    if len(positive) == 1:
        return 0.0
    return float(-np.sum(positive * np.log2(positive)) / 4)


def simulate_null(n_reads: int, methylated_sites: int, n_sims: int,
                  rng: np.random.Generator) -> np.ndarray:
    """Shuffle a fixed methylated-site count across four-CpG reads."""
    if n_reads < 1 or n_sims < 1 or not 0 <= methylated_sites <= 4 * n_reads:
        raise ValueError("invalid read depth, methylated-site count, or simulation count")
    sites = np.r_[np.ones(methylated_sites, dtype=np.uint8),
                  np.zeros(4 * n_reads - methylated_sites, dtype=np.uint8)]
    entropies = np.empty(n_sims, dtype=float)
    weights = np.array([8, 4, 2, 1], dtype=np.uint8)
    for index in range(n_sims):
        patterns = rng.permutation(sites).reshape(n_reads, 4) @ weights
        entropies[index] = entropy_from_counts(np.bincount(patterns, minlength=16))
    return np.sort(entropies)


def _validate_chunk(chunk: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Check that RLM's summaries agree with its pattern counts."""
    for name in ("start", "end", "entropy", "mean_methylation", "coverage", *PATTERNS):
        chunk[name] = pd.to_numeric(chunk[name], errors="raise")
    values = chunk.loc[:, PATTERNS].to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any()
            or not np.equal(values, np.floor(values)).all()):
        raise ValueError("RLM pattern counts must be nonnegative integers")
    counts = values.astype(np.int64)
    depths = counts.sum(axis=1)
    coverage = chunk["coverage"].to_numpy(dtype=float)
    starts = chunk["start"].to_numpy(dtype=float)
    ends = chunk["end"].to_numpy(dtype=float)
    observed = chunk["entropy"].to_numpy(dtype=float)
    means = chunk["mean_methylation"].to_numpy(dtype=float)
    if (not np.isfinite(coverage).all() or not np.equal(coverage, depths).all()
            or (depths < 1).any()):
        raise ValueError("RLM coverage must equal the sum of the 16 pattern counts")
    if (not np.isfinite(starts).all() or not np.isfinite(ends).all()
            or not np.equal(starts, np.floor(starts)).all()
            or not np.equal(ends, np.floor(ends)).all()
            or (starts < 0).any() or (ends <= starts).any()):
        raise ValueError("RLM coordinates must be nonnegative, ordered integers")
    methylated = counts @ PATTERN_WEIGHTS
    exact_means = methylated / (4 * depths)
    exact_entropy = np.array([entropy_from_counts(row) for row in counts])
    if (not np.isfinite(means).all() or not np.isfinite(observed).all()
            or not np.allclose(means, exact_means, rtol=0, atol=1e-5)
            or not np.allclose(observed, exact_entropy, rtol=0, atol=1e-5)):
        raise ValueError("RLM entropy or mean methylation disagrees with pattern counts")
    return depths, methylated, exact_entropy


def process_rlm_file(input_file: Union[str, Path], output_file: Union[str, Path],
                     n_sims: int = 1000, seed: int = 42,
                     reference_nreads: Optional[int] = 128,
                     chunk_size: int = 10000) -> None:
    """Write matched-null results directly from RLM's entropy BED output."""
    if n_sims < 1 or chunk_size < 1 or seed < 0:
        raise ValueError("n_sims and chunk_size must be positive; seed must be nonnegative")
    if reference_nreads is not None and reference_nreads < 1:
        raise ValueError("reference_nreads must be positive")
    source, destination = Path(input_file), Path(output_file)
    if source.resolve() == destination.resolve():
        raise ValueError("input and output must be different files")
    with source.open(encoding="utf-8") as handle:
        header = tuple(handle.readline().rstrip("\r\n").split("\t"))
    if header != REQUIRED_COLUMNS:
        raise ValueError("input must be the unmodified 23-column RLM entropy output")
    cache: dict[tuple[int, int], tuple[np.ndarray, float, float, float]] = {}

    def distribution(n_reads: int, methylated: int) -> tuple[np.ndarray, float, float, float]:
        key = (n_reads, methylated)
        if key not in cache:
            rng = np.random.default_rng(np.random.SeedSequence([seed, n_reads, methylated]))
            draws = simulate_null(n_reads, methylated, n_sims, rng)
            cache[key] = (draws, float(np.mean(draws)),
                          float(np.quantile(draws, 0.05)),
                          float(np.quantile(draws, 0.95)))
        return cache[key]

    with pd.read_csv(source, sep="\t", chunksize=chunk_size) as reader, \
            destination.open("w", encoding="utf-8") as handle:
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(handle, sep="\t", index=False)
        for chunk in reader:
            depths, methylated, observed = _validate_chunk(chunk)
            results = []
            chroms = chunk["#chr"].to_numpy()
            starts = chunk["start"].to_numpy()
            ends = chunk["end"].to_numpy()
            for index, (n_reads, meth_count, entropy) in enumerate(
                    zip(depths, methylated, observed)):
                n_reads, meth_count = int(n_reads), int(meth_count)
                methylation = meth_count / (4 * n_reads)
                null, null_mean, null_p05, null_p95 = distribution(n_reads, meth_count)
                row = {"chrom": chroms[index],
                       "start": int(starts[index]),
                       "end": int(ends[index]),
                       "observed_entropy": entropy,
                       "methylation": methylation,
                       "coverage": n_reads,
                       "minimum_possible_entropy": minimum_possible_entropy(methylation),
                       "matched_null_mean": null_mean,
                       "matched_null_p05": null_p05,
                       "matched_null_p95": null_p95,
                       "entropy_minus_matched_null": float(entropy - null_mean),
                       "empirical_p_lower": float((1 + np.searchsorted(null, entropy + 1e-12,
                                                                       side="right")) / (n_sims + 1)),
                       "empirical_p_upper": float((1 + n_sims - np.searchsorted(
                           null, entropy - 1e-12, side="left")) / (n_sims + 1))}
                if reference_nreads is not None:
                    reference_methylated = round(4 * reference_nreads * methylation)
                    _, reference_mean, reference_p05, reference_p95 = distribution(
                        reference_nreads, reference_methylated)
                    row.update(reference_nreads=reference_nreads,
                               reference_mean=reference_mean,
                               reference_p05=reference_p05,
                               reference_p95=reference_p95)
                results.append(row)
            pd.DataFrame(results, columns=OUTPUT_COLUMNS).to_csv(
                handle, sep="\t", index=False, header=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare RLM four-CpG entropy with a read-depth-matched null")
    parser.add_argument("--input", required=True, help="Unmodified RLM output_entropy.bed")
    parser.add_argument("--output", required=True)
    parser.add_argument("--nsims", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--reference-nreads", type=int, default=128,
                        help="Optional fixed-depth descriptive reference (default: 128)")
    parser.add_argument("--no-reference", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=10000)
    args = parser.parse_args()
    process_rlm_file(args.input, args.output, args.nsims, args.seed,
                     None if args.no_reference else args.reference_nreads,
                     args.chunk_size)


if __name__ == "__main__":
    main()
