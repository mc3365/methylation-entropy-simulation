"""Small synthetic calibration check for the CpG-matched coordination test.

This is a diagnostic, not a substitute for library-specific technical controls
or biological replicate validation. Run with --replicates 1000 for a more stable
estimate than the quick default.
"""

import argparse

import numpy as np

from methylation_entropy.rlm_entropy_v3 import entropy_from_counts
from methylation_entropy.rlm_entropy_v4 import simulate_independent_null


WEIGHTS = np.array([8, 4, 2, 1], dtype=np.uint8)
MARGINAL_RATES = np.array([0.2, 0.4, 0.6, 0.8])


def evaluate(n_reads, n_sims, n_replicates, seed, coordinated_fraction):
    rng = np.random.default_rng(seed)
    if not 0 <= coordinated_fraction <= min(MARGINAL_RATES):
        raise ValueError("coordinated fraction must be between 0 and 0.2")
    baseline_rates = (MARGINAL_RATES - coordinated_fraction) / (1 - coordinated_fraction)
    p_values = []
    for _ in range(n_replicates):
        reads = rng.random((n_reads, 4)) < baseline_rates
        if coordinated_fraction:
            coordinated = rng.random(n_reads) < coordinated_fraction
            reads[coordinated, :] = True
        patterns = reads.astype(np.uint8) @ WEIGHTS
        observed = entropy_from_counts(np.bincount(patterns, minlength=16))
        marginal_counts = tuple(int(x) for x in reads.sum(axis=0))
        draws = simulate_independent_null(n_reads, marginal_counts, n_sims, rng)
        p_values.append((1 + np.searchsorted(draws, observed + 1e-12,
                                              side="right")) / (n_sims + 1))
    return float(np.mean(np.asarray(p_values) <= 0.05))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=int, default=40)
    parser.add_argument("--nsims", type=int, default=999)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=71)
    args = parser.parse_args()
    if min(args.reads, args.nsims, args.replicates) < 1 or args.seed < 0:
        parser.error("counts must be positive and seed must be nonnegative")
    for label, fraction in (("independent", 0.0), ("coordinated mixture", 0.15)):
        rate = evaluate(args.reads, args.nsims, args.replicates, args.seed,
                        fraction)
        print(f"{label}: {rate:.3f} fraction with coordination p <= 0.05 "
              f"({args.replicates} datasets; {args.reads} reads each)")


if __name__ == "__main__":
    main()
