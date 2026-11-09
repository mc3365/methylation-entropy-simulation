"""Checks for the coverage- and CpG-marginal-matched RLM workflow."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from methylation_entropy.rlm_entropy_v3 import (
    PATTERNS, PATTERN_WEIGHTS, REQUIRED_COLUMNS, entropy_from_counts,
    process_rlm_file as process_pooled_file,
)
from methylation_entropy.rlm_entropy_v4 import (
    PATTERN_BITS, _draw_independent_patterns, binary_entropy,
    process_rlm_file, simulate_independent_null,
)


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "rlm_entropy_demo.bed"


class CpgMatchedNullTests(unittest.TestCase):
    def test_every_draw_preserves_all_four_cpg_counts(self):
        counts = (2, 5, 7, 9)
        patterns = _draw_independent_patterns(
            12, counts, 50, np.random.default_rng(11))
        self.assertEqual(patterns.shape, (50, 12))
        for index, count in enumerate(counts):
            observed = ((patterns >> (3 - index)) & 1).sum(axis=1)
            np.testing.assert_array_equal(observed, np.full(50, count))

    def test_entropy_and_pattern_conventions(self):
        np.testing.assert_array_equal(PATTERN_BITS[0], [0, 0, 0, 0])
        np.testing.assert_array_equal(PATTERN_BITS[-1], [1, 1, 1, 1])
        np.testing.assert_array_equal(binary_entropy([0, 0.5, 1]), [0, 1, 0])
        null = simulate_independent_null(
            8, (0, 0, 8, 8), 50, np.random.default_rng(3))
        np.testing.assert_array_equal(null, np.zeros(50))

    def test_simulation_matches_exact_small_null_probabilities(self):
        # With three reads and one methylated call at each of two CpGs, the
        # calls land on the same molecule with probability 1/3 and on distinct
        # molecules with probability 2/3. These cases have different entropies.
        null = simulate_independent_null(
            3, (1, 1, 0, 0), 30000, np.random.default_rng(23))
        low_entropy = entropy_from_counts(np.array([
            2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]))
        observed_low_fraction = float(np.mean(np.isclose(null, low_entropy)))
        self.assertAlmostEqual(observed_low_fraction, 1 / 3, delta=0.015)

    def test_same_pooled_mean_has_distinct_coordination_interpretations(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result.tsv"
            process_rlm_file(EXAMPLE, output, n_sims=2000, seed=9)
            result = pd.read_csv(output, sep="\t")
        self.assertEqual(result["coverage"].tolist(), [16, 16, 16, 10, 10])
        self.assertEqual(result["methylation"].tolist(), [0.5] * 5)
        self.assertEqual(result.loc[0, "independence_max_entropy"], 0)
        self.assertEqual(result.loc[0, "independent_null_mean"], 0)
        self.assertEqual(result.loc[0, "independent_null_sd"], 0)
        self.assertEqual(result.loc[0, "informative_cpg_count"], 0)
        self.assertEqual(result.loc[0, "empirical_p_coordination"], 1)
        self.assertEqual(result.loc[0, "independent_null_percentile"], 0.5)
        self.assertEqual(result.loc[1, "independence_max_entropy"], 1)
        self.assertEqual(result.loc[1, "informative_cpg_count"], 4)
        self.assertLess(result.loc[1, "empirical_p_coordination"], 0.01)
        self.assertLess(result.loc[1, "independent_null_percentile"], 0.01)
        self.assertLess(result.loc[2, "empirical_p_dispersion"], 0.01)
        self.assertGreater(result.loc[2, "independent_null_percentile"], 0.99)
        self.assertGreater(result.loc[4, "empirical_p_coordination"], 0.05)
        self.assertTrue((result["observed_total_correlation"] >= 0).all())
        self.assertFalse(any("reference" in column or "allelic" in column
                             for column in result.columns))

    def test_per_cpg_null_differs_from_pooled_null_when_cpg_rates_differ(self):
        counts = np.zeros(16, dtype=int)
        for pattern in ("gggg", "gggG", "ggGG", "gGGG"):
            counts[PATTERNS.index(pattern)] = 4
        n_reads = int(counts.sum())
        row = {"#chr": "chr1", "start": 100, "end": 102,
               "entropy": entropy_from_counts(counts),
               "epipolymorphism": 1 - float(np.sum((counts / n_reads) ** 2)),
               "mean_methylation": int(counts @ PATTERN_WEIGHTS) / (4 * n_reads),
               "coverage": n_reads}
        row.update(zip(PATTERNS, counts))
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "in.bed"
            matched = Path(directory) / "matched.tsv"
            pooled = Path(directory) / "pooled.tsv"
            pd.DataFrame([row], columns=REQUIRED_COLUMNS).to_csv(
                source, sep="\t", index=False)
            process_rlm_file(source, matched, n_sims=1000, seed=17)
            process_pooled_file(source, pooled, n_sims=1000, seed=17,
                                reference_nreads=None)
            matched_result = pd.read_csv(matched, sep="\t").iloc[0]
            pooled_result = pd.read_csv(pooled, sep="\t").iloc[0]
        self.assertEqual(
            [matched_result[f"cpg{j}_methylated_count"] for j in range(1, 5)],
            [0, 4, 8, 12])
        self.assertLess(matched_result["independent_null_mean"],
                        pooled_result["matched_null_mean"] - 0.1)

    def test_reproducible_across_chunk_sizes(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "a.tsv", Path(directory) / "b.tsv"
            process_rlm_file(EXAMPLE, first, n_sims=100, seed=7, chunk_size=1)
            process_rlm_file(EXAMPLE, second, n_sims=100, seed=7, chunk_size=4)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_rejects_inconsistent_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "bad.bed", Path(directory) / "out.tsv"
            data = pd.read_csv(EXAMPLE, sep="\t")
            data.loc[0, "coverage"] = 15
            data.to_csv(source, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "coverage"):
                process_rlm_file(source, output, n_sims=10)


if __name__ == "__main__":
    unittest.main()
