"""Checks for the depth-matched RLM workflow."""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from methylation_entropy.rlm_entropy_v3 import (
    PATTERNS, minimum_possible_entropy, process_rlm_file,
)


EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "rlm_entropy_demo.bed"


class RlmEntropyTests(unittest.TestCase):
    def test_minimum_is_a_population_bound(self):
        self.assertEqual(minimum_possible_entropy(0.5), 0)
        self.assertEqual(minimum_possible_entropy(0), 0)
        self.assertGreater(minimum_possible_entropy(0.1), 0)
        with self.assertRaises(ValueError):
            minimum_possible_entropy(float("nan"))

    def test_rlm_example_matches_each_depth_and_keeps_128_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "out.tsv"
            process_rlm_file(EXAMPLE, output, n_sims=150, seed=9, chunk_size=2)
            result = pd.read_csv(output, sep="\t")
        self.assertEqual(result["coverage"].tolist(), [16, 16, 16, 10, 10])
        self.assertEqual(result["observed_entropy"].tolist()[:4], [0, 0.25, 1, 0.25])
        self.assertAlmostEqual(result.loc[4, "observed_entropy"], 0.7304820237)
        self.assertEqual(result["minimum_possible_entropy"].tolist(), [0] * 5)
        self.assertEqual(result["reference_nreads"].tolist(), [128] * 5)
        self.assertEqual(result["reference_mean"].nunique(), 1)
        self.assertLess(result.loc[3, "matched_null_mean"],
                        result.loc[0, "matched_null_mean"])
        self.assertAlmostEqual(
            result.loc[0, "entropy_minus_matched_null"],
            result.loc[0, "observed_entropy"] - result.loc[0, "matched_null_mean"])
        self.assertLess(result.loc[0, "empirical_p_lower"], 0.05)
        self.assertLess(result.loc[2, "empirical_p_upper"], 0.05)
        self.assertGreater(result.loc[4, "observed_entropy"], result.loc[4, "matched_null_p05"])
        self.assertLess(result.loc[4, "observed_entropy"], result.loc[4, "matched_null_p95"])
        self.assertLess(result.loc[4, "observed_entropy"], result.loc[4, "reference_p05"])

    def test_seed_does_not_depend_on_chunk_size(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "a.tsv", Path(directory) / "b.tsv"
            process_rlm_file(EXAMPLE, first, n_sims=40, seed=7, chunk_size=1)
            process_rlm_file(EXAMPLE, second, n_sims=40, seed=7, chunk_size=4)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_rejects_coverage_that_disagrees_with_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            bad, output = Path(directory) / "bad.bed", Path(directory) / "out.tsv"
            data = pd.read_csv(EXAMPLE, sep="\t")
            data.loc[0, "coverage"] = 15
            data.to_csv(bad, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "coverage"):
                process_rlm_file(bad, output, n_sims=10, reference_nreads=None)

    def test_example_has_exactly_sixteen_pattern_columns(self):
        data = pd.read_csv(EXAMPLE, sep="\t")
        self.assertEqual(len(PATTERNS), 16)
        self.assertTrue(all(pattern in data.columns for pattern in PATTERNS))


if __name__ == "__main__":
    unittest.main()
