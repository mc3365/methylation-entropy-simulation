import json

import numpy as np
import pandas as pd
import pytest

from methylation_entropy.entropy_simulation_v2 import (
    calculate_methylation_entropy, get_minimum_entropy, process_file,
    simulate_entropy_ranges,
)


def test_entropy_known_values():
    assert calculate_methylation_entropy({"0000": 10}, 4) == 0
    assert calculate_methylation_entropy({"0000": 1, "1111": 1}, 4) == pytest.approx(0.25)


def test_closest_lookup_value():
    assert get_minimum_entropy(0.49, {"0.000": 0, "0.500": 0.25}) == 0.25


def test_reproducible_simulation():
    first = simulate_entropy_ranges(0.5, 8, 4, 10, np.random.default_rng(7))
    second = simulate_entropy_ranges(0.5, 8, 4, 10, np.random.default_rng(7))
    assert first == second
    assert 0 <= first["stochastic_min"] <= first["stochastic_max"] <= 1


def test_one_row_input(tmp_path):
    source, output, lookup = tmp_path / "in.bed", tmp_path / "out.tsv", tmp_path / "lookup.json"
    source.write_text("chr1\t1\t5\t0.2\t0.5\n", encoding="utf-8")
    lookup.write_text(json.dumps({"0.500": 0.25}), encoding="utf-8")
    process_file(source, output, lookup, n_sims=5, n_reads=8, workers=1, seed=1)
    result = pd.read_csv(output, sep="\t")
    assert len(result) == 1
    assert result.loc[0, "allelic_entropy"] == 0.25
