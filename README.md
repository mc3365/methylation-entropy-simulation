# Methylation Entropy Simulation

An installable Python module and command-line workflow for comparing observed
methylation entropy with stochastic and minimum-entropy reference models.

## Install

Python 3.9 or newer is required.

```bash
python -m pip install -e .
```

For tests:

```bash
python -m pip install -e '.[dev]'
pytest
```

## Run after RLM (recommended)

Run RLM with `--score entropy`, then pass its **unmodified**
`output_entropy.bed` to this command. RLM's file has a `#chr` header, all 16
four-CpG pattern counts, mean methylation, and `coverage`; these are needed to
match the null to the number of molecules at each locus.

```bash
methylation-entropy-rlm \
  --input output_entropy.bed \
  --output entropy_comparison.tsv \
  --nsims 1000 \
  --seed 42 \
  --no-reference
```

The `matched_null_*` columns simulate the **observed coverage and exact total
methylated-CpG count** for each locus. `entropy_minus_matched_null` and the
empirical lower/upper tail probabilities locate the observed value in that
distribution. The separate `reference_*` columns retain the original 128-read
comparison as a descriptive high-depth reference; they should not be used to
classify a lower-depth locus. The command above uses `--no-reference` so only
the depth-matched simulation runs. Omit that flag to include the optional
128-read reference.
The `minimum_possible_entropy` column is a population mathematical bound at the
observed mean. It is not an allele-specific model, and may be unattainable at a
particular finite coverage.

Try the small, synthetic RLM-format example:

```bash
methylation-entropy-rlm \
  --input examples/rlm_entropy_demo.bed \
  --output entropy_demo.tsv \
  --nsims 1000 --seed 42
```

All five rows have 50% mean methylation. The first three have 16 reads and
represent one repeated pattern, a two-pattern mixture, and all 16 patterns.
The final two have 10 reads; the last has entropy near the matched null but
well below the 128-read reference. See [Methodology and usage](docs/methodology.md).

To time the workflow on synthetic RLM-format data:

```bash
python examples/benchmark_rlm_entropy.py
```

The script reports both row counts and distinct `(coverage, methylated-site
count)` pairs. Runtime depends mainly on the number of distinct pairs because
each pair gets its own simulated distribution. Its default is 1,000 simulations
per pair, with the optional 128-read reference enabled.

## Legacy input

The original `methylation-entropy` command accepts a headerless, tab-separated
BED-like file:

```text
chrom  start  end  entropy  methylation
```

Methylation must be between 0 and 1. Invalid or missing values are skipped.
Because this format has no read depth or pattern counts, its fixed-depth null is
descriptive only. The included `mef_entropy_p0.bed` uses this legacy format; it
cannot be upgraded to a depth-matched analysis without the original RLM output.

## Legacy run

```bash
methylation-entropy \
  --input mef_entropy_p0.bed \
  --output entropy_simulation.tsv \
  --lookup minimum_entropy_lookup.json \
  --nsims 1000 \
  --nreads 128 \
  --workers 4 \
  --seed 42
```

Rebuild the minimum-entropy lookup:

```bash
methylation-entropy-lookup --output-dir results
```

`--seed` makes a run reproducible. The original analysis used 128 simulated
reads and 4 CpGs per read. The bundled lookup was built with 16 reads and 4
CpGs per read.

## Files

- `src/methylation_entropy/entropy_simulation_v2.py`: packaged simulation implementation
- `src/methylation_entropy/rlm_entropy_v3.py`: direct RLM integration with a depth-matched null
- `src/methylation_entropy/create_lookup_v2.py`: packaged lookup generator
- `examples/rlm_entropy_demo.bed`: small synthetic RLM-format demonstration
- `examples/benchmark_rlm_entropy.py`: reproducible synthetic runtime check
- `tests/`: automated tests
- `docs/methodology.md`: scientific rationale, complete workflow, and interpretation guide
- `example_mef_metaplots.Rmd`: downstream visualization example (site-specific paths)
- `simulation_hpc_raw_mef.sh`: original SLURM example (site-specific paths)

See [Methodology and usage](docs/methodology.md) for the full explanation of the
simulation design and guidance for interpreting its output.

The version suffixes retain earlier implementations for provenance.

## Citation

Add the associated manuscript or project citation before publishing.
