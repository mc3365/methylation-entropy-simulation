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

## Input

The input is a headerless, tab-separated BED-like file:

```text
chrom  start  end  entropy  methylation
```

Methylation must be between 0 and 1. Invalid or missing values are skipped.

## Run

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
- `src/methylation_entropy/create_lookup_v2.py`: packaged lookup generator
- `tests/`: automated tests
- `docs/methodology.md`: scientific rationale, complete workflow, and interpretation guide
- `example_mef_metaplots.Rmd`: downstream visualization example (site-specific paths)
- `simulation_hpc_raw_mef.sh`: original SLURM example (site-specific paths)

See [Methodology and usage](docs/methodology.md) for the full explanation of the
simulation design and guidance for interpreting its output.

The version suffixes distinguish revised code from the original source. Future
revisions should increment the suffix (`v3`, `v4`, and so on) while retaining
earlier versions for provenance.

## Citation

Add the associated manuscript or project citation before publishing.
