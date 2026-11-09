# Methylation Entropy Simulation

An installable Python tool for testing whether four-CpG methylation patterns
are more coordinated than expected when the CpGs vary independently. The
recommended analysis preserves each locus's read depth and each CpG's observed
methylated-call count.

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
four-CpG pattern counts, mean methylation, and `coverage`. The 16 counts supply
both the number of molecules and the methylated-call count at each CpG.

```bash
methylation-entropy-rlm \
  --input output_entropy.bed \
  --output entropy_comparison.tsv \
  --nsims 1000 \
  --seed 42
```

For every simulated draw, the tool independently permutes the four CpG columns
across the observed number of molecules. This preserves all four CpG methylation
counts exactly while removing coordination among CpGs. The
`independent_null_*` columns describe the resulting entropy distribution;
`entropy_minus_independent_null` is negative when the observed patterns are more
coordinated. `independent_null_percentile` directly locates the observed entropy
within its matched null: values near 0 indicate coordination, values near 1
indicate dispersion, and ties use their midrank. `empirical_p_coordination` is
the inclusive lower-tail Monte Carlo
probability. A small value indicates departure from this independence null,
**not proof of a clone or biological regulation**. The 5th–95th percentile
range is a null distribution, not a confidence interval for a biological effect.

`independence_max_entropy` is the entropy of an independent population at the
observed CpG methylation rates. `observed_total_correlation` is the gap between
that value and observed entropy. Its raw value is positive even for many finite
independent samples, so use the simulated null for calibration. Check
`informative_cpg_count` and `independent_null_sd` before interpreting a tail
probability; a degenerate null cannot test coordination.

Try the small, synthetic RLM-format example:

```bash
methylation-entropy-rlm \
  --input examples/rlm_entropy_demo.bed \
  --output entropy_demo.tsv \
  --nsims 1000 --seed 42
```

All five rows have 50% mean methylation. The first (`ggGG` on every read) has
zero entropy because each CpG is constant; that is expected under the per-CpG
null. The second (`gggg`/`GGGG`) also has low entropy but is strongly coordinated
because every CpG varies across molecules in lockstep. See
[Methodology and usage](docs/methodology.md).

To time the workflow on synthetic RLM-format data:

```bash
python examples/benchmark_rlm_entropy.py
```

The script reports row counts and distinct `(coverage, four CpG counts)`
signatures. Runtime depends mainly on the number of distinct signatures because
each gets its own simulated distribution. Its default is 1,000 simulations per
signature. Add `--model pooled` to benchmark the previous pooled shuffle.

For a small synthetic calibration check, run
`python examples/calibrate_cpg_null.py`. Its output is a diagnostic of the test
at stated simulation settings, not evidence about a biological sample.

## Previous RLM pooled shuffle

`methylation-entropy-rlm-pooled` retains the v3 method. It matches read depth
but only preserves the **total** methylated-call count across four CpGs. By
default it also computes a separate descriptive 128-read reference; use
`--no-reference` to omit it. This command is retained for comparison with prior
results, not as the recommended independence test.

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
- `src/methylation_entropy/rlm_entropy_v4.py`: recommended depth- and CpG-matched null
- `src/methylation_entropy/rlm_entropy_v3.py`: previous depth-matched pooled shuffle
- `src/methylation_entropy/create_lookup_v2.py`: packaged lookup generator
- `examples/rlm_entropy_demo.bed`: small synthetic RLM-format demonstration
- `examples/benchmark_rlm_entropy.py`: reproducible synthetic runtime check
- `examples/calibrate_cpg_null.py`: synthetic null calibration check
- `examples/run_rlm_entropy_from_bam.sh`: generic RLM pilot for an aligned BAM
- `tests/`: automated tests
- `docs/methodology.md`: scientific rationale, complete workflow, and interpretation guide
- `example_mef_metaplots.Rmd`: historical v2 plot code; not valid for v4 output
- `simulation_hpc_raw_mef.sh`: original SLURM example (site-specific paths)

See [Methodology and usage](docs/methodology.md) for the full explanation of the
simulation design and guidance for interpreting its output. The
[validation plan](docs/validation_plan.md) lists the real-data and biological
checks needed before a publication-ready release.

See the [changelog](CHANGELOG.md) before upgrading from v3: the
`methylation-entropy-rlm` null model and output schema changed in v0.4.

The version suffixes retain earlier implementations for provenance.

## Citation

Add the associated manuscript or project citation before publishing.
