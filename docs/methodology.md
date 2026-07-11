# Methylation Level and Entropy Simulation

## Purpose

This pipeline analyzes methylation entropy from whole-genome bisulfite
sequencing (WGBS) data. Its goal is to separate entropy changes that arise
naturally from changes in mean methylation from entropy patterns that may
reflect biological regulation of methylation heterogeneity.

Methylation entropy tends to increase as methylation approaches 50%. Therefore,
an entropy value cannot be interpreted independently of its corresponding mean
methylation level. The pipeline compares observed entropy with two reference
models:

1. **Minimum (allelic) entropy:** the theoretical lower bound at a given
   methylation level.
2. **Stochastic entropy:** the range expected when methylated CpGs are randomly
   distributed among DNA molecules while the total methylation level is held
   fixed.

These comparisons help address the central question: are observed entropy
changes explained by mean methylation alone, or is methylation heterogeneity
higher or lower than expected by chance?

## Pipeline overview

The workflow has three computational stages:

1. Generate or load a minimum-entropy lookup table.
2. Add minimum-entropy and stochastic reference values to observed WGBS data.
3. Visualize the observed and reference profiles around genomic features.

## 1. Minimum-entropy lookup table

The `methylation-entropy-lookup` command calculates the minimum possible
methylation entropy for every possible number of methylated CpGs, given a fixed
number of reads and CpGs per read.

For each total methylated-CpG count, the algorithm searches possible allocations
of reads containing 0, 1, 2, 3, or 4 methylated CpGs. It retains the allocation
with the lowest entropy and stores that value by methylation level.

The bundled `minimum_entropy_lookup.json` was calculated using:

- 16 reads
- 4 CpGs per read
- 64 total CpG observations

The 16-read setting reduces the cost of the exhaustive search. The stochastic
simulation uses 128 reads by default, so this difference should be considered
when interpreting the minimum and stochastic references.

To rebuild a lookup table:

```bash
methylation-entropy-lookup \
  --output-dir results \
  --nreads 16 \
  --ncpgs 4
```

This produces:

- `minimum_entropy_lookup.json`: methylation levels mapped to minimum entropy
- `minimum_entropy_lookup.csv`: a table for convenient inspection
- `minimum_entropy_combinations.json`: the minimizing read allocation at each level

## 2. Main entropy simulation

The `methylation-entropy` command reads observed methylation and entropy values
and adds the two reference models for each genomic locus.

For every valid row, it:

1. Finds the closest methylation level in the minimum-entropy lookup.
2. Calculates the corresponding allelic minimum entropy.
3. Preserves the rounded total number of methylated CpG observations.
4. Randomly distributes those methylated sites across simulated reads.
5. Calculates entropy for each simulated distribution.
6. Reports the simulation mean and the 5th and 95th percentiles.

The default Monte Carlo settings are 1,000 simulations, 128 reads, and 4 CpGs
per read. Multiprocessing can distribute genomic loci across CPU cores.

### Input format

The input is a headerless, tab-separated BED-like file with five columns:

```text
chrom  start  end  entropy  methylation
```

For example:

```text
chr1    10000    10004    0.3125    0.500
chr1    10004    10008    0.1250    0.250
```

Methylation values must be between 0 and 1. Rows with missing or out-of-range
methylation values are skipped.

### Example command

```bash
methylation-entropy \
  --input mef_entropy_p0.bed \
  --output entropy_simulation.tsv \
  --lookup minimum_entropy_lookup.json \
  --nsims 1000 \
  --nreads 128 \
  --ncpgs 4 \
  --workers 4 \
  --seed 42
```

Set `--seed` for reproducible random simulations. Set `--workers 1` when
debugging or when deterministic execution order is useful.

### Output columns

The output is a tab-separated table containing:

| Column | Meaning |
| --- | --- |
| `chrom` | Chromosome or sequence name |
| `start` | Genomic start coordinate |
| `end` | Genomic end coordinate |
| `observed_entropy` | Entropy measured in the input data |
| `methylation` | Observed mean methylation level |
| `allelic_entropy` | Closest theoretical minimum from the lookup table |
| `stochastic_mean` | Mean entropy across random simulations |
| `stochastic_min` | 5th percentile of simulated entropy |
| `stochastic_max` | 95th percentile of simulated entropy |

## 3. Visualization around genomic features

The included `example_mef_metaplots.Rmd` illustrates how to compare entropy
profiles around genomic features such as DNA methylation valleys (DMVs) and CpG
islands. Its file paths are specific to the original HPC environment and must be
changed before use.

The example workflow:

1. Loads genomic coordinates for target features.
2. Converts observed and simulated output columns to `GRanges` objects.
3. Uses `EnrichedHeatmap::normalizeToMatrix` to create binned matrices around
   each target region.
4. Calculates the mean and interquartile range across features for each bin.
5. Uses `ggplot2` to plot observed, stochastic, and allelic entropy together
   with methylation level.

The `extend` parameter controls the genomic distance around each target, and
`w` controls bin width. Both should be selected for the biological question and
feature scale.

The example metaplot ribbons summarize the 25th to 75th percentiles across
genomic features. This differs from `stochastic_min` and `stochastic_max`, which
are the 5th and 95th percentiles across simulations at each individual locus.

## Interpretation guide

The example visualization uses the following conceptual references:

- **Observed:** entropy measured from the WGBS data
- **Stochastic:** entropy expected from a random arrangement at the same mean methylation
- **Allelic:** the theoretical minimum entropy
- **Methylation:** the mean methylation level

### Observed entropy within the stochastic range

If observed entropy falls within the stochastic range, the methylation
heterogeneity may be largely explained by random CpG arrangement at the observed
mean methylation level.

### Observed entropy above the stochastic range

Entropy consistently above the stochastic range indicates more heterogeneous
methylation patterns than the random model predicts. This may suggest biological
processes that promote epigenetic variability, but additional validation is
needed before assigning a mechanism.

### Observed entropy near the allelic minimum

Entropy near the theoretical minimum indicates a comparatively uniform
methylation configuration. This may be consistent with tightly constrained or
regulated patterns, although the minimum is a mathematical reference rather than
proof of a specific biological process.

## HPC use

`simulation_hpc_raw_mef.sh` is retained as an example SLURM submission script.
Before using it, update the account-specific Conda environment, input directory,
output directory, email address, and cluster resource requests. The portable CLI
commands above are preferred for new workflows.

## Important considerations

- Entropy must be interpreted together with mean methylation.
- The lookup and simulation read depths differ by default (16 versus 128 reads).
- Monte Carlo percentiles describe simulation uncertainty at each locus.
- Metaplot percentiles describe variation across genomic features.
- The stochastic and minimum models are reference baselines, not direct evidence
  for a particular molecular mechanism.
