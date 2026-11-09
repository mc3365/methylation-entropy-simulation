# Methylation Level and Entropy Simulation

## Purpose

This pipeline analyzes four-CpG methylation patterns from whole-genome bisulfite
sequencing (WGBS). Its recommended question is whether CpGs on the same
molecules are more coordinated than expected from their individual methylation
levels and the number of observed molecules. It does not directly identify cell
clones, alleles, or regulatory mechanisms.

The original workflow used a pooled methylation shuffle at a fixed 128-read
depth and a separate mathematical minimum. The v3 workflow matched read depth
but still pooled the four CpGs. The recommended v4 workflow matches depth and
**each CpG's methylated-call count**. Earlier methods remain available for
comparison, but their outputs answer different questions.

## Recommended RLM workflow (v4)

RLM's complete 23-column entropy output contains the 16 four-CpG pattern counts
and coverage. For each window, v4 derives four CpG-specific methylated-call
counts from those pattern counts. Each null draw independently permutes each
CpG column across the observed number of molecules. Thus every draw preserves
the exact coverage and four observed CpG marginals while destroying
within-molecule association. No rounding of a decimal mean is needed.

```bash
methylation-entropy-rlm \
  --input sample_entropy.bed \
  --output sample_cpg_null.tsv \
  --nsims 1000 --seed 42
```

RLM's entropy is Shannon entropy of the 16 pattern frequencies divided by four.
The v4 output reports the same observed value, four marginal counts and rates,
and the conditional independence null's mean and 5th/95th percentiles. A
negative `entropy_minus_independent_null` means more coordinated patterns than
the null mean. `independent_null_percentile` is the observed entropy's midrank
within the simulated null, providing a direct location for plots: values near
zero indicate coordination, values near one indicate dispersion, and 0.5 is
central or entirely tied. `empirical_p_coordination` counts simulations with entropy at or
below the observed entropy; `empirical_p_dispersion` counts those at or above
it. Both use `(1 + count)/(1 + nsims)` and include ties. At 1,000 draws the
smallest value is `1/1001`; deeper tails need more simulations. The percentile
range describes simulated samples under the null, **not** uncertainty in a
biological effect estimate.

| Output column(s) | Meaning |
| --- | --- |
| `chrom`, `start`, `end`, `coverage` | RLM locus and number of complete four-CpG molecules |
| `methylation`, `observed_entropy` | RLM summaries recomputed and checked from pattern counts |
| `informative_cpg_count` | Number of CpGs that vary across the observed molecules; zero or one cannot identify inter-CpG dependence |
| `cpg1_methylated_count` through `cpg4_methylated_count`; `cpg1_methylation` through `cpg4_methylation` | Per-position counts and fractions; CpG 1 is the leftmost pattern character |
| `independence_max_entropy` | Population entropy ceiling at the observed CpG-specific rates |
| `observed_total_correlation` | Ceiling minus observed entropy; a descriptive, finite-sample-biased dependence measure |
| `independent_null_mean`, `independent_null_sd`, `independent_null_p05`, `independent_null_p95` | Conditional independence entropy distribution at observed depth and marginals |
| `entropy_minus_independent_null` | Observed minus null mean; negative suggests excess coordination |
| `independent_null_percentile` | Midrank location of observed entropy in its matched null, from 0 (low) to 1 (high) |
| `empirical_p_coordination`, `empirical_p_dispersion` | Inclusive one-sided Monte Carlo tail probabilities |
| `n_sims` | Number of null draws used for this row |

The population entropy ceiling at the observed marginal rates is
`independence_max_entropy = sum_j h(p_j)/4`, where `h` is binary entropy.
`observed_total_correlation` is that ceiling minus observed entropy. Its raw
value can be positive under independence at finite depth because empirical
joint entropy is downward-biased. Compare it with the matched simulation rather
than treating zero as a finite-sample cutoff. The ceiling need not be attainable
by a finite set of reads.

At 50% pooled methylation, 16 identical `ggGG` reads have entropy zero but each
CpG is constant across reads; v4 correctly gives a degenerate null at zero.
Its coordination p-value is one because there is no variation to test, not
because the data demonstrate independence.
Eight `gggg` plus eight `GGGG` reads also have 50% pooled methylation, but all
four CpGs vary in lockstep, yielding a small coordination-tail probability.
This distinction cannot be made from the five-column legacy output.

The null assumes molecules are exchangeable after conditioning on the four
observed CpG counts. Bisulfite conversion errors, PCR duplicates, mapping bias,
cell-type mixtures, allele-specific methylation, and heterogeneous molecule
sampling can affect interpretation. Overlapping four-CpG windows share CpGs and
often molecules, so window p-values are not independent biological replicates.
Low entropy can arise from positive concordance, anticorrelation, or other
dependence among CpGs; this omnibus entropy test does not identify a specific
pattern such as fully methylated molecules.

For DMV questions, summarize window-level effects within each DMV and compare
those summaries across biological replicates. Keep simulation percentiles
separate from between-DMV and between-replicate variation in figures.

The tool has been checked against synthetic RLM-format data. Full real RLM
outputs are needed to validate input compatibility, coverage distributions,
simulation calibration, and biological conclusions before publication.

The reproducible `examples/calibrate_cpg_null.py` diagnostic generates
independent reads and a coordinated mixture with the same population CpG
methylation rates. With its quick defaults (200 datasets, 40 reads per dataset,
999 null draws), 5.0% of independent datasets and 75.0% of coordinated datasets
had `empirical_p_coordination <= 0.05` at seed 71. These rates depend on the
specified depth, rates, mixture, and seed; they are a code calibration check,
not a measured biological effect or complete operating-characteristic study.

## Previous pooled RLM workflow (v3)

The `methylation-entropy-rlm-pooled` command accepts RLM's complete
`output_entropy.bed` without dropping columns. RLM reports entropy, the counts
of all 16 four-CpG patterns, pooled mean methylation, and coverage for each
four-CpG locus. Its entropy is Shannon entropy of the pattern frequencies
divided by four, the same convention used here. The command verifies that the
counts sum to coverage and reproduce RLM's entropy and mean (within the
rounding of its text output). It reads the file in chunks and caches simulations
for repeated `(coverage, methylated-site count)` pairs.

For a locus with `n` complete reads and `k` methylated CpG calls, the matched
null randomly permutes exactly `k` ones among `4n` positions, reshapes them
into `n` four-CpG patterns, and calculates entropy. Its mean and 5th/95th
percentiles describe the distribution of a **new sample at that depth under
this null**, not uncertainty in the observed entropy estimate. The input's
16 pattern counts make `n` and `k` exact, without rounding a decimal mean.

The output contains:

| Column(s) | Interpretation |
| --- | --- |
| `coverage`, `methylation`, `observed_entropy` | RLM summaries checked against pattern counts |
| `matched_null_mean`, `matched_null_p05`, `matched_null_p95` | Null distribution at observed depth and methylated-call count |
| `entropy_minus_matched_null` | Observed minus null mean; negative means more concentrated patterns than this null |
| `empirical_p_lower`, `empirical_p_upper` | Inclusive one-sided tails, `(1 + number of equally or more extreme simulations)/(1 + nsims)` |
| `minimum_possible_entropy` | Population lower bound at pooled mean; it may be unattainable at finite depth |
| `reference_nreads`, `reference_mean`, `reference_p05`, `reference_p95` | Optional fixed-depth descriptive reference, 128 reads by default |

The minimum bound is obtained by considering mixtures of at most two four-CpG
patterns with the required pooled mean. It is not a model of alleles. At 50%
methylation, for example, every molecule can have the same `ggGG` pattern,
giving entropy zero. A 50/50 mixture of `gggg` and `GGGG` instead has entropy
0.25. The synthetic `examples/rlm_entropy_demo.bed` includes both and a third
50%-methylated case with all 16 patterns once (entropy 1). Its final row has
10 reads and entropy about 0.73: ordinary under a 10-read random null at 50%
methylation, yet far below the descriptive 128-read reference.

Increasing `--nsims` makes the estimated null percentiles and tail
probabilities more precise; it does not change the assumed read depth. With
1,000 simulations, the smallest reported empirical tail probability is
`1/1001`. The optional 128-read reference is useful for showing the high-depth
stochastic expectation, but the matched columns should be used when asking
whether a measured locus is unusual.

A small tail probability indicates departure from **this pooled-mean shuffle
null**. It does not by itself establish biological regulation: different
CpGs can have different methylation rates, and cell or haplotype mixtures can
also produce structure. The recommended v4 null above addresses the unequal
CpG-rate limitation. Overlapping four-CpG windows are still not independent
tests.

## Legacy five-column pipeline overview (v2)

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

## Historical metaplot example (v2; not suitable for v4)

The included `example_mef_metaplots.Rmd` illustrates how to compare entropy
profiles around genomic features such as DNA methylation valleys (DMVs) and CpG
islands. Its file paths are specific to the original HPC environment and must be
changed before use. It consumes v2 columns and must **not** be used to plot v4
results. Its `stochastic_mean` trace is built from the old `stochastic_min`
column, so historical figures made with this script need to be audited before
scientific interpretation.

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

The example's observed ribbon summarizes the 25th to 75th percentiles across
genomic features. Its stochastic ribbon combines across-feature quartiles of
the *per-locus* simulated 5th and 95th percentiles. Those are different kinds
of variation and their overlap is not a calibrated test.

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
- Monte Carlo percentiles describe a null distribution at each locus, not
  uncertainty in a biological effect.
- Metaplot percentiles describe variation across genomic features.
- The stochastic and minimum models are reference baselines, not direct evidence
  for a particular molecular mechanism.
