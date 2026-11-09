# Changelog

## 0.4.0 — CpG-matched independence null

- `methylation-entropy-rlm` now permutes each of the four CpGs independently
  across the observed molecules, preserving coverage and each CpG's exact
  methylated-call count. It reports CpG marginals, an entropy ceiling, a
  descriptive total-correlation value, null summaries, and one-sided Monte
  Carlo tail probabilities.
- The previous pooled-count RLM method remains available as
  `methylation-entropy-rlm-pooled`. Its output columns and optional 128-read
  reference are unchanged. Users who need to reproduce v3 results should use
  this command explicitly.
- Added synthetic calibration, performance benchmarking, a generic RLM BAM
  pilot wrapper, and a real-data validation plan. This release has not yet been
  validated on an actual full RLM output file.

## 0.3.0 — Depth-matched RLM workflow

- Added direct ingestion of full RLM entropy output and a null that matches
  each window's read depth and total methylated-call count.
- Kept the 128-read stochastic comparison as an optional descriptive reference.
