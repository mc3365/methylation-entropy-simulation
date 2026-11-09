# Validation plan for an RLM-based release

The v4 per-CpG null is implemented and checked on synthetic RLM-format data.
The following real-data and biological checks are still required before calling
the tool publication-ready.

## 1. Check source files before rerunning alignment

- On the cluster, inspect the full human cancer RLM entropy files, preferably
  under `.../wgbs/human_cancer/RLM/entropy/`. Confirm that a candidate has the
  exact 23-column header, 16 pattern counts, mean methylation, and coverage.
  Reduced five-column BED files and bedGraphs cannot recover those counts.
- Inventory `hc_mef2_bam.tar.gz` without extracting it first. Identify BAMs,
  BAM indexes, samples, biological replicates, alignment reference, and library
  mode. The HC assay is hybrid capture around DMVs; it is appropriate for a DMV
  pilot but not a genome-wide coverage assessment.
- Do not start from the two raw WGBS archives unless the aligned BAMs prove
  unusable or the full MEF WGBS analysis is needed later.

## 2. Run one HC MEF pilot through RLM

- Extract or stage one aligned BAM and its index on cluster scratch. Confirm
  BAM integrity, sorting/indexing, reference contig names, and the metadata or
  tags required by the installed RLM build.
- Run RLM with the matching mouse reference and correct single- or paired-end
  mode. Existing project-specific cluster scripts can supply the local
  executable, reference path, and `-e` entropy-output option; verify these
  paths and flags against the selected BAM before running. The generic
  `examples/run_rlm_entropy_from_bam.sh` wrapper takes the BAM, reference,
  sequencing mode, and aligner explicitly. Upstream RLM defaults to BSMAP and
  a minimum coverage of 10; do not rely on those defaults unless they match
  this library.
- Keep the **full** RLM entropy BED. Feed it to `methylation-entropy-rlm` and
  retain the original RLM file, command, version, seed, and simulation count.
- Check the RLM pattern counts against coverage, entropy, and mean methylation.
  The CLI rejects rows that fail these checks rather than silently guessing.

## 3. Evaluate statistical behavior on real inputs

- Tabulate coverage, four CpG rates, missingness, and the number of informative
  windows in DMVs for each sample. Stratify summaries by depth and methylation.
  RLM's requirement that one molecule spans all four CpGs can favor CpG-dense
  regions, and the hybrid-capture pilot only covers targeted regions.
- Repeat synthetic null calibration across the depths and CpG-rate profiles
  seen in the pilot. Check false-positive calibration, tail resolution, and
  sensitivity to `--nsims`. The included `examples/calibrate_cpg_null.py` is a
  starting diagnostic, not a complete validation study.
- Check library-specific bisulfite conversion, duplication, and mapping QC.
  Add an explicit error model only when error rates can be estimated for the
  relevant libraries.
- Benchmark the CLI on a full cancer RLM file and the HC pilot. The v4
  simulation is batched and caches up to 256 distinct depth/marginal signatures
  to bound memory use; runtime will depend on how often signatures repeat.

## 4. Answer the DMV biological question

- Keep biological replicates separate. Define the DMV set and comparison rule
  before looking at coordination results. Avoid treating overlapping four-CpG
  windows as independent replicates.
- Summarize coverage-matched coordination effects per DMV and per biological
  replicate, then compare passages using replicate-level uncertainty.
- Analyze human cancer samples as a separate context; tumor purity and cell
  mixture are additional explanations for coordinated patterns.
- Treat a fitted fully methylated mixture component, if added later, as a
  model-dependent pattern proportion. Four-CpG bulk data alone cannot establish
  clone fraction or selection. Longer molecules, allele information, or lineage
  evidence would strengthen that interpretation.

## 5. Release gates

- Publish a small de-identified example from an actual RLM output if permitted,
  plus a reproducible end-to-end command and expected output.
- Check package installation, CLI behavior, and tests on supported Python
  versions; keep the synthetic calibration and real-file benchmark results.
- Add a license, citation, changelog, and method description with precise null
  assumptions and limitations before a public release or manuscript claim.
