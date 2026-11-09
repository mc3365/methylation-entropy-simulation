# Bouchet scratch pilot

These two Slurm examples keep extraction, analysis outputs, and Slurm logs on
scratch. They never delete files. They require an existing RLM build, a matching
reference FASTA, and an indexed BAM member inside a tar archive. Set all paths
for the particular project before submission; do not commit private data paths
or sample inventories to this public repository.

## 1. RLM on one HC BAM

From a checkout of this repository, set the paths in the shell. `BAM_MEMBER`
must exactly match a `.bam` member inside the archive; the script also extracts
its `.bam.bai`. The example is configured for paired-end BSMAP and minimum
coverage 10. Check the selected BAM's header before using those settings.

```bash
export ARCHIVE=/path/to/input_bams.tar.gz
export BAM_MEMBER=mef/example.bam
export REFERENCE=/path/to/matching/genome.fa
export RLM_BIN=/path/to/existing/RLM
export RLM_CONDA_ENV=/path/to/existing/rlm_environment
export OUT_ROOT=/nfs/roberts/scratch/your_group/your_user/rlm_pilot
mkdir -p "$OUT_ROOT"
sbatch --output="$OUT_ROOT/rlm-%j.log" examples/bouchet_hc_rlm_pilot.sbatch
```

The job creates a unique `sample_JOBID` directory under `OUT_ROOT`, extracts
the BAM and index there, checks the BAM, and writes the full RLM entropy and
single-read BEDs there. `--version-check false` avoids the interactive update
question during a batch job. The archive and reference are read only.

## 2. v4 on a small RLM sample

Find an existing Conda environment with both NumPy and pandas before this job.
Activate it to verify the imports; the RLM build environment may not include
these Python packages.

```bash
source /apps/software/system/software/miniconda/26.9.4/etc/profile.d/conda.sh
conda activate /path/to/existing/python_environment
python -c 'import numpy, pandas; print(numpy.__version__, pandas.__version__)'
```

Submit the v4 pilot from a checkout located on scratch. `INPUT_BED` is the full
RLM entropy BED from step 1. The job copies its header and first 200 windows
to a new scratch directory, then runs 999 conditional simulations per unique
depth and four-CpG methylation signature. It does not change the source BED.

```bash
export INPUT_BED=/nfs/roberts/scratch/your_group/your_user/rlm_pilot/sample_JOBID/output_entropy.bed
export PYTHON_CONDA_ENV=/path/to/existing/python_environment
export OUT_ROOT=/nfs/roberts/scratch/your_group/your_user/rlm_pilot
sbatch --output="$OUT_ROOT/v4-%j.log" examples/bouchet_v4_pilot.sbatch
```

Inspect the RLM and v4 Slurm logs and the number of output rows before
expanding to a full file or the other biological replicates. Keep replicates
separate. The first-200-row pilot checks compatibility and speed, not genome
wide calibration or a DMV biological result.
