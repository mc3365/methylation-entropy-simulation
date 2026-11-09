#!/usr/bin/env bash
# Extract one indexed HC BAM and run RLM. Submit with all paths supplied in the
# environment and an explicit --output path inside scratch (see docs).
#SBATCH --job-name=rlm-hc-pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00
set -euo pipefail

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${ARCHIVE:?Set the existing BAM tar.gz path}"
: "${BAM_MEMBER:?Set the BAM member path inside the archive}"
: "${REFERENCE:?Set the matching reference FASTA path}"
: "${RLM_BIN:?Set the existing RLM executable path}"
: "${RLM_CONDA_ENV:?Set the existing RLM Conda environment path}"
: "${OUT_ROOT:?Set a scratch output directory}"

[[ -d "$OUT_ROOT" ]] || { echo "Create OUT_ROOT on scratch before sbatch" >&2; exit 2; }
OUT_ROOT=$(realpath -e -- "$OUT_ROOT")
[[ "$OUT_ROOT" == /nfs/roberts/scratch/* ||
   "$OUT_ROOT" == /gpfs/gibbs/scratch/* ]] || {
    echo "OUT_ROOT must resolve inside scratch" >&2; exit 2;
}
[[ "$BAM_MEMBER" == *.bam && "$BAM_MEMBER" != /* &&
   "$BAM_MEMBER" != *../* && "$BAM_MEMBER" != ../* ]] || {
    echo "BAM_MEMBER must be a relative .bam path without .." >&2; exit 2;
}
[[ -f "$ARCHIVE" && -f "$REFERENCE" && -x "$RLM_BIN" ]] || {
    echo "Archive, reference, or RLM executable is missing" >&2; exit 2;
}
conda_sh=${CONDA_SH:-/apps/software/system/software/miniconda/26.9.4/etc/profile.d/conda.sh}
[[ -f "$conda_sh" ]] || { echo "Conda initialization script is missing" >&2; exit 2; }

module load SAMtools/1.21-GCC-13.3.0
source "$conda_sh"
conda activate "$RLM_CONDA_ENV"
command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 2; }

sample=${BAM_MEMBER##*/}
sample=${sample%.bam}
run_dir="${OUT_ROOT%/}/${sample}_${SLURM_JOB_ID}"
mkdir "$run_dir"  # Fails instead of overwriting an earlier run.
echo "Run directory: $run_dir"

tar -xzf "$ARCHIVE" -C "$run_dir" "$BAM_MEMBER" "${BAM_MEMBER}.bai"
bam_file="$run_dir/$BAM_MEMBER"
[[ -f "$bam_file" && -f "${bam_file}.bai" ]] || {
    echo "Expected BAM and index were not extracted" >&2; exit 2;
}
samtools quickcheck -v "$bam_file"

"$RLM_BIN" \
    --version-check false \
    --bam "$bam_file" \
    --reference "$REFERENCE" \
    --mode PE \
    --aligner bsmap \
    --score entropy \
    --coverage 10 \
    --output_entropy "$run_dir/output_entropy.bed" \
    --output_single_read "$run_dir/output_single_read_info.bed"

[[ -s "$run_dir/output_entropy.bed" ]] || {
    echo "RLM did not produce a nonempty entropy BED" >&2; exit 1;
}
echo "RLM entropy BED: $run_dir/output_entropy.bed"
