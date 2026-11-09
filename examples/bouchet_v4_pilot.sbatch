#!/usr/bin/env bash
# Run v4 on the first 200 rows of a full RLM BED, inside scratch.
#SBATCH --job-name=rlm-v4-pilot
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -euo pipefail

: "${SLURM_JOB_ID:?Submit this script with sbatch}"
: "${INPUT_BED:?Set the full RLM entropy BED path}"
: "${PYTHON_CONDA_ENV:?Set a Conda environment with NumPy and pandas}"
: "${OUT_ROOT:?Set a scratch output directory}"
[[ -d "$OUT_ROOT" ]] || { echo "Create OUT_ROOT on scratch before sbatch" >&2; exit 2; }
OUT_ROOT=$(realpath -e -- "$OUT_ROOT")
[[ "$OUT_ROOT" == /nfs/roberts/scratch/* ||
   "$OUT_ROOT" == /gpfs/gibbs/scratch/* ]] || {
    echo "OUT_ROOT must resolve inside scratch" >&2; exit 2;
}
[[ -f "$INPUT_BED" ]] || { echo "INPUT_BED is missing" >&2; exit 2; }
conda_sh=${CONDA_SH:-/apps/software/system/software/miniconda/26.9.4/etc/profile.d/conda.sh}
[[ -f "$conda_sh" ]] || { echo "Conda initialization script is missing" >&2; exit 2; }

source "$conda_sh"
conda activate "$PYTHON_CONDA_ENV"
python -c 'import numpy, pandas; print("NumPy", numpy.__version__, "pandas", pandas.__version__)'

tool_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
run_dir="${OUT_ROOT%/}/v4_${SLURM_JOB_ID}"
mkdir "$run_dir"  # Fails instead of overwriting an earlier run.
head -n 201 "$INPUT_BED" > "$run_dir/rlm_first200.bed"
[[ $(wc -l < "$run_dir/rlm_first200.bed") -gt 1 ]] || {
    echo "RLM BED has no data rows" >&2; exit 1;
}

PYTHONPATH="$tool_root/src${PYTHONPATH:+:$PYTHONPATH}" \
    python -m methylation_entropy.rlm_entropy_v4 \
        --input "$run_dir/rlm_first200.bed" \
        --output "$run_dir/entropy_cpg_null.tsv" \
        --nsims 999 --seed 42
echo "v4 pilot output: $run_dir/entropy_cpg_null.tsv"
