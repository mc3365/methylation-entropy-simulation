#!/usr/bin/env bash
# Run RLM on one already aligned BAM. Extract .bam files from archives first.
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
    echo "Usage: bash examples/run_rlm_entropy_from_bam.sh RLM_BIN BAM REFERENCE.fa MODE ALIGNER OUT_DIR [MIN_COVERAGE]" >&2
    echo "MODE: SE or PE; ALIGNER: bsmap, bismark, segemehl, or gem" >&2
    exit 2
fi

rlm_bin=$1
bam_file=$2
reference=$3
mode=$4
aligner=$5
out_dir=$6
min_coverage=${7:-10}

[[ -x "$rlm_bin" ]] || { echo "RLM executable not found: $rlm_bin" >&2; exit 2; }
[[ -f "$bam_file" ]] || { echo "BAM not found: $bam_file" >&2; exit 2; }
[[ -f "$reference" ]] || { echo "Reference not found: $reference" >&2; exit 2; }
[[ "$mode" == "SE" || "$mode" == "PE" ]] || { echo "MODE must be SE or PE" >&2; exit 2; }
case "$aligner" in
    bsmap|bismark|segemehl|gem) ;;
    *) echo "Unsupported ALIGNER: $aligner" >&2; exit 2 ;;
esac
[[ "$min_coverage" =~ ^[0-9]+$ ]] && (( min_coverage >= 1 && min_coverage <= 1000 )) || {
    echo "MIN_COVERAGE must be an integer from 1 to 1000" >&2
    exit 2
}

mkdir -p "$out_dir"
entropy_file="$out_dir/output_entropy.bed"
read_file="$out_dir/output_single_read_info.bed"
[[ ! -e "$entropy_file" && ! -e "$read_file" ]] || {
    echo "Output exists; choose a new OUT_DIR to avoid overwriting: $out_dir" >&2
    exit 2
}

if command -v samtools >/dev/null 2>&1; then
    samtools quickcheck -v "$bam_file"
fi

"$rlm_bin" \
    --bam "$bam_file" \
    --reference "$reference" \
    --mode "$mode" \
    --score entropy \
    --aligner "$aligner" \
    --coverage "$min_coverage" \
    --output_entropy "$entropy_file" \
    --output_single_read "$read_file"

echo "RLM entropy output: $entropy_file"
echo "Check the header, then run: methylation-entropy-rlm --input '$entropy_file' --output '$out_dir/entropy_cpg_null.tsv'"
