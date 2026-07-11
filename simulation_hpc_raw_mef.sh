#!/bin/bash
#SBATCH --job-name=Simulation_raw
#SBATCH --out="slurm-%j.out"
#SBATCH --time=1-00:00:00
#SBATCH --partition=bigmem
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem-per-cpu=128G
#SBATCH --mail-user=meiying.cui@yale.edu
#SBATCH --mail-type=ALL


#this version use lookup table for the minimal value

module purge
ml miniconda
conda activate /vast/palmer/pi/smith/mc3365/env/python_env

input_dir="/gpfs/gibbs/project/zsmith/mc3365/wgbs/rlm/mef_pdr"
main_dir="/vast/palmer/scratch/zsmith/mc3365/rlm_analysis/simulation/data"
code_dir="/vast/palmer/scratch/zsmith/mc3365/rlm_analysis/simulation/code/"


# Process each passage
for passage in _1 0 5; do
    echo "Processing passage ${passage}"
    
    python ${code_dir}/entropy_simulation_raw_multi_4.py \
        --input ${input_dir}/mef_entropy_p${passage}.bed \
        --output ${main_dir}/entropy_raw_simulation_p${passage}.txt \
        --lookup ${main_dir}/minimum_entropy_lookup.json \
        --nsims 1000
done
