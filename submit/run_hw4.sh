#!/bin/bash
#SBATCH --job-name=hw4_mix
#SBATCH --output=logs/hw4_%A_%a.out
#SBATCH --error=logs/hw4_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-3

module load R

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  N=100
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  N=1000
else
  N=10000
fi

Rscript R/run_scenario.R $N