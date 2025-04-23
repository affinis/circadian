#!/bin/bash

#SBATCH --job-name=MetaCycle-array
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --array=1-26
#SBATCH --exclusive
#SBATCH --output=log/MetaCycle_array_%A_%a.out
#SBATCH --error=log/MetaCycle_array_%A_%a.err

. ~/.bashrc_mini

conda activate Seurat5

Rscript --vanilla /lustre/home/acct-medll/medll/script/Rscript/runMetacyleByPseudoBulk.R ${SLURM_ARRAY_TASK_ID}
