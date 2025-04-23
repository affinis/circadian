#!/bin/bash

#SBATCH --job-name=MetaCycle-array
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --array=1-26
#SBATCH --exclusive
#SBATCH --output=log/MetaCycle_array_%a.out
#SBATCH --error=log/MetaCycle_array_%a.err
#SBATCH --open-mode=append

. ~/.bashrc_mini

conda activate Seurat5

Rscript --vanilla /lustre/home/acct-medll/medll/script/Rscript/runMetacyleByPseudoBulk_regular.R ${SLURM_ARRAY_TASK_ID}
