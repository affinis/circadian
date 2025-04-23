#!/bin/bash

#SBATCH --job-name=SCTransform
#SBATCH --partition=192c6t
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --output=log/SCTransform_%a.out
#SBATCH --error=log/SCTransform_%a.err

. ~/.bashrc_mini

conda activate Seurat5

Rscript --vanilla /lustre/home/acct-medll/medll/script/Rscript/runSCTransform.R
