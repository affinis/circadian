#!/bin/bash

#SBATCH --job-name=dorothea
#SBATCH --partition=192c6t
#SBATCH -N 1
#SBATCH -n 96
#SBATCH --output=log/dorothea.out
#SBATCH --error=log/dorothea.err
#SBATCH --open-mode=append

echo "------start----------`date`----------start-------"

. ~/.bashrc_mini

conda activate Seurat5

Rscript --vanilla ~/script/Rscript/dorothea.R

echo "------done----------`date`----------done-----------"
