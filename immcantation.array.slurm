#!/bin/bash

#SBATCH --job-name=immcantation-array
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --array=1-36
#SBATCH --output=log/immcantation_array_%A_%a.out
#SBATCH --error=log/immcantation_array_%A_%a.err

. ~/.bashrc_mini
#conda activate immcantation

echo "`which samtools`"

WRAPPER=/lustre/home/acct-medll/medll/data/immcantation_bkp/001.scripts/011.Wrapper_BT_CRout_single.sh
ID=$(cat id.immcantation.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)
CR_OUTS=/lustre/home/acct-medll/medll/data/cellranger_out/${ID}
echo $CR_OUTS

bash $WRAPPER $CR_OUTS
