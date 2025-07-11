#!/bin/bash

#SBATCH --job-name=SeaCells
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=1-7
#SBATCH --output=log/seacells_%a.out
#SBATCH --error=log/seacells_%a.err

echo "------start----------`date`----------start-------"

#. ~/.bashrc_mini

ID=$(cat id.seacells.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)

#conda init
source activate seacells

FILE_PATH=/dssg/home/acct-medll/medll/analysis/circadian/R/preparation_seacell_EmptyDrop
SEACELL_PY=/dssg/home/acct-medll/medll/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
ANNOTATION=/dssg/home/acct-medll/medll/analysis/circadian/R/droplet.annotation.tsv
proportion=0.08
min_gene=40

echo "`which python3`"

cd $FILE_PATH
#cd $FILE_PATH

if [ ! -d ${ID}_filtered_bc_matrix/seacells_$proportion ]
then
	mkdir ${ID}_filtered_bc_matrix/seacells_$proportion
fi

python3 $SEACELL_PY -i ${ID}_filtered_bc_matrix -a $ANNOTATION -o ${ID}_filtered_bc_matrix/seacells_$proportion -p $ID -r $proportion -m $min_gene


echo "------done----------`date`----------done-----------"
