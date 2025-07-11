#!/bin/bash

#SBATCH --job-name=SeaCells
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=1-18
#SBATCH --output=log/seacells_%a.out
#SBATCH --error=log/seacells_%a.err

echo "------start----------`date`----------start-------"

#. ~/.bashrc_mini

ID=$(cat id.seacells.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)

#conda init
source activate seacells

FILE_PATH=/dssg/home/acct-medll/medll/analysis/circadian/R/preparation_seacell_EmptyDrop
SEACELL_PY=/dssg/home/acct-medll/medll/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
ANNOTATION=/dssg/home/acct-medll/medll/analysis/circadian/R/cell.annotations.manual1_2_NI.predicted2.modforSeacells.tsv
proportion=0.08
min_gene=200

echo "`which python3`"

cd $FILE_PATH/$ID/outs/per_sample_outs/$ID/count/
#cd $FILE_PATH

if [ ! -f sample_filtered_feature_bc_matrix/barcodes.tsv.gz ]
then
	gzip sample_filtered_feature_bc_matrix/barcodes.tsv
fi

if [ `zcat sample_filtered_feature_bc_matrix/barcodes.tsv.gz | grep '^TF' | wc -l` -eq 0 ]
then
	cp sample_filtered_feature_bc_matrix/barcodes.tsv.gz sample_filtered_feature_bc_matrix/barcodes.tsv.gz.backup
	zcat sample_filtered_feature_bc_matrix/barcodes.tsv.gz | sed "s/^/TF_${ID}_/g" > sample_filtered_feature_bc_matrix/barcodes.tsv
	gzip -f sample_filtered_feature_bc_matrix/barcodes.tsv
fi

if [ ! -d seacells_$proportion ]
then
	mkdir seacells_$proportion
fi

python3 $SEACELL_PY -i sample_filtered_feature_bc_matrix -a $ANNOTATION -o seacells_$proportion -p $ID -r $proportion -m $min_gene


echo "------done----------`date`----------done-----------"
