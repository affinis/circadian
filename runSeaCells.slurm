#!/bin/bash

#SBATCH --job-name=SeaCells
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --output=log/seacells.out
#SBATCH --error=log/seacells.err
#SBATCH --open-mode=append

echo "------start----------`date`----------start-------"

. ~/.bashrc_mini

conda activate seacells

for sample in `ls /lustre/home/acct-medll/medll/data/cellranger_out | grep 'TF'`
do
	mkdir /lustre/home/acct-medll/medll/data/cellranger_out/$sample/outs/per_sample_outs/$sample/count/seacells
	cd /lustre/home/acct-medll/medll/data/cellranger_out/$sample/outs/per_sample_outs/$sample/count/sample_filtered_feature_bc_matrix
	if [ ! -f barcodes.tsv.gz ]
	then
		gzip barcodes.tsv
	fi
	cd ..
	prefix=`sed -e 's/_1/_CT9/g' -e 's/_2/_CT13/g' -e 's/_3/_CT17/g' -e 's/_4/_CT21/g' -e 's/_5/_CT1/g' -e 's/_6/_CT5/g' <<< $sample`
	python3 /lustre/home/acct-medll/medll/script/pySC_core/calculateMetacellsWithSeacells.py -i sample_filtered_feature_bc_matrix -a /lustre/home/acct-medll/medll/data/analysis/cell.annotations.predicted.and.NI.tsv -o seacells -p $prefix
	cd
done

echo "------done----------`date`----------done-----------"
