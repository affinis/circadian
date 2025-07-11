#!/bin/bash

#SBATCH --job-name=SeaCells
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --output=log/seacells.out
#SBATCH --error=log/seacells.err
#SBATCH --open-mode=append

echo "------start----------`date`----------start-------"

#. ~/.bashrc_mini

conda init
conda activate seacells

FILE_PATH=/dssg/home/acct-medll/medll/analysis/circadian/R/preparation_seacell
SEACELL_PY=/dssg/home/acct-medll/medll/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
ANNOTATION=/dssg/home/acct-medll/medll/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv

for sample in `ls $FILE_PATH | grep 'CT'`
do

	if [ ! -d $FILE_PATH/$sample/seacells ]
	then
		mkdir $FILE_PATH/$sample/seacells
	fi
	
	cd $FILE_PATH/$sample
	
	if [ -f $sample/genes.tsv.gz ]
	then
		mv $sample/genes.tsv.gz $sample/features.tsv.gz
	fi

	if [ -f $sample/features.tsv ]
	then
		gzip $sample/features.tsv
	fi

	if [ `zgrep "Gene Expression" $sample/features.tsv.gz | wc -l` -eq 0 ]
	then
		zcat $sample/features.tsv.gz | sed 's/$/\tGene Expression/g' > $sample/features.tsv
		mv $sample/features.tsv.gz $sample/features.tsv.gz.backup
		gzip $sample/features.tsv
	fi
	#if [ ! -f barcodes.tsv.gz ]
	#then
	#		gzip barcodes.tsv
	#fi
	#cd ..
	#prefix=`sed -e 's/_1/_CT9/g' -e 's/_2/_CT13/g' -e 's/_3/_CT17/g' -e 's/_4/_CT21/g' -e 's/_5/_CT1/g' -e 's/_6/_CT5/g' <<< $sample`
	prefix=$sample
	python3 $SEACELL_PY -i $sample -a $ANNOTATION -o seacells -p $prefix -r 0.04
	cd
done

echo "------done----------`date`----------done-----------"
