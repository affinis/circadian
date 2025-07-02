#!/bin/bash

# upstream: <Rscript>generateSeacellInputForEmptyDrop.R
# downstream:
# dependency: calculateMetacellsWithSeacells.modforSrtConverted.py
# caller:

WRAPPER_PATH=/tmpdata/LyuLin/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell_EmptyDrop

ANNOTATION_FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/droplet.annotation.tsv

echo "------start----------`date`----------start-------"

conda init
conda activate seacells

for sample in `ls $FILE_PATH`
do
	if [ ! -d $FILE_PATH/$sample/seacells ]
	then
		mkdir $FILE_PATH/$sample/seacells
	fi
	cd $FILE_PATH/
	prefix=$sample
	python3 $WRAPPER_PATH -i $sample -a $ANNOTATION_FILE_PATH -o $sample/seacells -p $prefix -m 40 -r 0.04
	cd ..
done

echo "------done----------`date`----------done-----------"
