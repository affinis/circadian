#!/bin/bash

# upstream: <Rscript>generateSeacellInputForEmptyDrop.R
# downstream:
# dependency: calculateMetacellsWithSeacells.modforSrtConverted.py
# caller:

WRAPPER_PATH=/tmpdata/LyuLin/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell_EmptyDrop

ANNOTATION_FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/droplet.annotation.tsv

echo "------start----------`date`----------start-------"

conda activate seacells

for sample in `ls $FILE_PATH`
do
	mkdir $FILE_PATH/$sample/seacells
	cd $FILE_PATH/$sample/
	prefix=$sample
	python3 $WRAPPER_PATH -i ./ -a $ANNOTATION_FILE_PATH -o seacells -p $prefix
	cd ..
done

echo "------done----------`date`----------done-----------"
