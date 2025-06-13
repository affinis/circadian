#!/bin/bash

# upstream: <Rscript>generateSeacellInput.R <R>circadian_core.R:generateAnnotationFile
# downstream:
# dependency: calculateMetacellsWithSeacells.modforSrtConverted.py
# caller:

WRAPPER_PATH=/tmpdata/LyuLin/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell

ANNOTATION_FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv

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
