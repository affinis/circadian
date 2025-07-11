#!/bin/bash

# upstream: <Rscript>generateSeacellInput.R <R>circadian_core.R:generateAnnotationFile
# downstream: <R>circadian_core.R:metacell2srt <R>circadian_core.R:readAllMetaCell
# dependency: calculateMetacellsWithSeacells.modforSrtConverted.py
# caller: NSF

WRAPPER_PATH=/tmpdata/LyuLin/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/preparation_seacell_byType

ANNOTATION_FILE_PATH=/tmpdata/LyuLin/analysis/circadian/R/cell.annotation.clean.sct.Azimuth.tsv
proportion=0.08
min_gene=200

echo "------start----------`date`----------start-------"

source activate seacells
n_sample=`ls $FILE_PATH | wc -l`
i=1
for sample in `ls $FILE_PATH`
do
	if [ ! -d $FILE_PATH/$sample/seacells_$proportion ]
	then
		mkdir $FILE_PATH/$sample/seacells_$proportion
	fi
	cd $FILE_PATH
	prefix=$sample
	python3 $WRAPPER_PATH -i $sample -a $ANNOTATION_FILE_PATH -o $sample/seacells_$proportion -p $prefix -r $proportion -m $min_gene
	echo "$i/${n_sample} completed !!"
	let i=i+1
done

echo "------done----------`date`----------done-----------"
