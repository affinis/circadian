#!/bin/bash

#SBATCH --job-name=SeaCells
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=1
#SBATCH --output=log/seacells_%a.out
#SBATCH --error=log/seacells_%a.err

echo "------start----------`date`----------start-------"


#conda init
source activate seacells

FILE_PATH=/dssg/home/acct-medll/medll/analysis/cellranger
SEACELL_PY=/dssg/home/acct-medll/medll/script/pySC_core/calculateMetacellsWithSeacells.modforSrtConverted.py
ANNOTATION=/dssg/home/acct-medll/medll/analysis/circadian/R/cell.annotations.manual1_2_NI.predicted2.modforSeacells.tsv
proportion=0.08
min_gene=200

echo "`which python3`"

for ID in `cat id.seacells.lst`
do
	cd $FILE_PATH/$ID/outs/per_sample_outs/$ID/count/preparation_seacells
	for type_dir in `ls`
	do
		cd $FILE_PATH/$ID/outs/per_sample_outs/$ID/count/preparation_seacells/$type_dir
		if [ ! -d seacells_$proportion ]
		then
			mkdir seacells_$proportion
		fi
		cd ..
		prefix=`sed 's/_filtered_bc_matrix//g' <<< $type_dir`
		python3 $SEACELL_PY -i $type_dir -a $ANNOTATION -o $type_dir/seacells_$proportion -p $prefix -r $proportion -m $min_gene
	done
done


echo "------done----------`date`----------done-----------"
