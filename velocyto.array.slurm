#!/bin/bash

#SBATCH --job-name=velocyto-array
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --array=1-6
#SBATCH --exclusive
#SBATCH --output=log/velocyto_array_%a.out
#SBATCH --error=log/velocyto_array_%a.err

#echo "`which samtools`"

#. ~/.bashrc_mini
conda init
conda activate velocyto

export PATH=/dssg/home/acct-medll/medll/miniconda3/envs/velocyto/bin:$PATH

VELOCYTO=~/miniconda3/envs/velocyto/bin/velocyto
ID=$(cat id.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)
CR_OUTS=~/analysis/cellranger/${ID}/outs/per_sample_outs/${ID}/count
REF=/dssg/home/acct-medll/medll/ref/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf
MASK=/dssg/home/acct-medll/medll/ref/hg38_rmsk.gtf
echo $CR_OUTS

if [ ! -f $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv ]
then
	gzip -d $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
fi

#module load samtools

ulimit -n 2048

$VELOCYTO run -b $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv -o $CR_OUTS/velocyto/ -m $MASK $CR_OUTS/sample_alignments.bam --samtools-memory 2048 $REF 
