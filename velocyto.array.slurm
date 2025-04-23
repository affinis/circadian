#!/bin/bash

#SBATCH --job-name=velocyto-array
#SBATCH --partition=cpu
#SBATCH -N 2
#SBATCH -n 80
#SBATCH --array=1-2
#SBATCH --ntasks-per-node=40
#SBATCH --output=log/velocyto_array_%A_%a.out
#SBATCH --error=log/velocyto_array_%A_%a.err

#echo "`which samtools`"

. ~/.bashrc_mini
conda activate velocyto

VELOCYTO=/lustre/home/miniconda3/envs/velocyto/bin/velocyto
ID=$(cat id.lst | sed -n ${SLURM_ARRAY_TASK_ID}p)
CR_OUTS=/lustre/home/acct-medll/medll/data/cellranger_out/${ID}/outs/per_sample_outs/${ID}/count
echo $CR_OUTS

if [ ! -f $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv ]
then
	gzip -d $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv.gz
fi

#module load samtools

ulimit -n 2048

velocyto run -b $CR_OUTS/sample_filtered_feature_bc_matrix/barcodes.tsv -o $CR_OUTS/velocyto/ -m /lustre/home/acct-medll/medll/ref/repeat_mask/hg38_rmsk.gtf $CR_OUTS/sample_alignments.bam /lustre/home/acct-medll/medll/ref/cellranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf --samtools-memory 800 --samtools-threads 32
