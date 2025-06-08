#!/bin/bash

set -e
set -u

INPUT=/tmpdata/LyuLin/data/RNASeq
OUTPUT=/tmpdata/LyuLin/analysis/circadian/RNASeq
ref=/tmpdata/LyuLin/ref/bwa/human/genome.fa

for sample in `ls $INPUT| cut -d '.' -f 1 |sort | uniq`
do
	mkdir $OUTPUT/$sample
	bwa mem -t 32 $ref $INPUT/$sample.R1.fq.gz $INPUT/$sample.R2.fq.gz | samtools view -bh - > $OUTPUT/$sample/$sample.bam
        samtools sort $OUTPUT/$sample/$sample.bam -o $OUTPUT/$sample/$sample.sorted.bam
	samtools index $OUTPUT/$sample/$sample.sorted.bam
done	
