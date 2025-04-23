#!/bin/bash

set -e
set -u

cp circadian_sample_summary.tsv circadian_sample_summary.tsv.backup

for sample in `ls /lustre/home/acct-medll/medll/data/cellranger_out`
do
	if [ `grep $sample circadian_sample_summary.tsv | wc -l` -eq 0 ]
	then
		echo "Adding sample: $sample"
		cat /lustre/home/acct-medll/medll/data/cellranger_out/$sample/outs/per_sample_outs/$sample/web_summary.html | grep -o -f html.pat | sed -e 's/{"table":{"header"://g' -e 's/}//g' | sed 's/"rows"://g' | sed 's/,\[/\t\[/' | cut -f 2 | sed -e 's/\[//g' -e 's/\]//g' | sed 's/","/\n/g' | sed 's/"//g' | sed "1i\\$sample" > log/${sample}.values
		paste circadian_sample_summary.tsv log/${sample}.values > circadian_sample_summary.new.tsv
		mv circadian_sample_summary.new.tsv circadian_sample_summary.tsv
	else
		continue
	fi
done
echo "Finished!"

#cat web_summary.html | grep -o -f ../../../../../../script/html.pat |sed -e 's/{"table":{"header"://g' -e 's/}//g' | sed 's/"rows"://g'
#cat material | sed 's/,\[/\t\[/' | cut -f 2 | sed -e 's/\[//g' -e 's/\]//g' | sed 's/","/\n/g' | sed 's/"//g'
