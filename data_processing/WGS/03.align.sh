#!/bin/bash

# bwa-mem2 2.1
# GATK 4.4.0.0

refPre="/path/to/GRCh38.fa"

sample=$1
inpath=$2
outpath=$3
readgroup=$4
finish=$5
error=$6

[[ ! -d $outpath$sample ]] && mkdir -p $outpath$sample

ls ${inpath}/${sample}/*_1_paired.fq.gz | while read fastq
do
	#get inupt file
	fastq=${fastq%_1_paired.fq.gz}
	fq1=${fastq}_1_paired.fq.gz
	fq2=${fastq}_2_paired.fq.gz
	name=$(basename $fastq)
	#extract read group from readgroup.csv
	rg=$(grep $name $readgroup | cut -d ',' -f 6)
	#If the job is cancelled by cluster, remove tmp files from last run
	if [[ -e "${outpath}/${sample}/${name}.namesort.bam.tmp.0000.bam" ]]; then
		rm "${outpath}/${sample}/${name}.namesort.bam.tmp"*
	fi
	#MarkDuplicatesSpark required input bam sorted by read name not by coordinates
	#MarkDuplicatesSpark will sort bam by coordinates and mark duplicates
	bwa-mem2 mem -t 32 -M -R "$rg" $refPre $fq1 $fq2 | samtools view -uhS - | \
	samtools sort -n -m 4G -@ 6 -o "${outpath}/${sample}/${name}.namesort.bam"
done

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi
