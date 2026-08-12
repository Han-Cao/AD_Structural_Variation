#!/bin/bash

# GATK 4.4.0.0

refPre="/path/to/GRCh38.fa"

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

[[ ! -d $outpath$sample ]] && mkdir -p $outpath$sample

ulimit -c unlimited

# get input bam files
INPUT_ARGUMENT=$(ls ${outpath}/${sample}/*namesort.bam | xargs -I {} echo '-I {}' | tr '\n' ' ')

#If the job is cancelled by cluster, remove tmp files from last run
if [[ -d "${outpath}/${sample}/${sample}.markdup.sort.bam.parts" ]]; then
	rm -r "${outpath}/${sample}/${sample}.markdup.sort.bam.parts"
fi

#markduplicates
gatk --java-options "-Xmx160G" MarkDuplicatesSpark \
${INPUT_ARGUMENT} \
-O "${outpath}/${sample}/${sample}.markdup.sort.bam" 

# conevert to cram
samtools view -C -T $refPre "${outpath}/${sample}/${sample}.markdup.sort.bam"  -o "${outpath}/${sample}.markdup.sort.cram"
samtools index "${outpath}/${sample}.markdup.sort.cram"

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
	# clean intermediate files
	rm -r "${outpath}/${sample}"
else
    touch "${error}/${sample}.error"
fi
