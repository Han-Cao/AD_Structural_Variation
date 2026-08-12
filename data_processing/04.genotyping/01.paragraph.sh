#!/bin/bash

vcf_sv="/path/to/All.genotype.filter.population.AC0.clean.vcf.gz"
ref="/path/to/GRCh38.fa"
tmp_dir="/tmp/"
depth_all="/path/to/sample_depth.txt"
threads=16

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

bam_raw=${inpath}/${sample}.markdup.sort.cram
bam_dedup=${tmp_dir}/${sample}.dedup.sort.cram

# Paragraph doesn't remove duplicated reads, so we have to do it manually
if [[ ! -e ${bam_dedup}.crai ]]; then 
	samtools view -@ 7 -F 1024 -h -C -T $ref -o $bam_dedup $bam_raw
	samtools index -@ $threads $bam_dedup
fi

# prepare manifest
depth=$(grep "^${sample}" $depth_all | cut -f 2)
max_depth=$(($depth * 20))
echo -e "id\tpath\tread length\tdepth" > input/manifest/${sample}.manifest.txt
echo -e "${sample}\t${bam_dedup}\t150\t${depth}" >> input/manifest/${sample}.manifest.txt

# run paragraph
/path/paragraph-v2.4a/bin/multigrmpy.py \
-i $vcf_sv \
-m input/manifest/${sample}.manifest.txt \
-r $ref \
-o ${outpath}/${sample} \
-t $threads \
-M $max_depth

