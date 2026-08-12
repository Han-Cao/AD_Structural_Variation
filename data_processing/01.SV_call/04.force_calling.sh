#!/bin/bash

ref_vcf="/scratch/PI/boip/Han/project/SV/represent/All_sample_pipeline_merged.representative.all.vcf.gz"
outpath="/scratch/PI/boip/Han/project/SV/force_calling/"

# minimap2
ls /scratch/PI/boip/Han/project/SV/align/minimap2/*bam | while read bam;
do
    name=$(basename $bam)
    name=${name%.bam}

    # sniffles
    sniffles \
    -i $bam \
    -v ${outpath}/${name}.sniffles.vcf \
    --reference "/scratch/PI/boip/Reference/Human_genome/GRCh38/GRCh38.fa" \
    -t 40 \
    --tandem-repeats "/scratch/PI/boip/Reference/Human_genome/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed" \
    --genotype-vcf $ref_vcf

    # cuteSV
    mkdir -p ${outpath}/${name}.cuteSV.vcf_wd
    cuteSV \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 200 \
    --diff_ratio_merging_DEL 0.5 \
    --genotype \
    -t 40 \
    -L -1 \
    -Ivcf $ref_vcf \
    $bam \
    "/scratch/PI/boip/Reference/Human_genome/GRCh38/GRCh38.fa" \
    ${outpath}/${name}.cuteSV.vcf \
    ${outpath}/${name}.cuteSV.vcf_wd

done


# NGMLR
ls /scratch/PI/boip/Han/project/SV/align/NGMLR/*bam | while read bam;
do
    name=$(basename $bam)
    name=${name%.bam}

    # sniffles
    sniffles \
    -i $bam \
    -v ${outpath}/${name}.sniffles.vcf \
    --reference "/scratch/PI/boip/Reference/Human_genome/GRCh38/GRCh38.fa" \
    -t 40 \
    --tandem-repeats "/scratch/PI/boip/Reference/Human_genome/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed" \
    --genotype-vcf $ref_vcf

    # cuteSV
    mkdir -p ${outpath}/${name}.cuteSV.vcf_wd
    cuteSV \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 200 \
    --diff_ratio_merging_DEL 0.5 \
    --genotype \
    -t 40 \
    -L -1 \
    -Ivcf $ref_vcf \
    $bam \
    "/scratch/PI/boip/Reference/Human_genome/GRCh38/GRCh38.fa" \
    ${outpath}/${name}.cuteSV.vcf \
    ${outpath}/${name}.cuteSV.vcf_wd

done