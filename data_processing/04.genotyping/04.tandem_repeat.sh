#!/bin/bash

# This code genotype FAM193B-CCG repeat

ulimit -S -n 100000

outpath="/path/to/tandem_repeat/"
ref="/path/to/GRCh38.fa"

# HipSTR
[[ ! -d ${outpath}/HipSTR ]] && mkdir -p ${outpath}/HipSTR

~/Software/HipSTR/HipSTR \
	--bam-files input/HipSTR_bam_files.txt \
	--fasta $ref \
	--regions input/FAM193B_STR.HipSTR.bed \
	--str-vcf  ${outpath}/HipSTR/HipSTR.HK.vcf.gz \
	--viz-out ${outpath}/HipSTR/HipSTR.HK.aln_viz.gz \
	--output-pls \
	--output-filters > ${outpath}/HipSTR/HipSTR.HK.log 2>&1

# GangSTR
[[ ! -d ${outpath}/GangSTR/HK ]] && mkdir -p ${outpath}/GangSTR/HK

parallel -j 64 "

sample=\$(basename {1})
sample=\${sample%%.markdup.sort.cram}

GangSTR \
--bam {1} \
--ref $ref \
--regions input/FAM193B_STR.GangSTR.bed \
--grid-threshold 250 \
--readlength 150 \
--bam-samps \${sample} \
--out ${outpath}/GangSTR/HK/\${sample} \
> ${outpath}/GangSTR/HK/\${sample}.log 2>&1

bgzip -f ${outpath}/GangSTR/HK/\${sample}.vcf
tabix -f ${outpath}/GangSTR/HK/\${sample}.vcf.gz

" :::: input/HipSTR_bam_files.txt

ls ${outpath}/GangSTR/HK/*.vcf.gz > ${outpath}/GangSTR.HK.vcf_list.txt

mergeSTR \
--vcfs-list ${outpath}/GangSTR.HK.vcf_list.txt \
--vcftype gangstr \
--out ${outpath}/GangSTR/GangSTR.HK \
> ${outpath}/GangSTR/GangSTR.HK.log 2>&1

bgzip -f ${outpath}/GangSTR/GangSTR.HK.vcf
tabix -f ${outpath}/GangSTR/GangSTR.HK.vcf.gz

# ExpansionHunter
[[ ! -d ${outpath}/ExpansionHunter/HK ]] && mkdir -p ${outpath}/ExpansionHunter/HK

parallel -j 64 "

sample=\$(basename {1})
sample=\${sample%%.markdup.sort.cram}

ExpansionHunter \
--reads {1} \
--reference $ref \
--variant-catalog input/FAM193B_STR.ExpansionHunter.json \
--output-prefix ${outpath}/ExpansionHunter/HK/\${sample} \
--threads 1 \
> ${outpath}/ExpansionHunter/HK/\${sample}.log 2>&1

bgzip -f ${outpath}/ExpansionHunter/HK/\${sample}.vcf
tabix -f ${outpath}/ExpansionHunter/HK/\${sample}.vcf.gz

" :::: input/HipSTR_bam_files.txt

ls ${outpath}/ExpansionHunter/HK/*.vcf.gz > ${outpath}/ExpansionHunter.HK.vcf_list.txt

mergeSTR \
--vcfs-list ${outpath}/ExpansionHunter.HK.vcf_list.txt \
--vcftype eh \
--out ${outpath}/ExpansionHunter/ExpansionHunter.HK \
> ${outpath}/ExpansionHunter/ExpansionHunter.HK.log 2>&1

bgzip ${outpath}/ExpansionHunter/ExpansionHunter.HK.vcf
tabix -f ${outpath}/ExpansionHunter/ExpansionHunter.HK.vcf.gz

# Ensembl TR
EnsembleTR \
--vcfs ${outpath}/HipSTR/HipSTR.HK.vcf.gz,${outpath}/GangSTR/GangSTR.HK.vcf.gz,${outpath}/ExpansionHunter/ExpansionHunter.HK.vcf.gz \
--ref $ref \
--out ${outpath}/EnsembleTR.FAM193B_STR.HK.vcf

bgzip -f ${outpath}/EnsembleTR.FAM193B_STR.HK.vcf
tabix -f ${outpath}/EnsembleTR.FAM193B_STR.HK.vcf.gz

# extract copy number for analysis
bcftools query -f '[%SAMPLE\t%NCOPY\t%SCORE\n]' ${outpath}/EnsembleTR.FAM193B_STR.HK.vcf.gz | \
sort -k1,1 > ${outpath}/EnsembleTR.FAM193B_STR.HK.genotype.txt