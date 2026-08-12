#!/bin/bash

# GATK 4.4.0.0

refPre="/path/to/GRCh38.fa"
dbsnp="/path/to/dbsnp_153.GRCh38.vcf.gz"
omni="/path/to/1000G_omni2.5.hg38.PASS.vcf.gz"
hapmap="/path/to/hapmap_3.3.hg38.vcf.gz"
g1000="/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5


finish=finish/BQSR/

[[ ! -d $outpath ]] && mkdir -p ${outpath}

gatk --java-options "-Xmx10G" BaseRecalibrator \
-VS LENIENT \
--known-sites $g1000 \
--known-sites $mills \
--known-sites $dbsnp \
-R $refPre \
-O "${outpath}/${sample}.markdup.sort.cram.recalibration_report.grp" \
-I "${inpath}/${sample}.markdup.sort.cram"

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi