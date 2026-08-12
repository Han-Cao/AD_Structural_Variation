#!/bin/bash
#SBATCH -p cpu -N 1 -n 40 -o log/VQSR/clean_VQSR.log

refPre="/path/to/GRCh38.fa"
vcf="/path/to/VQSR/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.vcf.gz"

outpath="/path/to/VQSR/"

# split to biallelic, remove duplicate variants, norm INDEL
# IMPORTANT: always use As_FilterStatus to keep high quality variants
bcftools norm -Ou -f $refPre --do-not-normalize -m - $vcf | \
bcftools norm -Ou -f $refPre -d exact -Ob -o ${outpath}/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.bcf
bcftools index ${outpath}/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.bcf

# only keep PASS variants
bcftools view \
-i 'AS_FilterStatus == "PASS"' \
-Oz -o ${outpath}/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.PASS.vcf.gz \
${outpath}/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.bcf
tabix -p vcf ${outpath}/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.PASS.vcf.gz