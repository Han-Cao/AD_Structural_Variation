#!/bin/bash

#SBATCH -p cpu -N 1 -n 40 -o log/sample_qc/make_bfile.log

vcf="/path/to/VQSR/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.norm.dedup.PASS.vcf.gz"
vcf_raw="/path/to/VQSR/HK_WGS_update.hardfilter.VQSR_SNP_INDEL.vcf.gz"
outpath="/path/to/sample_qc/"

[[ ! -d $outpath ]] && mkdir $outpath
[[ ! -d ${outpath}/stats ]] && mkdir ${outpath}/stats


# generate raw bfile
plink2 --vcf $vcf \
--make-bed \
--out ${outpath}/HK_WGS_update.PASS.raw

# extract autosome biallelic SNPs following HWE for QC
# these filters are for QC steps only, we have additional filters in downstream analysis
plink2 --vcf $vcf_raw \
--max-alleles 2 \
--var-filter \
--autosome \
--snps-only \
--hwe 5e-8 \
--make-bed \
--out ${outpath}/HK_WGS_update.PASS.biallele.snp.hwe.autosome

# filter by genotyping rate and MAF
plink2 --bfile ${outpath}/HK_WGS_update.PASS.biallele.snp.hwe.autosome \
--geno 0.01 \
--maf 0.01 \
--make-bed \
--out ${outpath}/HK_WGS_update.PASS.biallele.snp.hwe.autosome.geno1.maf1

# generate summary
plink2 --bfile ${outpath}/HK_WGS_update.PASS.raw \
--sample-counts \
--het \
--out ${outpath}/stats/HK_WGS_update.PASS.biallele.snp.hwe.autosome.stats
