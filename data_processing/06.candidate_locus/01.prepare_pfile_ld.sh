#!/bin/bash

bfile=/path/to/hk/bfile
outpath=/path/to/pfile_LD

# Estimate LD between FAM193B-CCG and nearby SNPs
vcf_small_var="/path/to/small_variants.vcf.gz"
outpath="/path/to/LD/"

# 1. prepare pfile for association analysis
# add DS to VCF for repeat
python add_ds_for_regression.py \
-i ${vcf_small_var} \
-o ${outpath}/repeat.dosage.regression.vcf.gz\
-d /path/to/dosage.txt

# concat repeat with nearby small variants
bcftools concat $vcf_small_var ${outpath}/repeat.dosage.regression.vcf.gz | \
bcftools view -f PASS -Oz -o ${outpath}/repeat.dosage.concat.regression.vcf.gz

# convert to pfile with dosage
plink2 --vcf ${outpath}/repeat.dosage.concat.regression.vcf.gz dosage=DS \
--keep /path/to/clean_samples.txt \
--double-id \
--sort-vars \
--make-pgen \
--out ${outpath}/repeat.dosage.concat.regression

# 2. prepare pfile for LD analysis
# add DS to VCF for repeat
python add_ds_for_ld.py \
-i ${vcf_small_var} \
-o ${outpath}/repeat.dosage.ld.vcf.gz\
-d /path/to/dosage.txt \
--max-cap ${QUANTILE_97.5} \
--min-cap ${QUANTILE_2.5}

# concat repeat with nearby small variants
bcftools concat $vcf_small_var ${outpath}/repeat.dosage.ld.vcf.gz | \
bcftools view -f PASS -Oz -o ${outpath}/repeat.dosage.concat.ld.vcf.gz

# convert to pfile with dosage
plink2 --vcf ${outpath}/repeat.dosage.concat.ld.vcf.gz dosage=DS \
--keep /path/to/clean_samples.unrelated.txt \
--double-id \
--sort-vars \
--make-pgen \
--out ${outpath}/repeat.dosage.concat.ld

# LD analysis
plink2 --pfile ${outpath}/repeat.dosage.concat.ld \
--r-unphased ref-based \
--ld-snp-list /path/to/FAM193B_CCG_var_list.txt \
--ld-window-r2 0 \
--out ${outpath}/repeat.dosage.concat.LD.r

plink2 --pfile ${outpath}/repeat.dosage.concat.ld \
--r2-unphased ref-based \
--ld-snp-list /path/to/FAM193B_CCG_var_list.txt \
--ld-window-r2 0 \
--out ${outpath}/repeat.dosage.concat.LD.r2