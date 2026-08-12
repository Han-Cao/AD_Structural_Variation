#!/bin/bash

bgen_array="/path/to/UKB.array.raw"
pfile_wgs="/path/to/repeat.dosage.concat.regression"
outpath="/path/to/association/"

# 1. prepare snp genotype files for regenie step1
plink2 --bgen ${bgen_array} ref-last \
--make-pgen \
--out ${outpath}/bfile_clean/UKB.array.raw
# this bgen has FID_IID format IID, we manually change the output psam to IID
mv ${outpath}/bfile_clean/UKB.array.raw.psam ${outpath}/bfile_clean/UKB.array.raw.psam.bkp
cut -f 2 -d '_' ${outpath}/bfile_clean/UKB.array.raw.psam.bkp > ${outpath}/bfile_clean/UKB.array.raw.psam

# only keep QC-pass samples and LD prune variants
plink2 --pfile ${outpath}/bfile_clean/UKB.array.raw \
--keep /path/to/clean_samples.txt \
--maf 0.01 --mac 200 --geno 0.01 --hwe 1e-15 \
--indep-pairwise 1000 100 0.9 \
--out ${outpath}/bfile_clean/UKB.array.prune

plink2 --pfile ${outpath}/bfile_clean/UKB.array.raw \
--keep /path/to/clean_samples.txt \
--extract ${outpath}/bfile_clean/UKB.array.prune.prune.in \
--make-pgen fid \
--out ${outpath}/bfile_clean/UKB.array.clean


# 2. regenie step 1
[[ ! -d ${outpath}/regenie_step1/ ]] && mkdir -p ${outpath}/regenie_step1/
regenie \
--step 1 \
--pgen ${pfile_array} \
--keep /path/to/UKB.clinical_case_control.samples.txt \
--phenoFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--phenoCol AD \
--covarFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--covarColList age,sex,center{1:3},pc_{1:20} \
--bt \
--bsize 1000 \
--threads 80 \
--out ${outpath}/regenie_step1/UKB.clinical_case_control

# 3. regenie step 2 
[[ ! -d ${outpath}/regenie_step2/ ]] && mkdir -p ${outpath}/regenie_step2/

# AD associations
regenie \
--step 2 \
--pgen /path/to/repeat.dosage.concat.regression \
--keep /path/to/UKB.clinical_case_control.samples.txt \
--phenoFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--phenoCol AD \
--covarFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--covarColList age,sex,center{1:3},pc_{1:20} \
--bt \
--pred ${outpath}/regenie_step1/UKB.clinical_case_control_pred.list \
--bsize 400 \
--threads 15 \
--out ${outpath}/regenie_step2/UKB.clinical_case_control

# condition on FAM193B-CCG
regenie \
--step 2 \
--pgen /path/to/repeat.dosage.concat.regression \
--keep /path/to/UKB.clinical_case_control.samples.txt \
--phenoFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--phenoCol AD \
--covarFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--covarColList age,sex,center{1:3},pc_{1:20} \
--condition-list /path/to/UKB.condition_tr.txt \
--bt \
--pred ${outpath}/regenie_step1/UKB.clinical_case_control_pred.list \
--bsize 400 \
--threads 15 \
--out ${outpath}/regenie_step2/UKB.clinical_case_control.condtion_tr

# condition on SNPs
regenie \
--step 2 \
--pgen /path/to/repeat.dosage.concat.regression \
--keep /path/to/UKB.clinical_case_control.samples.txt \
--phenoFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--phenoCol AD \
--covarFile /path/to/UKB.clinical_case_control.pheno_covar.txt \
--covarColList age,sex,center{1:3},pc_{1:20} \
--condition-list /path/to/UKB.condition_snp.txt \
--bt \
--pred ${outpath}/regenie_step1/UKB.clinical_case_control_pred.list \
--bsize 400 \
--threads 15 \
--out ${outpath}/regenie_step2/UKB.clinical_case_control.condtion_snp

# 4. local ancestry analysis
# since this include additional mixed ancestry samples, we need to redo step1

# regenie step 1
regenie \
--step 1 \
--pgen ${pfile_array} \
--keep /path/to/UKB.clinical_case_control_mixed_ancestry.samples.txt \
--phenoFile /path/to/UKB.clinical_case_control_mixed_ancestry.pheno_covar.txt \
--phenoCol AD \
--covarFile /path/to/UKB.clinical_case_control_mixed_ancestry.pheno_covar.txt \
--covarColList age,sex,center{1:3},pc_{1:20} \
--bt \
--bsize 1000 \
--threads 80 \
--out ${outpath}/regenie_step1/UKB.clinical_case_control_mixed_ancestry

# regenie step 2 per LA regions
[[ ! -d ${outpath}/regenie_step2/local_ancestry ]] && mkdir -p ${outpath}/regenie_step2/local_ancestry

# loop over 8 windows
for i in {1..8}; do
    regenie \
    --step 2 \
    --pgen  ${outpath}/UKB.clinical_case_control_mixed_ancestry.repeat \
    --keep /path/to/UKB.clinical_case_control_mixed_ancestry.samples.txt \
    --phenoFile /path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar_la.txt \
    --phenoCol AD \
    --covarFile /path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar_la.txt \
    --covarColList age,sex,center{1:3},pc_{1:20},EUR_${i},EAS_${i},SAS_${i},AFR_${i},AMR_${i} \
    --interaction EUR_${i} \
    --bt \
    --pred ${outpath}/regenie_step1/UKB.clinical_case_control_mixed_ancestry_pred.list \
    --bsize 400 \
    --threads 15 \
    --out ${outpath}/regenie_step2/local_ancestry/UKB.clinical_case_control_mixed_ancestry.EUR_${i}

    regenie \
    --step 2 \
    --pgen  ${outpath}/UKB.clinical_case_control_mixed_ancestry.repeat \
    --keep /path/to/UKB.clinical_case_control_mixed_ancestry.samples.txt \
    --phenoFile /path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar_la.txt \
    --phenoCol AD \
    --covarFile /path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar_la.txt \
    --covarColList age,sex,center{1:3},pc_{1:20},EUR_${i},EAS_${i},SAS_${i},AFR_${i},AMR_${i} \
    --interaction EAS_${i} \
    --bt \
    --pred ${outpath}/regenie_step1/UKB.clinical_case_control_mixed_ancestry_pred.list \
    --bsize 400 \
    --threads 15 \
    --out ${outpath}/regenie_step2/local_ancestry/UKB.clinical_case_control_mixed_ancestry.EAS_${i}
done

# extract results
head -n 1 ${outpath}/regenie_step2/local_ancestry/UKB.clinical_case_control_mixed_ancestry.EAS_1_AD.regenie \
> ${outpath}/regenie_step2/local_ancestry/local_ancestry_summary.txt

cat ${outpath}/regenie_step2/local_ancestry/*regenie | grep 'FAM193B-CCG-ADD' | grep 'ADD-INT_SNPx' \
>> ${outpath}/regenie_step2/local_ancestry/local_ancestry_summary.txt

cat ${outpath}/regenie_step2/local_ancestry/*regenie | grep 'FAM193B-CCG-DOM' | grep 'ADD-INT_SNPx' \
>> ${outpath}/regenie_step2/local_ancestry/local_ancestry_summary.txt