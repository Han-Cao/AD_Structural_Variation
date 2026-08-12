#!/bin/bash

# locus-wise association analysis
# chr5:177322587-177754013 in GRCh38

pfile_hk=/path/to/hk.pfile
outpath=/path/to/association/

[[ ! -d ${outpath}/ ]] && mkdir -p ${outpath}/

# 1. AD associaition
# local assocation on AD
plink2 --pfile $pfile_hk \
--maf 0.05 \
--ci 0.95 \
--logistic hide-covar omit-ref \
--covar /path/to/pheno_covar.txt \
--covar-name Sex,Age,PC1,PC2,PC3 \
--pheno /path/to/pheno_covar.txt \
--pheno-name AD \
--no-input-missing-phenotype \
--freq \
--out ${outpath}/HK_association.AD

# condition on FAM193B-CCG
plink2 --pfile $pfile_hk \
--maf 0.05 \
--ci 0.95 \
--logistic hide-covar omit-ref \
--condition Chinese_SV.INS.chr5_4392 \
--covar /path/to/pheno_covar.txt \
--covar-name Sex,Age,PC1,PC2,PC3 \
--pheno /path/to/pheno_covar.txt \
--pheno-name AD \
--no-input-missing-phenotype \
--freq \
--out ${outpath}/HK_association.AD.CCG_condition

# condition on top SNP
plink2 --pfile $pfile_hk \
--maf 0.05 \
--ci 0.95 \
--logistic hide-covar omit-ref \
--condition 5:177464920_T_C \
--covar /path/to/pheno_covar.txt \
--covar-name Sex,Age,PC1,PC2,PC3 \
--pheno /path/to/pheno_covar.txt \
--pheno-name AD \
--no-input-missing-phenotype \
--freq \
--out ${outpath}/HK_association.AD.SNP_condition


# 2. Blood eQTL analysis

# FAM193B eQTL
plink2 --pfile $pfile_hk \
--linear hide-covar omit-ref \
--maf 0.05 \
--ci 0.95 \
--covar /path/to/eQTL.covar.txt \
--covar-variance-standardize \
--pheno /path/to/FAM193B.txt \
--pheno-name FAM193B_rank \
--no-input-missing-phenotype \
--out ${outpath}/HK_association.FAM193B_rank

# exclude AD samples for SMR analysis
plink2 --pfile $pfile_hk \
--remove /path/to/AD.txt \
--linear hide-covar omit-ref \
--maf 0.05 \
--ci 0.95 \
--covar /path/to/eQTL.covar.txt \
--covar-variance-standardize \
--pheno /path/to/FAM193B.txt \
--pheno-name FAM193B_rank \
--no-input-missing-phenotype \
--out ${outpath}/HK_association.FAM193B_rank.exclude_AD

