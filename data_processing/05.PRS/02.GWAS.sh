#!/bin/bash

bfile_sv=/path/to/bfile/SRS.SV.autosome.MAF1.HWE1E-6.GENO10
bfile_snp=/path/to/bfile/SRS.SNP.autosome.MAF1.HWE1E-6.GENO1
outpath=/path/to/GWAS/association/

[[ ! -d $outpath ]] && mkdir -p $outpath

# SV GWAS on AD
plink2 --bfile $bfile_sv \
--maf 0.05 \
--exclude input/qc_fail_sv_bqsr.txt \
--allow-no-sex \
--keep-allele-order \
--logistic hide-covar omit-ref \
--ci 0.95 \
--covar input/SRS.pheno_covar.txt \
--covar-name Sex,Age,PC1,PC2,PC3 \
--pheno input/SRS.pheno_covar.txt \
--pheno-name AD \
--freq \
--out ${outpath}/SRS.SV.AD

# SNP GWAS on AD
plink2 --bfile $bfile_snp \
--maf 0.05 \
--allow-no-sex \
--keep-allele-order \
--logistic hide-covar omit-ref \
--ci 0.95 \
--covar input/SRS.pheno_covar.txt \
--covar-name Sex,Age,PC1,PC2,PC3 \
--pheno input/SRS.pheno_covar.txt \
--pheno-name AD \
--freq \
--out ${outpath}/SRS.SNP.AD
