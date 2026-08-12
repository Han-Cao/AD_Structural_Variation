#!/bin/bash

vcf_snp=/path/to/clean/SRS.SNP.vcf.gz
vcf_sv=/path/to/clean/SRS.SV.vcf.gz

outpath="/path/to/bfile/"

[[ ! -d $outpath ]] && mkdir -p $outpath

# extrat analysis samples
plink2 --vcf $vcf_snp \
--keep input/qc_pass_samples.txt \
--mac 1 \
--make-bed \
--set-all-var-ids "@:#_\$r_\$a" \
--out ${outpath}/SRS.SNP.autosome.raw

plink2 --vcf $vcf_sv \
--keep input/qc_pass_samples.txt \
--mac 1 \
--make-bed \
--autosome \
--out ${outpath}/SRS.SV.autosome.raw

# filter by maf, geno, hwe
plink2 --bfile ${outpath}/SRS.SNP.autosome.raw \
--maf 0.01 \
--hwe 1e-6 \
--geno 0.01 \
--make-bed \
--out ${outpath}/SRS.SNP.autosome.MAF1.HWE1E-6.GENO1

plink2 --bfile ${outpath}/SRS.SV.autosome.raw \
--hwe 1e-6 \
--geno 0.1 \
--make-bed \
--out ${outpath}/SRS.SV.autosome.MAF1.HWE1E-6.GENO10

# for PRS, we additionally apply --maf 0.05
plink2 --bfile ${outpath}/SRS.SNP.autosome.MAF1.HWE1E-6.GENO1 \
--maf 0.05 \
--make-bed \
--out ${outpath}/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1

plink2 --bfile ${outpath}/SRS.SV.autosome.MAF1.HWE1E-6.GENO10 \
--maf 0.05 \
--make-bed \
--out ${outpath}/SRS.SV.autosome.MAF5.HWE1E-6.GENO10

# merge SNP and SV
plink --bfile ${outpath}/SRS.SNP.autosome.MAF1.HWE1E-6.GENO1 \
--keep-allele-order \
--bmerge ${outpath}/SRS.SV.autosome.MAF1.HWE1E-6.GENO10 \
--make-bed \
--out ${outpath}/SRS.SNP.SV.autosome.MAF1.HWE1E-6.GENO1

plink --bfile ${outpath}/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1 \
--keep-allele-order \
--bmerge ${outpath}/SRS.SV.autosome.MAF5.HWE1E-6.GENO10 \
--make-bed \
--out ${outpath}/SRS.SNP.SV.autosome.MAF5.HWE1E-6.GENO1