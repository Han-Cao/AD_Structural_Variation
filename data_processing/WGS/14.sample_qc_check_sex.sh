#!/bin/bash
#SBATCH -p cpu -N 1 -n 1 -o log/sample_QC/check_sex.log

bfile="/path/to/HK_WGS_update.PASS.raw"
outpath="/path/to/check_sex/"

[[ ! -d $outpath ]] && mkdir -p $outpath

#split X chromosome (hg38)
plink --bfile $bfile \
--chr 23-25 \
--make-bed \
--split-x 2781479 155701383 \
--maf 0.05 \
--out "${outpath}/HK_WGS_update.PASS.chrX.MAF5"

# pruning
plink --bfile "${outpath}/HK_WGS_update.PASS.chrX.MAF5" \
--indep-pairphase 20000 2000 0.5 \
--out "${outpath}/HK_WGS_update.PASS.chrX.MAF5.prune"

#check sex
plink --bfile "${outpath}/HK_WGS_update.PASS.chrX.MAF5" \
--extract "${outpath}/HK_WGS_update.PASS.chrX.MAF5.prune.prune.in" \
--check-sex \
--out "${outpath}/HK_WGS_update.PASS.check_sex"