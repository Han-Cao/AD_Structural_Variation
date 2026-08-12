#!/bin/bash

outpath="/scratch/PI/boip/Han/project/SV/PRS/clump/"
[[ ! -d $outpath ]] && mkdir -p $outpath

# LD clumping to select SNPs
ls input/fold_AD/SRS.analysis_sample.fold*.fam | while read fam;
do
   fold=$(basename $fam)
   fold=${fold%.fam}
   fold=${fold#SRS.analysis_sample.}

   ############# SV analysis #############
   # clumping
   plink2 --bfile "/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10" \
   --clump-p1 0.1 \
   --clump-r2 0.1 \
   --clump-kb 500 \
   --clump "/path/to/PRS/discovery/SRS.clean.${fold}.SV.AD.glm.logistic" \
   --clump-snp-field ID \
   --clump-field P \
   --out "${outpath}/SRS.clean.AD.${fold}.SV"
   
   # extract clumped SVs
   awk "NR!=1{print \$3}" "${outpath}/SRS.clean.AD.${fold}.SV.clumps" \
   > "${outpath}/SRS.clean.AD.${fold}.SV.clumps.id"


   ############# SNP analysis #############
   # clumping
   plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
   --clump-p1 0.1 \
   --clump-r2 0.1 \
   --clump-kb 500 \
   --clump "/path/to/PRS/discovery/SRS.clean.${fold}.SNP.AD.glm.logistic" \
   --clump-snp-field ID \
   --clump-field P \
   --out "${outpath}/SRS.clean.AD.${fold}.SNP"

   # extract clumped SNPs
   awk "NR!=1{print \$3}" "${outpath}/SRS.clean.AD.${fold}.SNP.clumps" \
   > "${outpath}/SRS.clean.AD.${fold}.SNP.clumps.id"
done