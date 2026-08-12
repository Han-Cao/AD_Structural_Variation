#!/bin/bash

outpath="/path/to/PRS/score/"
[[ ! -d $outpath ]] && mkdir -p $outpath

ls input/fold_AD/SRS.analysis_sample.fold*.fam | while read fam;
do
   fold=$(basename $fam)
   fold=${fold%.fam}
   fold=${fold#SRS.analysis_sample.}

   ############# SV analysis #############

   result_sv=/path/to/PRS/discovery/SRS.clean.${fold}.SV.AD.glm.logistic

   plink2 --bfile "/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10" \
   --keep-allele-order \
   --score $result_sv 1 2 3 header center \
   --q-score-range input/q_score_range.txt $result_sv 1 4 \
   --keep $fam \
   --extract "/path/to/PRS/clump/SRS.clean.AD.${fold}.SV.clumps.id" \
   --out "${outpath}/SRS.clean.AD.${fold}.SV"


   ############# SNP analysis #############

   result_snp=/path/to/PRS/discovery/SRS.clean.${fold}.SNP.AD.glm.logistic
   # rs429358 is the only SNP within the APOE region after clumping across all folds, so we only need to exclude it
   plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
   --exclude-snp 19:44908684_T_C \
   --window 500 \
   --keep-allele-order \
   --score $result_snp 1 2 3 header center \
   --q-score-range input/q_score_range.txt $result_snp 1 4 \
   --keep $fam \
   --extract "/path/to/PRS/clump/SRS.clean.AD.${fold}.SNP.clumps.id" \
   --out "${outpath}/SRS.clean.AD.${fold}.SNP"
done
