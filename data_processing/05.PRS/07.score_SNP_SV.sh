#!/bin/bash

outpath="/path/to/PRS/SNP_SV/"
[[ ! -d $outpath ]] && mkdir -p $outpath

# this is the best cutoffs for SNP and SV
snp_cutoff=0.01
sv_cutoff=0.05

# clump 10 fold samples
ls input/fold_AD/SRS.analysis_sample.fold*.fam | while read fam;
do
   fold=$(basename $fam)
   fold=${fold%.fam}
   fold=${fold#SRS.analysis_sample.}

   clumps_snp="/path/to/PRS/clump/SRS.clean.AD.${fold}.SNP.clumps"
   clumps_sv="/path/to/PRS/clump/SRS.clean.AD.${fold}.SV.clumps"

   beta_sv=/path/to/PRS/discovery/SRS.clean.${fold}.SV.AD.glm.logistic
   beta_snp=/path/to/PRS/discovery/SRS.clean.${fold}.SNP.AD.glm.logistic

   ############# clump SV + SNP #############
   # prepare P value-only input for clumping
   echo "ID P" > "${outpath}/SRS.clean.AD.${fold}.SNP_SV.input.txt"
   awk "{if(\$4 < $sv_cutoff) print \$3,\$4}" $clumps_sv >> "${outpath}/SRS.clean.AD.${fold}.SNP_SV.input.txt"
   awk "{if(\$4 < $snp_cutoff) print \$3,\$4}" $clumps_snp >> "${outpath}/SRS.clean.AD.${fold}.SNP_SV.input.txt"

   plink2 --bfile "/path/to/bfile/SRS.SNP.SV.autosome.MAF5.HWE1E-6.GENO1" \
   --clump-p1 1 \
   --clump-r2 0.1 \
   --clump-kb 500 \
   --clump "${outpath}/SRS.clean.AD.${fold}.SNP_SV.input.txt" \
   --clump-snp-field ID \
   --clump-field P \
   --out "${outpath}/SRS.clean.AD.${fold}.SNP_SV"
   
   # extract clumped variants
   awk "NR!=1{print \$3}" "${outpath}/SRS.clean.AD.${fold}.SNP_SV.clumps" \
   > "${outpath}/SRS.clean.AD.${fold}.SNP_SV.clumps.id"


   ############# PRS #############
   # calcualte PRS after clumping for SNP+SV model
   plink2 --bfile "/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10" \
   --keep-allele-order \
   --score $beta_sv 1 2 3 header center \
   --keep $fam \
   --extract "${outpath}/SRS.clean.AD.${fold}.SNP_SV.clumps.id" \
   --out "${outpath}/SRS.clean.AD.${fold}.SV.score"

   plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
   --keep-allele-order \
   --score $beta_snp 1 2 3 header center \
   --keep $fam \
   --extract "${outpath}/SRS.clean.AD.${fold}.SNP_SV.clumps.id" \
   --out "${outpath}/SRS.clean.AD.${fold}.SNP.score"

done

# apply the representative model to all samples
fold_rep=fold62

beta_sv=/path/to/PRS/discovery/SRS.clean.${fold_rep}.SV.AD.glm.logistic
beta_snp=/path/to/PRS/discovery/SRS.clean.${fold_rep}.SNP.AD.glm.logistic

plink2 --bfile "/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10" \
--keep-allele-order \
--score $beta_sv 1 2 3 header center \
--remove /path/to/code/PRS/input/fold_AD_62.train.fam \
--extract "${outpath}/SRS.clean.AD.${fold_rep}.SNP_SV.clumps.id" \
--out "${outpath}/../represent/SRS.clean.AD.${fold_rep}.SV.score"

plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
--keep-allele-order \
--score $beta_snp 1 2 3 header center \
--remove /path/to/code/PRS/input/fold_AD_62.train.fam \
--extract "${outpath}/SRS.clean.AD.${fold_rep}.SNP_SV.clumps.id" \
--out "${outpath}/../represent/SRS.clean.AD.${fold_rep}.SNP.score"