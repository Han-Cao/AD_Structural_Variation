#!/bin/bash

outpath="/path/to/PRS/discovery/"
[[ ! -d $outpath ]] && mkdir -p $outpath

# re-run GWAS using training samples only for each fold
ls input/fold_AD/SRS.analysis_sample.fold*.fam | while read fam;
do
   fold=$(basename $fam)
   fold=${fold%.fam}
   fold=${fold#SRS.analysis_sample.}

   ############# SV analysis #############
   # GWAS
   plink2 --bfile "/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10" \
   --allow-no-sex \
   --keep-allele-order \
   --remove $fam \
   --logistic hide-covar omit-ref no-firth \
   --covar input/SRS.pheno_covar.txt \
   --covar-name Sex,Age,PC1,PC2,PC3 \
   --pheno input/SRS.pheno_covar.txt \
   --pheno-name AD \
   --out ${outpath}/SRS.clean.${fold}.SV
   
   # only keep ID, A1, beta and P for PRS
   echo -e "ID\tA1\tBETA\tP" > ${outpath}/SRS.clean.${fold}.SV.AD.glm.logistic.tmp
   awk "FNR>1{if(\$15<=0.1) print \$3, \$7, log(\$12), \$15}" ${outpath}/SRS.clean.${fold}.SV.AD.glm.logistic \
   >> ${outpath}/SRS.clean.${fold}.SV.AD.glm.logistic.tmp
   mv ${outpath}/SRS.clean.${fold}.SV.AD.glm.logistic.tmp ${outpath}/SRS.clean.${fold}.SV.AD.glm.logistic


   ############# SNP analysis #############
   # GWAS
   plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
   --allow-no-sex \
   --keep-allele-order \
   --remove $fam \
   --logistic hide-covar omit-ref no-firth \
   --covar input/SRS.pheno_covar.txt \
   --covar-name Sex,Age,PC1,PC2,PC3 \
   --pheno input/SRS.pheno_covar.txt \
   --pheno-name AD \
   --out ${outpath}/SRS.clean.${fold}.SNP

   # only keep ID, A1, beta and P for PRS
   echo -e "ID\tA1\tBETA\tP" > ${outpath}/SRS.clean.${fold}.SNP.AD.glm.logistic.tmp
   awk "FNR>1{if(\$15<=0.1) print \$3, \$7, log(\$12), \$15}" ${outpath}/SRS.clean.${fold}.SNP.AD.glm.logistic \
   >> ${outpath}/SRS.clean.${fold}.SNP.AD.glm.logistic.tmp
   mv ${outpath}/SRS.clean.${fold}.SNP.AD.glm.logistic.tmp ${outpath}/SRS.clean.${fold}.SNP.AD.glm.logistic
done