#!/bin/bash

# examine whether SNP + SV is due to LD by condition anlaysis on SVs

bim="/path/to/bfile/SRS.SV.autosome.MAF5.HWE1E-6.GENO10.bim"
inpath="/path/to/PRS/"
outpath="/path/to/PRS/condition/"
[[ ! -d $outpath ]] && mkdir -p $outpath

# prepare SV position file
awk 'OFS="\t"{print $2,$1,$4}' $bim | sort -k1,1 > "${outpath}/SRS.clean.AD.SV.pos"

ls input/fold_AD/SRS.analysis_sample.fold*.fam | while read fam;
do
    fold=$(basename $fam)
    fold=${fold%.fam}
    fold=${fold#SRS.analysis_sample.}

    # 1. extract lead SVs to adjust, join with their position
    clumped_vars="${inpath}/SNP_SV/SRS.clean.AD.${fold}.SNP_SV.clumps.id"
    grep 'Chinese_SV' $clumped_vars | sort | \
    join -11 -21 ${outpath}/SRS.clean.AD.SV.pos - | sort -k2,2n -k3,3n > "${outpath}/SRS.clean.AD.${fold}.lead_sv.txt"

    # 2. condition on lead SVs for SNPs within 1 MB
    mkdir -p "${outpath}/tmp"

    parallel -a ${outpath}/SRS.clean.AD.${fold}.lead_sv.txt --colsep " " -j 10 "

    # calcualte 1 MB window
    if [ {3} -lt 500000 ]; then
        start=0
    else
        start=\$(({3} - 500000))
    fi
    end=\$(({3} + 500000))

    # conditional analysis
    plink2 --bfile /path/to/bfile/SRS.SNP.SV.autosome.MAF5.HWE1E-6.GENO1 \
    --allow-no-sex \
    --keep-allele-order \
    --remove $fam \
    --logistic hide-covar omit-ref no-firth \
    --covar input/SRS.pheno_covar.txt \
    --covar-name Sex,Age,PC1,PC2,PC3 \
    --pheno input/SRS.pheno_covar.txt \
    --pheno-name AD \
    --chr {2} \
    --from-bp \${start} \
    --to-bp \${end} \
    --condition {1} \
    --out ${outpath}/tmp/SRS.clean.${fold}.SNP.condition_{1}
    "

    # merge all log files
    cat ${outpath}/tmp/*.log > ${outpath}/SRS.clean.${fold}.SNP.condition.log

    # merge all association files
    # prepare header
    echo -e "ID\tA1\tBETA\tP" > ${outpath}/SRS.clean.${fold}.SNP.condition.glm.logistic

    # sort SNPs by P value, remove SNPs without a P value, and keep duplicated SNP with largest P value
    cat ${outpath}/tmp/*.glm.logistic | cut -f 3,7,12,15 | grep -v 'Chinese_SV' | grep -v 'OR' | sort -k15,15gr | \
    awk '$4 != "NA" ' - | awk '!seen[$1]++ {print $1, $2, log($3), $4}' \
    >> ${outpath}/SRS.clean.${fold}.SNP.condition.glm.logistic

    # merge condition results with original summ stat
    # for duplicated SNPs, use the one in condition results
    result_snp_raw="${inpath}/discovery/SRS.clean.${fold}.SNP.AD.glm.logistic"

    cat ${outpath}/SRS.clean.${fold}.SNP.condition.glm.logistic $result_snp_raw | \
    awk '!seen[$1]++' > ${outpath}/SRS.clean.${fold}.SNP.merge.glm.logistic

    # 3. clump on merged P values
    # since we didn't pre-filter the P values in the merged GWAS results, we need to filter it by --clump-p1 0.01 (best cutoff) here
    plink2 --bfile /path/to/bfile/SRS.SNP.SV.autosome.MAF5.HWE1E-6.GENO1 \
    --clump-p1 0.01 \
    --clump-r2 0.1 \
    --clump-kb 500 \
    --clump ${outpath}/SRS.clean.${fold}.SNP.merge.glm.logistic \
    --clump-snp-field ID \
    --clump-field P \
    --out ${outpath}/SRS.clean.${fold}.SNP.merge
    
    # extract clumped vairnats ID
    awk "NR!=1{print \$3}" ${outpath}/SRS.clean.${fold}.SNP.merge.clumps \
    > ${outpath}/SRS.clean.${fold}.SNP.merge.clumps.id

    # 4. calculate PRS for SNPadj
    plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
    --keep-allele-order \
    --score ${outpath}/SRS.clean.${fold}.SNP.merge.glm.logistic 1 2 3 header center \
    --keep $fam \
    --extract ${outpath}/SRS.clean.${fold}.SNP.merge.clumps.id \
    --out ${outpath}/SRS.clean.${fold}.SNP.merge.score

    # 5. clean up
    rm -r "${outpath}/tmp"
done


# apply the representative model to all samples
fold_rep=fold62

plink2 --bfile "/path/to/bfile/SRS.SNP.autosome.MAF5.HWE1E-6.GENO1" \
--keep-allele-order \
--score ${outpath}/SRS.clean.${fold_rep}.SNP.merge.glm.logistic 1 2 3 header center \
--remove input/fold_AD_62.train.fam \
--extract ${outpath}/SRS.clean.${fold_rep}.SNP.merge.clumps.id \
--out "/path/to/PRS/SNP_SV/represent/SRS.clean.AD.${fold_rep}.SNP_condition.score"