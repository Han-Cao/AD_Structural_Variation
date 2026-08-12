#!/bin/bash
#SBATCH -p cpu -N 1 -n 20 -o log/VQSR/VQSR_INDEL.log

refPre="/path/to/GRCh38.fa"
dbsnp="/path/to/dbsnp_153.GRCh38.vcf.gz"
omni="/path/to/1000G_omni2.5.hg38.PASS.vcf.gz"
hapmap="/path/to/hapmap_3.3.hg38.vcf.gz"
g1000="/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

vcf="path/to/VQSR/HK_WGS_update.hardfilter.VQSR_SNP.vcf.gz"
outpath="path/to/VQSR/"

name=$(basename $vcf)
name=${name%.VQSR_SNP.vcf.gz}

[[ ! -d $outpath ]] && mkdir $outpath

gatk --java-options "-Xmx120G" VariantRecalibrator \
-R $refPre \
-V $vcf \
--resource:mills,known=true,training=true,truth=true,prior=12.0 $mills \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
-AS \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--max-gaussians 4 \
-mode INDEL \
-tranche 100.0 -tranche 99.95 -tranche 99.90 -tranche 99.80 -tranche 99.60 \
-tranche 99.50 -tranche 99.40 -tranche 99.30 -tranche 99.00 -tranche 98.00 \
-tranche 97.00 -tranche 96.00 -tranche 95.00 -tranche 94.00 -tranche 93.00 \
-tranche 92.00 -tranche 91.00 -tranche 90.00 -tranche 85.00 \
-O "${outpath}${name}.VQSR_INDEL.recal" \
--tranches-file "${outpath}${name}.VQSR_INDEL.tranches" \
--rscript-file "${outpath}${name}.VQSR_INDEL.R"