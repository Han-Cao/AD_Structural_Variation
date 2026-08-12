#!/bin/bash
#SBATCH -p cpu -N 1 -n 20 -o log/VQSR/VQSR_SNP.log

refPre="/path/to/GRCh38.fa"
dbsnp="/path/to/dbsnp_153.GRCh38.vcf.gz"
omni="/path/to/1000G_omni2.5.hg38.PASS.vcf.gz"
hapmap="/path/to/hapmap_3.3.hg38.vcf.gz"
g1000="/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

vcf="/path/to/Joint-call/HK_WGS_update.hardfilter.vcf.gz"
outpath="/path/to/VQSR/"

[[ ! -d $outpath ]] && mkdir $outpath

name=$(basename $vcf)
name=${name%.vcf.gz}

# GATK 4.4.0.0
# Remove MQ annotation if low sd
gatk --java-options "-Xmx120G" VariantRecalibrator \
-R $refPre \
-V $vcf \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
--resource:omni,known=false,training=true,truth=true,prior=12.0 $omni \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 $g1000 \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
-AS \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
--max-gaussians 8 \
-mode SNP \
-tranche 100.0 -tranche 99.95 -tranche 99.90 -tranche 99.80 -tranche 99.60 \
-tranche 99.50 -tranche 99.40 -tranche 99.30 -tranche 99.00 -tranche 98.00 \
-tranche 97.00 -tranche 96.00 -tranche 95.00 -tranche 94.00 -tranche 93.00 \
-tranche 92.00 -tranche 91.00 -tranche 90.00 -tranche 85.00 \
-O "${outpath}${name}.VQSR_SNP.recal" \
--tranches-file "${outpath}${name}.VQSR_SNP.tranches" \
--rscript-file "${outpath}${name}.VQSR_SNP.R"