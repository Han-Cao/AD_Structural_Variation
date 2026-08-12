#!/bin/bash

invcf="/path/to/Joint_call/genotype_sample/All.genotype.vcf.gz"
outpath="/path/to/Joint_call/filterSV/"

# model parameters are determined based on 
# the built-in k-fold cross validation in filterSV.py
# i.e., run with default paramters

# filter insertion
python /path/to/code/harmoniSV/filterSV.py \
-i $invcf \
-o ${outpath}/All.genotype.RF_INS \
--sv-type INS \
--feature SVLEN,MEAN_AF,STD_AF,DP_MINIMAP2_CUTESV,AF_MINIMAP2_CUTESV,AF_MINIMAP2_SVIM,DP_NGMLR_CUTESV,AF_NGMLR_CUTESV,AF_NGMLR_SVIM \
--min-support-method 6 \
--merge /path/to/Joint_call/merge_with_ref/Chinese_SV.merge_with_ref.txt \
--train-size 1 \
--target-prefix Chinese_SV \
--min-support-set 4 \
--min-re 4 \
--max-features 4 \
--max-depth 12 \
--min-samples-leaf 1 \
--bench-vcf /path/to/validation/SV_call/genotype_sample/HG002.genotype.vcf.gz \
--bench-sites /path/to/validation/bench/validation_set.txt \
--model-thread 20

# filter deletion
python /path/to/code/harmoniSV/filterSV.py \
-i ${outpath}/All.genotype.RF_INS.vcf.gz \
-o ${outpath}/All.genotype.RF_INS_DEL \
--sv-type DEL \
--feature SVLEN,MEAN_AF,STD_AF,DP_MINIMAP2_SNIFFLES,AF_MINIMAP2_SNIFFLES,AF_NGMLR_SVIM,DP_NGMLR_CUTESV,AF_NGMLR_CUTESV \
--min-support-method 6 \
--merge /path/to/Joint_call/merge_with_ref/Chinese_SV.merge_with_ref.txt \
--train-size 1 \
--target-prefix Chinese_SV \
--min-support-set 4 \
--min-re 4 \
--max-features 5 \
--max-depth 12 \
--min-samples-leaf 5 \
--bench-vcf /path/to/validation/SV_call/genotype_sample/HG002.genotype.vcf.gz \
--bench-sites /path/to/validation/bench/validation_set.txt \
--model-thread 20


# apply filter
# INS/DEL: use RF model by cutoff of precision = 0.95
# Only a small number of DUP/INV are identified, so we used hard filter
# DUP: hard filter by consensus call from Sniffles and cuteSV (SVIM don't call DUP)
# INV: hard filter by consensus call from Sniffles, cuteSV, and SVIM
bcftools filter -e '(SVTYPE=="DUP" & SUPP_METHOD<4) | (SVTYPE=="INV" & SUPP_METHOD<6)' \
-s CONSENSUS -Ou ${outpath}/All.genotype.RF_INS_DEL.vcf.gz | \
bcftools filter -e '(SVTYPE=="INS" & RF_SCORE < 0.779)' -s RF -Ou | \
bcftools filter -e '(SVTYPE=="DEL" & RF_SCORE < 0.135)' -s RF -Oz \
-o ${outpath}/All.genotype.filter.vcf.gz
tabix ${outpath}/All.genotype.filter.vcf.gz


# concat bench VCF
bcftools view -i 'SVTYPE="INS"' -Oz -o ${outpath}HG002.genotype.RF_INS.vcf.gz \
${outpath}All.genotype.RF_INS.bench.vcf.gz
tabix ${outpath}HG002.genotype.RF_INS.vcf.gz

bcftools view -i 'SVTYPE="DEL"' -Oz -o ${outpath}HG002.genotype.RF_DEL.vcf.gz \
${outpath}All.genotype.RF_INS_DEL.bench.vcf.gz
tabix ${outpath}HG002.genotype.RF_DEL.vcf.gz

bcftools concat -a -Oz -o ${outpath}HG002.genotype.RF_INS_DEL.vcf.gz \
${outpath}HG002.genotype.RF_INS.vcf.gz \
${outpath}HG002.genotype.RF_DEL.vcf.gz
tabix ${outpath}HG002.genotype.RF_INS_DEL.vcf.gz