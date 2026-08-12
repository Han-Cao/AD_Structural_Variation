#!/bin/bash

# INFO merging criteria
# Sum of INFO: depth, supporting reads
# Average of INFO: number of supporting methods
# Min of INFO: max number of supporting reads among all methods
# Max of INFO: RF score

# Merge VCF
python /path/to/code/harmoniSV/sample2pop.py \
-i /path/to/Joint_call/filterSV/All.genotype.filter.vcf.gz \
-o /path/to/Joint_call/population/All.genotype.filter.population.vcf.gz \
--filter-GT \
--info-first SVTYPE,SVLEN,END \
--info-sum DP*,RE* \
--info-avg SUPP_* \
--info-min MAX_RE \
--info-max RF_SCORE \
--info-to-format RF_SCORE \
--keep-format GT

tabix /path/to/Joint_call/population/All.genotype.filter.population.vcf.gz

# Remove SVs with AC=0
bcftools view -i 'AC>0' \
-Oz -o /path/to/Joint_call/population/All.genotype.filter.population.AC0.vcf.gz \
/path/to/Joint_call/population/All.genotype.filter.population.vcf.gz
tabix /path/to/Joint_call/population/All.genotype.filter.population.AC0.vcf.gz
