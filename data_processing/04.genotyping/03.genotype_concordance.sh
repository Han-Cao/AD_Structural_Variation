#!/bin/bash

truth_vcf="/path/to/All.genotype.filter.population.AC0.clean.vcf.gz"
outpath=/path/to/genotyping/bench/

harmonisv concordance \
-i "/path/to/SRS.overlap_samples.vcf.gz" \
-o "${outpath}/paragraph_concordance.35X.vcf.gz" \
-r $truth_vcf

# extract concordant SVs
bcftools query -i 'CON_GT=1' -f '%ID\n' "${outpath}/paragraph_concordance.35X.vcf.gz" > ${outpath}/concordance_sv.txt
