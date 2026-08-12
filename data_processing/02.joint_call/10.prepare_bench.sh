#!/bin/bash

invcf_path="/path/to/validation/bench/input/"
output_path="/path/to/validation/bench/output/"
sv_types="INS DEL"

# In this script, all input VCFs are generated 
# using the same code in 01.SV_call and 02.joint_call
# for the benchmark sample HG002

# prepare per-SVTYPE vcfs
for type in $sv_types; do
    # truth
    bcftools view -i "SVTYPE='$type'" \
    -Oz \
    -o /path/to/validation/SV_truth/HG002_SVs_Tier1_v0.6.PASS.${type}.vcf.gz \
    /path/to/validation/SV_truth/HG002_SVs_Tier1_v0.6.PASS.vcf.gz
    tabix /path/to/validation/SV_truth/HG002_SVs_Tier1_v0.6.PASS.${type}.vcf.gz

    # single method
    ls /path/to/validation/SV_call/Joint_call/single_method_clean/*vcf.gz | while read vcf;
    do
        name=$(basename $vcf | cut -f 1 -d '.')
        bcftools view -i "SVTYPE='$type'" -Oz -o ${invcf_path}/${name}.${type}.vcf.gz $vcf
        tabix ${invcf_path}/${name}.${type}.vcf.gz
    done

    # ensemble call
    bcftools view -i "SVTYPE='$type'" -Oz -o ${invcf_path}/ensemble.${type}.vcf.gz \
    /path/to/Joint_call/filterSV/HG002.genotype.RF_INS_DEL.vcf.gz
    tabix ${invcf_path}/ensemble.${type}.vcf.gz
done
