#!/bin/bash

# prepare input VCFs
bcftools view -o /path/to/Joint_call/merge_with_ref/Chinese_SV.vcf \
/path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.all.vcf.gz

echo '/path/to/Joint_call/merge_with_ref/Chinese_SV.vcf
/path/to/reference/jasmine_input/1KGP.vcf
/path/to/reference/jasmine_input/Chinese_405.vcf
/path/to/reference/jasmine_input/gnomAD.vcf
/path/to/reference/jasmine_input/HGDP.vcf
/path/to/reference/jasmine_input/HGSVC.vcf
/path/to/reference/jasmine_input/Icelander.vcf
/path/to/reference/jasmine_input/LRS15.vcf' > \
/path/to/Joint_call/merge_with_ref/merge_vcf.txt

# SV merging
jasmine \
file_list=/path/to/Joint_call/merge_with_ref/merge_vcf.txt \
out_file=/path/to/Joint_call/merge_with_ref/Chinese_SV.merge_with_ref.vcf \
threads=10 \
--keep_var_ids \
--ignore_strand 

# extract SV merging results
bcftools query -f '%IDLIST\n' /path/to/Joint_call/merge_with_ref/Chinese_SV.merge_with_ref.vcf | \
grep 'Chinese_SV' > /path/to/Joint_call/merge_with_ref/Chinese_SV.merge_with_ref.txt