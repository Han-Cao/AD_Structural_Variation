#!/bin/bash

invcf_path="/path/to/Joint_call/preprocess/"
outpath="/path/to/Joint_call/merge_sv/"

[[ ! -d $outpath ]] && mkdir -p $outpath

ls ${invcf_path}/*vcf > ${outpath}/SV_merge_vcf_list.txt

# merge across samples
jasmine \
file_list=${outpath}/SV_merge_vcf_list.txt \
out_file=${outpath}/All_sample_pipeline_merged.vcf \
threads=10 \
--keep_var_ids \
--ignore_strand 

# harmonize header
bash /path/to/code/harmoniSV/harmonize_vcf_header.sh \
-f ${outpath}/SV_merge_vcf_list.txt \
-o ${outpath}/All_sample_pipeline_merged.header \
-v ${outpath}/All_sample_pipeline_merged.vcf

bcftools reheader -h ${outpath}/All_sample_pipeline_merged.header \
-o ${outpath}/All_sample_pipeline_merged.reheader.vcf \
${outpath}/All_sample_pipeline_merged.vcf

bcftools sort -Oz -o ${outpath}/All_sample_pipeline_merged.vcf.gz ${outpath}/All_sample_pipeline_merged.reheader.vcf
rm ${outpath}/All_sample_pipeline_merged.vcf ${outpath}/All_sample_pipeline_merged.reheader.vcf

# write merged SV list
bcftools query -f '%IDLIST\n' ${outpath}/All_sample_pipeline_merged.vcf.gz \
> ${outpath}/All_sample_pipeline_merged.txt