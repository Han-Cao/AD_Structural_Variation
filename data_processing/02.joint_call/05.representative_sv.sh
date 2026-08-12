#!/bin/bash

inpath="/path/to/Joint_call/merge_sv/"

# identify representative SV by most frequent POS and SVLEN
parallel -j 20 '
python /path/to/code/harmoniSV/representSV.py \
-i /path/to/SV_call/harmonized/All_samples_pipelines.vcf.gz \
-o /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.chr{1}.vcf.gz \
-r chr{1} \
--merge /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.txt \
--sv-info /path/to/SV_call/harmonized/table/All_samples_pipelines.txt.gz \
--by-freq \
--id-prefix Chinese_SV \
--save-id \
--remove-hom-ref \
--min-len 50
' ::: {1..22} X

# concat all chr results
for i in {1..22} X; do
    echo /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.chr${i}.vcf.gz
done > /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.vcf_list.txt

bcftools concat \
-f /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.vcf_list.txt \
-n \
-Oz \
-o /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.all.vcf.gz
tabix /path/to/Joint_call/merge_sv/All_sample_pipeline_merged.representative.all.vcf.gz