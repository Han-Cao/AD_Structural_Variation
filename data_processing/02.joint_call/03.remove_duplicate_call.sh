#!/bin/bash

invcf_path="/path/to/SV_call/harmonized/"
jasmine_path="/path/to/SV_call/dup_call/"
outpath="/path/to/Joint_call/preprocess/"

[[ ! -d $outpath ]] && mkdir -p $outpath

ls ${invcf_path}/*harmonized.vcf.gz | while read invcf
do
    name=$(basename $invcf)
    name=${name%.vcf.gz}

    invcf_info="${invcf_path}/table/${name}.txt.gz"
    outvcf="${outpath}/${name}.dedup.vcf"
    jasmine_vcf="${jasmine_path}/${name}.preprocess.vcf"
    jasmine_merge_list=${jasmine_vcf%.vcf}.merge_list.txt

    # select the representative SV among duplicated calls by maximum supporting reads
    python /path/to/code/harmoniSV/representSV.py \
    -i $invcf \
    -o $outvcf \
    --merge $jasmine_merge_list \
    --sv-info $invcf_info \
    --by-max RE
done

