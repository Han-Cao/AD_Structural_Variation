#!/bin/bash

# GATK 4.4.0.0

refPre="/path/to/GRCh38.fa"

name=$1
interval_list=$2
outpath=$3
sample_map=$4
finish=$5
error=$6

threads=10

# run GenomicsDBImport
if [[ ! -d ${outpath}/${name}/ ]]; then
    # first run, create new GenomicsDB
    gatk --java-options "-Xmx40G" GenomicsDBImport \
    --genomicsdb-workspace-path "${outpath}/${name}" \
    -L $interval_list \
    --interval-padding 1000 \
    --sample-name-map $sample_map \
    -R $refPre \
    --genomicsdb-shared-posixfs-optimizations \
    --batch-size 50 \
    --reader-threads $threads
else
    # not first run, append to existing GenomicsDB
    # -L is not needed as it is already specified in the existing GenomicsDB
    gatk --java-options "-Xmx40G" GenomicsDBImport \
    --genomicsdb-update-workspace-path "${outpath}/${name}" \
    --sample-name-map $sample_map \
    -R $refPre \
    --genomicsdb-shared-posixfs-optimizations \
    --batch-size 50 \
    --reader-threads $threads
fi


# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${name}.error" ]] && mv "${error}/${name}.error" "${error}/${name}.error_fixed"
    touch "${finish}/${name}.done"
else
    touch "${error}/${name}.error"
fi