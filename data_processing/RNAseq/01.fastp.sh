#!/bin/bash

set -euxo pipefail

# fastp 0.23.4

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

[[ ! -d "${outpath}/${sample}" ]] && mkdir -p "${outpath}/${sample}"

# run fastp for all fastq files
ls ${inpath}${sample}/*_1.fq.gz | while read fastq
do
        fastq=${fastq%_1.fq.gz}
        fq1=${fastq}_1.fq.gz
        fq2=${fastq}_2.fq.gz
        name=$(basename $fastq)

        fastp --thread 10 \
        --in1 $fq1 --in2 $fq2 \
        --out1 "${outpath}/${sample}/${name}_1_paired.fq.gz" \
        --out2 "${outpath}/${sample}/${name}_2_paired.fq.gz" \
        --unpaired1 "${outpath}/${sample}/${name}_1_unpaired.fq.gz" \
        --unpaired2 "${outpath}/${sample}/${name}_2_unpaired.fq.gz" \
        --length_required 40 \
        --verbose \
        --json "${outpath}/${sample}/${name}_fastp.json" \
        --html "${outpath}/${sample}/${name}_fastp.html"
done

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi