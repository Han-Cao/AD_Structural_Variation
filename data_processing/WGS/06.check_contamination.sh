#!/bin/bash

# VerifyBamID 2.0.1

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

svd="/path/to/ref/VerifyBamID/1000g.phase3.100k.b38.vcf.gz.dat"
ref="/path/to/ref/Human_genome/GRCh38/GRCh38.fa"


# check contamination
# set wd to output path, the output will always to wd
[[ ! -d ${outpath}/${sample} ]] && mkdir -p ${outpath}/${sample}
cd ${outpath}/${sample}

verifybamid2 \
--SVDPrefix $svd \
--Reference $ref \
--BamFile ${inpath}/${sample}.markdup.sort.cram

# mark finish or error
if [[ $? -eq 0 ]]; then 
    # move output files
    mv ${outpath}/${sample}/result.selfSM ${outpath}/${sample}.selfSM
    mv ${outpath}/${sample}/result.Ancestry ${outpath}/${sample}.Ancestry
    rm -r ${outpath}/${sample}
    
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi
