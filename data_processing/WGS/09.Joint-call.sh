#!/bin/bash

refPre="/path/to/GRCh38.fa"
dbsnp="/path/to/dbsnp_153.GRCh38.vcf.gz"
omni="/path/to/1000G_omni2.5.hg38.PASS.vcf.gz"
hapmap="/path/to/hapmap_3.3.hg38.vcf.gz"
g1000="/path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

[[ ! -d $outpath ]] && mkdir $outpath

gatk --java-options "-Xmx80G" GenotypeGVCFs \
-R $refPre \
-V gendb://"${inpath}/${sample}" \
-O $outpath"/${sample}.raw.vcf.gz" \
-dbsnp $dbsnp \
-G StandardAnnotation \
-G AS_StandardAnnotation \
--genomicsdb-shared-posixfs-optimizations

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi