#!/bin/bash

inpath="/path/to/paragraph/output/"
outpath="/path/to/population/output/"

[[ ! -d ${outpath}/tmp ]] && mkdir -p ${outpath}/tmp

find $inpath -name 'genotypes.vcf.gz' | sort > ${outpath}/SRS.file_list.txt

# prepare GT-only tsv
parallel -j 10 '
vcf={1}
tmp_out={2}

name=$(dirname $vcf)
name=$(basename $name)

bcftools query -l $vcf > ${tmp_out}/${name}.txt
bcftools query -f "[%GT]\n" $vcf >> ${tmp_out}/${name}.txt

' :::: ${outpath}/SRS.file_list.txt ::: ${outpath}/tmp/

# generate VCF infos
vcf=$(head -n 1 ${outpath}/SRS.file_list.txt)

bcftools view -h $vcf | tail -n 1 | cut -f 1-9 > ${outpath}/merge_vcf_info.txt
bcftools view -H $vcf | cut -f 1-9 >> ${outpath}/merge_vcf_info.txt

# merge all VCFs
ulimit -n 2000

paste ${outpath}/merge_vcf_info.txt ${outpath}/tmp/* > ${outpath}/SRS.vcf

# concat vcf header
bcftools view -h $vcf | grep -v '^#CHROM' > ${outpath}/SRS.vcf.header
cat ${outpath}/SRS.vcf.header ${outpath}/SRS.vcf | bgzip > ${outpath}/SRS.SV.vcf.gz 
tabix ${outpath}/SRS.SV.vcf.gz 
