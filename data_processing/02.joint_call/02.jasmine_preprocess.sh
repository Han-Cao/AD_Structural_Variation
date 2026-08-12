#!/bin/bash

inpath="/path/to/SV_call/harmonized/"
outpath="/path/to/SV_call/dup_call/"
ref="/path/to/Reference/Human/GRCh38/GRCh38.fa"

[[ ! -d $outpath ]] && mkdir -p $outpath

# sniffles, with STRAND
ls ${inpath}/*sniffles*vcf.gz | while read vcf;
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}

    # convert to vcf format
    bcftools view $vcf > ${outpath}/${name}.input.vcf

    jasmine file_list=${outpath}/${name}.input.vcf \
    out_file=${outpath}/${name}.preprocess.vcf \
    genome_file=$ref \
    --comma_filelist \
    max_dist=200 \
    --allow_intrasample \
    --nonlinear_dist \
    --keep_var_ids \
    threads=10

    # write merged SV list
    bcftools query -f '%INTRASAMPLE_IDLIST\n' ${outpath}/${name}.preprocess.vcf \
    > ${outpath}/${name}.preprocess.merge_list.txt

    rm ${outpath}/${name}.input.vcf
done

# svim, without STRAND
ls ${inpath}/*svim*vcf.gz | while read vcf;
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}

    # convert to vcf format
    bcftools view $vcf > ${outpath}/${name}.input.vcf

    jasmine file_list=${outpath}/${name}.input.vcf \
    out_file=${outpath}/${name}.preprocess.vcf \
    genome_file=$ref \
    --comma_filelist \
    max_dist=200 \
    --allow_intrasample \
    --nonlinear_dist \
    --ignore_strand \
    --keep_var_ids \
    threads=10

    # write merged SV list
    bcftools query -f '%INTRASAMPLE_IDLIST\n' ${outpath}/${name}.preprocess.vcf \
    > ${outpath}/${name}.preprocess.merge_list.txt

    rm ${outpath}/${name}.input.vcf
done

# cuteSV, without STRAND
ls ${inpath}/*cuteSV*vcf.gz | while read vcf;
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}

    # convert to vcf format
    bcftools view $vcf > ${outpath}/${name}.input.vcf

    jasmine file_list=${outpath}/${name}.input.vcf \
    out_file=${outpath}/${name}.preprocess.vcf \
    genome_file=$ref \
    --comma_filelist \
    max_dist=200 \
    --allow_intrasample \
    --nonlinear_dist \
    --ignore_strand \
    --keep_var_ids \
    threads=10

    # write merged SV list
    bcftools query -f '%INTRASAMPLE_IDLIST\n' ${outpath}/${name}.preprocess.vcf \
    > ${outpath}/${name}.preprocess.merge_list.txt

    rm ${outpath}/${name}.input.vcf
done