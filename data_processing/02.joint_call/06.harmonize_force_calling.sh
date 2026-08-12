#!/bin/bash

inpath="/path/to/Joint_call/force_call/raw/"
outpath="/path/to/Joint_call/force_call/harmonized/"

[[ ! -d $outpath ]] && mkdir -p $outpath

# sniffles
# keep basic SV INFO, harmonize STRAND
# calculate DP and RE
ls ${inpath}/*sniffles.vcf.gz | while read vcf;
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}

    python /path/to/code/harmoniSV/harmonizeVCF.py \
    -i $vcf \
    -o ${outpath}/${name}.harmonized.vcf.gz \
    --info SVTYPE,SVLEN,END \
    --format-to-info DP=DR,DP=DV,RE=DV \
    --sum \
    --header /path/to/SV_call/header/harmonized_header.txt \
    --id-prefix ${name} 
done

# cuteSV
# keep basic SV INFO, STDEV of position and length
# calculate DP and RE
ls ${inpath}/*cuteSV.vcf.gz | while read vcf;
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}

    python /path/to/code/harmoniSV/harmonizeVCF.py \
    -i $vcf \
    -o ${outpath}/${name}.harmonized.vcf.gz \
    --info SVTYPE,SVLEN,END,RE \
    --format-to-info DP=DR,DP=DV \
    --sum \
    --header /path/to/SV_call/header/harmonized_header.txt \
    --id-prefix ${name} 
done


ls ${outpath}/*harmonized.vcf.gz > ${outpath}/All_samples_pipelines.txt
bcftools concat -a -f ${outpath}/All_samples_pipelines.txt -Oz -o ${outpath}/All_samples_pipelines.vcf.gz
tabix ${outpath}/All_samples_pipelines.vcf.gz

# genearte table for post-processing
[[ ! -d ${outpath}/table/ ]] && mkdir -p ${outpath}/table/

ls ${outpath}/*.vcf.gz | while read vcf
do
    name=$(basename $vcf)
    name=${name%.vcf.gz}.txt
    echo -e "ID\tCHR\tPOS\tSVTYPE\tSVLEN\tAC\tDP\tRE" > ${outpath}/table/${name}
    bcftools query -f '%ID\t%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%AC\t%DP\t%RE\n' $vcf >> ${outpath}/table/${name}
    gzip -f ${outpath}/table/${name}
done