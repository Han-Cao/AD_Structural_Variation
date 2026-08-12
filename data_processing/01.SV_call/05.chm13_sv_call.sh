#!/bin/bash

ref_grch38="/path/to/Reference/Human_genome/GRCh38/GRCh38.fa"
ref_chm13="/path/to/Reference/Human_genome/CHM13/chm13v2.0_maskedY_rCRS.fa"
outpath="/path/to/chm13/"
[[ ! -d ${outpath} ]] && mkdir -p ${outpath}

# minimap2
minimap2 -a -x asm5 --cs -r2k -t 36 $ref_grch38 $ref_chm13 | \
samtools sort --write-index -m4G -@4 -o ${outpath}/chm13.sorted.bam -

# SVIM-asm call SVs
[[ ! -d ${outpath}/working_dir ]] && mkdir -p ${outpath}/working_dir
svim-asm haploid ${outpath}/working_dir ${outpath}/chm13.sorted.bam $ref_grch38
mv ${outpath}/working_dir/variants.vcf ${outpath}/chm13.svim-asm.vcf
bgzip ${outpath}/chm13.svim-asm.vcf
tabix ${outpath}/chm13.svim-asm.vcf.gz