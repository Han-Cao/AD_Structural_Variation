#!/bin/bash

# overlap with GWAS loci within 10kb window
awk '{OFS="\t"} {print $1,$2-10000,$3+10000,$6}' input/Chinese_SV.bed > input/Chinese_SV.10kb.bed

bedtools intersect \
-a input/Chinese_SV.10kb.bed \
-b input/GWAS_catalog_Gene.hg38.bed \
-names GWAS_catalog \
-wao \
> /path/to/GWAS/overlap/Chinese_SV.overlap_GWAS.10kb.bed
