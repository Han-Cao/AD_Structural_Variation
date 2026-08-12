#!/bin/bash

# overlap with different features
bedtools intersect \
-a input/Chinese_SV.bed \
-b input/GRCh38_gene_all.clean.merge.sort.bed \
   input/GRCh38_gene_protein_coding.clean.merge.sort.bed \
   input/GRCh38_gene_lncRNA.clean.merge.sort.bed \
   input/GRCh38_gene_small_RNA.clean.merge.sort.bed \
   input/GRCh38_gene_IG_TR.clean.merge.sort.bed \
   input/GRCh38_gene_pseudogene.clean.merge.sort.bed \
   input/GRCh38_gene_intergenic.clean.merge.sort.bed \
   input/GRCh38_MANE_coding_exon.clean.merge.sort.bed \
   input/GRCh38_MANE_intron.clean.merge.sort.bed \
   input/GRCh38_MANE_UTR5.clean.merge.sort.bed \
   input/GRCh38_MANE_UTR3.clean.merge.sort.bed \
   input/GRCh38_PLS.clean.merge.sort.bed \
   input/GRCh38_pELS.clean.merge.sort.bed \
   input/GRCh38_dELS.clean.merge.sort.bed \
   input/GRCh38_CTCF.clean.merge.sort.bed \
-names gene_all protein_coding lncRNA small_RNA IG_TR pseudogene intergenic exon intron UTR5 UTR3 PLS pELS dELS CTCF \
-sorted \
-wao \
> /path/to/Chinese_SV.overlap_feature.bed
