#!/bin/bash

# run VEP to annotate SVs
# the output file is given in output/Chinese_SV.VEP_output.txt.gz

inpath="/path/to/annotation/VEP_input/"
outpath="/path/to/annotation/VEP_output/"

ccRE="/path/to/reference/annotation/GRCh38-cCREs_v3.vep.bed.gz"

vep --cache --dir_cache /path/to/Reference/vep/ \
--offline \
--force_overwrite \
--canonical \
--mane \
--tsl \
--custom ${ccRE},ccRE,bed,overlap \
--symbol \
--biotype \
--fasta "/path/to/Reference/Human_genome/GRCh38/GRCh38.fa" \
--overlaps \
--fork 4 \
-i "${inpath}/Chinese_SV.VEP_input.txt" \
-o "${outpath}/Chinese_SV.VEP_output.txt"
