#!/bin/bash

input=$1
output=$2

sniffles \
-i $input \
-v $output \
--reference "/path/to/Reference/Human_genome/GRCh38/GRCh38.fa" \
-t 40 \
--tandem-repeats "/path/to/Reference/Human_genome/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed"

