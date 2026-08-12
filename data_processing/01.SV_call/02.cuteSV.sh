#!/bin/bash

input=$1
output=$2

[[ ! -d ${output}_wd ]] && mkdir -p ${output}_wd

cuteSV \
--max_cluster_bias_INS 100 \
--diff_ratio_merging_INS 0.3 \
--max_cluster_bias_DEL 200 \
--diff_ratio_merging_DEL 0.5 \
--genotype \
-t 20 \
$input \
"/path/to/Reference/Human_genome/GRCh38/GRCh38.fa" \
$output \
${output}_wd