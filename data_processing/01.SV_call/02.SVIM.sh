#!/bin/bash

input=$1
outpath=$2

ref="/path/to/Reference/Human_genome/GRCh38/GRCh38.fa"

svim alignment \
$outpath \
$input \
$ref \
--max_sv_size 1000000
