#!/bin/bash

input=$1
output=$2

ref="/path/to/Reference/Human_genome/GRCh38/GRCh38.fa"

minimap2 -L -t 36 --MD -a -x map-pb $ref $input | samtools sort -@ 4 -o $output
samtools index $output
