#!/bin/bash

input=$1
output=$2

ref="/path/to/Reference/Human_genome/GRCh38/GRCh38.fa"

ngmlr --bam-fix -t 40 -x pacbio \
-r $ref \
-q $input \
-o ${output%.bam}.sam

# not compatable with newer samtools version
/home/hcaoad/Software/samtools-1.9/samtools sort -@ 10 \
-o ${output} ${output%.bam}.sam

/home/hcaoad/Software/samtools-1.9/samtools index ${output}