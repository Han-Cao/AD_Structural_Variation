#!/bin/bash

# star 2.7.10b
# rsem 1.3.3

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

# set no. of thread used (should be the same as SBATCH -c)
thread=10

[[ ! -d "${outpath}/${sample}" ]] && mkdir -p "${outpath}/${sample}"

# set reference
star_ref="/path/to/star_ref"
rsem_ref="/path/to/rsem_ref"

# extract all input fastq
read1=$(ls ${inpath}/${sample}/*_1_paired.fq.gz | paste -s -d ',')
read2=$(echo $read1 | sed 's/_1_paired.fq.gz/_2_paired.fq.gz/g')

# STAR
STAR --genomeDir $star_ref \
--outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20  --outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000 \
--alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1 \
--sjdbScore 1  --runThreadN $thread  --genomeLoad NoSharedMemory  --outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM  --outSAMheaderHD @HD VN:1.4 SO:unsorted \
--outFileNamePrefix "${outpath}/${sample}/${sample}.STAR."  --readFilesCommand zcat \
--readFilesIn $read1 $read2

# rsem
rsem-calculate-expression --bam --no-bam-output -p $thread \
--paired-end --forward-prob 0.5 \
"${outpath}/${sample}/${sample}.STAR.Aligned.toTranscriptome.out.bam" \
$rsem_ref \
"${outpath}/${sample}/${sample}.rsem"

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi
