#!/bin/bash

sample=$1
inpath=$2
outpath=$3
finish=$4
error=$5

refPre="/path/to/ref/GRCh38.fa"
dbsnp="/path/to/ref/dbsnp_153.GRCh38.vcf.gz"
g1000="/path/to/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="/path/to/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# interval list is available at https://github.com/Han-Cao/Genome_interval_list
singularity exec --nv \
-B "/path/to/ref/" \
-B "/path/to/work_dir/" \
--env TCMALLOC_MAX_TOTAL_THREAD_CACHE_BYTES=268435456 \
"/path/to/parabrick.4.0.1.sif" \
        pbrun haplotypecaller \
        --ref $refPre \
        --in-bam ${inpath}/${sample}.markdup.sort.cram \
        --in-recal-file ${inpath}/${sample}.markdup.sort.cram.recalibration_report.grp \
        --interval-file "/path/to/work_dir/code/input/wgs_calling_regions.hg38.interval_list" \
        --gvcf \
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 \
        -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        --out-variants ${outpath}/${sample}.g.vcf.gz \
        --num-gpus 1 

# mark finish or error
if [[ $? -eq 0 ]]; then 
    [[ -e "${error}/${sample}.error" ]] && mv "${error}/${sample}.error" "${error}/${sample}.error_fixed"
    touch "${finish}/${sample}.done"
else
    touch "${error}/${sample}.error"
fi
