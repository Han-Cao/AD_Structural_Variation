#!/bin/bash

# Use the PacBio HiFi data to benchmark the SV genotyping results in SRS

ref_vcf="/path/to/represent/All.genotype.filter.population.AC0.clean.concordance.siteonly.vcf.gz"
paragraph_vcf="/path/to/Paragraph.SV.vcf.gz"
bam="/path/to/hifi_sample.minimap2.GRCh38.cram"
outpath="/path/to/genotyping/benchmark/"

[[ ! -d $outpath ]] && mkdir -p $outpath

# 1. Re-genotyping by force-calling
# sniffles
sniffles \
-i $bam \
-v ${outpath}/hifi_sample.hifi.sniffles.vcf \
--reference "/path/to/Human_genome/GRCh38/GRCh38.fa" \
-t 40 \
--tandem-repeats "/path/to/Human_genome/GRCh38/human_GRCh38_no_alt_analysis_set.trf.bed" \
--genotype-vcf $ref_vcf

# cuteSV
mkdir -p ${outpath}/hifi_sample.cuteSV.vcf_wd
cuteSV \
--max_cluster_bias_INS 1000 \
--diff_ratio_merging_INS	0.9 \
--max_cluster_bias_DEL	1000 \
--diff_ratio_merging_DEL	0.5 \
--genotype \
-t 40 \
-L -1 \
-Ivcf $ref_vcf \
$bam \
"/path/to/Human_genome/GRCh38/GRCh38.fa" \
${outpath}/hifi_sample.hifi.cuteSV.vcf \
${outpath}/hifi_sample.cuteSV.vcf_wd
rm -r ${outpath}/hifi_sample.cuteSV.vcf_wd

# extract genotypes
bcftools query -f '%ID\t%FILTER\t[%GT\t%DR\t%DV]\n' ${outpath}/hifi_sample.hifi.sniffles.vcf > ${outpath}/hifi_sample.hifi.sniffles.genotype.txt
bcftools query -f '%ID\t%FILTER\t[%GT\t%DR\t%DV]\n' ${outpath}/hifi_sample.hifi.cuteSV.vcf > ${outpath}/hifi_sample.hifi.cuteSV.genotype.txt
bcftools view -s hifi_sample ${paragraph_vcf} | bcftools query -f '%ID\t[%GT]\n' > ${outpath}/hifi_sample.paragraph.genotype.txt


# 2. Benchmark
cd ${outpath}
Rscript --vanilla -e '
library(readr)
library(dplyr)

df_sniffles <- read_tsv("hifi_sample.hifi.sniffles.genotype.txt", col_names=c("ID", "FILTER", "GT", "DR", "DV"), na=".")
df_cuteSV <- read_tsv("hifi_sample.hifi.cuteSV.genotype.txt", col_names=c("ID", "FILTER", "GT", "DR", "DV"), na=".")
df_paragraph <- read_tsv("hifi_sample.paragraph.genotype.txt", col_names=c("ID", "GT"), na=".")

# hard filter hifi results by DP >= 10
# merge het and hom calls for benchmark (same as truvari)
df_sniffles <- df_sniffles %>% mutate(GT_sniffles = ifelse( (DR+DV ) >= 10, if_else(GT %in% c("0/1", "1/0", "1/1"), 1, 0), NA))
df_cuteSV <- df_cuteSV %>% mutate(GT_cuteSV = ifelse( (DR+DV ) >= 10, if_else(GT %in% c("0/1", "1/0", "1/1"), 1, 0), NA))
# paragraph may have NA genotypes, we need to handle it
df_paragraph <- df_paragraph %>% mutate(GT_paragraph = if_else(is.na(GT), NA, if_else(GT %in% c("0/1", "1/0", "1/1"), 1, 0)))

# only consider the concordant genotype between sniffles and cuteSV as the truth
df_truth <- inner_join(select(df_sniffles, ID, GT_sniffles), 
                       select(df_cuteSV, ID, GT_cuteSV), by="ID")
df_truth <- mutate(df_truth, GT_truth = ifelse(GT_sniffles == GT_cuteSV, GT_sniffles, NA))

# benchmark paragraph results
df_compare <- inner_join(df_truth, select(df_paragraph, ID, GT_paragraph), by="ID")
df_compare <- mutate(df_compare, SVTYPE = gsub("Chinese_SV\\.(\\w+)\\.chr.+", "\\1", ID))

df_benchmark <- group_by(df_compare, SVTYPE) %>%
              summarise(TP = sum(GT_truth == GT_paragraph & GT_paragraph == 1, na.rm=TRUE),
                        FP = sum(GT_truth != GT_paragraph & GT_paragraph == 1, na.rm=TRUE),
                        TN = sum(GT_truth == GT_paragraph & GT_paragraph == 0, na.rm=TRUE),
                        FN = sum(GT_truth != GT_paragraph & GT_paragraph == 0, na.rm=TRUE),
                        Precision = TP/(TP+FP),
                        Recall = TP/(TP+FN),
                        F1 = 2*Precision*Recall/(Precision+Recall))
'