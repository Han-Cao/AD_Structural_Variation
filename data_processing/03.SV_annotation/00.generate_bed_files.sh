#!/bin/bash

# download and prepare annotation bed files
# run this script in annotation directory, the results have been provided in the input/ folder

# prepare gene annotation by gene biotypes
wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.chr.gtf.gz
Rscript -e '
library(readr)
library(dplyr)

# read gtf
df_gtf <- read_tsv("Homo_sapiens.GRCh38.107.chr.gtf.gz", col_names=FALSE, comment = "#")
colnames(df_gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "info")

# extract gene and generate biotype
df_gene <- filter(df_gtf, feature == "gene")
df_gene$biotype <- gsub(".+gene_biotype \"(.+)\";", "\\1", df_gene$info)
df_gene <- mutate(df_gene,
                  chr = paste0("chr", chr),
                  biotype_key = ifelse(grepl("pseudogene", biotype), "pseudogene",
                                ifelse(grepl("protein_coding", biotype), "protein_coding",
                                ifelse(grepl("(mi|sn|sno|sca|s)RNA", biotype), "small_RNA",
                                ifelse(grepl("lncRNA", biotype), "lncRNA",
                                ifelse(grepl("^(IG|TR)", biotype), "IG_TR", "others"))))))
# convert to 0 based bed file
df_gene$start <- df_gene$start - 1
df_gene <- filter(df_gene, biotype_key != "others") %>% select(chr, start, end, biotype_key)
write_tsv(df_gene, "GRCh38_gene_all.bed", col_names=FALSE)

gene_types <- unique(df_gene$biotype_key)
for (type in gene_types){
    df_type <- filter(df_gene, biotype_key == type)
    write_tsv(df_type, paste0("GRCh38_gene_", type, ".bed"), col_names=FALSE)
}
'

# candidate cis-regulatory elements
# grouped by types (PLS, pELS, dELS, CTCF)
wget https://api.wenglab.org/screen_v13/fdownloads/V3/GRCh38-cCREs.bed
mv GRCh38-cCREs.bed GRCh38-cCREs_v3.bed
Rscript --vanilla -e '
library(readr)
library(dplyr)

df_ccRE <- read_tsv("GRCh38-cCREs_v3.bed", col_names=FALSE)
ccRE_list <- c("PLS", "pELS", "dELS", "CTCF")
for (ccRE in ccRE_list){
    df_subset <- filter(df_ccRE, grepl(ccRE, X6))
    write_tsv(select(df_subset, X1:X3, X6), paste0("GRCh38_", ccRE, ".bed"), col_names=FALSE)
}
'

# GRCh38_MANE_*.bed:
# exon, intron, UTR annotations of MANE select
# download from UCSC table browser, MANE v1.0


# extract autosome features, merge overlaped intervals, and sort
# write chr1-22 X to /path/to/reference/annotation/clean_pattern.txt
ls /path/to/reference/annotation/*bed | while read file;
do
    if [[ ! $file == *".clean.merge.sort.bed" ]]; then
        name=${file%.bed}.clean.merge.sort.bed
        bedtools sort -i $file | \
        bedtools merge -i - -c 4 -o distinct | \
        cut -f 1-4 | \
        grep -f /path/to/reference/annotation/clean_pattern.txt > $name
    fi
done

# extract intergenic reigon
bedtools complement -L \
-i /path/to/reference/annotation/GRCh38_gene_all.clean.merge.sort.bed \
-g /path/to/reference/annotation/GRCh38.genome \
> /path/to/reference/annotation/GRCh38_gene_intergenic.clean.merge.sort.bed

awk '{OFS="\t"}{print $0, "intergenic"}' \
/path/to/reference/annotation/GRCh38_gene_intergenic.clean.merge.sort.bed \
> /path/to/reference/annotation/GRCh38_gene_intergenic.clean.merge.sort.bed.tmp

mv /path/to/reference/annotation/GRCh38_gene_intergenic.clean.merge.sort.bed.tmp \
/path/to/reference/annotation/GRCh38_gene_intergenic.clean.merge.sort.bed