library(VariantAnnotation)
library(dplyr)
library(readr)

parse_vcf <- function(vcf){
    df_range <- rowRanges(vcf) %>% as_tibble() %>% select(seqnames, start, FILTER)
    df_range$ID <- rownames(vcf)
    df_range$strand <- "+"
    df_info <- info(vcf) %>% as_tibble() %>% select(SVTYPE, SVLEN)
    df_vcf <- bind_cols(df_range, df_info)
    df_vcf <- mutate(df_vcf, end = ifelse(SVTYPE=="INS", start, start+abs(SVLEN)-1))
}

vcf <- readVcf("/path/to/Joint_call/population/All.genotype.filter.population.AC0.vcf.gz", "hg38")
df_vcf <- parse_vcf(vcf)
df_vcf <- select(df_vcf, -FILTER, -SVLEN) %>%
          rename(chr=seqnames, type=SVTYPE) %>%
          mutate(chr=gsub("chr", "", chr)) %>%
          relocate(chr, start, end, type, strand, ID)
write_tsv(df_vcf, "/path/to/annotation/VEP_input/Chinese_SV.VEP_input.txt", col_names=FALSE)
