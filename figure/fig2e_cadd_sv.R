library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# LRS.SV.txt is generated using bcftools query from the population VCF
# Chinese_SV.CADD.txt is generated suing CADD-SV's web service

# read SV INFO
df_sv <- read_tsv("input/fig2e_cadd_sv/LRS.SV.txt", col_names=FALSE)
colnames(df_sv) <- c("ID", "SVLEN", "SVTYPE", "AC", "AN")
df_sv <- mutate(df_sv, AF=AC/AN, MAF = ifelse(AF>0.5, 1-AF, AF))

# read CADD
df_cadd <- read_tsv("input/fig2e_cadd_sv/Chinese_SV.CADD.txt")
df_cadd <- left_join(df_sv, select(df_cadd, ID:raw_score_flank), by="ID") %>% 
           filter(!is.na(CADD), MAF > 0)
df_cadd <- mutate(df_cadd, 
                  CADD_interval = cut(CADD, breaks=c(seq(0,20,2), 50), include.lowest=TRUE))

# annotate LoF SVs
df_lof <- read_tsv('input/fig2e_cadd_sv/SV_annotation_protein_coding_MANE_LoF.txt')
# remove those likely gain of function
df_lof <- filter(df_lof, ! Consequence %in% c("coding_sequence_variant,feature_elongation", "coding_sequence_variant,feature_elongation"))
df_cadd$lof <- df_cadd$ID %in% df_lof$ID

# read GWAS overlap results (code similar to fig.2c-d)
df_gwas <- read_tsv("input/fig2cd_gwas_overlap/Chinese_SV.overlap_GWAS.10kb.bed", 
                    col_names=FALSE, na=c(".", "NA"))
colnames(df_gwas) <- c("chr", "start", "end", "ID", "overlap_chr", "overlap_start", "overlap_end", "overlap_ID", "overlap_len")
df_gwas <- filter(df_gwas, !is.na(overlap_ID))
df_gwas$phenotype <- gsub('.+:(.+)', '\\1', df_gwas$overlap_ID)
df_gwas$gene <- gsub('(.+?):.+', '\\1', df_gwas$overlap_ID)
disease_list <- c("Cancer", "Cardiovascular_disease", "Digestive_system_disorder",
                  "Immune_system_disorder", "Immune_system_disorder", "Metabolic_disorder",
                  "Neurological_disorder", "Other_disease")
df_gwas <- filter(df_gwas, phenotype %in% disease_list)

# group SVs by GWAS loci and LoF
df_cadd <- mutate(df_cadd, 
                  gwas = ID %in% df_gwas$ID, 
                  gwas_lof = ID %in% filter(df_lof, SYMBCOL %in% df_gwas$gene)$ID,
                  risk_group=ifelse(!gwas, "Non-GWAS", ifelse(gwas_lof, "GWAS-LoF", "GWAS-non-LoF")))

# summarise LoF and GWAS by CADD interval
df_cadd_summary <- filter(df_cadd) %>% group_by(CADD_interval) %>% 
                summarise(lof=mean(lof), gwas=mean(gwas), gwas_lof=mean(gwas_lof), n = n()) %>%
                mutate(lof_non_gwas = lof - gwas_lof)
df_cadd_plot <- select(df_cadd_summary, CADD_interval, gwas, lof_non_gwas) %>%
            mutate(other = 1 - gwas - lof_non_gwas) %>%
            pivot_longer(-CADD_interval, names_to="group", values_to="freq") %>%
            mutate(group = factor(group, levels=c("other", "gwas", "lof_non_gwas")))

if (!dir.exists("output/fig2e_cadd_sv")) dir.create("output/fig2e_cadd_sv", recursive = T)
pdf("output/fig2e_cadd_sv/fig2e.cadd_sv.pdf")
ggplot(df_cadd_plot, aes(CADD_interval, freq, fill=group)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#999999", "#fbb040", "#ed1c24"))
dev.off()