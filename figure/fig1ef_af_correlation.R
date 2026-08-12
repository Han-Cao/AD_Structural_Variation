library(readr)
library(dplyr)
library(ggplot2)

theme_set(theme_classic() + 
          theme(text = element_text(size = 10, color = "black", face="plain"), 
                axis.line = element_line(linewidth = 0.2), axis.ticks = element_line(linewidth = 0.2)))

# these input files are generated using bcftools query on the corresponding VCFs 
# we keep those identified in both LRS and SRS with call rate > 80% for correlation analysis
df_af_lrs <- read_tsv("input/fig1ef_af_correlation/LRS_AF.SV.call_rate80.txt")
df_af_srs <- read_tsv("input/fig1ef_af_correlation/SRS_AF.SV.call_rate80.txt")

df_merge <- select(df_af_lrs, ID, SVTYPE, LRS_AF=AF) %>%
    left_join(select(df_af_srs, ID, SRS_AF=AF, concordant), by = "ID")

# AF correlation for raw genotyping results
cor_ins_raw <- cor.test(filter(df_merge, SVTYPE=="INS")$LRS_AF, filter(df_merge, SVTYPE=="INS")$SRS_AF)
cor_del_raw <- cor.test(filter(df_merge, SVTYPE=="DEL")$LRS_AF, filter(df_merge, SVTYPE=="DEL")$SRS_AF)

# AF correlation for concordant genotyping results
cor_ins_clean <- cor.test(filter(df_merge, concordant, SVTYPE=="INS")$LRS_AF, filter(df_merge, concordant, SVTYPE=="INS")$SRS_AF)
cor_del_clean <- cor.test(filter(df_merge, concordant, SVTYPE=="DEL")$LRS_AF, filter(df_merge, concordant, SVTYPE=="DEL")$SRS_AF)

# AF concordance
if (! dir.exists("output/fig1ef_af_correlation")) dir.create("output/fig1ef_af_correlation", recursive = T)

pdf("output/fig1ef_af_correlation/fig1e.AF_concordance_INS_raw.pdf", height = 2.2, width = 2.6)
ggplot(filter(df_merge, SVTYPE=="INS")) + 
    geom_hex(aes(x = LRS_AF, y = SRS_AF), binwidth = c(0.06, 0.04)) + 
    scale_fill_gradient(low = "white", high = "red", trans = "log10", limits=c(1,500), na.value = "red") +
    coord_fixed()
dev.off()

pdf("output/fig1ef_af_correlation/fig1e.AF_concordance_INS_clean.pdf", height = 2.2, width = 2.6)
ggplot(filter(df_merge, SVTYPE=="INS", concordant)) + 
    geom_hex(aes(x = LRS_AF, y = SRS_AF), binwidth = c(0.06, 0.05)) + 
    scale_fill_gradient(low = "white", high = "red", trans = "log10", limits=c(1,500), na.value = "red") +
    coord_fixed()
dev.off()

pdf("output/fig1ef_af_correlation/fig1f.AF_concordance_DEL_raw.pdf", height = 2.2, width = 2.6)
ggplot(filter(df_merge, SVTYPE=="DEL")) + 
    geom_hex(aes(x = LRS_AF, y = SRS_AF), binwidth = c(0.06, 0.04)) + 
    scale_fill_gradient(low = "white", high = "red", trans = "log10", limits=c(1,500), na.value = "red") +
    coord_fixed()
dev.off()

pdf("output/fig1ef_af_correlation/fig1f.AF_concordance_DEL_clean.pdf", height = 2.2, width = 2.6)
ggplot(filter(df_merge, SVTYPE=="DEL", concordant)) + 
    geom_hex(aes(x = LRS_AF, y = SRS_AF), binwidth = c(0.06, 0.05)) + 
    scale_fill_gradient(low = "white", high = "red", trans = "log10", limits=c(1,500), na.value = "red") +
    coord_fixed()
dev.off()

# report statistics
cat("####### Statistics #######\n")
cat("Fig1e: AF correlation for insertions (r)\n")
cat("raw: ", cor_ins_raw$estimate, "\n")
cat("filtered: ", cor_ins_clean$estimate, "\n\n")

cat("Fig1f: AF correlation for deletions (r)\n")
cat("raw: ", cor_del_raw$estimate, "\n")
cat("filtered: ", cor_del_clean$estimate, "\n")
