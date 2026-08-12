#!/bin/bash

[[ ! -d output/fig4f_meta_analysis ]] && mkdir -p output/fig4f_meta_analysis

# association by repeat dosage
# this reproduce the results in Supplemental Table 13
java -jar input/fig4f_meta_analysis/Metasoft.jar \
-input input/fig4f_meta_analysis/dosage_effect.txt \
-output output/fig4f_meta_analysis/dosage_effect.meta.txt \
-pvalue_table input/fig4f_meta_analysis/HanEskinPvalueTable.txt \
-log output/fig4f_meta_analysis/dosage_effect.meta.log

# association by repeat cutoff to estiamte odds ratios
# this reproduce the results in Fig. 4f
java -jar input/fig4f_meta_analysis/Metasoft.jar \
-input input/fig4f_meta_analysis/cutoff_effect.txt \
-output output/fig4f_meta_analysis/cutoff_effect.meta.txt \
-pvalue_table input/fig4f_meta_analysis/HanEskinPvalueTable.txt \
-log output/fig4f_meta_analysis/cutoff_effect.meta.log

# visualization

Rscript --vanilla -e '
library(readr)
library(dplyr)
library(ggplot2)

df_or_hk <- read_tsv("input/fig4f_meta_analysis/cutoff_effect.txt", 
                     col_names=c("cutoff", "beta_hk", "se_hk", "beta_mainland", "se_mainland"))
df_or_meta <- read_tsv("output/fig4f_meta_analysis/cutoff_effect.meta.txt")

# convert beta to odds ratio
df_or_hk <- mutate(df_or_hk, 
                   or = exp(beta_hk),
                   or_025 = exp(beta_hk - 1.96 * se_hk),
                   or_975 = exp(beta_hk + 1.96 * se_hk))

df_or_meta <- mutate(df_or_meta,
                     or = exp(BETA_FE),
                     or_025 = exp(BETA_FE - 1.96 * STD_FE),
                     or_975 = exp(BETA_FE + 1.96 * STD_FE))

df_plot <- bind_rows(filter(df_or_meta, !is.na(or)) %>% select(cutoff=RSID, or, or_025, or_975),
                     filter(df_or_hk, is.na(beta_mainland)) %>% select(cutoff, or, or_025, or_975))

# visualize

pdf("output/fig4f_meta_analysis/fig4f.meta_analysis.pdf")
ggplot(df_plot) + 
    geom_point(aes(cutoff, or)) +
    geom_linerange(aes(cutoff, ymin=or_025, ymax=or_975))
dev.off()
'