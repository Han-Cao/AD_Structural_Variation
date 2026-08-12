#!/bin/bash

# This script reproduce the SMR and HEIDI test results in Fig. 5e & Supplementary Table 20
# In the original analysis, we also performed multi-variant SMR test by --smr-multi and generated SMR plot data by --plot
# However, these additonal analysis require the UKB PLINK bfile which cannot be publicly shared

[[ ! -d 'output/fig5e_smr/' ]] && mkdir -p 'output/fig5e_smr/'

# Blood SMR
# AD assocaition: data_processing/06.candidate_locus/04.association_plink.sh
# Blood eQTL: data_processing/06.candidate_locus/04.association_plink.sh
smr \
--bld input/fig5e_smr/HK \
--gwas-summary input/fig5e_smr/HK.AD.ma \
--beqtl-summary input/fig5e_smr/HK.Blood_eQTL \
--target-snp FAM193B-CCG-ADD \
--ld-lower-limit 0.01 \
--out output/fig5e_smr/Blood_eQTL.FAM193B-CCG

# Brain SMR
# AD assocaition: data_processing/06.candidate_locus/05.association_regenie.sh
# Brain bulk eQTL: 10.1038/s41588-023-01300-6, 10.1038/s41588-024-02057-2
# Brain snRNA eQTL: fig5c-d_snrna_analysis.R

# Note: SMR and SMR-multi can produce slightly different HEIDI test results due to different SNP selection procedure
# In the original SMR-multi test, HEIDI test in brain is P = 0.259 (reported in the paper),
# while in this single-variant SMR, HEIDI test in brain is P = 0.245
smr \
--bld input/fig5e_smr/UKB \
--gwas-summary input/fig5e_smr/UKB.AD.ma \
--beqtl-summary input/fig5e_smr/EUR.Brain_eQTL \
--target-snp FAM193B-CCG-ADD \
--ld-lower-limit 0.01 \
--out output/fig5e_smr/Brain_eQTL.FAM193B-CCG

# visualization
Rscript --vanilla -e '
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

df_blood_smr <- read_tsv("output/fig5e_smr/Blood_eQTL.FAM193B-CCG.smr")
df_brain_smr <- read_tsv("output/fig5e_smr/Brain_eQTL.FAM193B-CCG.smr")

# merge SMR results
df_smr <- bind_rows(mutate(df_blood_smr, probeID="FAM193B-blood"),
                    df_brain_smr) %>%
        select(tissue=probeID, p_SMR, p_HEIDI)

# as stated above, SMR-multi cannot be done without individual-level PLINK bfile
# we use the existing results for visualization
tissue_list <- c("FAM193B-blood", "FAM193B-brain", 
                 "FAM193B-neuron", "FAM193B-astro", 
                 "FAM193B-oligo")
df_smr_multi <- tibble(tissue=tissue_list,
                       p_SMRmulti=c(2.3e-2, 5.8e-4, 3.4e-3, NA, NA))
df_smr <- inner_join(df_smr, df_smr_multi, by="tissue")


# plot
df_plot <- pivot_longer(df_smr, -tissue, names_to="test", values_to="p")
df_plot <- mutate(df_plot,
                  tissue = factor(tissue, levels=rev(tissue_list)),
                  test = factor(test, levels=c("p_HEIDI", "p_SMRmulti", "p_SMR")))

pdf("output/fig5e_smr/fig5e.smr_pvalues.pdf")
ggplot(df_plot) +
    geom_bar(aes(-log10(p), tissue, fill=test), stat="identity", position="dodge", color="black") +
    geom_vline(xintercept = -log10(0.05), linetype="dashed") +
    geom_vline(xintercept = -log10(0.05/5), linetype="dashed")
dev.off()
'
