library(readr)
library(dplyr)
library(ggplot2)

# due to research ethics, we cannot publicly share the raw individual-level genotype data
# so we generate a synthetic data using the synthpop package
# this can provide similar results as reported in the manuscript, but NOT the same
df_eqtl <- read_tsv("input/fig5a_blood_eqtl/Blood.gene_expression_geno.txt")

# this generates source data similar as the one used in in Fig. 5a
df_plot <- select(df_eqtl, FAM193B_norm, total_nCCG) %>% 
        mutate(mean_nCCG = total_nCCG / 2,
               CCG_interval = cut(mean_nCCG, breaks=c(7, 10, 14, 18, 22, Inf), right=FALSE))

# visualization
if(!dir.exists("output/fig5a_blood_eqtl")) dir.create("output/fig5a_blood_eqtl", recursive = T)
pdf("output/fig5a_blood_eqtl/fig5a.blood_eqtl.synthpop.pdf")
ggplot(df_plot, aes(CCG_interval, FAM193B_norm)) + 
    geom_boxplot(outlier.shape = NA, color="black") +
    geom_point(position=position_jitter(0.2))
dev.off()


# report statistics
# please note that given the real association is too significant, synthpop does not work well to generate similar P value
summ <- lm(FAM193B_norm ~ ., data=df_eqtl) %>% summary()
cat("####### Statistics from synthsized data ########\n")
print(summ$coefficients['total_nCCG',])
