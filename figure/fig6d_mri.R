library(readr)
library(dplyr)
library(robustbase)
library(RNOmni)
library(ggplot2)

# due to research ethics, we cannot publicly share the raw individual-level genotype data
# so we generate a synthetic data using the synthpop package
# this can provide similar results as reported in the manuscript, but NOT the same
df_gray_matter <- read_tsv('input/fig6d_mri/gray_matter_geno.txt')

# this is used to generate similar plots as Fig. 6d
df_gray_matter <- mutate(df_gray_matter,
                         CCG_interval = cut(max_nCCG, breaks=c(11, 14, 18, 22, Inf), right=FALSE))

# visualization
if (!dir.exists("output/fig6d_mri")) dir.create("output/fig6d_mri", recursive = T)
pdf("output/fig6d_mri/fig6d.gray_matter.synthpop.pdf")
ggplot(df_gray_matter, aes(CCG_interval, RankNorm(Gray_Matter))) + 
    geom_boxplot(outlier.shape = NA, color="black") +
    geom_point(position=position_jitter(0.2))
dev.off()

# report statistics
summ <- lmrob(RankNorm(Gray_Matter) ~ ., data=select(df_gray_matter, -CCG_interval)) %>% summary()
cat("####### Statistics from synthsized data ########\n")
print(summ$coefficients['max_nCCG',])