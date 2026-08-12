library(readr)
library(dplyr)
library(robustbase)
library(ggplot2)

# this is the source data used in Fig. 5b
df_blood_ad <- read_tsv("input/fig5b_blood_ad/Blood.gene_expression_ad.txt")

# visualization
if(!dir.exists("output/fig5b_blood_ad")) dir.create("output/fig5b_blood_ad", recursive = T)
pdf("output/fig5b_blood_ad/fig5b.blood_ad.pdf")
ggplot(df_blood_ad, aes(AD, FAM193B_norm, group=AD)) +
    geom_boxplot(outlier.shape = NA, color="black") +
    geom_point(position=position_jitter(0.2))
dev.off()

# report statistics
# this reproduce the statistics reported in the manuscript
summ <- lmrob(FAM193B_norm ~ ., data=df_blood_ad) %>% summary()
cat("####### Statistics ########\n")
print(summ$coefficients['AD',])