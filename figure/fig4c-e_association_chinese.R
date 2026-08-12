library(readr)
library(dplyr)
library(ggplot2)

# due to research ethics, we cannot publicly share the raw individual-level genotype data
# so we generate a synthetic data using the synthpop package
# this can provide similar results as reported in the manuscript, but NOT the same
df_hk <- read_tsv('input/fig4c-e_association_chinese/FAM193B_pheno_geno.HK.synthpop.txt')
df_mainland <- read_tsv('input/fig4c-e_association_chinese/FAM193B_pheno_geno.Chinese_Mainland.synthpop.txt')

# 1. logistic regression with all samples
# this generates similar results as Fig. 4c and Supplementary Table 13
# total_nCCG for additive model, max_nCCG for dominant model

# Hong Kong
summ_hk_add <- glm(AD ~ total_nCCG + Age + Sex + PC1 + PC2 + PC3, data=df_hk, family = "binomial") %>% summary()
summ_hk_dom <- glm(AD ~ max_nCCG + Age + Sex + PC1 + PC2 + PC3, data=df_hk, family = "binomial") %>% summary()

# Chinese Mainland
summ_mainland_add <- glm(AD ~ total_nCCG + Age + Sex + PC1 + PC2 + PC3, data=df_mainland, family = "binomial") %>% summary()
summ_mainland_dom <- glm(AD ~ max_nCCG + Age + Sex + PC1 + PC2 + PC3, data=df_mainland, family = "binomial") %>% summary()

# 2. association by cutoffs to estiamte odds ratios
# this generates similar results as Fig. 4d-e and Supplementary Table 15
assoc_by_cutoff <- function(data, cutoff_range, reg, target){
  df_result <- tibble()
  for (x in cutoff_range) {
    # analysis by cutoff
    df_work <-  mutate(data, 
                       high_nCCG = if_else(max_nCCG <= 11, FALSE,
                                           if_else(max_nCCG >= x, TRUE, NA)))
    cutoff_fit <- glm(reg, data = df_work, family = "binomial")
    result <- summary(cutoff_fit)$coef[target,]
    new_row <- tibble(Cutoff=x, Beta = result[1], SE=result[2], z=result[3], P=result[4])
    df_result <- bind_rows(df_result, new_row)
  }
  
  df_result <- mutate(df_result,
                      OR = exp(Beta), CI_025 = exp(Beta - 1.96 * SE), CI_975 = exp(Beta + 1.96 * SE))
  return(df_result)
}

df_result_cutoff_hk <- assoc_by_cutoff(data = df_hk, 
                                       cutoff_range = seq(12, 28),
                                       reg = 'AD ~ high_nCCG + Age + Sex + PC1 + PC2 + PC3',
                                       target = "high_nCCGTRUE")
df_result_cutoff_mainland <- assoc_by_cutoff(data = df_mainland, 
                                             cutoff_range = seq(12, 23),
                                             reg = 'AD ~ high_nCCG + Age + Sex + PC1 + PC2 + PC3',
                                             target = "high_nCCGTRUE")

# visualization
# This generates similar statistics reported in Fig. 4d-e
if (!dir.exists("output/fig4c-e_association_chinese/")) dir.create("output/fig4c-e_association_chinese/", recursive = T)
pdf("output/fig4c-e_association_chinese/fig4d.odds_ratio_hong_kong.synthpop.pdf")
ggplot(df_result_cutoff_hk) +
    geom_point(aes(x=Cutoff, y=OR)) +
    geom_linerange(aes(x=Cutoff, ymin=CI_025, ymax=CI_975))
dev.off()

pdf("output/fig4c-e_association_chinese/fig4e.odds_ratio_chinese_mainland.synthpop.pdf")
ggplot(df_result_cutoff_mainland) +
    geom_point(aes(x=Cutoff, y=OR)) +
    geom_linerange(aes(x=Cutoff, ymin=CI_025, ymax=CI_975))
dev.off()

# report statistics
cat("######## Statistical results in synthsized data ########\n")
cat("# Hong Kong (additive model)\n")
print(summ_hk_add$coefficients['total_nCCG',])

cat("# Hong Kong (dominant model)\n")
print(summ_hk_dom$coefficients['max_nCCG',])

cat("# Chinese Mainland (additive model)\n")
print(summ_mainland_add$coefficients['total_nCCG',])

cat("# Chinese Mainland (dominant model)\n")
print(summ_mainland_dom$coefficients['max_nCCG',])
