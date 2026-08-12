library(readr)
library(dplyr)
library(ggplot2)

# due to research ethics, we cannot publicly share the raw individual-level genotype data
# so we generate a synthetic data using the synthpop package
# since the synthetic data can not be used by Regenie for whole-genome regression
# we here provide the code to do basic logistic regression (same as the Chinese cohorts)
# this can provide similar results as reported in the manuscript, but NOT the same

# In real data, the differences between the statistics from Regenie and GLM are very small:
# Method	Model	Beta	SE	Z	P
# GLM	Additive	0.011	0.004	2.922	3.5E-03
# GLM	Dominant	0.013	0.005	2.524	1.2E-02
# Regenie	Additive	0.015	0.005	2.981	2.9E-03
# Regenie	Dominant	0.012	0.004	2.604	9.2E-03


df_ukb <- read_tsv('input/fig4g_association_ukb/FAM193B_pheno_geno.UKB.synthpop.txt')

# logistic regression with all samples
# this generates similar results as Supplementary Table 13

# additive model
summ_add <- glm(AD ~ ., data=select(df_ukb, -max_nCCG), family = "binomial") %>% summary()
# dominant model
summ_dom <- glm(AD ~ ., data=select(df_ukb, -total_nCCG), family = "binomial") %>% summary()


# association by cutoffs to estiamte odds ratios
# this generates similar results as Fig. 4g and Supplementary Table 16
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

reg_str <- paste0('AD ~ high_nCCG + age + sex +', 
                  paste0("center", 1:3, collapse = " + "),
                  "+",
                  paste0("pc_", 1:20, collapse = " + "))
df_result_cutoff_ukb <- assoc_by_cutoff(data = df_ukb, 
                                        cutoff_range = seq(12, 33),
                                        reg = reg_str,
                                        target = "high_nCCGTRUE")

if (!dir.exists("output/fig4g_association_ukb/")) dir.create("output/fig4g_association_ukb/", recursive = T)
# synthetic data does not work well for rare events, so the statistics of very long cutoffs can differ from the real data
pdf("output/fig4g_association_ukb/fig4g.odds_ratio_ukb.synthpop.pdf")
ggplot(df_result_cutoff_ukb) +
    geom_point(aes(x=Cutoff, y=OR)) +
    geom_linerange(aes(x=Cutoff, ymin=CI_025, ymax=CI_975))
dev.off()

# report statistics
cat("######## Statistical results in synthsized data ########\n")
cat("# UKB (additive model)\n")
print(summ_add$coefficients['total_nCCG',])

cat("# UKB (dominant model)\n")
print(summ_dom$coefficients['max_nCCG',])