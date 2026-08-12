library(readr)
library(dplyr)
library(tidyr)

min_max_norm <- function(x, scale = 2){
  x_min <- min(x, na.rm = T)
  x_max <- max(x, na.rm = T)
  
  return ((x - x_min) / (x_max - x_min) * scale)
}


# read per snp local ancestry
df_ukb_snp_ancestry <- read_table("/path/to/UKB.FAM193B.local_ancestry.prob_full.txt.gz")
df_ukb_snp_ancestry <- mutate(df_ukb_snp_ancestry, 
                              eid = gsub("_[01]", "", sample) %>% as.numeric(), 
                              hap = gsub("\\d+_", "", sample) %>% as.numeric())
# extract per snp locanl ancestry on each haplotype
df_ukb_snp_ancestry_hap1 <- filter(df_ukb_snp_ancestry, hap==0)
df_ukb_snp_ancestry_hap2 <- filter(df_ukb_snp_ancestry, hap==1)
# check all sample and SNPs are well aligned
all(filter(df_ukb_snp_ancestry, hap==0)$eid == filter(df_ukb_snp_ancestry, hap==1)$eid)
all(filter(df_ukb_snp_ancestry, hap==0)$SNPidx == filter(df_ukb_snp_ancestry, hap==1)$SNPidx)
# calculate per-snp average ancestry (i.e., mean of 2 haps) and
df_ukb_snp_ancestry_merge <- mutate(df_ukb_snp_ancestry_hap1,
                                    EUR = (df_ukb_snp_ancestry_hap1$EUR + df_ukb_snp_ancestry_hap2$EUR) / 2,
                                    AFR = (df_ukb_snp_ancestry_hap1$AFR + df_ukb_snp_ancestry_hap2$AFR) / 2,
                                    EAS = (df_ukb_snp_ancestry_hap1$EAS + df_ukb_snp_ancestry_hap2$EAS) / 2,
                                    AMR = (df_ukb_snp_ancestry_hap1$AMR + df_ukb_snp_ancestry_hap2$AMR) / 2,
                                    SAS = (df_ukb_snp_ancestry_hap1$SAS + df_ukb_snp_ancestry_hap2$SAS) / 2,
) %>% select(-hap)

# map SNPidx to genetic positon and calcualte average ancestry by window
df_ukb_snp_map <- read_tsv("/path/to/UKB.FAM193B.map.txt")
df_ukb_snp_map <- tibble(SNPidx=seq_along(df_ukb_snp_map$`position(base)`),
                         pos=df_ukb_snp_map$`position(base)`)
df_ukb_snp_ancestry_merge <- left_join(df_ukb_snp_ancestry_merge, df_ukb_snp_map, by="SNPidx")

# 50kb windows in 200kb flanking region
window_breaks <- c(177554490-2e5, 177554490-1.5e5, 177554490-1e5, 177554490-0.5e5, 177554490,
                   177554490+0.5e5, 177554490+1e5, 177554490+1.5e5, 177554490+2e5)
df_ukb_region_ancestry <- mutate(df_ukb_snp_ancestry_merge,
                                 window = cut(pos, breaks = window_breaks, include.lowest = TRUE))
df_ukb_region_ancestry <- filter(df_ukb_region_ancestry, !is.na(window)) %>% group_by(eid, window) %>%
  summarise(EUR = mean(EUR), EAS = mean(EAS), AFR = mean(AFR), 
            AMR = mean(AMR), SAS = mean(SAS)) %>% ungroup()


# get clean samples for analysis

# prepare phenotype for regenie
df_regenie_pheno_covar <- read_tsv('/path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar.txt')

# add local ancestry to covar
levels(df_ukb_region_ancestry$window) <- seq(1,8)
df_ukb_region_ancestry <- pivot_wider(df_ukb_region_ancestry, 
                                      id_cols = eid, names_from = window, values_from = EUR:SAS)
df_regenie_pheno_covar_la <- inner_join(df_regenie_pheno_covar, df_ukb_region_ancestry, by=c("IID"="eid"))
write_tsv(df_regenie_pheno_covar_la, "/path/to/UKB.clinical_case_control.include_mixed_ancestry.pheno_covar_la.txt")
