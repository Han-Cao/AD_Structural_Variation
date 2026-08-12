library(readr)
library(dplyr)

# read AD-related data
# UKB_phenotypes.AD_record.txt.gz is generated following UKB's guide (https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/gwas-using-alzheimers-disease)
# note that the original guide has a bug in the following code:
# father_ad_risk = 1 if 10 in row['illnesses_of_father'] else np.maximum(0.32, (100 - row['father_age'])/100)
# mother_ad_risk = 1 if 10 in row['illnesses_of_mother'] else np.maximum(0.32, (100 - row['mother_age'])/100)
# we need to change np.max to np.min to match the algorithm described in Jansen et al. 2019
df_ukb_pheno <- read_tsv("/path/to/UKB_phenotypes.AD_record.txt.gz")
# QC on ethnic backgroud
df_ethnic <- read_csv("/path/to/ethnic_bakckground.csv")
keep_ethnic <- c("British", "Irish", "White")
# for local ancestry analysis
# keep_ethnic = c("British", "Irish", "White", "White and Asian", "Any other white background")
df_ethnic_keep <- filter(df_ethnic, 
                         (p21000_i0 %in% keep_ethnic) | 
                           (is.na(p21000_i0) & 
                              ((p21000_i1 %in% keep_ethnic) | 
                               (p21000_i2 %in% keep_ethnic) |
                               (p21000_i3 %in% keep_ethnic))))
qc_fail_non_eur <- filter(df_ukb_pheno, ! (eid %in% df_ethnic_keep$eid)) %>% pull(eid)
# QC on proxy score
ukb_proxy_eid <- filter(df_ukb_pheno, 
                        father_age > 0, mother_age > 0,
                        ! grepl('-11', illnesses_of_father), ! grepl('-13', illnesses_of_father),
                        ! grepl('-11', illnesses_of_mother), ! grepl('-13', illnesses_of_mother)) %>% pull(eid)
# all genetic QCs
df_ukb_pheno <- filter(df_ukb_pheno, 
                        sex == sex_genetic, 
                        is.na(sex_chr_aneuploidy),
                        kinship_code != 10,
                        is.na(outlier_het_missing),
                        !is.na(pc_1),
                        !(eid %in% qc_fail_non_eur))

df_ukb_pheno <- select(df_ukb_pheno, eid:sex, has_ad_icd10, ad_risk_by_proxy, pc_1:pc_20)
# in the setting of proxy AD, clinical diagnosed AD == 2, we convert it to 1 for logistic regression
df_ukb_pheno <- mutate(df_ukb_pheno, AD = if_else(has_ad_icd10==2, 1,0))

# get latest age
df_age <- read_csv("/path/to/demographic.csv.gz", col_types = cols(.default=col_double()))
df_age <- df_age %>% select(eid = eid,
                            age_i0 = p21003_i0,
                            age_i1 = p21003_i1,
                            age_i2 = p21003_i2,
                            age_i3 = p21003_i3,
                            year_birth = p34)
df_age <- df_age %>% rowwise %>% mutate(age_latest = max(age_i0, age_i1, age_i2, age_i3, na.rm=TRUE))
df_age <- ungroup(df_age) %>% select(eid, age = age_latest)

# include WGS batch, convert number to factor
df_wgs <- read_tsv("/path/to/WGS_QC_metrics_merge.txt")
df_wgs <- select(df_wgs, eid, sequencing_center) %>% mutate(sequencing_center = paste0("c", sequencing_center))

# exclude NC samples with specific neurological disorders
# see Supplemental Methods for details
ukb_exclude_samples <- read_tsv("/path/to/UKB.exclude_NC_diseases.samples.txt", col_names = FALSE)

# read genotype for QC
df_ukb_str_geno <- read_tsv("genotype/EnsembleTR.FAM193B_STR.UKB.genotype.txt", na = ".", col_names = FALSE)
colnames(df_ukb_str_geno) <- c("IID", "GT", "Score")
df_ukb_str_geno <- separate(df_ukb_str_geno, GT, c("a1_nCCG", "a2_nCCG"), sep=",")
# normalize ncopy
# 1. round to int
# 2. -3 as we don't count the 5' CGCCGCCG in ref. In most samples, they are CGCCGCG
df_ukb_str_geno <- mutate(df_ukb_str_geno, 
                      a1_nCCG = round(as.numeric(a1_nCCG)) - 3,
                      a2_nCCG = round(as.numeric(a2_nCCG)) - 3) %>% rowwise() %>%
                mutate(max_nCCG = max(a1_nCCG, a2_nCCG),
                       mean_nCCG = mean(a1_nCCG, a2_nCCG),
                       total_nCCG = sum(a1_nCCG, a2_nCCG))

# join all tables
df_ukb_pheno <- inner_join(df_ukb_pheno, df_ukb_str_geno, by=c("eid"="IID"))
df_ukb_pheno <- left_join(df_ukb_pheno, df_wgs, by="eid")
df_ukb_pheno <- left_join(df_ukb_pheno, df_age, by="eid")
# exclude samples without proxy score or with specific neurological disorders
df_ukb_pheno_clean <- filter(df_ukb_pheno, Score >= 0.4,
                              ! (AD==0 & eid %in% ukb_exclude_samples$X1),
                              eid %in% ukb_proxy_eid) %>% distinct()
# exclude outliers
df_ukb_pheno_clean <- filter(df_ukb_pheno_clean, total_nCCG <= 58)
df_ukb_pheno_clean_ad <- filter(df_ukb_pheno_clean, AD==1, age >= 60)
df_ukb_pheno_clean_nc <- filter(df_ukb_pheno_clean, AD==0, age >= 60, ad_risk_by_proxy < 0.3)
df_ukb_pheno_clean_case_control <- bind_rows(df_ukb_pheno_clean_ad, df_ukb_pheno_clean_nc)

# prepare files for regenie analysis
# convert factor to dummy variables
df_regenie_pheno_covar <- mutate(df_ukb_pheno_clean_case_control,
                                 center1 = if_else(sequencing_center == "c1", 1, 0),
                                 center2 = if_else(sequencing_center == "c2", 1, 0),
                                 center3 = if_else(sequencing_center == "c3", 1, 0)) %>%
                          select(FID=eid, IID=eid, AD, age, sex, center1:center3, pc_1:pc_20)
# convert repeat count to DS (range 0-2)
df_regenie_dosage <- select(df_ukb_pheno_clean_case_control, eid, a1_nCCG, a2_nCCG, max_nCCG, total_nCCG) %>%
                        mutate(DS_add = min_max_norm(total_nCCG), DS_dom = min_max_norm(max_nCCG))
write_tsv(df_regenie_pheno_covar, "/path/to/UKB.clinical_case_control.pheno_covar.txt")
write_tsv(df_regenie_dosage, "/path/to/UKB.clinical_case_control.dosage.txt")
write_tsv(select(df_regenie_pheno_covar, FID,IID), "/path/to/UKB.clinical_case_control.samples.txt")
