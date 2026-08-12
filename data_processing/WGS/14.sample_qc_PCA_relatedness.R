# run this with srun -p cpu -n 1 -c 40 --pty bash
library(bigsnpr)
library(ggplot2)
library(dplyr)

options(bigstatsr.check.parallel.blas = FALSE)

bedfile <- "/path/to/sample_qc/HK_WGS_update.PASS.biallele.snp.hwe.autosome.geno1.maf1.bed"
rel_bedfile <- "/path/to/sample_qc/HK_WGS_update.PASS.raw.bed"
famfile <- "/path/to/sample_qc/HK_WGS_update.PASS.biallele.snp.hwe.autosome.geno1.maf1.fam"
plink2 <- "/path/to/plink2"
outpath <- "/path/to/Sample_QC/PCA_relatedness/"

# relateness
file_relatedness <- paste0(outpath, "QC_failed_relatedness.txt")
if (file.exists(file_relatedness)){
  remove_relatedness <- read.table(file_relatedness)[,1]
} else {

  rel <- snp_plinkKINGQC(
    plink2.path = plink2,
    bedfile.in = rel_bedfile,
    thr.king = 2^-3.5,
    make.bed = FALSE,
    ncores = nb_cores()
  )

  write.table(rel, paste0(outpath, "QC_failed_relatedness_table.txt"), quote=FALSE, row.names=FALSE)

  df_priority <- read.table("input/sample_QC_priority.txt", header=TRUE)
  df_sample_count <- data.frame(ID=c(rel$IID1, rel$IID2)) %>% count(ID) %>% arrange(desc(n))
  df_priority <- merge(df_priority, df_sample_count)

  remove_relatedness <- c()
  for (id in df_sample_count$ID){
    related_df <- filter(rel, IID1 == id | IID2 == id)
    related_df_remain <- filter(related_df, ! IID1 %in% remove_relatedness, ! IID2 %in% remove_relatedness)

    while (nrow(related_df_remain)>0){
      # rank by relation count and priority, remove top and loop
      remain_samples <- unique(c(related_df_remain$IID1, related_df_remain$IID2))
      remove_samples_df <- filter(df_priority, ID %in% remain_samples)
      remove_samples_df <- arrange(remove_samples_df, desc(n), priority)
      drop_sample <- remove_samples_df[1, "ID"]
      related_df_remain <- filter(related_df_remain, IID1 != drop_sample , IID2 != drop_sample)
      remove_relatedness <- c(remove_relatedness, drop_sample)
    }
  }

  write.table(remove_relatedness, file_relatedness, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# remove PCA outliers
obj.bed <- bed(bedfile)
sample_ids <- obj.bed$fam$sample.ID
ind.rel <- match(remove_relatedness, sample_ids)  # /!\ use $ID1 instead with old PLINK

ind.norel <- rows_along(obj.bed)[-ind.rel]

obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.norel, k = 20,
                       ncores = nb_cores())

# Outliers
prob <- bigutilsr::prob_dist(obj.svd$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

# histogram
pdf(paste0(outpath, "PCA_outliers_run1.pdf"))
ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(limits = c(0, NA), breaks=seq(0,1,0.1)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)")
dev.off()

# select cutoff by histgram
cutoff <- 0.7

# plot S by PCs
pdf(paste0(outpath, "PCA_outlier_PC_run1.pdf"))
plot_grid(plotlist = lapply(1:5, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)

# plot S by PCs applying cutoff based on histogram
plot_grid(plotlist = lapply(1:5, function(k) {
  plot(obj.svd, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S > cutoff) + 
    scale_colour_viridis_d()
}), scale = 0.95)
dev.off()


# Outliers iteration 2
# PCA without outliers
ind.row <- ind.norel[S < cutoff]
ind.col <- attr(obj.svd, "subset")
obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,
                        ind.col = ind.col, thr.r2 = NA,
                        k = 20, ncores = nb_cores())

prob <- bigutilsr::prob_dist(obj.svd2$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

# histogram
pdf(paste0(outpath, "PCA_outliers_run2.pdf"))
ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(limits = c(0, NA), breaks=seq(0,1,0.1)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)")
dev.off()

# select cutoff by histgram
cutoff <- 0.5

# plot S by PCs
pdf(paste0(outpath, "PCA_outlier_PC_run2.pdf"))
plot_grid(plotlist = lapply(1:5, function(k) {
  plot(obj.svd2, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S) +
    scale_colour_viridis_c()
}), scale = 0.95)

# plot S by PCs applying cutoff based on histogram
plot_grid(plotlist = lapply(1:5, function(k) {
  plot(obj.svd2, type = "scores", scores = 2 * k - 1:0, coeff = 0.6) +
    aes(color = S > cutoff) + 
    scale_colour_viridis_d()
}), scale = 0.95)
dev.off()

# Outliers iteration 3
ind.row <- ind.row[S < cutoff]
ind.col <- attr(obj.svd2, "subset")
obj.svd2 <- bed_autoSVD(obj.bed, ind.row = ind.row,
                        ind.col = ind.col, thr.r2 = NA,
                        k = 20, ncores = nb_cores())

prob <- bigutilsr::prob_dist(obj.svd2$u, ncores = nb_cores())
S <- prob$dist.self / sqrt(prob$dist.nn)

# histogram
pdf(paste0(outpath, "PCA_outliers_run3.pdf"))
ggplot() +
  geom_histogram(aes(S), color = "#000000", fill = "#000000", alpha = 0.5) +
  scale_x_continuous(limits = c(0, NA), breaks=seq(0,1,0.1)) +
  scale_y_sqrt(breaks = c(10, 100, 500)) +
  theme_bigstatsr() +
  labs(x = "Statistic of outlierness", y = "Frequency (sqrt-scale)")
dev.off()

# no outliers
# validation
pdf(paste0(outpath, "PCA_validation.pdf"))
plot(obj.svd2)
plot(obj.svd2, type = "loadings", loadings = 1:20, coef = 0.4)
plot(obj.svd2, type = "scores", scores = 1:20, coef = 0.4)
dev.off()

# write PCA
PCs <- predict(obj.svd2)
df_PCA <- data.frame(ID=sample_ids[ind.row])
df_PCA <- cbind(df_PCA, as.data.frame(PCs))
write.table(df_PCA, file.path(outpath, "PCA_vector.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

# write outliers
QC_outliers <- sample_ids[ind.norel[! ind.norel %in% ind.row]]
write.table(QC_outliers, file.path(outpath, "QC_failed_outliers.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

# write clean samples
QC_sex <- read.table(file.path("input/sex_error.fam"),header = FALSE)[,1]

clean_samples <- sample_ids[ind.row]
clean_samples <- clean_samples[! clean_samples %in% QC_sex]
write.table(clean_samples, file.path(outpath, "clean_samples.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
