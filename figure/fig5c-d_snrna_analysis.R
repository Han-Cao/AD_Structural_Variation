# FAM193B analysis
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lme4)
library(MAST)
library(Seurat)

# 0. read data
# the individual-level raw data is not publicly available, we only provide code here
# read snRNA-seq data
inhouse_raw <- readRDS("/path/to/annotated_brain_S102.rds")

# read FAM193B-CCG repeat
df_fam193b_dosage <- read_tsv("Input/genotype/EnsembleTR.FAM193B_STR.brain.genotype.txt", col_names=FALSE)
colnames(df_fam193b_dosage) <- c("Sample", "GT", "Score")
df_fam193b_dosage <- separate(df_fam193b_dosage, GT, c("a1_nCCG", "a2_nCCG"), sep=",")
# normalize ncopy
# 1. round to int
# 2. -3 as we don't count the 5' CGCCGCCG in ref. In most samples, they are CGCCGCG
df_fam193b_dosage <- mutate(df_fam193b_dosage, 
                          a1_nCCG = round(as.numeric(a1_nCCG)) - 3,
                          a2_nCCG = round(as.numeric(a2_nCCG)) - 3) %>% rowwise() %>%
  mutate(max_nCCG = max(a1_nCCG, a2_nCCG),
         mean_nCCG = mean(a1_nCCG, a2_nCCG),
         total_nCCG = sum(a1_nCCG, a2_nCCG))
df_fam193b_dosage <- left_join(df_fam193b_dosage, df_pca, by=c("Sample"="ID"))
df_fam193b_dosage <- filter(df_fam193b_dosage, Score >= 0.4)

# read other genetic data
df_fam193b_snp <- read_tsv("Input/genotype/FAM193B.200kb.txt")
df_pca <- read_tsv("Input/genotype/UKBBN_PCA.txt")

# 2. QC
# only keep protein coding genes
protein_coding_gene <- read_tsv("Input/HGNC_protein_coding.txt")
all_gene <- rownames(inhouse_raw)
keep_gene <- all_gene[all_gene %in% protein_coding_gene$symbol]
inhouse_raw <- inhouse_raw[keep_gene,]

# get expression
df_cell <- as_tibble(inhouse_raw@meta.data, rownames = "cell_id")
df_cell$FAM193B <- GetAssayData(inhouse_raw, "RNA", layer = "data")["FAM193B",] %>% as.numeric()

# merge neurons
df_cell[df_cell$celltype %in% c("excit", "inhibit"), "celltype"] <- "neuron"
df_cell <- left_join(df_cell, select(df_fam193b_dosage, Sample, max_nCCG, total_nCCG, PC1, PC2, PC3), by=c("stim" = "Sample"))

# annotate genotypes and merged celltypes
inhouse_raw$celltype <- df_cell$celltype
inhouse_raw$max_nCCG <- df_cell$max_nCCG
inhouse_raw$total_nCCG <- df_cell$total_nCCG

# get gene detection rate (as suggested by MAST)
expr_matrix <- GetAssayData(inhouse_raw, layer = "data", assay = "RNA")
detection_rate <- colSums(expr_matrix>0) / nrow(expr_matrix)
detection_rate <- scale(detection_rate) %>% as.numeric
df_cell$detection_rate <- detection_rate

# normalize covariates
df_cell$Age <- scale(df_cell$Age) %>% as.numeric()
df_cell$PMD <- scale(df_cell$PMD) %>% as.numeric()
df_cell$PC1 <- scale(df_cell$PC1) %>% as.numeric()
df_cell$PC2 <- scale(df_cell$PC2) %>% as.numeric()
df_cell$PC3 <- scale(df_cell$PC3) %>% as.numeric()


# 3. Subpopulation analysis
# to increase resolution, we redo clustering within different cell types
# we identify neuron of different layer using marker genes from https://www.science.org/doi/10.1126/science.adh1938
# for other cell types, we just name them by number
DefaultAssay(inhouse_raw) <- "integrated"

# function to extract statistics from MAST results
MAST_summary <- function(MAST_fit, doLRT){
  MAST_summary <- summary(MAST_fit, doLRT = doLRT)
  MAST_table <- MAST_summary$datatable
  merge(MAST_table[contrast==doLRT & component=="logFC", .(primerid, coef, ci.hi, ci.lo)],
        MAST_table[contrast==doLRT & component=="H", .(primerid, `Pr(>Chisq)`)])
}

# 3.1 neuron
# within neuron clustering
inhouse_neuron <- subset(inhouse_raw, celltype == "neuron")
inhouse_neuron <- FindNeighbors(inhouse_neuron, reduction = "pca", dims = 1:30)
inhouse_neuron <- FindClusters(inhouse_neuron, resolution = 0.05)

df_markers_neuron <- FindAllMarkers(inhouse_neuron, assay = "RNA", min.pct = 0.1, only.pos = TRUE)

# compare with known neuron markers from spatial data (https://www.science.org/doi/10.1126/science.adh1938)
df_neuron_spatial_marker <- read_tsv("Table/FAM193B/neuron_markers.txt")
df_markers_neuron <- left_join(df_markers_neuron, df_neuron_spatial_marker, by=c("gene"="Gene"))
# manually check top marker genes and assign subtypes
# Mixed cells are small cluster of outliers with multiple subtype markers
df_markers_neuron %>% filter(!is.na(Celltype), avg_log2FC > 1, pct.1 > pct.2)
df_neuron_cluster <- tibble(seurat_clusters=factor(c(0:12)),
                            subtype=c("Excit_L3",
                                      "Excit_L3/4/5",
                                      "Inhibit",
                                      "Inhibit",
                                      "Excit_L4",
                                      "Excit_L5",
                                      "Neuron_Mixed",
                                      "Excit_L5/6",
                                      "Excit_L6",
                                      "Neuron_Mixed",
                                      "Neuron_Mixed",
                                      "Excit_L5/6",
                                      "Neuron_Mixed"))
df_neuron_subtype <- as_tibble(inhouse_neuron@meta.data, rownames = "cell_id") %>% select(cell_id, seurat_clusters)
df_neuron_subtype <- left_join(df_neuron_subtype, df_neuron_cluster, by=c("seurat_clusters"))

inhouse_neuron$subtype <- df_neuron_subtype$subtype

# QTL
df_result_subtype_neuron <- tibble()
for (x in unique(df_neuron_subtype$subtype)) {
  df_work <- filter(df_neuron_subtype, subtype==x)
  df_work <- inner_join(df_cell, df_work, by="cell_id")
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_QTL_total <- mutate(res_QTL_total, subtype = x, n_cell = nrow(df_work), n_ind_cell50 = sum(table(df_work$stim) > 50))
  df_result_subtype_neuron <- bind_rows(df_result_subtype_neuron, res_QTL_total)
}

# 3.2 Astrocytes
inhouse_astro <- subset(inhouse_raw, celltype == "astro")
# within astrocyte clustering
inhouse_astro <- FindNeighbors(inhouse_astro, reduction = "pca", dims = 1:30)
inhouse_astro <- FindClusters(inhouse_astro, resolution = 0.05)

df_astro_cluster <- tibble(seurat_clusters=factor(c(0:1)), subtype=c("Ast_1", "Ast_2"))
df_astro_subtype <- as_tibble(inhouse_astro@meta.data, rownames = "cell_id") %>% select(cell_id, seurat_clusters)
df_astro_subtype <- left_join(df_astro_subtype, df_astro_cluster, by=c("seurat_clusters"))

inhouse_astro$subtype <- df_astro_subtype$subtype

# QTL
df_result_subtype_astro <- tibble()
for (x in unique(df_astro_subtype$subtype)) {
  df_work <- filter(df_astro_subtype, subtype==x)
  df_work <- inner_join(df_cell, df_work, by="cell_id")
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_QTL_total <- mutate(res_QTL_total, subtype = x, n_cell = nrow(df_work), n_ind_cell50 = sum(table(df_work$stim) > 50))
  df_result_subtype_astro <- bind_rows(df_result_subtype_astro, res_QTL_total)
}

# 3.3 Microglia
inhouse_micro <- subset(inhouse_raw, celltype == "mic")
# within microglia clustering
inhouse_micro <- FindNeighbors(inhouse_micro, reduction = "pca", dims = 1:50)
inhouse_micro <- FindClusters(inhouse_micro, resolution = 0.1)
# "Mixed" are small clusters of outliers
df_micro_cluster <- tibble(seurat_clusters=factor(c(0:3)), 
                           subtype=c("Mic_1", "Mic_2", "Mic_Mixed", "Mic_Mixed"))
df_micro_subtype <- as_tibble(inhouse_micro@meta.data, rownames = "cell_id") %>% select(cell_id, seurat_clusters)
df_micro_subtype <- left_join(df_micro_subtype, df_micro_cluster, by=c("seurat_clusters"))

inhouse_micro$subtype <- df_micro_subtype$subtype

# QTL
df_result_subtype_micro <- tibble()
for (x in unique(df_micro_subtype$subtype)) {
  df_work <- filter(df_micro_subtype, subtype==x)
  df_work <- inner_join(df_cell, df_work, by="cell_id")
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_QTL_total <- mutate(res_QTL_total, subtype = x, n_cell = nrow(df_work), n_ind_cell50 = sum(table(df_work$stim) > 50))
  df_result_subtype_micro <- bind_rows(df_result_subtype_micro, res_QTL_total)
}

# 3.4 Oligodendrocytes
inhouse_oligo <- subset(inhouse_raw, celltype == "oligo")
# within oligodendrocyte clustering
inhouse_oligo <- FindNeighbors(inhouse_oligo, reduction = "pca", dims = 1:30)
inhouse_oligo <- FindClusters(inhouse_oligo, resolution = 0.05)
# "Mixed" are small clusters of outliers
df_oligo_cluster <- tibble(seurat_clusters=factor(c(0:2)), 
                           subtype=c("Oligo_1", "Oligo_2", "Oligo_Mixed"))
df_oligo_subtype <- as_tibble(inhouse_oligo@meta.data, rownames = "cell_id") %>% select(cell_id, seurat_clusters)
df_oligo_subtype <- left_join(df_oligo_subtype, df_oligo_cluster, by=c("seurat_clusters"))

inhouse_oligo$subtype <- df_oligo_subtype$subtype

# QTL
df_result_subtype_oligo <- tibble()
for (x in unique(df_oligo_subtype$subtype)) {
  df_work <- filter(df_oligo_subtype, subtype==x)
  df_work <- inner_join(df_cell, df_work, by="cell_id")
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_QTL_total <- mutate(res_QTL_total, subtype = x, n_cell = nrow(df_work), n_ind_cell50 = sum(table(df_work$stim) > 50))
  df_result_subtype_oligo <- bind_rows(df_result_subtype_oligo, res_QTL_total)
}

# 3.5 OPC
inhouse_OPC <- subset(inhouse_raw, celltype == "OPC")
# within OPC clustering
inhouse_OPC <- FindNeighbors(inhouse_OPC, reduction = "pca", dims = 1:30)
inhouse_OPC <- FindClusters(inhouse_OPC, resolution = 0.05)
# "Mixed" are small clusters of outliers
df_OPC_cluster <- tibble(seurat_clusters=factor(c(0:3)), 
                           subtype=c("OPC_1", "OPC_2", "OPC_Mixed", "OPC_Mixed"))
df_OPC_subtype <- as_tibble(inhouse_OPC@meta.data, rownames = "cell_id") %>% select(cell_id, seurat_clusters)
df_OPC_subtype <- left_join(df_OPC_subtype, df_OPC_cluster, by=c("seurat_clusters"))

inhouse_OPC$subtype <- df_OPC_subtype$subtype

# QTL
df_result_subtype_OPC <- tibble()
for (x in unique(df_OPC_subtype$subtype)) {
  df_work <- filter(df_OPC_subtype, subtype==x)
  df_work <- inner_join(df_cell, df_work, by="cell_id")
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_QTL_total <- mutate(res_QTL_total, subtype = x, n_cell = nrow(df_work), n_ind_cell50 = sum(table(df_work$stim) > 50))
  df_result_subtype_OPC <- bind_rows(df_result_subtype_OPC, res_QTL_total)
}

# merge all subtypes
# Extended Data Fig. 7: subpopulation FAM193B QTL results
df_result_subtype <- bind_rows(df_result_subtype_astro, df_result_subtype_micro, df_result_subtype_neuron,
                               df_result_subtype_oligo, df_result_subtype_OPC)

write_tsv(df_result_subtype, "Table/FAM193B/subpopulation_FAM193B_CCG_QTL.txt")

# prepare clean snRNA dataset
df_cell_subtype <- bind_rows(df_astro_subtype, df_micro_subtype, df_neuron_subtype,
                             df_oligo_subtype, df_OPC_subtype) %>%
  filter(! subtype %in% c("Mic_Mixed", "Oligo_Mixed", "OPC_Mixed", "Neuron_Mixed"))

df_cell <- left_join(df_cell, select(df_cell_subtype, cell_id, subtype), by=c("cell_id"))
inhouse_raw$subtype <- df_cell$subtype
inhouse_clean <- inhouse_raw[, !is.na(inhouse_raw$subtype) | inhouse_raw$celltype == "endo"]

# we swap the order of OPC_2 and Oligo1 for better color visualization
inhouse_clean$subtype <- factor(inhouse_clean$subtype,
                                levels = c("Ast_1", "Ast_2", 
                                           "Excit_L3", "Excit_L3/4/5", "Excit_L4", "Excit_L5", 
                                           "Excit_L5/6", "Excit_L6", "Inhibit", 
                                           "Mic_1", "Mic_2", "OPC_2", "Oligo_2", 
                                           "OPC_1", "Oligo_1"))

# Dim plot
# major cell types
tiff("Plot/FAM193B/Dimplot_clean_celltype_nolabel.tiff", height = 2, width = 2, units = "in", res = 300)
DimPlot(inhouse_clean, group.by = "celltype", reduction="umap", raster = TRUE) + guides(color="none") +
  theme(plot.title = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())
dev.off()
# subtypes
tiff("Plot/FAM193B/Dimplot_clean_subtype_nolabel.tiff", height = 2, width = 2, units = "in", res = 300)
DimPlot(inhouse_clean, group.by = "subtype", reduction="umap", raster = TRUE) + guides(color="none") +
  theme(plot.title = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())
dev.off()
# Feature plot
tiff("Plot/FAM193B/Featureplot_clean_celltype_nolegend.tiff", height = 2, width = 2, units = "in", res = 300)
FeaturePlot(inhouse_clean, "FAM193B", reduction="umap", raster = TRUE, order = T, max.cutoff = 2) + guides(color="none") +
  theme(plot.title = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())
dev.off()

# 4. Major cell type eQTL

# Fig. 5c: association between FAM193B-CCG and FAM193B expression
celltypes <- df_cell$celltype %>% unique
df_result_cell <- tibble()
for (x in celltypes) {
  df_work <- filter(df_cell, celltype == x, !is.na(PC1))
  # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
  mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
  rownames(mat_work) <- c("FAM193B", "tmp")
  colnames(mat_work) <- rownames(df_work)
  mast_work <- FromMatrix(mat_work, cData = df_work)
  # run MAST
  res_QTL_total <- zlm(~total_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                       sca = mast_work[,!is.na(mast_work$total_nCCG)]) %>% MAST_summary(doLRT = "total_nCCG") %>% .[1,]
  res_new <- bind_rows(res_AD, res_QTL_total)
  res_new$celltype <- x
  df_result_cell <- bind_rows(df_result_cell, res_new)
}

write_tsv(df_result_cell, "Table/FAM193B/major_cell_type_FAM193B_CCG_eQTL.txt")

# Extended Data Fig. 6 & Supplementary Data 1: FAM193B cis eQTL in major cell types
snps <- colnames(df_fam193b_snp)[-1]
df_result_eQTL <- tibble()
for (x in celltypes) {
  for (y in snps) {
    df_work <- filter(df_cell, celltype == x, !is.na(PC1))
    # merge SNP genotype
    df_work <- left_join(df_work, select(df_fam193b_snp, ID, all_of(y)), by=c("stim"="ID"))
    colnames(df_work)[ncol(df_work)] <- "SNP"
    # MAST require at least 2 genes, we duplicate FAM193B expression to do analysis
    mat_work <- as.matrix(cbind(df_work$FAM193B, df_work$FAM193B)) %>% t
    rownames(mat_work) <- c("FAM193B", "tmp")
    colnames(mat_work) <- rownames(df_work)
    mast_work <- FromMatrix(mat_work, cData = df_work)
    # run MAST
    res_QTL_total <- zlm(~SNP + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                         sca = mast_work[,!is.na(mast_work$SNP)]) %>% MAST_summary(doLRT = "SNP") %>% .[1,]
    res_QTL_total$celltype <- x
    res_QTL_total$variable <- y
    df_result_eQTL <- bind_rows(df_result_eQTL, res_QTL_total)
  }
}
colnames(df_result_eQTL)[5] <- "P"
write_tsv(df_result_eQTL, "Table/FAM193B/major_cell_type_cis_eQTL.txt")



# 5. DEG analysis
df_result_deg <- tibble()
celltypes <- inhouse_raw$celltype %>% unique
for (x in celltypes) {
  mat_work <- GetAssayData(subset(inhouse_raw, celltype == x), layer = "data", assay = "RNA")
  covar_work <- tibble(cell_id = colnames(subset(inhouse_raw, celltype == x)))
  covar_work <- left_join(covar_work, df_cell, by="cell_id")
  # filter by percent expression
  mat_work <- mat_work[rowMeans(mat_work>0) > 0.05, ]
  mast_work <- FromMatrix(as.matrix(mat_work), cData = covar_work)
  mast_work <- mast_work[, !is.na(mast_work$max_nCCG)]
  # AD vs NC (statistics for Fig. 6b)
  res_deg_ADNC <- zlm(~Diagnosis + Age + Gender + PMD + detection_rate,
                      sca = mast_work) %>% MAST_summary(doLRT = "DiagnosisAD") %>% as_tibble()
  res_deg_ADNC$variable <- "AD vs NC"
  # association with FAM193B-CCG in NC (Fig. 5d, Fig. 6a-b)
  res_deg_maxCCG_NC_dose <- zlm(~max_nCCG + Age + Gender + PMD + detection_rate + PC1 + PC2 + PC3,
                           sca = mast_work[,mast_work$Diagnosis == "CONTROL" ]) %>% MAST_summary(doLRT = "max_nCCG") %>% as_tibble()
  res_deg_maxCCG_NC_dose$variable <- "max_nCCG_in_NC"
  
  res_new <- bind_rows(res_deg_ADNC, res_deg_maxCCG_NC_dose)
  res_new$celltype <- x
  df_result_deg <- bind_rows(df_result_deg, res_new)
  gc()
}

# multiple test correction within each cell type
colnames(df_result_deg)[5] <- "P"
df_result_deg <- df_result_deg %>% group_by(variable, celltype) %>% mutate(FDR = p.adjust(P, method = "BH")) %>% ungroup()
write_tsv(df_result_deg, "Table/FAM193B/DEG_CCG_dosage_effect.txt")
