library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(gridExtra)
library(tximport)
library(edgeR)
library(RNOmni)

inpath <- "/path/to/rsem/output/"
outpath <- "/path/to/output/"
outpath <- paste0(outpath, gsub("-", "", Sys.Date()), "/")
outpath_plot <- paste0(outpath, "plot/")

if(! dir.exists(outpath)){
    dir.create(outpath, recursive = TRUE)
}

if(! dir.exists(outpath_plot)){
    dir.create(outpath_plot, recursive = TRUE)
}

# set input
dir_rsem_abundance <- paste0(inpath, "rsem_abundance/")
dir_rsem_stat <- paste0(inpath, "rsem_stat/")

# extract sample ID from file path
extract_name <- function(x, patten_suffix){
    y <- basename(x)
    gsub(patten_suffix, "", y)
}

# convert expression matrix to tibble and map gene ID
convert_expr_mat <- function(x, type, df_map){
    if(type == "Gene"){
        df_new <- as_tibble(x, rownames="Gene")
        df_new <- right_join(df_map, df_new, by = "Gene")
    }
    else if (type=="Isoform") {
       df_new <- as_tibble(x, rownames="Isoform")
       df_new <- right_join(df_map, df_new, by="Isoform")
    }
    return(df_new)
}

# gene ID map
df_map_isoform <- read_tsv("/path/to/Homo_sapiens.GRCh38.109.id_map.txt", col_names = FALSE)
colnames(df_map_isoform) <- c("Gene", "Isoform", "Symbol")
df_map_gene <- select(df_map_isoform, Gene, Symbol) %>% filter(!duplicated(Gene))

############### rsem quantification
# gene-level
cat("Processing RSEM gene-level quantificaitons\n")
file_rsem_gene <- list.files(dir_rsem_abundance, ".+genes\\.results", full.names = TRUE)
names(file_rsem_gene) <- extract_name(file_rsem_gene, "\\.rsem\\.genes\\.results")
if(file.exists(paste0(outpath, "RSEM_gene_tximport.rds"))){
    txi_rsem_gene <- readRDS(paste0(outpath, "RSEM_gene_tximport.rds"))
} else {
    txi_rsem_gene <- tximport(file_rsem_gene, type = "rsem", txIn = FALSE, txOut = FALSE)
    saveRDS(txi_rsem_gene, paste0(outpath, "RSEM_gene_tximport.rds"))
}
df_rsem_gene_counts_raw <- convert_expr_mat(txi_rsem_gene$counts, "Gene", df_map_gene)
df_rsem_gene_tpm_raw <- convert_expr_mat(txi_rsem_gene$abundance, "Gene", df_map_gene)

write_tsv(df_rsem_gene_counts_raw, paste0(outpath, "RSEM_gene_count.raw.txt"))
write_tsv(df_rsem_gene_tpm_raw, paste0(outpath, "RSEM_gene_TPM.raw.txt"))

# isoform-level
cat("Processing RSEM isoform-level quantificaitons\n")
file_rsem_isoform <- list.files(dir_rsem_abundance, ".+isoforms\\.results", full.names = TRUE)
names(file_rsem_isoform) <- extract_name(file_rsem_isoform, "\\.rsem\\.isoforms\\.results")
if(file.exists(paste0(outpath, "RSEM_isoform_tximport.rds"))){
    txi_rsem_isoform <- readRDS(paste0(outpath, "RSEM_isoform_tximport.rds"))
} else {
    txi_rsem_isoform <- tximport(file_rsem_isoform, type = "rsem", txIn = TRUE, txOut = TRUE)
    saveRDS(txi_rsem_isoform, paste0(outpath, "RSEM_isoform_tximport.rds"))
}
df_rsem_isoform_counts_raw <- convert_expr_mat(txi_rsem_isoform$counts, "Isoform", df_map_isoform)
df_rsem_isoform_tpm_raw <- convert_expr_mat(txi_rsem_isoform$abundance, "Isoform", df_map_isoform)

write_tsv(df_rsem_isoform_counts_raw, paste0(outpath, "RSEM_isoform_count.raw.txt"))
write_tsv(df_rsem_isoform_tpm_raw, paste0(outpath, "RSEM_isoform_TPM.raw.txt"))

############### rsem stat
cat("Processing RSEM stat\n")
file_rsem_stat <- list.files(dir_rsem_stat, ".+rsem\\.cnt", full.names = TRUE)
names(file_rsem_stat) <- extract_name(file_rsem_stat, "\\.rsem\\.cnt")

if(file.exists(paste0(outpath, "RSEM_stat.txt"))){
    df_rsem_stat <- read_tsv(paste0(outpath, "RSEM_stat.txt"))
} else {
    list_rsem_stat <- lapply(file_rsem_stat, scan, what="numeric", nlines=3)
    df_rsem_stat <- do.call(rbind, list_rsem_stat) %>% as_tibble(rownames="ID") 
    colnames(df_rsem_stat) <- c("ID", "nUnalign", "nAligned", "nFilter", "nTotal",
                                "nUnique", "nMulti", "nUncertain", 
                                "nHits", "read_type")
    df_rsem_stat <- mutate(df_rsem_stat, across(-ID, as.integer))

    write_tsv(df_rsem_stat, paste0(outpath, "RSEM_stat.txt"))
    rm(list_rsem_stat)
}

############### normalization
norm_cpm <- function(txi, keep_rows){
        cts <- txi$counts
        normMat <- txi$length

        # rsem output 0 effective length if <=1, reset to 1 for library estimation
        normMat[normMat==0] <- 1

        # Obtaining per-observation scaling factors for length, adjusted to avoid
        # changing the magnitude of the counts.
        normMat <- normMat/exp(rowMeans(log(normMat)))
        normCts <- cts/normMat

        # RSEM produce 0 effective length

        # Computing effective library sizes from scaled counts, to account for
        # composition biases between samples.
        eff.lib <- calcNormFactors(normCts) * colSums(normCts)

        # Combining effective library sizes with the length factors, and calculating
        # offsets for a log-link GLM.
        normMat <- sweep(normMat, 2, eff.lib, "*")
        normMat <- log(normMat)

        # Creating a DGEList object for use in edgeR.
        y <- DGEList(cts)
        y <- scaleOffset(y, normMat)

        y <- y[keep_rows, ]
        cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)
        return(cpms)
}

# txi_rsem_gene <- readRDS(paste0(outpath, "RSEM_gene_tximport.rds"))
keep_gene_tpm_idx <- txi_rsem_gene$abundance %>% 
                     apply(., 1, function(x) (sum(x > 0.1) / length(x)) > 0.2)
keep_gene_count_idx <- txi_rsem_gene$counts %>% 
                       apply(., 1, function(x) (sum(x >= 6) / length(x)) > 0.2)
keep_gene_idx <- keep_gene_tpm_idx & keep_gene_count_idx

rsem_gene_tmm <- norm_cpm(txi_rsem_gene, keep_gene_idx)
rsem_gene_tmm_ranknorm <- apply(rsem_gene_tmm, 1, RankNorm) %>% t

df_rsem_gene_clean_tmm <- convert_expr_mat(rsem_gene_tmm, "Gene", df_map_gene)
df_rsem_gene_clean_tmm_ranknorm <- convert_expr_mat(rsem_gene_tmm_ranknorm, "Gene", df_map_gene)
write_tsv(df_rsem_gene_clean_tmm, paste0(outpath, "RSEM_gene_TMM.clean.txt"))
write_tsv(df_rsem_gene_clean_tmm_ranknorm, paste0(outpath, "RSEM_gene_TMM.clean.ranknorm.txt"))
