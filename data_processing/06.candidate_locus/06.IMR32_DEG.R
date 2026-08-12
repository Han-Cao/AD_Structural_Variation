library(tximport)
library(edgeR)
library(apeglm)
library(DESeq2)
library(dplyr)

###### SETTINGS ######
# set the input and output directories
sample_tab <- "input/RNAseq_samples.txt"
id_map <- "input/Homo_sapiens.GRCh38.109.id_map.txt"

###### ANALYSIS ######
# read tables for samples and genes
df_sample <- read_tsv(sample_tab)
df_map_isoform <- read_tsv(id_map, col_names = FALSE)
colnames(df_map_isoform) <- c("Gene", "Isoform", "Symbol")
df_map_gene <- dplyr::select(df_map_isoform, Gene, Symbol) %>% dplyr::filter(!duplicated(Gene)) # remove duplicated records at gene level

# load gene expression matrix
rsem_files <- df_sample$file
names(rsem_files) <- df_sample$sample
txi <- tximport(rsem_files, type = "rsem", txIn = FALSE, txOut = FALSE)
# technically, RSEM produce effect gene length of 0 when gene length < read length
# it can be converted to 1 to allow DESeq2 estimate library size
txi$length[txi$length == 0] <- 1

# convert expression matrix to table and map gene ID
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

norm_cpm <- function(txi){
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

        cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)
        return(cpms)
}

df_rsem_gene_counts <- convert_expr_mat(txi$counts, "Gene", df_map_gene)
df_rsem_gene_tpm <- convert_expr_mat(txi$abundance, "Gene", df_map_gene)
df_rsem_gene_cpm <- norm_cpm(txi) %>% convert_expr_mat("Gene", df_map_gene)

write_tsv(df_rsem_gene_counts, "/path/to/gene_expression_matrix_counts.txt")
write_tsv(df_rsem_gene_tpm, "/path/to/gene_expression_matrix_tpm.txt")
write_tsv(df_rsem_gene_cpm, "/path/to/gene_expression_matrix_cpm.txt")

# extract genes with at least half sample with TPM > 1
expr_mat <- select(df_rsem_gene_tpm, OE_1:PC_6) %>% as.matrix()
keep_idx <- rowMeans(expr_mat > 1) >= 0.5
keep_genes <- df_rsem_gene_tpm$Gene[keep_idx]


# DEG analysis
# force control as reference
df_sample$group <- factor(df_sample$group, levels = c("PC", "OE"))
# load txi to DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = df_sample, design = ~ group + batch)
# run analysis
dds <- DESeq(dds)
# DEG results
contrast <- c("group", "OE", "PC")
res <- results(dds, contrast = contrast, alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "group_OE_vs_PC", type = "apeglm")

# convert the results to table and map gene names
df_res <- as_tibble(res, rownames = "Gene")
df_resLFC <- as_tibble(resLFC, rownames = "Gene")

df_res <- right_join(df_map_gene, df_res, by=c("Gene"))
df_resLFC <- right_join(df_map_gene, df_resLFC, by=c("Gene"))

# filter out genes with very low expression
df_res <- filter(df_res, Gene %in% keep_genes)
df_resLFC <- filter(df_resLFC, Gene %in% keep_genes)

df_res$padj <- p.adjust(df_res$pvalue, method = "BH")
df_resLFC$padj <- p.adjust(df_resLFC$pvalue, method = "BH")

write_tsv(df_res, "/path/to/DEG_result.txt")
write_tsv(df_resLFC, "/path/to/DEG_result_lfcShrink.txt")
