library(readr)
library(tidyr)
library(dplyr)

# read data
df_vep <- read_tsv("output/Chinese_SV.VEP_output.txt.gz", 
                   comment = "##", col_types = list(Location = "c"))

# parse vep annotation file
grep_anno <- function(anno, tag){
    re_expr <- paste0("(.+;|^)", tag, "=(.+?)($|;.*)")
    ifelse(grepl(paste0(tag, "="), anno), gsub(re_expr, "\\2", anno), NA)
}

# extract vep annotation
colnames(df_vep)[1] <- "ID"
df_vep <- mutate(df_vep,
                 IMPACT = grep_anno(Extra, "IMPACT"),
                 SYMBCOL = grep_anno(Extra, "SYMBOL"),
                 BIOTYPE = grep_anno(Extra, "BIOTYPE"),
                 CANONICAL = grep_anno(Extra, "CANONICAL"),
                 MANE = grep_anno(Extra, "MANE_SELECT"),
                 TSL = grep_anno(Extra, "TSL") %>% as.numeric(),
                 Overlap = grep_anno(Extra, "OverlapBP") %>% as.numeric(),
                 ccRE = grep_anno(Extra, "ccRE"))

# summarise transcirpt BIOTYPE
df_vep <- mutate(df_vep,
                   protein_coding = grepl("protein_coding", BIOTYPE),
                   small_RNA = grepl("(mi|sn|sno|sca|s)RNA", BIOTYPE),
                   lncRNA = grepl("lncRNA", BIOTYPE),
                   pseudogene = grepl("pseudogene", BIOTYPE),
                   BIOTYPE_key = ifelse(protein_coding, "protein_coding",
                                 ifelse(small_RNA, "small_RNA",
                                 ifelse(lncRNA, "lncRNA",
                                 ifelse(pseudogene, "pseudogene", NA)))))

# protein coding SVs
df_protein_coding <- filter(df_vep, BIOTYPE_key == "protein_coding", !is.na(MANE))

# identify predicted loss of function SVs
df_LOF <- filter(df_protein_coding)
df_LOF <- mutate(df_LOF,
                 inframe_insertion = Allele %in% c("insertion", "duplication") & Protein_position != '-',
                 inframe_deletion = Allele %in% c("deletion", "INV") & Protein_position != '-',
                 stop_lost = grepl("stop_lost", Consequence),
                 transcript_ablation = grepl("transcript_ablation", Consequence),
                 transcript_amplification = grepl("transcript_amplification", Consequence))
df_LOF <- filter(df_LOF, inframe_insertion | inframe_deletion | stop_lost | transcript_ablation | transcript_amplification)

write_tsv(df_LOF, "output/SV_annotation_protein_coding_MANE_LoF.txt")
