options(stringsAsFactors = FALSE)

platform <- "ILLUMINA"
center <- "Novogene"

inpath <- "/path/to/fastp"
output <- "input/read_group.csv"

if(!dir.exists("input")){
    dir.create("input")
}


fq_files <- list.files(inpath, "_1_paired.fq.gz", full.names = TRUE, recursive = TRUE)

df_info <- data.frame(file=fq_files)
df_info$sample <- basename(dirname(df_info$file))
df_info$name <- gsub("_1_paired.fq.gz", "", basename(df_info$file))

name_group <- strsplit(df_info$name, "_")

extract_PU <- function(x){
    n_units <- length(x)
    return(paste(x[n_units-1], x[n_units], sep="."))
}

extract_LIB <- function(x){
    n_units <- length(x)
    return(x[n_units-2])
}

df_info$library <- sapply(name_group, extract_LIB)
df_info$platform_unit <- sapply(name_group, extract_PU)
df_info$RG <- paste0("@RG\\tID:", df_info$name, "\\tPU:", df_info$platform_unit, "\\tSM:", df_info$sample, "\\tLB:", df_info$library, "\\tPL:", platform, "\\tCN:", center)
df_info$updated_on <- as.character(Sys.Date())

if(file.exists(output)){
    df_info_old <- read.csv(output, stringsAsFactors = FALSE)
    df_info <- df_info[!df_info$name %in% df_info_old$name,]
    df_info <- rbind(df_info_old, df_info)
}

write.csv(df_info, output, row.names = FALSE, quote=FALSE)
