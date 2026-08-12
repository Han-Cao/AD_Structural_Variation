library(readr)
library(ggplot2)

# read TR genotypes, capped at 30 for visualization
df_tr <- read_tsv("input/fig4b_tr_distribution/tr_genotypes.txt", col_names="TR")

if(! dir.exists("output/fig4b_tr_distribution/")) dir.create("output/fig4b_tr_distribution/", recursive = T)
pdf("output/fig4b_tr_distribution/fig4b.tr_distribution.pdf")
ggplot(df_tr, aes(x=TR)) + 
    geom_histogram(binwidth=1, fill='#6799d1', color='black')
dev.off()