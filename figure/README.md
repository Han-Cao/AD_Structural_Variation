# Code to reproduce Figures

This directory contains code and related source data to reproduce the figures and statistical analyses in the manuscript.

## Requirements
The following tools are required to run the code:

```
bcftools=1.15.1
truvari=3.5.0
locuszoom
METASOFT
```

For R scripts, the following packages are required:

```
readr
dplyr
tidyr
ggplot2
ggridges
boot
boot.pval
RNOmni
robustbase
IRanges
pROC
rcompanion
nricens
TeachingDemos
magick
enrichplot

# Note: database version may affect results
# see script headers for the used version of the following databases
org.Hs.eg.db
msigdbr
clusterProfiler
```

## Reproduce Figures

Here, we provide scripts and associated source data (e.g., VCF, BED, summary statistics, and individual-level data if possible, see `input/` for details) to reproduce figures in the manuscripts. If all dependencies are correctly installed, the provided scripts can be executed in this working directory without any modifications. The reproduced figures will be saved in the `output/` directory, along with key statistical results printed to the terminal. We have tested all scripts on Linux and the results are provided in the `output_expected/` directory.

While all box, bar, and scatter plots in the manuscript were generated using [GraphPad Prism 8](https://www.graphpad.com/) to maintain visual consistency, we provide `ggplot2` code here to reproduce equivalent visualizations. Note that aesthetic details such as colors, spacing, and labeling can differ between the the results in `output/` and those in the manuscript.

For analyses requiring restricted individual-level genotype data (see Data Availability section in the manuscript), we supply synthetic example data. These data were synthesized from the original dataset using the R package `synthpop` and can be used to run the analyses and get similar results. Plots generated from synthetic data will be marked with the `.synthpop.pdf` suffix.

For wet lab experiments, we provide the source data here. The statistical analysis and visualizations can be reproduced using [GraphPad Prism 8](https://www.graphpad.com/).

More methodological details are provided in each script.