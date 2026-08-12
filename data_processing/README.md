# Code for genetic data processing

This folder contains the data processing code to generate key files for data analysis and figures. 

The final output files are used as the source data for figures, which are provided in the `figure/input/` folder.

Many scripts in this folder require raw FASTQ/BAM/VCF/PLINK files with sensitive genetic information, which cannot be provided in this repository. Please refer to the data availability section in the manuscript for details on data access.

SV discovery and genotyping:
- `01.SV_call`: SV calling on individual samples using long-read sequencing data
- `02.joint_call`: Joint SV calling
- `04.genotyping`: SV gentoyping using short-read sequencing data
- `harmoniSV`: Toolkit for joint SV calling

SV characerization and annotation:
- `03.SV_annotation`: Annotate SVs with genomic features

Genetic association analysis:
- `05.PRS`: GWAS and polygenic risk score analysis
- `06.candidate_locus`: Genetic analysis of candidate AD risk locus

Sequencing data analysis pipeline
- `WGS`: WGS data analysis
- `RNAseq`: RNA-seq data analysis