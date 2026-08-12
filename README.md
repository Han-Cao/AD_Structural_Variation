# Analysis of common structural variants in Alzheimer's disease


## Table of Contents
- [data_processing](./data_processing): Scripts for processing raw FASTQ/BAM/VCF/PLINK data
- [figure](./figure): Scripts for reproducing figures

Code details and instructions are provided in the `README.md` file within each directory.

## Directory structure

```
.
├── data_processing
│   ├── 01.SV_call
│   ├── 02.joint_call
│   ├── 03.SV_annotation
│   ├── 04.genotyping
│   ├── 05.PRS
│   ├── 06.candidate_locus
│   ├── RNAseq
│   ├── WGS
│   └── harmoniSV
└── figure
    ├── input (source data)
    └── output_expected (reproduced figures)
```

## Summary statistics

Due to GitHub file size limits, VCF files and GWAS summary statistics are available only in the [Zenodo release](https://doi.org/10.5281/zenodo.14851280).

## License
[MIT License](./LICENSE)