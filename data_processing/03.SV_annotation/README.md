# SV annotation

Annotate genomic features
- `00.generate_bed_files.sh`: download and generate genomic annotation bed files
- `01.overlap_feature.sh`: intersect SVs with genomic features using bedtools
- `02.summarise_feature_overlap.R`: count the length of features and overlap between features and SVs
- `03.prepare_VEP_input.R`: convert VCF to VEP input format
- `04.VEP_run.sh`: run VEP on SVs
- `05.summarise_VEP.R`: count SVs overlapping with genes grouped by biotypes, identify predicted LoF SVs

Annotate GWAS loci
- `02.overlap_GWAS.sh`: intersect GWAS loci with SVs using bedtools
