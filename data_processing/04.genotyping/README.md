# Genotyping of short-read sequencing data

SV genotyping
- `01.paragraph.sh`: per-sample SV genotyping using Paragraph
- `02.merge_vcf.sh`: combine per-sample VCF to population VCF
- `03.genotype_concordance.sh`: analyze genotype concordance of overlapping samples between LRS and SRS data

TR genotyping
- `04.tandem_repeat.sh`: tandem repeat genotyping

Benchmark
- `05.benchmark.sh`: benchmark SV genotyped by Paragraph using 30x PacBio HiFi data of one sample