# RNA-seq data analysis

Pre-processing:
- `01.fastp.sh`: quality control and reads trimming/filtering

Quantification:
- `02.star_rsem.sh`: mapping and quantification using STAR and RSEM
- `03.combine_results.R`: combine per-sample quantification into matrix

Post-processing:
- `04.PEER.sh`: calculate PEER factors for eQTL analysis