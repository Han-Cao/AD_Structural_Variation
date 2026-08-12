# Illumina WGS data analysis

Pre-processing:
- `01.fastp.sh`: quality control and reads trimming/filtering
- `02.extract_read_group.R`: prepare read group from fastq files

Variant calling:
- `03.align.sh`: align reads to GRCh38
- `04.markdup.sh`: mark duplicates
- `05.BQSR.sh`: base quality recalibration
- `07.HaplotypeCaller_applyBQSR.sh`: apply BQSR and call variants
- `08.GDBimport.sh`: import GVCFs to GenomicsDB
- `09.Joint-call.sh`: joint genotyping

Sample QC:
- `06.check_contaminaiton.sh`: check cross-sample contamination
- `14.sample_qc_make_bfile.sh`: basic sample QC and convert VCF to bfile
- `14.sample_qc_check_sex.sh`: compare genetic sex with phenotypic sex
- `14.sample_qc_PCA_relatedness.R`: check relatedness and run PCA

Variant QC:
- `10.hardfilter_concat.sh`: hard filter by ExcessHet and concat VCFs from genomic intervals
- `11.VQSR_SNP.sh`: variant quality score recalibration for SNPs
- `12.VQSR_INDEL.sh`: variant quality score recalibration for INDELs
- `13.clean_VQSR.sh`: apply VQSR filtering