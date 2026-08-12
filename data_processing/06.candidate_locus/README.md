# Genetic analysis for candidate locus

Prepare data
- `00.ukb_sample_qc.R`: sample QC for UK Biobank, prepare input files for Regenie
- `01.prepare_pfile_ld.sh`: prepare PLINK pfile with repeat dosage, calculate LD between repeat and SNPs
- `02.estiamte_local_ancestry.sh`: local ancestry estimation for UK Biobank
- `03.summarize_local_ancestry.R`: post-process local ancestry estimats for association analysis

Association analysis
- `04.association_plink.sh`: AD associaiton, conditional analysis, and blood eQTL for HK data
- `05.association_regenie.sh`: AD associaiton, conditional analysis, and local ancestry interaction analysis for UK Biobank
- `06.IMR32_DEG.R`: differential expression analysis in IMR32 cells with FAM193B overexpression

Scripts
- `add_ds_for_ld.py`: add repeat dosage to VCF FORMAT/DS for LD analysis
- `add_ds_for_regression.py`: add repeat dosage to VCF FORMAT/DS for regression analysis