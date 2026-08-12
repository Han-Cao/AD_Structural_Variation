# Polygenic risk score analysis

Prepare data
- `01.make_bfile.sh`: convert VCF to PLINK bfile and perform quality control
- `02.GWAS.sh`: run GWAS, this is not used for PRS
- `03.kfold_split.py`: split samples for 10 times repeated 10-fold cross validation

PRS
- `04.discovery_effect.sh`: run GWAS within each fold's training sample to estimate beta
- `05.clumping.sh`: LD clumping within each fold's training sample to select variants
- `06.score.sh`: calculate SNP and SV PRS for each fold's validation sample
- `07.score_SNP_SV.sh`: calculate SNP+SV PRS for each fold's validation sample
- `08.score_SNP_SV_condition.sh`: calculate conditional SNP+SV PRS for each fold's validation sample
