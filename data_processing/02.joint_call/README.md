# Joint SV calling

Preprocessing:
- `01.harmonize.sh`: harmonize VCFs from different SV callers to a standard format
- `02.jasmine_preprocess.sh`: identify duplicated SV calls within the same sample using Jasmine
- `03.remove_duplicate_call.sh`: remove duplicated SV calls identified by keeping the SV with the largest number of supporting reads
- `06.harmonize_force_calling.sh`: similar to `01.harmonize.sh`, but harmonize the force calling VCFs from `01.SV_call/04.force_calling.sh`

SV merging:
- `04.merge_sv`: merge SV calls from different callers and samples using Jasmine
- `05.representative_sv.sh`: identify the representative SV among merged SVs by selecting the one with most frequency POS and SVLEN
- `08.merge_with_ref.sh`: merge SVs with public SV dataset to classify TP and FP calls for model training
- `12.merge_filtered_with_ref.sh`: merge high-confidence SVs with public SV dataset and CHM13 SVs to evaluate robustness, this include additional SV dataset suggested by the reviewers

Postprocessing:
- `07.genotype_sample.sh`: per-sample ensemble genotyping using the results from all SV callers
- `09.filterSV.sh`: random forest model to filter high-confidence SVs
- `10.prepare_bench.sh`: prepare SV call set of HG002 for benchmark
- `11.merge_sample.sh`: convert sample-level VCF to population-level VCF