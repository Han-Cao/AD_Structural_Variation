# SV calling

Read alignment:
- `01.align_minimap2.sh`: align long-read sequencing data to GRCh38 using minimap2
- `01.align_NGMLR.sh`: align long-read sequencing data to GRCh38 using NGMLR

SV calling:
- `02.cuteSV.sh`: call SVs using cuteSV
- `02.Sniffles.sh`: call SVs using Sniffles
- `02.SVIM.sh`: call SVs using SVIM
- `03.filter_SV.py`: script to filter raw SV calls by length and quality
- `03.filter_SV_run.sh`: script to run filter_SV.py for all SV calling methods
- `05.CHM13_sv_call.sh`: call SVs on CHM13 comapred to GRCh38 using SVIM-asm

SV force calling:
- `04.force_calling.sh`: SV force calling using cuteSV and Sniffles, this is run after merging all SVs (`02.joint_call/04.merge_sv.sh`) and identifying representative SVs (`02.joint_call/05.representative_sv.sh`)
