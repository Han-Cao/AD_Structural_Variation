#!/bin/bash
#SBATCH -p himem-share -N 1 -n 11 -c 2

# per-sample SV genotyping
parallel -j $SLURM_NTASKS '
sample={1}
outpath="/path/to/Joint_call/genotype_sample/"

srun -N 1 -n 1 -o log/genotype/${sample}.log \
python /path/to/code/harmoniSV/genotypeSV.py \
-i /path/to/represent/All_sample_pipeline_merged.representative.all.vcf.gz \
-o ${outpath}/${sample}.genotype.vcf.gz \
--sample ${sample} \
--method-table input/method_table.txt \
--sv-info input/discovery_table.txt.gz \
--force-call-info input/force_call_table.txt.gz \
--genotyping-method force_call
' :::: input/samples.txt

# concat all samples
cat input/samples.txt | while read sample
do
    echo "/path/to/Joint_call/genotype_sample/${sample}.genotype.vcf.gz"
done > /path/to/Joint_call/genotype_sample/All_vcf.txt

bcftools concat -a -f /path/to/Joint_call/genotype_sample/All_vcf.txt \
-Oz -o /path/to/Joint_call/genotype_sample/All.genotype.vcf.gz
tabix /path/to/Joint_call/genotype_sample/All.genotype.vcf.gz
