#!/bin/bash
#SBATCH -p q1 -N 1 -n 20 -c 2 -o log/hardfilter_concat.log

fai="/path/to/GRCh38.fa.fai"

inpath="/path/to/Joint-call/"
# interval list is available at https://github.com/Han-Cao/Genome_interval_list
interval_folder="input/interval_list_bcftools/"

# extract variants by POS within the interval
parallel -j $SLURM_NTASKS '
vcf={1}
interval_folder={2}
output=${vcf%.raw.vcf.gz}.hardfilter.tmp.vcf.gz

name=$(basename $vcf)
name=${name%.raw.vcf.gz}

interval_list=${interval_folder}/${name}.interval_list

bcftools filter \
-R $interval_list --regions-overlap 0 \
-s "ExcessHet" -e "INFO/ExcessHet > 54.69" \
-O z -o $output $vcf
tabix $output
' ::: ${inpath}/*raw.vcf.gz ::: $interval_folder

for x in {000..222}
do
	echo ${inpath}/interval_${x}.hardfilter.tmp.vcf.gz
done > ${inpath}/vcf-merge.list

bcftools concat -n -Oz -o ${inpath}/HK_WGS_update.hardfilter.vcf.gz -f ${inpath}/vcf-merge.list
bcftools reheader -f $fai -o ${inpath}/HK_WGS_update.hardfilter.reheader.vcf.gz ${inpath}/HK_WGS_update.hardfilter.vcf.gz
mv ${inpath}/HK_WGS_update.hardfilter.reheader.vcf.gz ${inpath}/HK_WGS_update.hardfilter.vcf.gz
tabix -p vcf ${inpath}/HK_WGS_update.hardfilter.vcf.gz