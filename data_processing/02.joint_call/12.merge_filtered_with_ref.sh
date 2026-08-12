#!/bin/bash

# compare with external reference SV call sets
# as suggested by reviewers, we additionally added CCDG and 1KGP_ONT dataset

# 1. merge all SVs
bcftools view -o /path/to/Joint_call/benchmark/Chinese_SV.filter.vcf \
/path/to/Joint_call/population/All.genotype.filter.population.AC0.vcf.gz

echo '/path/to/Joint_call/benchmark/Chinese_SV.filter.vcf
/path/to/reference/jasmine_input/1KGP.vcf
/path/to/reference/jasmine_input/Chinese_405.vcf
/path/to/reference/jasmine_input/gnomAD.vcf
/path/to/reference/jasmine_input/HGDP.vcf
/path/to/reference/jasmine_input/HGSVC.vcf
/path/to/reference/jasmine_input/Icelander.vcf
/path/to/reference/jasmine_input/LRS15.vcf
/path/to/reference/jasmine_input/CCDG.normalize.vcf
/path/to/reference/jasmine_input/chm13.svim-asm.set_id.vcf' > \
/path/to/Joint_call/benchmark/merge_vcf.txt

jasmine \
file_list=/path/to/Joint_call/benchmark/merge_vcf.txt \
out_file=/path/to/Joint_call/benchmark/Chinese_SV.filter.merge_with_ref.vcf \
threads=40 \
--keep_var_ids \
--ignore_strand 

# 2. Extract SVs and frequencies of this study
bcftools query -f '%IDLIST\n' /path/to/Joint_call/benchmark/Chinese_SV.filter.merge_with_ref.vcf | \
grep 'Chinese_SV' > /path/to/Joint_call/benchmark/Chinese_SV.filter.merge_with_ref.txt

bcftools +fill-tags /path/to/Joint_call/population/All.genotype.filter.population.AC0.vcf.gz -- -t AC,AN | \
bcftools query -f '%ID\t%AC\n'  > /path/to/Joint_call/benchmark/Chinese_SV.filter.AC.txt