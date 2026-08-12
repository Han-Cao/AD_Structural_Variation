#!/bin/bash

#SBATCH -p cpu -N 1 -n 40

inpath="/path/to/SV_call/"
outpath="/path/to/SV_call_hardfilter/"

# sniffles
# keep SVs <= 1MB
ls ${inpath}/sniffles/*vcf | parallel -j 10 '
in_vcf={1}
outpath={2}

name=$(basename $in_vcf)
out_vcf=${outpath}/${name}.gz
prefix=${name%.sniffles.vcf}

python 03.filter_SV.py \
-i $in_vcf \
-o $out_vcf \
--prefix $prefix \
--mode Sniffles \
--maxlen 1000000 

' :::: - ::: $outpath

# svim
# keep by QUAL>=4 (corresponding to RE>=4), SVs <= 1MB
ls ${inpath}/svim/ | parallel -j 10 '
sample={1}
inpath={2}
outpath={3}

prefix=$(basename $sample)
in_vcf=${inpath}/svim/${prefix}/variants.vcf
out_vcf=${outpath}/${prefix}.svim.vcf.gz

python 03.filter_SV.py \
-i $in_vcf \
-o $out_vcf \
--prefix $prefix \
--mode SVIM \
--maxlen 1000000 \
--qual 4

' :::: - ::: $inpath ::: $outpath


# cuteSV
# keep SVs <= 1MB
ls ${inpath}/cuteSV/*vcf | parallel -j 10 '
in_vcf={1}
outpath={2}

name=$(basename $in_vcf)
out_vcf=${outpath}/${name}.gz
prefix=${name%.cuteSV.vcf}

python 03.filter_SV.py \
-i $in_vcf \
-o $out_vcf \
--prefix $prefix \
--mode cuteSV \
--maxlen 1000000 

' :::: - ::: $outpath
