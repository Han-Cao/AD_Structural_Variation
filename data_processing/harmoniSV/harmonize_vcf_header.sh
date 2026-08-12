#!/bin/bash

if [[ $# -eq 0 ]] ; then
	echo "Harmonize vcf INFO/FORMAT/FILTER headers. Other headers will use the first input vcf header."
	echo ""
	echo "Usage:"
	echo "-f input vcf list file"
	echo "-o output header"
  echo "-v reference vcf (optional)"
  exit 0
fi

while getopts ":f:o:v:" opt; do
  case $opt in
    f) input_file="$OPTARG"
    ;;
    o) output="$OPTARG"
	  ;;
    v) reference="$OPTARG"
  esac
done


if [[ ! -z $reference ]]; then
  bcftools view -h $reference |  grep '^##INFO=' > ${output}.tmp.INFO
  bcftools view -h $reference |  grep '^##FORMAT=' > ${output}.tmp.FORMAT
  bcftools view -h $reference |  grep '^##FILTER=' > ${output}.tmp.FILTER
fi

cat $input_file | while read vcf;
do
    bcftools view -h $vcf |  grep '^##INFO=' >> ${output}.tmp.INFO
    bcftools view -h $vcf |  grep '^##FORMAT=' >> ${output}.tmp.FORMAT
    bcftools view -h $vcf |  grep '^##FILTER=' >> ${output}.tmp.FILTER
done

if [[ ! -z $reference ]]; then
  VCF1=$reference
else
  VCF1=$(head -n 1 $input_file)
fi
bcftools view -h $VCF1 | head -n -1 | zgrep -v -E '^##(INFO|FORMAT|FILTER)=' > ${output}
sort ${output}.tmp.INFO | uniq >> ${output}
sort ${output}.tmp.FORMAT | uniq >> ${output}
sort ${output}.tmp.FILTER | uniq >> ${output}
bcftools view -h $VCF1 | tail -n 1 >> ${output}

[[ -e ${output}.tmp.INFO ]] && rm ${output}.tmp.INFO
[[ -e ${output}.tmp.FORMAT ]] && rm ${output}.tmp.FORMAT
[[ -e ${output}.tmp.FILTER ]] && rm ${output}.tmp.FILTER

