#!/bin/bash
#
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -o log/PEER.log

inpath=$1
input_matrix=${inpath}/RSEM_gene_TMM.clean.ranknorm.txt
outpath=${inpath}/PEER/

[[ ! -d $outpath ]] && mkdir -p $outpath

# transpose matrix
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $input_matrix | tr ' ' '\t' | tail -n +3 | cut -f 2- > ${outpath}/RSEM_gene_TMM.clean.ranknorm.tab

# run peertool
peertool \
-f ${outpath}/RSEM_gene_TMM.clean.ranknorm.tab \
-n 60 \
-o $outpath \
-i 1000
