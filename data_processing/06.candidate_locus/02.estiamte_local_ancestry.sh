#!/bin/bash

vcf_1kb=/path/to/1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
vcf_ukb="/path/to/small_variants.vcf.gz"
sample_1kg=/path/to/1kGP.unrelated_2504.txt

outpath="/path/to/local_ancestry/"

# prepare 1000G reference file
# keep only SNPs following HWE, set unique ID
bcftools view -r chr5:177322587-177754013 -v snps -S $sample_1kg -Ou $vcf_1kb | \
bcftools +fill-tags -Ou -- -t HWE | \
bcftools view -i 'HWE > 1e-15' -Ou | \
bcftools annotate --set-id '%CHROM:%POS\_%REF\_%ALT' -Oz -o ${outpath}/1000G.FAM193B.vcf.gz
tabix ${outpath}/1000G.FAM193B.vcf.gz
bcftools query -f '%ID\n' ${outpath}/1000G.FAM193B.vcf.gz > ${outpath}/1000G.FAM193B.var_id.txt

# prepare UK Biobank VCF
bcftools view -f PASS -S /path/to/UKB.clean.all_ancestry_samples.txt -Ou $vcf_ukb | \
bcftools +fill-tags -Ou -- -t MAF,HWE | \
bcftools view -i 'MAF > 0.01 & HWE > 1e-15' -Ou | \
bcftools annotate --set-id '%CHROM:%POS\_%REF\_%ALT' -Oz -o ${outpath}/UKB.FAM193B.vcf.gz
tabix ${outpath}/UKB.FAM193B.vcf.gz
bcftools query -f '%ID\n' ${outpath}/UKB.FAM193B.vcf.gz > ${outpath}/UKB.FAM193B.var_id.txt

# only keep intersecting variants
bcftools view -i "ID=@${outpath}/1000G.FAM193B.var_id.txt" \
-Oz -o ${outpath}/UKB.FAM193B.intersect_1kg.vcf.gz ${outpath}/UKB.FAM193B.vcf.gz
bcftools view -i "ID=@${outpath}/UKB.FAM193B.var_id.txt" \
-Oz -o ${outpath}/1000G.FAM193B.intersect_ukb.vcf.gz ${outpath}/1000G.FAM193B.vcf.gz

# phasing with beagle5
java -jar ~/bin/beagle.17Dec24.224.jar \
gt=${outpath}/UKB.FAM193B.intersect_1kg.vcf.gz \
ref=${outpath}/1000G.FAM193B.intersect_ukb.vcf.gz \
out=${outpath}/UKB.FAM193B.phased.beagle \
map=/path/to/beagle/plink.chr5.GRCh38.map \
chrom=chr5 \
burnin=5 \
iterations=40 \
nthreads=40

# prepare genetic map file from 1000G VCF
echo -e "position(base)\tgenetic_distance(centiMorgan)" > ${outpath}/UKB.FAM193B.map.txt
bcftools query -f '%POS\t%INFO/CM\n' ${outpath}/1000G.FAM193B.intersect_ukb.vcf.gz >> ${outpath}/UKB.FAM193B.map.txt

# prepare sample name file in chunks
[[ ! -d ${outpath}/UKB.FAM193B.local_ancestry_chunks ]] && mkdir -p ${outpath}/UKB.FAM193B.local_ancestry_chunks
bcftools query -l ${outpath}/UKB.FAM193B.phased.beagle.vcf.gz | \
split -l 1000 - -d --additional-suffix .txt ${outpath}/UKB.FAM193B.local_ancestry_chunks/chunk_

# local ancestry estiamte
cd /path/to/SparsePainter
ls ${outpath}/UKB.FAM193B.local_ancestry_chunks/chunk_*.txt | while read chunk;
do
    chunk_id=$(basename $chunk)
    chunk_id=${chunk_id%.txt}

    ./SparsePainter \
    -reffile ${outpath}/1000G.FAM193B.intersect_ukb.vcf.gz \
    -targetfile ${outpath}/UKB.FAM193B.phased.beagle.vcf.gz \
    -popfile /path/to/1000G.sample_pop.txt \
    -mapfile ${outpath}/UKB.FAM193B.map.txt \
    -namefile ${chunk} \
    -out ${outpath}/UKB.FAM193B.local_ancestry_chunks/${chunk_id}.local_ancestry \
    -prob -chunklength -chunkcount -aveSNP -aveind -sample \
    > ${outpath}/UKB.FAM193B.local_ancestry_chunks/${chunk_id}.local_ancestry.log 2>&1

    # postprocess prob results
    # the code provided by SparsePainter not work well... we use awk to do this
    zcat ${outpath}/UKB.FAM193B.local_ancestry_chunks/${chunk_id}.local_ancestry_prob.txt.gz | tail -n +2 | \
    awk '
    NR==1 {print "sample", "SNPidx", $3, $4, $5, $6, $7}
    NF==1 {sample=$1; next} 
    NF>=3 {
        for(snp=$1; snp<=$2; snp++) {
            printf "%s %d", sample, snp;
            for(i=3; i<=NF; i++) printf " %s", $i;
            print ""
        }
    }' | gzip -c > ${outpath}/UKB.FAM193B.local_ancestry_chunks/${chunk_id}.local_ancestry_prob_full.txt.gz
done

# merge results
# per-sample local ancestry
echo "individual_name EUR AFR EAS AMR SAS" > ${outpath}/UKB.FAM193B.local_ancestry.aveindprob.txt
cat ${outpath}/UKB.FAM193B.local_ancestry_chunks/*aveindprob.txt | \
grep -v individual_name >> ${outpath}/UKB.FAM193B.local_ancestry.aveindprob.txt

# per-SNP local ancestry
echo "sample SNPidx EUR AFR EAS AMR SAS" > ${outpath}/UKB.FAM193B.local_ancestry.prob_full.txt
zcat ${outpath}/UKB.FAM193B.local_ancestry_chunks/*prob_full.txt.gz | \
grep -v sample >> ${outpath}/UKB.FAM193B.local_ancestry.prob_full.txt
gzip ${outpath}/UKB.FAM193B.local_ancestry.prob_full.txt

# prepare map file used in downstream analysis
echo -e "position(base)\tgenetic_distance(centiMorgan)" > ${outpath}/UKB.FAM193B.map.txt
bcftools query -f '%POS\t%INFO/CM\n' ${outpath}/1000G.FAM193B.intersect_ukb.vcf.gz >> ${outpath}/UKB.FAM193B.map.txt
