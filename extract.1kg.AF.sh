#!/bin/bash

set -o nounset

# set number of processes
nproc=15

# get AF of the 5 superpopulations from the 1kg vcfs

# for each chromosome, extract the AF information from the vcf
# only do this for the sites in GTEx

genomes_1000_dir="$1"
variant_calling_dir="$2"


frq_file=$variant_calling_dir"merge_sample_names_SNPs.frq"

for i in $(seq 1 22); do echo $i; done | \
parallel --jobs $nproc "vcftools --gzvcf ${genomes_1000_dir}ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out ${variant_calling_dir}chr{}.SNPs --remove-indels --positions $frq_file --get-INFO EAS_AF --get-INFO AMR_AF --get-INFO AFR_AF --get-INFO EUR_AF --get-INFO SAS_AF"


cat ${variant_calling_dir}*SNPs.INFO > ${variant_calling_dir}concatenated_SNPs.INFO

python process.1kg.AF.py ${variant_calling_dir}concatenated_SNPs.INFO | sort -k1,1 -k2,2n > ${variant_calling_dir}SNPs.1kg.AF.bed


rm ${variant_calling_dir}*.INFO
gzip ${variant_calling_dir}*.bed
